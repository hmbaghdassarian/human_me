#!/usr/bin/env python
# coding: utf-8

import copy
import gc
import logging
import math
import multiprocessing
import os
import pathlib
import pickle
from typing import Dict, List, Optional, SupportsFloat, Tuple, Union
import warnings

import cobra
import numpy as np
import pandas as pd
import sympy
from cobra.core.dictlist import DictList
from six import iteritems, string_types
from sympy import lambdify
from tqdm import tqdm

from human_me.core.gene import ExpressedGene
from human_me.core.macromolecules.complex import Complex, RibosomalComplex
from human_me.core.macromolecules.macromolecule import Macromolecule
from human_me.core.reaction import (BiomassReaction,
                                    ComplexDegradationReaction,
                                    ExpressionReaction, MetabolicReaction)
from human_me.me_solver import solve_me
from human_me.preprocess import parse_complex
from human_me.utils import parameters as params
from human_me.utils.functions import flatten_list
from human_me.utils.machinery import rbps

logger = logging.getLogger(__name__)


class ME_Model(cobra.Model):
    # rewritten methods------------------------------------------------------------------------------------------------
    def __init__(self, id_or_model, name: Optional[str] = None, m_model: Optional[cobra.Model] = None,
                 n_cores: int = os.cpu_count(),
                 non_machinery: Dict[str, List[str]] = None, knock_out: Optional[List[str]] = None, additional_ko: List[str] = None):
        """[summary]

        Parameters
        ----------
        id_or_model : [type]
            [description]
        name : str, optional
            model name, by default None
        m_model : cobra.Model, optional
            The cobrapy model object that the ME_Model was built from. Only needed for checking model (.check_model), by default None
        n_cores : int, optional
            # of cores to parallelize on, by default os.cpu_count()
        non_machinery : Dict[str, List[str]], optional
            keys are HGNC IDs for non-machinery genes to be expressed, values are a list of compartments within the model for the gene to be transported to, by default None
            Exceptions are ubiquitin genes (HGNC:12468, HGNC:12463) and ribosomal genes
        knock_out : List[str], optional
            each element is the HGNC ID of a gene expressed in the model which should be knocked out

            *Note: you may want to instead knock-out during building if setting minimal_proteome = True and knocking out a
            gene that participates in a OR GPR rule (in case it is the one that is selected); otherwise
            me_model.knock_out() method should suffice, by default None
        additional_ko : List[str], optional
            a list of HGNC IDs for genes that were not explicitly knocked-out, but were only involved in catalysis of
            reactions catalyzed by a complex which contains another gene that was knocked-out, by default None
            this list is generated in build_me_model/me_builder

        Returns
        -------
        Nothing, but initializes the following variables for use later:
        self.m_model: cobra.Model
            the metabolic model from which the ME Model was built (added by builder class)
        self.S: Union[pd.DataFrame, np.array]
            model stoichiometric matrix
        self.solver_
            the LP solver for the model
        self.reaction_types: Dict[str, str]
            keys are category of the reaction, values are the reaction ID

        """
        super().__init__(id_or_model, name)
        if m_model is not None:
            if type(m_model) == cobra.Model:
                self.m_model = m_model.copy()
            else:
                raise ValueError('m_model must be None or a cobra.Model')
        self.S = None
        self.solver_ = None
        self.reaction_types = dict()
        self.n_cores = n_cores
        if self.n_cores in [0, 1, None]:
            self._par = False
        else:
            self._par = True

        if knock_out is None:
            self.knock_out = list()
        else:
            self.knock_out = knock_out

        if additional_ko is None:
            self.additional_ko = list()
        else:
            self.additional_ko = additional_ko

        if non_machinery is None:
            self.non_machinery = dict()
        else:
            self.non_machinery = non_machinery

    def _add_reactions(self, reaction_list: List[cobra.Reaction]):
        """Add reactions to the model
        Reactions with identifiers identical to a reaction already in the
        model are ignored.

        Parameters
        ----------
        reaction_list : List[cobra.Reaction]
            List of reactions to add
        """

        def existing_filter(rxn):
            if rxn.id in self.reactions:
                logger.warning(
                    "Ignoring reaction '%s' since it already exists.", rxn.id)
                return False
            return True

        # First check whether the reactions exist in the model.
        pruned = DictList(filter(existing_filter, reaction_list))

        # Add reactions. Also take care of genes and metabolites in the loop.
        for reaction in pruned:
            reaction._model = self
            # Build a `list()` because the dict will be modified in the loop.
            for metabolite in list(reaction.metabolites):
                # TODO: Should we add a copy of the metabolite instead?
                if metabolite not in self.metabolites:
                    self.add_metabolites(metabolite)
                # A copy of the metabolite exists in the model, the reaction
                # needs to point to the metabolite in the model.
                else:
                    # FIXME: Modifying 'private' attributes is horrible.
                    #                     stoichiometry = reaction._metabolites.pop(metabolite)
                    model_metabolite = self.metabolites.get_by_id(
                        metabolite.id)
                    #                     reaction._metabolites[model_metabolite] = stoichiometry
                    model_metabolite._reaction.add(reaction)

            for gene in list(reaction._genes):
                # If the gene is not in the model, add it
                if not self.genes.has_id(gene.id):
                    self.genes += [gene]
                    gene._model = self

                # Otherwise, make the gene point to the one in the model
                else:
                    model_gene = self.genes.get_by_id(gene.id)
                    if model_gene is not gene:
                        reaction._dissociate_gene(gene)
                        reaction._associate_gene(model_gene)

        self.reactions += pruned

        # make sure index and metabolite/list order are the same
        for idx, r_id in enumerate(self.reactions):
            if idx != self.reactions.index(r_id):
                raise ValueError('Indexing should be changed')

        for idx, m_id in enumerate(self.metabolites):
            if idx != self.metabolites.index(m_id):
                raise ValueError('Indexing should be changed')

    def _assign_reaction_types(self):
        """Organize reactions into their various categories. There will be overlap between the lists"""

        self.reaction_types['biomass'] = [r.id for r in self.reactions if isinstance(r, BiomassReaction)]
        self.reaction_types['metabolism'] = [r.id for r in self.reactions if isinstance(r, MetabolicReaction)]

        # get and initialize expression reactions
        expression_reactions = [r for r in self.reactions if isinstance(r, ExpressionReaction)]
        self.reaction_types['expression'] = {'all': [r.id for r in expression_reactions]}
        self.reaction_types['expression']['synthesis'] = {'protein': [], 'mRNA': [], 'complex': []}
        self.reaction_types['expression']['translation'] = list()
        self.reaction_types['expression']['sink'] = {'protein': [], 'mRNA': [], 'complex': []}
        self.reaction_types['expression']['ribosome_biogenesis'] = list()
        self.reaction_types['expression']['ubiquitin_biogenesis'] = list()

        # categorize expression reactions
        for r in expression_reactions:
            if r.synthesis:
                self.reaction_types['expression']['synthesis'][r.synthesis_type] += [r.id]
            if r.sink:
                self.reaction_types['expression']['sink'][r.sink_type] += [r.id]
            if hasattr(r, 'translation') and r.translation:
                self.reaction_types['expression']['translation'] += [r.id]
            if hasattr(r, 'ribosome_biogenesis') and r.ribosome_biogenesis:
                self.reaction_types['expression']['ribosome_biogenesis'] += [r.id]
            if hasattr(r, 'ubiquitin_biogenesis') and r.ubiquitin_biogenesis:
                self.reaction_types['expression']['ubiquitin_biogenesis'] += [r.id]

        self.reaction_types['coupled'] = [r.id for r in self.reactions if
                                          hasattr(r, 'coupled_metabolites') and r.coupled_metabolites != dict()]

    def _get_enzymes(self, reaction_list: List[cobra.Reaction]):
        """Track all .enzyme attributes. Run before _add_reactions.

        Parameters
        ----------
        reaction_list : List[cobra.Reaction]
            List of reactions to add
        """
        self._enzymes = list()
        for r in reaction_list:
            for m in r.metabolites:
                if hasattr(m, 'enzyme') and m.enzyme:
                    self._enzymes.append(m.id)

    def _map_enzymes(self):
        """When running _add_reactions, for dummy protein, seems to change .enzyme attribute to False.
        Generic correction in case it happens with other enzymes. Related to _get_enzymes() method"""
        for m_id in self._enzymes:
            self.metabolites.get_by_id(m_id).enzyme = True

    def _map_coupled_metabolites(self):
        """Reassigns metabolite object from r.metabolites to the .coupled_metabolites attribute of the reaction
        to ensure that the metabolite object is the most up to date version (prevents multiple copies from existing)"""

        for r in self.reactions:
            if hasattr(r, 'coupled_metabolites'):
                r._map_coupled_metabolites()

    def _map_metabolite_reactions_and_coupling(self):
        """Fixes error in which metabolites do not have all associated reactions or coupling constraints
        in the .reactions and .coupling_coefficients attributes respectively"""

        metab_reaction_map = {m.id: list() for m in self.metabolites}
        metab_coupling_map = {m.id: dict() for m in self.metabolites if hasattr(m, 'coupling_coefficient')}
        for r in self.reactions:
            for m in r.metabolites:
                metab_reaction_map[m.id] += [r.id]

                if hasattr(m, 'coupling_coefficient') and m.coupling_coefficient is not None:
                    for k, c in m.coupling_coefficient.items():
                        if k in metab_coupling_map[m.id]:
                            if c != metab_coupling_map[m.id][k]:
                                raise ValueError(
                                    'Disagreement in coupling coefficient for the same metabolite: ' + m.id)
                        else:
                            metab_coupling_map[m.id][k] = c
        for m_id, val in metab_coupling_map.items():
            if len(val) == 0:
                metab_coupling_map[m_id] = None

        for m_id, r_list in metab_reaction_map.items():
            metab = self.metabolites.get_by_id(m_id)
            metab._reaction = metab._reaction.union([self.reactions.get_by_id(r_id) for r_id in r_list])
            if m_id in metab_coupling_map:
                metab.coupling_coefficient = metab_coupling_map[m_id]

    def _map_reaction_metabolites(self):
        """Fixes error in which reactions do not have the most up to date metabolites"""

        for r in self.reactions:
            for m in r.metabolites:
                m_ = self.metabolites.get_by_id(m.id)
                if m != m_:
                    r._metabolites[m_] = r._metabolites.pop(m)

    def _clean_metabolites(self):
        """Remove or correct reactions assigned to metabolites which are not in the model"""
        rxn_ids = [r.id for r in self.reactions]
        for m in self.metabolites:
            for r in m.reactions:
                if r.id not in rxn_ids:
                    m._reaction.remove(r)
                elif r not in self.reactions:
                    m._reaction.remove(r)
                    m._reaction.add(self.reactions.get_by_id(r.id))

    def correct_object_tracking(self):
        """Resolves inconsistencies b/w metabolite.reactions and reaction.metabolites or reaction.coupled_metabolites"""
        # update metabolites
        self._clean_metabolites()
        self._map_enzymes()
        self._map_metabolite_reactions_and_coupling()
        # update reactions
        self._map_reaction_metabolites()
        self._map_coupled_metabolites()

    def add_reactions(self, reaction_list: List[cobra.Reaction]):
        """Add reactions to the model

        Parameters
        ----------
        reaction_list : List[cobra.Reaction]
            List of reactions to add
        """
        self._get_enzymes(reaction_list)
        self._add_reactions(reaction_list)
        self.correct_object_tracking()
        self._assign_reaction_types()

    def remove_reactions(self, reactions: List[Union[cobra.Reaction, str]], remove_orphans: bool = True):
        """Remove reactions from the model.

        Parameters
        ----------
        reactions : List[cobra.Reaction, str]
            A list with reactions (`cobra.Reaction`), or their id's, to remove
        remove_orphans : bool
            Remove orphaned genes and metabolites from the model as well, , by default True
        """
        if isinstance(reactions, string_types) or hasattr(reactions, "id"):
            warnings.warn("need to pass in a list")
            reactions = [reactions]

        for reaction in reactions:
            # Make sure the reaction is in the model
            try:
                reaction = self.reactions[self.reactions.index(reaction)]
            except ValueError:
                warnings.warn(reaction.id + 'not in model')

            else:
                self.reactions.remove(reaction)
                reaction._model = None

                for met in reaction._metabolites:
                    if reaction in met._reaction:
                        met._reaction.remove(reaction)
                        if remove_orphans and len(met._reaction) == 0:
                            self.remove_metabolites(met)

                for gene in reaction._genes:
                    if reaction in gene._reaction:
                        gene._reaction.remove(reaction)
                        if remove_orphans and len(gene._reaction) == 0:
                            self.genes.remove(gene)

                # remove reference to the reaction in all groups
                associated_groups = self.get_associated_groups(reaction)
                for group in associated_groups:
                    group.remove_members(reaction)

    # new methods------------------------------------------------------------------------------------------------

    def create_stoichiometric_matrix(self, array_type: str = 'numpy', mu_val: Optional[SupportsFloat] = None, inplace: bool = True):
        """Generate a stoichiometric array representation of the model.
        Adapted from cobra.util.array.create_stoichiometric_matrix to handle sympy.Expr objects.

        Parameters
        ----------
        array_type : str, optional
            Specifies the type of the stoichiometric matrix to be return (options ['numpy', 'pandas', 'sympy']), by default 'numpy'
        mu_val : Optional[SupportsFloat], optional
             A value for mu to replace in generating stoichiometric matrix. If None, will use sympy expressions, by default None
        inplace : bool, optional
            whether to update self.S (True) or return a new array (False), by default True

        Returns
        -------
        Union[pd.DataFrame, np.array, sympy.matrices.dense.MutableDenseMatrix]
            The the columns represent the reactions and rows represent
            metabolites. S[i,j] therefore contains the quantity of metabolite `i`
            produced (negative for consumed) by reaction `j`.
        """

        if array_type not in ['sympy', 'numpy', 'pandas']:
            raise TypeError('Incorrect array type specified')
        if array_type != 'sympy' and mu_val is None:
            raise ValueError('Must specify a mu_val for non-sympy matrices')
        if array_type == 'sympy' and mu_val is not None:
            warnings.warn(
                'Sympy array type will generate expression entries, mu_val will be disregarded. Use .replace_S_mu() to generate a numpy matrix with a specific mu_val')

        n_metabolites = len(self.metabolites)
        n_reactions = len(self.reactions)
        # initialize empty matrix
        array = np.zeros((n_metabolites, n_reactions))
        if array_type == 'sympy':
            array = sympy.Matrix(array)

        m_ind = self.metabolites.index
        r_ind = self.reactions.index

        if array_type == 'sympy':
            for reaction in self.reactions:
                for metabolite, stoich in iteritems(reaction.metabolites):
                    array[m_ind(metabolite), r_ind(reaction)] = stoich
        else:
            for reaction in self.reactions:
                #             reaction_type = isinstance(reaction, core.reaction.ME_Reaction)
                for metabolite, stoich in iteritems(reaction.metabolites):
                    if isinstance(stoich, sympy.Expr):
                        array[m_ind(metabolite.id), r_ind(reaction)] = float(stoich.subs(params.mu, mu_val))
                    else:
                        array[m_ind(metabolite.id), r_ind(reaction)] = stoich

        if array_type == 'pandas':
            metabolite_ids = [met.id for met in self.metabolites]
            reaction_ids = [rxn.id for rxn in self.reactions]
            array = pd.DataFrame(array, index=metabolite_ids, columns=reaction_ids)
        else:
            if array_type == 'sympy':
                self.replace_S_mu = lambdify(params.mu, array, modules='numpy')

        if inplace:
            self.S = array
        else:
            return array

    def initialize_solver(self, solver_type: str = 'qminos', precision: str = 'quad'):
        """Initialize the ME Model solver

        solver_type : str
            The solver to use for the linear programs (options ['qminos']), by default "qminos"
        precision: str
            The precision for the qminos solver (options ['double', 'quad', 'dq', 'dqq']), by default 'quad'
        """

        if solver_type == 'qminos':
            self.solver_ = solve_me.qminosSolver(precision=precision)
            self.solver_type = 'qminos'
            self.solver_precision = 'quad'
        else:
            raise ValueError('Only the qMINOS solver is currently implemented')

    def solve_lp(self, mu_val: SupportsFloat, objective: Optional[Dict[str, int]] = None, tolerance: SupportsFloat = 0) -> Tuple[np.array, int]:
        """Solves the linear program for a specified objective at a specified growth rate.

         Parameters
         ----------
         mu_val : SupportsFloat
             The growth value for which to solve the linear program [hr^-1]
         objective : Dict[str, int], optional
             The objective function to optimize. Dictionary represent a linear combination of reactions to optimize,
             with reaction ids as keys and the coefficient of the linear combination as the values.
             Values must either be 1 for maximization or -1 for minimization, by default {'biomass_dilution': 1}
         tolerance : float, optional
             Threshold below which expected sensitivity of solver is too low to detect infeasibility, by default 0

         Returns
         -------
         sln: np.array
             1D vector of fluxes in the optimal solution
         stat: int
             the solver status
                 0     Optimal solution found.
                 1     The problem is infeasible.
                 2     The problem is unbounded (or badly scaled).
                 3     Too many iterations.
                 4     Apparent stall.  The solution has not changed
                       for a large number of iterations (e.g. 1000).
         hsq:
             optimal basis (see qminospy.solver.QMINOS)
         """
        if objective is None:
            objective = {'biomass_dilution': 1}
        if self.solver_ is None:
            warnings.warn(
                'Solver is not initialized with ME_Model.intialize_solver, intializing with default parameters')
            self.initialize_solver()
        else:
            self.initialize_solver(solver_type=self.solver_type, precision=self.solver_precision)

        sln, stat, hs = self.solver_.solve_lp(me_model=self, mu_val=mu_val, objective=objective, tolerance=tolerance)
        return sln, stat, hs

    def maximize_growth(self, min_mu: SupportsFloat = 0, max_mu: SupportsFloat = 0.05,
                        mu_accuracy: SupportsFloat = 1e-10, increment: SupportsFloat = 0.02,
                        tolerance: SupportsFloat = 0, verbose: bool = True):
        """Binary search to find the maximum feasible growth rate.

        Parameters
        ----------
        me_model : human_me.core.model.ME_Model
            ME Model to solve
        min_mu : SupportsFloat, optional
            Expected minimum feasible growth rate, by default 0
        max_mu : SupportsFloat, optional
            Expected minimum infeasible growth rate (i.e., just above expected maximum feasible growth rate), by default 0.05
        mu_accuracy : SupportsFloat, optional
            The maximum error in mu after the binary search, by default 1e-10
        increment : SupportsFloat, optional
            The amount to increase growth by when searching for maximum infeasible growth rate from max_mu, by default 0.02
        tolerance : SupportsFloat, optional
            Threshold below which expected sensitivity of solver is too low to detect infeasibility, by default 0
        verbose : bool, optional
            Prints information about each linear program iteration, by default True

        Returns
        -------
        mu_max: int
            the maximum feasible growth value (in hours)
        res: Dict[float, Tuple[np.array, int]]]
            keys are all attempted growth values, values are dictionaries with keys as output from self.solve_lp
        """
        if self.solver_ is None:
            warnings.warn(
                'Solver is not initialized with ME_Model.intialize_solver, intializing with default parameters')
            self.initialize_solver()
        else:
            self.initialize_solver(solver_type=self.solver_type, precision=self.solver_precision)

        mu_max, res = self.solver_.maximize_growth(me_model=self,
                                                   min_mu=min_mu, max_mu=max_mu,
                                                   mu_accuracy=mu_accuracy, increment=increment,
                                                   tolerance=tolerance,
                                                   verbose=verbose)
        return mu_max, res

    def optimize(self, objective: Dict[str, int], mu_max: SupportsFloat, n_points: int = 10,
                 tolerance: SupportsFloat = 0, n_cores: Optional[int] = None, visualize: bool = True, fig_name: str = None):
        """General optimization of any non-growth objective

        Parameters
        ----------
        objective : Dict[str, int]
            The objective function to optimize. Dictionary represent a linear combination of reactions to optimize,
            with reaction ids as keys and the coefficient of the linear combination as the values.
            Values must either be 1 for maximization or -1 for minimization.
        mu_max : SupportsFloat
            the maximum growth value at which the model is feasible [hr^-1]; use .maximize_growth() method to identify (should be <= mu_max output)
            if using an experimental value, make sure it is feasible using the .solve_lp() method 
        n_points : int, optional
            # of growth values to consider between 0 and mu_max, by default 10
        tolerance : SupportsFloat, optional
            Threshold below which expected sensitivity of solver is too low to detect infeasibility, by default 0
        n_cores : Optional[int], optional
            the number of workers to use for parallelization, by default None
        visualize : bool, optional
            plot the relationship between growth and the objective function of interest, by default True
        fig_name : Optional[str], optional
            save the plotted figure to 'path/to/filename.ext', by default None

        Returns
        ----------
        sln: Tuple[float]
            first element is the growth value at which the non-growth objective is optimized
            second element is the optimized non-growth objective value
        predicted: pandas.DataFrame
            1000 growth values between 0 and mu_max, with corresponding interpolated objective values
        interp_fit: scipy.interpolate.interp1d
            a function to interpolate objective values from growth values, used to generated predicted
        optimal_vals: collections.OrderedDict
            keys are n_points growth values between 0 and mu_max, values are the objective value optimized at
            the corresponding growth value
        res: Dict[float, Tuple[Any]]
            keys are n_points growth values between 0 and mu_max, values are the output of .solve_lp at
            corresponding growth values with the objective set to the non-growth objective input
        """

        if self.solver_ is None:
            warnings.warn(
                'Solver is not initialized with ME_Model.intialize_solver, intializing with default parameters')
            self.initialize_solver()
        else:
            self.initialize_solver(solver_type=self.solver_type, precision=self.solver_precision)

        sln, predicted, interp_fit, optimal_vals, res = self.solver_.optimize(me_model=self, objective=objective,
                                                                              mu_max=mu_max, n_points=n_points,
                                                                              tolerance=tolerance,
                                                                              n_cores=n_cores, visualize=visualize,
                                                                              fig_name=fig_name)
        return sln, predicted, interp_fit, optimal_vals, res

    def infeasible_reactions(self, mu_val: SupportsFloat, sln, stat, tolerance: SupportsFloat = 1e-19) -> Dict[str, SupportsFloat]:
        """Binary search to find the maximum feasible growth rate

        Parameters
        ----------
        mu_val : SupportsFloat
            input growth value to ME_Model.solve_lp
        sln, stat : outputs of ME_Model.solve_lp
            Expected minimum feasible growth rate (~0)
        tolerance : SupportsFloat
            Threshold below which expected sensitivity of solver is too low to detect infeasibility

        Returns
        ----------
        ir : Dict[str, SupportsFloat]
            for reactions that cause infeasibility, keys are reaction ids for infeasible reactions and values are
            difference by which reaction flux is infeasible
        """
        if tolerance < 0:
            tolerance = abs(tolerance)

        ir = dict()
        for r in self.reactions:
            flux = sln[self.reactions.index(r.id)]

            ub = copy.copy(r.upper_bound)
            lb = copy.copy(r.lower_bound)

            if isinstance(ub, sympy.Expr):
                ub = float(ub.subs(params.mu, mu_val))
            if isinstance(lb, sympy.Expr):
                lb = float(lb.subs(params.mu, mu_val))

            if math.isnan(flux):
                ir[r.id] = flux
            elif flux > ub + tolerance:
                ir[r.id] = abs(flux - ub)
            elif flux < lb - tolerance:
                ir[r.id] = abs(lb - flux)

        if (len(ir) > 0 and stat == 0) or (len(ir) == 0 and stat != 0):
            warnings.warn(
                'There is a discrepancy between the solver status and reactions that violate bound constraints')
        return ir

    def check_me_mass_balance(self):
        """Checks that all reactions in ME Model are mass balance. Use after self.add_reactions"""
        print('Check reaction mass balances')

        for r in self.reactions:
            if len(r.check_mass_balance()) > 0:

                # account for instances when the original metabolic reaction is unbalances
                bool1 = isinstance(r, MetabolicReaction)
                if bool1:
                    # account for instances when the ME_Reaction version of reversible reactions
                    if self.m_model.reactions.get_by_id(r.cobra_id).reversibility and r.id.endswith(
                            '_R'):  # account for reversible reactions
                        bool1 = self.m_model.reactions.get_by_id(r.cobra_id).check_mass_balance() == {e: -c for e, c in
                                                                                                      r.check_mass_balance().items()}
                    else:
                        bool1 = r.check_mass_balance() == self.m_model.reactions.get_by_id(
                            r.cobra_id).check_mass_balance()
                else:
                    bool1 = False

                bool2 = isinstance(r, BiomassReaction)
                if not (bool1 or bool2):
                    raise ValueError('Atleast one reaction is not mass balanced')

    def _check_complete_reactions(self):
        """Checks that all the original metabolic model reactions have been included in the ME-Model"""
        if len(set([r.id for r in self.m_model.reactions]).difference(
                [self.reactions.get_by_id(r_id).cobra_id for r_id in self.reaction_types['metabolism']])) > 0:
            raise ValueError('Not all the original metabolic model reactions have been included in the ME-Model')

    def _check_hgncs(self):
        """Checks that all reactions and macromolecules that are expected to have an assigned hgnc_id, do"""

        # reactions
        bool1 = len([r for r in self.reactions if not hasattr(r, 'hgnc_id') and not (isinstance(r, (BiomassReaction, MetabolicReaction)))]) > 0

        no_hgnc = [r for r in self.reactions if type(r) == ExpressionReaction and r.hgnc_id is None]
        no_hgnc = [r for r in no_hgnc if (r.subsystem not in ['tRNA_Biogenesis', 'rRNA_expression', 'Complex_Formation',
                                                              'Complex_Degradation']) and (
            hasattr(r, 'ubiquitin_biogenesis') and not r.ubiquitin_biogenesis)]
        bool2 = len(no_hgnc) > 0

        if bool1 or bool2:
            raise ValueError('An expression reaction does not have an hgnc_id')

        # macromolecules
        fragments = ['3_trailer', '5_leader', 'ets', 'its']
        exceptions = ['ubiquitin_monomer_protein_c', 'cleaved_polyubiquitin_moiety_protein_c',
                      'ubiquitin_monomer_protein_n', 'cleaved_polyubiquitin_moiety_protein_n']
        hgnc_id_metabs = [m for m in self.metabolites if
                          isinstance(m, Macromolecule) and m.hgnc_id is None and m.type not in ['trna', 'rrna',
                                                                                                'complex'] and not [
                              hasattr(m, 'amt') and m._amt == 'complex'] and not (hasattr(m,
                                                                                          'fragment_type') and m.fragment_type in fragments) and m.id not in exceptions]
        if len(hgnc_id_metabs) > 0:
            raise ValueError('Some macromolecules did not get an HGNC ID assigned')

    def _check_coupling(self):
        """Checks that all reactions have received appropriate machinery (compares coupled metabolites to GPR)"""
        print('Make sure all reactions received correct coupled machinery')

        # set arguments
        if 'orphan' not in self.reaction_types:
            raise ValueError('Must specify a list of orphan reaction IDs')
        orphan = [self.reactions.get_by_id(r_id) for r_id in self.reaction_types['orphan']]

        test_reactions = [r for r in self.reactions if not (r.id in orphan) and not isinstance(r, BiomassReaction)]
        for r in tqdm(test_reactions):
            if isinstance(r, MetabolicReaction):
                r_ = self.m_model.reactions.get_by_id(r.cobra_id).copy()
            else:  # ExpressionReaction
                r_ = r.copy()

            if 'or' in r_.gene_reaction_rule and 'and' in r_.gene_reaction_rule:
                expected_machinery = parse_complex.eval_complex(r_.gene_reaction_rule)
            elif len(r_.genes) == 0:
                expected_machinery = []
            else:
                expected_machinery = [g.id for g in r_.genes]
                if 'and' in r_.gene_reaction_rule:
                    expected_machinery = [expected_machinery]

            ko = False

            expected_machinery2 = list()
            for em in expected_machinery:
                if type(em) == list:
                    expected_machinery2 += em
                else:
                    expected_machinery2.append(em)

            if len(set(expected_machinery2).difference(self.knock_out + self.additional_ko)) == 0:
                ko = True

            if (not hasattr(r, 'coupled_metabolites') or 'catalysis' not in r.coupled_metabolites.values()) and not ko:
                raise ValueError('Reaction does not have a record of coupled machinery: ' + r.id)
            actual_machinery = [m for m, v in r.coupled_metabolites.items() if v == 'catalysis']

            if len(actual_machinery) > 0:
                translation = isinstance(actual_machinery[0], RibosomalComplex)
                ribosomal_degradation = (len(actual_machinery) == 2)
            else:
                translation = False
                ribosomal_degradation = False
            if not ko:
                if not (translation or ribosomal_degradation):
                    # get complex or protein machinery information
                    cplx = False
                    if isinstance(actual_machinery[0], Complex):
                        cplx = True
                        am = list()
                        for p in actual_machinery[0].decompose_complex():
                            if p.type != 'protein':
                                raise ValueError('Non-proteins in complex machinery for ' + r.id)
                            am.append(p.id.split('_')[0])
                        actual_machinery = sorted(am)

                    # check that machinery matches
                    if len(expected_machinery) > 0:  # non-dummy
                        err = True
                        for rm in expected_machinery:
                            if type(rm) != list:
                                if not cplx and actual_machinery[0].id.split('_')[0] == rm:
                                    err = False
                            else:
                                if cplx and sorted(rm) == actual_machinery:
                                    err = False
                        if err:
                            if expected_machinery != self.knock_out:
                                raise ValueError('Machinery mismatch for ' + r.id)
                    else:
                        if len(actual_machinery) > 1 or not actual_machinery[0].dummy:  # dummy
                            raise ValueError('Non-dummy protein coupled to deorphaned reaction')
                elif ribosomal_degradation:
                    am = list()
                    for am_ in actual_machinery:
                        for p in am_.decompose_complex():
                            if p.type != 'protein':
                                raise ValueError('Non-proteins in complex machinery for ' + r.id)
                            am.append(p.id.split('_')[0])
                    actual_machinery = sorted(am)
                    if sorted(expected_machinery[0]) != actual_machinery:
                        raise ValueError('Incorrect machinery for ribosomal degradation: ' + r.id)
                elif translation:
                    actual_machinery = sorted(
                        [p.id.split('_')[0] for p in actual_machinery[0].decompose_complex() if p.type == 'protein'])
                    expected_machinery = sorted([p for p in expected_machinery[0] if p != 'ribosome'] + rbps)
                    if actual_machinery != expected_machinery:
                        raise ValueError('Incorrect machinery for translation: ' + r.id)
                else:
                    raise ValueError('Unaccounted for reaction criteria')

    def check_enzymes(self):
        """Makes sure all genes being expressed participate in a catalysis reaction (no unecessary expression reactions)"""

        proteins, complexes = [], []
        active_proteins, active_complexes = [], []
        for m in self.metabolites:
            if hasattr(m, 'type'):
                if m.type == 'protein':
                    proteins.append(m.id)
                    if m.enzyme:
                        active_proteins.append(m.id)
                elif m.type == 'complex':
                    complexes.append(m)
                    if m.enzyme:
                        active_complexes.append(m)
        # check complexes
        complexes = [m for m in complexes if not (
            isinstance(m, RibosomalComplex) or '_polyub_complex_' in m.id or np.all(
                    [isinstance(r_, ComplexDegradationReaction) for r_ in m.reactions]))]
        if len(list(set(complexes).difference(active_complexes))) > 0:
            raise ValueError('Unexpected inclusion of inactive complexes')

        # check monomers
        active_proteins += flatten_list([[p.id for p in m.decompose_complex()] for m in active_complexes])
        active_proteins = list(set([i.split('_')[0] if i.startswith('HGNC') else i for i in active_proteins]))
        proteins = set([i.split('_')[0] if 'HGNC' in i else i for i in proteins])
        proteins = [i for i in proteins if 'ubiquitin' not in i]
        ub_genes = ['HGNC:12463', 'HGNC:12468']
        if len(set(proteins).difference(
                active_proteins + self.additional_ko + ub_genes + list(self.non_machinery.keys()))) > 0:
            raise ValueError('Unexpected inclusion of inactive protein monomers')

    def check(self):
        """Check ME Model built correctly

        """
        self._check_complete_reactions()
        self._check_hgncs()
        self.check_me_mass_balance()
        self._check_coupling()
        self.check_enzymes()

    @staticmethod
    def create_expressed_gene(hgnc_id: str, relat_objects) -> ExpressedGene:
        g = ExpressedGene(hgnc_id)
        for m in relat_objects['macromolecules']:
            g.add_macromolecule(m)
        for r in relat_objects['reactions']:
            g.add_reaction(r)
        g.check()

        return g

    def _generate_expressed_genes(self):

        # map all hgnc_ids to all associated macromolecules and reactions
        hgnc_ids = set([m.hgnc_id for m in self.metabolites if isinstance(m, Macromolecule) and m.hgnc_id is not None])
        hgnc_ids = {hgnc_id: {'macromolecules': []} for hgnc_id in hgnc_ids}

        for m in self.metabolites:
            if isinstance(m, Macromolecule) and m.type in ['premrna', 'fragment_rna', 'mrna', 'protein', 'complex',
                                                           'proxy']:
                if m.hgnc_id is not None:
                    hgnc_ids[m.hgnc_id]['macromolecules'] += [m.id]
                else:
                    if m.type == 'complex':
                        c_hids = [p.hgnc_id for p in m.decompose_complex() if p.hgnc_id is not None]
                        for c_hid in c_hids:
                            hgnc_ids[c_hid]['macromolecules'] += [m.id]
                    elif m.type == 'proxy':
                        for c_hid in m._complex_hgnc_ids:
                            hgnc_ids[c_hid]['macromolecules'] += [m.id]

        #         for hgnc_id, v in hgnc_ids.items():
        #             v['reactions'] = list(set(flatten_list([[r.id for r in self.metabolites.get_by_id(m_id).reactions] for m_id in v['macromolecules']])))

        for v in hgnc_ids.values():
            v['reactions'] = list(set(flatten_list(
                [[r.id for r in self.metabolites.get_by_id(m_id).reactions] for m_id in v['macromolecules']])))
            v['macromolecules'] = [self.metabolites.get_by_id(m_id) for m_id in v['macromolecules']]
            v['reactions'] = [self.reactions.get_by_id(r_id) for r_id in v['reactions']]

        print('Add gene objects')
        if self._par:
            pool = multiprocessing.Pool(processes=self.n_cores)
            try:
                expressed_genes = pool.starmap(self.create_expressed_gene,
                                               zip(hgnc_ids.keys(), hgnc_ids.values()))
                pool.close()
                pool.join()
                gc.collect()
                self.expressed_genes = {g.hgnc_id: g for g in expressed_genes}
            except:
                pool.close()
                pool.join()
                gc.collect()
                raise ValueError('Parallelization failed while generating expressed genes list')
        else:
            self.expressed_genes = {hgnc_id: self.create_expressed_gene(hgnc_id, relat_objects) for hgnc_id, relat_objects in
                                    tqdm(hgnc_ids.items())}

    def knock_out_gene(self, hgnc_id: str, inplace: bool = False):
        """Knocks out a gene by blocking flux through synthesis of the associated mRNA molecule

        Parameters
        ----------
        hgnc_id : str
            gene id in HGNC format (HGNC:####); must be present in the model
        inplace : bool
            whether to return a separate model object with the knocked out reaction (False) or to
            change the reaction bounds in change, by default False
        """

        if hgnc_id not in self.expressed_genes:
            raise ValueError('The specified hgnc_id is not present in the ME Model')
        self.knock_out.append(hgnc_id)
        r_id = self.expressed_genes[hgnc_id].reactions['ExpressionReactions']['mrna']['synthesis']

        if inplace:
            self.reactions.get_by_id(r_id)._lower_bound = 0
            self.reactions.get_by_id(r_id)._upper_bound = 0
        else:
            me_model_copy = copy.deepcopy(self)
            me_model_copy.reactions.get_by_id(r_id)._lower_bound = 0
            me_model_copy.reactions.get_by_id(r_id)._upper_bound = 0
            return me_model_copy

    def pickle(self, file=os.path.join(os.path.abspath(os.getcwd()), 'me_model.pickle')):
        """Save ME_Model as a pickled object

        Parameters
        ----------
        file: str 
            will save to file = "full/path/to/filename.pickle", default './me_model.pickle'

        """
        if '.' in file:
            p = pathlib.Path(file)
            extensions = "".join(p.suffixes)
            file = str(p).replace(extensions, '.pickle')
        else:
            file = file + '.pickle'
        with open(file, 'wb') as handle:
            pickle.dump(self, handle)

    @staticmethod
    def load_pickled_model(file_name: str):
        """Loads a pickled me_model. Saved from me_model.pickle

        Parameters
        ----------
        file_name : str
            'full/path/to/me_model.pickle'

        Returns
        -------
        ME_Model
            ME model object
        """

        with open(file_name, 'rb') as handle:
            me_model = pickle.load(handle)
        me_model.correct_object_tracking()  # lost in pickling/loadings
        return me_model

    def copy(self):
        return copy.deepcopy(self)
