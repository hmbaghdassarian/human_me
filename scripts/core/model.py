#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle
import os
import pathlib
import copy

import cobra
from cobra.core.dictlist import DictList
from six import iteritems, string_types

import sympy
from sympy import lambdify
import numpy as np
import pandas as pd
import warnings
import copy
import math
from tqdm import tqdm

import sys
sys.path.insert(1, '../../scripts/')
from preprocess import parse_complex
from utils import parameters as params
from macromolecules.macromolecule import Macromolecule
from macromolecules.complex import Complex
from macromolecules.protein import Protein
from macromolecules.complex import Ribosomal_Complex

import core
from core.reaction import Biomass_Reaction, Expression_Reaction, Metabolic_Reaction, Complex_Degradation_Reaction

from me_solver import solve_me
# obj_type = type


# In[ ]:


def flatten_list(t):
    #https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    return [item for sublist in t for item in sublist]


# In[2]:


class ME_Model(cobra.Model):
# rewritten methods------------------------------------------------------------------------------------------------
    def __init__(self,  id_or_model, name = None, m_model = None):
        '''
        A simple object with an identifier
    
        Parameters
        ----------
        m_model: cobra.Model
            The cobrapy model object that the ME_Model was built from. Only needed for checking model (.check_model) 
        id_or_model: None or a string
            the identifier to associate with the object
        
        Returns
        -------
        Nothing, but initializes the following variables for use later:
        self.m_model: cobra.Model
            the metabolic model from which the ME Model was built (added by builder class)
        self.S:
            model stoichiometric matrix
        self.solver_:
            the LP solver for the model
        self.orphan: list
            a list of reactions that remain orphaned when dummy is being incorporated (added by builder class)
        self.deorphaned: list
            a list of reactions that are deoprhaned whend dummy protein is being incorporated (added by builder class)
            
        '''
        
        super().__init__(id_or_model, name)
        if m_model is not None:
            if type(m_model) == cobra.Model:
                self.m_model = m_model.copy()
            else:
                raise ValueError('m_model must be None or a cobra.Model')
        self.S = None
        self.solver_ = None
        self.reaction_types = dict()
        

    def _add_reactions(self, reaction_list):
        """Add reactions to the model.

        Reactions with identifiers identical to a reaction already in the
        model are ignored.

        Parameters
        ----------
        reaction_list : list
            A list of `cobra.Reaction` objects
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
        
        self._clean_metabolites()
        
    def _map_coupled_metabolites(self, verbose = False):
        '''Reassigns metabolite object from r.metabolites to the .coupled_metabolites attribute of the reaction
        to ensure that the metabolite object is the most up to date version (prevents multiple copies from existing)'''
        
        if verbose:
            print('Reassign .coupled_metabolites attribute')
        for r in self.reactions:
            if hasattr(r, 'coupled_metabolites'):
                r._map_coupled_metabolites()
    
    def _map_metabolite_reactions(self):
        '''Fixes strange error in which metabolites do not have all associated reactions 
        in the .reactions attribute'''
    
        metab_reaction_map = {m.id: list() for m in self.metabolites}
        for r in me_model.reactions:
            for m in r.metabolites:
                metab_reaction_map[m.id] += [r.id]

        for m_id, r_list in metab_reaction_map.items():
            metab = self.metabolites.get_by_id(m_id)
            metab._reaction = metab._reaction.union([self.reactions.get_by_id(r_id) for r_id in r_list])    
    
    def _assign_reaction_types(self):
        '''Organize reactions into their various categories. There will be overlap between the lists'''
        
        self.reaction_types['biomass'] = [r.id for r in self.reactions if isinstance(r, Biomass_Reaction)]
        self.reaction_types['metabolism'] = [r.id for r in self.reactions if isinstance(r, Metabolic_Reaction)]
        

        # get and initialize expression reactions
        expression_reactions = [r for r in self.reactions if isinstance(r, Expression_Reaction)]
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
                
        self.reaction_types['coupled'] = [r.id for r in self.reactions if hasattr(r, 'coupled_metabolites') and r.coupled_metabolites != dict()]
                
        
    def add_reactions(self, reaction_list, verbose = False):
        self._add_reactions(reaction_list)
        self._map_metabolite_reactions()
        self._map_coupled_metabolites(verbose = verbose)
        self._assign_reaction_types()
        
    
    def remove_reactions(self, reactions, remove_orphans=True):
        """Remove reactions from the model.

        Parameters
        ----------
        reactions : list
            A list with reactions (`cobra.Reaction`), or their id's, to remove

        remove_orphans : bool, default True
            Remove orphaned genes and metabolites from the model as well

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
    def _clean_metabolites(self):
        '''Remove reactions assigned to metabolites which are not in the model'''
        for m in self.metabolites:
            for r in m.reactions:
                if r not in self.reactions:
                    m._reaction.remove(r)
                    
    def create_stoichiometric_matrix(self, array_type = 'numpy', mu_val = None, inplace = True):

        """
        Adapted from cobra.util.array.create_stoichiometric_matrix to take in sympy.Expr objects

        Return a stoichiometric array representation of the given model.

        The the columns represent the reactions and rows represent
        metabolites. S[i,j] therefore contains the quantity of metabolite `i`
        produced (negative for consumed) by reaction `j`.

        Parameters
        -------
        array_type: one of ['numpy', 'pandas', 'sympy']
            Specifies the type of the stoichiometric matrix to be return

        Returns
        -------
        matrix of class `sympy.matrices.dense.MutableDenseMatrix` by default
            The stoichiometric matrix for the given model.
        """

        if array_type not in ['sympy', 'numpy', 'pandas']:
                raise ValueError('Incorrect array type specified')
        if array_type != 'sympy' and mu_val is None:
            raise ValueError('Must specify a mu_val for non-sympy matrices')
        if array_type == 'sympy' and mu_val is not None:
            warnings.warn('Sympy array type will generate expression entries, mu_val will be disregarded. Use .replace_S_mu() to generate a numpy matrix with a specific mu_val')

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

    def initialize_solver(self, solver_type = 'qminos', precision = 'quad'):
        '''Initialize the ME Model solver
        
        solver_type: string, default "qminos"
            The solver to use for the linear programs (no other options currently )
        precision: string, default "quad"
            The precision for the qminos solver (options ['double', 'quad', 'dq', 'dqq'])
        
        '''
        
        if solver_type == 'qminos':
            self.solver_ = solve_me.qminos_solver(precision = precision)
            self.solver_type = 'qminos'
            self.solver_precision = 'quad'
        else:
            raise ValueError('Only the qMINOS solver is currently implemented')
    
    def solve_lp(self, mu_val, objective = {'biomass_dilution': 1}, tolerance = 0):
        '''Solves the linear program for a specified objective at a specified growth rate

        Parameters
        ----------
        mu_val: float
            The growth value for which to solve the linear program
        objective: dict, default {'bimoass_dilution': 1}
            The objective function to optimize. Dictionary represent a linear combination of reactions to optimize,
            with reaction ids as keys and the coefficient of the lin. comb. as the values. 
            Values must either be 1 for maximization or -1 for minimization.

            Example: simplest case, to maximize reaction with id 'A', objective = {'A': 1}
        tolerance: float; default 0
            Threshold below which expected sensitivity of solver is too low to detect infeasibility

        Returns
        ----------
        sln: 1D np.array
            the vector of fluxes in the optimal solution
        stat: 
            the solver status 
                0     Optimal solution found.
                1     The problem is infeasible.
                2     The problem is unbounded (or badly scaled).
                3     Too many iterations.
                4     Apparent stall.  The solution has not changed
                      for a large number of iterations (e.g. 1000).
        hsq: 
            optimal basis (see qminospy.solver.QMINOS)
        '''
        if self.solver_ is None:
            warnings.warn('Solver is not initialized with ME_Model.intialize_solver, intializing with default parameters')
            self.initialize_solver()
        else:
            self.initialize_solver(solver_type = self.solver_type, precision = self.solver_precision)
            
        sln, stat, hs = self.solver_.solve_lp(me_model = self, mu_val = mu_val, objective = objective, tolerance = tolerance)
        return sln, stat, hs
    
    def maximize_growth(self, min_mu=0, max_mu=0.05, mu_accuracy=1e-10, increment = 0.02, 
                        tolerance = 0, verbose=True):
    
        '''Binary search to find the maximum feasible growth rate

        Parameters
        ----------
        min_mu: float, default 0
            Expected minimum feasible growth rate (~0)
        max_mu: float, default 0.05
            Expected minimum infeasible growth rate (i.e., just above expected maximum feasible growth rate)
        mu_accuracy: float, default 1e-4
            The maximum error in mu after the binary search
        increment: float, default 1
            The amount to increase growth by when searching for maximum infeasible growth rate from max_mu
        verbose: bool, default True
            Prints information about each linear program iteration

        Returns
        ----------
        mu_max: int
            the maximum feasible growth value (in hours)
        res: dict
            keys are all attempted growth values, values are dictionaries with keys as output from self.solve_lp
        '''
        if self.solver_ is None:
            warnings.warn('Solver is not initialized with ME_Model.intialize_solver, intializing with default parameters')
            self.initialize_solver()
        else:
            self.initialize_solver(solver_type = self.solver_type, precision = self.solver_precision)
        
        mu_max, res = self.solver_.maximize_growth(me_model = self, 
                                                     min_mu=min_mu, max_mu=max_mu, 
                                                     mu_accuracy=mu_accuracy, increment = increment,
                                                  tolerance = tolerance,
                                                     verbose=verbose)
        return mu_max, res

    def optimize(self, objective, mu_max, n_points = 10, 
                 tolerance = 0, n_cores = None, graph = True, fig_name = None):
        '''General optimization of any non-growth objective
        
        Parameters
        ----------
        objective: dict
            The objective function to optimize. Dictionary represent a linear combination of reactions to optimize,
            with reaction ids as keys and the coefficient of the lin. comb. as the values. 
            Values can only be all 1 for maximization, or all -1 for minimization. 

            Example: simplest case, to maximize reaction with id 'A', objective = {'A': 1}
        mu_max: float
            the maximum growth value at which the model is feasible; use .maximize_growth() method to identify
            (should be <= mu_max output of .maxmimize_growth() method)
        tolerance: float; default 0
            Threshold below which expected sensitivity of solver is too low to detect infeasibility
        n_cores: int, default None
            the number of workers to use for parallelization
        graph: bool; default True
            plot the relationship between growth and the objective function of interest
        fig_name: str; default None
            save the plotted figure to 
            
        Returns
        ----------
        sln: tuple 
            first element is the growth value at which the non-growth objective is optimized
            second element is the optimized non-growth objective value
        predicted: pd.DataFrame
            1000 growth values between 0 and mu_max, with corresponding interpolated objective values
        interp_fit: output of scipy.interpolate.interp1d
            a function to interpolate objective values from growth values, used to generated predicted
        optimal_vals: collections.OrderedDict
            keys are n_points growth values between 0 and mu_max, values are the objective value optimized at 
            the corresponding growth value
        res: dict
            keys are n_points growth values between 0 and mu_max, values are the output of .solve_lp at 
            corresponding growth values with the objective set to the non-growth objective input
        '''
        
        if self.solver_ is None:
            warnings.warn('Solver is not initialized with ME_Model.intialize_solver, intializing with default parameters')
            self.initialize_solver()
        else:
            self.initialize_solver(solver_type = self.solver_type, precision = self.solver_precision)
        
        sln, predicted, interp_fit, optimal_vals, res = self.solver_.optimize(me_model = self, objective = objective, 
                                                        mu_max = mu_max, n_points = n_points, tolerance = tolerance, 
                                                        n_cores = n_cores, graph = graph, fig_name = fig_name)
        return sln, predicted, interp_fit, optimal_vals, res
    
    def infeasible_reactions(self, mu_val, sln, stat, tolerance = 1e-19):
        '''Binary search to find the maximum feasible growth rate

        Parameters
        ----------
        mu_val: float
            input growth value to ME_Model.solve_lp
        sln, stat: outputs of ME_Model.solve_lp
            Expected minimum feasible growth rate (~0)
        tolerance: float
            Threshold below which expected sensitivity of solver is too low to detect infeasibility

        Returns
        ----------
        ir: dict
            for reactions that cause infeasibility, keys are reaction ids for infeasible reactions and values are 
            difference by which reaction flux is infeasible
        '''
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
            
            counter = 0
            if math.isnan(flux): 
                ir[r.id] = flux
            elif (flux > ub + tolerance):
                ir[r.id] = abs(flux - ub)
            elif (flux < lb - tolerance):
                ir[r.id] = abs(lb - flux)
                
        if (len(ir)>0 and stat == 0) or (len(ir)==0 and stat != 0):
            warnings.warn('There is a discrepancy between the solver status and reactions that violate bound constraints')
        return ir
    
    
    def check_me_mass_balance(self):
        '''Checks that all reactions in ME Model are mass balance. Use after self.add_reactions'''
        print('Check reaction mass balances')

        for r in self.reactions:
            if len(r.check_mass_balance())>0:

                # account for instances when the original metabolic reaction is unbalances 
                bool1 = isinstance(r, Metabolic_Reaction)
                if bool1:
                    # account for instances when the ME_Reaction version of reversible reactions
                    if self.m_model.reactions.get_by_id(r.cobra_id).reversibility and r.id.endswith('_R'): # account for reversible reactions
                        bool1 = self.m_model.reactions.get_by_id(r.cobra_id).check_mass_balance() == {e: -c for e,c in r.check_mass_balance().items()}
                    else:
                        bool1 = r.check_mass_balance() == self.m_model.reactions.get_by_id(r.cobra_id).check_mass_balance()
                else:
                    bool1 = False

                bool2 = isinstance(r, Biomass_Reaction)
                if not (bool1 or bool2):
                    raise ValueError('Atleast one reaction is not mass balanced')
    
    def _check_complete_reactions(self):
        '''Checks that all the original metabolic model reactions have been included in the ME-Model'''
        if len(set([r.id for r in self.m_model.reactions]).difference([self.reactions.get_by_id(r_id).cobra_id for r_id in self.reaction_types['metabolism']]))>0:
            raise ValueError('Not all the original metabolic model reactions have been included in the ME-Model')


    def _check_hgncs(self):
        '''Checks that all reactions and macromolecules that are expected to have an assigned hgnc_id, do'''
        
        # reactions
        bool1 = len([r for r in self.reactions if not hasattr(r, 'hgnc_id') and not (isinstance(r, Biomass_Reaction) or isinstance(r, Metabolic_Reaction))]) > 0

        no_hgnc = [r for r in self.reactions if type(r) == Expression_Reaction and r.hgnc_id is None]
        no_hgnc = [r for r in no_hgnc if (r.subsystem not in ['tRNA_Biogenesis', 'rRNA_expression', 'Complex_Formation', 
                    'Complex_Degradation']) and (hasattr(r, 'ubiquitin_biogenesis') and not r.ubiquitin_biogenesis)]
        bool2 = len(no_hgnc)>0

        if bool1 or bool2:
            raise ValueError('An expression reaction does not have an hgnc_id')
        
        # macromolecules
        fragments = ['3_trailer', '5_leader', 'ets', 'its']
        exceptions = ['ubiquitin_monomer_protein_c', 'cleaved_polyubiquitin_moiety_protein_c', 
                      'ubiquitin_monomer_protein_n', 'cleaved_polyubiquitin_moiety_protein_n']
        hgnc_id_metabs = [m for m in self.metabolites if isinstance(m, Macromolecule) and m.hgnc_id is None and                 m.type not in ['trna', 'rrna', 'complex'] and                   not m.id.endswith('COMPLEX_enzyme_deg_proxy') and                  not (hasattr(m, 'fragment_type') and m.fragment_type in fragments) and                  not m.id in exceptions]
        if len(hgnc_id_metabs)>0:
            raise ValueError('Some macromolecules did not get an HGNC ID assigned')
            
    def _check_coupling(self, orphan = None, knock_out = list(), additional_ko = list()):
        '''Checks that all reactions have received appropriate machinery (compares coupled metabolites to GPR)
        
        Parameters
        ----------
        orphan: list
            List of reaction IDs in model for reactions that are not expected to have any machinery
            Defaults to self.reaction_types['orphan']
        knock_out: list
            List of HGNC IDs of knocked out genes
        additional_ko: list
             a list of HGNC IDs for genes that were not explicitly knocked-out, but were only involved in catalysis of 
            reactions catalyzed by a complex which contains another gene that was knocked-out
            this list is generated in build_me_model/me_builder
        
        '''
        print('Make sure all reactions received correct coupled machinery')
        
        # set arguments
        if orphan is None:
            if not 'orphan' in self.reaction_types:
                raise ValueError('Must specify a list of orphan reaction IDs')
            orphan = [self.reactions.get_by_id(r_id) for r_id in self.reaction_types['orphan']]
        if knock_out is None:
            knock_out = list()
        if additional_ko is None:
            additional_ko = list()
        
        # define list of ribosomal protein hgnc ids
        rbps = ['HGNC:10404', 'HGNC:10420', 'HGNC:10421', 'HGNC:10424', 'HGNC:10425', 'HGNC:18501', 'HGNC:10426', 
        'HGNC:10429', 'HGNC:10440', 'HGNC:10441', 'HGNC:10442', 'HGNC:10383', 'HGNC:10384', 'HGNC:10385', 
        'HGNC:10386', 'HGNC:10387', 'HGNC:10388', 'HGNC:10389', 'HGNC:10396', 'HGNC:10397', 'HGNC:10401', 
        'HGNC:10402', 'HGNC:10405', 'HGNC:10409', 'HGNC:10410', 'HGNC:10411', 'HGNC:10413', 'HGNC:10414', 
        'HGNC:10416', 'HGNC:18476', 'HGNC:10418', 'HGNC:10419', 'HGNC:3597', 'HGNC:10417', 'HGNC:10304', 
        'HGNC:10306', 'HGNC:10368', 'HGNC:10307', 'HGNC:10330', 'HGNC:10305', 'HGNC:10369', 'HGNC:10302', 
        'HGNC:10359', 'HGNC:10311', 'HGNC:10298', 'HGNC:10371', 'HGNC:10299', 'HGNC:10349', 'HGNC:10372', 
        'HGNC:10364', 'HGNC:10350', 'HGNC:21370', 'HGNC:10312', 'HGNC:10331', 'HGNC:10315', 'HGNC:10313', 
        'HGNC:10325', 'HGNC:10348', 'HGNC:10332', 'HGNC:10327', 'HGNC:10354', 'HGNC:10317', 'HGNC:10301', 
        'HGNC:10333', 'HGNC:10351', 'HGNC:12458', 'HGNC:17050', 'HGNC:10334', 'HGNC:17976', 'HGNC:10328', 
        'HGNC:10340', 'HGNC:10360', 'HGNC:17094', 'HGNC:10316', 'HGNC:10377', 'HGNC:10345', 'HGNC:13631', 
        'HGNC:10362', 'HGNC:10329', 'HGNC:10346', 'HGNC:10344', 'HGNC:10363', 'HGNC:10336', 'HGNC:10347', 
        'HGNC:10353']


        test_reactions = [r for r in self.reactions if not r.id in orphan]
        for r in tqdm(test_reactions):
            if isinstance(r, Metabolic_Reaction):
                r_ = self.m_model.reactions.get_by_id(r.cobra_id).copy()
            else: # Expression_Reaction
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

            if len(set(expected_machinery2).difference(knock_out + additional_ko)) == 0:
                ko = True

            if (not hasattr(r, 'coupled_metabolites') or 'catalysis' not in r.coupled_metabolites.values()) and not ko:
                raise ValueError('Reaction does not have a record of coupled machinery: ' + r.id)
            actual_machinery = [m for m,v in r.coupled_metabolites.items() if v == 'catalysis']


            if len(actual_machinery) > 0:
                translation = isinstance(actual_machinery[0], Ribosomal_Complex)
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
                            else:
                                am.append(p.id.split('_')[0])
                        actual_machinery = sorted(am)

                    # check that machinery matches
                    if len(expected_machinery) > 0: # non-dummy
                        err = True
                        for rm in expected_machinery:
                            if type(rm) != list: 
                                if not cplx and actual_machinery[0].id.split('_')[0] == rm:
                                    err = False
                            else:
                                if cplx and sorted(rm) == actual_machinery:
                                    err = False
                        if err:
                            if expected_machinery != knock_out:
                                raise ValueError('Machinery mismatch for ' + r.id)
                    else:
                        if len(actual_machinery) > 1 or actual_machinery[0].type != 'dummy_protein': # dummy
                            raise ValueError('Non-dummy protein coupled to deorphaned reaction')
                elif ribosomal_degradation:
                    am = list()
                    for am_ in actual_machinery:
                        for p in am_.decompose_complex():
                            if p.type != 'protein':
                                raise ValueError('Non-proteins in complex machinery for ' + r.id)
                            else:
                                am.append(p.id.split('_')[0])
                    actual_machinery = sorted(am)
                    if sorted(expected_machinery[0]) != actual_machinery:
                        raise ValueError('Incorrect machinery for ribosomal degradation: ' + r.id)
                elif translation:
                    actual_machinery = sorted([p.id.split('_')[0] for p in actual_machinery[0].decompose_complex() if p.type == 'protein'])
                    expected_machinery = sorted([p for p in expected_machinery[0] if p != 'ribosome'] + rbps)
                    if actual_machinery != expected_machinery:
                        raise ValueError('Incorrect machinery for translation: ' + r.id)
                else:
                    raise ValueError('Unaccounted for reaction criteria')
    



    def check_enzymes(self, _additional_ko = list()):
        '''Makes sure all genes being expressed participate in a catalysis reaction (no unecessary expression reactions)
        
        Paramaters
        ----------
        additional_ko: list
            a list of HGNC IDs for genes that were not explicitly knocked-out, but were only involved in catalysis of 
            reactions catalyzed by a complex which contains another gene that was knocked-out
            this list is generated in build_me_model/me_builder
        '''

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
        complexes = [m for m in complexes if not (isinstance(m, Ribosomal_Complex) or '_polyub_complex_' in m.id                                          or np.all([isinstance(r_, Complex_Degradation_Reaction) for r_ in m.reactions]))]
        if len(list(set(complexes).difference(active_complexes)))>0:
            raise ValueError('Unexpected inclusion of inactive complexes')


        # check monomers
        active_proteins += flatten_list([[p.id for p in m.decompose_complex()] for m in active_complexes])
        active_proteins = list(set([i.split('_')[0] if i.startswith('HGNC') else i for i in active_proteins]))
        proteins = set([i.split('_')[0] if 'HGNC' in i else i for i in proteins])
        proteins = [i for i in proteins if 'ubiquitin' not in i]
        if len(set(proteins).difference(active_proteins + _additional_ko).difference({'HGNC:12463', 'HGNC:12468'}))>0:
            raise ValueError('Unexpected inclusion of inactive protein monomers')
    
    def check(self, orphan = None, knock_out = list(), _additional_ko = list()):
        '''Check reaction coupling and mass balance
        
        Parameters
        ----------
        orphan: list
            List of reaction IDs in model for reactions that are not expected to have any machinery
            Defaults to self.reaction_types['orphan']
        knock_out: list
            List of HGNC IDs of knocked out genes
        additional_ko: list
             a list of HGNC IDs for genes that were not explicitly knocked-out, but were only involved in catalysis of 
            reactions catalyzed by a complex which contains another gene that was knocked-out
            this list is generated in build_me_model/me_builder
        '''
        self._check_complete_reactions()
        self._check_hgncs()
        self.check_me_mass_balance()
        self._check_coupling(orphan = orphan, knock_out = knock_out, additional_ko = _additional_ko)
        self.check_enzymes(_additional_ko = _additional_ko)

    def pickle(self, file = os.path.join(os.path.abspath(os.getcwd()), 'me_model.pickle')):
        '''Save ME_Model as a pickled object
        
        Parameters
        ----------
        file: str, default saves to current directory
            will save to file = "full/path/to/filename.pickle"
        
        '''
        if '.' in file:
            p = pathlib.Path(file)
            extensions = "".join(p.suffixes)
            file = str(p).replace(extensions, '.pickle')
        else:
            file = file + '.pickle'
        with open(file, 'wb') as handle:
            pickle.dump(self, handle)
    def copy(self):
        return copy.deepcopy(self)
    
    

