#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle
import os
import pathlib

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
from utils import parameters as params
from macromolecules.complex import Complex
from macromolecules.protein import Protein

from core.reaction import ME_Reaction
from me_solver import solve_me


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
        self.orphan = None
        self.deorphaned = None
        

    def add_reactions(self, reaction_list):
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
                    stoichiometry = reaction._metabolites.pop(metabolite)
                    model_metabolite = self.metabolites.get_by_id(
                        metabolite.id)
                    reaction._metabolites[model_metabolite] = stoichiometry
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
                reaction_type = isinstance(reaction, ME_Reaction) and reaction.type != ['biomass']
                for metabolite, stoich in iteritems(reaction.metabolites):
                    if reaction_type and isinstance(stoich, sympy.Expr):
                        array[m_ind(metabolite), r_ind(reaction)] = float(stoich.subs(params.mu, mu_val))
                    else:
                        array[m_ind(metabolite), r_ind(reaction)] = stoich

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
    
    def solve_lp(self, mu_val, objective = {'biomass_dilution': 1}):
        '''Solves the linear program for a specified objective at a specified growth rate

        Parameters
        ----------
        mu_val: float
            The growth value for which to solve the linear program
        objective: dict, default {'bimoass_dilution': 1}
            The objective function to optimize. Dictionary represent a linear combination of reactions to optimize,
            with reaction ids as keys and the coefficient of the lin. comb. as the values. 

            Example: simplest case, to maximize reaction with id 'A', objective = {'A': 1}

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
            
        sln, stat, hs = self.solver_.solve_lp(me_model = self, mu_val = mu_val, objective = objective)
        return sln, stat, hs
    
    def maximize_growth(self, min_mu=0, max_mu=0.05, mu_accuracy=1e-4, increment = 1, verbose=True):
    
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
        mu_max,res = self.solver_.maximize_growth(me_model = self, 
                                                     min_mu=min_mu, max_mu=max_mu, 
                                                     mu_accuracy=mu_accuracy, increment = increment, 
                                                     verbose=verbose)
        return mu_max, res

    
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
            for reactions that cause feasibility, keys are reaction ids for infeasible reactions and values are optimizes reaction fluxes
        '''

        ir = dict()
        for r in self.reactions:
            flux = sln[self.reactions.index(r.id)]

            ub = copy.copy(r.upper_bound)
            lb = copy.copy(r.lower_bound)

            if isinstance(ub, sympy.Expr):
                ub = float(ub.subs(params.mu, mu_val))
            if isinstance(lb, sympy.Expr):
                lb = float(lb.subs(params.mu, mu_val))
            if math.isnan(flux) or ((flux > ub + tolerance) or (flux < lb - tolerance)):
                ir[r.id] = flux
                
        if (len(ir)>0 and stat == 0) or (len(ir)==0 and stat != 0):
            warnings.warn('There is a discrepancy between the solver status and reactions that violate bound constraints')
        return ir
    
    
    def check_me_mass_balance(self):
        '''Checks that all reactions in ME Model are mass balance. Use after self.add_reactions'''
        print('Check reaction mass balances')
        metabolic_reactions = [r.id for r in self.m_model.reactions]

        # strange exception ------
        exception = 'HGNC:9479_DEGRADATIONm' 
        check = [r for r in self.reactions if r.id != exception]
        if len([r for r in self.reactions if r.id == exception][0].check_mass_balance(tol = 1e-14))>0:
            err = True
        #---------------

        err = False
        for r in tqdm(check):
            if r.subsystem != 'Ribosome_Biogenesis' and r.subsystem != '':
                if isinstance(r, ME_Reaction):
                    if r.cobra_id is None and len(r.check_mass_balance())>0 and r.type != ['biomass']:
                        err = True
                        break
                    elif r.cobra_id is not None:
                        ogr = self.m_model.reactions.get_by_id(r.cobra_id).copy()
                        if (len([k for k in ogr.metabolites.keys() if k.elements is None]) == 0) and (r.check_mass_balance() != ogr.check_mass_balance()):
                            err = True
                            break
                else:
                    if r.id in metabolic_reactions:
                        ogr = self.m_model.reactions.get_by_id(r.id).copy()
                        if (len([k for k in ogr.metabolites.keys() if k.elements is None]) == 0) and (r.check_mass_balance() != ogr.check_mass_balance()):
                            err = True
                    elif len(r.check_mass_balance())>0:
                        err = True
            if err:
                raise ValueError('Not all expression module reactions are mass balanced') 
            
    def check_coupling(self):
        '''Checks that all reactions in ME Model received appropriate machinery. Use after self.add_reactions'''
        print('Check correct coupling of metabolic machinery')
        dummy = 'HGNC:DUMMY_folded_protein_c' in [m.id for m in self.metabolites]

        for r in tqdm(self.reactions):
            if isinstance(r, ME_Reaction):
                if r.type != ['biomass']:
                    if r.cobra_id is not None: 
                        r_ = self.m_model.reactions.get_by_id(r.cobra_id) # metabolic reactions
                    else:
                        r_ = r # non-metabolic reactions

                    machinery = [m for m,v in r.coupled_metabolites.items() if v == 'catalysis']
                    if len(machinery) > 1 and not r.id.endswith('_COMPLEX_PROTEASOMAL_DEGRADATIONc'):
                        raise ValueError('Unexpected coupling of multiple machinery: ' + r.id)
                    if not dummy:
                        if (len(r.genes) == 0 and len(machinery) > 0) or (len(r.genes) > 0 and len(machinery) == 0):
                            raise ValueError('Unexpected addition/lack of machinery: ' + r.id)
                    else:
                        if len(machinery) != 1 and not r.id.endswith('_COMPLEX_PROTEASOMAL_DEGRADATIONc'):
                            raise ValueError('Unexpected number of machinery coupled: ' + r.id)
                        machinery = machinery[0]
                        if ((len(r_.genes) == 0) and (machinery.id != 'HGNC:DUMMY_folded_protein_c')) or ((len(r_.genes) > 0) and (machinery.id == 'HGNC:DUMMY_folded_protein_c')):
                            raise ValueError('(A) Incorrect machinery coupled: ' + r.id)

                    if len(r.genes)>0:
                        if type(machinery) == list:
                            if len(machinery) != 1 and not r.id.endswith('_COMPLEX_PROTEASOMAL_DEGRADATIONc'):
                                raise ValueError('Unexpected number of machinery coupled: ' + r.id)
                            else:
                                machinery = machinery[0]

                        mach_me = [p for p in machinery.decompose_complex() if p.type == 'protein' or p.type == 'dummy_protein'] if isinstance(machinery, Complex) else [machinery]
                        mach_me = [m.id.split('_')[0] for m in mach_me]
                        mach_m = [g.id for g in list(r_.genes)]
                        if 'ribosome' in mach_m:
                            mach_m.remove('ribosome')
                            mach_m += [p.id.split('_')[0] for p in self.metabolites.get_by_id('mature_ribosome_complex_c').decompose_complex() if isinstance(p, Protein)]

                        if (len(set(mach_me).difference(mach_m)) > 0) or (len(set(mach_me).intersection(mach_m)) == 0):
                            raise ValueError('(B) Incorrect machinery coupled: ' + r.id)
                else:
                    pass
            else:
                if len(r.genes) > 0:
                    raise ValueError('Uncoupled reaction: ' + r.id)
                else:
                    if dummy and r.id not in self.orphan:
                        raise ValueError('Dummy protein did not couple: ' + r.id)
                    else:
                        pass

    def check(self):
        '''Check reaction coupling and mass balance'''
        self.check_me_mass_balance()
        self.check_coupling()

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
    
    

