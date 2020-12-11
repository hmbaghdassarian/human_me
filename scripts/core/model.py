#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
from cobra.core.dictlist import DictList
from cobra.util.context import get_context
from six import iteritems

import sympy
from sympy import lambdify
import warnings
import copy
import numpy as np

import sys
sys.path.insert(1, '../../scripts/')
from utils import parameters as params
from core.reaction import ME_Reaction
from me_solver import solve_me


# In[2]:


class ME_Model(cobra.Model):
    def __init__(self,  id_or_model, name = None):
        '''
        A simple object with an identifier
    
        Parameters
        ----------
        id: None or a string
            the identifier to associate with the object
            
        '''
        
        super().__init__(id_or_model, name)
        self.S = None

    def add_reactions(self, reaction_list):
        """Add reactions to the model.

        Reactions with identifiers identical to a reaction already in the
        model are ignored.

        The change is reverted upon exit when using the model as a context.

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

        context = get_context(self)

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
                    if context:
                        context(partial(
                            model_metabolite._reaction.remove, reaction))
        for gene in list(reaction._genes):
            # If the gene is not in the model, add it
            if not self.genes.has_id(gene.id):
                self.genes += [gene]
                gene._model = self

                if context:
                    # Remove the gene later
                    context(partial(self.genes.__isub__, [gene]))
                    context(partial(setattr, gene, '_model', None))

            # Otherwise, make the gene point to the one in the model
            else:
                model_gene = self.genes.get_by_id(gene.id)
                if model_gene is not gene:
                    reaction._dissociate_gene(gene)
                    reaction._associate_gene(model_gene)

        self.reactions += pruned

        if context:
            context(partial(self.reactions.__isub__, pruned))
        
        
        # make sure index and metabolite/list order are the same
        for idx, r_id in enumerate(self.reactions):
            if idx != self.reactions.index(r_id):
                raise ValueError('Indexing should be changed')

        for idx, m_id in enumerate(self.metabolites):
            if idx != self.metabolites.index(m_id):
                raise ValueError('Indexing should be changed')
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
        
    
    def solve_lp(self, mu_val, close_biomass_dilution = True, 
                 objective = {'biomass_dilution': 1}, solver_type = 'qminos', precision = 'quad'):
        
        '''
        
        mu_val is the growth value at which to optimize. 
        objective is a dictionary with keys as reaction ids to maximize as some linear combination and values as the coefficient for the linear objective
        Solver is a string, options of [qminos] - must have solveME and qMINOS installed
        
        Returns same outputs as qminospy.solver.solvelp:
        sln: optimal solution (reactions, metabolites, +1)
        stat: status
        hs: optimal basis
        
        
        stat:
        0     Optimal solution found.
        1     The problem is infeasible.
        2     The problem is unbounded (or badly scaled).
        3     Too many iterations.
        4     Apparent stall.  The solution has not changed
              for a large number of iterations (e.g. 1000).
        
        '''
        sln, stat, hs = solve_me.solve_lp(me_model = self, mu_val = mu_val, objective = objective, 
                                 close_biomass_dilution = close_biomass_dilution,
                                 solver_type = solver_type, precision = precision)
        
        return sln, stat, hs
    def infeasible_reactions(self, mu_val, sln, stat):
        '''
        Should only use for infeasible models to identify reactions that cause infeasibility.

        Inputs:
        sln: A solution output from ME_Model.solve_lp
        mu_val: The mu_value that was used in ME_Model.solve_lp

        Returns: a dictionary with keys as reaction ids and values as fluxes for reactions that caused infeasibility    
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
            if (flux > ub) or (flux < lb):
                ir[r.id] = flux
                
        if (len(ir)>0 and stat == 0) or (len(ir)==0 and stat != 0):
            warnings.warn('There is a discrepancy between the solver status and reactions that violate bound constraints')
        return ir

