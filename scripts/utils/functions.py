#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
from cobra.core.gene import parse_gpr

import pandas as pd
import itertools
import ast

import os
import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
# from utils.load_environmental_variables import *
from utils import metabolites as metab
from utils import parameters as params


# # Functions

# In[2]:


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


# In[3]:


def get_reaction_compartment(reaction):
    '''Input is a cobra.Reaction, output is a singular compartment. This function maps reactions to a particular 
    compartment according to some rules'''
    
    compartments_ = list(reaction.compartments.copy())
    if len(compartments_) > 1: # for reactions that occur in more than one compartment
        if 'c' in compartments_ and len(compartments_) == 2: # remove cytoplasmic compartment between the two for machinery
            compartments_.remove('c')
        else: # choose most common compartment 
            compartments_ = [max(compartments_, key = compartments_.count)]
    if len(compartments_) > 1:
        raise ValueError('Failed to map reaction to a singular compartment')
    elif compartments_[0] not in params.compartments.keys():
        raise ValueError('Mapped reaction to a compartment that is not allowed in ME model')
    else:
        return compartments_[0]


# In[4]:


def hydrolyze_atp(rxn, n_atp, compartment):
    '''
    Rxn is a dict for the cobra.Reaction.add_metabolite function.
    n_atp is the # of atp to hydrolyze
    compartment is the compartment for hydrolysis
    
    '''
    n_atp = round(n_atp)
    
    if metab.atp_compartments[compartment] in rxn.keys():
        rxn[metab.atp_compartments[compartment]] -= n_atp 
    else:
        rxn[metab.atp_compartments[compartment]] = -n_atp 

    if metab.h2o_compartments[compartment] in rxn.keys():
        rxn[metab.h2o_compartments[compartment]] -= n_atp 
    else:
        rxn[metab.h2o_compartments[compartment]] = -n_atp 

    if metab.adp_compartments[compartment] in rxn.keys():
        rxn[metab.adp_compartments[compartment]] += n_atp 
    else:
        rxn[metab.adp_compartments[compartment]] = n_atp

    if metab.pi_compartments[compartment] in rxn.keys():
        rxn[metab.pi_compartments[compartment]] += n_atp 
    else:
        rxn[metab.pi_compartments[compartment]] = n_atp

    if metab.h_compartments[compartment] in rxn.keys():
        rxn[metab.h_compartments[compartment]] += n_atp 
    else:
        rxn[metab.h_compartments[compartment]] = n_atp
    
    return rxn


# In[ ]:


def get_base_counts_and_elements(seq, triphosphate = True):
    '''
    
    Inputs:
    1) Seq is a Bio.Seq object or a string representing an RNA sequence. 
    2) triphosphate is a boolean. If true, assume the 5' end has a triphosphate, otherwise assume monophosphate. 
   
   Outputs:
    1) base_counts is a dictionary with keys as one-letter RNA base characters and values as the # of 
    occurences of that base in the RNA sequence
    2) elements is a dictionary emulating cobra.Metabolite.elements
   
   '''
    base_counts = dict()
    for base_letter in metab.seq_element_map.keys():
        base_counts[base_letter] = seq.count(base_letter)
        
    elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0}
    for base_letter in metab.seq_element_map.keys():
        for element in elements.keys():
            elements[element] += base_counts[base_letter]* metab.seq_element_map[base_letter][element]   
    
    #3' OH end
    elements['H'] += 1 
    elements['O'] += 1
    
    # 5' end
    if triphosphate:
        elements['P'] += 2
        elements['O'] += 6
    else:
        elements['H'] += 1
      
        
    return base_counts, elements


# In[8]:


def make_protein_metabolite(id_, amino_acid_counts, L_protein, compartment):
    '''
    ID is a string to name the protein metabolite. 
    Amino acid counts is a dictionary with keys as the aa one letter code and counts as the number of occurences of that amino acid in the protein sequence
    L_protein is the length of the protein
    Compartment is the location of the protein (one letter string, corresponds to Recon2.2s compartments)
    
    Will return a cobra.Metabolite object with relevant charge and elements.
    
    '''
    
    protein_metabolite = cobra.Metabolite(id_ + '_protein[' + compartment + ']')
    protein_metabolite.compartment = compartment
    
    elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}
    if compartment in metab.seq_amino_acid_map_compartments.keys():
        for aa_code, aa_count in amino_acid_counts.items():
            aa_elements = metab.seq_amino_acid_map_compartments[compartment][aa_code].elements
            for element in aa_elements:
                elements[element] += aa_count*aa_elements[element]
    else:
        raise ValueError('Must add this compartment to make_protein_metabolite function')

    # peptide bond formation
    elements['H'] -= 2*(L_protein-1)
    elements['O'] -= 1*(L_protein-1)

    protein_metabolite.elements = elements
    # assume charge of amino acid is the ssame regardless of metabolite
    protein_metabolite.charge = sum([metab.seq_amino_acid_map_compartments[compartment][aa_code].charge*aa_count for aa_code, aa_count in amino_acid_counts.items()])
    return protein_metabolite


# In[12]:


def parse_me_reaction_id(x):
    if 'HGNC' in x.split('_')[0]:
        return '_'.join(x.split('_')[1:])
    else:
        return x


# In[13]:


def eval_complex_recur_full(expr):
    '''
    
    Recursive parsing of gprs into lists of complexes. Input expr is a cobra.parse_gpr(gpr_string), 
    output is a list of lists, each entry of which is a complex joined by 'AND'; netsted lists are joined 
    with each other by 'OR'. 
    
    Inspired by corda source code, should cite them.
    
    '''
    
    # corda: https://github.com/resendislab/corda/blob/master/corda/util.py
    if isinstance(expr, ast.Expression):
        return eval_complex_recur_full(expr.body) # Here!
    elif isinstance(expr, ast.Name):
        return cobra.core.gene.ast2str(expr)
    elif isinstance(expr, ast.BoolOp):
        op = expr.op
        if isinstance(op, ast.Or):
            return [eval_complex_recur_full(i) for i in expr.values] # Here!
        elif isinstance(op, ast.And):
            or_op = False
            names = []
            bool_ors = []
            bool_ands = []
            for e in expr.values:
                if isinstance(e, ast.BoolOp):
                    if isinstance(e.op, ast.Or):
                        or_op = True
                        bool_ors.append(e)
                    else:
                        bool_ands.append(e)
                elif isinstance(e, ast.Name):
                    names.append(eval_complex_recur_full(e))
            
            if or_op:
                product = []
                if len(bool_ors) > 1:
                    product = list(itertools.product(*[eval_complex_recur_full(i) for i in bool_ors]))
                    
                    
                ba_lists = []
                for ba in bool_ands:
                    ba_list = [eval_complex_recur_full(j) for i in bool_ands for j in i.values]
                    ba_lists.append(ba_list)
                
                result = []    
                if len(ba_lists) > 0:
                    bal_results = [eval_complex_recur_full(bal) for bal in ba_lists]
                    for br in bal_results:
                        if len(product) == 0:
                            result += [[eval_complex_recur_full(j)] + names + br for i in bool_ors for j in i.values] # Here!
                        else:
                            result += [list(p) + names + br for p in product]
                else:
                    if len(product) == 0:
                        result += [[eval_complex_recur_full(j)] + names for i in bool_ors for j in i.values] # Here!
                    else:
                        result += [list(p) + names for p in product]
                return result
                    
            else:
                return [eval_complex_recur_full(i) for i in expr.values] # Here!

def unnested_list(nested_list):
    lists = []
    non_lists = []
    for l in nested_list:
        if isinstance(l, list):
            lists.append(l)
        else:
            non_lists.append(l)
    size = len(lists)
    if size == len(nested_list):
        new_list = []
        for l in lists:
            new_list += unnested_list(l)
        return new_list
    elif size == 0:
        return [nested_list]
    else:
        return [item for sublist in [unnested_list(l) for l in lists] for item in sublist] + non_lists
    
def eval_complex(expr):
    '''Input is a gene reaction rule, output is a list. If the length of the list is longer than 1, these are
    due to the presence of OR. If the entry itself is a list longer than one, this entry is a complex.'''
    
    complexes = unnested_list(eval_complex_recur_full(parse_gpr(expr)[0]))
    
    # make sure machinery is always output in the same order
    for i in range(len(complexes)):
        if isinstance(complexes[i],list):
            complexes[i] = sorted(complexes[i])
    return complexes


# In[14]:


def SASA(mw):
    '''Estimate the protein solvent-accessible surface area from the molecular weight'''
    return mw**(0.75)


# In[15]:


# def check_me_mass_balance(r, metabolic_model = params.human_model):
#     '''r is a cobra.Reaction object'''
#     if len(r.genes) == 0:
#         return r.check_mass_balance()
#     else:
    
#         metabolic_reaction_names = [r.name for r in metabolic_model.reactions if len(r.genes)>0]
#         if r.name in metabolic_reaction_names: # metabolic reactions
#             # remove coupling constraint
#             rxn = {m:c for m,c in r.metabolites.items()  if ('protein' not in m.id) and ('complex' not in m.id)}
#         else: # expression reactions
#             raise ValueError('Do not currently have code base to get mass balance of expression reactions')
# #             if 'TRANSLATION' in r.id:
# #                 rxn = {m:c for m,c in r[0].metabolites.items()  if ('protein' not in m.id) and ('complex' not in m.id)}
# #                 rxn = {m:c for m,c in rxn.items()  if ('mrna' not in m.id) or ('proxy' in m.id)}
#     r_ = cobra.Reaction(' ')
#     r_.add_metabolites(rxn) 
    
#     return r_.check_mass_balance()


# In[16]:


import sympy
from operator import attrgetter
from collections import defaultdict
import warnings
from math import isinf
from six import iteritems


class ME_Reaction(cobra.Reaction):
    '''Inherited from cobra.Reaction. Allows stochiometric coefficient and reaction bounds to be a function of mu.'''

    def __init__(self,  id, type_, name='', subsystem='', lower_bound=0.0, upper_bound=None, 
                cobra_id = None):
        '''Helps distinguish between these reactions, which have mu in bounds, and coupling reactions, which 
        have mu in stochiometric coefficient.
        
        Cobra ID can be input to keep track of original cobra model id from which this reaction was derived.
        
        '''
        
        if len(set(type_).difference(['biomass', 'translation', 'catalysis'])) > 0:
            raise ValueError('Specified reaction type is not known')
        
        super().__init__(id, name, subsystem, lower_bound, upper_bound)
        self.type = type_
        self.cobra_id = cobra_id
    
    def check_me_bounds(self, lb, ub):
        if self.type == ['biomass']:
            if isinstance(lb, sympy.Expr) or isinstance(ub, sympy.Expr):
                if not params.mu in lb.free_symbols and params.mu in ub.free_symbols:
                    raise ValueError('Currently, if reaction bounds are a function of mu, they must be for both the upper and lower bound')
        else:
            if isinstance(lb, sympy.Expr) or isinstance(ub, sympy.Expr):
                raise ValueError('Reaction bounds can only be a function of mu for reactions of type biomass')

            
    def replace_bound_mu(self, mu_val = 1, values = None, _ = True, inplace = False):
        '''
        Assumes growth is always > 0. Gives numeric values to bounds for certain methods.
        
        '''
        
        if _:
            lb, ub = self._lower_bound, self._upper_bound
        else:
            lb, ub = self.lower_bound, self.upper_bound
            
        self.check_me_bounds(lb,ub)    
        
        if isinstance(lb, sympy.Expr):  # check_me_bounds makes sure both lb and ub are symp.Expr objects
            # replace growth with 1 (assuming growth always > 0)
            lb,ub = float(lb.subs(params.mu,mu_val)), float(ub.subs(params.mu,mu_val)) 
        else:
            if self.type == ['biomass']:
                warnings.warn('Bounds do not have a mu value')
        
        if values == None:
            if not inplace:
                return lb, ub
            else:
                self._lower_bound, self._upper_bound = lb,ub
        else:
            if not isinstance(values, list):
                raise TypeError('values must a list')
                
            for i in range(len(values)):
                if isinstance(values[i], sympy.Expr): # assumes the sympy expression always containts mu
                    values[i] = values[i].subs(params.mu, mu_val)
            if not inplace:
                return lb, ub, values
            else:
                raise ValueError('Either values must be None or inplace False')
    
    @property
    def reversibility(self):
        """Whether the reaction can proceed in both directions (reversible)

        This is computed from the current upper and lower bounds.

        """
        lb,ub = self.replace_bound_mu() 
        if not (isinstance(lb,sympy.Expr) and isinstance(ub,sympy.Expr)):
            return lb < 0 < ub
#         else: # if mu is just in one bound 
#             raise ValueError('For now, mu must be in both reaction bounds if boundaries are a function of mu')
    
    def build_reaction_string(self, use_metabolite_names=False):
        """Generate a human readable reaction string"""

        def format(number):
            return "" if number == 1 else str(number).rstrip(".") + " "

        id_type = 'id'
        if use_metabolite_names:
            id_type = 'name'
        reactant_bits = []
        product_bits = []
        for met in sorted(self._metabolites, key=attrgetter("id")):
            coefficient = self._metabolites[met]
            name = str(getattr(met, id_type))
            
            if not isinstance(coefficient, sympy.Expr):
                if coefficient >= 0:
                    product_bits.append(format(coefficient) + name)
                else:
                    reactant_bits.append(format(abs(coefficient)) + name)
            else:
                if float(coefficient.subs(params.mu, 1)) >= 0:
                    product_bits.append(format(coefficient) + name)
                else:
                    reactant_bits.append(format(coefficient).replace('-', '') + name)

        reaction_string = ' + '.join(reactant_bits)
        if not self.reversibility:
            lb,ub = self.replace_bound_mu(_ = False)
            if lb < 0 and ub <= 0:
                reaction_string += ' <-- '
            else:
                reaction_string += ' --> '
        else:
            reaction_string += ' <=> '
        reaction_string += ' + '.join(product_bits)
        return reaction_string 

    def check_mass_balance(self):
        """Compute mass and charge balance for the reaction

        returns a dict of {element: amount} for unbalanced elements.
        "charge" is treated as an element in this dict
        This should be empty for balanced reactions.
        """
        reaction_element_dict = defaultdict(int)
        for metabolite, coefficient in iteritems(self._metabolites):
            if not isinstance(coefficient, sympy.Expr): # don't include coupled metabolites
                if metabolite.charge is not None:
                    reaction_element_dict["charge"] +=                         coefficient * metabolite.charge
                if metabolite.elements is None:
                    raise ValueError("No elements found in metabolite %s"
                                     % metabolite.id)
                for element, amount in iteritems(metabolite.elements):
                    reaction_element_dict[element] += coefficient * amount
            else:
                if len(set(self.type).difference(['translation', 'catalysis'])) > 0:
                    raise ValueError('Mu can only be a coefficient in translation and catalysis reactions')
                else:
                    # for the exceptional situaiton in which a machinery is catalyzing its own expression reaction
                    # assume it is a reactant with coefficient -1...not robust
                    
                    # should only happen with peroxisomal protein degradation and 
                    if len(self.genes) == 1 and self.genes == metabolite.id.split('_')[0]:
                        for element, amount in iteritems(metabolite.elements):
                            reaction_element_dict[element] = -1*amount
                    
        # filter out 0 values
        return {k: v for k, v in iteritems(reaction_element_dict) if v != 0}
    
    def replace_coefficient_mu(self, mu_val):
        if len(set(self.type).difference(['translation', 'catalysis'])) > 0:
            raise ValueError('Mu can only be a coefficient in translation and catalysis reactions')
        
        if not (mu_val > 0):
            raise ValueError('Mu must be > 0')
        
#         # this method is quicker then looping through all metabolites - 
#         # specifically coded for current coupling format so not very robust
#         new_rxn = self.metabolites.copy()
#         if 'translation' in self.type:
#             mtr = [m for m in self.reactants if ('mrna_deg_proxy' in m.id) or ('mrna[c]' in m.id)]
#             for m in mtr:
#                 new_val = float(new_rxn[m].subs(params.mu, mu_val))
#                 new_rxn[m] = new_val
#         if 'catalysis' in self.type:
#             raise ValueError('Not yet encoded')
#         return self.add_metabolites(new_rxn, combine = False)
            
            
        new_rxn = self.metabolites.copy()
#         counter = 0
        for met, coeff in self.metabolites.items():
            if isinstance(coeff, sympy.Expr):
                new_rxn[met] = float(coeff.subs(params.mu, mu_val))
#                 counter += 1
        self.add_metabolites(new_rxn, combine = False)
        
#         if counter < 1:
#             warnings.warn('No mu values to replace')
        
    
    @property
    def reactants(self):
        """Return a list of reactants for the reaction."""
        
        reactants_ = list()
        for k,v in iteritems(self._metabolites):
            if not isinstance(v, sympy.Expr) and v < 0:
                reactants_.append(k)
            elif isinstance(v, sympy.Expr):
                if float(v.subs(params.mu, 1)) < 0:
                    reactants_.append(k)
        return reactants_

    @property
    def products(self):
        """Return a list of products for the reaction"""
        products_ = list()
        for k,v in iteritems(self._metabolites):
            if not isinstance(v, sympy.Expr) and v >= 0:
                products_.append(k)
            elif isinstance(v, sympy.Expr):
                if float(v.subs(params.mu, 1)) >= 0:
                    products_.append(k)
        return products_
    
#     def update_variable_bounds(self): ## needs to be modified for cobra.Model.add_reactions to work
#         if self.model is None:
#             return
#         # We know that `lb <= ub`.
#         # assume bounds as function of mu always > 0 so can add to the first conditional statemnt
#         if isinstance(self._lower_bound, sympy.Expr) or (self._lower_bound > 0): 
#             self.forward_variable.set_bounds(
#                 lb=self._lower_bound if (isinstance(self._lower_bound, sympy.Expr) or not isinf(self._lower_bound)) else None,
#                 ub=self._upper_bound if (isinstance(self._upper_bound, sympy.Expr) or not isinf(self._upper_bound)) else None
#             )
#             self.reverse_variable.set_bounds(lb=0, ub=0)
#         elif self._upper_bound < 0:
#             self.forward_variable.set_bounds(lb=0, ub=0)
#             self.reverse_variable.set_bounds(
#                 lb=None if isinf(self._upper_bound) else -self._upper_bound,
#                 ub=None if isinf(self._lower_bound) else -self._lower_bound
#             )
#         else:
#             self.forward_variable.set_bounds(
#                 lb=0,
#                 ub=None if isinf(self._upper_bound) else self._upper_bound
#             )
#             self.reverse_variable.set_bounds(
#                 lb=0,
#                 ub=None if isinf(self._lower_bound) else -self._lower_bound
#             )


# In[17]:


from cobra.core.dictlist import DictList
from cobra.util.context import get_context
from sympy import lambdify
import warnings
# import scipy
import numpy as np
import sys
sys.path.insert(1, '../../scripts/')
from me_solver import solve_me



# from sympy.parsing.sympy_parser import parse_expr
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
                        array[m_ind(metabolite), r_ind(reaction)] = stoich.subs(params.mu, mu_val)
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
        x: optimal solution
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
        
        return solve_me.solve_lp(me_model = self, mu_val = mu_val, objective = objective, 
                                 close_biomass_dilution = close_biomass_dilution,
                                 solver_type = solver_type, precision = precision)


# In[18]:


# # test LP
# test_model = ME_Model(cobra.io.load_json_model('/data2/hratch/Software/qminos_solver/solvemepy/examples/models/iJO1366.json'))
# # test_model.constraints[0].lb = -5
# # test_model.constraints[0].ub = 2

# # test_model.constraints[20].ub = 50
# # test_model.constraints[20].lb = 30
# xq,statq,hsq = test_model.solve_lp(mu_val = 0.03, objective = {'BIOMASS_Ec_iJO1366_core_53p95M': 1}, 
#                                   precision = 'quad')


# In[19]:


# test_model = cobra.io.load_json_model('/data2/hratch/Software/qminos_solver/solvemepy/examples/models/iJO1366.json')
# test_model.objective = {test_model.reactions.BIOMASS_Ec_iJO1366_core_53p95M: 1}


# In[20]:


# test_reaction = ME_Reaction('test', type_ = ['catalysis'])


# A, B, C = cobra.Metabolite('mA'), cobra.Metabolite('mB'), cobra.Metabolite('mC')
# rxn_A, rxn_B = ME_Reaction('rA', type_ = ['translation']), ME_Reaction('rB', type_ = ['translation'])
# rxn_C = ME_Reaction('rC', type_ = ['biomass'])
# rxn_A.add_metabolites({A: -2*params.mu, B: 3*params.mu, C: 2})
# rxn_B.add_metabolites({C: -5*params.mu, B: 10})
# rxn_C._lower_bound, rxn_C._upper_bound = params.mu, params.mu

# test_model = ME_Model('test')
# test_model.add_reactions([rxn_A, rxn_B, rxn_C])

# xq,statq,hsq = test_model.solve_lp(mu_val = 0.03, objective = {'rA': 1})

