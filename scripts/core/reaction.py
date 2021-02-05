#!/usr/bin/env python
# coding: utf-8

# In[2]:


import cobra

import sympy
from operator import attrgetter
from collections import defaultdict
import warnings
from math import isinf
from six import iteritems

import sys
sys.path.insert(1, '../../scripts/')
from utils import parameters as params
from utils import machinery as mach


# In[2]:


class ME_Reaction(cobra.Reaction):
    '''
    
    Inherited from cobra.Reaction. Allows stochiometric coefficient and reaction bounds to be a function of mu.
    Note: not used in all expression module reactions, just those requiring mu as a parameter.
    
    '''

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
        self.coupled_metabolites = None
        self._protein_deg_proxy = False
    
    def _couple(self, metabolite, type):
        '''Add coupling coefficient and associated metadata to reaction for a coupled metabolite
        
        Parameters
        ----------
        metabolite: macromolecules.macromolecule.Macromolecule
            a macromolecule with an associated coupling coefficient
        type: str
            the type of reactions that are being coupled (one of ['mrna_degradation', 'mrna_formation', 'catalysis', 
            'enzyme_degradation'])
        
        '''
        
        if metabolite.coupling_coefficient is None:
            raise ValueError('Cannot add coupling metadata to reaction for a metabolite without coupling coefficient metadata')
        if type not in metabolite.coupling_coefficient: # this also checks correct coupling types defined
            raise ValueError('Incorrect coupling coefficient type specified for this metabolite')

        if self.coupled_metabolites is None:
            self.coupled_metabolites = {metabolite: type}
        else:
            self.coupled_metabolites[metabolite] = type

        self.add_metabolites({metabolite: metabolite.coupling_coefficient[type]}, combine = True)
    
    def couple(self, metabolites, types):
        '''Add coupling coefficient and associated metadata to reaction for a coupled metabolites
        
        Parameters
        ----------
        metabolites: macromolecules.macromolecule.Macromolecule or list 
            list of macromolecules or single macromolecule with associated coupling coefficients
        types: str or list
            the type of reactions that are being coupled (options: ['mrna_degradation', 'mrna_formation', 'catalysis'])
        
        '''

        if isinstance(metabolites, list):
            for metabolite, type in dict(zip(metabolites,types)).items():
                self._couple(metabolite,type)
        else:
            self._couple(metabolites,types)
        
            
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
            # replace growth with input mu val (assuming growth always > 0)
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
                    values[i] = float(values[i].subs(params.mu, mu_val))
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

    def check_mass_balance(self, tol = 0):
        """Compute mass and charge balance for the reaction

        returns a dict of {element: amount} for unbalanced elements.
        "charge" is treated as an element in this dict
        This should be empty for balanced reactions.
        """
        
        reaction_element_dict = defaultdict(int)
        md = self._metabolites.copy()
        if self.coupled_metabolites is not None:
            for metabolite, type in self.coupled_metabolites.items():
                md[metabolite] -= metabolite.coupling_coefficient[type] # coupling not part of mass balance
        for metabolite, coefficient in iteritems(md):    
            if metabolite.charge is not None:
                reaction_element_dict["charge"] += coefficient * metabolite.charge
            for element, amount in iteritems(metabolite.elements):
                reaction_element_dict[element] += coefficient * amount

        return {k: v for k, v in iteritems(reaction_element_dict) if abs(v) > tol}


    
    def replace_coefficient_mu(self, mu_val):
        if len(set(self.type).difference(['translation', 'catalysis'])) > 0:
            raise ValueError('Mu can only be a coefficient in translation and catalysis reactions')
        
        if not (mu_val > 0):
            raise ValueError('Mu must be > 0')
            
        new_rxn = self.metabolites.copy()
        for met, coeff in self.metabolites.items():
            if isinstance(coeff, sympy.Expr):
                new_rxn[met] = float(coeff.subs(params.mu, mu_val))
        self.add_metabolites(new_rxn, combine = False)
        
    
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
    def _add_protein_deg_proxy(self, protein_deg_proxy):
        '''Proxy metabolite for protein degradation'''
        if not protein_deg_proxy.type == 'proxy':
            raise ValueError('Expected proxy macromolecule')
        if self._protein_deg_proxy:
            raise ValueError('Protein degradation proxy already added')
            
        self.add_metabolites({protein_deg_proxy: 1})
        self._protein_deg_proxy = True
        self.protein_deg_proxy = protein_deg_proxy


# In[ ]:


class Protein_Degradation_Reaction(cobra.Reaction):
    def __init__(self, id=None, name='', subsystem='', lower_bound=0.0, upper_bound=None):
        cobra.Reaction.__init__(self, id=id, name=name, subsystem=subsystem, lower_bound=lower_bound, 
                                upper_bound=upper_bound)
        self._macromolecules = [] # list of macromolecule ids associated with this degradation reaction
        self._enzymes = None # list of enzyme ids associated with this degradation reaction
        self.sink = False # whether the reaction is the "final" (to amino acids) degradation reaction
        self.subsystem = 'Protein_Degradation'
        self._ribosomal_degradation = False # see complex_degradation_reaction for details
        
    def _update_tracking(self, macromolecules):
        '''Mutual tracking of degradation reactions associated with a macromolecule and vice-versa'''
        if type(macromolecules) != list:
            macromolecules._degradation_reactions.append(self.id)
            self._macromolecules.append(macromolecules)
        else:
            for macromolecule in macromolecules:
                macromolecule._degradation_reactions.append(self.id)
                self._macromolecules.append(macromolecule)
    def _consolidate_macromolecules(self):
        '''Remove redundant macromolecules'''
        for m in self._macromolecules:
            m._consolidate_degradation_rxns()
        self._macromolecules = list(set(self._macromolecules))
        
    def _update_enzymes(self):
        '''Update enzymes list to include macromolecules that are classified as enzymes'''
        self._enzymes = [m for m in self._macromolecules if m.enzyme]
        for m in self._enzymes:
            if self.id not in m._degradation_reactions:
                raise ValueErorr('Improper tracking of degradation reactions and associated macromolecules')
    def _set_proteasomal_degradation(self, **kwargs):
        '''For code consistency, mainly for Complex_Degradation_Reaction, see that method'''
        
        self.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)

class Complex_Degradation_Reaction(cobra.Reaction):
    def __init__(self, id=None, name='', subsystem='', lower_bound=0.0, upper_bound=None):
        cobra.Reaction.__init__(self, id=id, name=name, subsystem=subsystem, lower_bound=lower_bound, 
                                upper_bound=upper_bound)
        self._macromolecules = [] # list of macromolecule ids associated with this degradation reaction
        self._enzymes = None # list of enzyme ids associated with this degradation reaction
        self.sink = False # whether the reaction is the "final" (to amino acids) degradation reaction
        self.subsystem = 'Complex_Degradation'
        self._ribosomal_degradation = False
        
    def _update_tracking(self, macromolecules):
        '''Mutual tracking of degradation reactions associated with a macromolecule and vice-versa'''
        if type(macromolecules) != list:
            macromolecules._degradation_reactions.append(self.id)
            self._macromolecules.append(macromolecules)
        else:
            for macromolecule in macromolecules:
                macromolecule._degradation_reactions.append(self.id)
                self._macromolecules.append(macromolecule)
    def _consolidate_macromolecules(self):
        '''Remove redundant macromolecules'''
        for m in self._macromolecules:
            m._consolidate_degradation_rxns()
        self._macromolecules = list(set(self._macromolecules))
        
    def _update_enzymes(self):
        '''Update enzymes list to include macromolecules that are classified as enzymes'''
        self._enzymes = [m for m in self._macromolecules if m.enzyme]
        for m in self._enzymes:
            if self.id not in m._degradation_reactions:
                raise ValueErorr('Improper tracking of degradation reactions and associated macromolecules')
    def _set_proteasomal_degradation(self, macromolecule, ribosomal_complex):
        '''Quick addition of attribute for build_me script, since current format has the machinery for
        the proteosomal degradation different than standard complexes (to degrad rRNAs as well). 
        Change in machinery hard-coded into degradation.degrade script and double-checked in build_me script.
        
        macromolecule: Protein or Complex instance
        ribosomal_complex: bool
            whether or not the macromolecule is a ribosomal complex
        '''
        
        if not ribosomal_complex:
            machinery_ = mach.proteasome_machinery
        else:
            self._ribosomal_degradation = True
            # hard-coded
            machinery_ = list(set(mach.proteasome_machinery + mach.exosome['HGNC ID (gene)'].tolist()))

            # this is more flexible, but hardcoded check in degradation.proteasomal_degradation
            # renders this unecessary (keep for future iterations)
            
#             # add machinery
#             counter = 0
#             machinery_ = list()
#             if len(set([m for m in mdc if m.type == 'protein'])) > 0:
#                 machinery_ += mach.proteasome_machinery
#                 counter += 1
#             if len(rm) > 0:
#                 machinery_ += mach.exosome['HGNC ID (gene)'].tolist() #rrna degradation machinery
#                 counter += 1
        
#             if counter != 2:
#                 err = 'Internal: Only expect mature ribosome complex with rRNA and protein to be degraded. 
#                 err += 'Should work with current code, but double check'
#                 raise ValueError(err)
       
        self.gene_reaction_rule = ' and '.join(machinery_)
        

