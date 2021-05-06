#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra

import sympy
from operator import attrgetter
from collections import defaultdict
import warnings
from math import isinf
import numpy as np
from six import iteritems
import copy

import sys
sys.path.insert(1, '../../scripts/')
from utils import parameters as params
from utils import machinery as mach


# In[2]:


class ME_Reaction(cobra.Reaction):
    '''Inherited from cobra.Reaction. Allows stochiometric coefficient to be a function of mu.
    
    '''

    def __init__(self,  id, name='', subsystem='', lower_bound=0.0, upper_bound=None):
        
        '''
        Helps distinguish between these reactions, which have mu in bounds, and coupling reactions, which 
        have mu in stochiometric coefficient.
        '''
        
        super().__init__(id, name, subsystem, lower_bound, upper_bound)
        self.coupled_metabolites = dict()
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
        
        
        self.add_metabolites({metabolite: metabolite.coupling_coefficient[type]}, combine = True)
        mmap = {m.id: m for m in self._metabolites}
        if metabolite.id in mmap: # maintain the coupling attributes this way
            self._metabolites[metabolite] = self._metabolites.pop(mmap[metabolite.id])
            metabolite._reaction.add(self)
            
        # maintain the coupling attributes
#         if mmap[metabolite.id].coupling_coefficient is None:
#             mmap[metabolite.id].coupling_coefficient = metabolite.coupling_coefficient
#         else:
#             for type,coeff in metabolite.coupling_coefficient.items():
#                 if type not in mmap[metabolite.id].coupling_coefficient:
#                     mmap[metabolite.id].coupling_coefficient[type] = coeff
#                 elif  mmap[metabolite.id].coupling_coefficient[type] != coeff:
#                     raise ValueError('Mismatch in coupling coefficient calculation')

#         mmap[metabolite.id].enzyme = True
        
#         for attr in ['keff', ]
#         if mmap[metabolite.id].keff is None:
#             mmap[metabolite.id].keff = metabolite.keff

        if self.coupled_metabolites == dict():
            self.coupled_metabolites = {metabolite: type}
        else:
            self.coupled_metabolites[metabolite] = type

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

            
    def _check_me_bounds(self, lb, ub):
        if isinstance(lb, sympy.Expr) or isinstance(ub, sympy.Expr):
            raise ValueError('Reaction bounds can only be a function of mu for reactions of type biomass')
    
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
            if self.lower_bound < 0 and self.upper_bound <= 0:
                reaction_string += ' <-- '
            else:
                reaction_string += ' --> '
        else:
            reaction_string += ' <=> '
        reaction_string += ' + '.join(product_bits)
        return reaction_string

    def _map_coupled_metabolites(self):
        mmap = {m.id: m for m in self.metabolites}
        cm = dict()
        for md,type_ in self.coupled_metabolites.items():
#             mmap[md.id].coupling_coefficient = md.coupling_coefficient
            cm[mmap[md.id]] = type_
        self.coupled_metabolites = cm 
            
    def check_mass_balance(self, tol = 0, sympy_tol = 1e-15):
        """Compute mass and charge balance for the reaction

        returns a dict of {element: amount} for unbalanced elements.
        "charge" is treated as an element in this dict
        This should be empty for balanced reactions.
        
        sympy_tol: float
            sympy.expr conversions may result in some error, account for this when getting rid of the 
            coupling coefficient values
        """
        
        reaction_element_dict = defaultdict(int)
        mmap = {m.id: m for m in self._metabolites}
        
        md = dict()
        for m,c in self.metabolites.items():
            if m.id not in md:
                md[m.id] = c
            else:
                raise ValueError('Same metabolite id, different objects in reaction') #md[m.id] += c
        
        # deal with coupled metabolites (also required id mapping above)
        for metabolite, type in self.coupled_metabolites.items():
            md[metabolite.id] -= metabolite.coupling_coefficient[type] # coupling not part of mass balance
            md[metabolite.id] = float(md[metabolite.id])
            for val in [1,0]: 
                if abs((np.sign(md[metabolite.id])*val) - md[metabolite.id]) < 1e-15:
                    md[metabolite.id] = np.sign(md[metabolite.id])*val

        for m_id, coefficient in iteritems(md):   
            metabolite = mmap[m_id]
            if metabolite.charge is not None:
                reaction_element_dict["charge"] += coefficient * metabolite.charge
            for element, amount in iteritems(metabolite.elements):
                reaction_element_dict[element] += coefficient * amount

        return {k: v for k, v in iteritems(reaction_element_dict) if abs(v) > tol}
    
    def replace_coefficient_mu(self, mu_val, inplace = True):
        if not (mu_val > 0):
            raise ValueError('Mu must be > 0')
            
        new_rxn = self.metabolites.copy()
        for met, coeff in self.metabolites.items():
            if isinstance(coeff, sympy.Expr):
                new_rxn[met] = float(coeff.subs(params.mu, mu_val))
        
        if inplace:
            self.add_metabolites(new_rxn, combine = False)
        else: 
            reaction = copy.deepcopy(self)
            reaction.add_metabolites(new_rxn, combine = False)
            return reaction
        
    
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
        
class Metabolic_Reaction(ME_Reaction):
    '''Inherited from ME_Reaction, specifies the metabolic reactions in the model'''
    
    def __init__(self,  id, cobra_id, name='', subsystem='', lower_bound=0.0, upper_bound=None):
        '''cobra_id specifies the original reaction name in the M-Model'''

        super().__init__(id, name, subsystem, lower_bound, upper_bound)
        self.cobra_id = cobra_id


class Expression_Reaction(ME_Reaction):
    '''Inherited from ME_Reaction, specifies the expression reactions in the model'''
    
    def __init__(self,  id, subsystem, name='', lower_bound=0.0, upper_bound=None, 
                hgnc_id = None, 
                 synthesis = False, synthesis_type = None, sink = False, sink_type = None,
                 ubiquitin_biogenesis = False, ribosome_biogenesis = False):
        '''
        
        Parameters
        ----------
        subystem: str
            one of 'tRNA_Biogenesis', 'rRNA_expression', 'mRNA_expression', 'Protein_Expression', 'Protein_Degradation', 'Complex_Formation', 'Complex_Degradation'
        synthesis: bool
            whether the reaction represents the "main" synthesis/production for the macromolecule
            intended for use with genes (reactions with an associated hgnc id, and complexes)
        synthesis_type: str
            one of ['mRNA', 'protein', 'complex']
            if synthesis is True, the type of macromolecule being synthesized should also be specified
                *for mRNA, the synthesis reaction is coupled to its respective protein translation reaction
                *for proteins and complexes, the synthesis reaction is the final reaction producing the enzyme which 
                will be coupled to the metabolic catalysis reaction
        sink: bool 
            whether the reaction represents the "main" sink/degradation for the macromolecule
            intended for use with genes (reactions with an associated hgnc id, and complexes)
            
                *for mRNA, the degradation reaction will be coupled to the protein translation reaction
                *for proteins and complexes, the degradation reaction will be coupled to the respective
                metabolic catalysis reaction
                *exceptions are synthesis and sink of reactions in ubiquitin_biogenesis (True); these are 
                assigned as synthesis and sink to track ubiquitin, but are not themselves coupled to anything
        sink_type: str
            one of ['mRNA', 'protein', 'complex']
            if sink True, the type of macromolecule being degraded should also be specified
        ubiquitin_biogenesis: bool
            whether the Expression_Reaction is part of ubiquitin_biogenesis reactions, only used to ignore hgnc_id is None
        ribosome_biogenesis: bool
            whether the Expression_Reaction is part of ribosome_biogenesis reactions, only used to ignore hgnc_id is None
        '''
        
        if subsystem not in ['tRNA_Biogenesis', 'rRNA_expression', 'mRNA_expression', 'Protein_Expression', 
                             'Protein_Degradation', 'Complex_Formation', 'Complex_Degradation']:
            raise ValueError('Must specify an appropriate expression subsystem')
            
        super().__init__(id, name, subsystem, lower_bound, upper_bound)
        
        self.ubiquitin_biogenesis = ubiquitin_biogenesis
        if (not (self.subsystem in ['tRNA_Biogenesis', 'rRNA_expression', 'Complex_Formation', 
                                   'Complex_Degradation'] or self.ubiquitin_biogenesis)) and (hgnc_id is None):
            raise ValueError('Must specify hgnc_id of the gene being expressed')
        
        self.hgnc_id = hgnc_id
        self.synthesis = synthesis
        if self.synthesis and synthesis_type not in ['mRNA', 'protein', 'complex']: 
            raise ValueError('The synthesis type must be specified')
        else:
            self.synthesis_type = synthesis_type
    
        self.sink = sink
        if self.sink and sink_type not in ['mRNA', 'protein', 'complex']: 
            raise ValueError('The synthesis type must be specified')
        else:
            self.sink_type = sink_type
            
        self.ribosome_biogenesis = ribosome_biogenesis

class Protein_Expression_Reaction(Expression_Reaction):
    '''Inherited from Expression_Reaction, specifies the protein expression reactions in the model'''
    
    def __init__(self,  id, name='', lower_bound=0.0, upper_bound=None, 
                hgnc_id = None, translation = False, synthesis = False, 
                 ubiquitin_biogenesis = False, ribosome_biogenesis = False):
        '''
        
        Parameters
        ----------
        translation: bool
            whether the reaction represents the "main" synthesis/production for the protein 
            represents coupling of initial protein product to mRNA (mRNA-->protein coupling)
        synthesis: bool
            whether the reaction represents the "main" synthesis/production of an enzyme 
            that will be coupled to a metabolic reaction as a monomer 
        ubiquitin_biogenesis: bool
            whether the Expression_Reaction is part of ubiquitin_biogenesis reactions, only used to ignore hgnc_id is None
        ribosome_biogenesis: bool
            whether the Expression_Reaction is part of ribosome_biogenesis reactions, only used to ignore hgnc_id is None
        '''

        synthesis_type = None
        if synthesis:
            synthesis_type = 'protein'
        super().__init__(id=id, 
                         name=name, lower_bound=lower_bound, upper_bound=upper_bound, hgnc_id = hgnc_id, 
                         synthesis = synthesis, synthesis_type = synthesis_type,
                         ubiquitin_biogenesis = ubiquitin_biogenesis, ribosome_biogenesis = ribosome_biogenesis,
                         sink = False, sink_type = None, subsystem='Protein_Expression')
        self.translation = translation
        
        


class Protein_Degradation_Reaction(Expression_Reaction):
    def __init__(self, id, hgnc_id, sink = False, sink_type = None, name='', lower_bound=0.0, upper_bound=None):
        '''
        sink: bool
            whether the reaction represents the "main" sink/degradation for the macromolecule
            intended for use with genes (reactions with an associated hgnc id, and complexes)
        '''
        super().__init__(id=id, subsystem='Protein_Degradation', sink = sink, sink_type = sink_type, 
                         name=name, lower_bound=lower_bound, upper_bound=upper_bound, hgnc_id = hgnc_id)
        self._macromolecules = [] # list of macromolecule ids associated with this degradation reaction
        self._enzymes = None # list of enzyme ids associated with this degradation reaction
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

class Complex_Degradation_Reaction(Expression_Reaction):
    def __init__(self, id=None, sink = False, sink_type = None, 
                 name='', lower_bound=0.0, upper_bound=None, hgnc_id = None):
        '''
        
        sink: bool
            whether the reaction represents the "main" sink/degradation for the macromolecule
            intended for use with genes (reactions with an associated hgnc id, and complexes)
        hgnc_id: None
            always None, for internal use with expression/protein_expression/degradation script
        '''
        super().__init__(id=id, subsystem='Complex_Degradation', sink = sink, sink_type = sink_type, 
                         name=name, lower_bound=lower_bound, upper_bound=upper_bound)
        self._macromolecules = [] # list of macromolecule ids associated with this degradation reaction
        self._enzymes = None # list of enzyme ids associated with this degradation reaction
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
            
            # Option 1: degrade rRNA with ribosomal degradation - see also expression.protein_expression.degradation.proteasomal_degradation
            machinery_ = list(set(mach.proteasome_machinery + mach.exosome['HGNC ID (gene)'].tolist()))

#             # # Option 2: degrade proteins with ribosomal degradation, releasing rRNA as intact - see also expression.protein_expression.degradation.proteasomal_degradation
#             machinery_ = mach.proteasome_machinery

            
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
        


# In[3]:


def to_metabolic_reaction(reaction, id = None):
    '''Convert from cobra.Reaction to human_me.core.Metabolic_Reaction'''
    
    if id is None:
        id = reaction.id
    rxn = Metabolic_Reaction(id = id, cobra_id = reaction.id, name = reaction.name, subsystem = reaction.subsystem,
                             lower_bound = reaction.lower_bound, upper_bound = reaction.upper_bound)
    for k in set(reaction.__dict__.keys()).difference(['_id', 'name', 'subsystem', 
                                                       '_lower_bound', '_upper_bound', 
                                                      '_model']):
        rxn.__dict__[k] = copy.deepcopy(reaction.__dict__[k])
    return rxn
        


# In[4]:


class Biomass_Reaction(cobra.Reaction):
    '''Specifies biomass reactions in the model, allowing reaction bounds to be a function of mu'''
    
    def __init__(self,  id, name='', subsystem='', lower_bound=0.0, upper_bound=None):
        super().__init__(id, name, subsystem, lower_bound, upper_bound)

    def _check_me_bounds(self, lb, ub):
        if isinstance(lb, sympy.Expr) or isinstance(ub, sympy.Expr):
            if not params.mu in lb.free_symbols and params.mu in ub.free_symbols:
                raise ValueError('Currently, if reaction bounds are a function of mu, they must be for both the upper and lower bound')
    def replace_bound_mu(self, mu_val = 1, values = None, inplace = False, _ = True):
        '''
        Assumes growth is always > 0. Gives numeric values to bounds for certain methods.
        
        Parameters
        ----------
        mu_val: float
            The value for mu to replace the bounds that contain a mu expression with
        values: list or None
            Each entry is an expression containing mu to be replaced by mu_val; inplace must be False
        inplace: bool; default False
            Whether to replace the bounds inplace on the reaction object (True), or return the bounds
        _: bool
            internal use, whether to use cobra.Reaction._upper_bound or cobra.Reaction.upper_bound
        
        
        '''
        
        if _:
            lb, ub = copy.copy(self._lower_bound), copy.copy(self._upper_bound)
        else:
            lb, ub = copy.copy(self.lower_bound), copy.copy(self.upper_bound)
            
        self._check_me_bounds(lb,ub)    
        
        if isinstance(lb, sympy.Expr):  # _check_me_bounds makes sure both lb and ub are symp.Expr objects
            # replace growth with input mu val (assuming growth always > 0)
            lb,ub = float(lb.subs(params.mu,mu_val)), float(ub.subs(params.mu,mu_val)) 
#         else:
#             warnings.warn('Bounds do not have a mu value')
        
        if values == None:
            if not inplace:
                return lb, ub
            else:
                self._lower_bound, self._upper_bound = lb,ub
        else:
            if not isinstance(values, list):
                raise TypeError('values must a list')
                
            for i in range(len(values)):
                if isinstance(values[i], sympy.Expr): # assumes the sympy expression always contains mu
                    values[i] = float(values[i].subs(params.mu, mu_val))
            if not inplace:
                return lb, ub, values
            else:
                raise ValueError('Either values must be None or inplace False')
    
    @property
    def reversibility(self):
        """
        Whether the reaction can proceed in both directions (reversible)

        This is computed from the current upper and lower bounds.

        """
        lb,ub = self.replace_bound_mu() 
        return lb < 0 < ub
    
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
            if coefficient >= 0:
                product_bits.append(format(coefficient) + name)
            else:
                reactant_bits.append(format(abs(coefficient)) + name)

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

