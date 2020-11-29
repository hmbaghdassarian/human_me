#!/usr/bin/env python
# coding: utf-8

# In[2]:


import cobra
import itertools
import uuid
from macromolecules.macromolecule import Macromolecule


# In[3]:


def flatten_list(list_):
    return [item for sublist in list_ for item in sublist]
    
class Complex(Macromolecule):
    def __init__(self, metabolites, complex_id = None):
        '''
        Inputs:
        metabolites is a list of Macromolecule objects (protein or RNA or complex, not generic metabolites)
        complex_id is a string for the id of the complex metabolite, otherwise will form a random id
        Output:
        A Macromolecule object representing the complex formed between macromolecules
        
        '''
        # checks
        if type(metabolites) != list or len(metabolites) == 0:
            raise ValueError('Must provide a list of macromolecules to form complex')
        # cobra metabolite not set up, check for bc they don't have the attribute .type
        if len([m for m in metabolites if not(isinstance(m, Macromolecule))]) > 0:
            raise ValueError('Generic cobra.Metabolite cannot form complexes with macromolecules currently')
        
        self.type = 'complex'
        self.components = {m: metabolites.count(m) for m in metabolites}
        # parse compartment    
        compartments = list(set([m.compartment for m in self.components.keys()]))
        if len(compartments) == 1:
            compartment = compartments[0]
#         # exception of ribosome complex
#         elif (len(compartments) == 2) and ('c' in compartments) and ('mature_ribosome_complex_c' in [m.id for m in self.components.keys()]):
#             compartment = 'c'
        else:
            raise ValueError('Metabolites are not in the same compartment')
        
        
        # parse metabolite id
        if complex_id == None:
            self.temp_id = str(uuid.uuid4().fields[0])
#             self.get_temp_id()        
        else: 
            self.temp_id = complex_id
        
        elements = dict()
        for m,count in self.components.items():
            for k,v in m.elements.items():
                if k in elements.keys():
                    elements[k] += v*count
                else:
                    elements[k] = v*count
        
        # make the metabolite
        Macromolecule.__init__(self, id = self.temp_id + '_complex_' + compartment, compartment = compartment,
                                  charge = sum([m.charge*count for m, count in self.components.items()]), 
                              elements = elements)
#         self.mass = self.get_complex_biomass() # to avoid mastking with non-complex macromolecules
    def form_complex(self, reaction_id = None):

        '''
        Output: A cobra.Reaction object representing the complex formation between metabolites stored in self.complex_formation

        '''
        
        if reaction_id == None:
            reaction_id = self.temp_id + '_COMPLEX_FORMATION' + self.compartment
        else:
            reaction_id = reaction_id + '_COMPLEX_FORMATION' + self.compartment        
        
        
        complex_formation = cobra.Reaction(reaction_id)
        rxn = {m: -count for m,count in self.components.items()}
        rxn[self] = 1
        complex_formation.add_metabolites(rxn)
        complex_formation.lower_bound = -1000 # reversible
        
        return complex_formation
    
    def decompose_complex(self, decomposed_complex = None):
        '''Recursive method to get the complex by its individual components, including nested complexes'''
        
        if decomposed_complex == None:
            all_metab = flatten_list([[m]*count for m, count in self.components.items()])
            decomposed_complex = Complex(metabolites = all_metab, complex_id = 'ignore')
        
        if 'complex' not in [m.type for m in decomposed_complex.components.keys()]:
            return decomposed_complex.components
        else:
            metabolites_ = flatten_list([[m]*count for m, count in decomposed_complex.components.items() if m.type != 'complex'])
            metabolites_ += flatten_list(flatten_list([[[m_]*count_ for m_, count_ in m.components.items()]*count for m, count in decomposed_complex.components.items() if m.type == 'complex']))
            return self.decompose_complex(decomposed_complex = Complex(metabolites = metabolites_, complex_id = 'ignore'))

#     def get_temp_id(self):
#         # to name complex -- does not currently name homodimers correctly (will give same name as dimer)
#         temp_id = '_'.join(sorted(set([k.id.split('_')[0] if 'HGNC:' in k.id else '_'.join(k.id.split('_')[:-1]) for k in self.decompose_complex().keys()])))
        
#         # deal with homo-oligomers to have unique IDs - likely won't come up 
#         all_metab = flatten_list([[m]*count for m,count in self.components.items()])
#         cpx_only = flatten_list([[m]*count for m,count in self.components.items() if type(m) == Complex])
#         if all_metab == cpx_only:
#             print('woo')
#             combs = itertools.combinations(test.components.keys(),2)
#             counter = 0
#             for comb in combs:
#                 if comb[0] == comb[1]:
#                     counter += 1
#             if counter == len(list(combs)):
#                 temp_id += '_' + str(uuid.uuid4().fields[-1])[:5] # add random id
        
#         self.temp_id = temp_id

    def get_complex_biomass(self):
        '''Returns a dictionary of the complex biomass by its individual component types'''

        biomass_by_type = dict()
        for m, count in self.decompose_complex().items():
            if m.type in biomass_by_type.keys():
                biomass_by_type[m.type] += count*m.mass
            else:
                biomass_by_type[m.type] = count*m.mass

        return biomass_by_type


# In[ ]:


# def get_complex_biomass_change(complex_products, complex_reactants):
#     '''Input is two lists of type COMPLEX, one representing those on the product side, one representing those on the reactant side
#     output is a dictionary of biomass change for each respective biomass type.'''
    
#     product_biomass = dict()
#     for cp in complex_products:
#         if type(cp)!= Complex:
#             raise TypeError('All complex products must be a COMPLEX object')
#         for bt, mw in cp.get_complex_biomass().items():
#             if bt in product_biomass.keys():
#                 product_biomass[bt] += mw
#             else:
#                 product_biomass[bt] = mw
    
#     reactant_biomass = dict()
#     for cr in complex_reactants:
#         if type(cr)!= Complex:
#             raise TypeError('All complex reactants must be a COMPLEX object')
#         for bt, mw in cr.get_complex_biomass().items():
#             if bt in reactant_biomass.keys():
#                 reactant_biomass[bt] += mw
#             else:
#                 reactant_biomass[bt] = mw
    
#     for bt in set(product_biomass.keys()).difference(reactant_biomass.keys()):
#         reactant_biomass[bt] = 0
#     for bt in set(reactant_biomass.keys()).difference(product_biomass.keys()):
#         product_biomass[bt] = 0    
    
#     return {bt: product_biomass[bt] - reactant_biomass[bt] for bt in product_biomass.keys() if product_biomass[bt] - reactant_biomass[bt] != 0}

