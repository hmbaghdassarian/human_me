
import cobra
import itertools
import uuid

import sys
sys.path.insert(1, '../../scripts/')
from macromolecules.macromolecule import Macromolecule
from uniform_processes.biomass import biomass_mapper

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
        compartments = list(set([m.compartment for m in self.components]))
        if len(compartments) == 1:
            compartment = compartments[0]
        # exception of ribosome complex
        elif (sorted(compartments) == ['c', 'r']) and ('mature_ribosome_complex_c' in [m.id for m in self.components]):
            compartment = 'c'
        else:
            raise ValueError('Metabolites forming a complex must all be in the same compartment')
        
        
        # parse metabolite id
        if complex_id == None:
            self.temp_id = str(uuid.uuid4().fields[0])
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
#         self.mass = str(self.mass) # to avoid mitaking with non-complex macromolecules
    def form_complex(self, reaction_id = None):

        '''
        Output: A cobra.Reaction object representing the complex formation between metabolites stored in self.complex_formation

        '''
        
        if reaction_id is None:
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

    def udate_id(self):
        '''In cases where complex id is too long (see build_me_model generate_complex_reactions method)'''
        self.temp_id = str(uuid.uuid4().fields[0])
        self.id = self.temp_id + '_complex_' + self.compartment

    def get_complex_biomass(self):
        '''Returns a dictionary of the complex biomass by its individual component types'''

        biomass_by_type = dict()
        for m, count in self.decompose_complex().items():
            if m.type in biomass_by_type.keys():
                biomass_by_type[m.type] += count*m.mass
            else:
                biomass_by_type[m.type] = count*m.mass

        return biomass_by_type

def add_biomass_change(reaction):
    '''
    
    Input: list of cobra.Reactions
    Output: dictionary delineating the change in biomass (products - substrates) for the different categories
    of biomass.
    
    '''
    biomass_change = dict()
    for m, count in reaction.metabolites.items():
        if isinstance(m, Macromolecule):
            if type(m) != Complex:
                if m.type in biomass_change:
                    biomass_change[m.type] += (count*m.mass)
                else:
                    biomass_change[m.type] = (count*m.mass)
            else:
                for type_, mass_ in m.get_complex_biomass().items():
                    if type_ in biomass_change:
                        biomass_change[type_] += (count*mass_)
                    else:
                        biomass_change[type_] = (count*mass_)
    biomass_change = {biomass_mapper[k]:v for k,v in biomass_change.items()}
    reaction.add_metabolites(biomass_change, combine = False)
