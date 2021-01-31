#!/usr/bin/env python
# coding: utf-8

# In[8]:


import cobra
import itertools
import collections
import uuid
import numpy as np

import sys
sys.path.insert(1, '../../scripts/')
from macromolecules.macromolecule import Macromolecule
from utils import machinery as mach


# In[9]:


cotransloc_ids = set([mid + '_folded_protein_r' for mid in mach.ctnm + mach.translation_efs]) 
def flatten_list(list_):
    return [item for sublist in list_ for item in sublist]


# In[91]:


class Complex(Macromolecule):
    '''Complexes formed by non-covalent interactions between macromolecules'''
    def __init__(self, metabolites, complex_id = None, ignore_compartment = False):
        '''
        Parameters
        ----------
        metabolites: list 
            each entry is a Macromolecule object (protein or RNA or complex, not generic metabolites)
        complex_id: str
            the id of the complex metabolite; if None, will generate a random id
        ignore_compartment: bool 
            whether to ignore the metabolite compartments, mainly for internal use
        
        Returns
        ----------
        A Complex object representing the complex formed between input macromolecules
        
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
        
        # test compartment consistency - exception of ribosome complexes
        comp_ids = [m.id for m in self.components]
        cotransloc_cond = (len(cotransloc_ids.difference(comp_ids))==0) or ('mature_ribosome_complex_c' in comp_ids)
        
        if len(compartments) == 1:
            compartment = compartments[0]
        elif (sorted(compartments) == ['c', 'r']) and cotransloc_cond:
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
        
        self.reaction_id = None # none before running form_complex(); this is used in update_id() method
        self._deg_initialized = False
        self.enzyme = False
        
    def update_id(self, new_id = None):
        '''In cases where complex id is too long (see build_me_model generate_complex_reactions method)'''

        if self.reaction_id is not None:
            raise ValueError('Reaction and complex IDs will be consistent since you are updating the id after forming the reaction.')
        if new_id is None:
            self.temp_id = str(uuid.uuid4().fields[0])
        else:
            self.temp_id = new_id
        self.id = self.temp_id + '_complex_' + self.compartment
    
    
    def form_complex(self, reaction_id = None, reversible = True):

        '''The reaction required to generate the Complex object
        Note: assumes non-covalent complex formation (in terms of elemental balance)

        Parameters
        -----------
        reaction_id: str
            ID to assign to the reaction
        reversible: bool
            Whether the reaction is reversible. Setting to True may make model more efficient (allows reuse of 
            self.components if involved in other reactions)

        Returns
        ----------
        complex_formation: cobra.Reaction 
            the complex formation between metabolites stored in self.components

        '''
        
        if reaction_id is None:
            self.reaction_id = self.temp_id + '_COMPLEX_FORMATION' + self.compartment
        else:
            self.reaction_id = reaction_id + '_COMPLEX_FORMATION' + self.compartment        
        
        
        complex_formation = cobra.Reaction(self.reaction_id)
        rxn = {m: -count for m,count in self.components.items()}
        rxn[self] = 1
        complex_formation.add_metabolites(rxn)
        if reversible:
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

    def get_complex_biomass(self):
        '''Returns a dictionary of the complex biomass by its individual component types'''

        biomass_by_type = dict()
        for m, count in self.decompose_complex().items():
            if m.type in biomass_by_type.keys():
                biomass_by_type[m.type] += count*(m.formula_weight/1000)
            else:
                biomass_by_type[m.type] = count*(m.formula_weight/1000)

        return biomass_by_type

    def _initialize_deg_params(self):
        '''Initialize attributes for creating degradation reactions'''
        
        self._deg_initialized = True
        dc = self.decompose_complex()
        complex_types = list(set([m.type for m in dc]))
        if len(complex_types) > 1 or complex_types[0] != 'protein':
            raise ValueError('Complexes can only be degraded if composed only of proteins')
        
        self._amino_acid_counts = collections.Counter()
        self._ptms = collections.Counter()
        self._L_protein = 0
        for p,c in dc.items():
            self._L_protein += p._L_protein*c
            for i in range(c):
                self._amino_acid_counts.update(p._amino_acid_counts)
                if hasattr(p, '_ptms'):
                    self._ptms.update(p._ptms)
        
        if len(self._ptms) > 0:
            raise ValueError('PTMs to Complexes is currently unaccounted for and will likely lead to imbalances in degradation reactions')
            
        self._deg_id = self.temp_id
        self._degradation_reactions = []
        
        del dc
        del complex_types
        
    def _consolidate_degradation_rxns(self):
        '''Remove redundant IDs'''
        self._degradation_reactions = list(set(self._degradation_reactions))
    
    def change_compartment(self, new_compartment):
        '''Returns a copy of the complex metabolite, but in new compartment'''
        return self.change_compartment_and_components(new_compartment)
    def change_compartment_and_components(self, new_compartment):
        '''Returns a copy of the complex metabolite, but in new compartment. 
        Recursive to change all components (nested complexes and their components)'''


        if new_compartment == self.compartment:
            raise ValueError('The macromolecule is already in this compartment')
        if new_compartment not in params.compartments.keys():
            err = 'Specified compartment is not considered in the ME Model. Please input one of the following: ' 
            err += ', '.join(list(params.compartments.keys()))

        metabolites = list()
        for m,c in self.components.items():
            if m.type != 'complex':
                metabolites += [m.change_compartment(new_compartment)]*c
            else:
                metabolites += [m.change_compartment_and_components(new_compartment)]*c

        new_complex = Complex(metabolites = metabolites, complex_id = complex.temp_id)
        if self._deg_initialized:
            new_complex._initialize_deg_params()
        return new_complex


# In[92]:


import sys
import random
import pandas as pd
from expression.gene_information import gene_information
import expression.build_mrna_expression_reactions as build_mrna
from expression.protein_expression import ubiquitin
from macromolecules.protein import Protein
from utils import parameters as params

proteins = []
for i in range(2):
    psim_toy = pd.DataFrame(columns = ['HGNC_ID', 'PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ', 'POLYA_LENGTH', 'TMD', 
                                   'SP', 'N_INTRONS', 'DSB', 'GPI', 'OG', 'LOCATION'])

    hgnc_id, premrna_seq = 'HGNC:TOY', ''.join(random.choices(['U', 'C', 'G', 'A'], k = 100))
    mrna_seq = premrna_seq[25:75]
    # note that there is no check that the protein_sequence corresponds to the mrna_sequence beyond checking for the length
    protein_seq = ''.join(random.choices(params.amino_acids, k = int(len(mrna_seq)/3)))
    polyA_length, tmd, sp, n_exons, dsb, gpi, og  = None, 1, True, None, 2, 2, 2
    ub_args = ubiquitin.express_ubiquitin(compress_mrna = False)

    import itertools
    reactions = list()

    location = list(params.compartments.keys())
    psim_toy.loc[0,:] = [hgnc_id, premrna_seq, mrna_seq, protein_seq, polyA_length, tmd, sp, n_exons, dsb, gpi, og, location]
    gene_info = gene_information(hgnc_id, premrna_seq, mrna_seq, protein_seq,
                     ptms = {}, tmd = tmd, sp = sp, polyA_length = polyA_length, 
                     n_exons = n_exons) 
    gene_info.get_final_locations(metabolic_model = cobra.Model(''), final_locations = location)
    proteins.append(Protein(id_ = 'a', compartment = 'c', gene_info = gene_info))
    
    proteins_ = list()
    for p in proteins:
        proteins_.append(p.change_compartment('g'))
#     proteins = proteins_

cplx = Complex(metabolites = proteins_+proteins_, complex_id = 'test')
cplx = Complex(metabolites = [cplx, cplx] + proteins_ + proteins_, complex_id = 'test')
cplx._initialize_deg_params()


# In[93]:


new_complex = cplx.change_compartment('c')


# In[95]:


new_complex.components


# In[97]:


cplx.components


# In[ ]:





# In[69]:


def change_compartment(complex, new_compartment):
    '''Returns a copy of the macromolecule metabolite, but in new compartment. 
    Recursive to change all components (nested complexes and their components)'''
    
    
    if new_compartment == complex.compartment:
        raise ValueError('The macromolecule is already in this compartment')
    if new_compartment not in params.compartments.keys():
        err = 'Specified compartment is not considered in the ME Model. Please input one of the following: ' 
        err += ', '.join(list(params.compartments.keys()))

    metabolites = list()
    for m,c in complex.components.items():
        if m.type != 'complex':
            metabolites += [m.change_compartment(new_compartment)]*c
        else:
            metabolites += [change_compartment(m, new_compartment)]*c
    
    new_complex = Complex(metabolites = metabolites, complex_id = complex.temp_id)
    if complex._deg_initialized:
        new_complex._initialize_deg_params()
    return new_complex

    
    


# In[82]:





# In[72]:


new_complex = change_compartment(complex = cplx, new_compartment = 'c')


# In[74]:


new_complex._L_protein


# In[57]:


cplx.components


# In[51]:


cplx._L_protein


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


def add_complex_metabolites(cplx, met_to_add, complex_id):
    '''Add a metabolite to an existing complex object, returning it as a new complex
    
    Parameters
    -----------
    cplx: macromolecules.complex.Complex
        The Complex object to be appended
    met_to_add: dict
        The metabolites to add to the complex. Keys are objects inherited from 
        macrmolecules.macromolecule.Macromolecule and values are the number of copies of this metabolite to add.
    complex_id: str
        The new id for the new complex
    
    Returns
    ----------
    cplx2: macromolecules.complex.Complex
        The new Complex object with the added metabolites
    
    '''
    
    
    mtblts = list()
    for m,c in cplx.components.items():
        mtblts += [m]*c
    for m,c in met_to_add.items():
        mtblts += [m]*c

    cplx2 = Complex(metabolites = mtblts,
            complex_id = complex_id)
    if cplx._deg_initialized:
        cplx2._initialize_deg_params()
        cplx2._degradation_reactions += cplx._degradation_reactions # inherit degradation reactions
        cplx2._consolidate_degradation_rxns()
    return cplx2

