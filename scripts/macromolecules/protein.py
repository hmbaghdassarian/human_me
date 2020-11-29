#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
import sys
sys.path.insert(1, '../../scripts/')
from utils import metabolites as metab
from macromolecules.macromolecule import Macromolecule


# In[2]:


class Protein(Macromolecule):
    def __init__(self, compartment, id_, gene_info = None, amino_acid_counts = None):
        '''
        
        Generates a Macromolecule in the compartment for a protein with either 1) gene_info (gene_information object) or all of 
        2) id (string)and amino_acid_counts (dictionary, keys as 1-letter amino acide code values as number of 
        occurences in the protein).
        
        If gene_info and id_ are both not None, will concatenate the two strings.
        
        '''
        if gene_info is not None and (amino_acid_counts is not None):
            raise ValueError('Please specify either gene_info only or amino_acid_counts only')
        elif gene_info is None and ((id_ is None) or (amino_acid_counts is None)):
            raise ValueError('Please specify either gene_info or id_/amino_acid_counts')
        if id_ is None:
            raise ValueError('Unaccounted for condition in protein id naming')
        
        if gene_info is not None:
            L_protein = gene_info.L_protein
            id_ = gene_info.hgnc_id + '_' + id_ + '_protein_' + compartment 
            amino_acid_counts = gene_info.amino_acid_counts
        else:
            L_protein = sum(amino_acid_counts.values())
            id_ = id_ + '_protein_' + compartment
        
        charge = sum([metab.seq_amino_acid_map_compartments[compartment][aa_code].charge*aa_count for aa_code, aa_count in amino_acid_counts.items()])
        
        elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}
        if compartment in metab.seq_amino_acid_map_compartments.keys():
            for aa_code, aa_count in amino_acid_counts.items():
                aa_elements = metab.seq_amino_acid_map_compartments[compartment][aa_code].elements
                for element in aa_elements:
                    elements[element] += aa_count*aa_elements[element]
        else:
            raise ValueError('Internal: Must add ' + compartment + ' compartment to amino acid map in metab')

        # peptide bond formation
        elements['H'] -= 2*(L_protein-1)
        elements['O'] -= 1*(L_protein-1)
        
        Macromolecule.__init__(self, id = id_, compartment = compartment, charge = charge, elements = elements)
        self.type = 'protein'

