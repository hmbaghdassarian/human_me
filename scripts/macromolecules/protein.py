#!/usr/bin/env python
# coding: utf-8

# In[31]:


import cobra
import sys
sys.path.insert(1, '../../scripts/')
from utils import parameters as params
from utils import metabolites as metab
# from utils import machinery as mach
# from uniform_processes.biomass import biomass_rna_mapper
# from utils import functions as func


# In[55]:


class Protein(cobra.Metabolite):
    def __init__(self, compartment, id_, gene_info = None, amino_acid_counts = None):
        '''
        
        Generates a cobra.Metabolite in the compartment for a protein with either 1) gene_info (gene_information object) or all of 
        2) id (string)and amino_acid_counts (dictionary, keys as 1-letter amino acide code values as number of 
        occurences in the protein).
        
        If gene_info and id_ are both not None, will concatenate the two strings.
        
        '''
        
        if compartment not in params.compartments.keys():
            err = 'Specified compartment is not considered in the ME Model. Please input one of the following: ' 
            err += ', '.join(list(params.compartments.keys()))
            raise ValueError(err)
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
        
        cobra.Metabolite.__init__(self, id = id_, compartment = compartment, charge = charge)
        self.elements = elements
        self.mass = self.formula_weight/1000
    def change_compartment(self, new_compartment):
        '''Returns a copy of the protein metabolite, but in new compartment'''
        
        if new_compartment == self.compartment:
            raise ValueError('The protein is already in this compartment')
        if new_compartment not in params.compartments.keys():
            err = 'Specified compartment is not considered in the ME Model. Please input one of the following: ' 
            err += ', '.join(list(params.compartments.keys()))
        
        new_protein = self.copy()
        new_protein.id = '_'.join(self.id.split('_')[:-1]) + '_' + new_compartment
        new_protein.compartment = new_compartment
        
        return new_protein


# In[61]:


# import random
# import cobra
# import pandas as pd
# from utils import parameters as params
# from utils import functions as func


# psim_toy = pd.DataFrame(columns = ['HGNC_ID', 'PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ', 'POLYA_LENGTH', 'TMD', 
#                                'SP', 'N_INTRONS', 'DSB', 'GPI', 'OG', 'LOCATION'])

# hgnc_id, premrna_seq = 'HGNC:TOY', ''.join(random.choices(['U', 'C', 'G', 'A'], k = 100))
# mrna_seq = premrna_seq[25:75]
# # note that there is no check that the protein_sequence corresponds to the mrna_sequence beyond checking for the length
# protein_seq = ''.join(random.choices(params.amino_acids, k = int(len(mrna_seq)/3)))
# polyA_length, tmd, sp, n_introns, dsb, gpi, og  = None, 1, True, 0, 2, 2, 2
# location = ['c'] # cytoplasm and golgi

# psim_toy.loc[0,:] = [hgnc_id, premrna_seq, mrna_seq, protein_seq, polyA_length, tmd, sp, n_introns, dsb, gpi, og, location]
# from expression.gene_information import gene_information
# gene_info = gene_information(hgnc_id, premrna_seq, mrna_seq, protein_seq,
#                  ptms = {'dsb': dsb, 'og': og, 'gpi': gpi}, tmd = tmd, sp = sp, polyA_length = polyA_length, 
#                  n_introns = n_introns)
# gene_info.get_final_locations(metabolic_model = cobra.Model(''), final_locations = location)

# a = func.make_protein_metabolite(id_ = gene_info.hgnc_id + '_unprocessed_folded', 
#                 amino_acid_counts = gene_info.amino_acid_counts, L_protein = gene_info.L_protein,
#                 compartment = 'r')
# b = Protein(compartment = 'r', id_ = 'unprocessed_folded', gene_info = gene_info)
# c = Protein(compartment = 'r', id_ = gene_info.hgnc_id + '_unprocessed_folded', 
#             amino_acid_counts = gene_info.amino_acid_counts)

