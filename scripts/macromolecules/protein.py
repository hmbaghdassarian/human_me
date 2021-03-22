#!/usr/bin/env python
# coding: utf-8

# In[31]:


import cobra
import sys
sys.path.insert(1, '../../scripts/')
from utils import metabolites as metab
from macromolecules.macromolecule import Macromolecule


# In[5]:


class Protein(Macromolecule):
    def __init__(self, compartment, id_, gene_info = None, amino_acid_counts = None, dummy = False):
        '''
        
        Generates a Macromolecule in the compartment for a protein with either 1) gene_info (gene_information object) or all of 
        2) id (string)and amino_acid_counts (dictionary, keys as 1-letter amino acide code values as number of 
        occurences in the protein).
        
        
        
        If gene_info and id_ are both not None, will concatenate the two strings.
        
        '''
        
        '''Inheritcs from Macromolecule. Class for Protein objects in ME-Model

        Parameters
        ----------
        compartment: str
            same as cobra.Metabolite.__init__
        id_: str
            same as cobra.Metabolite.__init__ (id)
        gene_info: gene_information object
        amino_acid_counts: dict
            keys are amino acids, values are the number of occurences in the protein sequence
        dummy: bool, default False
            whether the protein is a dummy protein for the unmodeled protein fraction of the ME-Model

        '''
        if gene_info is not None and (amino_acid_counts is not None):
            raise ValueError('Please specify either gene_info only or amino_acid_counts only')
        elif gene_info is None and ((id_ is None) or (amino_acid_counts is None)):
            raise ValueError('Please specify either gene_info or id_/amino_acid_counts')
        if id_ is None:
            raise ValueError('Unaccounted for condition in protein id naming')
        
        if gene_info is not None:
            # for degradation
            self.length = gene_info.L_protein
            self._amino_acid_counts = gene_info.amino_acid_counts
            self._ptms = gene_info.ptms
            self._deg_id = gene_info.hgnc_id
            # for rest of pipeline
            self.k_deg =  gene_info.coupling_params['alpha_p']
            id_ = gene_info.hgnc_id + '_' + id_ + '_protein_' + compartment 
            
        else:
            self._amino_acid_counts = amino_acid_counts
            self.length = sum(amino_acid_counts.values())
            id_ = id_ + '_protein_' + compartment            
        
        charge = sum([metab.seq_amino_acid_map_compartments[compartment][aa_code].charge*aa_count for aa_code, aa_count in self._amino_acid_counts.items()])
        
        elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}
        if compartment in metab.seq_amino_acid_map_compartments.keys():
            for aa_code, aa_count in self._amino_acid_counts.items():
                aa_elements = metab.seq_amino_acid_map_compartments[compartment][aa_code].elements
                for element in aa_elements:
                    elements[element] += aa_count*aa_elements[element]
        else:
            raise ValueError('Internal: Must add ' + compartment + ' compartment to amino acid map in metab')

        # peptide bond formation
        elements['H'] -= 2*(self.length-1)
        elements['O'] -= 1*(self.length-1)
        
        Macromolecule.__init__(self, id = id_, compartment = compartment, charge = charge, elements = elements)
        
        if not dummy:
            self.type = 'protein'
        else:
            self.type = 'dummy_protein'
        
        self.enzyme = False # whether the protein is involved in catalysis of a reaction
        self.keff = None
        self._degradation_reactions = [] # associated degradation reactions for protein monomer, if any
    
    def _consolidate_degradation_rxns(self):
        '''Remove redundant IDs'''
        self._degradation_reactions = list(set(self._degradation_reactions))


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

