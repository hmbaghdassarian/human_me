#!/usr/bin/env python
# coding: utf-8

# In[2]:


import cobra
import pandas as pd

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import molecular_weight as calculate_molecular_weight

import requests, sys, json, re, warnings

import sys
sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *
from utils import *


# In[44]:


class gene_information():
    '''This class compiles all the necessary information for a given protein to be expressed in the 
    ME model. 
    
    Notes: As of right now, machinery PTMs are not considered. Proteins processed through the secretory pathway includes 
    1) anything assigned to the Secreted Protein module, or 2) machinery with a final location
    in the following compartments: ['l', 'r', 'e', 'x', 'g', 'pm']. See get_final_locations() method
    for details.'''
    
    
    def __init__(self, metabolic_model, hgnc_id, 
                 premrna_seq = None, mrna_seq = None, protein_seq = None,
                 ptms = {}, tmd = 0, sp = False, keff = None, polyA_length = None, n_introns = None):
        '''
        
        1) Metabolic model is a cobrapy model - required. 
        
        2) HGNC ID is a string in the format HGNC:#### - required.
        
        3-5) Relevant string representing sequence - required
        
        6) PTMs is a dictionary with keys as the ptm and values as the number of that ptm for that gene.
        PTMs are not currently considered for machinery. - optional
        
        7) TMD is an integer indicating the number of transmembrane domains the protein has. This is only relevant
        for proteins processed into secretory pathway. - optional
        
        8) SP is a boolean indicating whether a protein has a signal peptide. 
        Not used in current format - unimplemented
        
        9) keff is a float representing the kinetic constant the enzyme in [units]. - optional
        10) polyA_length is an integer representing the length of the polyA tail. This information will be estimated
        if not provided. - optional
        
        11) n_introns is an integer representing the length of the polyA tail. This information will be estimated
        if not provided. Should be specific to the transcript isoform. - optional
        
        '''
        
        self.hgnc_id = hgnc_id
        
        # current structure assumes that a protein is either machinery (catalyzing a reaction) or
        # a secreted protein (processed through secretory pathway, does not catalyze reaction) but not both
        
        machinery = [g.id for g in metabolic_model.genes] # not super efficient to do this each time
        if hgnc_id in machinery:
            self.module = 'Machinery'
        else:
            self.module = 'Non-Machinery'
        
        # sequence check
        if premrna_seq == None or mrna_seq == None or protein_seq == None:
            raise ValueError('All of the sequence types (premrna, mrna, protein) must be provided')
        if len(set(premrna_seq).difference(['A', 'U', 'G', 'C'])) > 0:
            raise ValueError('The premrna sequence contains bases which are not allowed')
        if len(set(mrna_seq).difference(['A', 'U', 'G', 'C'])) > 0:
            raise ValueError('The mrna sequence contains bases which are not allowed')
        if len(set(protein_seq).difference(amino_acids)) > 0:
            raise ValueError('The protein sequence contains amino acids which are not allowed')
        
        if len(premrna_seq) < len(mrna_seq):
            raise ValueError('The premrna sequence provided is shorter than the mrna sequence provided')
        elif len(premrna_seq) == len(mrna_seq):
            if premrna_seq != mrna_seq:
                raise ValueError('Premrna and mrna sequences are the same length, but not the same sequence')
            if self.n_introns > 0:
                raise ValueError('Premrna and mrna sequences are the same length, but you have indicated this is not an intronless gene')
        else:
            premrna_base_counts, mrna_base_counts = dict(), dict()
            for base_letter in seq_element_map.keys():
                premrna_base_counts[base_letter] = premrna_seq.count(base_letter)
                mrna_base_counts[base_letter] = mrna_seq.count(base_letter)
            for k,v in mrna_base_counts.items():
                if v > premrna_base_counts[k]:
                    raise ValueError('Number of ' + k + ' bases in premrna sequence less than that of mrna sequence')
            
            self.premrna_base_counts = premrna_base_counts
            self.mrna_base_counts = mrna_base_counts
            
            
        if len(mrna_seq) < len(protein_seq)*3:
            warnings.warn('The mrna and protein sequence lengths are inconsistent')
        
        self.premrna_seq = premrna_seq
        self.mrna_seq = mrna_seq
        self.protein_seq = protein_seq
        
        self.protein_mass = calculate_molecular_weight(seq=self.protein_seq, seq_type='protein')

        self.ptms = ptms
        self.tmd = tmd
        self.sp = sp
        
        if polyA_length == None or pd.isna(polyA_length) or (polyA_length >= 0 and round(polyA_length) == polyA_length):
            self.polyA_length = polyA_length
        else:
            raise ValueError('polyA_length must either be an integer >= 0 or None/nan')
        
        if n_introns == None or pd.isna(n_introns) or (n_introns >= 0 or round(n_introns) == n_introns):
            self.n_introns = n_introns
        else:
            raise ValueError('n_introns must either be an integer >= 0 or None/nan')
        
        self.keff = keff
        if self.keff == None and self.module == 'Machinery':
            warnings.warn('No keff specified for this enzyme, will assume a value in model building')
        elif self.keff != None and self.module == 'Non-Machinery':
            warnings.warn('keff specified for a non-machinery protein, will not be used')
            
            
        self.final_locations = None
       
    def get_final_locations(self, metabolic_model, final_locations = None):
        '''Assigns a final compartment for proteins. For machinery, extracts this from the model. 
        For secreted proteins, final_locations should be specified by a list of strings
        within the allowable compartments. This method helps define necessary transport reactions.
        
        The final output will be a dictionary with keys as the final locations and values as the method of 
        synthesis (Traditional Expression, Mitochondrial Expression, Canonical Secretion, Non-Canonical Secretion) 
        depending on Boolean rules. Traditional are those that don't go through the secretory pathway.'''
        
        if self.module == 'Machinery':
            if final_locations != None:
                warnings.warn('Final location extacted from cobrapy model, will disregard user input.')

            rxns = list(metabolic_model.genes.get_by_id(self.hgnc_id).reactions)
            final_locations = []
            
            for r in rxns:
            # proteins can be associated with multiple locations due to multiple reactions, but for each reaction
            # we want that protein to be associated with one compartment
                compartments = r.compartments.copy()
                if len(compartments) == 1: # not needed but more efficient
                    final_locations += list(compartments)
                    pass
                elif len(compartments) == 2: # for reactions that occur in more than one compartment
                    if 'c' in compartments: # remove cytoplasmic compartment between the two for machinery
                        compartments.remove('c')
                    else: # choose compartment on reactant side if no cytoplasmic compartment
                        reactant_compartments = set([m.compartment for m in r.reactants])
                        if len(reactant_compartments) == 1:
                            compartments = reactant_compartments
                        else:
                            compartments = max(reactant_compartments, key = list(reactant_compartments).count)
                elif len(compartments) > 2: # hardcoded for ASPGLUm reaction
                    compartments = {'i'}
                
                final_locations += list(compartments)
            final_locations = sorted(set(final_locations)) # redundancy from multiple reactions

                 
        if self.module == 'Non-Machinery':
            if final_locations == None:
                raise ValueError('For non-machinery, must specify the final locations')
            if type(final_locations) != list:
                raise ValueError('Final locations must be a list of string')
            if len(set(final_locations).difference(compartments.keys())) > 0:
                error = 'At least one of the locations specified is not allowed in this model.'
                raise ValueError(error + ' Allowable comparments include: ' + ', '.join(list(compartments.keys())))

   
        # transport rules
        # assume location dictates transport pathway ind of sp;
        # assume all genes are transported to mitochondria 
        # thus, two modes of transport:
        # 1) cytosolic transport: cytosolic translation-->import to final compartment
        # 2) canonical secretion: transport/translation via secretory pathway to final compartment
        
        # can expand on these based on signal peptide and transmembrane domain logic in the future
        self.final_locations = {}    
        for loc in final_locations: # no signal peptide consideration
            if loc in ['n', 'c', 'x', 'm', 'i']: 
                # mitochondrial expression not considered
                self.final_locations[loc] = 'Cytosolic Tranport'
                if self.sp: 
                    warnings.warng('Signal peptides not considered for these compartments')
            else:
                self.final_locations[loc] = 'Canonical Secretion'
                if not self.sp:
                    # add non-canonical in future
                    
                    # current structure assumes signal peptide presence for multi-localizing proteins with atleast
                    # one compartment in secretory pathway. in the future, presence of signal peptide could be 
                    # conditional for each location, somewhat analogous to transcript isoforms
                    
                    warning_ = 'Final location is part of secretory pathway, but no signal peptide indicated.'
                    warning_ += 'Non canonical secretion is not considered currently. Changing sp to True'
                    warnings.warn(warning_)
                    self.sp = True

    def check_gene_information(self):
        if self.final_locations == None:
            raise ValueError('Must specify a final location for the gene. Use the get_final_locations() method')
        if len(self.ptms) > 0:
            if self.module == 'Machinery':
                # change in the future
                warnings.warn('PTMs are not considered for machinery proteins currently')
            elif len(set(self.ptms.keys()).difference(allowed_ptms.keys())) > 0:
                warnings.warn('Atleast one of the PTMs provided will not be considered in this model')
        print('No errors raised')


# Usage

# In[21]:


# psim_me = pd.read_csv(local_data_path + 'processed/psim_me.csv', index_col = 0)
# human_model = cobra.io.load_json_model(local_data_path + 'processed/corrected_recon2_2.json')
# sp_dict = {1: True, 0: False}
# ptm_cols = ['DSB', 'GPI', 'NG', 'OG']
# ptm_keys = list(allowed_ptms.keys())


# In[22]:


# psim_me.head()


# In[23]:


# gene1_id = human_model.genes[0].id

# idx  = psim_me[psim_me['HGNC_ID'] == gene1_id].index
# ptms_ = dict(zip(ptm_keys, psim_me.loc[idx, ptm_cols].iloc[0,:].tolist()))
# ptms_ = {k:v for k,v in ptms_.items() if v != 0 and not pd.isna(v)}
# fl = psim_me.loc[idx, 'Location'].tolist()[0]

# pm,m,p = psim_me.loc[idx, 'PREMRNA_SEQ'].tolist()[0], psim_me.loc[idx, 'MRNA_SEQ'].tolist()[0], psim_me.loc[idx, 'PROTEIN_SEQ'].tolist()[0]

# sp = psim_me.loc[idx, 'SP'].tolist()[0]
# if pd.isna(sp):
#     sp = 0
# sp = sp_dict[sp]
# tmd = psim_me.loc[idx,'TMD'].tolist()[0]
# if pd.isna(tmd):
#     tmd = 0
# polyA_length_ = psim_me.loc[idx, 'POLYA_LENGTH'].tolist()[0]


# In[42]:


# # initialize the gene class
# gene1 = gene_information(human_model, hgnc_id = gene1_id, 
#                          premrna_seq=pm, mrna_seq=m, protein_seq=p,
#                          ptms = ptms_, tmd = tmd, sp = sp, 
#                         keff = None, polyA_length = polyA_length_, n_introns = None)
# print(gene1.module)
# print(gene1.hgnc_id)
# print(gene1.sp)
# print(gene1.ptms)
# print(gene1.tmd)
# print(gene1.polyA_length)
# print(gene1.n_introns)
# print(gene1.protein_mass)


# In[38]:


# # get the gene's final locations, final_locations list does not need to be specified for machinery
# gene1.get_final_locations(metabolic_model = human_model, final_locations=fl)
# print(gene1.final_locations)


# In[39]:


# gene1.check_gene_information()


# In[ ]:





# In[ ]:




