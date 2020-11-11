#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
# from utils.load_environmental_variables import *
from utils import machinery as mach
from utils import parameters as params
from utils import functions as func

from uniform_processes.build_trna_expression_reactions import charged_trna_metabolites, modified_trna_transcript_c
from expression.gene_information import gene_information


# In[2]:


charged_trna_map = {v.id.split('_')[2]: v for v in charged_trna_metabolites}
charged_trna_map_mw = {aa_code: aa_metab.formula_weight/1000 for aa_code, aa_metab in charged_trna_map.items()}
modified_trna_transcript_c_mw = modified_trna_transcript_c.formula_weight/1000


# In[ ]:


sp_dict = {1: True, 0: False, float('nan'): False}
ptm_cols = ['DSB', 'GPI', 'NG', 'OG']
ptm_keys = list(params.allowed_ptms.keys())

def generate_geneinfo_object(hgnc_id, psim = params.psim_me, machinery_list = mach.metabolic_machinery, 
                             metabolic_model = params.human_model):
    '''Generates gene information object from PSIM'''
    
    idx = psim[psim.HGNC_ID == hgnc_id].index.tolist()
    if len(idx) == 0:
        raise ValueError(hgnc_id + ' is not in the PSIM')
    if len(idx) > 1:
        warnings.warn('More than one entry of this gene by HGNC ID in PSIM, taking the first')

    entries = psim.loc[idx[0],:]
    if type(entries['LOCATION']) == str:
        entries['LOCATION'] = list(entries['LOCATION'].split(']')[0].split('[')[1].split(','))

    gene_info = gene_information(hgnc_id = entries['HGNC_ID'], 
                    premrna_seq = entries['PREMRNA_SEQ'], mrna_seq = entries['MRNA_SEQ'], 
                    protein_seq = entries['PROTEIN_SEQ'], 
                    machinery_list = machinery_list,
                    ptms = dict(zip(['dsb', 'og', 'gpi'],[entries['DSB'], entries['OG'], entries['GPI']])),
                    tmd = entries['TMD'], sp = entries['SP'], polyA_length = entries['POLYA_LENGTH'], 
                    n_introns = entries['N_INTRONS'])
    gene_info.get_final_locations(metabolic_model = metabolic_model, 
                                  final_locations = entries['LOCATION'])
    gene_info.check_gene_information()
    return gene_info

