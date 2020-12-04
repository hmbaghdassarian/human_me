#!/usr/bin/env python
# coding: utf-8

# In[3]:


import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
# from utils.load_environmental_variables import *
from utils import machinery as mach
from utils import parameters as params
from utils import functions as func

from uniform_processes.build_trna_expression_reactions import charged_trna_metabolites, modified_trna_transcript_c
from expression.gene_information import gene_information


# In[4]:


charged_trna_map = {v.id.split('_')[2]: v for v in charged_trna_metabolites}


# In[23]:


ptm_cols = ['DSB', 'GPI', 'NG', 'OG']
ptm_keys = list(params.allowed_ptms.keys())
cp_keys = ['mrna_half_life', 'alpha_p', 'ptr', 'ptr_tissue', 'constant_ptr']

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
    
    cp_values = entries['MRNA_HALF_LIFE'], entries['ALPHA_P'], entries['PTR'], entries['PTR_TISSUE'], entries['CONSTANT_PTR']

    gene_info = gene_information(hgnc_id = entries['HGNC_ID'], 
                    premrna_seq = entries['PREMRNA_SEQ'], mrna_seq = entries['MRNA_SEQ'], 
                    protein_seq = entries['PROTEIN_SEQ'], 
                    machinery_list = machinery_list,
                    ptms = dict(zip(['dsb', 'og', 'gpi'],[entries['DSB'], entries['OG'], entries['GPI']])),
                    tmd = entries['TMD'], sp = entries['SP'], polyA_length = entries['POLYA_LENGTH'], 
                    n_introns = entries['N_INTRONS'], 
                    coupling_params = dict(zip(cp_keys, cp_values)))
    gene_info.get_final_locations(metabolic_model = metabolic_model, 
                                  final_locations = entries['LOCATION'])
    gene_info.check_gene_information()
    return gene_info

