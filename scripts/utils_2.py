#!/usr/bin/env python
# coding: utf-8

# In[2]:


# all scripts call utils, this is for utils that depend on earlier scripts
# can't be in utils
import sys
sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *
from utils import *
from build_trna_expression_reactions import charged_trna_metabolites, modified_trna_transcript_c
from gene_information import gene_information


# In[ ]:


charged_trna_map = {v.id.split('_')[2]: v for v in charged_trna_metabolites}


# In[ ]:


sp_dict = {1: True, 0: False, float('nan'): False}
ptm_cols = ['DSB', 'GPI', 'NG', 'OG']
ptm_keys = list(allowed_ptms.keys())

def generate_geneinfo_object(hgnc_id, psim = psim_me):
    '''Generates gene information object from PSIM'''
    
    idx = psim[psim.HGNC_ID == hgnc_id].index.tolist()
    if len(idx) == 0:
        raise ValueError(hgnc_id + ' is not in the PSIM')
    if len(idx) > 1:
        warnings.warn('More than one entry of this gene by HGNC ID in PSIM, taking the first')

    entries = psim.loc[idx[0],:]

    gene_info = gene_information(metabolic_model=human_model, hgnc_id = entries['HGNC_ID'], 
                    premrna_seq = entries['PREMRNA_SEQ'], mrna_seq = entries['MRNA_SEQ'], 
                    protein_seq = entries['PROTEIN_SEQ'], 
                    ptms = dict(zip(['dsb', 'og', 'gpi'],[entries['DSB'], entries['OG'], entries['GPI']])),
                    tmd = entries['TMD'], sp = entries['SP'], polyA_length = entries['POLYA_LENGTH'], 
                    n_introns = entries['N_INTRONS'])
    gene_info.get_final_locations(metabolic_model = human_model, final_locations = entries['LOCATION'])
    gene_info.check_gene_information()
    return gene_info

# def generate_geneinfo_object(hgnc_id, final_locations = [], psim = psim_me, keff = None, n_introns = None):
#     '''Generates gene information object from PSIM'''

#     idx  = psim[psim['HGNC_ID'] == hgnc_id].index
#     ptms_ = dict(zip(ptm_keys, psim.loc[idx, ptm_cols].iloc[0,:].tolist()))
#     ptms_ = {k:v for k,v in ptms_.items() if v != 0 and not pd.isna(v)}
#     fl = psim.loc[idx, 'Location'].tolist()[0]

#     pm,m,p = psim.loc[idx, 'PREMRNA_SEQ'].tolist()[0], psim.loc[idx, 'MRNA_SEQ'].tolist()[0], psim.loc[idx, 'PROTEIN_SEQ'].tolist()[0]

#     sp = psim.loc[idx, 'SP'].tolist()[0]
#     if pd.isna(sp):
#         sp = 0
#     sp = sp_dict[sp]
#     tmd = psim.loc[idx,'TMD'].tolist()[0]
#     if pd.isna(tmd):
#         tmd = 0
    
#     pa = psim_me.loc[idx, 'POLYA_LENGTH'].tolist()[0]
#     if pd.isna(pa):
#         polyA_length_ = pa
#     else:
#         polyA_length_ = round(float(pa))
#     gene_info = gene_information(metabolic_model = human_model, hgnc_id = hgnc_id, 
#                              premrna_seq=pm, mrna_seq=m, protein_seq=p,
#                              ptms = ptms_, tmd = tmd, sp = sp, 
#                             keff = keff, polyA_length = polyA_length_, n_introns= n_introns)
#     gene_info.get_final_locations(human_model, final_locations = final_locations)
#     gene_info.check_gene_information()
#     return gene_info

