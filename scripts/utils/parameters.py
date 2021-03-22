#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd 
import numpy as np
import cobra
from sympy.parsing.sympy_parser import parse_expr 

import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils.load_environmental_variables import processed_data_path, build_files_path


# In[2]:


psim_me = pd.read_hdf(processed_data_path + 'corrected_psim.h5', key = 'corrected')

psim_me['SP'] = psim_me['SP'].apply(lambda x: bool(x))
human_model = cobra.io.read_sbml_model(processed_data_path + 'corrected_model.xml')


# In[3]:


mu = parse_expr('mu')


# In[4]:


# compartments = human_model.compartments.copy()
# compartments['pm'] = 'plasma membrane'

compartments = {'c': 'cytoplasm', 'l': 'lysosome', 'r': 'endoplasmic reticulum', 'e': 'extracellular space', 'm': 'mitochondrion',
 'g': 'Golgi apparatus', 'n': 'nucleus', 'b': 'boundary', 'i': 'mitochondrial intermembrane space', 
 'x': 'peroxisome', 'pm': 'plasma membrane'}

allowed_ptms = {'dsb': 'disulfide bond formation', 'gpi': 'GPI Anchor', 'og': 'O-linked glycosylation'}#,
               #'ng': 'N-linked glycosylation'}

amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

allowed_trna_modifications = {}
#number_BiP = len(gene_info.protein_seq)/40


# In[5]:


# universal variables and inputs

rate_intron = 10/67000 # 10 introns / 67 kbp (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5199132/)
# L_polyA_n = 250 # https://www.nature.com/articles/s41592-019-0503-y
n_ub = 4 # see this - no. of ubiquitins to add to protein

transport_translocation_atp_cost = 0.5 # 1 ATP/2 residues
proteolysis_translocation_atp_cost = 0.5 # 1 ATP/2 residues

ptt_length = 160 # amino acid length greater than which co-translatioanl translocatoin occurs rather than post-translational
nuclear_diffusion_limit = 40 # 40 kDA and less proteins diffuse through nucleus
L_sp = 22 # secretory pathway signal peptide degradation
Kv = 0.7 # secretory pathway vesicle coat coefficients

membrane_diffusion_limit = 504 # 504 Da includes ATP, uncharged molecules at this diffusion limit are passive, no dummy


# In[6]:


# coupling parameters

# enzyme
keff_median = 3.983*3600 # units: hr^-1 (3.983 in s^-1)

# central dogma
alpha_m_median =  0.06108233261605428 # units: hours (Gregersen et al ) median value
alpha_p_median = 0.018342530808268292 # units: hours ^-1 (Cambridge et al 2011) median value
ptr_median = 65162.83940608428 # (Eraslan et al 2019) median value

ptr = pd.read_csv(build_files_path + 'PTR_Gagneur_processed.tsv', sep = '\t', index_col = 0)
# don't groupby hgnc ID median, because if tissue option is used, can include unmapped ids in calculation
ptr.drop(columns = ['ENSG_ID'], inplace = True)
ptr.columns = pd.Series(ptr.columns).apply(lambda x: x.split('_')[0] if '_PTR' in x else x).tolist()

alpha_p = pd.read_csv(build_files_path + 'protein_turnover.csv', index_col = 0)
alpha_p = alpha_p.groupby(alpha_p.HGNC_ID).median().kdeg # have true median stored above

alpha_m = pd.read_csv(build_files_path + 'Gregersen_mrna_turnover_processed.tsv', sep = '\t', index_col = 0)
alpha_m = alpha_m.groupby(alpha_m.HGNC_ID).median().median_turnover # have true median stored above

turnover = {'alpha_m': alpha_m, 'alpha_p': alpha_p, 
           'alpha_m_median': alpha_m_median, 'alpha_p_median': alpha_p_median}

# ribosome
rrna_degradation_constant = np.log(2)/72 # bioid 108025
# ribosomal_degradation_rate = np.log(2)/300 #bioid 110053 # unused


# In[7]:


# biomass

# constant fractions
dna_frac = 0.014
carb_frac = 0.071
lipid_frac = 0.097
# other_frac = 0.054

unmodeled_protein_frac = 1-0.12041534186261499

