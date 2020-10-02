#!/usr/bin/env python
# coding: utf-8

# In[63]:


import pandas as pd 
import numpy as np
import cobra

import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils.load_environmental_variables import *


# In[ ]:


# inputs
psim_me = pd.read_csv(local_data_path + 'processed/corrected_psim_me.csv', index_col = 0) #psim_me = pd.read_csv(root_path+'TRASH.csv', index_col = 0)
psim_me['SP'] = psim_me['SP'].map({1: True, 0: False})
human_model = cobra.io.load_json_model(local_data_path + 'processed/corrected_model.json')


# In[ ]:


compartments = {'c': 'cytosol',  'l': 'lysosome', 'm': 'mitochondria', 'r': 'endoplasmic reticulum', 
                'e': 'extracellular space', 'x': 'peroxisome/glyoxysome', 'n': 'nucleus', 'g': 'golgi apparatus',
                'i': 'inner mitochondrial compartment', 'pm': 'plasma membrane'}
allowed_ptms = {'dsb': 'disulfide bond formation', 'gpi': 'GPI Anchor', 'og': 'O-linked glycosylation'}#,
               #'ng': 'N-linked glycosylation'}

amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

allowed_trna_modifications = {}
#number_BiP = len(gene_info.protein_seq)/40


# In[6]:


# universal variables and inputs

rate_intron = 10/67000 # 10 introns / 67 kbp
n_ub = 4 # see this - no. of ubiquitins to add to protein

transport_translocation_atp_cost = 0.5 # 1 ATP/2 residues
proteolysis_translocation_atp_cost = 0.5 # 1 ATP/2 residues

ptt_length = 160 # amino acid length greater than which co-translatioanl translocatoin occurs rather than post-translational
nuclear_diffusion_limit = 40 # 40 kDA and less proteins diffuse through nucleus
L_sp = 22 # secretory pathway signal peptide degradation
Kv = 0.7 # secretory pathway vesicle coat coefficients


# In[65]:


# kinetic parameters

# enzyme
keff_median = 3.983*3600 # units: hr^-1 (3.983 in s^-1)


# central dogma
mrna_half_life = 10 #units: hours
alpha_m = np.log(2)/mrna_half_life
alpha_p = 0.02 # units: hours ^-1


# In[ ]:


# biomass

# constant fractions
dna_frac = 0.014
carb_frac = 0.071
lipid_frac = 0.097
other_frac = 0.054

