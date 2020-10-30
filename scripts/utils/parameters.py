#!/usr/bin/env python
# coding: utf-8

# In[63]:


import pandas as pd 
import numpy as np
import cobra
from sympy.parsing.sympy_parser import parse_expr 

import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils.load_environmental_variables import processed_data_path


# In[ ]:


psim_me = pd.read_csv(processed_data_path + 'corrected_psim_me.csv', index_col = 0) 
psim_me['SP'] = psim_me['SP'].map({1: True, 0: False})
human_model = cobra.io.load_json_model(processed_data_path + 'corrected_model.json')


# In[ ]:


mu = parse_expr('mu')


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
# L_polyA_n = 250 # https://www.nature.com/articles/s41592-019-0503-y
n_ub = 4 # see this - no. of ubiquitins to add to protein

transport_translocation_atp_cost = 0.5 # 1 ATP/2 residues
proteolysis_translocation_atp_cost = 0.5 # 1 ATP/2 residues

ptt_length = 160 # amino acid length greater than which co-translatioanl translocatoin occurs rather than post-translational
nuclear_diffusion_limit = 40 # 40 kDA and less proteins diffuse through nucleus
L_sp = 22 # secretory pathway signal peptide degradation
Kv = 0.7 # secretory pathway vesicle coat coefficients


# In[2]:


# kinetic parameters

# enzyme
keff_median = 3.983*3600 # units: hr^-1 (3.983 in s^-1)


# central dogma
mrna_half_life = 10 #units: hours
alpha_p = 0.02 # units: hours ^-1

coupling_params = {'mrna_half_life': mrna_half_life, 
                  'alpha_p': alpha_p}


# In[5]:


# biomass

# constant fractions
dna_frac = 0.014
carb_frac = 0.071
lipid_frac = 0.097
other_frac = 0.054

unmodeled_protein_frac = 1-0.12041534186261499

