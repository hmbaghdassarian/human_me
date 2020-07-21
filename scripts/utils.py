#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
import pandas as pd

import sys
sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *


# Necessary variables and functions shared across scripts

# Gene information variables

# In[2]:


# define necessary variables 
compartments = {'c': 'cytosol',  'l': 'lysosome', 'm': 'mitochondria', 'r': 'endoplasmic reticulum', 
                'e': 'extracellular space', 'x': 'peroxisome/glyoxysome', 'n': 'nucleus', 'g': 'golgi apparatus',
                'i': 'inner mitochondrial compartment', 'pm': 'plasma membrane'}
allowed_ptms = {'dsb': 'disulfide bond formation', 'gpi': 'GPI Anchor', 
               'ng': 'N-linked glycosylation', 'og': 'O-linked glycosylation'}

amino_acids = ['A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 
              'P', 'S', 'T', 'W', 'Y', 'V']

human_model = cobra.io.load_json_model(local_data_path + 'processed/corrected_recon2_2.json')
atp_n = human_model.metabolites.get_by_id('atp[n]')
gtp_n = human_model.metabolites.get_by_id('gtp[n]')
seq_metabolite_map = {human_model.metabolites.get_by_id('utp[n]'): 'U' , 
                      gtp_n: 'G',
                      human_model.metabolites.get_by_id('ctp[n]'): 'C',
                      atp_n: 'A'}

# RNA backbone elements
seq_element_map = dict()
for k,v in seq_metabolite_map.items():
    elements = k.elements
    elements['O'] = elements['O'] - 7 # lost from incoming ntp
    elements['P'] = elements['P'] - 2 # lost from incoming ntp
    elements['H'] = elements['H'] - 1 # lost from 3' end of growing strand
    seq_element_map[v] = elements


# mrna_expression variables

# In[3]:


ppi_n =  human_model.metabolites.get_by_id('ppi[n]')


# machinery

# elongation machinery------------------------------------------------------------
rnap = pd.read_csv(local_data_path + 'raw/RNAP_HUGO.csv', index_col = None, skiprows = [0])
rnap2 = rnap[rnap['Approved name'].isin([i for i in rnap['Approved name'] if ' II ' in i])]

tfiis, tfiif, ell = ['HGNC:11612', 'HGNC:11614'], ['HGNC:4652', 'HGNC:4653'], ['HGNC:23114', 'HGNC:17064', 'HGNC:23113']
elongin = pd.read_csv(local_data_path + 'raw/elongin.csv', index_col = None, skiprows = [0])
elongator = pd.read_csv(local_data_path + 'raw/elongator.csv', index_col = None, skiprows = [0])
fact = ['HGNC:11327', 'HGNC:11465']
ec = rnap2['HGNC ID (gene)'].tolist() + elongin['HGNC ID (gene)'].tolist() + elongator['HGNC ID (gene)'].tolist()
ec += tfiis + tfiif + ell + fact

# processing variables------------------------------------------------------------
# L_polyA_n = 250 # https://www.nature.com/articles/s41592-019-0503-y
from polyA_statistics import calculate_polyA_length#, min_polyA_mean, polyA_params, polyA_mod, 
pi_n =  human_model.metabolites.get_by_id('pi[n]')
h_n =  human_model.metabolites.get_by_id('h[n]')
h2o_n =  human_model.metabolites.get_by_id('h2o[n]')
gp = gtp_n.elements
gp['O'] -= 6
gp['P'] -= 2
amet_n =  human_model.metabolites.get_by_id('amet[n]')
ahcys_n =  human_model.metabolites.get_by_id('ahcys[n]')
adp_n = human_model.metabolites.get_by_id('adp[n]')
rate_intron = 10/67000 # 10 introns / 67 kbp

# processing machinery
cpsf = ['HGNC:2324', 'HGNC:2327', 'HGNC:19124', 'HGNC:2325', 'HGNC:2326', 'HGNC:25651', 'HGNC:13871']
cstf = ['HGNC:2483', 'HGNC:2484', 'HGNC:2485']
cfim, cfiim = ['HGNC:14981', 'HGNC:15970', 'HGNC:14982'], ['HGNC:30097', 'HGNC:16999']
polyA = cpsf + cstf + cfim + cfiim

nelf = ['HGNC:12768', 'HGNC:24324', 'HGNC:15934', 'HGNC:13974']
capping = nelf + ['HGNC:10073', 'HGNC:10075', 'HGNC:21077', 'HGNC:7658', 'HGNC:7659', 'HGNC:11467', 'HGNC:11469',
                 'HGNC:29200', 'HGNC:17970']

# ignore snRNA for now, can change in the future
spliceosome = pd.read_csv(local_data_path + 'raw/spliceosome.txt', index_col = None, sep = '\t')
spliceosome = spliceosome[spliceosome['Locus type'] != 'RNA, small nuclear']
spliceosome = sorted(set(spliceosome['HGNC ID'].tolist()))


# lariat degradataion------------------------------------------------------------
exosome = pd.read_csv(local_data_path + 'raw/exosome.csv', index_col = None, skiprows = [0])
lariat_machinery = {'Linearization': ['HGNC:15594'] ,
                    "5' Degradation": ['HGNC:12836'], 
               "Exosome": exosome.loc[:, 'HGNC ID (gene)'].tolist() + ['HGNC:29911'], 
               'NEXT Complex': ['HGNC:18734', 'HGNC:9904', 'HGNC:25265']}
lm1 = ' and '.join(lariat_machinery['Linearization'] + lariat_machinery["5' Degradation"])
lm2 = ' and '.join(lariat_machinery['Linearization'] + lariat_machinery["Exosome"] + lariat_machinery["NEXT Complex"])
lm_rule = '({})'.format(lm1) + ' or ' + '({})'.format(lm2)

nmp_map_n = {'C': human_model.metabolites.get_by_id('cmp[n]'), 
          'U': human_model.metabolites.get_by_id('ump[n]'), 
          'G': human_model.metabolites.get_by_id('gmp[n]'), 
          'A': human_model.metabolites.get_by_id('amp[n]')}
ntp_map_n = {v: k for k,v in seq_metabolite_map.items()}

# mrna export------------------------------------------------------------
tho = pd.read_csv(local_data_path + 'raw/tho.csv', index_col = None, skiprows = [0])
trex = tho.loc[:, 'HGNC ID (gene)'].tolist() + ['HGNC:17821', 'HGNC:25407', 'HGNC:24971', 'HGNC:24511', 'HGNC:24432',
                                               'HGNC:23782', 'HGNC:29093', 'HGNC:3447', 'HGNC:8071', 'HGNC:15913',
                                               'HGNC:25091', 'HGNC:24101', 'HGNC:7658']

# mrna degradation------------------------------------------------------------
h_c =  human_model.metabolites.get_by_id('h[c]')
h2o_c =  human_model.metabolites.get_by_id('h2o[c]')
pi_c = human_model.metabolites.get_by_id('pi[c]')
amet_c = human_model.metabolites.get_by_id('amet[c]')
ahcys_c = human_model.metabolites.get_by_id('ahcys[c]')
amp_c = human_model.metabolites.get_by_id('amp[c]')
nmp_map_c = {'C': human_model.metabolites.get_by_id('cmp[c]'), 
          'U': human_model.metabolites.get_by_id('ump[c]'), 
          'G': human_model.metabolites.get_by_id('gmp[c]'), 
          'A': amp_c}
ndp_map_c = {'C': human_model.metabolites.get_by_id('cdp[c]'), 
          'U': human_model.metabolites.get_by_id('udp[c]'), 
          'G': human_model.metabolites.get_by_id('gdp[c]'), 
          'A': human_model.metabolites.get_by_id('adp[c]')}

ccr4_not = pd.read_csv(local_data_path + 'raw/CCR4_NOT.csv', index_col = None, skiprows = [0])
# pabp3 isoform not included
deadenylation_machinery = {'CCR4_NOT Deadenylation': ccr4_not.loc[:, 'HGNC ID (gene)'].tolist(), 
                'PARN Deadenylation': ['HGNC:8609'], 
               'PABP Deadenylation': ['HGNC:20074', 'HGNC:29991', 'HGNC:8554']}
deadenylation_machinery = [item for sublist in [v for v in deadenylation_machinery.values()] for item in sublist]


mrna_degradation_machinery_1 = {"Exosome": exosome.loc[:, 'HGNC ID (gene)'].tolist(), 
               'Cap_Degradation': ['HGNC:29812']}
mrna_degradation_machinery_1 = [item for sublist in [v for v in mrna_degradation_machinery_1.values()] for item in sublist]
decapping_degradation_machinery = {'LSM1-7 Complex': ['HGNC:20472', 'HGNC:13940', 'HGNC:17874', 'HGNC:17259', 
                                   'HGNC:17162', 'HGNC:17017', 'HGNC:20470'],
                'Decapping': ['HGNC:18714', 'HGNC:24451', 'HGNC:24452', 'HGNC:17157'], 
               "5' Exonuclease": ['HGNC:30654']}
decapping_degradation_machinery = [item for sublist in [v for v in decapping_degradation_machinery.values()] for item in sublist]
degradation_rule1 = ' and '.join(deadenylation_machinery + mrna_degradation_machinery_1)
decapping_rule = ' and '.join(deadenylation_machinery + decapping_degradation_machinery)


# rrna_expression variables

# In[4]:


# variables needed
ntp_map_c = {'C': human_model.metabolites.get_by_id('ctp[c]'), 
             'U': human_model.metabolites.get_by_id('utp[c]'), 
             'G': human_model.metabolites.get_by_id('gtp[c]'), 
             'A': human_model.metabolites.get_by_id('atp[c]')}
ntp_map_n = {v:k for k,v in seq_metabolite_map.items()}

# 5s rrna
rnap3 = rnap[rnap['Approved name'].isin([i for i in rnap['Approved name'] if ' III ' in i])]
tfiiia, tfiiib = ['HGNC:4662'], ['HGNC:13652', 'HGNC:11551', 'HGNC:11588']
tfiiic = ['HGNC:4664', 'HGNC:4665', 'HGNC:4666', 'HGNC:4667', 'HGNC:4668', 'HGNC:20872']
rnap3_transcription_machinery = rnap3['HGNC ID (gene)'].tolist() + tfiiia + tfiiib + tfiiic
REXO5 = 'HGNC:24661'
xpo1 = ['HGNC:12825']

# other rrnas
rnap1 = rnap[rnap['Approved name'].isin([i for i in rnap['Approved name'] if ' I ' in i])]
taf = ['HGNC:11532', 'HGNC:11533', 'HGNC:11534']
ubf = ['HGNC:12511']
rnap1_tfs = taf + ubf
UTP10 = ['HGNC:25517'] # a' cleavage
RNASEN = ['HGNC:17904'] # assumed site 02 endonucleolytic cleavage
RMRP = ['HGNC:10031']
UTP23 = ['HGNC:28224']
UTP24 = ['HGNC:20220']
PARN = ['HGNC:8609']
PAPD5 = ['HGNC:30758'] #TENT4B
NOB1 = ['HGNC:29540']
LAS1 = ['HGNC:25726']
DIS3 = ['HGNC:20604']
ISG20L2 = ['HGNC:25745']
ERI1 = ['HGNC:23994']


# In[5]:


def get_base_counts_and_elements(seq, triphosphate = True):
    '''Seq is a Bio.Seq object. Triphosphate is a boolean indicated whether the 5 end has a triphosphate. 
    Otherwise assume it is a monophosphate'''
    base_counts = dict()
    for base_letter in seq_element_map.keys():
        base_counts[base_letter] = seq.count(base_letter)
        
    elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0}
    for base_letter in seq_element_map.keys():
        for element in elements.keys():
            elements[element] += base_counts[base_letter]* seq_element_map[base_letter][element]   
    
    #3' end
    elements['H'] += 1 
    elements['O'] += 1
    
    # 5' end
    if triphosphate:
        elements['P'] += 2
        elements['O'] += 6
    else:
        elements['H'] += 1
      
        
    return base_counts, elements

