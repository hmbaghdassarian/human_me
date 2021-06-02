#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils.parameters import human_model


# In[2]:


atp_n = human_model.metabolites.get_by_id('atp_n')
gtp_n = human_model.metabolites.get_by_id('gtp_n')
ppi_n =  human_model.metabolites.get_by_id('ppi_n')
ppi_c =  human_model.metabolites.get_by_id('ppi_c')

# mrna expression
# processing variables------------------------------------------------------------
 
pi_n =  human_model.metabolites.get_by_id('pi_n')
h_n =  human_model.metabolites.get_by_id('h_n')
h2o_n =  human_model.metabolites.get_by_id('h2o_n')
gp = gtp_n.elements
gp['O'] -= 6
gp['P'] -= 2
amet_n =  human_model.metabolites.get_by_id('amet_n')
ahcys_n =  human_model.metabolites.get_by_id('ahcys_n')
adp_n = human_model.metabolites.get_by_id('adp_n')

# mrna degradation------------------------------------------------------------
h_c =  human_model.metabolites.get_by_id('h_c')
h2o_c =  human_model.metabolites.get_by_id('h2o_c')
pi_c = human_model.metabolites.get_by_id('pi_c')
amet_c = human_model.metabolites.get_by_id('amet_c')
ahcys_c = human_model.metabolites.get_by_id('ahcys_c')
amp_c = human_model.metabolites.get_by_id('amp_c')

# rrna expression
atp_c = human_model.metabolites.get_by_id('atp_c')

# protein expression

# nucleus------------------------------------------------------------
gdp_n = human_model.metabolites.get_by_id('gdp_n')

# secretory pathway------------------------------------------------------------
o2_r = human_model.metabolites.get_by_id('o2_r')
h2o2_r = human_model.metabolites.get_by_id('h2o2_r')

o2_c = human_model.metabolites.get_by_id('o2_c')
h2o2_c = human_model.metabolites.get_by_id('h2o2_c')

hdca_r = human_model.metabolites.get_by_id('hdca_r')
gpi_hs_r = human_model.metabolites.get_by_id('gpi_hs_r')
# balanced_gpi = {'N': 6,'O': 42,'S': 1,'P': 5, 'C': 58, 'I': 2,'F': 1,'H': 107}
# m4ataer_0 = human_model.metabolites.get_by_id('gpi_sig_r')
# m4ataer_1 = human_model.metabolites.get_by_id('m_em_3gacpail_hs_r')
# m4ataer_2 = human_model.metabolites.get_by_id('m_em_3gacpail_prot_hs_r')
# m4ataer_3 = human_model.metabolites.get_by_id('pre_prot_r')
# M4ATAer = {m4ataer_0:1,m4ataer_1:-1,m4ataer_2:1, m4ataer_3:-1}
udpacgal_g = human_model.metabolites.get_by_id('udpacgal_g')
udpgal_g = human_model.metabolites.get_by_id('udpgal_g')
uacgam_g = human_model.metabolites.get_by_id('uacgam_g')
h_g = human_model.metabolites.get_by_id('h_g')
udp_g = human_model.metabolites.get_by_id('udp_g')

udpgal_r = human_model.metabolites.get_by_id('udpgal_r')
uacgam_r = human_model.metabolites.get_by_id('uacgam_r')
udpacgal_r = human_model.metabolites.get_by_id('udpacgal_r')
udp_r = human_model.metabolites.get_by_id('udp_r')

hdca_l = human_model.metabolites.get_by_id('hdca_l')
h_l = human_model.metabolites.get_by_id('h_l')
h2o_l = human_model.metabolites.get_by_id('h2o_l')
o2_l = human_model.metabolites.get_by_id('o2_l')
h2o2_l = human_model.metabolites.get_by_id('h2o2_l')

udpacgal_l = human_model.metabolites.get_by_id('udpacgal_l')
udp_l = human_model.metabolites.get_by_id('udp_l')


# In[3]:


# mrna expression
seq_metabolite_map = {human_model.metabolites.get_by_id('utp_n'): 'U' , 
                      gtp_n: 'G',
                      human_model.metabolites.get_by_id('ctp_n'): 'C',
                      atp_n: 'A'}

# RNA backbone elements
seq_element_map = dict()
for k,v in seq_metabolite_map.items():
    elements = k.elements
    elements['O'] = elements['O'] - 7 # lost from incoming ntp
    elements['P'] = elements['P'] - 2 # lost from incoming ntp
    elements['H'] = elements['H'] - 1 # lost from 3' end of growing strand
    seq_element_map[v] = elements

# lariat degradataion------------------------------------------------------------
nmp_map_n = {'C': human_model.metabolites.get_by_id('cmp_n'), 
          'U': human_model.metabolites.get_by_id('ump_n'), 
          'G': human_model.metabolites.get_by_id('gmp_n'), 
          'A': human_model.metabolites.get_by_id('amp_n')}
ntp_map_n = {v: k for k,v in seq_metabolite_map.items()}

# mrna degradation------------------------------------------------------------
nmp_map_c = {'C': human_model.metabolites.get_by_id('cmp_c'), 
          'U': human_model.metabolites.get_by_id('ump_c'), 
          'G': human_model.metabolites.get_by_id('gmp_c'), 
          'A': amp_c}
ndp_map_c = {'C': human_model.metabolites.get_by_id('cdp_c'), 
          'U': human_model.metabolites.get_by_id('udp_c'), 
          'G': human_model.metabolites.get_by_id('gdp_c'), 
          'A': human_model.metabolites.get_by_id('adp_c')}

# rrna expression
ntp_map_c = {'C': human_model.metabolites.get_by_id('ctp_c'), 
             'U': human_model.metabolites.get_by_id('utp_c'), 
             'G': human_model.metabolites.get_by_id('gtp_c'), 
             'A': atp_c}
ntp_map_n = {v:k for k,v in seq_metabolite_map.items()}

# trna expression
seq_amino_acid_map_c = {
    'A': human_model.metabolites.get_by_id('ala_L_c'),
    'R': human_model.metabolites.get_by_id('arg_L_c'),
    'N': human_model.metabolites.get_by_id('asn_L_c'),
    'D': human_model.metabolites.get_by_id('asp_L_c'),
    'C': human_model.metabolites.get_by_id('cys_L_c'),
    'E': human_model.metabolites.get_by_id('glu_L_c'),
    'Q': human_model.metabolites.get_by_id('gln_L_c'),
    'G': human_model.metabolites.get_by_id('gly_c'),
    'H': human_model.metabolites.get_by_id('his_L_c'),
    'I': human_model.metabolites.get_by_id('ile_L_c'),
    'L': human_model.metabolites.get_by_id('leu_L_c'),
    'K': human_model.metabolites.get_by_id('lys_L_c'),
    'M': human_model.metabolites.get_by_id('met_L_c'),
    'F': human_model.metabolites.get_by_id('phe_L_c'),
    'P': human_model.metabolites.get_by_id('pro_L_c'),
    'S': human_model.metabolites.get_by_id('ser_L_c'),
    'T': human_model.metabolites.get_by_id('thr_L_c'),
    'W': human_model.metabolites.get_by_id('trp_L_c'),
    'Y': human_model.metabolites.get_by_id('tyr_L_c'),
    'V': human_model.metabolites.get_by_id('val_L_c'), 
}

seq_amino_acid_map_m = {aa_code: human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_m') for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}
seq_amino_acid_map_l = {aa_code: human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_l') for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}
seq_amino_acid_map_x = {aa_code: human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_x') for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}
seq_amino_acid_map_n = {aa_code: human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_n') for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}
seq_amino_acid_map_r = {aa_code: human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_r') for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}

seq_amino_acid_map_compartments = {'c': seq_amino_acid_map_c, 'x': seq_amino_acid_map_x, 'r': seq_amino_acid_map_r,
                                  'm': seq_amino_acid_map_m, 'n': seq_amino_acid_map_n, 'l': seq_amino_acid_map_l}


# In[4]:



adp_c = ndp_map_c['A']

atp_m = human_model.metabolites.get_by_id('atp_m')
adp_m = human_model.metabolites.get_by_id('adp_m')
h_m = human_model.metabolites.get_by_id('h_m')
pi_m = human_model.metabolites.get_by_id('pi_m')
h2o_m = human_model.metabolites.get_by_id('h2o_m')
h_i = human_model.metabolites.get_by_id('h_i')

h_x = human_model.metabolites.get_by_id('h_x')
h2o_x = human_model.metabolites.get_by_id('h2o_x')
pi_x = human_model.metabolites.get_by_id('pi_x')
atp_x = human_model.metabolites.get_by_id('atp_x')
adp_x = human_model.metabolites.get_by_id('adp_x')

h_r = human_model.metabolites.get_by_id('h_r')
h2o_r = human_model.metabolites.get_by_id('h2o_r')
pi_r = human_model.metabolites.get_by_id('pi_r')
atp_r = human_model.metabolites.get_by_id('atp_r')
adp_r = human_model.metabolites.get_by_id('adp_r')

pi_l = human_model.metabolites.get_by_id('pi_l')
atp_l = human_model.metabolites.get_by_id('atp_l')
adp_l = human_model.metabolites.get_by_id('adp_l')



atp_compartments = {'c': atp_c, 'm': atp_m, 'i': atp_m, 'x': atp_x, 'n': atp_n, 'r': atp_r, 'l': atp_l}
adp_compartments = {'c': adp_c, 'm': adp_m, 'i': adp_m, 'x': adp_x, 'n': adp_n, 'r': adp_r, 'l': adp_l}
h2o_compartments = {'c': h2o_c, 'm': h2o_m, 'i': h2o_m, 'x': h2o_x, 'n': h2o_n, 'r': h2o_r, 'l': h2o_l}
pi_compartments = {'c': pi_c, 'm': pi_m, 'i': pi_m, 'x': pi_x, 'n': pi_n, 'r': pi_r, 'l': pi_l}
h_compartments = {'c': h_c, 'm': h_m, 'i': h_i, 'x': h_x, 'n': h_n, 'r': h_r, 'l': h_l}


# In[5]:


# biomass

# dna
datp_n = human_model.metabolites.get_by_id('datp_n')
dctp_n = human_model.metabolites.get_by_id('dctp_n')
dgtp_n = human_model.metabolites.get_by_id('dgtp_n')
dttp_n = human_model.metabolites.get_by_id('dttp_n')

# carbohydrate
g6p_c = human_model.metabolites.get_by_id('g6p_c')

# lipid
chsterol_c = human_model.metabolites.get_by_id('chsterol_c')
clpn_hs_c = human_model.metabolites.get_by_id('clpn_hs_c')
pail_hs_c = human_model.metabolites.get_by_id('pail_hs_c')
pchol_hs_c = human_model.metabolites.get_by_id('pchol_hs_c')
pe_hs_c = human_model.metabolites.get_by_id('pe_hs_c')
pglyc_hs_c = human_model.metabolites.get_by_id('pglyc_hs_c')
ps_hs_c = human_model.metabolites.get_by_id('ps_hs_c')
sphmyln_hs_c = human_model.metabolites.get_by_id('sphmyln_hs_c')


# In[6]:


# h2o = {'r': h2o_r, 'c': h2o_c, 'l': h2o_l, 'm': h2o_m, 'n': h2o_n, 'x': h2o_x}
nmp_map = {'n': nmp_map_n, 'c': nmp_map_c}
ntp_map = {'n': ntp_map_n, 'c': ntp_map_c}

