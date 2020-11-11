#!/usr/bin/env python
# coding: utf-8

# In[65]:


import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
# from utils.load_environmental_variables import *
from utils.parameters import human_model


# In[66]:


atp_n = human_model.metabolites.get_by_id('atp[n]')
gtp_n = human_model.metabolites.get_by_id('gtp[n]')
ppi_n =  human_model.metabolites.get_by_id('ppi[n]')
ppi_c =  human_model.metabolites.get_by_id('ppi[c]')

# mrna expression
# processing variables------------------------------------------------------------
 
pi_n =  human_model.metabolites.get_by_id('pi[n]')
h_n =  human_model.metabolites.get_by_id('h[n]')
h2o_n =  human_model.metabolites.get_by_id('h2o[n]')
gp = gtp_n.elements
gp['O'] -= 6
gp['P'] -= 2
amet_n =  human_model.metabolites.get_by_id('amet[n]')
ahcys_n =  human_model.metabolites.get_by_id('ahcys[n]')
adp_n = human_model.metabolites.get_by_id('adp[n]')

# mrna degradation------------------------------------------------------------
h_c =  human_model.metabolites.get_by_id('h[c]')
h2o_c =  human_model.metabolites.get_by_id('h2o[c]')
pi_c = human_model.metabolites.get_by_id('pi[c]')
amet_c = human_model.metabolites.get_by_id('amet[c]')
ahcys_c = human_model.metabolites.get_by_id('ahcys[c]')
amp_c = human_model.metabolites.get_by_id('amp[c]')

# rrna expression
atp_c = human_model.metabolites.get_by_id('atp[c]')

# protein expression

# nucleus------------------------------------------------------------
gdp_n = human_model.metabolites.get_by_id('gdp[n]')

# secretory pathway------------------------------------------------------------
o2_r = human_model.metabolites.get_by_id('o2[r]')
h2o2_r = human_model.metabolites.get_by_id('h2o2[r]')
hdca_r = human_model.metabolites.get_by_id('hdca[r]')
gpi_hs_r = human_model.metabolites.get_by_id('gpi_hs[r]')
# balanced_gpi = {'N': 6,'O': 42,'S': 1,'P': 5, 'C': 58, 'I': 2,'F': 1,'H': 107}
# m4ataer_0 = human_model.metabolites.get_by_id('gpi_sig[r]')
# m4ataer_1 = human_model.metabolites.get_by_id('m_em_3gacpail_hs[r]')
# m4ataer_2 = human_model.metabolites.get_by_id('m_em_3gacpail_prot_hs[r]')
# m4ataer_3 = human_model.metabolites.get_by_id('pre_prot[r]')
# M4ATAer = {m4ataer_0:1,m4ataer_1:-1,m4ataer_2:1, m4ataer_3:-1}
udpacgal_g = human_model.metabolites.get_by_id('udpacgal[g]')
udpgal_g = human_model.metabolites.get_by_id('udpgal[g]')
uacgam_g = human_model.metabolites.get_by_id('uacgam[g]')
h_g = human_model.metabolites.get_by_id('h[g]')
udp_g = human_model.metabolites.get_by_id('udp[g]')

udpgal_r = human_model.metabolites.get_by_id('udpgal[r]')
uacgam_r = human_model.metabolites.get_by_id('uacgam[r]')
udpacgal_r = human_model.metabolites.get_by_id('udpacgal[r]')
udp_r = human_model.metabolites.get_by_id('udp[r]')

hdca_l = human_model.metabolites.get_by_id('hdca[l]')
h_l = human_model.metabolites.get_by_id('h[l]')
h2o_l = human_model.metabolites.get_by_id('h2o[l]')
o2_l = human_model.metabolites.get_by_id('o2[l]')
h2o2_l = human_model.metabolites.get_by_id('h2o2[l]')

udpacgal_l = human_model.metabolites.get_by_id('udpacgal[l]')
udp_l = human_model.metabolites.get_by_id('udp[l]')


# In[67]:


# mrna expression
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

# lariat degradataion------------------------------------------------------------
nmp_map_n = {'C': human_model.metabolites.get_by_id('cmp[n]'), 
          'U': human_model.metabolites.get_by_id('ump[n]'), 
          'G': human_model.metabolites.get_by_id('gmp[n]'), 
          'A': human_model.metabolites.get_by_id('amp[n]')}
ntp_map_n = {v: k for k,v in seq_metabolite_map.items()}

# mrna degradation------------------------------------------------------------
nmp_map_c = {'C': human_model.metabolites.get_by_id('cmp[c]'), 
          'U': human_model.metabolites.get_by_id('ump[c]'), 
          'G': human_model.metabolites.get_by_id('gmp[c]'), 
          'A': amp_c}
ndp_map_c = {'C': human_model.metabolites.get_by_id('cdp[c]'), 
          'U': human_model.metabolites.get_by_id('udp[c]'), 
          'G': human_model.metabolites.get_by_id('gdp[c]'), 
          'A': human_model.metabolites.get_by_id('adp[c]')}

# rrna expression
ntp_map_c = {'C': human_model.metabolites.get_by_id('ctp[c]'), 
             'U': human_model.metabolites.get_by_id('utp[c]'), 
             'G': human_model.metabolites.get_by_id('gtp[c]'), 
             'A': atp_c}
ntp_map_n = {v:k for k,v in seq_metabolite_map.items()}

# trna expression
seq_amino_acid_map_c = {
    'A': human_model.metabolites.get_by_id('ala_L[c]'),
    'R': human_model.metabolites.get_by_id('arg_L[c]'),
    'N': human_model.metabolites.get_by_id('asn_L[c]'),
    'D': human_model.metabolites.get_by_id('asp_L[c]'),
    'C': human_model.metabolites.get_by_id('cys_L[c]'),
    'E': human_model.metabolites.get_by_id('glu_L[c]'),
    'Q': human_model.metabolites.get_by_id('gln_L[c]'),
    'G': human_model.metabolites.get_by_id('gly[c]'),
    'H': human_model.metabolites.get_by_id('his_L[c]'),
    'I': human_model.metabolites.get_by_id('ile_L[c]'),
    'L': human_model.metabolites.get_by_id('leu_L[c]'),
    'K': human_model.metabolites.get_by_id('lys_L[c]'),
    'M': human_model.metabolites.get_by_id('met_L[c]'),
    'F': human_model.metabolites.get_by_id('phe_L[c]'),
    'P': human_model.metabolites.get_by_id('pro_L[c]'),
    'S': human_model.metabolites.get_by_id('ser_L[c]'),
    'T': human_model.metabolites.get_by_id('thr_L[c]'),
    'W': human_model.metabolites.get_by_id('trp_L[c]'),
    'Y': human_model.metabolites.get_by_id('tyr_L[c]'),
    'V': human_model.metabolites.get_by_id('val_L[c]'), 
}

seq_amino_acid_map_m = {aa_code: human_model.metabolites.get_by_id(met_obj.id.replace('[c]', '[m]')) for aa_code, met_obj in seq_amino_acid_map_c.items()}
seq_amino_acid_map_l = {aa_code: human_model.metabolites.get_by_id(aa_metabolite.id.replace('[c]', '[l]')) for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}
seq_amino_acid_map_x = {aa_code: human_model.metabolites.get_by_id(met_obj.id.replace('[c]', '[x]')) for aa_code, met_obj in seq_amino_acid_map_c.items()}
seq_amino_acid_map_n = {aa_code: human_model.metabolites.get_by_id(aa_metabolite.id.replace('[c]', '[n]')) for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}
seq_amino_acid_map_r = {aa_code: human_model.metabolites.get_by_id(aa_metabolite.id.replace('[c]', '[r]')) for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}

seq_amino_acid_map_compartments = {'c': seq_amino_acid_map_c, 'x': seq_amino_acid_map_x, 'r': seq_amino_acid_map_r,
                                  'm': seq_amino_acid_map_m, 'n': seq_amino_acid_map_n, 'l': seq_amino_acid_map_l}



# In[74]:



adp_c = ndp_map_c['A']

atp_m = human_model.metabolites.get_by_id('atp[m]')
adp_m = human_model.metabolites.get_by_id('adp[m]')
h_m = human_model.metabolites.get_by_id('h[m]')
pi_m = human_model.metabolites.get_by_id('pi[m]')
h2o_m = human_model.metabolites.get_by_id('h2o[m]')
h_i = human_model.metabolites.get_by_id('h[i]')

h_x = human_model.metabolites.get_by_id('h[x]')
h2o_x = human_model.metabolites.get_by_id('h2o[x]')
pi_x = human_model.metabolites.get_by_id('pi[x]')
atp_x = human_model.metabolites.get_by_id('atp[x]')
adp_x = human_model.metabolites.get_by_id('adp[x]')

h_r = human_model.metabolites.get_by_id('h[r]')
h2o_r = human_model.metabolites.get_by_id('h2o[r]')
pi_r = human_model.metabolites.get_by_id('pi[r]')
atp_r = human_model.metabolites.get_by_id('atp[r]')
adp_r = human_model.metabolites.get_by_id('adp[r]')

pi_l = human_model.metabolites.get_by_id('pi[l]')
atp_l = human_model.metabolites.get_by_id('atp[l]')
adp_l = human_model.metabolites.get_by_id('adp[l]')



atp_compartments = {'c': atp_c, 'm': atp_m, 'i': atp_m, 'x': atp_x, 'n': atp_n, 'r': atp_r, 'l': atp_l}
adp_compartments = {'c': adp_c, 'm': adp_m, 'i': adp_m, 'x': adp_x, 'n': adp_n, 'r': adp_r, 'l': adp_l}
h2o_compartments = {'c': h2o_c, 'm': h2o_m, 'i': h2o_m, 'x': h2o_x, 'n': h2o_n, 'r': h2o_r, 'l': h2o_l}
pi_compartments = {'c': pi_c, 'm': pi_m, 'i': pi_m, 'x': pi_x, 'n': pi_n, 'r': pi_r, 'l': pi_l}
h_compartments = {'c': h_c, 'm': h_m, 'i': h_i, 'x': h_x, 'n': h_n, 'r': h_r, 'l': h_l}


# In[69]:


# biomass

# dna
datp_n = human_model.metabolites.get_by_id('datp[n]')
dctp_n = human_model.metabolites.get_by_id('dctp[n]')
dgtp_n = human_model.metabolites.get_by_id('dgtp[n]')
dttp_n = human_model.metabolites.get_by_id('dttp[n]')

# carbohydrate
g6p_c = human_model.metabolites.get_by_id('g6p[c]')

# lipid
chsterol_c = human_model.metabolites.get_by_id('chsterol[c]')
clpn_hs_c = human_model.metabolites.get_by_id('clpn_hs[c]')
pail_hs_c = human_model.metabolites.get_by_id('pail_hs[c]')
pchol_hs_c = human_model.metabolites.get_by_id('pchol_hs[c]')
pe_hs_c = human_model.metabolites.get_by_id('pe_hs[c]')
pglyc_hs_c = human_model.metabolites.get_by_id('pglyc_hs[c]')
ps_hs_c = human_model.metabolites.get_by_id('ps_hs[c]')
sphmyln_hs_c = human_model.metabolites.get_by_id('sphmyln_hs[c]')


# In[75]:


h2o = {'r': h2o_r, 'c': h2o_c, 'l': h2o_l, 'm': h2o_m, 'n': h2o_n, 'x': h2o_x}
nmp_map = {'n': nmp_map_n, 'c': nmp_map_c}
ntp_map = {'n': ntp_map_n, 'c': ntp_map_c}

