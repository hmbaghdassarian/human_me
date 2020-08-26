#!/usr/bin/env python
# coding: utf-8

# In[19]:


import cobra
import pandas as pd

import sys
sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *


# In[4]:


# universal variables and inputs
psim_me = pd.read_csv(local_data_path + 'processed/corrected_psim_me.csv', index_col = 0) #psim_me = pd.read_csv(root_path+'TRASH.csv', index_col = 0)
psim_me['SP'] = psim_me['SP'].map({1: True, 0: False})
human_model = cobra.io.load_json_model(local_data_path + 'processed/corrected_model.json')
ptt_length = 160 # amino acid length greater than which co-translatioanl translocatoin occurs rather than post-translational
nuclear_diffusion_limit = 40000 # 40 kDA and less proteins diffuse through nucleus
rate_intron = 10/67000 # 10 introns / 67 kbp
L_sp = 22 # secretory pathway signal peptide degradation
Kv = 0.7 # secretory pathway vesicle coat coefficients


# Gene information class variables

# In[48]:


# define necessary variables 
rs = pd.read_csv(local_data_path + 'raw/small_ribosomal_protein.csv', index_col = None, skiprows = [0])
rl = pd.read_csv(local_data_path + 'raw/large_ribosomal_protein.csv', index_col = None, skiprows = [0])


expression_model = cobra.io.json.load_json_model(root_path + 'expression_module_model.json')
# to work with gene_information class
expression_model_2 = expression_model.copy()
for r in expression_model_2.genes.get_by_id('ribosome').reactions:
    r.gene_reaction_rule = r.gene_reaction_rule.replace('ribosome', ' and '.join(rs['HGNC ID (gene)'].tolist() + rl['HGNC ID (gene)'].tolist()))


expression_psim = pd.read_csv(root_path + 'expression_module_psim.csv', index_col = 0)
metabolic_machinery = [g.id for g in human_model.genes] # not efficient to do this each time
expression_machinery = expression_psim.HGNC_ID.tolist()
all_machinery = metabolic_machinery + expression_machinery


compartments = {'c': 'cytosol',  'l': 'lysosome', 'm': 'mitochondria', 'r': 'endoplasmic reticulum', 
                'e': 'extracellular space', 'x': 'peroxisome/glyoxysome', 'n': 'nucleus', 'g': 'golgi apparatus',
                'i': 'inner mitochondrial compartment', 'pm': 'plasma membrane'}
allowed_ptms = {'dsb': 'disulfide bond formation', 'gpi': 'GPI Anchor', 'og': 'O-linked glycosylation'}#,
               #'ng': 'N-linked glycosylation'}

amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

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

# In[49]:


ppi_n =  human_model.metabolites.get_by_id('ppi[n]')


# machinery

# elongation machinery------------------------------------------------------------
rnap = pd.read_csv(local_data_path + 'raw/RNAP_HUGO.csv', index_col = None, skiprows = [0])
rnap2 = rnap[rnap['Approved name'].isin([i for i in rnap['Approved name'] if ' II ' in i])]

tfiis, tfiif, ell = ['HGNC:11612', 'HGNC:11614'], ['HGNC:4652', 'HGNC:4653'], ['HGNC:23114', 'HGNC:17064', 'HGNC:23113']
elongin = pd.read_csv(local_data_path + 'raw/elongin.csv', index_col = None, skiprows = [0])
elongin.drop(index = elongin[elongin['HGNC ID (gene)'] == 'HGNC:24617'].index, inplace = True)
elongator = pd.read_csv(local_data_path + 'raw/elongator.csv', index_col = None, skiprows = [0])
fact = ['HGNC:11327', 'HGNC:11465']
ec = rnap2['HGNC ID (gene)'].tolist() + elongin['HGNC ID (gene)'].tolist() + elongator['HGNC ID (gene)'].tolist()
ec += tfiis + tfiif + ell + fact

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
XRN1 = ['HGNC:30654']
decapping_degradation_machinery = {'LSM1-7 Complex': ['HGNC:20472', 'HGNC:13940', 'HGNC:17874', 'HGNC:17259', 
                                   'HGNC:17162', 'HGNC:17017', 'HGNC:20470'],
                'Decapping': ['HGNC:18714', 'HGNC:24451', 'HGNC:24452', 'HGNC:17157'], 
               "5' Exonuclease": XRN1}
decapping_degradation_machinery = [item for sublist in [v for v in decapping_degradation_machinery.values()] for item in sublist]
degradation_rule1 = ' and '.join(deadenylation_machinery + mrna_degradation_machinery_1)
decapping_rule = ' and '.join(deadenylation_machinery + decapping_degradation_machinery)


# rrna_expression variables

# In[50]:


# variables needed
atp_c = human_model.metabolites.get_by_id('atp[c]')
ntp_map_c = {'C': human_model.metabolites.get_by_id('ctp[c]'), 
             'U': human_model.metabolites.get_by_id('utp[c]'), 
             'G': human_model.metabolites.get_by_id('gtp[c]'), 
             'A': atp_c}
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
RAN = ['HGNC:9846']
XPO1 = ['HGNC:12825']


# trna variables

# In[51]:


# add to utils
allowed_trna_modifications = {}

TRNT1 = ['HGNC:17341']
RNASEP = ['HGNC:30129', 'HGNC:30329', 'HGNC:30081', 'HGNC:17689', 'HGNC:19949', 'HGNC:21300',
         'HGNC:17688', 'HGNC:20992', 'HGNC:30361', 'HGNC:30327']
RNASEZ = ['HGNC:14198']
trna_splicing_machinery = ['HGNC:28422', 'HGNC:16791', 'HGNC:15506', 'HGNC:27561']
XPOT = ['HGNC:12826'] # nuclear export of trna

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
ppi_c =  human_model.metabolites.get_by_id('ppi[c]')

classI_synthetase = pd.read_csv(local_data_path + 'raw/classI_aa_trna_synthetases.csv', index_col = None, 
                                skiprows = [0])
classII_synthetase = pd.read_csv(local_data_path + 'raw/classII_aa_trna_synthetases.csv', index_col = None, 
                                 skiprows = [0])
trna_synthetase = pd.concat([classI_synthetase, classII_synthetase], axis = 0)
trna_synthetase.reset_index(inplace = True, drop = True)
drop_idx = [i for i in trna_synthetase.index if ('mitochondria' in trna_synthetase.loc[i, 'Approved name']) or ('2' in trna_synthetase.loc[i, 'Approved name'])]
trna_synthetase.drop(index = drop_idx, inplace = True)
trna_synthetase.drop_duplicates(keep='first', inplace = True)
trna_synthetase.reset_index(inplace = True, drop = True)

seq_synthetase_map = {
    'A': trna_synthetase[trna_synthetase['Approved name'] == 'alanyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'R': trna_synthetase[trna_synthetase['Approved name'] == 'arginyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'N': trna_synthetase[trna_synthetase['Approved name'] == 'asparaginyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'D': trna_synthetase[trna_synthetase['Approved name'] == 'aspartyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'C': trna_synthetase[trna_synthetase['Approved name'] == 'cysteinyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'E': trna_synthetase[trna_synthetase['Approved name'] == 'glutamyl-prolyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'Q': trna_synthetase[trna_synthetase['Approved name'] == 'glutaminyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'G': trna_synthetase[trna_synthetase['Approved name'] == 'glycyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'H': trna_synthetase[trna_synthetase['Approved name'] == 'histidyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'I': trna_synthetase[trna_synthetase['Approved name'] == 'isoleucyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'L': trna_synthetase[trna_synthetase['Approved name'] == 'leucyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'K': trna_synthetase[trna_synthetase['Approved name'] == 'lysyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'M': trna_synthetase[trna_synthetase['Approved name'] == 'methionyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'F': trna_synthetase[trna_synthetase['Approved name'] == 'phenylalanyl-tRNA synthetase subunit alpha']['HGNC ID (gene)'].tolist() + trna_synthetase[trna_synthetase['Approved name'] == 'phenylalanyl-tRNA synthetase subunit beta']['HGNC ID (gene)'].tolist(),
    'P': trna_synthetase[trna_synthetase['Approved name'] == 'glutamyl-prolyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'S': trna_synthetase[trna_synthetase['Approved name'] == 'seryl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'T': trna_synthetase[trna_synthetase['Approved name'] == 'threonyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'W': trna_synthetase[trna_synthetase['Approved name'] == 'tryptophanyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'Y': trna_synthetase[trna_synthetase['Approved name'] == 'tyrosyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist(),
    'V': trna_synthetase[trna_synthetase['Approved name'] == 'valyl-tRNA synthetase 1']['HGNC ID (gene)'].tolist()
}


# In[52]:


# cytoplasmic transport expression
transport_translocation_atp_cost = 0.5 # 1 ATP/2 residues
proteolysis_translocation_atp_cost = 0.5 # 1 ATP/2 residues

translation_efs = ['HGNC:3189', 'HGNC:3214', 'HGNC:3208', 'HGNC:3300']
n_ub = 4 # see this - no. of ubiquitins to add to protein

# # DON'T DELETE------------------------------------------------
# import pandas as pd
# e3_ligase = pd.read_csv(local_data_path + 'raw/E3_HPA.tsv', sep = '\t')
# e3_ligase = e3_ligase[e3_ligase['RNA cell line specificity'] == 'Low cell line specificity']
# e3_ligase = e3_ligase.loc[e3_ligase[e3_ligase['Subcellular main location'].notna()].index,:]
# e3_ligase = e3_ligase.loc[[i for i in e3_ligase.index if ('Cytosol' in e3_ligase.loc[i, 'Subcellular main location'])], :]
# # e3_ligase = e3_ligase[e3_ligase['Subcellular main location'] == 'Cytosol']
# top_gene_idx = e3_ligase.iloc[:, e3_ligase.columns.tolist().index('Tissue RNA - adipose tissue [NX]'):].mean(axis = 1).sort_values(ascending = False).index.tolist()[3]
# e3_uniprot_id = e3_ligase.loc[top_gene_idx, 'Uniprot']

# e2_ligase = pd.read_csv(local_data_path + 'raw/E2_HPA.tsv', sep = '\t')
# e2_ligase = e2_ligase[e2_ligase['RNA cell line specificity'] == 'Low cell line specificity']
# e2_ligase = e2_ligase.loc[e2_ligase[e2_ligase['Subcellular main location'].notna()].index,:]
# e2_ligase = e2_ligase.loc[[i for i in e2_ligase.index if ('Cytosol' in e2_ligase.loc[i, 'Subcellular main location'])], :]
# # e2_ligase = e2_ligase[e2_ligase['Subcellular main location'] == 'Cytosol']
# top_gene_idx = e2_ligase.iloc[:, e2_ligase.columns.tolist().index('Tissue RNA - adipose tissue [NX]'):].mean(axis = 1).sort_values(ascending = False).index.tolist()[0]
# e2_uniprot_id = e2_ligase.loc[top_gene_idx, 'Uniprot']

# e3_ligase = pd.read_csv(local_data_path + 'raw/E3_HPA.tsv', sep = '\t')
# e3_ligase = e3_ligase[e3_ligase['RNA cell line specificity'] == 'Low cell line specificity']
# e3_ligase = e3_ligase.loc[e3_ligase[e3_ligase['Subcellular main location'].notna()].index,:]
# e3_ligase = e3_ligase.loc[[i for i in e3_ligase.index if ('Nucleoplasm' in e3_ligase.loc[i, 'Subcellular main location'])], :]
# # e3_ligase = e3_ligase[e3_ligase['Subcellular main location'] == 'Cytosol']
# top_gene_idx = e3_ligase.iloc[:, e3_ligase.columns.tolist().index('Tissue RNA - adipose tissue [NX]'):].mean(axis = 1).sort_values(ascending = False).index.tolist()[0]
# e3_uniprot_id = e3_ligase.loc[top_gene_idx, 'Uniprot']
# e2_ligase = pd.read_csv(local_data_path + 'raw/E2_HPA.tsv', sep = '\t')
# e2_ligase = e2_ligase[e2_ligase['RNA cell line specificity'] == 'Low cell line specificity']
# e2_ligase = e2_ligase.loc[e2_ligase[e2_ligase['Subcellular main location'].notna()].index,:]
# e2_ligase = e2_ligase.loc[[i for i in e2_ligase.index if ('Nucleoplasm' in e2_ligase.loc[i, 'Subcellular main location'])], :]
# # e2_ligase = e2_ligase[e2_ligase['Subcellular main location'] == 'Cytosol']
# top_gene_idx = e2_ligase.iloc[:, e2_ligase.columns.tolist().index('Tissue RNA - adipose tissue [NX]'):].mean(axis = 1).sort_values(ascending = False).index.tolist()[0]
# e2_uniprot_id = e2_ligase.loc[top_gene_idx, 'Uniprot']
# # DON'T DELETE------------------------------------------------

USP5, UBA1, UBE2D3, STUB1 = ['HGNC:12628'], ['HGNC:12469'], ['HGNC:12476'], ['HGNC:11427']
RNF181, UB2EV1 = ['HGNC:28037'], ['HGNC:12494']
UB_ligases_c = UBA1 + UBE2D3 + STUB1
UB_ligases_n = UBA1 + UB2EV1 + RNF181

proteasome_structural = ['HGNC:9554', 'HGNC:9560', 'HGNC:9557', 'HGNC:9556', 'HGNC:9564', 'HGNC:9565', 'HGNC:9558',
                        'HGNC:9566', 'HGNC:9567']
proteasome_ubiquitin = ['HGNC:9559', 'HGNC:15759', 'HGNC:16889', 'HGNC:9561', 'HGNC:12612', 'HGNC:19678']
proteasome_atpase = ['HGNC:9548', 'HGNC:9547', 'HGNC:9551', 'HGNC:9553', 'HGNC:9549', 'HGNC:9552']
proteasome_machinery = proteasome_structural + proteasome_ubiquitin + proteasome_atpase

seq_amino_acid_map_m = {aa_code: human_model.metabolites.get_by_id(met_obj.id.replace('[c]', '[m]')) for aa_code, met_obj in seq_amino_acid_map_c.items()}


# mitochondria

TOM = ['HGNC:31369', 'HGNC:34528', 'HGNC:21648', 'HGNC:20947', 'HGNC:18002', 'HGNC:18001', 'HGNC:11985']
# HGNC's tim23 already contains PAM
TIM23_PAM = pd.read_csv(local_data_path + 'raw/tim23_complex.csv',  index_col = None, skiprows = [0])['HGNC ID (gene)'].tolist()
HSP70_m = ['HGNC:5244'] # mitocondrial version
OXA = ['HGNC:8526'] # inner membrane transport

mLON, iAAA   = ['HGNC:9479'], ['HGNC:12843']# mitochondrial proteases
#mAAA =  ['HGNC:315', 'HGNC:11237'] 
HSP70_c, HSP40_c  = ['HGNC:5233'], ['HGNC:5229']

# peroxisome
seq_amino_acid_map_x = {aa_code: human_model.metabolites.get_by_id(met_obj.id.replace('[c]', '[x]')) for aa_code, met_obj in seq_amino_acid_map_c.items()}
PEX5, L_PEX5 = ['HGNC:9719'], 639 # Uniprot and PSIM_ME agree on this number
peroxins = ['HGNC:22965', 'HGNC:8859', 'HGNC:8850', 'HGNC:8856', 'HGNC:8855'] + PEX5
AWP1 = ['HGNC:30164']
LONP2 = ['HGNC:20598']



# nucleus
seq_amino_acid_map_n = {aa_code: human_model.metabolites.get_by_id(aa_metabolite.id.replace('[c]', '[n]')) for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}


gdp_n = human_model.metabolites.get_by_id('gdp[n]')
importins = ['HGNC:6400', 'HGNC:6394']


# In[53]:


# ribosome biogenesis
UCHL3 = ['HGNC:12515']
pre40s_rbfs = ['HGNC:25542', 'HGNC:21173', 'HGNC:32790', 'HGNC:29100']
pre60s_rbfs = ['HGNC:18477', 'HGNC:25789', 'HGNC:19440', 'HGNC:20870', 'HGNC:17083', 'HGNC:4333']

eif1, eif2 = ['HGNC:3249', 'HGNC:3250'], ['HGNC:3265', 'HGNC:3266', 'HGNC:3267']
eif3 = pd.read_csv(local_data_path + 'raw/eif3.csv', index_col = None, skiprows = [0])['HGNC ID (gene)'].tolist()
eif4f = ['HGNC:3282', 'HGNC:3284', 'HGNC:3287', 'HGNC:3296']
eif5 = ['HGNC:3299', 'HGNC:30793']
eifs = eif1 + eif2 + eif3 + eif4f + eif5 + ['HGNC:8554']


# In[54]:


# secretory pathway

ASNA1, WRB  = ['HGNC:752'], ['HGNC:12790']
ptnm = ['HGNC:20090', 'HGNC:5238', 'HGNC:7670', 'HGNC:18276', 'HGNC:16993', 
                            'HGNC:18277', 'HGNC:11846', 'HGNC:21082', 'HGNC:10759', 'HGNC:11323', 'HGNC:11324', 
                            'HGNC:11325', 'HGNC:11326'] + ASNA1
ctnm = ['HGNC:20090', 'HGNC:5238', 'HGNC:16931', 'HGNC:7670', 'HGNC:10448', 'HGNC:17718', 'HGNC:23400', 
        'HGNC:18276', 'HGNC:16993', 'HGNC:18277', 'HGNC:11846', 'HGNC:21082', 'HGNC:10759', 'HGNC:24624', 
        'HGNC:23401', 'HGNC:28962', 'HGNC:26212', 'HGNC:11299', 'HGNC:11300', 'HGNC:11301', 'HGNC:11302', 
        'HGNC:11303', 'HGNC:11307', 'HGNC:11323', 'HGNC:11324', 'HGNC:11325', 'HGNC:11326', 'HGNC:11740']
sp_map = {'5682': 'HGNC:9530', '5683': 'HGNC:9531', '5684': 'HGNC:9532', '5685': 'HGNC:9533', '5686': 'HGNC:9534', 
          '5688': 'HGNC:9536', '143471': 'HGNC:22985', '5689': 'HGNC:9537', '5690': 'HGNC:9539', 
          '5691': 'HGNC:9540', '5692': 'HGNC:9541', '5693': 'HGNC:9542', '5694': 'HGNC:9543', '5695': 'HGNC:9544', 
          '5696': 'HGNC:9545', '5698': 'HGNC:9546', '5699': 'HGNC:9538', '122706': 'HGNC:31963', 
          '5707': 'HGNC:9554', '5708': 'HGNC:9559', '5710': 'HGNC:9561', '5711': 'HGNC:9563', '9861': 'HGNC:9564', 
          '5713': 'HGNC:9565', '5714': 'HGNC:9566', '5715': 'HGNC:9567', '5716': 'HGNC:9555', '5717': 'HGNC:9556', 
          '5718': 'HGNC:9557', '5719': 'HGNC:9558', '10213': 'HGNC:16889', '8624': 'HGNC:3043', 
          '56984': 'HGNC:24929', '84262': 'HGNC:22420', '389362': 'HGNC:21108'}
sp_rule = '((5682) or (5683) or (5684) or (5685) or (5686) or (5688) or (143471)) and ((5689) or (5690) or (5691) or (5692) or (5693) or (5694) or (5695) or (5696) or (5698) or (5699) or (122706)) and (8624) and (56984) and (84262) and (389362) and (5707) and (5708) and (5711) and (5710) and (9861) and (5713) and (5714) and (5715) and (5716) and (5717) and (5718) and (5719) and (10213)'
for k,v in sp_map.items():
    if k != '5698':
        sp_rule = sp_rule.replace(k,v)
sp_rule = sp_rule.replace('5698', 'HGNC:9546')
    

o2_r = human_model.metabolites.get_by_id('o2[r]')
h2o2_r = human_model.metabolites.get_by_id('h2o2[r]')

P4HB = ['HGNC:8548']
gpi_machinery = ['HGNC:4446', 'HGNC:25712', 'HGNC:8965', 'HGNC:14937', 'HGNC:14938', 'HGNC:15791']
hdca_r = human_model.metabolites.get_by_id('hdca[r]')
gpi_hs_r = human_model.metabolites.get_by_id('gpi_hs[r]')
balanced_gpi = {'N': 6,'O': 42,'S': 1,'P': 5,'E': 1,'C': 58,'I': 2,'F': 1,'H': 107, 'R': 1}
# M4ATAer = human_model.reactions.get_by_id('M4ATAer').metabolites
m4ataer_0 = human_model.metabolites.get_by_id('gpi_sig[r]')
m4ataer_1 = human_model.metabolites.get_by_id('m_em_3gacpail_hs[r]')
m4ataer_2 = human_model.metabolites.get_by_id('m_em_3gacpail_prot_hs[r]')
m4ataer_3 = human_model.metabolites.get_by_id('pre_prot[r]')
M4ATAer = {m4ataer_0:1,m4ataer_1:-1,m4ataer_2:1, m4ataer_3:-1}
#number_BiP = len(gene_info.protein_seq)/40


# In[55]:



copii_r_m = ['HGNC:14562', 'HGNC:4430', 'HGNC:6632', 'HGNC:9758', 'HGNC:10535', 'HGNC:10697', 'HGNC:29006', 
             'HGNC:10700', 'HGNC:10701', 'HGNC:10703', 'HGNC:17052', 'HGNC:11440']
copii_gpi_m = ['HGNC:14562', 'HGNC:4430', 'HGNC:9758', 'HGNC:10535', 'HGNC:10697', 'HGNC:29006', 'HGNC:10700', 
               'HGNC:10701', 'HGNC:10703', 'HGNC:17052', 'HGNC:11440']
udpacgal_g = human_model.metabolites.get_by_id('udpacgal[g]')
udpgal_g = human_model.metabolites.get_by_id('udpgal[g]')
uacgam_g = human_model.metabolites.get_by_id('uacgam[g]')
h_g = human_model.metabolites.get_by_id('h[g]')
udp_g = human_model.metabolites.get_by_id('udp[g]')

og_rule = '(HGNC:16347 or HGNC:19873 or HGNC:4124 or HGNC:4127 or HGNC:4131 or HGNC:4125 or HGNC:4129 or HGNC:19875 or HGNC:4128 or HGNC:23242 or HGNC:4123 or HGNC:4130 or HGNC:4126) and HGNC:24337 and HGNC:24338 and HGNC:4205'
copi_m = ['HGNC:649', 'HGNC:14562', 'HGNC:2230', 'HGNC:2231', 'HGNC:2232', 'HGNC:2234', 'HGNC:2236', 'HGNC:2243', 'HGNC:19356', 
          'HGNC:9758', 'HGNC:10700', 'HGNC:11443', 'HGNC:15942', 'HGNC:25847']
clathrin_m = ['HGNC:652', 'HGNC:2090', 'HGNC:2091', 'HGNC:2092', 'HGNC:17842', 'HGNC:16064', 'HGNC:17079', 
              'HGNC:14902', 'HGNC:11441', 'HGNC:11442', 'HGNC:11430']


# In[56]:


retro_mach_glyco = ['HGNC:16695', 'HGNC:28454', 'HGNC:14236', 'HGNC:18261', 'HGNC:10717', 'HGNC:30396', 'HGNC:20738', 
           'HGNC:12520', 'HGNC:12666']
ERDJ5 = ['HGNC:24637']

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

HSP90AB1 = ['HGNC:5258']
escrt = pd.read_csv(local_data_path + 'raw/escrt_complexes.txt', sep = '\t')['HGNC ID'].tolist()
eps = ['HGNC:3419', 'HGNC:21604']
endocytic_machinery = sorted(set(proteasome_ubiquitin + escrt + eps + clathrin_m))

seq_amino_acid_map_l = {aa_code: human_model.metabolites.get_by_id(aa_metabolite.id.replace('[c]', '[l]')) for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}
cathepsins = ['HGNC:2527', 'HGNC:2529', 'HGNC:9251']


# # Functions

# In[ ]:


def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__


# In[ ]:


def get_reaction_compartment(reaction):
    '''Input is a cobra.Reaction, output is a singular compartment. This function maps reactions to a particular 
    compartment according to some rules'''
    
    compartments_ = list(reaction.compartments.copy())
    if len(compartments_) > 1: # for reactions that occur in more than one compartment
        if 'c' in compartments_ and len(compartments_) == 2: # remove cytoplasmic compartment between the two for machinery
            compartments_.remove('c')
        else: # choose most common compartment 
            compartments_ = [max(compartments_, key = compartments_.count)]
    if len(compartments_) > 1:
        raise ValueError('Failed to map reaction to a singular compartment')
    elif compartments_[0] not in compartments.keys():
        raise ValueError('Mapped reaction to a compartment that is not allowed in ME model')
    else:
        return compartments_[0]


# In[ ]:


seq_amino_acid_map_compartments = {'c': seq_amino_acid_map_c, 'x': seq_amino_acid_map_x, 
                                  'm': seq_amino_acid_map_m, 'n': seq_amino_acid_map_n, 'l': seq_amino_acid_map_l}
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

def hydrolyze_atp(rxn, n_atp, compartment):
    '''
    Rxn is a dict for the cobra.Reaction.add_metabolite function.
    n_atp is the # of atp to hydrolyze
    compartment is the compartment for hydrolysis
    
    '''
    n_atp = round(n_atp)
    
    if atp_compartments[compartment] in rxn.keys():
        rxn[atp_compartments[compartment]] -= n_atp 
    else:
        rxn[atp_compartments[compartment]] = -n_atp 

    if h2o_compartments[compartment] in rxn.keys():
        rxn[h2o_compartments[compartment]] -= n_atp 
    else:
        rxn[h2o_compartments[compartment]] = -n_atp 

    if adp_compartments[compartment] in rxn.keys():
        rxn[adp_compartments[compartment]] += n_atp 
    else:
        rxn[adp_compartments[compartment]] = n_atp

    if pi_compartments[compartment] in rxn.keys():
        rxn[pi_compartments[compartment]] += n_atp 
    else:
        rxn[pi_compartments[compartment]] = n_atp

    if h_compartments[compartment] in rxn.keys():
        rxn[h_compartments[compartment]] += n_atp 
    else:
        rxn[h_compartments[compartment]] = n_atp
    
    return rxn


# In[5]:


def get_base_counts_and_elements(seq, triphosphate = True):
    '''
    
    Inputs:
    1) Seq is a Bio.Seq object or a string representing an RNA sequence. 
    2) triphosphate is a boolean. If true, assume the 5' end has a triphosphate, otherwise assume monophosphate. 
   
   Outputs:
    1) base_counts is a dictionary with keys as one-letter RNA base characters and values as the # of 
    occurences of that base in the RNA sequence
    2) elements is a dictionary emulating cobra.Metabolite.elements
   
   '''
    base_counts = dict()
    for base_letter in seq_element_map.keys():
        base_counts[base_letter] = seq.count(base_letter)
        
    elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0}
    for base_letter in seq_element_map.keys():
        for element in elements.keys():
            elements[element] += base_counts[base_letter]* seq_element_map[base_letter][element]   
    
    #3' OH end
    elements['H'] += 1 
    elements['O'] += 1
    
    # 5' end
    if triphosphate:
        elements['P'] += 2
        elements['O'] += 6
    else:
        elements['H'] += 1
      
        
    return base_counts, elements


# In[ ]:


def make_rna_metabolite(metabolite_name, seq, molecule_type, compartment = 'n', triphosphate = True):
    
    '''
    Inputs:
    1) metabolite_name is the name of the RNA molecule (unique ID)
    2) seq is a string representing the one-letter sequence of the RNA molecule.
    3) molecule type = ['mrna', 'trna', 'rrna']
    4) compartment is the one-letter string representing the location of the RNA molecule (usually 'n' or 'c')
    5) triphosphate is a boolean. If true, assume the 5' end has a triphosphate, otherwise assume monophosphate.
    
    Outputs:
    1) rna_metabolite is an object of cobra.Metabolite representing teh RNA molecule
    2) base_counts is a dictionary with keys as one-letter RNA base characters and values as the # of 
    occurences of that base in the RNA sequence 
    
    '''
    
    if molecule_type not in ['mrna', 'trna', 'rrna']:
        raise ValueError('molecule_type must be mrna, trna, or rrna')

    rna_metabolite = cobra.Metabolite(metabolite_name + '_' + molecule_type + '[' + compartment + ']')
    rna_metabolite.compartment = compartment
    base_counts, elements = get_base_counts_and_elements(seq, triphosphate = triphosphate) # utils function

    rna_metabolite.elements = elements
    rna_metabolite.charge = -len(seq)
    
    if triphosphate:
        rna_metabolite.charge -= 3
    
    return rna_metabolite, base_counts


# In[ ]:


def rna_exonucleolytic_degradation(rna_metabolite, rna_base_counts, rna_sequence, reaction_name, 
                                   triphosphate = True, nucleus = True):
    ''' 
    
    Generates a reaction for exonucleolytic cleavage of an RNA molecule (RNA-->NMPs).
    Inputs:
    1) rna_metabolite is a cobra.Metabolite object representing the rna molecule to be degraded.
    2) rna_base_counts is a dictionary with keys as one-letter RNA base characters and values as the # of 
    occurences of that base in the RNA sequence.
    3) rna_sequence is the ordered (5'-->3') sequence of the RNA molecule, as a string of one-letter bases
    4) reaction_name is a string representing the name you want to give the reaction
    5) triphosphate is a boolean. If true, assume the 5' end has a triphosphate, otherwise assume monophosphate. 
    6) nucleus is a boolean. If true, the degradation reaction is taking place in the nucleus. Otherwise, assume
    it takes place in the cytoplasm (current iteration of model only degrades RNA in these two compartments).
    
    Output: a degradation reaction of type cobra.Reaction
    no GPRs or subsystems added to reaction
    
    '''
    # exonucleolytic cleavage of RNA reaction

    

    if nucleus: 
        rna_degradation = cobra.Reaction(reaction_name + '_DEGRADATIONn')
        rxn = dict()
        rxn[h2o_n] = -sum(rna_base_counts.values())+1
        rxn[rna_metabolite] = -1
        for k,v in nmp_map_n.items():
            rxn[v] = rna_base_counts[k]

        # triphosphate on 5' end
        if triphosphate:
            rxn[nmp_map_n[rna_sequence[0]]] -= 1
            rxn[ntp_map_n[rna_sequence[0]]]  = 1  
            rxn[h_n] = sum(rna_base_counts.values())-1
        else:
            rxn[h_n] = sum(rna_base_counts.values()) # extra H on 5' end <--unsure about this

        rna_degradation.add_metabolites(rxn)

        
    else:
        rna_degradation = cobra.Reaction(reaction_name + '_DEGRADATIONc')
        rxn = dict()
        rxn[h2o_c] = -sum(rna_base_counts.values())+1
        rxn[rna_metabolite] = -1
        for k,v in nmp_map_c.items():
            rxn[v] = rna_base_counts[k]

        # triphosphate on 5' end
        if triphosphate:
            rxn[nmp_map_c[rna_sequence[0]]] -= 1
            rxn[ntp_map_c[rna_sequence[0]]]  = 1  
            rxn[h_c] = sum(rna_base_counts.values())-1
        else:
            rxn[h_c] = sum(rna_base_counts.values()) # extra H on 5' end <--unsure about this

        rna_degradation.add_metabolites(rxn)
        
    return rna_degradation 


# In[ ]:


seq_amino_acid_map_r = {aa_code: human_model.metabolites.get_by_id(aa_metabolite.id.replace('[c]', '[r]')) for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}

seq_amino_acid_map_compartments = {'c': seq_amino_acid_map_c, 'x': seq_amino_acid_map_x, 'r': seq_amino_acid_map_r,
                                  'm': seq_amino_acid_map_m, 'n': seq_amino_acid_map_n}
def make_protein_metabolite(id_, amino_acid_counts, L_protein, compartment):
    '''
    ID is a string to name the protein metabolite. 
    Amino acid counts is a dictionary with keys as the aa one letter code and counts as the number of occurences of that amino acid in the protein sequence
    L_protein is the length of the protein
    Compartment is the location of the protein (one letter string, corresponds to Recon2.2s compartments)
    
    Will return a cobra.Metabolite object with relevant charge and elements.
    
    '''
    
    protein_metabolite = cobra.Metabolite(id_ + '_protein[' + compartment + ']')
    protein_metabolite.compartment = compartment
    elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}
    if compartment in seq_amino_acid_map_compartments.keys():
        for aa_code, aa_count in amino_acid_counts.items():
            aa_elements = seq_amino_acid_map_compartments[compartment][aa_code].elements
            for element in aa_elements:
                elements[element] += aa_count*aa_elements[element]
    else:
        raise ValueError('Must add this compartment to make_protein_metabolite function')

    # peptide bond formation
    elements['H'] -= 2*(L_protein-1)
    elements['O'] -= 1*(L_protein-1)

    protein_metabolite.elements = elements
    # assume charge of amino acid is the ssame regardless of metabolite
    protein_metabolite.charge = sum([seq_amino_acid_map_compartments[compartment][aa_code].charge*aa_count for aa_code, aa_count in amino_acid_counts.items()])
    return protein_metabolite


# In[1]:


def make_complex_metabolite(complex_id = None, **complex_info):# metabolites, *ids, *metabolite_types):
    '''
    
    Inputs:
    Complex info is a dictionary with three keys ['METABOLITES', 'IDS', 'METABOLITE_TYPES']
    Each value is a list:
        Metabolites is a list of cobra.Metabolite objects
        IDs is a list of string identifiers corresponding to each metabolite object
        Metabolite_types is a list of strings; possible values are ['protein', 'rrna', 'trna', 'mrna',  'metabolite']
        This means complexes can form between any of these species, including other complexes; metabolite is a M-model metabolite
    complex_id is a string for the id of the complex metabolite, otherwise will form one from metabolite ids
    Output:
    A cobra.Metabolite object representing the complex formed between metabolites
    
    '''
    if sorted(set(complex_info.keys())) != ['IDS', 'METABOLITES', 'METABOLITE_TYPES']:
        raise ValueError('Invalid complex information keys or insufficient complex information keys')
    
    metabolites, ids, metabolite_types = complex_info['METABOLITES'], complex_info['IDS'], complex_info['METABOLITE_TYPES']
    
    if len(set(metabolite_types).difference(['protein', 'rrna', 'trna', 'mrna',  'metabolite', 'complex']))>1:
        raise ValueError('At least one of the metabolite types is not considered in complex formation currently')
    
    
    compartments = list(set([m.compartment for m in metabolites]))
    if len(compartments) == 1:
        compartment = compartments[0]
    # exception of ribosome complex
    elif (len(compartments) == 2) and ('c' in compartments) and ('mature_ribosome_complex_complex[c]' in [m.id for m in metabolites]):
        compartment = 'c'
    else:
        raise ValueError('Metabolites are not in the same compartment')
    
    mt_type = '_'.join(list(set(metabolite_types)))
    
    ids_ = '_'.join(ids)
    
    if complex_id == None:
        id_ = ids_ + '_' + mt_type
    else: 
        id_ = complex_id + '_' + mt_type
    
    complex_id = id_ + '_complex' + '[' + compartment + ']'
    if len(complex_id)>(256-8-4-len(mt_type)): #-8 and -4 for _complex and compartment appended to end
        err_msg = 'Cobrapy requires metabolite ids to be less than 256 characters, please specify a '
        err_msg += 'shorter user-defined complex id'
        raise ValueError(err_msg)
    
        
    complex_metabolite = cobra.Metabolite(complex_id)
    complex_metabolite.compartment = compartment
    complex_metabolite.charge = sum([m.charge for m in metabolites])
    
    elements = dict()
    for m in metabolites:
        for k,v in m.elements.items():
            if k in elements.keys():
                elements[k] += v
            else:
                elements[k] = v
    complex_metabolite.elements = elements
    
    return complex_metabolite, id_

def form_complex(reaction_id = None, complex_id = None, **complex_info):
    
    '''
    
    Inputs:
    Complex info is a dictionary with three keys ['METABOLITES', 'IDS', 'METABOLITE_TYPES']
    Each value is a list:
        Metabolites is a list of cobra.Metabolite objects
        IDs is a list of string identifiers corresponding to each metabolite object
        Metabolite_types is a list of strings; possible values are ['protein', 'rrna', 'trna', 'mrna',  'metabolite', 'complex']
  
    Output:
    A cobra.Reaction object representing the complex formation between metabolites
    
    '''
    
    complex_metabolite, id_ = make_complex_metabolite(complex_id, **complex_info)
    metabolites = complex_info['METABOLITES']
    compartment = list(set([m.compartment for m in metabolites]))[0]

    if reaction_id == None:
        reaction_id = id_ + '_COMPLEX_FORMATION' + compartment
    else:
        reaction_id = reaction_id + '_COMPLEX_FORMATION' + compartment
    complex_formation = cobra.Reaction(reaction_id)
    
    rxn = {m: -1 for m in metabolites}
    rxn[complex_metabolite] = 1
    complex_formation.add_metabolites(rxn)
    complex_formation.lower_bound = -1000
    
    return complex_formation, complex_metabolite

