#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
# import escher
# from escher import Builder
import pandas as pd 

from Bio.Seq import Seq
from Bio.Alphabet import generic_rna

import numpy as np
import statsmodels.api as sm
import scipy.stats as st

import sys
sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *
from gene_information import human_model, atp_n, gtp_n, seq_metabolite_map, seq_element_map#, allowed_ptms
from gene_information import gene_information


# # start of actual script

# In[2]:


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


mrna_degradation_machinery_1 = {"Exosome": exosome.loc[:, 'HGNC ID (gene)'].tolist() + ['HGNC:29911'], 
               'Cap_Degradation': ['HGNC:29812']}
mrna_degradation_machinery_1 = [item for sublist in [v for v in mrna_degradation_machinery_1.values()] for item in sublist]
decapping_degradation_machinery = {'LSM1-7 Complex': ['HGNC:20472', 'HGNC:13940', 'HGNC:17874', 'HGNC:17259', 
                                   'HGNC:17162', 'HGNC:17017', 'HGNC:20470'],
                'Decapping': ['HGNC:18714', 'HGNC:24451', 'HGNC:24452', 'HGNC:17157'], 
               "5' Exonuclease": ['HGNC:30654']}
decapping_degradation_machinery = [item for sublist in [v for v in decapping_degradation_machinery.values()] for item in sublist]
degradation_rule1 = ' and '.join(deadenylation_machinery + mrna_degradation_machinery_1)
decapping_rule = ' and '.join(deadenylation_machinery + decapping_degradation_machinery)


# In[3]:


class Transcript():
    def __init__(self, gene_information):
        '''Input is an object of the gene_information class, output is all the information needed to 
        build the transcription reations.'''
        self.premrna_seq = Seq(gene_information.premrna_seq, generic_rna)
        self.mrna_seq = Seq(gene_information.mrna_seq, generic_rna) 
        self.id = gene_information.hgnc_id
        self.premrna_base_counts, self.mrna_base_counts = gene_information.premrna_base_counts, gene_information.mrna_base_counts
        
        self.polyA_length = calculate_polyA_length(gene_information.polyA_length)
        self.mrna_seq_length = len(self.mrna_seq)
        self.premrna_seq_length = len(self.premrna_seq)
        
        # metabolite output of transcriptional elongation and processing reactions----------------------
        # combined to save on compute time/for loops
        self.premrna_transcript, self.mrna_transcript_n = cobra.Metabolite(self.id + '_premrna_transcript[n]'), cobra.Metabolite(self.id + '_mrna_transcript[n]')
        self.premrna_transcript.compartment, self.mrna_transcript_n.compartment = 'n', 'n'
        
        elongated_elements, processed_elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0}, {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0}
        for base_letter in seq_element_map.keys():
            for element in elongated_elements.keys():
                elongated_elements[element] += self.premrna_base_counts[base_letter]* seq_element_map[base_letter][element]
                processed_elements[element] += (self.mrna_base_counts[base_letter]* seq_element_map[base_letter][element]) 
        
        #3 and 5' ends
        for dict_ in [elongated_elements, processed_elements]:
            dict_['P'] += 2
            dict_['O'] += 7
            dict_['H'] += 1
        
        self.premrna_transcript.elements = elongated_elements
        self.premrna_transcript.charge = -self.premrna_seq_length - 3 # -3 for 5' end triphosphate
        
        ### processed specific
        
        
        # lariats
        if self.premrna_seq_length > self.mrna_seq_length:
            self.lariats = cobra.Metabolite(self.id + '_lariats[n]')
            self.lariats.compartment = 'n'
            self.lariat_base_counts = {k: v - self.mrna_base_counts[k] for k,v in self.premrna_base_counts.items()}
            
            lariat_elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0}
            for base_letter in seq_element_map.keys():
                for element in lariat_elements.keys():
                    lariat_elements[element] += self.lariat_base_counts[base_letter]* seq_element_map[base_letter][element]
         
            self.lariats.elements = lariat_elements
            self.lariat_length = sum(self.lariat_base_counts.values())
            self.lariats.charge = -self.lariat_length
            
            if gene_information.n_introns == None or pd.isna(gene_information.n_introns):
                self.n_lariats = round(self.premrna_seq_length * rate_intron)
                if self.n_lariats < 1: # atleast one intron must be generated
                    self.n_lariats = 1
    #             n_splices = n_lariats * 2 # because on average, one more exon than intron
            else:
                self.n_lariats = gene_information.n_introns
        
        else:
            self.lariats = None
        
        #mrna
        
        for element in processed_elements.keys():
            processed_elements[element] += (self.polyA_length* seq_element_map['A'][element]) # polyA tail
            processed_elements[element] += gp[element] #5' cap rxn2 - addition of Gp
        
        # 5' cap 
        processed_elements['P'] -= 1 # rxn 1: lost of third triphosphate by RTPase
        processed_elements['O'] -= 4 # rxn 1: loss of third triophosphate by RTPase
        processed_elements['C'] += 2 # rxn 3-4: methyltransferase - cap0 and cap1 structure
        processed_elements['H'] += 5 # methyltransferase - cap0 and cap1 structure
    
        
        self.mrna_transcript_n.elements = processed_elements
        
        # -3 as for premrna_transcript, +2 for cap
        self.mrna_transcript_n.charge = -self.mrna_seq_length - 3 - self.polyA_length + 2 

        
    def build_transcript_elongation_reaction(self):
        '''Input is an object of the Transcript class. Output is reaction (cobra.Reaction object) for
        transcriptional elongation of that gene.'''
        
        # elongation reaction
        # https://www.google.com/search?q=rna+polymerization+reaction&source=lnms&tbm=isch&sa=X&ved=2ahUKEwiN_73Vk7rqAhXOsJ4KHW5lB4UQ_AUoAXoECA4QAw&biw=1920&bih=1001#imgrc=w7XH4mHmJglCuM
        self.transcript_elongation = cobra.Reaction(self.id + '_TRANSCRIPTION_ELONGATION')
        self.transcript_elongation.subsytem = 'Transcription'
        
        rxn = dict()
        for ntp, base_letter in seq_metabolite_map.items():
            rxn[ntp] = -1*self.premrna_base_counts[base_letter]
        # pyrophosphate released per base added, -1 for 3/5' ends
        rxn[ppi_n] = self.premrna_seq_length - 1
        rxn[self.premrna_transcript] = 1
        # ATP consumption due to PTMs of nucleosomes
        # https://www.pnas.org/content/pnas/suppl/2015/10/29/1514974112.DCSupplemental/pnas.1514974112.sapp.pdf
        # can perhaps add later

        self.transcript_elongation.add_metabolites(rxn)
        self.transcript_elongation.gene_reaction_rule = ' and '.join(ec) # GPRs
                
    def build_transcript_processing_reaction(self):
        '''Processing includes capping, splicing, and polyA tail.'''
        
        # combine in to one to not create too many reactions
        # capping itself is 4 reactions
        
        self.transcript_processing = cobra.Reaction(self.id + '_TRANSCRIPTION_PROCESSING')
        self.transcript_processing.subsytem = 'Transcription'
        rxn = dict()
        rxn[atp_n], rxn[ppi_n] = -self.polyA_length, self.polyA_length # polyA tail 
        
        # 5' cap: https://sites.google.com/site/learnorganicchem/organic-molecules/biomolecules/rna/rna-processing?tmpl=%2Fsystem%2Fapp%2Ftemplates%2Fprint%2F&showPrintDialog=1
        rxn[h2o_n], rxn[pi_n] = -1, 1 #rtpase
        rxn[gtp_n] = -1 #gp transfer
        rxn[ppi_n] += 1 # gp transfer
        rxn[amet_n], rxn[ahcys_n] = -2, 2 # methyltransferase - cap0 and cap1 structure
        rxn[h_n] = 1 # methyltransferase cap1
        
        #4 ATP consumed per capping reaction
        rxn[atp_n] -= 4 
        rxn[h2o_n] -= 4
        rxn[adp_n] = 4
        rxn[pi_n] += 4
        rxn[h_n] += 4

        # transcripts
        rxn[self.premrna_transcript] = -1
        rxn[self.mrna_transcript_n] = 1
        
        # splicing
        if self.lariats != None:
            rxn[self.lariats] = 1
    
            # 10 ATP consumed per intron during splicing
            rxn[atp_n] -= 10*self.n_lariats 
            rxn[h2o_n] -= 10*self.n_lariats
            rxn[adp_n] += 10*self.n_lariats
            rxn[pi_n] += 10*self.n_lariats
            rxn[h_n] += 10*self.n_lariats
            
            
        self.transcript_processing.add_metabolites(rxn)
        self.transcript_processing.gene_reaction_rule = ' and '.join(polyA + capping + spliceosome) # GPRs
    
    def build_lariat_degradation_reaction(self):
        # diagram: fig2 https://schoolbag.info/chemistry/chemical_biology/148.html
        if self.lariats != None:

            # approach 2: disregard lariat linearization as it is one phisphodiester bond break
            self.lariat_degradation = cobra.Reaction(self.id + '_LARIAT_DEGRADATION')
            self.lariat_degradation.subsytem = 'Transcription'
            rxn = dict()
            rxn[h2o_n] = -self.lariat_length # why doesn't this require -1?
            rxn[h_n] = self.lariat_length
            rxn[self.lariats] = -1
            for k,v in nmp_map_n.items():
                rxn[v] = self.lariat_base_counts[k]

            self.lariat_degradation.add_metabolites(rxn)
            self.lariat_degradation.gene_reaction_rule = lm_rule
            
        else:
            self.lariat_degradation = None
    def build_transcript_export_reaction(self):
        # make the cytosolic mrna metabolite
        self.mrna_transcript_c = self.mrna_transcript_n.copy()
        self.mrna_transcript_c.id = self.mrna_transcript_c.id.replace('[n]', '[c]')
        self.mrna_transcript_c.compartment = 'c'

        # make the transport reaction
        self.transcript_export = cobra.Reaction(self.id + '_TRANSCRIPTtn')
        self.transcript_export.name = 'mRNA nuclear export'
        self.transcript_export.subsytem = 'Transcription'
        rxn = dict()
        rxn[self.mrna_transcript_n], rxn[self.mrna_transcript_c] = -1, 1
        # 10 ATP consumer per transcript exported
        rxn[atp_n], rxn[h2o_n], rxn[adp_n], rxn[pi_n], rxn[h_n] = -10, -10, 10, 10, 10

        
        self.transcript_export.add_metabolites(rxn)
        # can change this GPR as an if statement in future based on following source:
        # https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.0010144.g005
        self.transcript_export.gene_reaction_rule = ' and '.join(trex)

    def build_transcript_degradation_reactions(self):
        self.transcript_degradation = cobra.Reaction(self.id + '_transcript_degradation')
        self.transcript_degradation.subsytem = 'Transcription'

        rxn = dict()
        rxn[h2o_c] = -(self.mrna_seq_length + self.polyA_length)+1
        rxn[h_c] = self.mrna_seq_length + self.polyA_length-1
        rxn[self.mrna_transcript_c] = -1
        for k,v in nmp_map_c.items():
            rxn[v] = self.mrna_base_counts[k]
        rxn[amp_c] += self.polyA_length

        # no m7g metabolite in recon2.2, so just reverse the methylation instead
        rxn[amet_c], rxn[ahcys_c] = 2, -2 # reverse methyltransferase - cap0 and cap1 structure

        # 3'-->5' degradation------------------------------------------------------------------------------------
        self.transcript_degradation_1 = cobra.Reaction(self.id + "_3'to5'_transcript_degradation")
        self.transcript_degradation_1.subsytem = 'Transcription'

        rxn_1 = rxn.copy()
        # fig 1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6565619/
        # 5' cap - degradation from 3'-->5' direction (scavenger - HIT mechanism)
        rxn_1[nmp_map_c[self.mrna_seq[0]]] -= 1
        rxn_1[ndp_map_c[self.mrna_seq[0]]] = 1
        rxn_1[h2o_c] -= 1
        rxn_1[h_c] += 1
        rxn_1[nmp_map_c['G']] += 1
        self.transcript_degradation_1.add_metabolites(rxn_1)
        self.transcript_degradation_1.gene_reaction_rule = degradation_rule1


        # 5'->3' degradation (decapping) ------------------------------------------------------------------------------------
        self.transcript_degradation_2_decapping = cobra.Reaction(self.id + "_decapping_transcript_degradation")
        self.transcript_degradation_2_decapping.subsytem = 'Transcription'

        rxn_2_decapping = rxn.copy()
        # 5' cap - from 5'-->3' direction (DCP1/DCP2 - NUDIX mechanism)
        rxn_2_decapping[h2o_c] -= 1
        rxn_2_decapping[h_c] += 1
        rxn_2_decapping[ndp_map_c['G']] = 1


        self.transcript_degradation_2_decapping.add_metabolites(rxn_2_decapping)
        self.transcript_degradation_2_decapping.gene_reaction_rule = decapping_rule



# In[4]:


#lariat degradation in two reactions
#             # approach 1: split up into two reactions
            
            
#             # rxn 1: lariat linearization by dbr1
#             self.linearized_lariats = self.lariats.copy()
#             self.linearized_lariats.id = self.linearized_lariats.id.split('_')[0] + '_linearized_lariats'
#             elements = self.linearized_lariats.elements.copy()
#             elements['H'] += self.n_lariats
#             elements['O'] += self.n_lariats
#             self.linearized_lariats.elements = elements
#             self.linearized_lariats.charge -= self.n_lariats


#             self.lariat_linearization = cobra.Reaction(self.id + '_lariat_linearization')
#             self.lariat_linearization.subsytem = 'Transcription'
#             rxn = dict()
#             rxn[h2o_n] = -self.n_lariats
#             rxn[h_n] = self.n_lariats
#             rxn[self.lariats] = -1
#             rxn[self.linearized_lariats] = 1
#             self.lariat_linearization.add_metabolites(rxn)
#             self.lariat_linearization.gene_reaction_rule = lariat_linearization[0]
            
#             # rxn 2: lariat degradation by exonucleases - not quite mass balanced, must fix if use this approach
#             self.lariat_degradation = cobra.Reaction(self.id + '_lariat_degradation')
#             self.lariat_degradation.subsytem = 'Transcription'
#             rxn = dict()
#             rxn[h2o_n] = -self.lariat_length
#             rxn[h_n] = self.lariat_length
#             rxn[self.linearized_lariats] = -1
#             for k,v in nmp_map.items():
#                 rxn[v] = self.lariat_base_counts[k]

#             self.lariat_degradation.add_metabolites(rxn)
#             self.lariat_degradation.gene_reaction_rule = ' and '.join(lariat_nucleases)


