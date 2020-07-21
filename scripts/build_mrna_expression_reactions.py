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
from gene_information import gene_information
from utils import *


# # start of actual script

# In[58]:


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
            rxn[h2o_n] = -self.lariat_length # no -1 bc of additional 2 to 5' bond 
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
#         self.transcript_degradation = cobra.Reaction(self.id + '_TRANSCRIPT_DEGRADATION')
#         self.transcript_degradation.subsytem = 'Transcription'

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
        self.transcript_degradation_1 = cobra.Reaction(self.id + "_3'to5'_TRANSCRIPT_DEGRADATION")
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
        self.transcript_degradation_2_decapping = cobra.Reaction(self.id + "_DECAPPING_TRANSCRIPT_DEGRADATION")
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


