#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
import pandas as pd 

from Bio.Seq import Seq
from Bio.Alphabet import generic_rna

import numpy as np
import statsmodels.api as sm
import scipy.stats as st

import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
# from utils.load_environmental_variables import *
from utils import machinery as mach
from utils import parameters as params
from utils import metabolites as metab
from utils import functions as func
from utils.polyA_statistics import calculate_polyA_length
from uniform_processes import biomass

from macromolecules.RNA import RNA_fragment, pre_mRNA, mRNA


# In[2]:


class express_mrna():
    def __init__(self, gene_info):
        self.gene_info = gene_info
        self.reactions = []
    
    def transcribe_premrna(self):
        # elongation reaction
        # https://www.google.com/search?q=rna+polymerization+reaction&source=lnms&tbm=isch&sa=X&ved=2ahUKEwiN_73Vk7rqAhXOsJ4KHW5lB4UQ_AUoAXoECA4QAw&biw=1920&bih=1001#imgrc=w7XH4mHmJglCuM
        
        self.premrna = pre_mRNA(self.gene_info)
        self.transcript_elongation = self.premrna.synthesize(id_ = self.gene_info.hgnc_id + '_TRANSCRIPTION_ELONGATION')
        self.reactions.append(self.transcript_elongation)   
        
    def process_mrna(self):
        '''Processing includes capping, splicing, and polyA tail.'''
        # combine in to one to not create too many reactions (capping itself is 4 reactions)
        # make mrna_n metabolite
        self.mrna_n = mRNA(self.gene_info, compartment = 'n')
        self.polyA_length = int(calculate_polyA_length(self.gene_info.polyA_length))
        self.mrna_n.update_metabolite(seq = ''.join(['A']*self.polyA_length), 
                                     append = True, append_to = '3_primed')
        
        #+2 for cap
        self.mrna_n.charge += 2#(-self.polyA_length + 2) 

        transcript_processing = cobra.Reaction(self.gene_info.hgnc_id + '_TRANSCRIPTION_PROCESSING')
        transcript_processing.subsytem = 'mRNA_expression'
        rxn = dict()
        
        rxn[metab.atp_n], rxn[metab.ppi_n] = -self.polyA_length, self.polyA_length # polyA tail 

        # 5' cap: https://sites.google.com/site/learnorganicchem/organic-molecules/biomolecules/rna/rna-processing?tmpl=%2Fsystem%2Fapp%2Ftemplates%2Fprint%2F&showPrintDialog=1
        rxn[metab.h2o_n], rxn[metab.pi_n] = -1, 1 #rtpase
        rxn[metab.gtp_n] = -1 #gp transfer
        rxn[metab.ppi_n] += 1 # gp transfer
        rxn[metab.amet_n], rxn[metab.ahcys_n] = -2, 2 # methyltransferase - cap0 and cap1 structure
        rxn[metab.h_n] = 1 # methyltransferase cap1

        #4 ATP consumed per capping reaction
        rxn = func.hydrolyze_atp(rxn, n_atp = 4, compartment = 'n')

        processed_elements = self.mrna_n.elements.copy()
        for element in processed_elements.keys():
#             processed_elements[element] += (self.polyA_length* metab.seq_element_map['A'][element]) # polyA tail
            processed_elements[element] += metab.gp[element] #5' cap rxn2 - addition of Gp

        # 5' cap 
        processed_elements['P'] -= 1 # rxn 1: lost of third triphosphate by RTPase
        processed_elements['O'] -= 4 # rxn 1: loss of third triophosphate by RTPase
        processed_elements['C'] += 2 # rxn 3-4: methyltransferase - cap0 and cap1 structure
        processed_elements['H'] += 5 # methyltransferase - cap0 and cap1 structure
        self.mrna_n.elements = processed_elements

        rxn[self.premrna], rxn[biomass.premrna_] = -1, -self.premrna.formula_weight/1000
        rxn[self.mrna_n], rxn[biomass.mrna_] = 1, self.mrna_n.formula_weight/1000

        # splicing
        if self.premrna.length > self.mrna_n.length - self.polyA_length: 
            lariat_seq = ''
            for nt in ['A', 'U', 'G', 'C']:

                if nt != 'A':
                    diff = self.premrna.sequence.count(nt) - self.mrna_n.sequence.count(nt)
                    lariat_seq += ''.join([nt]*diff)
                else:
                    diff = self.premrna.sequence.count(nt) - (self.mrna_n.sequence.count(nt)-self.polyA_length)
                    lariat_seq += ''.join([nt]*diff)
            
            self.lariat = RNA_fragment(metabolite_name = self.gene_info.hgnc_id, fragment_type='lariat', 
                                       seq = lariat_seq, triphosphate = False)

            if self.gene_info.n_introns == None:
                n_lariats = self.premrna.length * params.rate_intron # removed ROUND()
                if n_lariats < 1: # atleast one intron must be generated
                    n_lariats = 1
            #             n_splices = n_lariats * 2 # because on average, one more exon than intron
            else:
                n_lariats = self.gene_info.n_introns  

            rxn[self.lariat], rxn[biomass.other_rna_] = 1, self.lariat.formula_weight/1000
            rxn[metab.h2o_n] -= 1 # endonucleolytic cleavage
            # 10 ATP consumed per intron during splicing
            rxn = func.hydrolyze_atp(rxn, n_atp = 10*n_lariats, compartment = 'n')
            # lariat degradation - no linearization reaction (just one triphosphate consumption)
            lariat_degradation = self.lariat.exonucleolytic_degradation(reaction_name = self.gene_info.hgnc_id + '_lariats')
            lariat_degradation.subsystem = 'mRNA_expression'
            lariat_degradation.gene_reaction_rule = mach.lm_rule  
            if list(lariat_degradation.compartments) != ['n']:
                raise ValueError('Lariat degradation must be confined to nuclear compartment')
            else:
                self.reactions.append(lariat_degradation)
        else:
            lariat_degradation = None

        transcript_processing.add_metabolites(rxn)
        transcript_processing.gene_reaction_rule = ' and '.join(mach.polyA + mach.capping + mach.spliceosome) 
        if len(transcript_processing.check_mass_balance()) > 0:
            raise ValueError('Transcript processing for ' + self.gene_info.hgnc_id + ' is unbalanced')
        elif list(transcript_processing.compartments) != ['n']:
            raise ValueError('Transcript processing must be confined to nuclear compartment')
        else:
            self.transcript_processing = transcript_processing
            self.reactions.append(transcript_processing)
    def export_mrna(self):
        # make the cytosolic mrna metabolite
        self.mrna_c = self.mrna_n.copy()
        self.mrna_c.id = self.mrna_c.id.replace('[n]', '[c]')
        self.mrna_c.compartment = 'c'

        # make the transport reaction
        mrna_export = cobra.Reaction(self.gene_info.hgnc_id + '_mRNA_EXPORTtn')
        mrna_export.name = 'mRNA nuclear export'
        mrna_export.subsytem = 'mRNA_expression'
        rxn = dict()
        rxn[self.mrna_n], rxn[self.mrna_c] = -1, 1
        
#         # NEW
#         self.mrna_dilution_proxy = cobra.Metabolite(self.gene_info.hgnc_id + '_mrna_dilution_proxy')
#         rxn[self.mrna_dilution_proxy] = 1 
        
        # 10 ATP consumer per transcript exported
        rxn = func.hydrolyze_atp(rxn, n_atp = 10, compartment = 'n')

        mrna_export.add_metabolites(rxn)
        # can change this GPR as an if statement in future based on following source:
        # https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.0010144.g005
        mrna_export.gene_reaction_rule = ' and '.join(mach.trex)
        
        if len(mrna_export.check_mass_balance()) > 0:
            raise ValueError('mRNA export for ' + self.gene_info.hgnc_id + ' is unbalanced')
        else:
            self.mrna_export = mrna_export
            self.reactions.append(mrna_export)
            
    def degrade_mrna(self, decapping = True, three_to_five = False):
        '''

        Right now, only one of the two degradation pathways is included. We assume the 5' to 3' pathway is present.
        This is simply to limit the total number of reactions

        '''
        
        rxn = self.mrna_c.exonucleolytic_degradation(reaction_name = '', balanced = False)
        rxn = rxn.metabolites.copy()
        del rxn[[m for m in rxn.keys() if m.id == metab.ntp_map_c[self.gene_info.mrna_seq[0]].id][0]]

        # no m7g metabolite in recon2.2, so just reverse the methylation instead
        rxn[metab.amet_c], rxn[metab.ahcys_c] = 2, -2 # reverse methyltransferase - cap0 and cap1 structure

        # proxy metabolite for coupling mRNA degradation to protein synthesis flux
        self.mrna_deg_proxy = cobra.Metabolite(self.gene_info.hgnc_id + '_mrna_deg_proxy')
        rxn[self.mrna_deg_proxy] = 1 

        h2o_c = [m for m in rxn.keys() if m.id == 'h2o[c]'][0] # won't load directly from metab for some reason
        h_c = [m for m in rxn.keys() if m.id == 'h[c]'][0]

        if three_to_five: 
            transcript_degradation_1 = cobra.Reaction(self.gene_info.hgnc_id + "_3'to5'_mRNA_DEGRADATIONc")
            transcript_degradation_1.subsytem = 'mRNA_expression'
            rxn_1 = rxn.copy()

            rxn_1[metab.ndp_map_c[self.gene_info.mrna_seq[0]]] = 1

            gmp_c = [m for m in rxn.keys() if m.id == 'gmp[c]'][0]
            rxn_1[h2o_c] -= 1
            rxn_1[h_c] += 1
            rxn_1[gmp_c] += 1
            transcript_degradation_1.add_metabolites(rxn_1)
            transcript_degradation_1.gene_reaction_rule = mach.degradation_rule1

            if len(transcript_degradation_1.check_mass_balance()) > 0:
                raise ValueError('3 primed to 5 primed degradation for ' + self.gene_info.hgnc_id + ' is unbalanced')
            elif list(transcript_degradation_1.compartments) != ['c']:
                raise ValueError('Transcript degradation must be confined to cytosolic compartment')
            else:
                self.reactions.append(transcript_degradation_1)
        if decapping:
            transcript_degradation_2_decapping = cobra.Reaction(self.gene_info.hgnc_id + "_DECAPPING_mRNA_DEGRADATIONc")
            transcript_degradation_2_decapping.subsytem = 'mRNA_expression'

            rxn_2 = rxn.copy()
            rxn_2[[m for m in rxn_2.keys() if m.id == metab.nmp_map_c[self.gene_info.mrna_seq[0]].id][0]] += 1
            # 5' cap - from 5'-->3' direction (DCP1/DCP2 - NUDIX mechanism)
            rxn_2[h2o_c] -= 1
            rxn_2[h_c] += 1
            rxn_2[metab.ndp_map_c['G']] = 1


            transcript_degradation_2_decapping.add_metabolites(rxn_2)
            transcript_degradation_2_decapping.gene_reaction_rule = mach.decapping_rule
            if len(transcript_degradation_2_decapping.check_mass_balance()) > 0:
                raise ValueError('Decapping degradation for ' + self.gene_info.hgnc_id + ' is unbalanced')
            elif list(transcript_degradation_2_decapping.compartments) != ['c']:
                raise ValueError('Transcript degradation must be confined to cytosolic compartment')
            else:
                self.reactions.append(transcript_degradation_2_decapping)
    def compress_mrna_module(self):
        rxns_to_remove = [self.transcript_elongation, self.transcript_processing, self.mrna_export]
        rxn = dict()
        rxn_map = dict()
        for r in rxns_to_remove:
            for met, coeff in r.metabolites.items():
                if met.id in rxn.keys():
                    rxn[met.id] += coeff
                else:
                    rxn[met.id] = coeff
                    rxn_map[met.id] = met
        rxn = {rxn_map[k]: v for k,v in rxn.items()}

        transcription = cobra.Reaction(self.gene_info.hgnc_id + "_TRANSCRIPTION")
        transcription.subsytem = 'mRNA_expression'
        transcription.add_metabolites(rxn)

        transcription.gene_reaction_rule = ' and '.join(sorted(set([item.id for sublist in [list(r.genes) for r in rxns_to_remove] for item in sublist])))

        if len(transcription.check_mass_balance()) > 0:
            raise ValueError('Condensed transcription reaction for ' + self.gene_info.hgnc_id + ' is unbalanced')
        
        for r in rxns_to_remove:
            self.reactions.remove(r)
        self.reactions.append(transcription)


# In[3]:


def get_mrna_expression_reactions(gene_info, compress_mrna = False):
    '''
    gene_info is an object of the gene_information class. Returns all reactions associated with mRNA
    expression and necessary metabolites for other modules in the ME model. 
    
    compress_mrna is a boolean. If True, will make the transcription, processing, and export reactions one 
    single reaction
    
    '''
    self = express_mrna(gene_info)
    self.transcribe_premrna()
    self.process_mrna()
    self.export_mrna()
    self.degrade_mrna()
    if compress_mrna:
        self.compress_mrna_module()

#     return self.reactions, self.mrna_dilution_proxy, self.mrna_deg_proxy 
    return self.reactions, self.mrna_c, self.mrna_deg_proxy 

