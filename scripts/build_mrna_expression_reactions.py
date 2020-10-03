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
sys.path.insert(1, '../scripts/') # comment out in python script
# from utils.load_environmental_variables import *
from utils import machinery as mach
from utils import parameters as params
from utils import metabolites as metab
from utils import functions as func
from utils.polyA_statistics import calculate_polyA_length
from uniform_processes import biomass
# from gene_information import gene_information


# In[6]:


def transcribe_premrna(gene_info):
    # elongation reaction
    # https://www.google.com/search?q=rna+polymerization+reaction&source=lnms&tbm=isch&sa=X&ved=2ahUKEwiN_73Vk7rqAhXOsJ4KHW5lB4UQ_AUoAXoECA4QAw&biw=1920&bih=1001#imgrc=w7XH4mHmJglCuM
    premrna_transcript_n, premrna_base_counts = func.make_rna_metabolite(gene_info.hgnc_id + '_premrna', 
                                            gene_info.premrna_seq, molecule_type = 'mrna', compartment = 'n', 
                                            triphosphate = True)  
    premrna_mw = func.get_metabolite_mw(premrna_transcript_n)
        
    transcript_elongation = cobra.Reaction(gene_info.hgnc_id + '_TRANSCRIPTION_ELONGATION')
    transcript_elongation.subsytem = 'mRNA_expression'
        
    rxn = dict()
    for ntp, base_letter in metab.seq_metabolite_map.items():
        rxn[ntp] = -1*premrna_base_counts[base_letter]
    # pyrophosphate released per base added, -1 for 3/5' ends
    L_premrna = len(gene_info.premrna_seq)
    rxn[metab.ppi_n] = L_premrna - 1
    rxn[premrna_transcript_n] = 1
    rxn[biomass.premrna_] = premrna_mw # biomass reaction
    # ATP consumption due to PTMs of nucleosomes
    # https://www.pnas.org/content/pnas/suppl/2015/10/29/1514974112.DCSupplemental/pnas.1514974112.sapp.pdf
    # can perhaps add later

    transcript_elongation.add_metabolites(rxn)
    transcript_elongation.gene_reaction_rule = ' and '.join(mach.ec) # GPRs
    
    return transcript_elongation, premrna_transcript_n, L_premrna, premrna_base_counts, premrna_mw

def process_mrna(gene_info, premrna_transcript_n, L_premrna, premrna_base_counts, premrna_mw):
        '''Processing includes capping, splicing, and polyA tail.'''
        # combine in to one to not create too many reactions (capping itself is 4 reactions)
        

        transcript_processing = cobra.Reaction(gene_info.hgnc_id + '_TRANSCRIPTION_PROCESSING')
        transcript_processing.subsytem = 'mRNA_expression'
        rxn = dict()
        
        polyA_length = calculate_polyA_length(gene_info.polyA_length)
        rxn[metab.atp_n], rxn[metab.ppi_n] = -polyA_length, polyA_length # polyA tail 
        
        # 5' cap: https://sites.google.com/site/learnorganicchem/organic-molecules/biomolecules/rna/rna-processing?tmpl=%2Fsystem%2Fapp%2Ftemplates%2Fprint%2F&showPrintDialog=1
        rxn[metab.h2o_n], rxn[metab.pi_n] = -1, 1 #rtpase
        rxn[metab.gtp_n] = -1 #gp transfer
        rxn[metab.ppi_n] += 1 # gp transfer
        rxn[metab.amet_n], rxn[metab.ahcys_n] = -2, 2 # methyltransferase - cap0 and cap1 structure
        rxn[metab.h_n] = 1 # methyltransferase cap1
        
        #4 ATP consumed per capping reaction
        rxn = func.hydrolyze_atp(rxn, n_atp = 4, compartment = 'n')

        # transcripts
        mrna_transcript_n, mrna_base_counts = func.make_rna_metabolite(gene_info.hgnc_id, gene_info.mrna_seq, 
                                              molecule_type = 'mrna', compartment = 'n',triphosphate = True) 
        processed_elements = mrna_transcript_n.elements
        for element in processed_elements.keys():
            processed_elements[element] += (polyA_length* metab.seq_element_map['A'][element]) # polyA tail
            processed_elements[element] += metab.gp[element] #5' cap rxn2 - addition of Gp
        # 5' cap 
        processed_elements['P'] -= 1 # rxn 1: lost of third triphosphate by RTPase
        processed_elements['O'] -= 4 # rxn 1: loss of third triophosphate by RTPase
        processed_elements['C'] += 2 # rxn 3-4: methyltransferase - cap0 and cap1 structure
        processed_elements['H'] += 5 # methyltransferase - cap0 and cap1 structure
        mrna_transcript_n.elements = processed_elements
        
        #+2 for cap
        mrna_transcript_n.charge += (-polyA_length + 2) 
        
        mw_mrna = func.get_metabolite_mw(mrna_transcript_n)
        rxn[premrna_transcript_n], rxn[biomass.premrna_] = -1, -premrna_mw
        rxn[mrna_transcript_n], rxn[biomass.mrna_] = 1, mw_mrna
        
        processing_reactions = list()
        
        # splicing
        L_mrna = len(gene_info.mrna_seq)
        if L_premrna > L_mrna: 
            lariat_seq = ''.join([k*(v - mrna_base_counts[k]) for k,v in premrna_base_counts.items()])
            lariats_n, lariats_base_counts = func.make_rna_metabolite(gene_info.hgnc_id + '_lariats', 
                                             lariat_seq,molecule_type = 'mrna', compartment = 'n',
                                             triphosphate = False)
            
            if gene_info.n_introns == None:
                n_lariats = round(L_premrna * params.rate_intron)
                if n_lariats < 1: # atleast one intron must be generated
                    n_lariats = 1
    #             n_splices = n_lariats * 2 # because on average, one more exon than intron
            else:
                n_lariats = gene_info.n_introns                    
            
            mw_lariats = func.get_metabolite_mw(lariats_n)
            rxn[lariats_n], rxn[biomass.other_rna_] = 1, mw_lariats
            rxn[metab.h2o_n] -= 1 # endonucleolytic cleavage
            # 10 ATP consumed per intron during splicing
           
            rxn = func.hydrolyze_atp(rxn, n_atp = 10*n_lariats, compartment = 'n')
            
            # lariat degradation - no linearization reaction (just one triphosphate consumption)
            lariat_degradation = func.rna_exonucleolytic_degradation(lariats_n, lariats_base_counts, lariat_seq, 
                                  gene_info.hgnc_id + '_lariats',triphosphate = False, nucleus = True)
            lariat_degradation.add_metabolites({biomass.other_rna_: - mw_lariats})
            lariat_degradation.subsystem = 'mRNA_expression'
            lariat_degradation.gene_reaction_rule = mach.lm_rule
            processing_reactions += [lariat_degradation]
        else:
            lariat_degradation = None

        transcript_processing.add_metabolites(rxn)
        transcript_processing.gene_reaction_rule = ' and '.join(mach.polyA + mach.capping + mach.spliceosome) 
        processing_reactions += [transcript_processing]   
        
        return processing_reactions, mrna_transcript_n, polyA_length, L_mrna, mrna_base_counts, mw_mrna

def export_mrna(gene_info, mrna_transcript_n):
    # make the cytosolic mrna metabolite
    mrna_transcript_c = mrna_transcript_n.copy()
    mrna_transcript_c.id = mrna_transcript_c.id.replace('[n]', '[c]')
    mrna_transcript_c.compartment = 'c'

    # make the transport reaction
    mrna_export = cobra.Reaction(gene_info.hgnc_id + '_mRNA_EXPORTtn')
    mrna_export.name = 'mRNA nuclear export'
    mrna_export.subsytem = 'mRNA_expression'
    rxn = dict()
    rxn[mrna_transcript_n], rxn[mrna_transcript_c] = -1, 1
    # 10 ATP consumer per transcript exported
    rxn = func.hydrolyze_atp(rxn, n_atp = 10, compartment = 'n')

    mrna_export.add_metabolites(rxn)
    # can change this GPR as an if statement in future based on following source:
    # https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.0010144.g005
    mrna_export.gene_reaction_rule = ' and '.join(mach.trex)
    
    return mrna_export, mrna_transcript_c
        

def degrade_mrna(gene_info, mrna_transcript_c, polyA_length, L_mrna, mrna_base_counts, mw_mrna, decapping = True, 
                three_to_five = False):
    '''
    
    Right now, only one of the two degradation pathways is included. We assume the 5' to 3' pathway is present.
    This is simply to limit the total number of reactions
    
    '''
    
    degradation_reactions = list()
    
    
    rxn = dict()
    rxn[metab.h2o_c] = -(L_mrna + polyA_length)+1
    rxn[metab.h_c] = L_mrna + polyA_length-1
    rxn[mrna_transcript_c], rxn[biomass.mrna_] = -1, -mw_mrna
    
    for k,v in metab.nmp_map_c.items():
        rxn[v] = mrna_base_counts[k]
    rxn[metab.amp_c] += polyA_length

    # no m7g metabolite in recon2.2, so just reverse the methylation instead
    rxn[metab.amet_c], rxn[metab.ahcys_c] = 2, -2 # reverse methyltransferase - cap0 and cap1 structure
    
    # proxy metabolite for coupling mRNA degradation to protein synthesis flux
    mrna_deg_proxy = cobra.Metabolite(gene_info.hgnc_id + 'mrna_deg_proxy')
    rxn[mrna_deg_proxy] = 1 
    
    

    # 3'-->5' degradation------------------------------------------------------------------------------------
    if three_to_five:
        transcript_degradation_1 = cobra.Reaction(gene_info.hgnc_id + "_3'to5'_mRNA_DEGRADATIONc")
        transcript_degradation_1.subsytem = 'Transcription'

        rxn_1 = rxn.copy()
        # fig 1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6565619/
        # 5' cap - degradation from 3'-->5' direction (scavenger - HIT mechanism)
        rxn_1[metab.nmp_map_c[gene_info.mrna_seq[0]]] -= 1
        rxn_1[metab.ndp_map_c[gene_info.mrna_seq[0]]] = 1
        rxn_1[metab.h2o_c] -= 1
        rxn_1[metab.h_c] += 1
        rxn_1[metab.nmp_map_c['G']] += 1
        transcript_degradation_1.add_metabolites(rxn_1)
        transcript_degradation_1.gene_reaction_rule = mach.degradation_rule1
        degradation_reactions.append(transcript_degradation_1)


    # 5'->3' degradation (decapping) ------------------------------------------------------------------------------------
    if decapping:
        transcript_degradation_2_decapping = cobra.Reaction(gene_info.hgnc_id + "_DECAPPING_mRNA_DEGRADATIONc")
        transcript_degradation_2_decapping.subsytem = 'Transcription'

        rxn_2_decapping = rxn.copy()
        # 5' cap - from 5'-->3' direction (DCP1/DCP2 - NUDIX mechanism)
        rxn_2_decapping[metab.h2o_c] -= 1
        rxn_2_decapping[metab.h_c] += 1
        rxn_2_decapping[metab.ndp_map_c['G']] = 1


        transcript_degradation_2_decapping.add_metabolites(rxn_2_decapping)
        transcript_degradation_2_decapping.gene_reaction_rule = mach.decapping_rule
        degradation_reactions.append(transcript_degradation_2_decapping)
    
    return degradation_reactions, mrna_deg_proxy
    
    
def get_mrna_expression_reactions(gene_info):
    '''
    gene_info is an object of the gene_information class. Returns all reactions associated with mRNA
    expression and necessary metabolites for other modules in the ME model. 
    
    '''
    transcript_elongation, premrna_transcript_n, L_premrna, premrna_base_counts, premrna_mw = transcribe_premrna(gene_info)
    
    processing_res = process_mrna(gene_info, premrna_transcript_n, L_premrna, premrna_base_counts, premrna_mw)
    processing_reactions, mrna_transcript_n, polyA_length, L_mrna, mrna_base_counts, mw_mrna = processing_res
    
    mrna_export, mrna_transcript_c = export_mrna(gene_info, mrna_transcript_n)
    degradation_reactions, mrna_deg_proxy = degrade_mrna(gene_info, mrna_transcript_c, polyA_length, L_mrna, mrna_base_counts, mw_mrna)
    
    
    
    reactions = [transcript_elongation] + processing_reactions + [mrna_export] + degradation_reactions
    return reactions, mrna_transcript_c, mrna_deg_proxy 

