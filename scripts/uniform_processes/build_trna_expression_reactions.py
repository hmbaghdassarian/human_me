#!/usr/bin/env python
# coding: utf-8

# In[17]:


import warnings
import cobra

import pandas as pd

import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils.load_environmental_variables import *
from utils import machinery as mach
from utils import parameters as params
from utils import metabolites as metab
from utils import functions as func

from uniform_processes import biomass


# -To change to individual tRNA molecules rather than a generic one, start with the charge_trna function and also change the trna_biogenesis function
# 
# -To add modifications, will need to change the allowed_trna_modifications dictionary in utils, and edit the following functions: 1) modify_trna_nuclear, 2) modify_trna_cytosolic, 3) degrade_trna

# # TRNA Information Class

# In[18]:


class trna_information():
    def __init__(self, maturetrna_sequence , id_, three_trailer_seq = None , five_leader_seq = None, 
                 modifications = {}, intron_sequences = None):
        '''
        1) Mature trna sequence is the RNA sequence of the final processed tRNA, represented as a string from 5'
        to 3' end, including CCA.
        2) id_ should be map to isodecoder somehow
        2-3) five_leader_seq and three_trailer_seq are the RNA sequences of the 5' leader and 3' trailer sequences that are 
        excised. Represented as a string from 5' to 3' end.
        5) Modifications is a dictionary with keys as possible modifications (string) and values as the 
        number of modifications that occur (integer). See allowed_trna_modifications for all trna modifications
        incorporated in this model. 
        6) Intron sequences is a list of RNA sequence corresponding to each intron. 
        '''
        
        
        if maturetrna_sequence[-3:] != 'CCA':
            warnings.warn('CCA tail not present in provided mature sequence, adding to 3 primed end')
            maturetrna_sequence += 'CCA'
#         if len(maturetrna_sequence) < 73 or len(maturetrna_sequence) > 93:
#             # https://www.nature.com/articles/nrm.2017.77#Sec2
#             warnings.warn('Mature tRNA sequence not in the expected length range (76<=L<=93)')
            
#         if anticodon_sequence != None:
#             warnings.warn('Current iteration of ME-model synthesizes a generic tRNA, not a codon-specific one')
#             if anticodon_sequence not in maturetrna_sequence:
#                 raise ValueError('Anticodon sequence not in mature tRNA sequence')
                
        if len(set(modifications.keys()).difference(params.allowed_trna_modifications)) > 0:
            warning_ = 'At least one of the listed modifications is not currently considered in this model'
            warnings.warn(warning_)
            modifications = {k:v for k in modifications.keys() if k in params.allowed_trna_modifications.keys()}
        
        if intron_sequences != None and type(intron_sequences) != list:
            raise ValueError('intron_sequences must be a list of sequences, one for each intron')


        
        self.maturetrna_sequence = maturetrna_sequence
        self.id = id_
        self.three_trailer_seq = three_trailer_seq
        self.five_leader_seq = five_leader_seq
#         self.anticodon_sequence = anticodon_sequence
        self.modifications = modifications
        self.intron_sequences = intron_sequences
        
        if self.five_leader_seq != None:
            self.pretrna_sequence = self.five_leader_seq + self.maturetrna_sequence[:-3] 
        else:
            self.pretrna_sequence = self.maturetrna_sequence[:-3]
        
        if self.intron_sequences != None: # position doesn't matter, just need # of elements and mass balance
            self.pretrna_sequence += ''.join(self.intron_sequences)
            
        if self.three_trailer_seq != None:
            self.pretrna_sequence += self.three_trailer_seq   
        
        pretrna_base_counts, trna_base_counts = dict(), dict()
        for base_letter in metab.seq_element_map.keys():
            pretrna_base_counts[base_letter] = self.pretrna_sequence.count(base_letter)
            trna_base_counts[base_letter] = self.maturetrna_sequence.count(base_letter)
        for k,v in trna_base_counts.items():
            if v > pretrna_base_counts[k]:
                raise ValueError('Number of ' + k + ' bases in pretrna sequence less than that of trna sequence')

        self.pretrna_base_counts = pretrna_base_counts
        self.trna_base_counts = trna_base_counts


# # Reactions

# In[19]:


def update_trna_degradation(trna_degradation_reaction):
    trna_degradation_reaction.subsytem = 'tRNA_Biogenesis'
    trna_degradation_reaction.gene_reaction_rule = ' and '.join(mach.lariat_machinery['Exosome'])
    return trna_degradation_reaction


# In[45]:


def transcribe_pretrna(trna_info):
    pretrna_transcript_n, pretrna_base_counts = func.make_rna_metabolite(trna_info.id + '_pre', 
                                                trna_info.pretrna_sequence, molecule_type = 'trna',
                                                compartment = 'n', triphosphate = True)

    pretrna_transcription = cobra.Reaction('TRANSCRIPTION_PRE_TRNA_' + trna_info.id)
    pretrna_transcription.subsytem = 'tRNA_Biogenesis'
    rxn = dict()
    for ntp, base_letter in metab.seq_metabolite_map.items():
        rxn[ntp] = -1*pretrna_base_counts[base_letter]
    rxn[metab.ppi_n] = len(trna_info.pretrna_sequence) - 1
    pretrna_mw = func.get_metabolite_mw(pretrna_transcript_n)
    rxn[pretrna_transcript_n], rxn[biomass.trna_] = 1, pretrna_mw
    pretrna_transcription.add_metabolites(rxn)
    pretrna_transcription.gene_reaction_rule = ' and '.join(mach.rnap3_transcription_machinery)     
    return pretrna_transcription, pretrna_transcript_n, pretrna_mw

def process_trna(trna_info, pretrna_transcript_n, pretrna_mw):
    '''
    
    This reaction processes pre-tRNA into mature tRNA in the nucleus. 
    This includes: CCA synthesis, 5' leader and 3' trailer cleavage (and degradation as separate reactions),
    and splicing (and intron degradation as a separate reactions for each intron).
    
    '''
    # processing includes 5' leader and 3' trailer degradation, CCA synthesis
    # in the future, should inlcude splicing
    
    rxn = {pretrna_transcript_n: -1}
    # CCA synthesis
    rxn[metab.ntp_map_n['C']] = -2
    rxn[metab.ntp_map_n['A']] = -1
    rxn[metab.ppi_n] = 3
    trna_processing_machinery = mach.TRNT1.copy()
    
    # initialize
    reactions = list()
    rxn[metab.h2o_n] = 0 
    
    other_trna_biomass = 0 # initialize
    
    # 5' cleavage
    if trna_info.five_leader_seq != None: # if there is a 5' leader sequence
        five_frag_n, five_frag_base_counts = func.make_rna_metabolite(trna_info.id + "_5'_leader_fragment", 
                                             trna_info.five_leader_seq, compartment = 'n', molecule_type = 'trna',
                                             triphosphate = True)
        five_frag_mw = func.get_metabolite_mw(five_frag_n)
        other_trna_biomass += five_frag_mw
        rxn[five_frag_n] = 1
        rxn[metab.h2o_n] -= 1 #endonuclolytic cleavage (RNAse P)
        trna_processing_machinery += mach.RNASEP
        
        five_leader_degradation = func.rna_exonucleolytic_degradation(five_frag_n, five_frag_base_counts, trna_info.five_leader_seq, 
                              trna_info.id + "_5'_leader_fragment_tRNA", triphosphate = True, 
                                  nucleus = True)
        five_leader_degradation.add_metabolites({biomass.other_rna_: -five_frag_mw})
        five_leader_degradation = update_trna_degradation(five_leader_degradation)
        reactions += [five_leader_degradation]
        
        tp = False 
    else:
        tp = True

    # mature tRNA
    trna_transcript_n, trna_base_counts = func.make_rna_metabolite(trna_info.id, trna_info.maturetrna_sequence, 
                                          molecule_type = 'trna', compartment = 'n', triphosphate = tp)
    biomass_change = func.get_metabolite_mw(trna_transcript_n) - pretrna_mw
    rxn[trna_transcript_n], rxn[biomass.trna_] = 1, biomass_change

    # 3' cleavage
    if trna_info.three_trailer_seq != None:
        three_frag_n, three_frag_base_counts = func.make_rna_metabolite(trna_info.id + "_3'_trailer_fragment", 
                                               trna_info.three_trailer_seq, molecule_type = 'trna', 
                                               compartment = 'n',triphosphate = False)
        three_frag_mw = func.get_metabolite_mw(three_frag_n)
        other_trna_biomass += three_frag_mw
        rxn[three_frag_n] = 1
        rxn[metab.h2o_n] -= 1 #endonuclolytic cleavage (RNase Z)
        trna_processing_machinery += mach.RNASEZ
        
        
        three_trailer_degradation = func.rna_exonucleolytic_degradation(three_frag_n, three_frag_base_counts, 
                                    trna_info.three_trailer_seq,trna_info.id + "_3'_trailer_fragment_tRNA", 
                                    triphosphate = False,nucleus = True)
        three_trailer_degradation.add_metabolites({biomass.other_rna_: -three_frag_mw})
        three_trailer_degradation = update_trna_degradation(three_trailer_degradation)
        reactions += [three_trailer_degradation]

        

    # splicing of intron
    if trna_info.intron_sequences != None:
        n_introns = len(trna_info.intron_sequences)
        trna_introns_n = dict()
        for i in range(len(trna_info.intron_sequences)):
            trna_intron_n, trna_intron_base_counts = func.make_rna_metabolite(trna_info.id + "_intron_" + str(i), 
                                                     trna_info.intron_sequences[i], molecule_type = 'trna',
                                                     compartment = 'n',triphosphate = False)
            intron_mw = func.get_metabolite_mw(trna_intron_n)
            other_trna_biomass += intron_mw
            rxn[trna_intron_n] = 1
            trna_processing_machinery += mach.trna_splicing_machinery
            
            intron_degradation = func.rna_exonucleolytic_degradation(trna_intron_n, trna_intron_base_counts, 
                                                           trna_info.intron_sequences[i],
                                                           trna_info.id + "_intron_" + str(i) + '_tRNA', 
                                                           triphosphate = False, nucleus = True)
            intron_degradation.add_metabolites({biomass.other_rna_: -intron_mw})
            intron_degradation = update_trna_degradation(intron_degradation)
            reactions += [intron_degradation]

        rxn[h2o_n] -= n_introns
    
    if other_trna_biomass > 0:
        rxn[biomass.other_rna_] = other_trna_biomass
    trna_processing = cobra.Reaction('PROCESSING_TRNA_' + trna_info.id)
    trna_processing.subsytem = 'tRNA_Biogenesis'
    trna_processing.add_metabolites(rxn)
    trna_processing.gene_reaction_rule = ' and '.join(trna_processing_machinery)

    reactions += [trna_processing]
    
    return reactions, trna_transcript_n

def modify_trna_nuclear(trna_info, trna_transcript_n):
    '''Add to the if statement in the future. This is for nuclear modifications'''
    if len(trna_info.modifications) > 0:
        raise ValueError('Modifications are not currently considered')
#         modified_trna_transcript_n = trna_transcript_n.copy()
#         modified_trna_transcript_n.id = trna_info.id + '_modified_trna[n]'
        ####
    else:
        trna_modifications_nuclear = None
        modified_trna_transcript_n = trna_transcript_n
    return trna_modifications_nuclear, modified_trna_transcript_n

def primary_export_trna(trna_info, modified_trna_transcript_n):
    trna_transcript_c = modified_trna_transcript_n.copy()
    trna_transcript_c.id = trna_transcript_c.id.replace('[n]', '[c]')
    trna_transcript_c.compartment = 'c'

    trna_primary_export = cobra.Reaction(trna_info.id + 'PRIMARY_EXPORTtn')
    trna_primary_export.subsytem = 'tRNA_Biogenesis'
    trna_primary_export.name = 'trna nuclear export'
    
    export_rxn = {modified_trna_transcript_n: -1, trna_transcript_c: 1}
    # gtp hydrolysis on cytoplasmic side for export (see protein_expression nuclear_transport for details)
    export_rxn[metab.ntp_map_c['G']], export_rxn[metab.h2o_c], export_rxn[metab.ndp_map_c['G']], export_rxn[metab.pi_c], export_rxn[metab.h_c]  = -1, -1, 1, 1, 1
    trna_primary_export.add_metabolites(export_rxn)
    trna_primary_export.gene_reaction_rule = ' and '.join(mach.XPOT + mach.RAN)
    
    return trna_primary_export, trna_transcript_c


def modify_trna_cytosolic(trna_info, trna_transcript_c):
    '''Add to the if statement in the future. This is for cytosolic modifications'''
    if len(trna_info.modifications) > 0:
        raise ValueError('Modifications are not currently considered')
    else:
        trna_modifications_cytosolic = None
        modified_trna_transcript_c = trna_transcript_c
    
    modified_trna_transcript_c_mw = func.get_metabolite_mw(modified_trna_transcript_c)
    return trna_modifications_cytosolic, modified_trna_transcript_c, modified_trna_transcript_c_mw

def degrade_trna(trna_info, modified_trna_transcript_c, modified_trna_transcript_c_mw):
    # currently built on the assumption that there are no post-transcriptional modifications on trna 
    
    if trna_info.five_leader_seq != None: # if there is a 5' leader sequence
        tp = False 
    else:
        tp = True
    trna_degradation = func.rna_exonucleolytic_degradation(modified_trna_transcript_c, trna_info.trna_base_counts, 
                       trna_info.maturetrna_sequence, trna_info.id + '_tRNA', triphosphate = tp, nucleus = False)
    trna_degradation.add_metabolites({biomass.trna_: -modified_trna_transcript_c_mw})
    trna_degradation.subsytem = 'tRNA_Biogenesis'
    trna_degradation.gene_reaction_rule = mach.XRN1[0]
    return trna_degradation 
    

def charge_trna(trna_info, modified_trna_transcript_c, modified_trna_transcript_c_mw):
    '''tRNA charging reaction combines the activation and charging steps into one reaction.'''
    
    
    # in the future, this should take into account anticodon sequence, which should be in 
    # trna_info id
    
    # now, since just a generic charging reaction, will create one for each amino acid (for loop)
    
    # diagram: https://www.researchgate.net/figure/The-reaction-scheme-for-the-two-steps-of-aminoacylation-reaction-at-the-active-site-of_fig4_231225238
    trna_charging_reactions, charged_trna_metabolites = [], []
    for code, aa in metab.seq_amino_acid_map_c.items():
        elements = modified_trna_transcript_c.elements
        # attachment - loss of hydrogen from tRNA hydroxyl, and oxygen from amino acid carboxyl
        elements['H'] -= 1 
        elements['O'] -= 1
        for element, count in aa.elements.items():
            if element in elements.keys():
                elements[element] += count
            else:
                elements[element] = count

        charged_trna_c = cobra.Metabolite('charged_' + trna_info.id + '_' + code + '_trna[c]')
        charged_trna_c.compartment = 'c'
        charged_trna_c.elements = elements
        # +1 for loss of negative charge on oxygen of amino acid
        charged_trna_c.charge = modified_trna_transcript_c.charge + aa.charge + 1 

        trna_charging = cobra.Reaction('CHARGING_TRNA_' + trna_info.id + '_' + code)
        trna_charging.subsytem = 'tRNA_Biogenesis'
        biomass_change = func.get_metabolite_mw(charged_trna_c) - modified_trna_transcript_c_mw
        rxn = {modified_trna_transcript_c: -1, aa: -1, charged_trna_c: 1, metab.atp_c: -1, metab.ppi_c: 1, 
               metab.amp_c: 1, biomass.trna_: biomass_change}
        trna_charging.add_metabolites(rxn)
        # add gprs
        genes = mach.seq_synthetase_map[code]
        if len(genes) == 1:
            trna_charging.gene_reaction_rule = genes[0]
        else:
            trna_charging.gene_reaction_rule = ' and '.join(genes)
        

        trna_charging_reactions.append(trna_charging)
        charged_trna_metabolites.append(charged_trna_c)
    return trna_charging_reactions, charged_trna_metabolites
    

def trna_biogenesis(trna_info):
    '''trna_info is an object of class trna_information'''
    pretrna_transcription, pretrna_transcript_n, pretrna_mw = transcribe_pretrna(trna_info)
    trna_processing_reactions, trna_transcript_n = process_trna(trna_info, pretrna_transcript_n, pretrna_mw)
    # right now, no modification reaction and modified and non-modified are same metabolite
    trna_modifications_nuclear, modified_trna_transcript_n = modify_trna_nuclear(trna_info, trna_transcript_n) 
    trna_primary_export, trna_transcript_c = primary_export_trna(trna_info, modified_trna_transcript_n)
    trna_modifications_cytosolic, modified_trna_transcript_c, modified_trna_transcript_c_mw  = modify_trna_cytosolic(trna_info, trna_transcript_c)
    trna_degradation = degrade_trna(trna_info, modified_trna_transcript_c, modified_trna_transcript_c_mw)
    trna_charging_reactions, charged_trna_metabolites = charge_trna(trna_info, modified_trna_transcript_c, modified_trna_transcript_c_mw)
    #---------------------------------------------------------------------------------------------------
    
    reactions = [pretrna_transcription] + trna_processing_reactions + [trna_primary_export, trna_degradation] 
    reactions += trna_charging_reactions
    if trna_modifications_nuclear != None:
        reactions += [trna_modifications_nuclear]
    if trna_modifications_cytosolic != None:
        reactions += [trna_modifications_cytosolic]
    
    # reactiosn will be added to model, both charged and uncharged (modified) trna will be involved in translationr reactions
    return reactions, charged_trna_metabolites, modified_trna_transcript_c, modified_trna_transcript_c_mw
     


# # Consensus Sequences
# 
# positionally-independent (position won't effect .elements of cobra.Metabolite in cobra.Reaction)

# In[39]:


def get_base_frequency(seq_col, L):
    base_counts = {'T': 0, 'C': 0, 'G': 0, 'A': 0}
    for seq in trna_data[seq_col]:
        for base in base_counts.keys():
            base_counts[base] += seq.count(base)
    total = sum(base_counts.values())   
    base_frequencies = {base: (counts/total) for base,counts in base_counts.items()}
    base_frequencies['U'] = base_frequencies['T']
    base_frequencies.pop('T')

    final_seq = ''
    for base, frequency in base_frequencies.items():
        final_seq += base*round(frequency*L)

    # this assumes that if lengths are off, they are only off by one (accurate assumption for this dataset)
    if len(final_seq) != L:
        base_counts = {base: final_seq.count(base) for base in base_frequencies.keys()}
        if len(final_seq) > L:
            final_seq = ''
            to_remove = [k for k,v in base_frequencies.items() if v == min(base_frequencies.values())][0]
            base_counts[to_remove] -= 1
            for base in base_counts.keys():
                final_seq += base*base_counts[base]
        else:
            final_seq = ''
            to_add = [k for k,v in base_frequencies.items() if v == max(base_frequencies.values())][0]
            base_counts[to_add] += 1
            for base in base_counts.keys():
                final_seq += base*base_counts[base]

    return final_seq

trna_data = pd.read_excel(local_data_path + 'raw/trna_leaders_and_trailers.xlsx')
trna_data['Mature_Length'] = trna_data['mature seq'].apply(lambda x: len(x))

L_mature = trna_data.Mature_Length.value_counts()[trna_data.Mature_Length.value_counts() == trna_data.Mature_Length.value_counts().max()].index.tolist()[0]
L_leader = trna_data[' leader length'].value_counts()[trna_data[' leader length'].value_counts() == trna_data[' leader length'].value_counts().max()].index.tolist()[0]
L_trailer = trna_data[' trailer length'].value_counts()[trna_data[' trailer length'].value_counts() == trna_data[' trailer length'].value_counts().max()].index.tolist()[0]

mature_seq = get_base_frequency('mature seq', L_mature)
leader_seq, trailer_seq = get_base_frequency('leader seq', L_leader), get_base_frequency('trailer seq', L_trailer)

# CCA tail (position matters in the code)
base_counts = {base: mature_seq.count(base) for base in ['U', 'C', 'G', 'A']}
base_counts['C'] -= 2
base_counts['A'] -= 1
mature_seq = ''
for base in base_counts.keys():
    mature_seq += base*base_counts[base]
mature_seq += 'CCA'


# # Generate reactions

# In[46]:


trna_info = trna_information(maturetrna_sequence = mature_seq , id_ = 'generic', three_trailer_seq = trailer_seq, 
                             five_leader_seq = leader_seq, modifications = {},
                             intron_sequences = None)
trna_biogenesis_reactions, charged_trna_metabolites, modified_trna_transcript_c, modified_trna_transcript_c_mw  = trna_biogenesis(trna_info)


# In[24]:


# trna_model = cobra.Model('trna_biogenesis')
# trna_model.add_reactions(trna_biogenesis_reactions)
# cobra.io.save_json_model(trna_model, local_data_path + 'interim/trna_biogenesis.json')


# In[25]:


# trna_model


# In[26]:


# import escher
# from escher import Builder
# builder = escher.Builder(map_json=local_data_path + 'figures/trna_biogenesis_map.json')
# builder

