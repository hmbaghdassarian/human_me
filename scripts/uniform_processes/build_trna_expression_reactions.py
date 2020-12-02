#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings
import cobra

import pandas as pd

import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils.load_environmental_variables import raw_data_path
from utils import machinery as mach
from utils import parameters as params
from utils import metabolites as metab
# from utils import functions as func
from macromolecules.RNA import tRNA, RNA_fragment
from macromolecules.complex import add_biomass_change


# -To change to individual tRNA molecules rather than a generic one, start with the charge_trna function and also change the trna_biogenesis function
# 
# -To add modifications, will need to change the allowed_trna_modifications dictionary in utils, and edit the following functions: 1) modify_trna_nuclear, 2) modify_trna_cytosolic, 3) degrade_trna

# # TRNA Information Class

# In[2]:


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

# In[3]:


# def update_trna_degradation(trna_degradation_reaction):
#     trna_degradation_reaction.subsytem = 'tRNA_Biogenesis'
#     trna_degradation_reaction.gene_reaction_rule = ' and '.join(mach.lariat_machinery['Exosome'])
#     return trna_degradation_reaction


# In[4]:


class express_trna():
    def __init__(self, trna_info):
        self.trna_info = trna_info
        self.reactions = []
    def transcribe_pretrna(self):
        self.pretrna_n = tRNA(metabolite_name=self.trna_info.id + '_pre', seq=self.trna_info.pretrna_sequence, 
                       compartment = 'n', triphosphate = True)
        pretrna_transcription = self.pretrna_n.synthesize(id_ = 'TRANSCRIPTION_PRE_TRNA_' + self.trna_info.id)
        self.reactions.append(pretrna_transcription)
    def process_trna(self):
        '''

        This reaction processes pre-tRNA into mature tRNA in the nucleus. 
        This includes: CCA synthesis, 5' leader and 3' trailer cleavage (and degradation as separate reactions),
        and splicing (and intron degradation as a separate reactions for each intron).
        
        '''

        rxn = {self.pretrna_n: -1}
        # CCA synthesis
        rxn[metab.ntp_map_n['C']] = -2
        rxn[metab.ntp_map_n['A']] = -1
        rxn[metab.ppi_n] = 3
        trna_processing_machinery = mach.TRNT1.copy()

        # initialize
        reactions = list()
        rxn[metab.h2o_n] = 0 

        if self.trna_info.five_leader_seq != None: # if there is a 5' leader sequence
                five_leader = RNA_fragment(metabolite_name=self.trna_info.id, seq = self.trna_info.five_leader_seq, 
                    compartment = 'n', triphosphate = True, fragment_type = '5_leader')


                rxn[five_leader] = 1
                rxn[metab.h2o_n] -= 1 #endonuclolytic cleavage (RNAse P)
                trna_processing_machinery += mach.RNASEP

                five_leader_degradation = five_leader.exonucleolytic_degradation(reaction_name = self.trna_info.id + "_5'_leader_fragment_tRNA",
                                                                                update = True)
                self.reactions.append(five_leader_degradation)

                tp = False 
        else:
            tp = True

        self.trna_n = tRNA(metabolite_name=self.trna_info.id, seq = self.trna_info.maturetrna_sequence, 
                          compartment = 'n', triphosphate = tp)
        rxn[self.trna_n] = 1

        # 3' cleavage
        if self.trna_info.three_trailer_seq != None:
            three_trailer = RNA_fragment(metabolite_name = self.trna_info.id, seq = self.trna_info.three_trailer_seq, 
                                          fragment_type = '3_trailer',
                                         compartment = 'n',triphosphate = False)
            rxn[three_trailer] = 1
            rxn[metab.h2o_n] -= 1 #endonuclolytic cleavage (RNase Z)
            trna_processing_machinery += mach.RNASEZ

            three_trailer_degradation = three_trailer.exonucleolytic_degradation(reaction_name = self.trna_info.id + "_3'_trailer_fragment_tRNA")
            self.reactions.append(three_trailer_degradation)

        # splicing of intron
        if self.trna_info.intron_sequences != None:
            n_introns = len(self.trna_info.intron_sequences)
            trna_processing_machinery += mach.trna_splicing_machinery

            trna_introns_n = dict()
            for i in range(len(self.trna_info.intron_sequences)):
                trna_intron_n = RNA_fragment(metabolite_name = self.trna_info.id + '_' + str(i), 
                                            seq = self.trna_info.intron_sequences[i], fragment_type = 'trna_intron', 
                                            compartment = 'n', triphosphate = False)

                rxn[trna_intron_n] = 1

                intron_degradation = trna_intron_n.exonucleolytic_degradation(reaction_name = self.trna_info.id + "_intron_" + str(i) + '_tRNA', 
                                                                             update = True)
                self.reactions.append(intron_degradation)

            rxn[h2o_n] -= n_introns

        trna_processing = cobra.Reaction('PROCESSING_TRNA_' + self.trna_info.id)
        trna_processing.subsytem = 'tRNA_Biogenesis'
        trna_processing.add_metabolites(rxn)
        trna_processing.gene_reaction_rule = ' and '.join(trna_processing_machinery)

        if len(trna_processing.check_mass_balance()) > 0:
            raise ValueError('tRNA processing for ' + self.trna_info.id + ' is unbalanced')
        elif list(trna_processing.compartments) != ['n']:
            raise ValueError('tRNA processing must be confined to nuclear compartment')
        else:
            self.reactions.append(trna_processing)
    def modify_trna_nuclear(self):
        '''Add to the if statement in the future. This is for nuclear modifications.'''

        if len(self.trna_info.modifications) > 0:
            raise ValueError('tRNA modifications are not currently considered')
    #         modified_trna_transcript_n = trna_transcript_n.copy()
    #         modified_trna_transcript_n.id = trna_info.id + '_modified_trna[n]'
            ####
        else:
            trna_modifications_nuclear = None
            self.modified_trna_n = self.trna_n # not copy, same object
    def primary_export_trna(self):
        self.trna_c = self.modified_trna_n.change_compartment('c')

        trna_primary_export = cobra.Reaction(self.trna_info.id + 'PRIMARY_EXPORTtn')
        trna_primary_export.subsytem = 'tRNA_Biogenesis'
        trna_primary_export.name = 'trna nuclear export'

        export_rxn = {self.modified_trna_n: -1, self.trna_c: 1}
        # gtp hydrolysis on cytoplasmic side for export (see protein_expression nuclear_transport for details)
        export_rxn[metab.ntp_map_c['G']], export_rxn[metab.h2o_c], export_rxn[metab.ndp_map_c['G']], export_rxn[metab.pi_c], export_rxn[metab.h_c]  = -1, -1, 1, 1, 1
        trna_primary_export.add_metabolites(export_rxn)
        trna_primary_export.gene_reaction_rule = ' and '.join(mach.XPOT + mach.RAN)

        self.reactions.append(trna_primary_export)
    def modify_trna_cytosolic(self):
        '''Add to the if statement in the future. This is for cytosolic modifications'''
        if len(self.trna_info.modifications) > 0:
            raise ValueError('Modifications are not currently considered')
        else:
            trna_modifications_cytosolic = None
            self.modified_trna_c = self.trna_c
    def degrade_trna(self):
        # currently built on the assumption that there are no post-transcriptional modifications on trna 

        if self.trna_info.five_leader_seq != None: # if there is a 5' leader sequence
            tp = False 
        else:
            tp = True
        
        trna_degradation = self.modified_trna_c.exonucleolytic_degradation(reaction_name = self.trna_info.id + '_tRNA')
        trna_degradation.subsytem = 'tRNA_Biogenesis'
        trna_degradation.gene_reaction_rule = mach.XRN1[0]
        self.reactions.append(trna_degradation) 
    
    def charge_trna(self):
        '''tRNA charging reaction combines the activation and charging steps into one reaction.'''


        # in the future, this should take into account anticodon sequence, which should be in 
        # trna_info id

        # now, since just a generic charging reaction, will create one for each amino acid (for loop)

        # diagram: https://www.researchgate.net/figure/The-reaction-scheme-for-the-two-steps-of-aminoacylation-reaction-at-the-active-site-of_fig4_231225238
        trna_charging_reactions, charged_trna_metabolites = [], []
        for code, aa in metab.seq_amino_acid_map_c.items():
            elements = self.modified_trna_c.elements.copy()
            # attachment - loss of hydrogen from tRNA hydroxyl, and oxygen from amino acid carboxyl
            elements['H'] -= 1 
            elements['O'] -= 1
            for element, count in aa.elements.items():
                if element in elements.keys():
                    elements[element] += count
                else:
                    elements[element] = count

            charged_trna_c = tRNA(metabolite_name='charged_' + self.trna_info.id + '_' + code, 
                                 seq = '', compartment = 'c', triphosphate = self.modified_trna_c.triphosphate)
            charged_trna_c.compartment = 'c'
            charged_trna_c.elements = elements
            # +1 for loss of negative charge on oxygen of amino acid
            charged_trna_c.charge = self.modified_trna_c.charge + aa.charge + 1 

            trna_charging = cobra.Reaction('CHARGING_TRNA_' + self.trna_info.id + '_' + code)
            rxn = {self.modified_trna_c: -1, aa: -1, charged_trna_c: 1, metab.atp_c: -1, metab.ppi_c: 1, 
                   metab.amp_c: 1}
            trna_charging.add_metabolites(rxn)
            # add gprs
            genes = mach.seq_synthetase_map[code]
            if len(genes) == 1:
                trna_charging.gene_reaction_rule = genes[0]
            else:
                trna_charging.gene_reaction_rule = ' and '.join(genes)
            
            if len(trna_charging.check_mass_balance()) > 0:
                raise ValueError('tRNA charging for ' + self.trna_info.id + '_' + code + ' is unbalanced')
            elif list(trna_charging.compartments) != ['c']:
                raise ValueError('tRNA charging must be confined to cytosolic compartment')
            else:            
                trna_charging_reactions.append(trna_charging)
                charged_trna_metabolites.append(charged_trna_c)
        
        self.reactions += trna_charging_reactions
        self.charged_trna_metabolites = charged_trna_metabolites
    def add_biomass(self):
        for r in self.reactions:
            r.subsystem = 'tRNA_Biogenesis'
            add_biomass_change(r)


# In[5]:


def trna_biogenesis(trna_info):
    tb = express_trna(trna_info)
    tb.transcribe_pretrna()
    tb.process_trna()
    tb.modify_trna_nuclear()
    tb.primary_export_trna()
    tb.modify_trna_cytosolic()
    tb.degrade_trna()
    tb.charge_trna()
    tb.add_biomass()
    
    return tb.reactions, tb.charged_trna_metabolites, tb.modified_trna_c


# # Consensus Sequences
# 
# positionally-independent (position won't effect .elements of cobra.Metabolite in cobra.Reaction)

# In[6]:


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

trna_data = pd.read_excel(raw_data_path + 'trna_leaders_and_trailers.xlsx')
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

# In[7]:


trna_info = trna_information(maturetrna_sequence = mature_seq , id_ = 'generic', three_trailer_seq = trailer_seq, 
                             five_leader_seq = leader_seq, modifications = {},
                             intron_sequences = None)
trna_biogenesis_reactions, charged_trna_metabolites, modified_trna_transcript_c  = trna_biogenesis(trna_info)

