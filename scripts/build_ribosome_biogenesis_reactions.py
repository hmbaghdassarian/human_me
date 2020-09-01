#!/usr/bin/env python
# coding: utf-8

# In[1]:


# special for ribosome only, make sure to include code in building me model to not make ribosomes the same way 
# as other proteins


# In[1]:


import cobra

import pandas as pd
import numpy as np

from Bio.Seq import Seq
from Bio import SeqIO

import sys
sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *
from utils import *
from utils_2 import *

import build_mrna_expression_reactions as build_mrna
import build_protein_expression_reactions as build_protein


# # rRNA

# In[2]:


# rrna sequences
# assume the ncbi 45s is actually 47s...see notes for details
rrna_47s_seq = SeqIO.read(local_data_path + 'raw/45s_rrna_seq.txt', "fasta").seq.transcribe()
rrna_18s_seq = SeqIO.read(local_data_path + 'raw/18s_rrna_seq.txt', "fasta").seq.transcribe()
rrna_28s_seq = SeqIO.read(local_data_path + 'raw/28s_rrna_seq.txt', "fasta").seq.transcribe()
rrna_5_8s_seq = SeqIO.read(local_data_path + 'raw/5_8s_rrna_seq.txt', "fasta").seq.transcribe()
ets_5_seq = rrna_47s_seq[:rrna_47s_seq.index(rrna_18s_seq)]
its_1_seq = rrna_47s_seq[rrna_47s_seq.index(rrna_18s_seq) + len(rrna_18s_seq):rrna_47s_seq.index(rrna_5_8s_seq)]
its_2_seq = rrna_47s_seq[rrna_47s_seq.index(rrna_5_8s_seq) + len(rrna_5_8s_seq):rrna_47s_seq.index(rrna_28s_seq)]
ets_3_seq = rrna_47s_seq[rrna_47s_seq.index(rrna_28s_seq) + len(rrna_28s_seq):]

pre_rrna_5s_seq = SeqIO.read(local_data_path + 'raw/5s_rrna_seq.txt', "fasta").seq.transcribe()
rrna_5s_seq = pre_rrna_5s_seq[:120] # 120 is length of mature 5s_rrna



# rrna cut sites
# cut site indexes, relative to how far right (3' end is right) they are of a certain feature
# Fig. 3b https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4361047/
# scaled from length in figure to length of my sequence
A_prime_index = int(round(np.median([414,420])*len(ets_5_seq)/3657)) # location from 5' end of 47s
A_0_index = int(round(1642*len(ets_5_seq)/3657)) - A_prime_index # how far to the right of 45s is the A_0 site
site_2_index = int(round((6470-5527)*len(its_1_seq)/(6623-5527))) # how far right of end of 18s
site_4_index = int(round((7570-6779)*(len(its_2_seq)/(7935-6779)))) # how far to the right of the end of 5.8s/how far into ITS2 is the site 4 cut location
e_index = int(round((np.median([5606,5609])-5527)*len(its_1_seq)/(6623-5527)))# bp to right of end of 18s/start of its_1
conserved_stall_idx = int(round((np.median([6117,6192])-5527)*len(its_1_seq)/(6623-5527))) # how far right of end of 18s does RRP6 stall to form 21S-C
# Fig. 6b https://www.sciencedirect.com/science/article/pii/S1097276513005844?via%3Dihub
seven_s_idx = 190 
five_eight_plus_forty_idx = 40
six_s_index = 1 #https://www.nature.com/articles/s41594-019-0234-x?draft=collection

# # original: from 2+ difference soures
# # 420 and numerator from figure 3A https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3964915/
# A_prime_index = round(420*(len(ets_5_seq)/(1800 + 2000 + 420)))
# # Fig1D: https://www.researchgate.net/figure/Mapping-the-cleavages-in-human-ITS1-A-Alternative-processing-pathways-of-human-rRNA_fig1_235729322
# site_2_index = int(round(np.median([6396 -5520,6508-5520]) *(len(its_1_seq)/(6603-5520))))
# #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3964915/
# A_0_index = round(1800*(len(ets_5_frag2_seq)/(1800 + 2000)))  
# # supplementary Fig. 2 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3632142/
# conserved_stall_idx = int(round(np.median([590, 635])))
# # Fig1D: https://www.researchgate.net/figure/Mapping-the-cleavages-in-human-ITS1-A-Alternative-processing-pathways-of-human-rRNA_fig1_235729322
# site_4_index = int(round((7564-6773)*(len(its_2_seq)/(7891 - 6773))))
# e_index = 80 # fig 3c:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3632142/
# yeast_its2_length = 420 # https://www.microbiologyresearch.org/docserver/fulltext/jmm/66/2/126_jmm000426.pdf?expires=1594927728&id=id&accname=guest&checksum=7C6B2DF6CE3C3080E28605D15B99DF1E
# yeast_c2 = 140 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4361047/
# seven_s_idx = int(round((yeast_c2/yeast_its2_length)*len(its_2_seq))) # how far to the right of 5.8s does sequence extend to form the 7s rrrna



# In[3]:


psim_rib = psim_me.copy()
def format_location(x):
    return ['n', 'c']
psim_rib.LOCATION = psim_rib.LOCATION.apply(lambda x: format_location(x))


# In[4]:


def build_ribosome_protein_expression_reactions():
    '''Reactions associated with transcription and translation of ribosomal proteins'''
    
    # RPS27A and RPL40 ubiquitin-fusion in the future
    
    rs_ids = rs['HGNC ID (gene)'].tolist()
    rs_expression_reactions, rs_protein_metabolites = list(), list()
    for i in rs_ids:
        gene_info = generate_geneinfo_object(hgnc_id = i, psim = psim_rib, 
                    machinery_list = list(), metabolic_model = cobra.Model())
        gene_info.final_locations = {'c': 'Cytosolic Tranport', 'n': 'Cytosolic Tranport'}
        mrna_expression_reactions, mrna_transcript_c = build_mrna.mrna_expression(gene_info)
        protein_expression_reactions, protein_metabolites = build_protein.get_protein_expression_reactions(gene_info)
        protein_expression_reactions = protein_expression_reactions[:-3] # no nuclear degradation
        rs_expression_reactions += mrna_expression_reactions + protein_expression_reactions
        rs_protein_metabolites += protein_metabolites

    rl_ids = rl['HGNC ID (gene)'].tolist()
    RPL40_HGNC = 'HGNC:12458'
    rl_ids.remove(RPL40_HGNC) # RPL40 is a ubiquitin monomer
    rl_expression_reactions, rl_protein_metabolites = list(), list()
    for i in rl_ids:
        gene_info = generate_geneinfo_object(hgnc_id = i, psim = psim_rib, 
                    machinery_list = list(), metabolic_model = cobra.Model())
        gene_info.final_locations = {'c': 'Cytosolic Tranport', 'n': 'Cytosolic Tranport'}
        mrna_expression_reactions, mrna_transcript_c = build_mrna.mrna_expression(gene_info)
        protein_expression_reactions, protein_metabolites = build_protein.get_protein_expression_reactions(gene_info)
        protein_expression_reactions = protein_expression_reactions[:-3] # no nuclear degradation
        rl_expression_reactions += mrna_expression_reactions + protein_expression_reactions
        rl_protein_metabolites += protein_metabolites

    # RPL40-UB FUSION----------------------------------------------------------------------
    gene_info = generate_geneinfo_object(hgnc_id = RPL40_HGNC, psim = psim_rib, 
                    machinery_list = list(), metabolic_model = cobra.Model())
    gene_info.final_locations = {'n': 'Cytosolic Tranport'}
    mrna_expression_reactions, mrna_transcript_c = build_mrna.mrna_expression(gene_info)
    translation_elongation_c, unfolded_protein_c = build_protein.translate_protein_cytosolic(gene_info)
    protein_expression_reactions.append(translation_elongation_c)

    # cleaved protein sequence, gene_info object, and cobra.Metabolite
    processed_seq = gene_info.protein_seq[:gene_info.protein_seq.index(build_protein.single_ubiquitin_sequence)] + gene_info.protein_seq[gene_info.protein_seq.index(build_protein.single_ubiquitin_sequence) + len(build_protein.single_ubiquitin_sequence):]
    psim_temp = psim_rib.copy()
    psim_temp.loc[psim_temp[psim_temp.HGNC_ID == RPL40_HGNC].index, 'PROTEIN_SEQ'] = processed_seq
    gene_info = generate_geneinfo_object(hgnc_id = RPL40_HGNC, psim = psim_temp, 
                    machinery_list = list(), metabolic_model = cobra.Model())
    
    gene_info.final_locations = {'n': 'Cytosolic Tranport'}

    processed_unfolded_protein_c = make_protein_metabolite(id_ = RPL40_HGNC + '_processed_unfolded',
                                    amino_acid_counts = gene_info.amino_acid_counts, 
                                                           L_protein = len(processed_seq), compartment = 'c')
    ub_cleavage = cobra.Reaction(gene_info.hgnc_id + '_UBIQUITIN_CLEAVAGEc')
    ub_cleavage.subsytem = 'Protein_Expression'
    ub_cleavage.add_metabolites({unfolded_protein_c:-1, h2o_c: -1, 
                                 build_protein.ub_c: 1, processed_unfolded_protein_c: 1})
    ub_cleavage.gene_reaction_rule = UCHL3[0]

    protein_folding_cytosolic, folded_protein_c = build_protein.fold_protein_cytosolic(gene_info, 
                                                                                       processed_unfolded_protein_c)
    nuclear_import, folded_protein_n = build_protein.transport_nuclear_protein(gene_info, folded_protein_c)

    rl_expression_reactions += mrna_expression_reactions + [translation_elongation_c, ub_cleavage, protein_folding_cytosolic, nuclear_import] + build_protein.degrade_cytosolic_protein(gene_info, folded_protein_c)
    rl_protein_metabolites += [folded_protein_c, folded_protein_n]
    
    return rs_expression_reactions, rs_protein_metabolites, rl_expression_reactions, rl_protein_metabolites


# In[6]:


def update_rrna_degradation(rrna_degradation_reaction, nucleus = True):
    rrna_degradation_reaction.subsytem = 'Ribosome_Biogenesis'
    if nucleus:
        rule_part2 = ' and '.join(lariat_machinery['Exosome'] + lariat_machinery['NEXT Complex']) + ')'
        rrna_degradation_reaction.gene_reaction_rule = lariat_machinery["5' Degradation"][0] + ' or (' + rule_part2
    else:
        rrna_degradation_reaction.gene_reaction_rule = ' and '.join(exosome['HGNC ID (gene)'].tolist())
    return rrna_degradation_reaction

def build_rrna5s_reactions(rpl5_n, rpl11_n):
    
    # TRANSCRIPTION - basically emulates Transcript.transcript_elongation reaction
    pre_rrna5s_n, pre_base_counts5s = make_rna_metabolite('pre_5s', pre_rrna_5s_seq, molecule_type = 'rrna', compartment = 'n')

    rrna5s_transcription = cobra.Reaction('TRANSCRIPTION_PRE_RRNA5s')
    rrna5s_transcription.subsytem = 'Ribosome_Biogenesis'
    rxn = dict()
    for ntp, base_letter in seq_metabolite_map.items():
        rxn[ntp] = -1*pre_base_counts5s[base_letter]
    rxn[ppi_n] = len(pre_rrna_5s_seq) - 1
    rxn[pre_rrna5s_n] = 1
    rrna5s_transcription.add_metabolites(rxn)
    rrna5s_transcription.gene_reaction_rule = ' and '.join(rnap3_transcription_machinery) 
    
    # PROCESSING - mature rrna (3->5' exonucleolytic cleave of last 24 bases) and complex formation with RPL5/RPL11
    rrna5s_processing = cobra.Reaction('PROCESSING_RRNA5s')
    rrna5s_processing.subsytem = 'Ribosome_Biogenesis'
    rrna5s_n, base_counts5s = make_rna_metabolite('5s', rrna_5s_seq, molecule_type = 'rrna', compartment = 'n') 
    deg_base_counts = dict()
    for k,v in pre_base_counts5s.items():
        deg_base_counts[k] = v - base_counts5s[k]

    rxn = dict()
    rxn[h2o_n] = -sum(deg_base_counts.values()) # no -1 because all bonds 5'-most bond cleave
    rxn[h_n] = sum(deg_base_counts.values())
    for k,v in nmp_map_n.items():
        rxn[v] = deg_base_counts[k]
    
    metabolites = [rrna5s_n, rpl5_n, rpl11_n]
    complex_info = {'METABOLITES': metabolites, 'IDS': [m.id.split('_')[0] for m in metabolites], 
                                   'METABOLITE_TYPES': [m.id.split('_')[-1].split('[')[0] for m in metabolites]}
    rrna5s_complex_n, rrna5s_complex_n_id  = make_complex_metabolite(**complex_info)


    rxn[pre_rrna5s_n], rxn[rpl5_n], rxn[rpl11_n]  = -1, -1, -1
    rxn[rrna5s_complex_n] = 1

    rrna5s_processing.add_metabolites(rxn)
    rrna5s_processing.gene_reaction_rule = REXO5

    
    # TRANSPORT - will be transported as pre60s later, but make an rrna5s cytoplasmic for degradation, as
    # ribosome dissociates in cytoplasm
    # must add nucleocytoplasmic export via ran gtp: https://www.sciencedirect.com/science/article/pii/S0171933504702575?via%3Dihub
    rrna5s_c = rrna5s_n.copy()
    rrna5s_c.id = rrna5s_c.id.replace('[n]', '[c]')
    rrna5s_c.compartment = 'c'
    
    # Degradation
    rrna5s_degradation = rna_exonucleolytic_degradation(rrna5s_c, base_counts5s, rrna_5s_seq, reaction_name = '5s_rRNA')
    rrna5s_degradation.subsytem = 'Ribosome_Biogenesis'
    rrna5s_degradation.gene_reaction_rule = ' and '.join(exosome['HGNC ID (gene)'].tolist())
    
    rrna5s_reactions = [rrna5s_transcription, rrna5s_processing, rrna5s_degradation]
    
    
    return rrna5s_reactions, rrna5s_complex_n, rrna5s_c


# In[7]:


# ets_5_frag1 is from 5' end of 47s to A' site
# ets_5_frag2 is from A' to 18s
# ets_5_frag3 is from A' to A0
# ets_5_frag4 is from A0 site to site 1 (start of 18s)
# its_1_frag1_seq is between site E and teh conserved stall location of RRP6
# its_1_frag2_seq is less than E site (some degradation) + a polyA/U tail

def build_other_rrna_reactions(rrna5s_complex_n, rs_protein_metabolites, rl_protein_metabolites, 
                              rpl5_n, rpl11_n, rrna5s_c):
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6315592/ figure 2
    
    # 47s transcription------------------------------------------------------------------------------------
    rrna_47s_n, base_counts_47s = make_rna_metabolite('47s', rrna_47s_seq, molecule_type = 'rrna', compartment = 'n')
    rrna_47s_transcription = cobra.Reaction('TRANSCRIPTION_RRNA_47s')
    rrna_47s_transcription.subsytem = 'Ribosome_Biogenesis'
    rxn = dict()
    for ntp, base_letter in seq_metabolite_map.items():
        rxn[ntp] = -1*base_counts_47s[base_letter]
    rxn[ppi_n] = len(rrna_47s_seq) - 1
    rxn[rrna_47s_n] = 1
    rrna_47s_transcription.add_metabolites(rxn)
    rrna_47s_transcription.gene_reaction_rule = ' and '.join(rnap1['HGNC ID (gene)'].tolist() + rnap1_tfs)  
    
    # 45s formation------------------------------------------------------------------------------------
    ets_5_frag1_seq = ets_5_seq[:A_prime_index] 
    ets_5_frag2_seq = ets_5_seq[A_prime_index:]
    rrna_45s_seq = rrna_47s_seq[A_prime_index:rrna_47s_seq.index(rrna_28s_seq) + len(rrna_28s_seq)]

    ets_5_frag1_n, base_counts_ets_5_frag1 = make_rna_metabolite('ets_5_frag1', ets_5_frag1_seq, 
                                                                 molecule_type = 'rrna', compartment = 'n')
    rrna_45s_n, base_counts_rrna_45s = make_rna_metabolite('45s', rrna_45s_seq, molecule_type = 'rrna', 
                                                           compartment = 'n', triphosphate=False)
    ets_3_n, base_counts_ets_3 = make_rna_metabolite('ets_3', ets_3_seq, molecule_type = 'rrna', 
                                                     compartment = 'n', triphosphate = False)

    rrna_45s_formation = cobra.Reaction('FORMATION_RRNA_45s')
    rrna_45s_formation.subsytem = 'Ribosome_Biogenesis'
    rxn = dict()
    rxn[h2o_n] = -2 # 2 endonuclolytic cleavage events to go from 47s to 45s
    rxn[rrna_47s_n], rxn[rrna_45s_n], rxn[ets_3_n], rxn[ets_5_frag1_n] = -1, 1, 1, 1
    rrna_45s_formation.add_metabolites(rxn)
    rrna_45s_formation.gene_reaction_rule = ' and '.join(UTP10 + RNASEN)
    
    ets_3_degradation = rna_exonucleolytic_degradation(ets_3_n, base_counts_ets_3, ets_3_seq, reaction_name = 'ets_3_rRNA', triphosphate = False)
    ets_3_degradation = update_rrna_degradation(ets_3_degradation)
    ets_5_frag1_degradation = rna_exonucleolytic_degradation(ets_5_frag1_n, base_counts_ets_5_frag1, ets_5_frag1_seq, reaction_name = 'ets_5_frag1_rRNA', triphosphate = True)
    ets_5_frag1_degradation = update_rrna_degradation(ets_5_frag1_degradation)
    
    #45S-->30S + 32.5S------------------------------------------------------------------------------------
    idx_30s = rrna_45s_seq.index(rrna_18s_seq) + len(rrna_18s_seq)
    rrna_30s_seq = ets_5_frag2_seq + rrna_18s_seq + rrna_45s_seq[idx_30s:idx_30s + site_2_index]
    rrna_32_5s_seq = rrna_45s_seq[idx_30s + site_2_index:]

    rrna_30s_n, base_counts_rrna_30s = make_rna_metabolite('30s', rrna_30s_seq, molecule_type = 'rrna', 
                                                           compartment = 'n', triphosphate=False)
    rrna_32_5s_n, base_counts_rrna_32_5s = make_rna_metabolite('32_5s', rrna_32_5s_seq, molecule_type = 'rrna', 
                                                               compartment = 'n', triphosphate=False)

    rrna_30s_formation = cobra.Reaction('FORMATION_RRNA_30s_32_5s')
    rrna_30s_formation.subsytem = 'Ribosome_Biogenesis'
    rxn = dict()
    rxn[h2o_n] = -1 # endonuclolytic cleavage event at site 2
    rxn[rrna_45s_n], rxn[rrna_30s_n], rxn[rrna_32_5s_n] = -1, 1, 1
    rrna_30s_formation.add_metabolites(rxn)
    rrna_30s_formation.gene_reaction_rule = RMRP[0]
    
    #26s formation------------------------------------------------------------------------------------
    rrna_26s_seq = ets_5_frag2_seq[A_0_index:] + rrna_18s_seq + its_1_seq[:site_2_index]
    ets_5_frag3_seq = ets_5_frag2_seq[:A_0_index]
    rrna_26s_n, base_counts_rrna_26s = make_rna_metabolite('26s', rrna_26s_seq, molecule_type = 'rrna', 
                                                           compartment = 'n', triphosphate=False)
    ets_5_frag3_n, base_counts_ets_5_frag3 = make_rna_metabolite('ets_5_frag3', ets_5_frag3_seq, 
                                                                 molecule_type = 'rrna', compartment = 'n', 
                                                                 triphosphate=False)

    rrna_26s_formation = cobra.Reaction('FORMATION_RRNA_26s')
    rrna_26s_formation.subsytem = 'Ribosome_Biogenesis'
    rxn = dict()
    rxn[h2o_n] = -1 # endonuclolytic cleavage event at site 2
    rxn[rrna_30s_n], rxn[rrna_26s_n], rxn[ets_5_frag3_n] = -1, 1, 1
    rrna_26s_formation.add_metabolites(rxn)
    rrna_26s_formation.gene_reaction_rule = UTP23[0]

    ets_5_frag3_degradation = rna_exonucleolytic_degradation(ets_5_frag3_n, base_counts_ets_5_frag3, ets_5_frag3_seq, reaction_name = 'ets_5_frag3_rRNA', triphosphate = False)
    ets_5_frag3_degradation = update_rrna_degradation(ets_5_frag3_degradation)
    
    # 21S formation------------------------------------------------------------------------------------
    rrna_21s_seq = rrna_18s_seq + its_1_seq[:site_2_index]
    ets_5_frag4_seq = ets_5_frag2_seq[A_0_index:]
    rrna_21s_n, base_counts_rrna_21s = make_rna_metabolite('21s', rrna_21s_seq, molecule_type = 'rrna', compartment = 'n',
                                                           triphosphate=False)
    ets_5_frag4_n, base_counts_ets_5_frag4 = make_rna_metabolite('ets_5_frag4', ets_5_frag4_seq, molecule_type = 'rrna', compartment = 'n',
                                                                 triphosphate=False)

    rrna_21s_formation = cobra.Reaction('FORMATION_RRNA_21s')
    rrna_21s_formation.subsytem = 'Ribosome_Biogenesis'
    rxn = dict()
    rxn[h2o_n] = -1 # endonuclolytic cleavage event at site 1
    rxn[rrna_26s_n], rxn[rrna_21s_n], rxn[ets_5_frag4_n] = -1, 1, 1
    rrna_21s_formation.add_metabolites(rxn)
    rrna_21s_formation.gene_reaction_rule = UTP24[0]
    ets_5_frag4_degradation = rna_exonucleolytic_degradation(ets_5_frag4_n, base_counts_ets_5_frag4, ets_5_frag4_seq, reaction_name = 'ets_5_frag4_rRNA', triphosphate = False)
    ets_5_frag4_degradation = update_rrna_degradation(ets_5_frag4_degradation)
    # 21SC formation------------------------------------------------------------------------------------

    rrna_21sc_seq = rrna_18s_seq + its_1_seq[:conserved_stall_idx]
    rrna_21sc_n, base_counts_rrna_21sc = make_rna_metabolite('21sc', rrna_21sc_seq, molecule_type = 'rrna', 
                                                             compartment = 'n', triphosphate=False)

    deg_seq = its_1_seq[conserved_stall_idx: site_2_index]
    base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

    rrna_21sc_formation = cobra.Reaction('FORMATION_RRNA_21sc')
    rrna_21sc_formation.subsytem = 'Ribosome_Biogenesis'
    rxn = dict()
    rxn[rrna_21s_n], rxn[rrna_21sc_n] = -1,1
    # exonucleolytic cleavage
    for k,v in nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[h2o_n] = -len(deg_seq)
    rxn[h_n] = len(deg_seq)
    rrna_21sc_formation.add_metabolites(rxn)
    rrna_21sc_formation.gene_reaction_rule = exosome[exosome['Approved symbol'] == 'EXOSC10']['HGNC ID (gene)'].tolist()[0]

    # 18se formation------------------------------------------------------------------------------------
    rrna_18se_seq = rrna_18s_seq + its_1_seq[:e_index]
    its_1_frag1_seq = its_1_seq[e_index:conserved_stall_idx]
    rrna_18se_n, base_counts_rrna_18se = make_rna_metabolite('18se', rrna_18se_seq, molecule_type = 'rrna', 
                                                             compartment = 'n', triphosphate=False)
    its_1_frag1_n, base_counts_its_1_frag1 = make_rna_metabolite('its_1_frag1', its_1_frag1_seq, molecule_type = 'rrna', compartment = 'n', triphosphate=False)

    rrna_18se_formation = cobra.Reaction('FORMATION_RRNA_18se')
    rrna_18se_formation.subsytem = 'Ribosome_Biogenesis'
    # endonuclolytic cleavage event at site E
    rrna_18se_formation.add_metabolites({h2o_n: -1, rrna_21sc_n: -1, rrna_18se_n: 1, its_1_frag1_n: 1})
    rrna_18se_formation.gene_reaction_rule = UTP24[0]
    its_1_frag1_degradation = rna_exonucleolytic_degradation(its_1_frag1_n, base_counts_its_1_frag1, its_1_frag1_seq, reaction_name = 'its_1_frag1_rRNA', triphosphate = False)
    its_1_frag1_degradation = update_rrna_degradation(its_1_frag1_degradation)
    
    # 18se nuclear processing------------------------------------------------------------------------------------
    rrna_18se_processed_seq = rrna_18se_seq[:-int(0.75*e_index)] # degradation of 60/80 bps of ITS1 by PARN
    rrna_18se_processed_seq += 'U'*int(0.125*e_index)+'A'*int(0.125*e_index) # polyU by PAPD5
    deg_seq = rrna_18se_seq[-int(0.75*e_index):]

    rrna_18se_processed_n, base_counts_rrna_18se_processed = make_rna_metabolite('18se_processed', rrna_18se_processed_seq, molecule_type = 'rrna', compartment = 'n', triphosphate=False)
    base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

    rrna_18se_processing = cobra.Reaction('PROCESSING_RRNA_18se')
    rrna_18se_processing.subsytem = 'Ribosome_Biogenesis'
    rxn = dict()
    rxn[rrna_18se_n], rxn[rrna_18se_processed_n] = -1,1
    # exonucleolytic cleavage
    for k,v in nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[h2o_n] = -len(deg_seq)
    rxn[h_n] = len(deg_seq)
    # polyU/A tail synthesis<--unsure why ppi_n has -1 but it mass balances
    rxn[atp_n], rxn[ntp_map_n['U']], rxn[ppi_n] = -int(0.125*e_index), -int(0.125*e_index), int(0.25*e_index)-1

    rrna_18se_processing.add_metabolites(rxn)
    rrna_18se_processing.gene_reaction_rule = ' and '.join(PARN + PAPD5)
    
    # pre40s complex------------------------------------------------------------------------------------
    metabolites = [m for m in rs_protein_metabolites if m.compartment == 'n'] + [rrna_18se_processed_n]
    complex_info = {'METABOLITES': metabolites, 'IDS': [m.id.split('_')[0] for m in metabolites], 
                                   'METABOLITE_TYPES': [m.id.split('_')[-1].split('[')[0] for m in metabolites]}
    pre40s_complex_formation, pre40s_complex_n = form_complex(complex_id = 'pre40s', **complex_info)
    pre40s_complex_formation.id = 'pre40s_COMPLEX_FORMATIONn'
    pre40s_complex_formation.lower_bound = 0
    pre40s_complex_formation.gene_reaction_rule = ' and '.join(pre40s_rbfs)

    # pre40s nucleocytoplasmic export-----------------------------------------------------------------------
    pre40s_complex_c = pre40s_complex_n.copy()
    pre40s_complex_c.id = pre40s_complex_c.id.replace('[n]', '[c]')
    pre40s_complex_c.compartment = 'c'

    pre40s_transport = cobra.Reaction('pre40s_NUCLEAR_EXPORTtn')
    pre40s_transport.subsytem = 'Ribosome_Biogenesis'
    pre40s_transport.name = 'pre40s nuclear export'
    rxn = {pre40s_complex_n: -1, pre40s_complex_c: 1}
    rxn[ntp_map_c['G']], rxn[h2o_c], rxn[ndp_map_c['G']], rxn[pi_c], rxn[h_c]  = -1, -1, 1, 1, 1
    pre40s_transport.add_metabolites(rxn)
    pre40s_transport.gene_reaction_rule = ' and '.join(tfiiia + RAN + XPO1)
    
    # 18s/mature 40s formation------------------------------------------------------------------------------------
    its_1_frag2_seq = its_1_seq[:int(0.25*e_index)+1] + 'U'*int(0.125*e_index)+'A'*int(0.125*e_index)
    rrna_18s_c, base_counts_rrna_18s = make_rna_metabolite('18s', rrna_18s_seq, molecule_type = 'rrna', 
                                                           compartment = 'c', triphosphate=False)
    its_1_frag2_c, base_counts_its_1_frag2 = make_rna_metabolite('its_1_frag2', its_1_frag2_seq, 
                                                                 molecule_type = 'rrna', compartment='c',
                                                                 triphosphate=False)

    rrna_18s_formation = cobra.Reaction('40s_MATURATION')
    rrna_18s_formation.subsytem = 'Ribosome_Biogenesis'
    # endonuclolytic cleavage event at site 3
    metabolites = [m for m in rs_protein_metabolites if m.compartment == 'c'] + [rrna_18s_c]
    complex_info = {'METABOLITES': metabolites, 'IDS': [m.id.split('_')[0] for m in metabolites], 
                                   'METABOLITE_TYPES': [m.id.split('_')[-1].split('[')[0] for m in metabolites]}
    forty_s_complex_c, forty_s_id = make_complex_metabolite(complex_id = '40s', **complex_info)
    rrna_18s_formation.add_metabolites({h2o_n: -1, pre40s_complex_c: -1, forty_s_complex_c: 1, its_1_frag2_c: 1})
    rrna_18s_formation.gene_reaction_rule = NOB1[0]

    its_1_frag2_degradation = rna_exonucleolytic_degradation(its_1_frag2_c, base_counts_its_1_frag2, its_1_frag2_seq, reaction_name = 'its_1_frag2_rRNA', triphosphate = False, nucleus = False)
    its_1_frag2_degradation = update_rrna_degradation(its_1_frag2_degradation, nucleus = False)
    
    
    #18s degradation------------------------------------------------------------------------------------
    rrna_18s_degradation = rna_exonucleolytic_degradation(rrna_18s_c, base_counts_rrna_18s, rrna_18s_seq, 
                                                     reaction_name = '18s_rrna_degradation', 
                                                     triphosphate = False)
    rrna_18s_degradation = update_rrna_degradation(rrna_18s_degradation, nucleus = False)
    
    
    # 32S formation------------------------------------------------------------------------------------
    deg_seq = its_1_seq[site_2_index:]
    rrna_32s_seq = rrna_32_5s_seq[len(deg_seq):]

    rrna_32s_n, base_counts_rrna_32s = make_rna_metabolite('32s', rrna_32s_seq, molecule_type = 'rrna', 
                                                           compartment = 'n',triphosphate=False)
    base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

    rrna_32s_formation = cobra.Reaction('FORMATION_RRNA_32s')
    rrna_32s_formation.subsytem = 'Ribosome_Biogenesis'

    rxn = dict()
    rxn[rrna_32_5s_n], rxn[rrna_32s_n] = -1,1
    # exonucleolytic cleavage
    for k,v in nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[h2o_n] = -len(deg_seq)
    rxn[h_n] = len(deg_seq)

    rrna_32s_formation.add_metabolites(rxn)
    rrna_32s_formation.gene_reaction_rule = lariat_machinery["5' Degradation"][0]
    
    #32s-->12s + 28.5s------------------------------------------------------------------------------------
 
    rrna_12s_seq = rrna_5_8s_seq + its_2_seq[:site_4_index]
    rrna_28_5s_seq = its_2_seq[site_4_index:] + rrna_28s_seq
    rrna_12s_n, base_counts_rrna_12s = make_rna_metabolite('12s', rrna_12s_seq, molecule_type = 'rrna', compartment = 'n', 
                                                           triphosphate=False)
    rrna_28_5s_n, base_counts_rrna_28_5s = make_rna_metabolite('28_5s', rrna_28_5s_seq, molecule_type = 'rrna', compartment = 'n', 
                                                               triphosphate=False)

    rrna_12s_28_5s_formation = cobra.Reaction('FORMATION_RRNA_12s_28_5s')
    rrna_12s_28_5s_formation.subsytem = 'Ribosome_Biogenesis'
    rxn = dict()
    rxn[h2o_n] = -1 # endonuclolytic cleavage event at site 4
    rxn[rrna_32s_n], rxn[rrna_12s_n], rxn[rrna_28_5s_n] = -1,1,1
    rrna_12s_28_5s_formation.add_metabolites(rxn)
    rrna_12s_28_5s_formation.gene_reaction_rule = LAS1[0]
    
    #28s formation------------------------------------------------------------------------------------
    deg_seq = rrna_28_5s_seq[:rrna_28_5s_seq.index(rrna_28s_seq)]

    rrna_28s_n, base_counts_rrna_28s = make_rna_metabolite('28s', rrna_28s_seq, molecule_type = 'rrna', compartment = 'n', 
                                                           triphosphate=False)
    base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

    rrna_28s_formation = cobra.Reaction('FORMATION_RRNA_28s')
    rrna_28s_formation.subsytem = 'Ribosome_Biogenesis'

    rxn = dict()
    rxn[rrna_28_5s_n], rxn[rrna_28s_n] = -1,1
    # exonucleolytic cleavage
    for k,v in nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[h2o_n] = -len(deg_seq)
    rxn[h_n] = len(deg_seq)

    rrna_28s_formation.add_metabolites(rxn)
    rrna_28s_formation.gene_reaction_rule = lariat_machinery["5' Degradation"][0]
    
    #7s formation------------------------------------------------------------------------------------
    deg_seq = its_2_seq[seven_s_idx: site_4_index]
    rrna_7s_seq = rrna_5_8s_seq + its_2_seq[:seven_s_idx]
    rrna_7s_n, base_counts_rrna_7s = make_rna_metabolite('7s', rrna_7s_seq, molecule_type = 'rrna', compartment = 'n', triphosphate=False)
    base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

    rrna_7s_formation = cobra.Reaction('FORMATION_RRNA_7s')
    rrna_7s_formation.subsytem = 'Ribosome_Biogenesis'

    rxn = dict()
    rxn[rrna_12s_n], rxn[rrna_7s_n] = -1,1
    # exonucleolytic cleavage
    for k,v in nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[h2o_n] = -len(deg_seq)
    rxn[h_n] = len(deg_seq)

    rrna_7s_formation.add_metabolites(rxn)
    rrna_7s_formation.gene_reaction_rule = ' and '.join(DIS3 + ISG20L2)
    
    #5.8s+40 formation------------------------------------------------------------------------------------
    deg_seq = its_2_seq[five_eight_plus_forty_idx: seven_s_idx]
    rrna_5_8s_plus_40_seq = rrna_5_8s_seq + its_2_seq[:five_eight_plus_forty_idx]
    rrna_5_8s_plus_40_n, base_counts_rrna_5_8s_plus_40 = make_rna_metabolite('5_8s_plus_40', rrna_5_8s_plus_40_seq, molecule_type = 'rrna', compartment = 'n', triphosphate=False)
    base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

    rrna_5_8s_plus_40_formation = cobra.Reaction('FORMATION_RRNA_5_8s_plus_40')
    rrna_5_8s_plus_40_formation.subsytem = 'Ribosome_Biogenesis'

    rxn = dict()
    rxn[rrna_7s_n], rxn[rrna_5_8s_plus_40_n] = -1,1
    # exonucleolytic cleavage
    for k,v in nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[h2o_n] = -len(deg_seq)
    rxn[h_n] = len(deg_seq)

    rrna_5_8s_plus_40_formation.add_metabolites(rxn)
    rrna_5_8s_plus_40_formation.gene_reaction_rule = ' and '.join(DIS3 + ISG20L2)
    
    #6s formation------------------------------------------------------------------------------------

    deg_seq = its_2_seq[six_s_index: five_eight_plus_forty_idx]
    rrna_6s_seq = rrna_5_8s_seq + its_2_seq[:six_s_index]
    rrna_6s_n, base_counts_rrna_6s = make_rna_metabolite('6s', rrna_6s_seq, molecule_type = 'rrna', compartment = 'n', 
                                                         triphosphate=False)
    base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

    rrna_6s_formation = cobra.Reaction('FORMATION_RRNA_6s')
    rrna_6s_formation.subsytem = 'Ribosome_Biogenesis'

    rxn = dict()
    rxn[rrna_5_8s_plus_40_n], rxn[rrna_6s_n] = -1,1
    # exonucleolytic cleavage
    for k,v in nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[h2o_n] = -len(deg_seq)
    rxn[h_n] = len(deg_seq)

    rrna_6s_formation.add_metabolites(rxn)
    rrna_6s_formation.gene_reaction_rule = exosome[exosome['Approved symbol'] == 'EXOSC10']['HGNC ID (gene)'].tolist()[0]

    # pre60s complex formation------------------------------------------------------------------------------------
    rl_2 = list(set(rl_protein_metabolites).difference([rpl5_n, rpl11_n]))
    metabolites = [m for m in rl_2 if m.compartment == 'n'] + [rrna_28s_n, rrna_6s_n, rrna5s_complex_n]
    complex_info = {'METABOLITES': metabolites, 'IDS': [m.id.split('_')[0] for m in metabolites], 
                                   'METABOLITE_TYPES': [m.id.split('_')[-1].split('[')[0] for m in metabolites]}
    pre60s_complex_formation, pre60s_complex_n = form_complex(complex_id = 'pre60s', **complex_info)
    # add a gtp hydrolysis to the complex formation: https://www.embopress.org/doi/full/10.15252/embj.2018100278
    rxn = pre60s_complex_formation.metabolites.copy()
    rxn[ntp_map_n['G']], rxn[h2o_n], rxn[gdp_n], rxn[pi_n], rxn[h_n]  = -1, -1, 1, 1, 1
    pre60s_complex_formation.add_metabolites(rxn)

    pre60s_complex_formation.id = 'pre60s_COMPLEX_FORMATIONn'
    pre60s_complex_formation.lower_bound = 0
    pre60s_complex_formation.gene_reaction_rule = ' and '.join(pre60s_rbfs)
    
    # pre60s nucleocytoplasmic export-----------------------------------------------------------------------
    pre60s_complex_c = pre60s_complex_n.copy()
    pre60s_complex_c.id = pre60s_complex_c.id.replace('[n]', '[c]')
    pre60s_complex_c.compartment = 'c'

    pre60s_transport = cobra.Reaction('pre60s_NUCLEAR_EXPORTtn')
    pre60s_transport.subsytem = 'Ribosome_Biogenesis'
    pre60s_transport.name = 'pre60s nuclear export'
    rxn = {pre60s_complex_n: -1, pre60s_complex_c: 1}
    rxn[ntp_map_c['G']], rxn[h2o_c], rxn[ndp_map_c['G']], rxn[pi_c], rxn[h_c]  = -1, -1, 1, 1, 1
    pre60s_transport.add_metabolites(rxn)
    pre60s_transport.gene_reaction_rule = ' and '.join(tfiiia + RAN + XPO1)
    
    
    #5.8s/mature 60s formation------------------------------------------------------------------------------------
    deg_seq = its_2_seq[:six_s_index]
    rrna_5_8s_c, base_counts_rrna_5_8s = make_rna_metabolite('5_8s', rrna_5_8s_seq, molecule_type = 'rrna', compartment = 'c', 
                                                             triphosphate=False)
    base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)
    
    rrna_28s_c = rrna_28s_n.copy()
    rrna_28s_c.id = rrna_28s_c.id.replace('[n]', '[c]')
    rrna_28s_c.compartment = 'c'
    
    metabolites = [m for m in rl_2 if m.compartment == 'c'] + [rrna_28s_c, rrna_5_8s_c, rrna5s_c]
    complex_info = {'METABOLITES': metabolites, 'IDS': [m.id.split('_')[0] for m in metabolites], 
                                   'METABOLITE_TYPES': [m.id.split('_')[-1].split('[')[0] for m in metabolites]}
    sixty_s_complex_c, sixty_s_id = make_complex_metabolite(complex_id = '60s', **complex_info)

    rrna_5_8s_formation = cobra.Reaction('60s_maturation')
    rrna_5_8s_formation.subsytem = 'Ribosome_Biogenesis'
    rxn = dict()
    
    rxn[pre60s_complex_c], rxn[sixty_s_complex_c] = -1,1
    # exonucleolytic cleavage
    for k,v in nmp_map_c.items():
        rxn[v] = base_counts_deg[k]
    rxn[h2o_c] = -len(deg_seq)
    rxn[h_c] = len(deg_seq)

    rrna_5_8s_formation.add_metabolites(rxn)
    rrna_5_8s_formation.gene_reaction_rule = ERI1[0]
    
    
    #5.8s and 28s degradation------------------------------------------------------------------------------------
    rrna_28s_degradation = rna_exonucleolytic_degradation(rrna_28s_c, base_counts_rrna_28s, rrna_28s_seq, 
                                                     reaction_name = '28s_rrna', 
                                                     triphosphate = False)
    rrna_28s_degradation = update_rrna_degradation(rrna_28s_degradation, nucleus = False)
    
    rrna_5_8s_degradation = rna_exonucleolytic_degradation(rrna_5_8s_c, base_counts_rrna_5_8s, rrna_5_8s_seq, 
                                                 reaction_name = '5_8s_rrna', 
                                                 triphosphate = False)
    rrna_5_8s_degradation = update_rrna_degradation(rrna_5_8s_degradation, nucleus = False)
    
    
    #------------------------------------------------------------------------------------
    all_reactions = [rrna_47s_transcription, rrna_45s_formation, ets_3_degradation, ets_5_frag1_degradation, 
                     rrna_30s_formation, rrna_26s_formation, ets_5_frag3_degradation, rrna_21s_formation, 
                     ets_5_frag4_degradation, rrna_21sc_formation, rrna_18se_formation, its_1_frag1_degradation, 
                     rrna_18se_processing, pre40s_complex_formation, pre40s_transport, 
                     rrna_18s_formation, rrna_18s_degradation, its_1_frag2_degradation, 
                     rrna_32s_formation, rrna_12s_28_5s_formation, rrna_28s_formation, 
                    rrna_7s_formation, rrna_5_8s_plus_40_formation, rrna_6s_formation, 
                    rrna_5_8s_formation, rrna_28s_degradation, rrna_5_8s_degradation]
    
    for i in range(len(all_reactions)): # because degradation reactions don't have proper subsystem assigned
        r = all_reactions[i]
        r.subsystem = 'Ribosome_Biogenesis'
        all_reactions[i] = r
        
    mature_ribosomal_precomplexes = [forty_s_complex_c, sixty_s_complex_c]
    mature_rrna_metabolites = [rrna_5_8s_c, rrna_28s_c, rrna_18s_c]
    

    return all_reactions, mature_ribosomal_precomplexes, mature_rrna_metabolites


# In[8]:


def build_ribosome():
    rs_expression_reactions, rs_protein_metabolites, rl_expression_reactions, rl_protein_metabolites = build_ribosome_protein_expression_reactions()
    rpl5_n = [m for m in rl_protein_metabolites if m.id == 'HGNC:10360_folded_protein[n]'][0]
    rpl11_n = [m for m in rl_protein_metabolites if m.id == 'HGNC:10301_folded_protein[n]'][0]
    rrna5s_reactions, rrna5s_complex_n, rrna5s_c = build_rrna5s_reactions(rpl5_n, rpl11_n)
    other_rrna_reactions, mature_ribosomal_precomplexes, mature_rrna_metabolites = build_other_rrna_reactions(rrna5s_complex_n, rs_protein_metabolites, rl_protein_metabolites, rpl5_n, rpl11_n, rrna5s_c)

    # ribosome complex formation
    complex_info = {'METABOLITES': mature_ribosomal_precomplexes, 
                    'IDS': ['mature', 'ribosome'], 
                    'METABOLITE_TYPES': ['complex', 'complex']}
    ribosome_complex_formation, ribosome_complex_c = form_complex(**complex_info)
    # add a gtp hydrolysis to the complex formation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5861459/
    rxn = ribosome_complex_formation.metabolites.copy()
    rxn[ntp_map_c['G']], rxn[h2o_c], rxn[ndp_map_c['G']], rxn[pi_c], rxn[h_c]  = -1, -1, 1, 1, 1
    ribosome_complex_formation.add_metabolites(rxn)
    
    ribosome_complex_formation.id = 'RIBOSOME_COMPLEX_FORMATIONc'
    ribosome_complex_formation.lower_bound = 0
    ribosome_complex_formation.gene_reaction_rule = ' and '.join(eifs)

    # ribosome complex dissociation
    ribosome_complex_dissociation = cobra.Reaction('RIBOSOME_COMPLEX_DISSOCIATIONc')
    ind_mets = rl_protein_metabolites + rs_protein_metabolites + mature_rrna_metabolites + [rrna5s_c]
    rxn = {m: -1 for m in (ind_mets) if m.compartment == 'c'}
    rxn[ribosome_complex_c] = 1
    ribosome_complex_dissociation.add_metabolites(rxn)

    all_reactions = rrna5s_reactions + other_rrna_reactions + rs_expression_reactions +  rl_expression_reactions
    all_reactions += [ribosome_complex_formation, ribosome_complex_dissociation]
    return  all_reactions, ribosome_complex_c


# In[12]:


ribosomal_reactions, ribosome_complex_c = build_ribosome()
del psim_rib


# # Original build_rrna_reactions
# orignal version, where all RRNA processing was done independent of complex formation with ribosomal proteins
# based on fig. 2: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6315592/#B64-biomolecules-08-00123

# In[10]:


# def update_rrna_degradation(rrna_degradation_reaction, nucleus = True):
#     rrna_degradation_reaction.subsytem = 'Ribosome_Biogenesis'
#     if nucleus:
#         rule_part2 = ' and '.join(lariat_machinery['Exosome'] + lariat_machinery['NEXT Complex']) + ')'
#         rrna_degradation_reaction.gene_reaction_rule = lariat_machinery["5' Degradation"][0] + ' or (' + rule_part2
#     else:
#         rrna_degradation_reaction.gene_reaction_rule = ' and '.join(exosome['HGNC ID (gene)'].tolist())
#     return rrna_degradation_reaction

# def build_rrna5s_reactions():
    
#     # TRANSCRIPTION - basically emulates Transcript.transcript_elongation reaction
#     pre_rrna5s_n, pre_base_counts5s = make_rna_metabolite('pre_5s', pre_rrna_5s_seq, molecule_type = 'rrna', compartment = 'n')

#     rrna5s_transcription = cobra.Reaction('TRANSCRIPTION_PRE_RRNA5s')
#     rrna5s_transcription.subsytem = 'Ribosome_Biogenesis'
#     rxn = dict()
#     for ntp, base_letter in seq_metabolite_map.items():
#         rxn[ntp] = -1*pre_base_counts5s[base_letter]
#     rxn[ppi_n] = len(pre_rrna_5s_seq) - 1
#     rxn[pre_rrna5s_n] = 1
#     rrna5s_transcription.add_metabolites(rxn)
#     rrna5s_transcription.gene_reaction_rule = ' and '.join(rnap3_transcription_machinery) 
    
#     # PROCESSING - mature rrna (3->5' exonucleolytic cleave of last 24 bases)
#     rrna5s_processing = cobra.Reaction('PROCESSING_RRNA5s')
#     rrna5s_processing.subsytem = 'Ribosome_Biogenesis'
#     rrna5s_n, base_counts5s = make_rna_metabolite('5s', rrna_5s_seq, molecule_type = 'rrna', compartment = 'n') 
#     deg_base_counts = dict()
#     for k,v in pre_base_counts5s.items():
#         deg_base_counts[k] = v - base_counts5s[k]

#     rxn = dict()
#     rxn[h2o_n] = -sum(deg_base_counts.values()) # no -1 because all bonds 5'-most bond cleave
#     rxn[h_n] = sum(deg_base_counts.values())
#     for k,v in nmp_map_n.items():
#         rxn[v] = deg_base_counts[k]


#     rxn[pre_rrna5s_n] = -1
#     rxn[rrna5s_n] = 1


#     rrna5s_processing.add_metabolites(rxn)
#     rrna5s_processing.gene_reaction_rule = REXO5

    
#     # TRANSPORT 
#     # must add nucleocytoplasmic export via ran gtp: https://www.sciencedirect.com/science/article/pii/S0171933504702575?via%3Dihub
#     rrna5s_c = rrna5s_n.copy()
#     rrna5s_c.id = rrna5s_c.id.replace('[n]', '[c]')
#     rrna5s_c.compartment = 'c'
    
#     rrna5s_transport = cobra.Reaction('RRNA5s_NUCLEAR_EXPORTtn')
#     rrna5s_transport.name = 'rRNA5s nuclear export'
#     rrna5s_transport.subsytem = 'Ribosome_Biogenesis'
#     rxn = {rrna5s_n: -1, rrna5s_c: 1}
#     # gtp hydrolysis on cytoplasmic side for export (see protein_expression nuclear_transport for details)
#     rxn[ntp_map_c['G']], rxn[h2o_c], rxn[ndp_map_c['G']], rxn[pi_c], rxn[h_c]  = -1, -1, 1, 1, 1
#     rrna5s_transport.add_metabolites(rxn)
#     rrna5s_transport.gene_reaction_rule = ' and '.join(tfiiia + RAN + XPO1)
    
#     # Degradation
#     rrna5s_degradation = rna_exonucleolytic_degradation(rrna5s_c, base_counts5s, rrna_5s_seq, reaction_name = '5s_rRNA')
#     rrna5s_degradation.subsytem = 'Ribosome_Biogenesis'
#     rrna5s_degradation.gene_reaction_rule = ' and '.join(exosome['HGNC ID (gene)'].tolist())
    
#     rrna5s_reactions = [rrna5s_transcription, rrna5s_processing, rrna5s_transport, rrna5s_degradation]
    
    
#     return rrna5s_reactions, rrna5s_c



# # ets_5_frag1 is from 5' end of 47s to A' site
# # ets_5_frag2 is from A' to 18s
# # ets_5_frag3 is from A' to A0
# # ets_5_frag4 is from A0 site to site 1 (start of 18s)
# # its_1_frag1_seq is between site E and teh conserved stall location of RRP6
# # its_1_frag2_seq is less than E site (some degradation) + a polyA/U tail

# def build_other_rrna_reactions():
#     # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6315592/ figure 2
    
#     # 47s transcription------------------------------------------------------------------------------------
#     rrna_47s_n, base_counts_47s = make_rna_metabolite('47s', rrna_47s_seq, molecule_type = 'rrna', compartment = 'n')
#     rrna_47s_transcription = cobra.Reaction('TRANSCRIPTION_RRNA_47s')
#     rrna_47s_transcription.subsytem = 'Ribosome_Biogenesis'
#     rxn = dict()
#     for ntp, base_letter in seq_metabolite_map.items():
#         rxn[ntp] = -1*base_counts_47s[base_letter]
#     rxn[ppi_n] = len(rrna_47s_seq) - 1
#     rxn[rrna_47s_n] = 1
#     rrna_47s_transcription.add_metabolites(rxn)
#     rrna_47s_transcription.gene_reaction_rule = ' and '.join(rnap1['HGNC ID (gene)'].tolist() + rnap1_tfs)  
    
#     # 45s formation------------------------------------------------------------------------------------
#     ets_5_frag1_seq = ets_5_seq[:A_prime_index] 
#     ets_5_frag2_seq = ets_5_seq[A_prime_index:]
#     rrna_45s_seq = rrna_47s_seq[A_prime_index:rrna_47s_seq.index(rrna_28s_seq) + len(rrna_28s_seq)]

#     ets_5_frag1_n, base_counts_ets_5_frag1 = make_rna_metabolite('ets_5_frag1', ets_5_frag1_seq, 
#                                                                  molecule_type = 'rrna', compartment = 'n')
#     rrna_45s_n, base_counts_rrna_45s = make_rna_metabolite('45s', rrna_45s_seq, molecule_type = 'rrna', 
#                                                            compartment = 'n', triphosphate=False)
#     ets_3_n, base_counts_ets_3 = make_rna_metabolite('ets_3', ets_3_seq, molecule_type = 'rrna', 
#                                                      compartment = 'n', triphosphate = False)

#     rrna_45s_formation = cobra.Reaction('FORMATION_RRNA_45s')
#     rrna_45s_formation.subsytem = 'Ribosome_Biogenesis'
#     rxn = dict()
#     rxn[h2o_n] = -2 # 2 endonuclolytic cleavage events to go from 47s to 45s
#     rxn[rrna_47s_n], rxn[rrna_45s_n], rxn[ets_3_n], rxn[ets_5_frag1_n] = -1, 1, 1, 1
#     rrna_45s_formation.add_metabolites(rxn)
#     rrna_45s_formation.gene_reaction_rule = ' and '.join(UTP10 + RNASEN)
    
#     ets_3_degradation = rna_exonucleolytic_degradation(ets_3_n, base_counts_ets_3, ets_3_seq, reaction_name = 'ets_3_rRNA', triphosphate = False)
#     ets_3_degradation = update_rrna_degradation(ets_3_degradation)
#     ets_5_frag1_degradation = rna_exonucleolytic_degradation(ets_5_frag1_n, base_counts_ets_5_frag1, ets_5_frag1_seq, reaction_name = 'ets_5_frag1_rRNA', triphosphate = True)
#     ets_5_frag1_degradation = update_rrna_degradation(ets_5_frag1_degradation)
    
#     #45S-->30S + 32.5S------------------------------------------------------------------------------------
#     idx_30s = rrna_45s_seq.index(rrna_18s_seq) + len(rrna_18s_seq)
#     rrna_30s_seq = ets_5_frag2_seq + rrna_18s_seq + rrna_45s_seq[idx_30s:idx_30s + site_2_index]
#     rrna_32_5s_seq = rrna_45s_seq[idx_30s + site_2_index:]

#     rrna_30s_n, base_counts_rrna_30s = make_rna_metabolite('30s', rrna_30s_seq, molecule_type = 'rrna', 
#                                                            compartment = 'n', triphosphate=False)
#     rrna_32_5s_n, base_counts_rrna_32_5s = make_rna_metabolite('32_5s', rrna_32_5s_seq, molecule_type = 'rrna', 
#                                                                compartment = 'n', triphosphate=False)

#     rrna_30s_formation = cobra.Reaction('FORMATION_RRNA_30s_32_5s')
#     rrna_30s_formation.subsytem = 'Ribosome_Biogenesis'
#     rxn = dict()
#     rxn[h2o_n] = -1 # endonuclolytic cleavage event at site 2
#     rxn[rrna_45s_n], rxn[rrna_30s_n], rxn[rrna_32_5s_n] = -1, 1, 1
#     rrna_30s_formation.add_metabolites(rxn)
#     rrna_30s_formation.gene_reaction_rule = RMRP[0]
    
#     #26s formation------------------------------------------------------------------------------------
#     rrna_26s_seq = ets_5_frag2_seq[A_0_index:] + rrna_18s_seq + its_1_seq[:site_2_index]
#     ets_5_frag3_seq = ets_5_frag2_seq[:A_0_index]
#     rrna_26s_n, base_counts_rrna_26s = make_rna_metabolite('26s', rrna_26s_seq, molecule_type = 'rrna', 
#                                                            compartment = 'n', triphosphate=False)
#     ets_5_frag3_n, base_counts_ets_5_frag3 = make_rna_metabolite('ets_5_frag3', ets_5_frag3_seq, 
#                                                                  molecule_type = 'rrna', compartment = 'n', 
#                                                                  triphosphate=False)

#     rrna_26s_formation = cobra.Reaction('FORMATION_RRNA_26s')
#     rrna_26s_formation.subsytem = 'Ribosome_Biogenesis'
#     rxn = dict()
#     rxn[h2o_n] = -1 # endonuclolytic cleavage event at site 2
#     rxn[rrna_30s_n], rxn[rrna_26s_n], rxn[ets_5_frag3_n] = -1, 1, 1
#     rrna_26s_formation.add_metabolites(rxn)
#     rrna_26s_formation.gene_reaction_rule = UTP23[0]

#     ets_5_frag3_degradation = rna_exonucleolytic_degradation(ets_5_frag3_n, base_counts_ets_5_frag3, ets_5_frag3_seq, reaction_name = 'ets_5_frag3_rRNA', triphosphate = False)
#     ets_5_frag3_degradation = update_rrna_degradation(ets_5_frag3_degradation)
    
#     # 21S formation------------------------------------------------------------------------------------
#     rrna_21s_seq = rrna_18s_seq + its_1_seq[:site_2_index]
#     ets_5_frag4_seq = ets_5_frag2_seq[A_0_index:]
#     rrna_21s_n, base_counts_rrna_21s = make_rna_metabolite('21s', rrna_21s_seq, molecule_type = 'rrna', compartment = 'n',
#                                                            triphosphate=False)
#     ets_5_frag4_n, base_counts_ets_5_frag4 = make_rna_metabolite('ets_5_frag4', ets_5_frag4_seq, molecule_type = 'rrna', compartment = 'n',
#                                                                  triphosphate=False)

#     rrna_21s_formation = cobra.Reaction('FORMATION_RRNA_21s')
#     rrna_21s_formation.subsytem = 'Ribosome_Biogenesis'
#     rxn = dict()
#     rxn[h2o_n] = -1 # endonuclolytic cleavage event at site 1
#     rxn[rrna_26s_n], rxn[rrna_21s_n], rxn[ets_5_frag4_n] = -1, 1, 1
#     rrna_21s_formation.add_metabolites(rxn)
#     rrna_21s_formation.gene_reaction_rule = UTP24[0]
#     ets_5_frag4_degradation = rna_exonucleolytic_degradation(ets_5_frag4_n, base_counts_ets_5_frag4, ets_5_frag4_seq, reaction_name = 'ets_5_frag4_rRNA', triphosphate = False)
#     ets_5_frag4_degradation = update_rrna_degradation(ets_5_frag4_degradation)
#     # 21SC formation------------------------------------------------------------------------------------

#     rrna_21sc_seq = rrna_18s_seq + its_1_seq[:conserved_stall_idx]
#     rrna_21sc_n, base_counts_rrna_21sc = make_rna_metabolite('21sc', rrna_21sc_seq, molecule_type = 'rrna', 
#                                                              compartment = 'n', triphosphate=False)

#     deg_seq = its_1_seq[conserved_stall_idx: site_2_index]
#     base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

#     rrna_21sc_formation = cobra.Reaction('FORMATION_RRNA_21sc')
#     rrna_21sc_formation.subsytem = 'Ribosome_Biogenesis'
#     rxn = dict()
#     rxn[rrna_21s_n], rxn[rrna_21sc_n] = -1,1
#     # exonucleolytic cleavage
#     for k,v in nmp_map_n.items():
#         rxn[v] = base_counts_deg[k]
#     rxn[h2o_n] = -len(deg_seq)
#     rxn[h_n] = len(deg_seq)
#     rrna_21sc_formation.add_metabolites(rxn)
#     rrna_21sc_formation.gene_reaction_rule = exosome[exosome['Approved symbol'] == 'EXOSC10']['HGNC ID (gene)'].tolist()[0]

#     # 18se formation------------------------------------------------------------------------------------
#     rrna_18se_seq = rrna_18s_seq + its_1_seq[:e_index]
#     its_1_frag1_seq = its_1_seq[e_index:conserved_stall_idx]
#     rrna_18se_n, base_counts_rrna_18se = make_rna_metabolite('18se', rrna_18se_seq, molecule_type = 'rrna', 
#                                                              compartment = 'n', triphosphate=False)
#     its_1_frag1_n, base_counts_its_1_frag1 = make_rna_metabolite('its_1_frag1', its_1_frag1_seq, molecule_type = 'rrna', compartment = 'n', triphosphate=False)

#     rrna_18se_formation = cobra.Reaction('FORMATION_RRNA_18se')
#     rrna_18se_formation.subsytem = 'Ribosome_Biogenesis'
#     # endonuclolytic cleavage event at site E
#     rrna_18se_formation.add_metabolites({h2o_n: -1, rrna_21sc_n: -1, rrna_18se_n: 1, its_1_frag1_n: 1})
#     rrna_18se_formation.gene_reaction_rule = UTP24[0]
#     its_1_frag1_degradation = rna_exonucleolytic_degradation(its_1_frag1_n, base_counts_its_1_frag1, its_1_frag1_seq, reaction_name = 'its_1_frag1_rRNA', triphosphate = False)
#     its_1_frag1_degradation = update_rrna_degradation(its_1_frag1_degradation)
#     # 18se nuclear processing------------------------------------------------------------------------------------
#     rrna_18se_processed_seq = rrna_18se_seq[:-int(0.75*e_index)] # degradation of 60/80 bps of ITS1 by PARN
#     rrna_18se_processed_seq += 'U'*int(0.125*e_index)+'A'*int(0.125*e_index) # polyU by PAPD5
#     deg_seq = rrna_18se_seq[-int(0.75*e_index):]

#     rrna_18se_processed_n, base_counts_rrna_18se_processed = make_rna_metabolite('18se_processed', rrna_18se_processed_seq, molecule_type = 'rrna', compartment = 'n', triphosphate=False)
#     base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

#     rrna_18se_processing = cobra.Reaction('PROCESSING_RRNA_18se')
#     rrna_18se_processing.subsytem = 'Ribosome_Biogenesis'
#     rxn = dict()
#     rxn[rrna_18se_n], rxn[rrna_18se_processed_n] = -1,1
#     # exonucleolytic cleavage
#     for k,v in nmp_map_n.items():
#         rxn[v] = base_counts_deg[k]
#     rxn[h2o_n] = -len(deg_seq)
#     rxn[h_n] = len(deg_seq)
#     # polyU/A tail synthesis<--unsure why ppi_n has -1 but it mass balances
#     rxn[atp_n], rxn[ntp_map_n['U']], rxn[ppi_n] = -int(0.125*e_index), -int(0.125*e_index), int(0.25*e_index)-1

#     rrna_18se_processing.add_metabolites(rxn)
#     rrna_18se_processing.gene_reaction_rule = ' and '.join(PARN + PAPD5)

#     # 18se nucleocytoplasmic export-----------------------------------------------------------------------
#     rrna_18se_processed_c = rrna_18se_processed_n.copy()
#     rrna_18se_processed_c.id = rrna_18se_processed_c.id.replace('[n]', '[c]')
#     rrna_18se_processed_c.compartment = 'c'

#     rrna_18se_transport = cobra.Reaction('RRNA_18se_tn')
#     rrna_18se_transport.subsytem = 'Ribosome_Biogenesis'
#     rrna_18se_transport.name = 'rRNA18se nuclear export'
#     rxn = {rrna_18se_processed_n: -1, rrna_18se_processed_c: 1}
#     rxn[ntp_map_c['G']], rxn[h2o_c], rxn[ndp_map_c['G']], rxn[pi_c], rxn[h_c]  = -1, -1, 1, 1, 1
#     rrna_18se_transport.add_metabolites(rxn)
#     rrna_18se_transport.gene_reaction_rule = ' and '.join(tfiiia + RAN + XPO1)
     
    
#     # 18s formation------------------------------------------------------------------------------------
#     its_1_frag2_seq = its_1_seq[:int(0.25*e_index)+1] + 'U'*int(0.125*e_index)+'A'*int(0.125*e_index)
#     rrna_18s_c, base_counts_rrna_18s = make_rna_metabolite('18s', rrna_18s_seq, molecule_type = 'rrna', 
#                                                            compartment = 'c', triphosphate=False)
#     its_1_frag2_c, base_counts_its_1_frag2 = make_rna_metabolite('its_1_frag2', its_1_frag2_seq, 
#                                                                  molecule_type = 'rrna', compartment='c',
#                                                                  triphosphate=False)

#     rrna_18s_formation = cobra.Reaction('FORMATION_RRNA_18s')
#     rrna_18s_formation.subsytem = 'Ribosome_Biogenesis'
#     # endonuclolytic cleavage event at site 3
#     rrna_18s_formation.add_metabolites({h2o_n: -1, rrna_18se_processed_c: -1, rrna_18s_c: 1, its_1_frag2_c: 1})
#     rrna_18s_formation.gene_reaction_rule = NOB1[0]

#     its_1_frag2_degradation = rna_exonucleolytic_degradation(its_1_frag2_c, base_counts_its_1_frag2, its_1_frag2_seq, reaction_name = 'its_1_frag2_rRNA', triphosphate = False, nucleus = False)
#     its_1_frag2_degradation = update_rrna_degradation(its_1_frag2_degradation, nucleus = False)
    
#     # 32S formation------------------------------------------------------------------------------------
#     deg_seq = its_1_seq[site_2_index:]
#     rrna_32s_seq = rrna_32_5s_seq[len(deg_seq):]

#     rrna_32s_n, base_counts_rrna_32s = make_rna_metabolite('32s', rrna_32s_seq, molecule_type = 'rrna', 
#                                                            compartment = 'n',triphosphate=False)
#     base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

#     rrna_32s_formation = cobra.Reaction('FORMATION_RRNA_32s')
#     rrna_32s_formation.subsytem = 'Ribosome_Biogenesis'

#     rxn = dict()
#     rxn[rrna_32_5s_n], rxn[rrna_32s_n] = -1,1
#     # exonucleolytic cleavage
#     for k,v in nmp_map_n.items():
#         rxn[v] = base_counts_deg[k]
#     rxn[h2o_n] = -len(deg_seq)
#     rxn[h_n] = len(deg_seq)

#     rrna_32s_formation.add_metabolites(rxn)
#     rrna_32s_formation.gene_reaction_rule = lariat_machinery["5' Degradation"][0]
    
#     #32s-->12s + 28.5s------------------------------------------------------------------------------------
 
#     rrna_12s_seq = rrna_5_8s_seq + its_2_seq[:site_4_index]
#     rrna_28_5s_seq = its_2_seq[site_4_index:] + rrna_28s_seq
#     rrna_12s_n, base_counts_rrna_12s = make_rna_metabolite('12s', rrna_12s_seq, molecule_type = 'rrna', compartment = 'n', 
#                                                            triphosphate=False)
#     rrna_28_5s_n, base_counts_rrna_28_5s = make_rna_metabolite('28_5s', rrna_28_5s_seq, molecule_type = 'rrna', compartment = 'n', 
#                                                                triphosphate=False)

#     rrna_12s_28_5s_formation = cobra.Reaction('FORMATION_RRNA_12s_28_5s')
#     rrna_12s_28_5s_formation.subsytem = 'Ribosome_Biogenesis'
#     rxn = dict()
#     rxn[h2o_n] = -1 # endonuclolytic cleavage event at site 4
#     rxn[rrna_32s_n], rxn[rrna_12s_n], rxn[rrna_28_5s_n] = -1,1,1
#     rrna_12s_28_5s_formation.add_metabolites(rxn)
#     rrna_12s_28_5s_formation.gene_reaction_rule = LAS1[0]
    
#     #28s formation------------------------------------------------------------------------------------
#     deg_seq = rrna_28_5s_seq[:rrna_28_5s_seq.index(rrna_28s_seq)]

#     rrna_28s_n, base_counts_rrna_28s = make_rna_metabolite('28s', rrna_28s_seq, molecule_type = 'rrna', compartment = 'n', 
#                                                            triphosphate=False)
#     base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

#     rrna_28s_formation = cobra.Reaction('FORMATION_RRNA_28s')
#     rrna_28s_formation.subsytem = 'Ribosome_Biogenesis'

#     rxn = dict()
#     rxn[rrna_28_5s_n], rxn[rrna_28s_n] = -1,1
#     # exonucleolytic cleavage
#     for k,v in nmp_map_n.items():
#         rxn[v] = base_counts_deg[k]
#     rxn[h2o_n] = -len(deg_seq)
#     rxn[h_n] = len(deg_seq)

#     rrna_28s_formation.add_metabolites(rxn)
#     rrna_28s_formation.gene_reaction_rule = lariat_machinery["5' Degradation"][0]
    
#     #28s export------------------------------------------------------------------------------------
    
#     rrna_28s_c = rrna_28s_n.copy()
#     rrna_28s_c.id = rrna_28s_c.id.replace('[n]', '[c]')
#     rrna_28s_c.compartment = 'c'

#     rrna_28s_transport = cobra.Reaction('RRNA_28s_tn')
#     rrna_28s_transport.subsytem = 'Ribosome_Biogenesis'
#     rrna_28s_transport.name = 'rRNA28s nuclear export'
#     rrna_28s_transport.add_metabolites({rrna_28s_n: -1, rrna_28s_c: 1})
#     rrna_28s_transport.gene_reaction_rule = 'xpo1_nucleocytoplasmic_export'
    
#     #7s formation------------------------------------------------------------------------------------
#     deg_seq = its_2_seq[seven_s_idx: site_4_index]
#     rrna_7s_seq = rrna_5_8s_seq + its_2_seq[:seven_s_idx]
#     rrna_7s_n, base_counts_rrna_7s = make_rna_metabolite('7s', rrna_7s_seq, molecule_type = 'rrna', compartment = 'n', triphosphate=False)
#     base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

#     rrna_7s_formation = cobra.Reaction('FORMATION_RRNA_7s')
#     rrna_7s_formation.subsytem = 'Ribosome_Biogenesis'

#     rxn = dict()
#     rxn[rrna_12s_n], rxn[rrna_7s_n] = -1,1
#     # exonucleolytic cleavage
#     for k,v in nmp_map_n.items():
#         rxn[v] = base_counts_deg[k]
#     rxn[h2o_n] = -len(deg_seq)
#     rxn[h_n] = len(deg_seq)

#     rrna_7s_formation.add_metabolites(rxn)
#     rrna_7s_formation.gene_reaction_rule = ' and '.join(DIS3 + ISG20L2)
    
#     #5.8s+40 formation------------------------------------------------------------------------------------
#     deg_seq = its_2_seq[five_eight_plus_forty_idx: seven_s_idx]
#     rrna_5_8s_plus_40_seq = rrna_5_8s_seq + its_2_seq[:five_eight_plus_forty_idx]
#     rrna_5_8s_plus_40_n, base_counts_rrna_5_8s_plus_40 = make_rna_metabolite('5_8s_plus_40', rrna_5_8s_plus_40_seq, molecule_type = 'rrna', compartment = 'n', triphosphate=False)
#     base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

#     rrna_5_8s_plus_40_formation = cobra.Reaction('FORMATION_RRNA_5_8s_plus_40')
#     rrna_5_8s_plus_40_formation.subsytem = 'Ribosome_Biogenesis'

#     rxn = dict()
#     rxn[rrna_7s_n], rxn[rrna_5_8s_plus_40_n] = -1,1
#     # exonucleolytic cleavage
#     for k,v in nmp_map_n.items():
#         rxn[v] = base_counts_deg[k]
#     rxn[h2o_n] = -len(deg_seq)
#     rxn[h_n] = len(deg_seq)

#     rrna_5_8s_plus_40_formation.add_metabolites(rxn)
#     rrna_5_8s_plus_40_formation.gene_reaction_rule = ' and '.join(DIS3 + ISG20L2)
    
#     #6s export------------------------------------------------------------------------------------

#     deg_seq = its_2_seq[six_s_index: five_eight_plus_forty_idx]
#     rrna_6s_seq = rrna_5_8s_seq + its_2_seq[:six_s_index]
#     rrna_6s_n, base_counts_rrna_6s = make_rna_metabolite('6s', rrna_6s_seq, molecule_type = 'rrna', compartment = 'n', 
#                                                          triphosphate=False)
#     base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

#     rrna_6s_formation = cobra.Reaction('FORMATION_RRNA_6s')
#     rrna_6s_formation.subsytem = 'Ribosome_Biogenesis'

#     rxn = dict()
#     rxn[rrna_5_8s_plus_40_n], rxn[rrna_6s_n] = -1,1
#     # exonucleolytic cleavage
#     for k,v in nmp_map_n.items():
#         rxn[v] = base_counts_deg[k]
#     rxn[h2o_n] = -len(deg_seq)
#     rxn[h_n] = len(deg_seq)

#     rrna_6s_formation.add_metabolites(rxn)
#     rrna_6s_formation.gene_reaction_rule = exosome[exosome['Approved symbol'] == 'EXOSC10']['HGNC ID (gene)'].tolist()[0]

#     #6s transport------------------------------------------------------------------------------------
#     rrna_6s_c = rrna_6s_n.copy()
#     rrna_6s_c.id = rrna_6s_c.id.replace('[n]', '[c]')
#     rrna_6s_c.compartment = 'c'

#     rrna_6s_transport = cobra.Reaction('RRNA_6s_tn')
#     rrna_6s_transport.subsytem = 'Ribosome_Biogenesis'
#     rrna_6s_transport.name = 'rRNA6s nuclear export'
#     rrna_6s_transport.add_metabolites({rrna_6s_n: -1, rrna_6s_c: 1})
#     rrna_6s_transport.gene_reaction_rule = 'xpo1_nucleocytoplasmic_export'
    
#     #5.8s formation------------------------------------------------------------------------------------

#     deg_seq = its_2_seq[:six_s_index]
#     rrna_5_8s_c, base_counts_rrna_5_8s = make_rna_metabolite('5_8s', rrna_5_8s_seq, molecule_type = 'rrna', compartment = 'c', 
#                                                              triphosphate=False)
#     base_counts_deg, elements_deg = get_base_counts_and_elements(deg_seq)

#     rrna_5_8s_formation = cobra.Reaction('FORMATION_RRNA_5_8s')
#     rrna_5_8s_formation.subsytem = 'Ribosome_Biogenesis'

#     rxn = dict()
#     rxn[rrna_6s_c], rxn[rrna_5_8s_c] = -1,1
#     # exonucleolytic cleavage
#     for k,v in nmp_map_c.items():
#         rxn[v] = base_counts_deg[k]
#     rxn[h2o_c] = -len(deg_seq)
#     rxn[h_c] = len(deg_seq)

#     rrna_5_8s_formation.add_metabolites(rxn)
#     rrna_5_8s_formation.gene_reaction_rule = ERI1[0]
    
    
#     #------------------------------------------------------------------------------------
#     all_reactions = [rrna_47s_transcription, rrna_45s_formation, ets_3_degradation, ets_5_frag1_degradation, 
#                      rrna_30s_formation, rrna_26s_formation, ets_5_frag3_degradation, rrna_21s_formation, 
#                      ets_5_frag4_degradation, rrna_21sc_formation, rrna_18se_formation, its_1_frag1_degradation, 
#                      rrna_18se_processing, rrna_18se_transport, 
#                      rrna_18s_formation, its_1_frag2_degradation, 
#                      rrna_32s_formation, rrna_12s_28_5s_formation, rrna_28s_formation, rrna_28s_transport, 
#                     rrna_7s_formation, rrna_5_8s_plus_40_formation, rrna_6s_formation, rrna_6s_transport, 
#                     rrna_5_8s_formation]
#     mature_rrna_metabolites = [rrna_5_8s_c, rrna_18s_c, rrna_28s_c]
    

#     return all_reactions, mature_rrna_metabolites


# In[11]:


# def build_rrna_reactions():
#     '''Reactions associated with ribosomal RNA biogenesis.'''
#     rrna5s_reactions, rrna5s_c = build_rrna5s_reactions()
#     other_rrna_reactions, other_mature_rrna_metabolites = build_other_rrna_reactions()
    
#     rrna_reactions = rrna5s_reactions + other_rrna_reactions
#     rrna_metabolites = [rrna5s_c] + other_mature_rrna_metabolites
    
#     return rrna_reactions, rrna_metabolites

