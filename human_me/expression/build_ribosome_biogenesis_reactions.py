#!/usr/bin/env python
# coding: utf-8
import random
from typing import Any, Dict, List

import pandas as pd
import numpy as np

import human_me.expression.gene_expression.build_mrna_expression_reactions as build_mrna
from human_me.core.macromolecules.complex import RibosomalComplex
from human_me.core.macromolecules.protein import Protein
from human_me.core.macromolecules.RNA import RNA_fragment, rRNA
from human_me.core.reaction import ExpressionReaction
from human_me.expression.gene_expression import gene_information
from human_me.expression.gene_expression.protein_expression import \
    build_protein_expression_reactions as build_protein
from human_me.expression.gene_expression.protein_expression import degradation
from human_me.io import HiddenPrints
from human_me.utils import functions as func
from human_me.utils import machinery as mach
from human_me.data.file_paths import build_files_url

# # rRNA
func.read_fasta_url(build_files_url + '45s_rrna_seq.txt')
# rrna sequences
# assume the ncbi 45s is actually 47s...see notes for details
rrna_47s_seq = func.read_fasta_url(build_files_url + '45s_rrna_seq.txt').seq.transcribe()
rrna_18s_seq = func.read_fasta_url(build_files_url + '18s_rrna_seq.txt').seq.transcribe()
rrna_28s_seq = func.read_fasta_url(build_files_url + '28s_rrna_seq.txt').seq.transcribe()
rrna_5_8s_seq = func.read_fasta_url(build_files_url + '5_8s_rrna_seq.txt').seq.transcribe()
ets_5_seq = rrna_47s_seq[:rrna_47s_seq.index(rrna_18s_seq)]
its_1_seq = rrna_47s_seq[rrna_47s_seq.index(rrna_18s_seq) + len(rrna_18s_seq):rrna_47s_seq.index(rrna_5_8s_seq)]
its_2_seq = rrna_47s_seq[rrna_47s_seq.index(rrna_5_8s_seq) + len(rrna_5_8s_seq):rrna_47s_seq.index(rrna_28s_seq)]
ets_3_seq = rrna_47s_seq[rrna_47s_seq.index(rrna_28s_seq) + len(rrna_28s_seq):]

pre_rrna_5s_seq = func.read_fasta_url(build_files_url + '5s_rrna_seq.txt').seq.transcribe()
rrna_5s_seq = pre_rrna_5s_seq[:120]  # 120 is length of mature 5s_rrna

# rrna cut sites
# cut site indexes, relative to how far right (3' end is right) they are of a certain feature
# Fig. 3b https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4361047/
# scaled from length in figure to length of my sequence
A_prime_index = int(round(np.median([414, 420]) * len(ets_5_seq) / 3657))  # location from 5' end of 47s
A_0_index = int(round(1642 * len(ets_5_seq) / 3657)) - A_prime_index  # how far to the right of 45s is the A_0 site
site_2_index = int(round((6470 - 5527) * len(its_1_seq) / (6623 - 5527)))  # how far right of end of 18s
site_4_index = int(round((7570 - 6779) * (len(its_2_seq) / (
    7935 - 6779))))  # how far to the right of the end of 5.8s/how far into ITS2 is the site 4 cut location
e_index = int(round(
    (np.median([5606, 5609]) - 5527) * len(its_1_seq) / (6623 - 5527)))  # bp to right of end of 18s/start of its_1
conserved_stall_idx = int(round((np.median([6117, 6192]) - 5527) * len(its_1_seq) / (
    6623 - 5527)))  # how far right of end of 18s does RRP6 stall to form 21S-C
# Fig. 6b https://www.sciencedirect.com/science/article/pii/S1097276513005844?via%3Dihub
SEVEN_S_IDX = 190
FIVE_EIGHT_PLUS_FORTY_IDX = 40
SIX_S_IDX = 1  # https://www.nature.com/articles/s41594-019-0234-x?draft=collection

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
# SEVEN_S_IDX = int(round((yeast_c2/yeast_its2_length)*len(its_2_seq))) # how far to the right of 5.8s does sequence extend to form the 7s rrrna


# In[3]:
exclude = ['POLYUBIQUITINATIONn', 'DEUBIQUITINATIONn', 'DEGRADATIONn']

def set_ribosomal_psim(psim_me: pd.DataFrame) -> pd.DataFrame:
    """Format psim information specifically for ribosomal biogenesis

    Parameters
    ----------
    psim_me : pd.DataFrame
        see PSIM_README.md for details

    Returns
    -------
    psim_rib : pd.DataFrame
        formatted dataframe
    """
    psim_rib = psim_me.copy()
    psim_rib.LOCATION = psim_rib.LOCATION.apply(lambda x: ['n', 'c'])
    return psim_rib
    


def cleave_ub(hgnc_id: str, model_metabolites, psim_rib, ub_args: Dict[str, Any], compress_mrna: bool, 
              charged_trna_map, modified_trna_transcript_c,
              stochastic: bool, seeds: List[int]):
    """Generates reactions specific for ubiquitin-protein fusions. RPL40 and RPS27A have ubiquitin fusions.

    Parameters
    ----------
    hgnc_id : str
        Gene HGNC ID
    me_input_model : cobra.Model
        'the corrected input metabolic model (as provided in preprocess.correct_inputs.correct_model)'
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin
    psim_rib : pd.DataFrame
        formatted PSIM for ribosomal biogenesis
    ub_args : Dict[str, Any]
        set of variables generated by this ubiquitin.express_ubiquitin to be used throughout model building
    compress_mrna : bool
        whether to condense elongation, processing, and nuclear export reactions into a single reaction
    charged_trna_map : Dict[str, macromolecules.RNA.tRNA]
        output of creat_trna()
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin    
    stochastic : bool
        Whether a potentially stochastic output should be stochastic, or choose a default behavior instead
    seeds : List[int]
        A list of seeds for if stochastic is set to True
    """
    gene_info = gene_information.generate_from_psim(hgnc_id=hgnc_id, psim=psim_rib,
                                                    machinery_list=list(), reactions=None, nonmachinery_locations=['n'],
                                                    stochastic=stochastic, seed=seeds[0])
    gene_info = func.convert_gi(gene_info, non_machinery=dict())

    mrna_expression_reactions, mrna_transcript_c, mrna_deg_proxy = build_mrna.get_mrna_expression_reactions(gene_info=gene_info,
                                                                                                            model_metabolites=model_metabolites,
                                                                                                            compress_mrna=compress_mrna)
    translation_elongation_c, unfolded_protein_c = build_protein.c_trln.translate_protein_cytosolic(gene_info, 
                                                                                                    mrna_transcript_c,
                                                                                                    mrna_deg_proxy, 
                                                                                                    charged_trna_map=charged_trna_map, 
                                                                                                    modified_trna_transcript_c=modified_trna_transcript_c,
                                                                                                    model_metabolites=model_metabolites,)
    translation_elongation_c.ubiquitin_biogenesis = True
    # cleaved protein sequence, gene_info object, and cobra.Metabolite
    processed_seq = gene_info.protein_seq[
        :gene_info.protein_seq.index(ub_args['SINGLE_UB_SEQ'])] + gene_info.protein_seq[
        gene_info.protein_seq.index(
            ub_args[
                'SINGLE_UB_SEQ']) + len(
            ub_args[
                'SINGLE_UB_SEQ']):]
    psim_temp = psim_rib.copy()
    psim_temp.loc[psim_temp[psim_temp.HGNC_ID == hgnc_id].index, 'PROTEIN_SEQ'] = processed_seq
    gene_info = gene_information.generate_from_psim(hgnc_id=hgnc_id, psim=psim_temp,
                                                    machinery_list=list(), reactions=None, nonmachinery_locations=['n'],
                                                    stochastic=stochastic, seed=seeds[1])
    del psim_temp
    gene_info = func.convert_gi(gene_info, non_machinery=dict())

    processed_unfolded_protein_c = Protein(id_='processed_unfolded', compartment='c', model_metabolites=model_metabolites,
                                           gene_info=gene_info)
    ub_cleavage = ExpressionReaction(id=gene_info.hgnc_id + '_UBIQUITIN_CLEAVAGEc', hgnc_id=gene_info.hgnc_id,
                                     subsystem='Protein_Expression', ribosome_biogenesis=True,
                                     ubiquitin_biogenesis=True)
    ub_cleavage.add_metabolites({unfolded_protein_c: -1, model_metabolites.h2o_c: -1,
                                 ub_args['ub_c']: 1, processed_unfolded_protein_c: 1})
    ub_cleavage.gene_reaction_rule = mach.UCHL3[0]

    protein_folding_cytosolic, folded_protein_c = build_protein.fold_protein_cytosolic(gene_info,
                                                                                       processed_unfolded_protein_c, 
                                                                                       model_metabolites=model_metabolites)
    folded_protein_c.alpha_p = gene_info.coupling_params['alpha_p']

    # for degradation------------------
    folded_protein_c._L_protein = gene_info.L_protein
    folded_protein_c._amino_acid_counts = gene_info.amino_acid_counts
    folded_protein_c._ptms = gene_info.ptms
    folded_protein_c._hgnc_id = gene_info.hgnc_id
    # ------------------
    nuclear_import, folded_protein_n = build_protein.transport_nuclear_protein(gene_info, folded_protein_c, model_metabolites=model_metabolites)
    dcp = degradation.degrade(folded_protein_c, model_metabolites=model_metabolites, **{'ub_args': ub_args})

    to_add = [translation_elongation_c, ub_cleavage, protein_folding_cytosolic,
              nuclear_import] + dcp + mrna_expression_reactions

    return to_add, folded_protein_c, folded_protein_n


def build_ribosome_protein_expression_reactions(model_metabolites, psim_rib, modified_trna_transcript_c, charged_trna_map,
                                                ub_args: Dict[str, Any], compress_mrna: bool, stochastic: bool, seeds: List[int]):
    """Reactions associated with transcription and translation of ribosomal proteins

    Parameters
    ----------
    me_input_model : cobra.Model
        the corrected input metabolic model (as provided in preprocess.correct_inputs.correct_model)
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin
    psim_rib : pd.DataFrame
        formatted PSIM for ribosomal biogenesis
    modified_trna_transcript_c : macromolecules.RNA.tRNA
        output of create_trna()
    charged_trna_map : Dict[str, macromolecules.RNA.tRNA]
        output of creat_trna()
    ub_args : Dict[str, Any]
        set of variables generated by this ubiquitin.express_ubiquitin to be used throughout model building
    compress_mrna : bool
        whether to condense elongation, processing, and nuclear export reactions into a single reaction
    stochastic : bool
        Whether a potentially stochastic output should be stochastic, or choose a default behavior instead
    seeds : List[int]
        A list of seeds for if stochastic is set to True
    """

    # small ribosome proteins--------------------------------------------------------------------------------
    rs_ids = mach.rs['HGNC ID (gene)'].tolist()
    RPS27A_HGNC = 'HGNC:10417'
    rs_ids.remove(RPS27A_HGNC)  # RPS27A contains a ubiquitin monomer
    rs_expression_reactions, rs_protein_metabolites = list(), list()
    seed_idx = 0
    for i in rs_ids:
        gene_info = gene_information.generate_from_psim(hgnc_id=i, psim=psim_rib,
                                                        machinery_list=list(), reactions=None,
                                                        nonmachinery_locations=['n', 'c'],
                                                        stochastic=stochastic, seed=seeds[seed_idx])
        gene_info = func.convert_gi(gene_info, non_machinery=dict())
        mrna_expression_reactions, mrna_transcript_c, mrna_deg_proxy = build_mrna.get_mrna_expression_reactions(
            gene_info=gene_info, model_metabolites=model_metabolites, compress_mrna=compress_mrna)
        protein_expression_reactions, protein_metabolites = build_protein.get_protein_expression_reactions(gene_info, 
                                                                                                           mrna_transcript_c,
                                                                                                           mrna_deg_proxy,
                                                                                                           ub_args, 
                                                                                                            modified_trna_transcript_c=modified_trna_transcript_c, 
                                                                                                            charged_trna_map=charged_trna_map, model_metabolites=model_metabolites)
        protein_expression_reactions = [r for r in protein_expression_reactions if
                                        r.id.split('_')[-1] not in exclude]  # no nuclear degradation
        rs_expression_reactions += mrna_expression_reactions + protein_expression_reactions
        rs_protein_metabolites += protein_metabolites
        seed_idx += 1

    # large ribosome proteins--------------------------------------------------------------------------------
    rl_ids = mach.rl['HGNC ID (gene)'].tolist()
    RPL40_HGNC = 'HGNC:12458'
    rl_ids.remove(RPL40_HGNC)  # RPL40 contains a ubiquitin monomer
    rl_expression_reactions, rl_protein_metabolites = list(), list()
    for i in rl_ids:
        gene_info = gene_information.generate_from_psim(hgnc_id=i, psim=psim_rib,
                                                        machinery_list=list(), reactions=None,
                                                        nonmachinery_locations=['n', 'c'],
                                                        stochastic=stochastic, seed=seeds[seed_idx])
        gene_info = func.convert_gi(gene_info, non_machinery=dict())
        mrna_expression_reactions, mrna_transcript_c, mrna_deg_proxy = build_mrna.get_mrna_expression_reactions(
            gene_info=gene_info, model_metabolites=model_metabolites, compress_mrna=compress_mrna)
        protein_expression_reactions, protein_metabolites = build_protein.get_protein_expression_reactions(gene_info, 
                                                                                                           mrna_transcript_c,
                                                                                                           mrna_deg_proxy,
                                                                                                           ub_args, 
                                                                                                           modified_trna_transcript_c=modified_trna_transcript_c, 
                                                                                                           charged_trna_map=charged_trna_map,
                                                                                                           model_metabolites=model_metabolites, 
                                                                                                           )
        protein_expression_reactions = [r for r in protein_expression_reactions if
                                        r.id.split('_')[-1] not in exclude]  # no nuclear degradation
        rl_expression_reactions += mrna_expression_reactions + protein_expression_reactions
        rl_protein_metabolites += protein_metabolites
        seed_idx += 1

    # ubiquitin fusions--------------------------------------------------------------------------------------
    # RPS27A
    to_add, folded_protein_c, folded_protein_n = cleave_ub(hgnc_id=RPS27A_HGNC, model_metabolites=model_metabolites, psim_rib=psim_rib, ub_args=ub_args,
                                                           compress_mrna=compress_mrna, charged_trna_map=charged_trna_map, modified_trna_transcript_c=modified_trna_transcript_c,  
                                                           stochastic=stochastic,
                                                           seeds=seeds[seed_idx: seed_idx + 2])
    seed_idx += 2
    rs_expression_reactions += to_add
    rs_protein_metabolites += [folded_protein_c, folded_protein_n]
    # RPL40
    to_add, folded_protein_c, folded_protein_n = cleave_ub(hgnc_id=RPL40_HGNC,  model_metabolites=model_metabolites, psim_rib=psim_rib, ub_args=ub_args,
                                                           compress_mrna=compress_mrna, charged_trna_map=charged_trna_map, modified_trna_transcript_c=modified_trna_transcript_c, 
                                                           stochastic=stochastic,
                                                           seeds=seeds[seed_idx: seed_idx + 2])
    rl_expression_reactions += to_add
    for r in rl_expression_reactions + rs_expression_reactions:
        r.ribosome_biogenesis = True
    rl_protein_metabolites += [folded_protein_c, folded_protein_n]

    return rs_expression_reactions, rs_protein_metabolites, rl_expression_reactions, rl_protein_metabolites


def build_rrna5s_reactions(rpl5_n, rpl11_n, model_metabolites):
    # TRANSCRIPTION - basically emulates Transcript.transcript_elongation reaction
    pre_rrna5s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='pre_5s', seq=pre_rrna_5s_seq, compartment='n')
    rrna5s_transcription = pre_rrna5s_n.synthesize(id_='TRANSCRIPTION_PRE_RRNA5s')
    rrna5s_transcription.gene_reaction_rule = ' and '.join(mach.rnap3_transcription_machinery)

    # PROCESSING - mature rrna (3->5' exonucleolytic cleave of last 24 bases) and complex formation with RPL5/RPL11
    rrna5s_processing = ExpressionReaction(subsystem='Complex_Formation', ribosome_biogenesis=True,
                                           id='PROCESSING_RRNA5s')
    rrna5s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='5s', seq=rrna_5s_seq, compartment='n')
    deg_base_counts = dict()
    for k, v in pre_rrna5s_n.base_counts.items():
        deg_base_counts[k] = v - rrna5s_n.base_counts[k]

    rxn = dict()
    rxn[model_metabolites.h2o_n] = -sum(deg_base_counts.values())  # no -1 because all bonds 5'-most bond cleave
    rxn[model_metabolites.h_n] = sum(deg_base_counts.values())
    for k, v in model_metabolites.nmp_map_n.items():
        rxn[v] = deg_base_counts[k]

    rrna5s_complex_n = RibosomalComplex(metabolites=[rrna5s_n, rpl5_n, rpl11_n], complex_id='rrna5s')

    rxn[pre_rrna5s_n], rxn[rpl5_n], rxn[rpl11_n] = -1, -1, -1
    rxn[rrna5s_complex_n] = 1
    rrna5s_processing.add_metabolites(rxn)
    rrna5s_processing.gene_reaction_rule = mach.REXO5

    # TRANSPORT - will be transported as pre60s later, but make an rrna5s cytoplasmic for degradation, as
    # ribosome dissociates in cytoplasm
    # must add nucleocytoplasmic export via ran gtp: https://www.sciencedirect.com/science/article/pii/S0171933504702575?via%3Dihub
    rrna5s_c = rrna5s_n.change_compartment('c')

    # Degradation
    rrna5s_degradation = rrna5s_c.exonucleolytic_degradation(reaction_name='5s_rRNA')
    rrna5s_degradation.gene_reaction_rule = ' and '.join(mach.exosome['HGNC ID (gene)'].tolist())

    rrna5s_reactions = [rrna5s_transcription, rrna5s_processing, rrna5s_degradation]

    return rrna5s_reactions, rrna5s_complex_n, rrna5s_c


# In[5]:


# ets_5_frag1 is from 5' end of 47s to A' site
# ets_5_frag2 is from A' to 18s
# ets_5_frag3 is from A' to A0
# ets_5_frag4 is from A0 site to site 1 (start of 18s)
# its_1_frag1_seq is between site E and teh conserved stall location of RRP6
# its_1_frag2_seq is less than E site (some degradation) + a polyA/U tail

def build_other_rrna_reactions(rrna5s_complex_n, rs_protein_metabolites, rl_protein_metabolites,
                               rpl5_n, rpl11_n, rrna5s_c, model_metabolites, reversible_complex_formation=False):
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6315592/ figure 2

    # 47s transcription------------------------------------------------------------------------------------
    rrna_47s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='47s', seq=rrna_47s_seq, compartment='n')
    rrna_47s_transcription = rrna_47s_n.synthesize(id_='TRANSCRIPTION_RRNA_47s')
    rrna_47s_transcription.gene_reaction_rule = ' and '.join(mach.rnap1['HGNC ID (gene)'].tolist() + mach.rnap1_tfs)

    # 45s formation------------------------------------------------------------------------------------
    ets_5_frag1_seq = ets_5_seq[:A_prime_index]
    ets_5_frag2_seq = ets_5_seq[A_prime_index:]
    rrna_45s_seq = rrna_47s_seq[A_prime_index:rrna_47s_seq.index(rrna_28s_seq) + len(rrna_28s_seq)]

    ets_5_frag1_n = RNA_fragment(model_metabolites=model_metabolites, metabolite_name='5_frag1', seq=ets_5_frag1_seq, compartment='n', fragment_type='ets')  # ets5
    ets_3_n = RNA_fragment(model_metabolites=model_metabolites, metabolite_name='3', seq=ets_3_seq, compartment='n', triphosphate=False, fragment_type='ets')
    rrna_45s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='45s', seq=rrna_45s_seq, compartment='n', triphosphate=False)

    rrna_45s_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                            id='FORMATION_RRNA_45s')

    rxn = dict()
    rxn[model_metabolites.h2o_n] = -2  # 2 endonuclolytic cleavage events to go from 47s to 45s
    rxn[rrna_47s_n], rxn[rrna_45s_n], rxn[ets_3_n], rxn[ets_5_frag1_n] = -1, 1, 1, 1
    rrna_45s_formation.add_metabolites(rxn)
    rrna_45s_formation.gene_reaction_rule = ' and '.join(mach.UTP10 + mach.RNASEN)

    ets_3_degradation = ets_3_n.exonucleolytic_degradation(reaction_name='ets_3_rRNA', update=True)
    ets_5_frag1_degradation = ets_5_frag1_n.exonucleolytic_degradation(reaction_name='ets_5_frag1_rRNA', update=True)
    # 45S-->30S + 32.5S------------------------------------------------------------------------------------
    idx_30s = rrna_45s_seq.index(rrna_18s_seq) + len(rrna_18s_seq)
    rrna_30s_seq = ets_5_frag2_seq + rrna_18s_seq + rrna_45s_seq[idx_30s:idx_30s + site_2_index]
    rrna_32_5s_seq = rrna_45s_seq[idx_30s + site_2_index:]

    rrna_30s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='30s', seq=rrna_30s_seq, compartment='n', triphosphate=False)
    rrna_32_5s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='32_5s', seq=rrna_32_5s_seq, compartment='n', triphosphate=False)

    rrna_30s_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                            id='FORMATION_RRNA_30s_32_5s')
    rxn = dict()
    rxn[model_metabolites.h2o_n] = -1  # endonuclolytic cleavage event at site 2
    rxn[rrna_45s_n], rxn[rrna_30s_n], rxn[rrna_32_5s_n] = -1, 1, 1
    rrna_30s_formation.add_metabolites(rxn)
    rrna_30s_formation.gene_reaction_rule = ''  # mach.RMRP[0]<--disregarded for now bc ribozyme

    # 26s formation------------------------------------------------------------------------------------
    rrna_26s_seq = ets_5_frag2_seq[A_0_index:] + rrna_18s_seq + its_1_seq[:site_2_index]
    ets_5_frag3_seq = ets_5_frag2_seq[:A_0_index]

    rrna_26s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='26s', seq=rrna_26s_seq, compartment='n', triphosphate=False)
    ets_5_frag3_n = RNA_fragment(model_metabolites=model_metabolites, metabolite_name='5_frag3', seq=ets_5_frag3_seq, compartment='n', triphosphate=False,
                                 fragment_type='ets')

    rrna_26s_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                            id='FORMATION_RRNA_26s')

    rxn = dict()
    rxn[model_metabolites.h2o_n] = -1  # endonuclolytic cleavage event at site 2
    rxn[rrna_30s_n], rxn[rrna_26s_n], rxn[ets_5_frag3_n] = -1, 1, 1
    rrna_26s_formation.add_metabolites(rxn)
    rrna_26s_formation.gene_reaction_rule = mach.UTP23[0]

    ets_5_frag3_degradation = ets_5_frag3_n.exonucleolytic_degradation(reaction_name='ets_5_frag3_rRNA', update=True)

    # 21S formation------------------------------------------------------------------------------------

    rrna_21s_seq = rrna_18s_seq + its_1_seq[:site_2_index]
    ets_5_frag4_seq = ets_5_frag2_seq[A_0_index:]

    rrna_21s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='21s', seq=rrna_21s_seq, compartment='n', triphosphate=False)
    ets_5_frag4_n = RNA_fragment(model_metabolites=model_metabolites, metabolite_name='5_frag4', seq=ets_5_frag4_seq, fragment_type='ets', compartment='n',
                                 triphosphate=False)

    rrna_21s_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                            id='FORMATION_RRNA_21s')
    rxn = dict()

    rxn[model_metabolites.h2o_n] = -1  # endonuclolytic cleavage event at site 1
    rxn[rrna_26s_n], rxn[rrna_21s_n], rxn[ets_5_frag4_n] = -1, 1, 1
    rrna_21s_formation.add_metabolites(rxn)
    rrna_21s_formation.gene_reaction_rule = mach.UTP24[0]

    ets_5_frag4_degradation = ets_5_frag4_n.exonucleolytic_degradation(reaction_name='ets_5_frag4_rRNA', update=True)

    # 21SC formation------------------------------------------------------------------------------------
    rrna_21sc_seq = rrna_18s_seq + its_1_seq[:conserved_stall_idx]
    rrna_21sc_n = rRNA(model_metabolites=model_metabolites, metabolite_name='21sc', seq=rrna_21sc_seq, compartment='n', triphosphate=False)

    deg_seq = its_1_seq[conserved_stall_idx: site_2_index]
    base_counts_deg, elements_deg = func.get_base_counts_and_elements(deg_seq, model_metabolites=model_metabolites)

    rrna_21sc_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                             id='FORMATION_RRNA_21sc')
    rxn = dict()
    rxn[rrna_21s_n], rxn[rrna_21sc_n] = -1, 1
    # exonucleolytic cleavage
    for k, v in model_metabolites.nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[model_metabolites.h2o_n] = -len(deg_seq)
    rxn[model_metabolites.h_n] = len(deg_seq)
    rrna_21sc_formation.add_metabolites(rxn)
    rrna_21sc_formation.gene_reaction_rule = \
        mach.exosome[mach.exosome['Approved symbol'] == 'EXOSC10']['HGNC ID (gene)'].tolist()[0]

    # 18se formation------------------------------------------------------------------------------------
    rrna_18se_seq = rrna_18s_seq + its_1_seq[:e_index]
    its_1_frag1_seq = its_1_seq[e_index:conserved_stall_idx]
    rrna_18se_n = rRNA(model_metabolites=model_metabolites, metabolite_name='18se', seq=rrna_18se_seq, compartment='n', triphosphate=False)
    its_1_frag1_n = RNA_fragment(model_metabolites=model_metabolites, metabolite_name='1_frag1', seq=its_1_frag1_seq, fragment_type='its', compartment='n',
                                 triphosphate=False)

    rrna_18se_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                             id='FORMATION_RRNA_18se')
    # endonuclolytic cleavage event at site E
    rrna_18se_formation.add_metabolites({model_metabolites.h2o_n: -1, rrna_21sc_n: -1, rrna_18se_n: 1, its_1_frag1_n: 1})
    rrna_18se_formation.gene_reaction_rule = mach.UTP24[0]

    its_1_frag1_degradation = its_1_frag1_n.exonucleolytic_degradation(reaction_name='its_1_frag1_rRNA', update=True)

    # 18se nuclear processing------------------------------------------------------------------------------------
    rrna_18se_processed_seq = rrna_18se_seq[:-int(0.75 * e_index)]  # degradation of 60/80 bps of ITS1 by PARN
    rrna_18se_processed_seq += 'U' * int(0.125 * e_index) + 'A' * int(0.125 * e_index)  # polyU by PAPD5
    deg_seq = rrna_18se_seq[-int(0.75 * e_index):]

    rrna_18se_processed_n = rRNA(model_metabolites=model_metabolites, metabolite_name='18se_processed', seq=rrna_18se_processed_seq, compartment='n', triphosphate=False)
    base_counts_deg, elements_deg = func.get_base_counts_and_elements(deg_seq, model_metabolites=model_metabolites)

    rrna_18se_processing = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                              id='PROCESSING_RRNA_18se')
    rxn = dict()
    rxn[rrna_18se_n], rxn[rrna_18se_processed_n] = -1, 1
    # exonucleolytic cleavage
    for k, v in model_metabolites.nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[model_metabolites.h2o_n] = -len(deg_seq)
    rxn[model_metabolites.h_n] = len(deg_seq)
    # polyU/A tail synthesis<--unsure why ppi_n has -1 but it mass balances
    rxn[model_metabolites.atp_n], rxn[model_metabolites.ntp_map_n['U']], rxn[model_metabolites.ppi_n] = -int(0.125 * e_index), -int(0.125 * e_index), int(
        0.25 * e_index) - 1

    rrna_18se_processing.add_metabolites(rxn)
    rrna_18se_processing.gene_reaction_rule = ' and '.join(mach.PARN + mach.PAPD5)

    # pre40s complex------------------------------------------------------------------------------------
    metabolites = [m for m in rs_protein_metabolites if m.compartment == 'n'] + [rrna_18se_processed_n]
    pre40s_complex_n = RibosomalComplex(metabolites=metabolites, complex_id='pre40s')
    pre40s_complex_formation = pre40s_complex_n.form_complex(reversible=reversible_complex_formation)
    pre40s_complex_formation.lower_bound = 0
    pre40s_complex_formation.gene_reaction_rule = ' and '.join(mach.pre40s_rbfs)

    # pre40s nucleocytoplasmic export-----------------------------------------------------------------------
    pre40s_complex_c = pre40s_complex_n.change_compartment('c')

    pre40s_transport = ExpressionReaction(subsystem='Complex_Formation', ribosome_biogenesis=True,
                                          id='pre40s_NUCLEAR_EXPORTtn')
    pre40s_transport.name = 'pre40s nuclear export'
    rxn = {pre40s_complex_n: -1, pre40s_complex_c: 1}
    # gtp hydrolysis
    rxn[model_metabolites.ntp_map_c['G']], rxn[model_metabolites.h2o_c], rxn[model_metabolites.ndp_map_c['G']], rxn[model_metabolites.pi_c], rxn[
        model_metabolites.h_c] = -1, -1, 1, 1, 1
    pre40s_transport.add_metabolites(rxn)
    pre40s_transport.gene_reaction_rule = ' and '.join(mach.tfiiia + mach.RAN + mach.XPO1)

    # 18s/mature 40s formation------------------------------------------------------------------------------------
    its_1_frag2_seq = its_1_seq[:int(0.25 * e_index) + 1] + 'U' * int(0.125 * e_index) + 'A' * int(0.125 * e_index)
    rrna_18s_c = rRNA(model_metabolites=model_metabolites, metabolite_name='18s', seq=rrna_18s_seq, compartment='c', triphosphate=False)
    its_1_frag2_c = RNA_fragment(model_metabolites=model_metabolites, metabolite_name='1_frag2', seq=its_1_frag2_seq, fragment_type='its', compartment='c',
                                 triphosphate=False)

    rrna_18s_formation = ExpressionReaction(subsystem='Complex_Formation', ribosome_biogenesis=True,
                                            id='40s_MATURATION')
    # endonuclolytic cleavage event at site 3
    metabolites = [m for m in rs_protein_metabolites if m.compartment == 'c'] + [rrna_18s_c]
    forty_s_complex_c = RibosomalComplex(metabolites=metabolites, complex_id='40s')

    rrna_18s_formation.add_metabolites({model_metabolites.h2o_n: -1, pre40s_complex_c: -1, forty_s_complex_c: 1,
                                        its_1_frag2_c: 1})
    rrna_18s_formation.gene_reaction_rule = mach.NOB1[0]

    its_1_frag2_degradation = its_1_frag2_c.exonucleolytic_degradation(reaction_name='its_1_frag2_rRNA', update=True)

    # 18s degradation------------------------------------------------------------------------------------
    rrna_18s_degradation = rrna_18s_c.exonucleolytic_degradation(reaction_name='18s_rrna_degradation', update=True)

    # 32S formation------------------------------------------------------------------------------------
    deg_seq = its_1_seq[site_2_index:]
    rrna_32s_seq = rrna_32_5s_seq[len(deg_seq):]

    rrna_32s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='32s', seq=rrna_32s_seq, compartment='n', triphosphate=False)
    base_counts_deg, elements_deg = func.get_base_counts_and_elements(deg_seq, model_metabolites=model_metabolites)

    rrna_32s_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                            id='FORMATION_RRNA_32s')
    rxn = dict()
    rxn[rrna_32_5s_n], rxn[rrna_32s_n] = -1, 1
    # exonucleolytic cleavage
    for k, v in model_metabolites.nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[model_metabolites.h2o_n] = -len(deg_seq)
    rxn[model_metabolites.h_n] = len(deg_seq)

    rrna_32s_formation.add_metabolites(rxn)
    rrna_32s_formation.gene_reaction_rule = mach.lariat_machinery["5' Degradation"][0]

    # 32s-->12s + 28.5s------------------------------------------------------------------------------------

    rrna_12s_seq = rrna_5_8s_seq + its_2_seq[:site_4_index]
    rrna_28_5s_seq = its_2_seq[site_4_index:] + rrna_28s_seq
    rrna_12s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='12s', seq=rrna_12s_seq, compartment='n', triphosphate=False)
    rrna_28_5s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='28_5s', seq=rrna_28_5s_seq, compartment='n', triphosphate=False)

    rrna_12s_28_5s_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                                  id='FORMATION_RRNA_12s_28_5s')
    rxn = dict()
    rxn[model_metabolites.h2o_n] = -1  # endonuclolytic cleavage event at site 4
    rxn[rrna_32s_n], rxn[rrna_12s_n], rxn[rrna_28_5s_n] = -1, 1, 1
    rrna_12s_28_5s_formation.add_metabolites(rxn)
    rrna_12s_28_5s_formation.gene_reaction_rule = mach.LAS1[0]

    # 28s formation------------------------------------------------------------------------------------
    deg_seq = rrna_28_5s_seq[:rrna_28_5s_seq.index(rrna_28s_seq)]

    rrna_28s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='28s', seq=rrna_28s_seq, compartment='n', triphosphate=False)
    base_counts_deg, elements_deg = func.get_base_counts_and_elements(deg_seq, model_metabolites=model_metabolites)

    rrna_28s_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                            id='FORMATION_RRNA_28s')

    rxn = dict()
    rxn[rrna_28_5s_n], rxn[rrna_28s_n] = -1, 1
    # exonucleolytic cleavage
    for k, v in model_metabolites.nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[model_metabolites.h2o_n] = -len(deg_seq)
    rxn[model_metabolites.h_n] = len(deg_seq)

    rrna_28s_formation.add_metabolites(rxn)
    rrna_28s_formation.gene_reaction_rule = mach.lariat_machinery["5' Degradation"][0]

    # 7s formation------------------------------------------------------------------------------------
    deg_seq = its_2_seq[SEVEN_S_IDX: site_4_index]
    rrna_7s_seq = rrna_5_8s_seq + its_2_seq[:SEVEN_S_IDX]
    rrna_7s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='7s', seq=rrna_7s_seq, compartment='n', triphosphate=False)
    base_counts_deg, elements_deg = func.get_base_counts_and_elements(deg_seq, model_metabolites=model_metabolites)

    rrna_7s_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                           id='FORMATION_RRNA_7s')

    rxn = dict()
    rxn[rrna_12s_n], rxn[rrna_7s_n] = -1, 1
    # exonucleolytic cleavage
    for k, v in model_metabolites.nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[model_metabolites.h2o_n] = -len(deg_seq)
    rxn[model_metabolites.h_n] = len(deg_seq)

    rrna_7s_formation.add_metabolites(rxn)
    rrna_7s_formation.gene_reaction_rule = ' and '.join(mach.DIS3 + mach.ISG20L2)

    # 5.8s+40 formation------------------------------------------------------------------------------------
    deg_seq = its_2_seq[FIVE_EIGHT_PLUS_FORTY_IDX: SEVEN_S_IDX]
    rrna_5_8s_plus_40_seq = rrna_5_8s_seq + its_2_seq[:FIVE_EIGHT_PLUS_FORTY_IDX]
    rrna_5_8s_plus_40_n = rRNA(model_metabolites=model_metabolites, metabolite_name='5_8s_plus_40', seq=rrna_5_8s_plus_40_seq, compartment='n', triphosphate=False)
    base_counts_deg, elements_deg = func.get_base_counts_and_elements(deg_seq, model_metabolites=model_metabolites)

    rrna_5_8s_plus_40_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                                     id='FORMATION_RRNA_5_8s_plus_40')

    rxn = dict()
    rxn[rrna_7s_n], rxn[rrna_5_8s_plus_40_n] = -1, 1
    # exonucleolytic cleavage
    for k, v in model_metabolites.nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[model_metabolites.h2o_n] = -len(deg_seq)
    rxn[model_metabolites.h_n] = len(deg_seq)

    rrna_5_8s_plus_40_formation.add_metabolites(rxn)
    rrna_5_8s_plus_40_formation.gene_reaction_rule = ' and '.join(mach.DIS3 + mach.ISG20L2)

    # 6s formation------------------------------------------------------------------------------------

    deg_seq = its_2_seq[SIX_S_IDX: FIVE_EIGHT_PLUS_FORTY_IDX]
    rrna_6s_seq = rrna_5_8s_seq + its_2_seq[:SIX_S_IDX]
    rrna_6s_n = rRNA(model_metabolites=model_metabolites, metabolite_name='6s', seq=rrna_6s_seq, compartment='n', triphosphate=False)
    base_counts_deg, elements_deg = func.get_base_counts_and_elements(deg_seq, model_metabolites=model_metabolites)

    rrna_6s_formation = ExpressionReaction(subsystem='rRNA_expression', ribosome_biogenesis=True,
                                           id='FORMATION_RRNA_6s')

    rxn = dict()
    rxn[rrna_5_8s_plus_40_n], rxn[rrna_6s_n] = -1, 1
    # exonucleolytic cleavage
    for k, v in model_metabolites.nmp_map_n.items():
        rxn[v] = base_counts_deg[k]
    rxn[model_metabolites.h2o_n] = -len(deg_seq)
    rxn[model_metabolites.h_n] = len(deg_seq)

    rrna_6s_formation.add_metabolites(rxn)
    rrna_6s_formation.gene_reaction_rule = \
        mach.exosome[mach.exosome['Approved symbol'] == 'EXOSC10']['HGNC ID (gene)'].tolist()[0]

    # pre60s complex formation------------------------------------------------------------------------------------
    rl_2 = list(set(rl_protein_metabolites).difference([rpl5_n, rpl11_n]))
    metabolites = [m for m in rl_2 if m.compartment == 'n'] + [rrna_28s_n, rrna_6s_n, rrna5s_complex_n]
    pre60s_complex_n = RibosomalComplex(metabolites=metabolites, complex_id='pre60s')
    pre60s_complex_formation = pre60s_complex_n.form_complex(reaction_id='pre60s',
                                                             reversible=reversible_complex_formation)
    # add a gtp hydrolysis to the complex formation: https://www.embopress.org/doi/full/10.15252/embj.2018100278
    rxn = dict()
    # gtp hydrolysis
    rxn[model_metabolites.ntp_map_n['G']], rxn[model_metabolites.h2o_n], rxn[model_metabolites.gdp_n], rxn[model_metabolites.pi_n], rxn[model_metabolites.h_n] = -1, -1, 1, 1, 1
    pre60s_complex_formation.add_metabolites(rxn)
    pre60s_complex_formation.lower_bound = 0
    pre60s_complex_formation.gene_reaction_rule = ' and '.join(mach.pre60s_rbfs)

    # pre60s nucleocytoplasmic export-----------------------------------------------------------------------
    pre60s_complex_c = pre60s_complex_n.change_compartment('c')

    pre60s_transport = ExpressionReaction(subsystem='Complex_Formation', ribosome_biogenesis=True,
                                          id='pre60s_NUCLEAR_EXPORTtn')
    pre60s_transport.name = 'pre60s nuclear export'
    rxn = {pre60s_complex_n: -1, pre60s_complex_c: 1}
    # gtp hydrolysis
    rxn[model_metabolites.ntp_map_c['G']], rxn[model_metabolites.h2o_c], rxn[model_metabolites.ndp_map_c['G']], rxn[model_metabolites.pi_c], rxn[
        model_metabolites.h_c] = -1, -1, 1, 1, 1
    pre60s_transport.add_metabolites(rxn)
    pre60s_transport.gene_reaction_rule = ' and '.join(mach.tfiiia + mach.RAN + mach.XPO1)

    # 5.8s/mature 60s formation------------------------------------------------------------------------------------
    deg_seq = its_2_seq[:SIX_S_IDX]
    rrna_5_8s_c = rRNA(model_metabolites=model_metabolites, metabolite_name='5_8s', seq=rrna_5_8s_seq, compartment='c', triphosphate=False)
    base_counts_deg, elements_deg = func.get_base_counts_and_elements(deg_seq, model_metabolites=model_metabolites)

    rrna_28s_c = rrna_28s_n.change_compartment('c')

    metabolites = [m for m in rl_2 if m.compartment == 'c'] + [rrna_28s_c, rrna_5_8s_c, rrna5s_c]
    sixty_s_complex_c = RibosomalComplex(metabolites=metabolites, complex_id='60s')
    rrna_5_8s_formation = ExpressionReaction(subsystem='Complex_Formation', ribosome_biogenesis=True,
                                             id='60s_maturation')
    rxn = dict()

    rxn[pre60s_complex_c], rxn[sixty_s_complex_c] = -1, 1
    # exonucleolytic cleavage
    for k, v in model_metabolites.nmp_map_c.items():
        rxn[v] = base_counts_deg[k]
    rxn[model_metabolites.h2o_c] = -len(deg_seq)
    rxn[model_metabolites.h_c] = len(deg_seq)
    rrna_5_8s_formation.add_metabolites(rxn)
    rrna_5_8s_formation.gene_reaction_rule = mach.ERI1[0]

    # 5.8s and 28s degradation------------------------------------------------------------------------------------
    rrna_28s_degradation = rrna_28s_c.exonucleolytic_degradation(reaction_name='28s_rrna', update=True)
    rrna_5_8s_degradation = rrna_5_8s_c.exonucleolytic_degradation(reaction_name='5_8s_rrna', update=True)

    # ------------------------------------------------------------------------------------
    all_reactions = [rrna_47s_transcription, rrna_45s_formation, ets_3_degradation, ets_5_frag1_degradation,
                     rrna_30s_formation, rrna_26s_formation, ets_5_frag3_degradation, rrna_21s_formation,
                     ets_5_frag4_degradation, rrna_21sc_formation, rrna_18se_formation, its_1_frag1_degradation,
                     rrna_18se_processing, pre40s_complex_formation, pre40s_transport,
                     rrna_18s_formation, rrna_18s_degradation, its_1_frag2_degradation,
                     rrna_32s_formation, rrna_12s_28_5s_formation, rrna_28s_formation,
                     rrna_7s_formation, rrna_5_8s_plus_40_formation, rrna_6s_formation, pre60s_complex_formation,
                     pre60s_transport, rrna_5_8s_formation, rrna_28s_degradation, rrna_5_8s_degradation]

    mature_ribosomal_precomplexes = [forty_s_complex_c, sixty_s_complex_c]
    mature_rrna_metabolites = [rrna_5_8s_c, rrna_28s_c, rrna_18s_c]

    return all_reactions, mature_ribosomal_precomplexes, mature_rrna_metabolites


def build_ribosome(model_metabolites, psim_rib, modified_trna_transcript_c, charged_trna_map, ub_args: Dict[str, Any], compress_mrna: bool, reversible_complex_formation: bool, stochastic: bool, seed: int):
    """Generate all ribosome biogenesis reactions

    Parameters
    ----------
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin
    psim_rib : pd.DataFrame
        formatted PSIM for ribosomal biogenesis
    modified_trna_transcript_c : macromolecules.RNA.tRNA
        output of create_trna()
    charged_trna_map : Dict[str, macromolecules.RNA.tRNA]
        output of creat_trna()
    ub_args : Dict[str, Any]
        set of variables generated by this ubiquitin.express_ubiquitin to be used throughout model building
    compress_mrna : bool
        whether to condense elongation, processing, and nuclear export reactions into a single reaction
    reversible_complex_formation : bool
        whether complex formation reactions are reversible
    stochastic : bool
        Whether a potentially stochastic output should be stochastic, or choose a default behavior instead
    seed : int
        A seed for if stochastic is set to True
    """
    with HiddenPrints():
        random.seed(seed)
        seeds = random.sample(range(0, int((2 ** 32 - 1))), k=87)
        rs_expression_reactions, rs_protein_metabolites, rl_expression_reactions, rl_protein_metabolites = build_ribosome_protein_expression_reactions(
            model_metabolites=model_metabolites, psim_rib=psim_rib, modified_trna_transcript_c=modified_trna_transcript_c, charged_trna_map=charged_trna_map,
            ub_args=ub_args, compress_mrna=compress_mrna, stochastic=stochastic, seeds=seeds)
    rpl5_n = [m for m in rl_protein_metabolites if m.id == 'HGNC:10360_folded_protein_n'][0]
    rpl11_n = [m for m in rl_protein_metabolites if m.id == 'HGNC:10301_folded_protein_n'][0]
    rrna5s_reactions, rrna5s_complex_n, rrna5s_c = build_rrna5s_reactions(rpl5_n, rpl11_n, model_metabolites)
    other_rrna_reactions, mature_ribosomal_precomplexes, mature_rrna_metabolites = build_other_rrna_reactions(
        rrna5s_complex_n, rs_protein_metabolites, rl_protein_metabolites,
        rpl5_n, rpl11_n, rrna5s_c, model_metabolites=model_metabolites,
        reversible_complex_formation=reversible_complex_formation)

    # ribosome complex formation
    ribosome_complex_c = RibosomalComplex(metabolites=mature_ribosomal_precomplexes, complex_id='mature_ribosome')
    ribosome_complex_formation = ribosome_complex_c.form_complex(reversible=reversible_complex_formation)
    # add a gtp hydrolysis to the complex formation: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5861459/
    rxn = dict()
    rxn[model_metabolites.ntp_map_c['G']], rxn[model_metabolites.h2o_c], rxn[model_metabolites.ndp_map_c['G']], rxn[model_metabolites.pi_c], rxn[
        model_metabolites.h_c] = -1, -1, 1, 1, 1
    ribosome_complex_formation.add_metabolites(rxn)

    #     ribosome_complex_formation.id = 'RIBOSOME_COMPLEX_FORMATIONc'
    ribosome_complex_formation.lower_bound = 0
    ribosome_complex_formation.gene_reaction_rule = ' and '.join(mach.eifs)

    all_reactions = rrna5s_reactions + other_rrna_reactions + [ribosome_complex_formation]
    # ribosome complex dissociation
    if reversible_complex_formation:
        ribosome_complex_dissociation = ExpressionReaction(subsystem='Complex_Degradation', ribosome_biogenesis=True,
                                                           id='RIBOSOME_COMPLEX_DISSOCIATIONc')
        ind_mets = rl_protein_metabolites + rs_protein_metabolites + mature_rrna_metabolites + [rrna5s_c]
        rxn = {m: -1 for m in ind_mets if m.compartment == 'c'}
        rxn[ribosome_complex_c] = 1
        ribosome_complex_dissociation.add_metabolites(rxn)
        all_reactions.append(ribosome_complex_dissociation)

    all_reactions += rs_expression_reactions + rl_expression_reactions

    return all_reactions, ribosome_complex_c
