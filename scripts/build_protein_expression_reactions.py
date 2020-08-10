#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *
from utils import *
from utils_2 import *
import build_mrna_expression_reactions as bm


# In[151]:


def translate_protein_cytosolic(gene_info):    

    # peptide bond formation: https://d1j63owfs0b5j3.cloudfront.net/pop-quiz/answerImage/Amino-Acid-1-popquiz.png
    # tRNA amino acide release: https://rnajournal.cshlp.org/content/14/8/1526/F1.expansion.html

    rxn = {charged_trna_map[aa_code]: -aa_count for aa_code, aa_count in gene_info.amino_acid_counts.items()} # tRNA consumption
    rxn[modified_trna_transcript_c] = gene_info.L_protein
    rxn[h2o_c] = -gene_info.L_protein # release of peptide from tRNA, addition of -OH to uncharged tRNA
    rxn[h_c] = gene_info.L_protein # release of peptide from tRNA, addition of -OH to uncharged tRNA
    rxn[h2o_c] += gene_info.L_protein - 1 # peptide bond formation (hydrolysis)
    
    # gtp hydrolysis per aa added
    rxn[ntp_map_c['G']] = -gene_info.L_protein 
    rxn[h2o_c] -= gene_info.L_protein
    rxn[ndp_map_c['G']] = gene_info.L_protein
    rxn[pi_c] = gene_info.L_protein
    rxn[h_c] += gene_info.L_protein
    

    rxn_c = rxn.copy()
    unfolded_protein_c = make_protein_metabolite(id_ = gene_info.hgnc_id + '_unfolded', 
                amino_acid_counts = gene_info.amino_acid_counts, L_protein = gene_info.L_protein,
                compartment = 'c')
    rxn_c[unfolded_protein_c] = 1
    
    translation_elongation = cobra.Reaction(gene_info.hgnc_id + '_CYTOSOLIC_TRANSLATION_ELONGATION')
    translation_elongation.subsytem = 'Protein_Expression'
    translation_elongation.add_metabolites(rxn_c)

    translation_elongation.gene_reaction_rule = ' and '.join(translation_efs + ['ribosome']) # GPRs

    return translation_elongation, unfolded_protein_c

def fold_protein_cytosolic(gene_info, unfolded_protein_c):
    # extending proteostasis network in the future would be good
    # will need to make sure inputs to each compartment-specific reactions are at the correct folding stage
    # e.g., mitochondria currently takes unfolded protein, and in future we may want it to take a partially folded
    
    folded_protein_c = unfolded_protein_c.copy()
    folded_protein_c.id = folded_protein_c.id.replace('unfolded', 'folded')
    rxn = {unfolded_protein_c: -1, folded_protein_c: 1}
    protein_folding = cobra.Reaction(gene_info.hgnc_id + '_CYTOSOLIC_PROTEIN_FOLDING')
    protein_folding.subsytem = 'Protein_Expression'
    
    if gene_info.L_protein > 100: #chaperone assisted for larger proteins - https://www.nature.com/articles/nature10317
        rxn = hydrolyze_atp(rxn, n_atp = gene_info.L_protein*proteolysis_translocation_atp_cost, compartment = 'c')
        protein_folding.gene_reaction_rule = ' and '.join(HSP40_c + HSP70_c) # GPRs
    
    
    protein_folding.add_metabolites(rxn)

    
    
    return protein_folding, folded_protein_c


# Ubiquitin expression

# In[83]:


# UBC
ubc_psim = psim_me[psim_me['HGNC_ID'] == 'HGNC:12468'] # UBC
ubc_psim['Location'] = 'c'
ubc_info = gene_information(metabolic_model = human_model, hgnc_id = ubc_psim['HGNC_ID'].values.tolist()[0], 
                         premrna_seq=ubc_psim['PREMRNA_SEQ'].values.tolist()[0], 
                            mrna_seq=ubc_psim['MRNA_SEQ'].values.tolist()[0], 
                            protein_seq=ubc_psim['PROTEIN_SEQ'].values.tolist()[0],
                            polyA_length = round(ubc_psim['POLYA_LENGTH'].values.tolist()[0]))
ubc_info.get_final_locations(human_model, final_locations=['c'])
ubc_mrna_expression_reactions = bm.mrna_expression(ubc_info)

# ubiquitin monomer
single_ubiquitin_sequence = ubc_info.protein_seq[:76]
monoub_aa_counts = {k: single_ubiquitin_sequence.count(k) for k in amino_acids}
L_monoub = len(single_ubiquitin_sequence)
n_ub_monomers = ubc_info.protein_seq.count(single_ubiquitin_sequence)
ub_c = make_protein_metabolite(id_ = 'ubiquitin_monomer', amino_acid_counts = monoub_aa_counts,
                               L_protein = L_monoub, compartment = 'c')

# monomerization from ubc polyub
# amino_acid_counts_ubc = {k: ubc_info.protein_seq.count(k) for k in amino_acids}
# L_ubc = len(ubc_info.protein_seq)

ubc_translation_reaction_cytosolic, ubc_c = translate_protein_cytosolic(ubc_info)

ubiquitin_monomerization_ubc = cobra.Reaction(ubc_info.hgnc_id + '_monomerization')
ubiquitin_monomerization_ubc.subsytem = 'Protein_Expression'
rxn = {ubc_c:-1, ub_c: n_ub_monomers, seq_amino_acid_map_c[ubc_info.protein_seq[n_ub_monomers*76:]]: 1, 
      h2o_c: -n_ub_monomers}
ubiquitin_monomerization_ubc.add_metabolites(rxn)
ubiquitin_monomerization_ubc.gene_reaction_rule = USP5[0]

# UBB
ubb_psim = psim_me[psim_me['HGNC_ID'] == 'HGNC:12463'] # UBB
ubb_psim['Location'] = 'c'
ubb_info = gene_information(metabolic_model = human_model, hgnc_id = ubb_psim['HGNC_ID'].values.tolist()[0], 
                         premrna_seq=ubb_psim['PREMRNA_SEQ'].values.tolist()[0], 
                            mrna_seq=ubb_psim['MRNA_SEQ'].values.tolist()[0], 
                            protein_seq=ubb_psim['PROTEIN_SEQ'].values.tolist()[0],
                            polyA_length = round(ubb_psim['POLYA_LENGTH'].values.tolist()[0]))
ubb_info.get_final_locations(human_model, final_locations=['c'])
ubb_mrna_expression_reactions = bm.mrna_expression(ubb_info)

# amino_acid_counts_ubb = {k: ubb_info.protein_seq.count(k) for k in amino_acids}
# L_ubb = len(ubb_info.protein_seq)

ubb_translation_reaction_cytosolic, ubb_c = translate_protein_cytosolic(ubb_info)

# monomerization from ubb polyub
n_ub_monomers = ubb_info.protein_seq.count(single_ubiquitin_sequence)
ubiquitin_monomerization_ubb = cobra.Reaction(ubb_info.hgnc_id + '_monomerization')
ubiquitin_monomerization_ubb.subsytem = 'Protein_Expression'
rxn = {ubb_c:-1, ub_c: n_ub_monomers, seq_amino_acid_map_c[ubb_info.protein_seq[n_ub_monomers*76:]]: 1, 
      h2o_c: -n_ub_monomers}
ubiquitin_monomerization_ubb.add_metabolites(rxn)
ubiquitin_monomerization_ubc.gene_reaction_rule = USP5[0]

# breakdown of the polyubiquitin cleaved from proteins in ubiquitin-proteasome pathway
polyub_aa_counts = {aa_code: aa_count*n_ub for aa_code, aa_count in monoub_aa_counts.items()}
polyub_c = make_protein_metabolite(id_ = 'cleaved_polyubiquitin_moiety', amino_acid_counts = polyub_aa_counts,
                               L_protein = L_monoub*n_ub, compartment = 'c')
ubiquitin_monomerization_polyub = cobra.Reaction('polyubiquitin_monomerization')
ubiquitin_monomerization_polyub.subsytem = 'Protein_Expression'
rxn = {polyub_c:-1, ub_c: n_ub, h2o_c: -(n_ub-1)}
ubiquitin_monomerization_polyub.add_metabolites(rxn)
ubiquitin_monomerization_polyub.gene_reaction_rule = USP5[0]

# nuclear import of ubiquitin
nuclear_import_ub_mono = cobra.Reaction('ubiquitin_monomer_nuclear_importtn')
nuclear_import_ub_mono.subsytem = 'Protein_Expression'
ub_n = ub_c.copy()
ub_n.id, ub_n.compartment = ub_n.id.replace('[c]', '[n]'), 'n'
nuclear_import_ub_mono.add_metabolites({ub_n: 1, ub_c: -1})
nuclear_import_ub_mono.lower_bound = -1000

# nuclear export of polyubiquitin moiety
nuclear_export_ub_poly = cobra.Reaction('polyubiquitin_moiety_nuclear_exporttn')
nuclear_export_ub_poly.subsytem = 'Protein_Expression'
polyub_n = polyub_c.copy()
polyub_n.id, polyub_n.compartment = polyub_n.id.replace('[c]', '[n]'), 'n'
nuclear_export_ub_poly.add_metabolites({polyub_n: -1, polyub_c: 1})
nuclear_export_ub_poly.lower_bound = -1000

# degradation
degradation_ub = cobra.Reaction('ubiquitin_monomer_degradation')
degradation_ub.subsytem = 'Protein_Expression'
rxn = {seq_amino_acid_map_c[aa_code]: aa_counts for aa_code, aa_counts in monoub_aa_counts.items()}
rxn[ub_c] = -1
rxn[h2o_c] =  -(L_monoub-1)
# atp hydrolysis for translocation/unfolding by 26S - known 1 ATP per 2 residues - https://www.nature.com/articles/s41586-018-0736-4
rxn = hydrolyze_atp(rxn, n_atp = L_monoub/2, compartment = 'c')

degradation_ub.add_metabolites(rxn)
degradation_ub.gene_reaction_rule = ' and '.join(proteasome_machinery)

ub_reactions = ubc_mrna_expression_reactions + [ubc_translation_reaction_cytosolic] 
ub_reactions += ubb_mrna_expression_reactions + [ubb_translation_reaction_cytosolic]
ub_reactions += [ubiquitin_monomerization_ubc, ubiquitin_monomerization_ubb, ubiquitin_monomerization_polyub, degradation_ub]
ub_reactions += [nuclear_import_ub_mono, nuclear_export_ub_poly]


# # Degradation (Ubiquitin-Proteasome)

# In[168]:


def protein_polyubiquitination(gene_info, protein_metabolite, compartment):
    
    polyu_protein_aa_counts = gene_info.amino_acid_counts.copy()
    for aa_code,aa_counts in monoub_aa_counts.items():
        if aa_code in polyu_protein_aa_counts:
            polyu_protein_aa_counts[aa_code] += aa_counts*n_ub
        else: 
            polyu_protein_aa_counts[aa_code] = aa_counts*n_ub
    
    if compartment == 'c':
        if protein_metabolite.compartment != 'c':
            raise ValueError('Compartment mismatch for polyubiquitination')
    
        polyubiquitinate_protein = cobra.Reaction(protein_metabolite.id + '_CYTOPLASMIC_POLYUBIQUITINATION')
        polyubiquitinate_protein.subsytem = 'Protein_Expression'

        polyub_protein_c = make_protein_metabolite(id_ = protein_metabolite.id + '_polyub', 
                           amino_acid_counts = polyu_protein_aa_counts, L_protein = gene_info.L_protein + (L_monoub*n_ub),
                           compartment = 'c') 
        rxn = {protein_metabolite: -1, ub_c: -n_ub, polyub_protein_c:1, h2o_c: n_ub}
        # 1 ATP hydrolysis per ubiquitin monomer added (https://link.springer.com/article/10.1007/s10637-020-00894-6)
        rxn = hydrolyze_atp(rxn, n_atp = n_ub, compartment = 'c')

        polyubiquitinate_protein.add_metabolites(rxn)
        polyubiquitinate_protein.gene_reaction_rule = ' and '.join(UB_ligases_c)
        return polyubiquitinate_protein, polyub_protein_c
    elif compartment == 'n':
        if protein_metabolite.compartment != 'n':
            raise ValueError('Compartment mismatch for polyubiquitination')

        polyubiquitinate_protein = cobra.Reaction(protein_metabolite.id + '_NUCLEAR_POLYUBIQUITINATION')
        polyubiquitinate_protein.subsytem = 'Protein_Expression'

        polyub_protein_n = make_protein_metabolite(id_ = protein_metabolite.id + '_polyub', 
                           amino_acid_counts = polyu_protein_aa_counts, L_protein = gene_info.L_protein + (L_monoub*n_ub),
                           compartment = 'n') 
        rxn = {protein_metabolite: -1, ub_n: -n_ub, polyub_protein_n: 1, h2o_n: n_ub}
        # 1 ATP hydrolysis per ubiquitin monomer added (https://link.springer.com/article/10.1007/s10637-020-00894-6)
        rxn = hydrolyze_atp(rxn, n_atp = n_ub, compartment = 'n')


        polyubiquitinate_protein.add_metabolites(rxn)
        polyubiquitinate_protein.gene_reaction_rule = ' and '.join(UB_ligases_n)
        
        return polyubiquitinate_protein, polyub_protein_n

    else:
        raise ValueError('Current compartment does not have polyubiquitination')

def proteasomal_degradation(gene_info, protein_metabolite, polyub_protein_metabolite, compartment):
    if compartment == 'c':
        if (protein_metabolite.compartment != 'c') or (polyub_protein_metabolite.compartment != 'c'):
            raise ValueError('Compartment mismatch for cytoplasmic proteasomal degradation')
        
        
        deubiquitination = cobra.Reaction(protein_metabolite.id + '_CYTOPLASMIC_DEUBIQUITINATION')
        deubiquitination.subsytem = 'Protein_Expression'
        deubiquitination.add_metabolites({polyub_protein_metabolite: -1, h2o_c: -1, protein_metabolite: 1, polyub_c: 1})
        deubiquitination.gene_reaction_rule = ' and '.join(proteasome_machinery)

        protein_degradation = cobra.Reaction(protein_metabolite.id + '_CYTOPLASMIC_PROTEASOMAL_DEGRADATION')
        protein_degradation.subsytem = 'Protein_Expression'
        rxn = {seq_amino_acid_map_c[aa_code]: aa_counts for aa_code, aa_counts in gene_info.amino_acid_counts.items()}
        rxn[polyub_protein_metabolite], rxn[h2o_c], rxn[polyub_c] = -1, -gene_info.L_protein, 1
        # atp hydrolysis for translocation/unfolding  - known 1 ATP per 2 residues - https://www.nature.com/articles/s41586-018-0736-4
        L_polub_protein = (gene_info.L_protein + (L_monoub*n_ub)) 
        rxn = hydrolyze_atp(rxn, n_atp = L_polub_protein/2, compartment = 'c')


        protein_degradation.add_metabolites(rxn)
        protein_degradation.gene_reaction_rule = ' and '.join(proteasome_machinery)

        protein_degradation_reactions = [deubiquitination, protein_degradation]

        return protein_degradation_reactions
    elif compartment == 'n':
        if (protein_metabolite.compartment != 'n') or (polyub_protein_metabolite.compartment != 'n'):
            raise ValueError('Compartment mismatch for nuclear proteasomal degradation')


        deubiquitination = cobra.Reaction(protein_metabolite.id + '_NUCLEAR_DEUBIQUITINATION')
        deubiquitination.subsytem = 'Protein_Expression'
        deubiquitination.add_metabolites({polyub_protein_metabolite: -1, h2o_n: -1, protein_metabolite: 1, polyub_n: 1})
        deubiquitination.gene_reaction_rule = ' and '.join(proteasome_machinery)

        protein_degradation = cobra.Reaction(protein_metabolite.id + '_NUCLEAR_PROTEASOMAL_DEGRADATION')
        protein_degradation.subsytem = 'Protein_Expression'
        rxn = {seq_amino_acid_map_n[aa_code]: aa_counts for aa_code, aa_counts in gene_info.amino_acid_counts.items()}
        rxn[polyub_protein_metabolite], rxn[h2o_n], rxn[polyub_n] = -1, -gene_info.L_protein, 1
        # atp hydrolysis for translocation/unfolding  - known 1 ATP per 2 residues - https://www.nature.com/articles/s41586-018-0736-4
        L_polub_protein = (gene_info.L_protein + (L_monoub*n_ub)) 
        rxn = hydrolyze_atp(rxn, n_atp = L_polub_protein/2, compartment = 'n')


        protein_degradation.add_metabolites(rxn)
        protein_degradation.gene_reaction_rule = ' and '.join(proteasome_machinery)

        protein_degradation_reactions = [deubiquitination, protein_degradation]

        return protein_degradation_reactions
    else:
        raise ValueError('Current compartment does not have proteasomal degradation')


# # Cytosolic Degradation

# In[193]:


def degrade_cytosolic_protein(gene_info, folded_protein_c):
    polyubiquitinate_folded_protein_c, polyub_protein_c = protein_polyubiquitination(gene_info, 
                                                          protein_metabolite = folded_protein_c, compartment = 'c') 
    cytosolic_proteasomal_degradation_reactions = proteasomal_degradation(gene_info, 
                                                  protein_metabolite = folded_protein_c, 
                                                  polyub_protein_metabolite = polyub_protein_c, compartment = 'c')
    
    
    return [polyubiquitinate_folded_protein_c] + cytosolic_proteasomal_degradation_reactions


# # Nuclear Reactions

# In[213]:


def transport_nuclear_protein(gene_info, folded_protein_c):

    folded_protein_n = folded_protein_c.copy()
    folded_protein_n.id, folded_protein_n.compartment = folded_protein_n.id.replace('[c]', '[n]'), 'n' 

    nuclear_import = cobra.Reaction(gene_info.hgnc_id + '_NUCLEAR_IMPORTtn')
    nuclear_import.subsytem = 'Protein_Expression'
#     nuclear_export = nuclear_import.copy()
#     nuclear_export.id = nuclear_export.id.replace('IMPORT', 'EXPORT')

    import_rxn = {folded_protein_c: -1, folded_protein_n: 1}
#     export_rxn = {folded_protein_c: 1, folded_protein_n: -1}

    if gene_info.protein_mass > nuclear_diffusion_limit:
        # gtp hydrolysis per import
        import_rxn[gtp_n], import_rxn[h2o_n], import_rxn[gdp_n], import_rxn[pi_n], import_rxn[h_n]  = -1, -1, 1, 1, 1
        nuclear_import.add_metabolites(import_rxn)
        nuclear_import.gene_reaction_rule = ' and '.join(importins + RAN)
#         export_rxn[gtp_c], export_rxn[h2o_c], export_rxn[gdp_c], export_rxn[pi_c], export_rxn[h_c]  = -1, -1, 1, 1, 1
#         nuclear_export.add_metabolites(export_rxn)
#         nuclear_export.gene_reaction_rule = ' and '.join(XPO1 + RAN)

    else: # diffusion
        nuclear_import.add_metabolites(import_rxn)
        nuclear_import.lower_bound = -1000
#         nuclear_export.add_metabolites(export_rxn)
        
    return nuclear_import, folded_protein_n

def degrade_nuclear_protein(gene_info, folded_protein_n):
    polyubiquitinate_folded_protein_n, polyub_protein_n = protein_polyubiquitination(gene_info, 
                                                          protein_metabolite = folded_protein_n, compartment = 'n')
    nuclear_proteasomal_degradation_reactions = proteasomal_degradation(gene_info, 
                                                protein_metabolite = folded_protein_n, 
                                                polyub_protein_metabolite = polyub_protein_n, compartment = 'n')
    return [polyubiquitinate_folded_protein_n] + nuclear_proteasomal_degradation_reactions

def get_nuclear_reactions(gene_info, folded_protein_c):
    nuclear_import, folded_protein_n = transport_nuclear_protein(gene_info, folded_protein_c)
    nuclear_degradation_reactions = degrade_nuclear_protein(gene_info, folded_protein_n)
    
    return [nuclear_import] + nuclear_degradation_reactions, folded_protein_n


# # Mitochondrial Reactions
# 

# In[212]:


# i is intermembrane space, but called inner in compartments BIGG
# stick to notation and use inner instead of inter in reaction naming

def transport_mitochondrial_matrix(gene_info, unfolded_protein_c):
    # transport and folding
    if unfolded_protein_c.compartment != 'c':
        raise ValueError('Only cytoplasmic proteins can be transported to mitochondrial matrix')
    
    mitochondrial_matrix_transport = cobra.Reaction(gene_info.hgnc_id + '_MITOCHONDRIAL_MATRIXtn')
    mitochondrial_matrix_transport.subsytem = 'Protein_Expression'
    pre_protein_m = unfolded_protein_c.copy()
    pre_protein_m.id = pre_protein_m.id.replace('[c]', '[m]')
    pre_protein_m.compartment = 'm'
    pre_protein_m.id = pre_protein_m.id.replace('unfolded', 'folded_pre')
    
    rxn = {unfolded_protein_c: -1, pre_protein_m: 1}
    # ATP hydrolysis for transport, assums 1 ATP consumed per 2 residues
    rxn = hydrolyze_atp(rxn, n_atp = gene_info.L_protein*transport_translocation_atp_cost, compartment = 'm')
    
    
    mitochondrial_matrix_transport.add_metabolites(rxn)
    mitochondrial_matrix_transport.gene_protein_rule = ' and '.join(TOM + TIM23_PAM + HSP70_m)
    
    return mitochondrial_matrix_transport, pre_protein_m

def mitochondrial_matrix_protein_processing(gene_info, pre_protein_m):
    # implement this in the future: cleavage of MTS (and degradation of MTS)
    processed_protein_m, aa_counts_processed_m, L_processed_protein_m = pre_protein_m, gene_info.amino_acid_counts.copy(), gene_info.L_protein
    process_mitochondrial_matrix_protein = None
    return process_mitochondrial_matrix_protein, processed_protein_m, aa_counts_processed_m, L_processed_protein_m


def degrade_mitochondrial_protein(gene_info, protein_metabolite, compartment, L_protein, amino_acid_counts):
    rxn = {seq_amino_acid_map_m[aa_code]: aa_counts for aa_code, aa_counts in amino_acid_counts.items()}
    rxn[protein_metabolite], rxn[h2o_m] = -1, -(L_protein-1)
    
    if compartment == 'm':
        mitochondrial_degradation = cobra.Reaction(gene_info.hgnc_id + '_MITOCHONDRIAL_MATRIX_DEGRADATION')
        mitochondrial_degradation.gene_protein_rule = mLON[0]
        
        # ATP hydrolysis by LON: 2 ATP per residue - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518814/
        rxn = hydrolyze_atp(rxn, n_atp = L_protein*2, compartment = 'm')
        

    elif compartment == 'i':
        mitochondrial_degradation = cobra.Reaction(gene_info.hgnc_id + '_INNER_MITOCHONDRIAL_DEGRADATION')
        mitochondrial_degradation.gene_protein_rule = iAAA[0]#' and '.join(mAAA + iAAA)
        
        # ATP hydrolysis by m/i-AAA: 1 ATP per 2 residues -- no source, assumes same as 26S proteasome
        rxn = hydrolyze_atp(rxn, n_atp = L_protein*proteolysis_translocation_atp_cost, compartment = 'i')
  
    mitochondrial_degradation.subsytem = 'Protein_Expression'
    mitochondrial_degradation.add_metabolites(rxn)
    
    return mitochondrial_degradation

def transport_mitochondrial_inter(gene_info, processed_protein_m):
    # upper left Fig 12-29 https://www.ncbi.nlm.nih.gov/books/NBK26828/ 
    # import to matrix then re-export to inter membrane space
    
    if processed_protein_m.compartment != 'm':
        raise ValueError('Only the mechanism of mitochondrial matrix import and re-export to inter membrane is considered')
    
    mitochondrial_inter_transport = cobra.Reaction(gene_info.hgnc_id + '_MITOCHONDRIAL_INNERtn')
    mitochondrial_inter_transport.subsytem = 'Protein_Expression'
    pre_protein_i = processed_protein_m.copy()
    pre_protein_i.id = processed_protein_m.id.replace('[m]', '[i]')
    pre_protein_i.compartment = 'i'
    
    rxn = {processed_protein_m: -1, pre_protein_i: 1}
    
    mitochondrial_inter_transport.add_metabolites(rxn)
    mitochondrial_inter_transport.gene_protein_rule = OXA[0]
    
    return mitochondrial_inter_transport, pre_protein_i

def mitochondrial_inter_protein_processing(gene_info, pre_protein_i):
    # implement this in the future: cleavage of secondary sequence (and degradation)
    processed_protein_i, aa_counts_processed_i, L_processed_protein_i = pre_protein_i, gene_info.amino_acid_counts.copy(), gene_info.L_protein
    process_mitochondrial_matrix_protein = None
    return process_mitochondrial_matrix_protein, processed_protein_i, aa_counts_processed_i, L_processed_protein_i

def get_mitochondrial_reactions(gene_info, unfolded_protein_c, compartments):
    mitochondrial_matrix_transport, pre_protein_m = transport_mitochondrial_matrix(gene_info, unfolded_protein_c)
    process_mitochondrial_matrix_protein, processed_protein_m, aa_counts_processed_m, L_processed_protein_m = mitochondrial_matrix_protein_processing(gene_info, pre_protein_m)
    
    mitochondrial_reactions = [mitochondrial_matrix_transport]
    mitochondrial_protein_metabolites = list()
    if process_mitochondrial_matrix_protein != None:
        mitochondrial_reactions += [process_mitochondrial_matrix_protein]
    
    if 'm' in compartments:
        mitochondrial_matrix_degradation = degrade_mitochondrial_protein(gene_info, protein_metabolite = processed_protein_m, compartment = 'm', L_protein = L_processed_protein_m, amino_acid_counts = aa_counts_processed_m)
        mitochondrial_reactions += [mitochondrial_matrix_degradation]
        mitochondrial_protein_metabolites += [processed_protein_m]
    if 'i' in compartments:
        mitochondrial_inter_transport, pre_protein_i = transport_mitochondrial_inter(gene_info, processed_protein_m)
        process_mitochondrial_inter_protein, processed_protein_i, aa_counts_processed_i, L_processed_protein_i = mitochondrial_inter_protein_processing(gene_info, pre_protein_i)
        if process_mitochondrial_matrix_protein != None:
                mitochondrial_reactions += [process_mitochondrial_inter_protein]        
        mitochondrial_inter_degradation = degrade_mitochondrial_protein(gene_info, protein_metabolite = processed_protein_i, compartment = 'i', L_protein = L_processed_protein_i, amino_acid_counts = aa_counts_processed_i)
        mitochondrial_reactions += [mitochondrial_inter_transport, mitochondrial_inter_degradation]
        mitochondrial_protein_metabolites += [processed_protein_i]

    return mitochondrial_reactions, mitochondrial_protein_metabolites


# # Peroxisomal

# In[211]:


def transport_peroxisome(gene_info, folded_protein_c):
    if folded_protein_c.compartment != 'c':
        raise ValueError('Only cytoplasmic proteins can be transported to mitochondrial matrix')
    
    peroxisomal_transport = cobra.Reaction(gene_info.hgnc_id + '_PEROXISOMEtn')
    peroxisomal_transport.subsytem = 'Protein_Expression'
    folded_protein_x = folded_protein_c.copy()
    folded_protein_x.id = folded_protein_x.id.replace('[c]', '[x]')
    folded_protein_x.compartment = 'x'
    
    rxn = {folded_protein_c: -1, folded_protein_x: 1}
    # ATP hydrolysis for transport--translocation of protein, export of PEX5S receptor
    rxn = hydrolyze_atp(rxn, n_atp = (gene_info.L_protein+L_PEX5)*transport_translocation_atp_cost, 
                        compartment = 'x')
    
    
    peroxisomal_transport.add_metabolites(rxn)
    peroxisomal_transport.gene_protein_rule = ' and '.join(peroxins + AWP1)
    
    return peroxisomal_transport, folded_protein_x

def degrade_peroxisomal_protein(gene_info, folded_protein_x):
    
    
    peroxisomal_degradation = cobra.Reaction(gene_info.hgnc_id + '_PEROXISOMAL_DEGRADATION')
    peroxisomal_degradation.subsytem = 'Protein_Expression'
    peroxisomal_degradation.gene_protein_rule = LONP2[0]

    rxn = {seq_amino_acid_map_x[aa_code]: aa_counts for aa_code, aa_counts in gene_info.amino_acid_counts.items()}
    rxn[folded_protein_x], rxn[h2o_x] = -1, -(gene_info.L_protein-1)
    # ATP hydrolysis by LON: 2 ATP per residue - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518814/
    rxn = hydrolyze_atp(rxn, n_atp = gene_info.L_protein*2, compartment = 'x')
    peroxisomal_degradation.add_metabolites(rxn)
    
    return peroxisomal_degradation

def get_peroxisomal_reactions(gene_info, folded_protein_c):
    peroxisomal_transport, folded_protein_x = transport_peroxisome(gene_info, folded_protein_c)
    peroxisomal_degradation = degrade_peroxisomal_protein(gene_info, folded_protein_x)
    
    return [peroxisomal_transport, peroxisomal_degradation], folded_protein_x


# In[214]:


def get_protein_expression_reactions(gene_info):
    # after transport, expand these to secretory pathways
    protein_expression_reactions, protein_metabolites = list(), list()
    
    # cytoplasmic translation 
    if 'Cytosolic Tranport' in gene_info.final_locations.values() or gene_info.L_protein <= 160: 
        translation_elongation_c, unfolded_protein_c = translate_protein_cytosolic(gene_info)
        protein_expression_reactions.append(translation_elongation_c)
    
    if 'c' in gene_info.final_locations.keys() or 'x' in gene_info.final_locations.keys() or 'n' in gene_info.final_locations.keys():
        protein_folding_cytosolic, folded_protein_c = fold_protein_cytosolic(gene_info, unfolded_protein_c)
        protein_expression_reactions += [protein_folding_cytosolic]

        
        if 'c' in gene_info.final_locations.keys() or 'x' in gene_info.final_locations.keys() or ('n' in gene_info.final_locations.keys() and gene_info.protein_mass <= nuclear_diffusion_limit):
           # cytoplasmic degradation of folded proteins: cytoplasmic proteins, peroxisomal proteins, or nuclear proteins undergoing passive diffusion
            protein_expression_reactions += degrade_cytosolic_protein(gene_info, folded_protein_c)
            
            if 'c' in gene_info.final_locations.keys():
                protein_metabolites += [folded_protein_c]
            
            if 'x' in gene_info.final_locations.keys():
                peroxisomal_reactions, folded_protein_x = get_peroxisomal_reactions(gene_info, folded_protein_c)
                protein_expression_reactions += peroxisomal_reactions
                protein_metabolites += [folded_protein_x]
        
        if 'n' in gene_info.final_locations.keys():
            nuclear_reactions, folded_protein_n = get_nuclear_reactions(gene_info, folded_protein_c)
            protein_expression_reactions += nuclear_reactions
            protein_metabolites += [folded_protein_n]
            
            
    if 'i' in gene_info.final_locations.keys(): # no folding for i, but cytoplasmic degradation
        protein_expression_reactions += degrade_cytosolic_protein(gene_info, unfolded_protein_c)
    
    # mitochondrial transport and degradation ('i' and 'm')
    if ('m' in gene_info.final_locations.keys()) or ('i' in gene_info.final_locations.keys()):
        if ('m' in gene_info.final_locations.keys()) and ('i' in gene_info.final_locations.keys()):
            mitochondrial_reactions, mitochondrial_protein_metabolites = get_mitochondrial_reactions(gene_info, unfolded_protein_c, compartments = ['m','i'])
        elif 'm' in gene_info.final_locations.keys():
            mitochondrial_reactions, mitochondrial_protein_metabolites = get_mitochondrial_reactions(gene_info, unfolded_protein_c, compartments = ['m'])
        elif 'i' in gene_info.final_locations.keys():
            mitochondrial_reactions, mitochondrial_protein_metabolites = get_mitochondrial_reactions(gene_info, unfolded_protein_c, compartments = ['i'])
        protein_expression_reactions += mitochondrial_reactions
        protein_metabolites += mitochondrial_protein_metabolites
    
        
    return protein_expression_reactions, protein_metabolites


# In[9]:


# gene_info.final_locations = {'c': 'Cytosolic Tranport', 'm': 'Cytosolic Tranport', 
#                              'i': 'Cytosolic Tranport'}
# mit_reactions = get_protein_expression_reactions(gene_info)
# mit_mod = cobra.Model('mitochondrial_expression')
# mit_mod.add_reactions(mit_reactions)
# import escher
# builder = escher.Builder(model = mit_mod)


# In[183]:


# sp_dict = {1: True, 0: False, float('nan'): False}
# ptm_cols = ['DSB', 'GPI', 'NG', 'OG']
# ptm_keys = list(allowed_ptms.keys())

# gene1_id = human_model.genes[0].id


# idx  = psim_me[psim_me['HGNC_ID'] == gene1_id].index
# ptms_ = dict(zip(ptm_keys, psim_me.loc[idx, ptm_cols].iloc[0,:].tolist()))
# ptms_ = {k:v for k,v in ptms_.items() if v != 0 and not pd.isna(v)}
# fl = psim_me.loc[idx, 'Location'].tolist()[0]

# pm,m,p = psim_me.loc[idx, 'PREMRNA_SEQ'].tolist()[0], psim_me.loc[idx, 'MRNA_SEQ'].tolist()[0], psim_me.loc[idx, 'PROTEIN_SEQ'].tolist()[0]

# sp = psim_me.loc[idx, 'SP'].tolist()[0]
# if pd.isna(sp):
#     sp = 0
# sp = sp_dict[sp]
# tmd = psim_me.loc[idx,'TMD'].tolist()[0]
# if pd.isna(tmd):
#     tmd = 0
# polyA_length_ = psim_me.loc[idx, 'POLYA_LENGTH'].tolist()[0]
# gene_info = gene_information(metabolic_model = human_model, hgnc_id = gene1_id, 
#                          premrna_seq=pm, mrna_seq=m, protein_seq=p,
#                          ptms = ptms_, tmd = tmd, sp = sp, 
#                         keff = None, polyA_length = polyA_length_, n_introns= None)
# gene_info.get_final_locations(human_model)
# mrna_expression_reactions = bm.mrna_expression(gene_info)


# In[215]:


# gene_info.final_locations = {'n': 'Cytosolic Tranport', 'c': 'Cytosolic Transport', 'i': 'Cytosolic Transport', 
#                             'm': 'Cytosolic Transport', 'x': 'Cytosolic Transport'}
# protein_expression_reactions, protein_metabolites = get_protein_expression_reactions(gene_info)


# In[ ]:





# In[ ]:





# In[ ]:


#ub_reactions

