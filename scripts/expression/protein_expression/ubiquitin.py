#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.insert(1, '../../../scripts/') # comment out in python script
from utils import machinery as mach
from utils import parameters as params
from utils import metabolites as metab
from utils import functions as func
from core.reaction import Expression_Reaction

from expression.gene_information import gene_information
import expression.build_mrna_expression_reactions as build_mrna
from expression.protein_expression import cytosolic_translation as c_trln

from macromolecules.protein import Protein
from uniform_processes import biomass


# Ubiquitin expression

# In[2]:


def express_ubiquitin(compress_mrna = False):
    # UBC
    ubc_psim = params.psim_me[params.psim_me['HGNC_ID'] == 'HGNC:12468'] # UBC
    ubc_psim['Location'] = 'c'
    ubc_info = gene_information(hgnc_id = ubc_psim['HGNC_ID'].values.tolist()[0], 
                             premrna_seq=ubc_psim['PREMRNA_SEQ'].values.tolist()[0], 
                                mrna_seq=ubc_psim['MRNA_SEQ'].values.tolist()[0], 
                                protein_seq=ubc_psim['PROTEIN_SEQ'].values.tolist()[0],
                                polyA_length = round(ubc_psim['POLYA_LENGTH'].values.tolist()[0]))
    ubc_info.get_final_locations(params.human_model, final_locations=['c'])
    ubc_mrna_expression_reactions, ubc_transcript_c, ubc_deg_proxy = build_mrna.get_mrna_expression_reactions(ubc_info, compress_mrna = compress_mrna)

    # ubiquitin monomer
    single_ubiquitin_sequence = ubc_info.protein_seq[:76]
    monoub_aa_counts = {k: single_ubiquitin_sequence.count(k) for k in params.amino_acids}
    L_monoub = len(single_ubiquitin_sequence)
    n_ub_monomers = ubc_info.protein_seq.count(single_ubiquitin_sequence)
    ub_c = Protein(compartment = 'c', id_ = 'ubiquitin_monomer', amino_acid_counts = monoub_aa_counts)


    ubc_translation_reaction_cytosolic, ubc_c = c_trln.translate_protein_cytosolic(ubc_info, ubc_transcript_c, ubc_deg_proxy)

    ubiquitin_monomerization_ubc = Expression_Reaction(ubc_info.hgnc_id + '_MONOMERIZATIONc',
                                   subsystem = 'Protein_Expression', ubiquitin_biogenesis = True, 
                                   hgnc_id = ubc_info.hgnc_id, synthesis = True, synthesis_type = 'protein')
    rxn = {ubc_c:-1, ub_c: n_ub_monomers, metab.seq_amino_acid_map_c[ubc_info.protein_seq[n_ub_monomers*76:]]: 1, 
          metab.h2o_c: -n_ub_monomers}
    ubiquitin_monomerization_ubc.add_metabolites(rxn)
    ubiquitin_monomerization_ubc.gene_reaction_rule = mach.USP5[0]

    # UBB
    ubb_psim = params.psim_me[params.psim_me['HGNC_ID'] == 'HGNC:12463'] # UBB
    ubb_psim['Location'] = 'c'
    ubb_info = gene_information(hgnc_id = ubb_psim['HGNC_ID'].values.tolist()[0], 
                             premrna_seq=ubb_psim['PREMRNA_SEQ'].values.tolist()[0], 
                                mrna_seq=ubb_psim['MRNA_SEQ'].values.tolist()[0], 
                                protein_seq=ubb_psim['PROTEIN_SEQ'].values.tolist()[0],
                                polyA_length = round(ubb_psim['POLYA_LENGTH'].values.tolist()[0]))
    ubb_info.get_final_locations(params.human_model, final_locations=['c'])
    ubb_mrna_expression_reactions, ubb_transcript_c, ubb_deg_proxy = build_mrna.get_mrna_expression_reactions(ubb_info, compress_mrna = compress_mrna)

    ubb_translation_reaction_cytosolic, ubb_c = c_trln.translate_protein_cytosolic(ubb_info, ubb_transcript_c, ubb_deg_proxy)

    # monomerization from ubb polyub
    n_ub_monomers = ubb_info.protein_seq.count(single_ubiquitin_sequence)
    ubiquitin_monomerization_ubb = Expression_Reaction(ubb_info.hgnc_id + '_MONOMERIZATIONc', 
                                   subsystem = 'Protein_Expression', ubiquitin_biogenesis = True, 
                                   hgnc_id = ubb_info.hgnc_id, synthesis = True, synthesis_type = 'protein')

    rxn = {ubb_c:-1, ub_c: n_ub_monomers, metab. seq_amino_acid_map_c[ubb_info.protein_seq[n_ub_monomers*76:]]: 1, 
          metab.h2o_c: -n_ub_monomers}
    ubiquitin_monomerization_ubb.add_metabolites(rxn)
    ubiquitin_monomerization_ubb.gene_reaction_rule = mach.USP5[0]

    # breakdown of the polyubiquitin cleaved from proteins in ubiquitin-proteasome pathway
    polyub_aa_counts = {aa_code: aa_count*params.n_ub for aa_code, aa_count in monoub_aa_counts.items()}
    polyub_c = Protein(id_ = 'cleaved_polyubiquitin_moiety', amino_acid_counts = polyub_aa_counts,
                        compartment = 'c')
    ubiquitin_monomerization_polyub = Expression_Reaction('POLYUBIQUITIN_MONOMERIZATIONc', 
                                                          subsystem = 'Protein_Expression', 
                                                          ubiquitin_biogenesis = True, hgnc_id = None)

    rxn = {polyub_c:-1, ub_c: params.n_ub, metab.h2o_c: -(params.n_ub-1)}
    ubiquitin_monomerization_polyub.add_metabolites(rxn)
    ubiquitin_monomerization_polyub.gene_reaction_rule = mach.USP5[0]

    # nuclear import of ubiquitin
    nuclear_import_ub_mono = Expression_Reaction('UBIQUITIN_MONOMER_IMPORTtn', subsystem = 'Protein_Expression', 
                                                 ubiquitin_biogenesis = True, hgnc_id = None)
    ub_n = ub_c.change_compartment('n')
    ub_n.id, ub_n.compartment = '_'.join(ub_n.id.split('_')[:-1]) + '_n', 'n'
    nuclear_import_ub_mono.add_metabolites({ub_n: 1, ub_c: -1})
    nuclear_import_ub_mono.lower_bound = -1000

    # nuclear export of polyubiquitin moiety
    nuclear_export_ub_poly = Expression_Reaction('POLYUBIQUITIN_MOIETY_EXPORTtn', 
                                                subsystem = 'Protein_Expression', ubiquitin_biogenesis = True, hgnc_id = None)
    polyub_n = polyub_c.change_compartment('n')
    
    if polyub_n.formula_weight/1000 > params.nuclear_diffusion_limit:
        raise ValueError('Unaccounted for non-passive nuclear export mechanism')
    
    polyub_n.id, polyub_n.compartment = '_'.join(polyub_n.id.split('_')[:-1]) + '_n', 'n'
    nuclear_export_ub_poly.add_metabolites({polyub_n: -1, polyub_c: 1})
    nuclear_export_ub_poly.lower_bound = -1000

    # degradation
    degradation_ub = Expression_Reaction('UBIQUITIN_MONOMER_DEGRADATIONc', 
                                        subsystem = 'Protein_Degradation', ubiquitin_biogenesis = True, hgnc_id = None, 
                                         sink = True, sink_type = 'protein')
    rxn = {metab.seq_amino_acid_map_c[aa_code]: aa_counts for aa_code, aa_counts in monoub_aa_counts.items()}
    rxn[ub_c] = -1
    rxn[metab.h2o_c] =  -(L_monoub-1)
    # atp hydrolysis for translocation/unfolding by 26S - known 1 ATP per 2 residues - https://www.nature.com/articles/s41586-018-0736-4
    rxn = func.hydrolyze_atp(rxn, n_atp = L_monoub/2, compartment = 'c')

    degradation_ub.add_metabolites(rxn)
    degradation_ub.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)

    ub_reactions = ubc_mrna_expression_reactions + [ubc_translation_reaction_cytosolic] 
    ub_reactions += ubb_mrna_expression_reactions + [ubb_translation_reaction_cytosolic]
    ub_reactions += [ubiquitin_monomerization_ubc, ubiquitin_monomerization_ubb, ubiquitin_monomerization_polyub, degradation_ub]
    ub_reactions += [nuclear_import_ub_mono, nuclear_export_ub_poly]
      
    ub_args = {'ub_reactions': ub_reactions, 'ub_c': ub_c, 'ub_n': ub_n, 'polyub_c': polyub_c, 
               'polyub_n': polyub_n, 'monoub_aa_counts': monoub_aa_counts, 'L_monoub': L_monoub, 
              'single_ubiquitin_sequence': single_ubiquitin_sequence}
    
    return ub_args

