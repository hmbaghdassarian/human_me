#!/usr/bin/env python
# coding: utf-8

from typing import Tuple

from human_me.core.macromolecules.macromolecule import Proxy
from human_me.core.macromolecules.protein import Protein
from human_me.core.reaction import ProteinExpressionReaction
from human_me.utils import machinery as mach


def translate_protein_cytosolic(gene_info, mrna_transcript_c, mrna_deg_proxy: Proxy, modified_trna_transcript_c, charged_trna_map, model_metabolites) -> Tuple[ProteinExpressionReaction, Protein]:
    """Generate cytosolic translation reaction

    Parameters
    ----------
    gene_info : GeneInformation
        representation of gene to be expressed
    me_input_model : cobra.Model
        the corrected input metabolic model (as provided in preprocess.correct_inputs.correct_model)
    mrna_transcript_c : mRNA
        the final, cytosolic mRNA transcript
    mrna_deg_proxy : Proxy
        proxy metabolite generated in mRNA degradation reaction for coupling
    modified_transcript_c : macromolecules.RNA.tRNA
        output of create_trna()
    charged_trna_map : Dict[str, macromolecules.RNA.tRNA]
        output of creat_trna()
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin

    Returns
    -------
    translation_elongation: ProteinExpressionReaction
        translation elongation reaction
    unfolded_protein_c: ProteinExpressionReaction
        cytosolic translation product
    """
    # peptide bond formation: https://d1j63owfs0b5j3.cloudfront.net/pop-quiz/answerImage/Amino-Acid-1-popquiz.png
    # tRNA amino acide release: https://rnajournal.cshlp.org/content/14/8/1526/F1.expansion.html

    rxn = dict()
    for aa_code, aa_count in gene_info.amino_acid_counts.items():
        rxn[charged_trna_map[aa_code]] = -aa_count  # tRNA consumption

    rxn[modified_trna_transcript_c] = gene_info.L_protein

    rxn[model_metabolites.h2o_c] = -gene_info.L_protein  # release of peptide from tRNA, addition of -OH to uncharged tRNA
    rxn[model_metabolites.h_c] = gene_info.L_protein  # release of peptide from tRNA, addition of -OH to uncharged tRNA
    rxn[model_metabolites.h2o_c] += gene_info.L_protein - 1  # peptide bond formation (hydrolysis)

    # gtp hydrolysis per aa added
    rxn[model_metabolites.ntp_map_c['G']] = -gene_info.L_protein
    rxn[model_metabolites.h2o_c] -= gene_info.L_protein
    rxn[model_metabolites.ndp_map_c['G']] = gene_info.L_protein
    rxn[model_metabolites.pi_c] = gene_info.L_protein
    rxn[model_metabolites.h_c] += gene_info.L_protein

    unfolded_protein_c = Protein(compartment='c', id_='unfolded', model_metabolites=model_metabolites, gene_info=gene_info)
    rxn[unfolded_protein_c] = 1

    translation_elongation = ProteinExpressionReaction(gene_info.hgnc_id + '_TRANSLATION_ELONGATIONc',
                                                       hgnc_id=gene_info.hgnc_id, translation=True)
    translation_elongation.gene_reaction_rule = ' and '.join(mach.translation_efs + ['ribosome'])  # GPRs

    translation_elongation.add_metabolites(rxn)
    # coupling
    mrna_deg_proxy.couple(value=-gene_info.coupling['mrna_degradation'])
    mrna_transcript_c.couple(type='mrna_formation', value=-gene_info.coupling['mrna_formation'])
    translation_elongation.couple(metabolites=[mrna_deg_proxy, mrna_transcript_c],
                                  types=['mrna_degradation', 'mrna_formation'])

    return translation_elongation, unfolded_protein_c
