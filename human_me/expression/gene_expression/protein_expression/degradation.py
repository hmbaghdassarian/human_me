#!/usr/bin/env python
# coding: utf-8

# This script provides a complete set of specific degradation reactions for complexes, based on their compartment. For proteins, conditional inclusion of various reactions based on gene features, etc is implemented in the build_protein_express_reactions script.
from typing import List, Optional, Tuple, Union

from human_me.utils import machinery as mach
from human_me.utils import parameters as params
from human_me.utils import functions as func

from human_me.core.macromolecules.protein import Protein
from human_me.core.macromolecules.complex import Complex, RibosomalComplex
from human_me.core.macromolecules.complex import add_complex_metabolites

from human_me.core.reaction import ProteinDegradationReaction, ComplexDegradationReaction

deg_reaction_map = {'protein': ProteinDegradationReaction, 'complex': ComplexDegradationReaction}

DegradedMacromoleculeType = Union[Protein, Complex, RibosomalComplex]
DegradationReactionType = Union[ProteinDegradationReaction, ComplexDegradationReaction]
# # Cytosol and Nucleus


def protein_polyubiquitination(macromolecule: DegradedMacromoleculeType, model_metabolites, **kwargs) -> Tuple[DegradationReactionType, DegradedMacromoleculeType]:
    """Adds 4 ubiquitins to the macromolecule

    Parameters
    ----------
    macromolecule : DegradedMacromoleculeType
        the macromolecule to be polyuibquitinated

    Returns
    -------
    polyubiquitinate_protein : DegradationReactionType
        the polyubiquitination reaction
    polyub_macromolecule : DegradedMacromoleculeType
        the polyubiquitinated macromolecule
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin
    """

    ub_args = kwargs['ub_args']

    r_id = macromolecule.id if macromolecule.type == 'protein' else macromolecule._deg_id
    polyubiquitinate_protein = deg_reaction_map[macromolecule.type](r_id + '_POLYUBIQUITINATION' + macromolecule.compartment,
                                                                           hgnc_id=macromolecule.hgnc_id)

    cmap = {'n': 'n', 'c': 'c', 'pm': 'c'}
    if macromolecule.compartment not in cmap:
        raise ValueError(macromolecule.id + ': Current compartment, ' + macromolecule.compartment + ' does not have polyubiquitination reactions available')

    if macromolecule.type == 'protein':
        polyu_protein_aa_counts = macromolecule._amino_acid_counts.copy()
        for aa_code, aa_counts in ub_args['monoub_aa_counts'].items():
            if aa_code in polyu_protein_aa_counts:
                polyu_protein_aa_counts[aa_code] += aa_counts * params.N_UB
            else:
                polyu_protein_aa_counts[aa_code] = aa_counts * params.N_UB
        polyub_macromolecule = Protein(id_='_'.join(macromolecule.id.split('_')[:-1]) + '_polyub',
                                       model_metabolites=model_metabolites,
                                       compartment=cmap[macromolecule.compartment],
                                       amino_acid_counts=polyu_protein_aa_counts)
        polyub_macromolecule._ptms = macromolecule._ptms
        polyub_macromolecule.hgnc_id = macromolecule.hgnc_id

        # deal with PTMs - mass balance elements appropriately, macromlecule input has the PTMs,
        # and they are removed in proceeding reactions, but polyub_macromolecule needs to have the appropriate ptms
        if macromolecule.compartment == 'pm' and len(polyub_macromolecule._ptms) > 0:
            elements = polyub_macromolecule.elements
            if 'dsb' in polyub_macromolecule._ptms:
                elements['H'] -= 2 * polyub_macromolecule._ptms['dsb']
            if 'og' in macromolecule._ptms.keys():
                number_Oglycans = macromolecule._ptms['og']
                balance_og = {'C': (8 + 6 + 8) * number_Oglycans,  # each 1/3 entry is for each 1/3 reactions in Jahir's model, in case want to separate in the future
                              'H': (13 + 10 + 13) * number_Oglycans,
                              'N': (1 + 0 + 1) * number_Oglycans,
                              'O': (5 + 5 + 5) * number_Oglycans}
                for e, c in balance_og.items():
                    elements[e] += c
            if 'gpi' in macromolecule._ptms.keys():
                for e, c in model_metabolites.balanced_gpi.items():
                    if e in elements:
                        elements[e] += c
                    else:
                        elements[e] = c

            polyub_macromolecule.elements = elements

        if macromolecule.compartment == 'pm':
            polyub_macromolecule = polyub_macromolecule.change_compartment('pm')
    else:
        fused_poly_ub = ub_args['polyub_' + cmap[macromolecule.compartment]].copy()
        fused_poly_ub.id = fused_poly_ub.id.replace('cleaved_', 'fused_')
        # peptide bond formation
        elements = fused_poly_ub.elements.copy()
        elements['H'] -= 2
        elements['O'] -= 1
        fused_poly_ub.elements = elements
        fused_poly_ub.compartment = macromolecule.compartment

        polyub_macromolecule = add_complex_metabolites(cplx=macromolecule,
                                                       met_to_add={fused_poly_ub: 1},
                                                       complex_id=macromolecule.temp_id + '_polyub')

    rxn = {macromolecule: -1, ub_args['ub_' + cmap[macromolecule.compartment]]: -params.N_UB,
           polyub_macromolecule: 1, model_metabolites.h2o_compartments[cmap[macromolecule.compartment]]: params.N_UB}
    # 1 ATP hydrolysis per ubiquitin monomer added (https://link.springer.com/article/10.1007/s10637-020-00894-6)
    rxn = func.hydrolyze_atp(rxn, n_atp=params.N_UB, compartment=cmap[macromolecule.compartment], model_metabolites=model_metabolites)
    polyubiquitinate_protein.add_metabolites(rxn)
    if macromolecule.compartment in ['c', 'n']:
        polyubiquitinate_protein.gene_reaction_rule = ' and '.join(mach.UB_ligases[macromolecule.compartment])
    else:
        polyubiquitinate_protein.gene_reaction_rule = ' and '.join(mach.UB_ligases_c + mach.HSP70_c + mach.HSP90AB1)

    polyubiquitinate_protein._update_tracking(macromolecules=[macromolecule, polyub_macromolecule])
    return polyubiquitinate_protein, polyub_macromolecule


def proteasomal_degradation(macromolecule: DegradedMacromoleculeType, model_metabolites, **kwargs) -> List[DegradationReactionType]:
    """Protoeasomal degradation of macromolecule.   

    Parameters
    ----------
    macromolecule : DegradedMacromoleculeType
        the macromolecule to be degraded
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin

    Returns
    -------
    List[DegradationReactionType]
        the deubiquitination and proteosomal degradation reaction
    """
    polyub_macromolecule = kwargs['polyub_macromolecule']
    ub_args = kwargs['ub_args']

    if macromolecule.compartment != polyub_macromolecule.compartment:
        raise ValueError('Polyubiquitinated and deubiquitinated protein compartment does not match')
    if macromolecule.compartment not in ['c', 'n']:
        raise ValueError(macromolecule.id + ': Only proteins/complexes in nucleus and cytosol are considered for proteasomal degradation')
    r_id = macromolecule.id if macromolecule.type == 'protein' else macromolecule._deg_id

    # ------------------------------deubiquitination------------------------------------
    deubiquitination = deg_reaction_map[macromolecule.type](r_id + '_DEUBIQUITINATION' + macromolecule.compartment,
                                                            hgnc_id=macromolecule.hgnc_id)
    deubiquitination.add_metabolites({polyub_macromolecule: -1,
                                      model_metabolites.h2o_compartments[macromolecule.compartment]: -1,
                                      macromolecule: 1, ub_args['polyub_' + macromolecule.compartment]: 1})

    deubiquitination.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)

    # ------------------------------degradation------------------------------------
    protein_degradation = deg_reaction_map[macromolecule.type](r_id + '_PROTEASOMAL_DEGRADATION' + macromolecule.compartment,
                                                               hgnc_id=macromolecule.hgnc_id, sink=True, sink_type=macromolecule.type)

    protein_length = macromolecule.length if type(macromolecule) != RibosomalComplex else macromolecule.length['protein']
    h2o_length = protein_length

    rcp = False
    if macromolecule.type == 'protein':
        aac = macromolecule._amino_acid_counts.copy()
    else:  # complexes
        if type(macromolecule) != RibosomalComplex:  # non ribosomal complexes
            h2o_length -= sum(macromolecule.decompose_complex().values())  # non-covalent bonds don't require h2o,
        else:  # ribosomal complexes
            rcp = True
            h2o_length -= sum([v for m, v in macromolecule.decompose_complex().items() if m.type == 'protein'])  # non-covalent bonds don't require h2o,
        aac = polyub_macromolecule._amino_acid_counts.copy()
        aac.subtract(ub_args['polyub_' + macromolecule.compartment]._amino_acid_counts)
        h2o_length += 1  # 1 required for fused polyub

    rxn = {model_metabolites.seq_amino_acid_map_compartments[macromolecule.compartment][aa_code]: aa_counts for aa_code, aa_counts in aac.items()}
    rxn[polyub_macromolecule] = -1
    rxn[model_metabolites.h2o_compartments[macromolecule.compartment]] = -h2o_length
    rxn[ub_args['polyub_' + macromolecule.compartment]] = 1
    # atp hydrolysis for translocation/unfolding  - known 1 ATP per 2 residues - https://www.nature.com/articles/s41586-018-0736-4
    rxn = func.hydrolyze_atp(rxn, n_atp=protein_length / 2,
                             compartment=macromolecule.compartment, model_metabolites=model_metabolites)

    # ------------
    # add attribute to indicate whether unique ribosomal complex machinery
    # also adds gene_reaction rules:
    protein_degradation._set_proteasomal_degradation(macromolecule=macromolecule, 
                                                     ribosomal_complex=rcp)
    # ------------
    protein_degradation.add_metabolites(rxn)

    if rcp:
        mdc = macromolecule.decompose_complex()
        rm = [m for m in mdc if m.type != 'protein']
        # hardcoded check
        if sorted([m.id for m in rm]) != ['18s_rrna_c', '28s_rrna_c', '5_8s_rrna_c', '5s_rrna_c'] or not len([m for m in mdc if m.type == 'protein']) > 0:
            err = 'Internal: Only expect rrna of mature ribosome complex to be degraded. Should work with current code'
            err += ' but double check as other RNA molecules are not expected and appropriate machinery is added.'
            raise ValueError(err)

        # 3 spots:
        # here, core.reaction.ComplexDegradationReaction._set_proteasomal_degradation, and build_me_model.generate_complex_reactions

        # Option 1: degrade rRNA with ribosomal degradation - see also
        for rm_ in rm:
            new_rxn = rm_.exonucleolytic_degradation(reaction_name='', update=True)
            rxn = {m: c for m, c in new_rxn.metabolites.items() if not hasattr(m, 'type') or m.type != 'rrna'}
            protein_degradation.add_metabolites(rxn, combine=True)
            # this causes multiple copies of same metabolite to be in reaction
#             for m,c in new_rxn.metabolites.items():
#                 if not hasattr(m, 'type') or m.type != 'rrna':
#                     if m in rxn.keys():
#                         rxn[m] += c
#                     else:
#                         rxn[m] = c
#         # Option 2: degrade proteins with ribosomal degradation, releasing rRNA as intact - see also core.reaction.ComplexDegradationReaction._set_proteasomal_degradation
#         for rm_ in rm:
#             if rm_ in rxn.keys():
#                 rxn[rm_] += mdc[rm_]
#             else:
#                 rxn[rm_] = mdc[rm_]

    # tracking
#     protein_degradation.sink = True
    for r in [deubiquitination, protein_degradation]:
        r._update_tracking(macromolecules=[macromolecule, polyub_macromolecule])

    return [deubiquitination, protein_degradation]


def degrade_cytosolic_nuclear_protein(macromolecule: DegradedMacromoleculeType, model_metabolites, **kwargs) -> List[DegradationReactionType]:
    """Degradation reactions for cytosolic or nuclear proteins.

    Parameters
    ----------
    macromolecule : DegradedMacromoleculeType
        The macromolecule to be degraded
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin

    Returns
    -------
    List[DegradationReactionType]
        All associated degradation reactions of the macromolecule
    """
    ub_args = kwargs['ub_args']
    if macromolecule.type == 'protein':
        polyubiquitinate_macromolecule, polyub_macromolecule = protein_polyubiquitination(macromolecule=macromolecule, model_metabolites=model_metabolites, ub_args = ub_args)
    else:
        polyubiquitinate_macromolecule, polyub_macromolecule = protein_polyubiquitination(macromolecule=macromolecule, model_metabolites=model_metabolites, ub_args=ub_args)
    proteasomal_degradation_rxns = proteasomal_degradation(macromolecule=macromolecule, model_metabolites=model_metabolites,
                                                           polyub_macromolecule = polyub_macromolecule,
                                                           ub_args = ub_args)

    return [polyubiquitinate_macromolecule] + proteasomal_degradation_rxns


# # Mitochondria and Intermembrane Space
def degrade_mitochondrial_protein(macromolecule: DegradedMacromoleculeType, model_metabolites) -> List[DegradationReactionType]:
    """Degradation reactions for proteins localized to mitochondria.

    Parameters
    ----------
    macromolecule : DegradedMacromoleculeType
        The macromolecule to be degraded
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin

    Returns
    -------
    List[DegradationReactionType]
        The degradation reaction
    """
    rxn = {model_metabolites.seq_amino_acid_map_m[aa_code]: aa_counts for aa_code, aa_counts in macromolecule._amino_acid_counts.items()}

    h2o_length = macromolecule.length - 1 if macromolecule.type == 'protein' else macromolecule.length - sum(macromolecule.decompose_complex().values())  # non-covalent bonds, +1,-1 cancel out
    rxn[macromolecule], rxn[model_metabolites.h2o_m] = -1, -h2o_length

    if macromolecule.compartment == 'm':
        mitochondrial_degradation = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_DEGRADATIONm',
                                                                         hgnc_id=macromolecule.hgnc_id, sink=True, sink_type=macromolecule.type)
        mitochondrial_degradation.gene_reaction_rule = mach.mLON[0]

        # ATP hydrolysis by LON: 2 ATP per residue - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518814/
        rxn = func.hydrolyze_atp(rxn, n_atp=macromolecule.length * 2, compartment='m', model_metabolites=model_metabolites)

    elif macromolecule.compartment == 'i':
        mitochondrial_degradation = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_DEGRADATIONi',
                                                                         hgnc_id=macromolecule.hgnc_id, sink=True, sink_type=macromolecule.type)
        mitochondrial_degradation.gene_reaction_rule = mach.iAAA[0]  # ' and '.join(mAAA + iAAA)
        # in the future, may want to add ubqituin-proteasome:

        # ATP hydrolysis by m/i-AAA: 1 ATP per 2 residues -- no source, assumes same as 26S proteasome
        rxn = func.hydrolyze_atp(rxn, n_atp=macromolecule.length * params.PROTEOLYSIS_TRANSLOCATION_ATP_COST, compartment='i', model_metabolites=model_metabolites)

    mitochondrial_degradation.add_metabolites(rxn)

#     mitochondrial_degradation.sink = True
    mitochondrial_degradation._update_tracking(macromolecule)

    return [mitochondrial_degradation]


def degrade_peroxisomal_protein(macromolecule: DegradedMacromoleculeType, model_metabolites) -> List[DegradationReactionType]:
    """Degradation rof proteins localized to peroxisome.

    Parameters
    ----------
    macromolecule : DegradedMacromoleculeType
        The macromolecule to be degraded
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin

    Returns
    -------
    List[DegradationReactionType]
        The degradation reactions
    """

    peroxisomal_degradation = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_DEGRADATIONx',
                                                                   hgnc_id=macromolecule.hgnc_id, sink=True, sink_type=macromolecule.type)
    peroxisomal_degradation.gene_reaction_rule = mach.LONP2[0]

    h2o_length = macromolecule.length - 1 if macromolecule.type == 'protein' else macromolecule.length - sum(macromolecule.decompose_complex().values())  # non-covalent bonds +1,-1 cancel out

    rxn = {model_metabolites.seq_amino_acid_map_x[aa_code]: aa_counts for aa_code, aa_counts in macromolecule._amino_acid_counts.items()}
    rxn[macromolecule], rxn[model_metabolites.h2o_x] = -1, -h2o_length
    # ATP hydrolysis by LON: 2 ATP per residue - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518814/
    rxn = func.hydrolyze_atp(rxn, n_atp=macromolecule.length * 2, compartment='x', model_metabolites=model_metabolites)
    peroxisomal_degradation.add_metabolites(rxn)

#     peroxisomal_degradation.sink = True
    peroxisomal_degradation._update_tracking(macromolecule)

    return [peroxisomal_degradation]


# # Secretory Pathway Degradation
def unfold_secretory_protein(macromolecule: DegradedMacromoleculeType, model_metabolites) -> Tuple[DegradedMacromoleculeType]:
    """Remove PTMs and unfold proteins for lysosomal and secretory compartments. For lysosomal degradation, 
    only remove PTMs, there is no unfolding/misfolding as in ERAD.

    Parameters
    ----------
    macromolecule : DegradedMacromoleculeType
        The macromolecule to be degraded
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin

    Returns
    -------
    Tuple[DegradedMacromoleculeType]
        Unfolded and PTM removed version of macromolecule

    Raises
    ------
    ValueError
        [description]
    """
    if macromolecule.compartment not in ['r', 'l']:
        raise ValueError('Protein metabolite does not have correct compartment')

    unfold_protein = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_UNFOLD' + macromolecule.compartment,
                                                          hgnc_id=macromolecule.hgnc_id)
    unfolded_protein = macromolecule.copy()
    if macromolecule.compartment == 'r':  # not "unfolding/misfolding" for lysosomal degradation
        if 'folded' in unfolded_protein.id:
            unfolded_protein.id = unfolded_protein.id.replace('folded', 'unfolded')
        elif 'HGNC:' in unfolded_protein.id:
            unfolded_protein.id = unfolded_protein.id[:unfolded_protein.id.index('_')] + '_unfolded' + unfolded_protein.id[unfolded_protein.id.index('_'):]
        else:
            unfolded_protein.id = 'unfolded_' + unfolded_protein.id

    rxn = dict()
    rxn[macromolecule] = -1

    unfold_mach = list()
    elements = unfolded_protein.elements.copy()
#     lysosomal_degradation_ptm_condition = 'gpi' in macromolecule._ptms.keys() and len(macromolecule._ptms.keys()) == 1

    # PTM removals HERE
    if 'ng' in macromolecule._ptms.keys():
        raise ValueError('N-glycosylation not yet incorporated')
#     if lysosomal_degradation_ptm_condition:
#         raise ValueError('GPI-anchored proteins with no other ptms should be degraded via lysosomal pathway')
    if 'gpi' in macromolecule._ptms.keys():
        unfolded_protein.id = unfolded_protein.id.replace('_GPI', '')

        for e, c in model_metabolites.balanced_gpi.items():
            elements[e] -= c

        rxn[model_metabolites.gpi_hs_r] = 1
        if macromolecule.compartment == 'r':
            rxn[model_metabolites.hdca_r], rxn[model_metabolites.h_r], rxn[model_metabolites.h2o_r] = -1, -1, 1
        elif macromolecule.compartment == 'l':
            rxn[model_metabolites.hdca_l], rxn[model_metabolites.h_l], rxn[model_metabolites.h2o_l] = -1, -1, 1

    if 'dsb' in macromolecule._ptms.keys():
        unfolded_protein.id = unfolded_protein.id.replace('_DSB', '')
        number_DSB = macromolecule._ptms['dsb']
        elements['H'] += 2 * number_DSB
        # incorporate exchange with reductase in future versions
        if macromolecule.compartment == 'r':
            rxn[model_metabolites.o2_r], rxn[model_metabolites.h2o2_r] = number_DSB, -number_DSB
        elif macromolecule.compartment == 'l':
            rxn[model_metabolites.o2_l], rxn[model_metabolites.h2o2_l] = number_DSB, -number_DSB
        unfold_mach += mach.ERDJ5
    if 'og' in macromolecule._ptms.keys():
        unfolded_protein.id = unfolded_protein.id.replace('_OG', '')
        number_Oglycans = macromolecule._ptms['og']
        balance_og = {'C': (8 + 6 + 8) * number_Oglycans,  # each 1/3 entry is for each 1/3 reactions in Jahir's model, in case want to separate in the future
                      'H': (13 + 10 + 13) * number_Oglycans,
                      'N': (1 + 0 + 1) * number_Oglycans,
                      'O': (5 + 5 + 5) * number_Oglycans}
        for e, c in balance_og.items():
            elements[e] -= c

        if macromolecule.compartment == 'r':
            rxn[model_metabolites.udpacgal_r], rxn[model_metabolites.udpgal_r], rxn[model_metabolites.uacgam_r] = number_Oglycans, number_Oglycans, number_Oglycans
            rxn[model_metabolites.udp_r] = -3 * number_Oglycans
            if model_metabolites.h_r in rxn.keys():
                rxn[model_metabolites.h_r] -= 3 * number_Oglycans
            else:
                rxn[model_metabolites.h_r] = -3 * number_Oglycans
        elif macromolecule.compartment == 'l':
            rxn[model_metabolites.udpacgal_l], rxn[model_metabolites.udpgal_g], rxn[model_metabolites.uacgam_g] = number_Oglycans, number_Oglycans, number_Oglycans
            rxn[model_metabolites.udp_l] = -3 * number_Oglycans
            if model_metabolites.h_l in rxn.keys():
                rxn[model_metabolites.h_l] -= 3 * number_Oglycans
            else:
                rxn[model_metabolites.h_l] = -3 * number_Oglycans

    ###########
    unfolded_protein._ptms = dict()
    unfolded_protein.elements = elements
    rxn[unfolded_protein] = 1

    unfold_protein.add_metabolites(rxn)
    if len(unfold_mach) > 1:
        unfold_protein.gene_reaction_rule = ' and '.join(unfold_mach)
    elif len(unfold_mach) == 1:
        unfold_protein.gene_reaction_rule = unfold_mach[0]

    # PTMs - merged with unfolding reaction for now  <--same structure as adding the PTMs
    #     if 'ng' in macromolecule._ptms.keys():
    #         raise ValueError('N-glycosylation not yet incorporated')
    #     if 'dsb' in macromolecule._ptms.keys():
    #     else:
    #         unmodified_protein_r = unfolded_protein
    #     if 'gpi' in macromolecule._ptms.keys():
    #     else:
    #         unmodified_protein_r = unmodified_protein_r
    #     if 'og' in macromolecule._ptms.keys()
    #     else:
    #         unmodified_protein_r = unmodified_protein_r
    #     if len(macromolecule._ptms) == 0:
    #         unmodified_protein_r = unfolded_protein

    unmodified_protein = unfolded_protein  # for adding PTMs as separate reactions in future, if want to
    unfold_protein._update_tracking(macromolecules=[macromolecule, unfolded_protein, unmodified_protein])

    return unfold_protein, unmodified_protein


def retrograde_er(macromolecule: DegradedMacromoleculeType, model_metabolites, retro_protein_r: Optional[DegradedMacromoleculeType] = None) -> Tuple[DegradationReactionType, DegradedMacromoleculeType]:
    """ER retrograde transport

    Parameters
    ----------
    macromolecule : DegradedMacromoleculeType
        Macromolecule undergoing retrograde transport
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin
    retro_protein_r : DegradedMacromoleculeType, optional
        If the ER version of this protein exists, provide it here, by default None

    Returns
    -------
    retrograde_transport : DegradationReactionType
        the retrograde transport reaction
    retro_protein_r : DegradedMacromoleculeType
        the retrograde transported macromolecule
    """
    if macromolecule.compartment != 'g':
        raise ValueError('ER retrograde transport can only occur for Golgi macromolecules')

    V = macromolecule.formula_weight / 1000 * 1.21 / 1000.0  # Protein Volume in nm^3
    copi_coeff = int(round(143793.19 * params.K_V / V))

    if retro_protein_r is None:
        retro_protein_r = macromolecule.change_compartment('r')

    rxn = {macromolecule: -copi_coeff, retro_protein_r: copi_coeff}
    # gtp hydrolysis
    rxn[model_metabolites.ntp_map_c['G']], rxn[model_metabolites.h2o_c], rxn[model_metabolites.ndp_map_c['G']], rxn[model_metabolites.pi_c], rxn[model_metabolites.h_c] = -127, -127, 127, 127, 127

    retrograde_transport = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_COPI_RETROtr', hgnc_id=macromolecule.hgnc_id)
    retrograde_transport.add_metabolites(rxn)
    retrograde_transport.gene_reaction_rule = ' and '.join(mach.copi_m)
    retrograde_transport._update_tracking([macromolecule, retro_protein_r])

    return retrograde_transport, retro_protein_r


def build_erad_reactions(macromolecule: DegradedMacromoleculeType, model_metabolites, **kwargs):
    """Generate all ERAD reactions

    Parameters
    ----------
    macromolecule : DegradedMacromoleculeType
        macromolecule to undergo ERAD
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin
    """
    if 'unfolded_protein_c' in kwargs:
        if macromolecule.type == 'complex':
            raise ValueError(macromolecule.id + ': Unexpected input argument for ERAD of complexes')
        unfolded_protein_c = kwargs['unfolded_protein_c']
    else:
        unfolded_protein_c = None
    
    if 'ub_args' in kwargs:
        ub_args = kwargs['ub_args']
    elif macromolecule.type == 'complex':
        raise ValueError('"ub_args" must be provided in kwargs for ERAD of complexes')

    macromolecule_r = macromolecule
    if macromolecule.compartment == 'g':
        if macromolecule.type == 'protein':
            err = 'Internal: Unexpected direct ERAD of Golgi protein (hard-coded seperate retrograde transport'
            err += 'reaction); this situation is only for complexes. Should work without this error, but just double check. '
            raise ValueError(err)
        retrograde_transport, macromolecule_r = retrograde_er(macromolecule, model_metabolites=model_metabolites)
    elif macromolecule.compartment == 'r':
        retrograde_transport = None
    else:
        raise ValueError('ERAD is only for Golgi or ER proteins in current model implementation')

    unfold_er_protein, unmodified_protein_r = unfold_secretory_protein(macromolecule_r, model_metabolites=model_metabolites)

    # Retro-translocation
    retrotranslocate_protein = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_RETROTRANSLOCATION', hgnc_id=macromolecule.hgnc_id)
    if unfolded_protein_c is None:  # those that underwent co-translational rather than post-translational
        unfolded_protein_c = unmodified_protein_r.change_compartment('c')
        # replace unfolded bc, if it underwent co-translational translocation
        # the signal peptide degradation step makes it such that it is not the same as any
        # cytosolically translated unfolded proteins (multi-localization), i.e., 'i' destined proteins
        # and their degradation
        unfolded_protein_c.id = unfolded_protein_c.id.replace('unfolded', 'retrotranslocated_unfolded')

    rxn = {unmodified_protein_r: -1, unfolded_protein_c: 1}

    # FUTURE: separate nonglycosylated and glycosylated ERAD
    # if 'ng' in macromolecule._ptms.keys() or 'og' in macromolecule._ptms.keys():

    # glyco based ERAD from Jahir
    rxn = func.hydrolyze_atp(rxn, n_atp=6, compartment='c', model_metabolites=model_metabolites)  # from Jahir retro_TRANSLOC_2
    retrotranslocate_protein.add_metabolites(rxn)
    retrotranslocate_protein.gene_reaction_rule = ' and '.join(mach.retro_mach_glyco)

    erad_reactions = [unfold_er_protein, retrotranslocate_protein]
    if retrograde_transport is not None:
        erad_reactions = [retrograde_transport] + erad_reactions

    # protein portion hard-coded into protein_expression script, complex degradation is added through build_me
    if macromolecule.type == 'protein':
        for r in erad_reactions:
            r._update_tracking([macromolecule, macromolecule_r, unmodified_protein_r, unfolded_protein_c])  # redundanciese dealt with in _consolidate_tracking
        return erad_reactions, unfolded_protein_c
    # this portion is hard-coded in protein_expression portion
    erad_reactions += degrade_cytosolic_nuclear_protein(macromolecule=unfolded_protein_c, model_metabolites=model_metabolites, ub_args = ub_args)
    # if macromolecule.type == 'protein':
    #     erad_reactions += degrade_cytosolic_nuclear_protein(macromolecule=unfolded_protein_c, model_metabolites=model_metabolites, ub_args = ub_args)
    # else: 
    #     erad_reactions += degrade_cytosolic_nuclear_protein(macromolecule=unfolded_protein_c, model_metabolites=model_metabolites, ub_args = ub_args)
    for r in erad_reactions:
        r._update_tracking([macromolecule, macromolecule_r, unmodified_protein_r, unfolded_protein_c])
    return erad_reactions


def build_endocytosis_reactions(macromolecule_pm: DegradedMacromoleculeType, model_metabolites, **kwargs) -> Tuple[List[DegradationReactionType], DegradedMacromoleculeType]:
    """Create the endocytosis reactions

    Parameters
    ----------
    macromolecule_pm : DegradedMacromoleculeType
        macromolecule to undergo endocytosis
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin

    Returns
    -------
    Tuple[List[DegradationReactionType], DegradedMacromoleculeType]
        Associated reactions and output macromolecule after endocytosis
    """
    ub_args = kwargs['ub_args']

    if 'macromolecule_l' in kwargs:
        macromolecule_l = kwargs['macromolecule_l']
        if macromolecule_l is not None and macromolecule_l.type == 'complex':
            raise ValueError(macromolecule_l.id + ': Unexpected provision of lysosomal complex in endocytosis (internal)')
    else:
        macromolecule_l = None

    # polyubiquitination for lysosomal targetting-------------------------------------
    polyubiquitinate_protein, polyub_macromolecule_pm = protein_polyubiquitination(macromolecule=macromolecule_pm, model_metabolites=model_metabolites, 
                                                                                   ub_args = ub_args)
    # endocytosis--------------------------------------------------------------------------
    # combine dequbiquitination with endocytosis https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3987138/
    if macromolecule_l is None:
        macromolecule_l = macromolecule_pm.change_compartment('l')
    else:
        if macromolecule_l.length != macromolecule_pm.length:
            raise ValueError('Endocytosis: lysosomal and plasma membrane protein lengths disagree')

    endocytosis = deg_reaction_map[macromolecule_pm.type](macromolecule_pm._deg_id + '_CLATHRIN_ENDOCYTOSIS',
                                                          hgnc_id=macromolecule_pm.hgnc_id)
    rxn = {polyub_macromolecule_pm: -1, macromolecule_l: 1, model_metabolites.h2o_c: -1, ub_args['polyub_c']: 1}

    # gtp hydrolysis for vesicle scission
    rxn[model_metabolites.ntp_map_c['G']] = -round(macromolecule_pm.length) * params.TRANSPORT_TRANSLOCATION_ATP_COST
    rxn[model_metabolites.h2o_c] -= round(macromolecule_pm.length) * params.TRANSPORT_TRANSLOCATION_ATP_COST
    rxn[model_metabolites.ndp_map_c['G']] = round(macromolecule_pm.length) * params.TRANSPORT_TRANSLOCATION_ATP_COST
    rxn[model_metabolites.pi_c] = round(macromolecule_pm.length) * params.TRANSPORT_TRANSLOCATION_ATP_COST
    rxn[model_metabolites.h_c] = round(macromolecule_pm.length) * params.TRANSPORT_TRANSLOCATION_ATP_COST

    endocytosis.add_metabolites(rxn)
    endocytosis.gene_reaction_rule = ' and '.join(mach.endocytic_machinery)

    for r in [polyubiquitinate_protein, endocytosis]:
        r._update_tracking([macromolecule_pm, macromolecule_l, polyub_macromolecule_pm])
    return [polyubiquitinate_protein, endocytosis], macromolecule_l


def lysosomal_degradation(macromolecule: DegradedMacromoleculeType, model_metabolites) -> List[DegradationReactionType]:
    """Generate reactions and associated macromolecules for lysosomal degradation

    Parameters
    ----------
    macromolecule : DegradedMacromoleculeType
        macromolecule to undergo lysosomal degradation
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin

    Returns
    -------
    List[DegradationReactionType]
        The associated lysosomal degradation reactions
    """
    if macromolecule.compartment != 'l':
        raise ValueError('This reaction only occurs in the lysosome')
    lysosomal_degradation_rxns = list()

    if len(macromolecule._ptms.keys()) > 0:
        # not unfolding lysosomal protein is just removing the PTMs, not an actual unfolding reaction
        unfold_lysosomal_protein, unmodified_macromolecule = unfold_secretory_protein(macromolecule, model_metabolites=model_metabolites)
        lysosomal_degradation_rxns += [unfold_lysosomal_protein]
    else:
        unmodified_macromolecule = macromolecule

    degrade_lysosomal_protein = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_LYSOSOMAL_DEGRADATIONl',
                                                                     hgnc_id=macromolecule.hgnc_id, sink=True, sink_type=macromolecule.type)

    h2o_length = macromolecule.length - 1 if macromolecule.type == 'protein' else macromolecule.length - sum(macromolecule.decompose_complex().values())  # non-covalent bonds, +1,-1 cancel out

    rxn = {model_metabolites.seq_amino_acid_map_l[aa_code]: aa_counts for aa_code, aa_counts in macromolecule._amino_acid_counts.items()}
    rxn[unmodified_macromolecule], rxn[model_metabolites.h2o_l] = -1, -h2o_length
    rxn = func.hydrolyze_atp(rxn, n_atp=macromolecule.length * params.PROTEOLYSIS_TRANSLOCATION_ATP_COST, compartment='l', model_metabolites=model_metabolites)

    degrade_lysosomal_protein.add_metabolites(rxn)
    degrade_lysosomal_protein.gene_reaction_rule = ' and '.join(mach.cathepsins)
#     degrade_lysosomal_protein.sink = True

    lysosomal_degradation_rxns += [degrade_lysosomal_protein]

    for r in lysosomal_degradation_rxns:
        r._update_tracking([macromolecule, unmodified_macromolecule])

    return lysosomal_degradation_rxns


def degrade_lysosomal_pm_protein(macromolecule: DegradationReactionType, model_metabolites, **kwargs) -> List[DegradationReactionType]:
    """Lysoomal degradation of membrane proteins

    Parameters
    ----------
    macromolecule : DegradationReactionType
        macromolecule to undergo degradation
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin

    Returns
    -------
    List[DegradationReactionType]
        associated degradation reactions
    """

    if 'ub_args' in kwargs:
        ub_args = kwargs['ub_args']
    else:
        if macromolecule.compartment == 'pm':
            raise ValueError('ub_args must be provided in **kwargs for degradation of plasma-membrane proteins')
        ub_args = None

    macromolecule_l = macromolecule
    if macromolecule.compartment == 'pm':
        if macromolecule.type == 'protein':
            err = 'Internal: Unexpected direct endocytosis of plasma-membrane protein when should be hardcoded in protein_expression script'
            err += ' this situation is expected only to be used for complexes. Should work without this error,'
            err += ' but just double check.'
            raise ValueError(err)
        endocytosis_reactions, macromolecule_l = build_endocytosis_reactions(macromolecule_pm=macromolecule, model_metabolites=model_metabolites,
                                                                             ub_args = ub_args)
    else:
        endocytosis_reactions = None
    if macromolecule_l.compartment != 'l':
        raise ValueError('Must be a lysosomal protein/complex being degraded')

    deg_reactions = lysosomal_degradation(macromolecule_l, model_metabolites)
    if endocytosis_reactions is not None:
        deg_reactions += endocytosis_reactions

    for r in deg_reactions:
        r._update_tracking([macromolecule, macromolecule_l])  # redundanciese dealt with in _consolidate_tracking
    return deg_reactions


degrade_reaction_map = {'c': degrade_cytosolic_nuclear_protein, 'n': degrade_cytosolic_nuclear_protein,
                        'm': degrade_mitochondrial_protein, 'i': degrade_mitochondrial_protein,
                        'x': degrade_peroxisomal_protein,
                        'r': build_erad_reactions, 'g': build_erad_reactions,
                        'l': degrade_lysosomal_pm_protein, 'pm': degrade_lysosomal_pm_protein}


def degrade(macromolecule: DegradedMacromoleculeType, model_metabolites,  **kwargs) -> List[DegradationReactionType]:
    """Compartment-specific degradation reactions for proteins or protein-protein complexes

    Parameters
    ----------
    macromolecule : DegradedMacromoleculeType
        macromolecule to be degraded
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin

    Returns
    -------
    List[DegradationReactionType]
        degradation reactions

    Raises
    ------
    ValueError
        [description]
    ValueError
        [description]
    """
    if macromolecule.type not in ['protein', 'complex']:
        raise ValueError('Macromolecule to degrade must be protein or complex')
    if macromolecule.type == 'complex' and not macromolecule._deg_initialized:
        macromolecule._initialize_deg_params()

    deg_reactions = degrade_reaction_map[macromolecule.compartment](macromolecule, model_metabolites, **kwargs)
    if isinstance(macromolecule, RibosomalComplex):
        for r in deg_reactions:
            r.ribosome_biogenesis = True

    if macromolecule.compartment not in ['r', 'g'] or macromolecule.type == 'complex':
        dr = deg_reactions
    else:
        dr = deg_reactions[0]

    err = False
    for r in dr:
        r._consolidate_macromolecules()
        if len(r.check_mass_balance()) > 0:
            err = True
    if err:
        raise ValueError(macromolecule.id + ': Degradation reactions are unbalanced')
    return deg_reactions
