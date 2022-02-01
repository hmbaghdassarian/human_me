#!/usr/bin/env python
# coding: utf-8
from collections import OrderedDict
from typing import Dict, List, Optional, Union

import cobra

from human_me.utils import parameters as params
from human_me.core.reaction import BiomassReaction
from human_me.core.macromolecules.macromolecule import Macromolecule


class Biomass(cobra.Metabolite):
    """An object of type Biomass inherited from cobra.Metabolite"""

    def __init__(self, id: str = None, compartment: Optional[str] = None):
        super().__init__(id=id, compartment=compartment)


# make the biomass metabolites
biomass_ = Biomass('biomass')
biomass_dilution = BiomassReaction('biomass_dilution')
biomass_dilution.add_metabolites({biomass_: -1})
biomass_dilution._lower_bound, biomass_dilution._upper_bound = params.mu, params.mu

biomass_reactions = [biomass_dilution]

# constant
dna_ = Biomass('biomass_DNA')
carb_ = Biomass('biomass_carbohydrate')
lipid_ = Biomass('biomass_lipid')
# other_ = Biomass('biomass_other')

# variable
protein_ = Biomass('biomass_protein')
unmodeled_protein_ = Biomass('biomass_unmodeled_protein')
trna_ = Biomass('biomass_tRNA')
rrna_ = Biomass('biomass_rRNA')
mrna_ = Biomass('biomass_mRNA')
premrna_ = Biomass('biomass_premRNA')
other_rna_ = Biomass('biomass_other_RNA')

biomass_mapper = {'rrna': rrna_, 'protein': protein_, 'dummy_protein': unmodeled_protein_,
                  'mrna': mrna_, 'trna': trna_, 'fragment_rna': other_rna_,
                  'premrna': premrna_}

# biomass formation reactions

biomass_metabolites = [dna_, carb_, lipid_, trna_, rrna_, mrna_, premrna_, other_rna_,
                       protein_, unmodeled_protein_]  # , other_]
for bm in biomass_metabolites:
    reaction_ = BiomassReaction('_'.join(bm.id.split('_')[1:]) + '_biomass_to_biomass')
    reaction_.add_metabolites({bm: -1, biomass_: 1})
    biomass_reactions.append(reaction_)

# # protein biomass with unmodeled protein
upb_reaction = biomass_reactions.pop(len(biomass_reactions) - 1)
# pb_reaction = cobra.Reaction('protein_biomass_to_biomass')
# pb_reaction.add_metabolites({protein_: -1, biomass_: 1})


# The following reactions convert the biomass components which are a constant proportion from the metabolic model formulation to the ME model formulation. Briefly, the coefficients of the precursor reactions must be scaled by their molecular weight, and the product must be equal to the constant proportion of that class of biomass, bounded by growth (flux through reaction = growth rate).

# constant biomass reactions
# TODO: make *_coef variables customizable by user
def create_biomass_reactions(model_metabolites, biomass_reactions: List[BiomassReaction] = biomass_reactions):
    # DNA------------------------------------------------------
    dna_reaction = BiomassReaction('DNA_biomass_formation')

    # coefs from original RECON2.2
    datp_coef = 0.941642857142857
    dctp_coef = 0.674428571428572
    dgtp_coef = 0.707
    dttp_coef = 0.935071428571429

    # original coefficient from DNA biomass formation reaction*metabolite molecular weight
    rxn = {model_metabolites.datp_n: -datp_coef * model_metabolites.datp_n.formula_weight / 1000,
        model_metabolites.dctp_n: -dctp_coef * model_metabolites.dctp_n.formula_weight / 1000,
        model_metabolites.dgtp_n: -dgtp_coef * model_metabolites.dgtp_n.formula_weight / 1000,
        model_metabolites.dttp_n: -dttp_coef * model_metabolites.dttp_n.formula_weight / 1000,
        dna_: params.DNA_FRAC}
    dna_reaction.add_metabolites(rxn)
    dna_reaction._lower_bound, dna_reaction._upper_bound = params.mu, params.mu

    # CARBOHYDRATE------------------------------------------------------
    g6p_coef = 3.87591549295775
    carbohydrate_reaction = BiomassReaction('carbohydrate_biomass_formation')
    rxn = {model_metabolites.g6p_c: -g6p_coef * model_metabolites.g6p_c.formula_weight / 1000,
        carb_: params.CARB_FRAC}
    carbohydrate_reaction.add_metabolites(rxn)
    carbohydrate_reaction._lower_bound, carbohydrate_reaction._upper_bound = params.mu, params.mu

    # LIPID------------------------------------------------------
    chsterol_coef = 0.210319587628866
    clpn_hs_coef = 0.120185567010309
    pail_hs_coef = 0.240360824742268
    pchol_hs_coef = 1.59237113402062
    pe_hs_coef = 0.570865979381443
    pglyc_hs_coef = 0.0300412371134021
    ps_hs_coef = 0.0600927835051546
    sphmyln_hs_coef = 0.180268041237113

    CLPN_HS_C_MW = 508.21930 / 1000  # ChEBI 28494
    PAIL_HS_C_MW = 387.211 / 1000  # ChEBI 57880
    PCHOL_HS_C_MW = 311.226 / 1000  # ChEBI 64482
    PE_HS_C_MW = 269.146 / 1000  # ChEBI 16038
    PGLYC_HS_C_MW = 299.14860 / 1000  # ChEB 60523
    PS_HS_C_MW = 312.14740 / 1000  # ChEBI 58436
    SPHMYLN_HS_C_MW = 492.630  # ChEBI 62490

    lipid_reaction = BiomassReaction('lipid_biomass_formation')
    rxn = {model_metabolites.chsterol_c: -chsterol_coef * model_metabolites.chsterol_c.formula_weight / 1000,
        model_metabolites.clpn_hs_c: -clpn_hs_coef * CLPN_HS_C_MW,
        model_metabolites.pail_hs_c: -pail_hs_coef * PAIL_HS_C_MW,
        model_metabolites.pchol_hs_c: -pchol_hs_coef * PCHOL_HS_C_MW,
        model_metabolites.pe_hs_c: -pe_hs_coef * PE_HS_C_MW,
        model_metabolites.pglyc_hs_c: -pglyc_hs_coef * PGLYC_HS_C_MW,
        model_metabolites.ps_hs_c: -ps_hs_coef * PS_HS_C_MW,
        model_metabolites.sphmyln_hs_c: -sphmyln_hs_coef * SPHMYLN_HS_C_MW,
        lipid_: params.LIPID_FRAC}
    lipid_reaction.add_metabolites(rxn)
    lipid_reaction._lower_bound, lipid_reaction._upper_bound = params.mu, params.mu

    biomass_reactions += [dna_reaction, carbohydrate_reaction, lipid_reaction]
    return biomass_reactions

def add_biomass_change(reaction: cobra.Reaction, inplace: bool = True) -> Union[None, Dict[str, float]]:
    """Calculate net biomass change in a reaction.

    Parameters
    ----------
    reaction : cobra.Reaction
        reaction uponw which to calculate biomass change
    inplace : bool, optional
        whether to modify the reaction (True) or return a dictionary representation of th new reaction.metabolites, by default True

    Returns
    -------
    None
        if inplace is True, updates the input reaction
    Dict[str, float]
        if inplace is False, the new reaction.metabolites representation is returned
    """

    biomass_change = dict()
    # must order for precision (order of adding masses effects final sum)
    md_ = reaction._metabolites.copy()
    md_map = {m.id: m for m in md_}
    md = OrderedDict({md_map[m_id]: md_[md_map[m_id]] for m_id in sorted(md_map)})

    reaction._map_coupled_metabolites()
    for metabolite, type in reaction.coupled_metabolites.items():
        md[metabolite] -= metabolite.coupling_coefficient[type]  # coupling not part of mass balance

    # extracellular proteins do not contribute to biomass
    #     md = {m:count for m,count in md.items() if m.compartment != 'e'}

    for m, count in md.items():
        if m.compartment != 'e':
            if isinstance(m, Macromolecule):
                if m.type != 'complex':  # includes ribosomes
                    if m.type in biomass_change:
                        biomass_change[m.type] += (count * m.formula_weight / 1000)
                    else:
                        biomass_change[m.type] = (count * m.formula_weight / 1000)
                else:  # complexes
                    for type_, mass_ in m.get_complex_biomass().items():
                        if type_ in biomass_change:
                            biomass_change[type_] += (count * mass_)
                        else:
                            biomass_change[type_] = (count * mass_)

    # exclude trna charging/uncharging from change in trna biomass
    # this removes tradeoffs between generating protein biomass and maintaining trna biomass
    if (hasattr(reaction, 'trna_charging') and reaction.trna_charging) or (
            hasattr(reaction, 'translation') and reaction.translation):
        del biomass_change['trna']

    # proxy metabolites do not contribute to bimoass
    if 'proxy' in biomass_change:
        del biomass_change['proxy']

    biomass_change = {biomass_mapper[k]: v for k, v in biomass_change.items()}

    if inplace:
        reaction.add_metabolites(biomass_change, combine=False)
    else:
        return biomass_change
