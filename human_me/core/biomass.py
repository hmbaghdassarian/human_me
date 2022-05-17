#!/usr/bin/env python
# coding: utf-8
from collections import OrderedDict
from typing import Dict, List, Optional, Union
import warnings

import cobra

from human_me.utils import parameters as params
from human_me.utils.parameters import biomass_parameters
from human_me.core.reaction import BiomassReaction
from human_me.core.macromolecules.macromolecule import Macromolecule


class Biomass(cobra.Metabolite):
    """An object of type Biomass inherited from cobra.Metabolite"""

    def __init__(self, id: str = None, compartment: Optional[str] = None):
        super().__init__(id=id, compartment=compartment)


# make the biomass metabolites
biomass_ = Biomass('biomass_total')
biomass_dilution = BiomassReaction('biomass_dilution')
biomass_dilution.add_metabolites({biomass_: -1})
biomass_dilution._lower_bound, biomass_dilution._upper_bound = params.mu, params.mu

biomass_reactions = [biomass_dilution]

# constant
dna_ = Biomass('biomass_DNA')
carb_ = Biomass('biomass_carbohydrate')
lipid_ = Biomass('biomass_lipid')
other_ = Biomass('biomass_other')

# variable
trna_ = Biomass('biomass_tRNA')
rrna_ = Biomass('biomass_rRNA')
mrna_ = Biomass('biomass_mRNA')
premrna_ = Biomass('biomass_premRNA')
other_rna_ = Biomass('biomass_other_RNA')
protein_ = Biomass('biomass_protein')
orphan_protein_ = Biomass('biomass_orphan_protein')

biomass_mapper = {'rrna': rrna_, 'trna': trna_, 'premrna': premrna_, 'mrna': mrna_, 'fragment_rna': other_rna_, 
                'protein': protein_, 'orphan_protein': orphan_protein_}

# biomass formation reactions
constant_biomass_metabolites = [dna_, carb_, lipid_, other_]
biomass_metabolites = constant_biomass_metabolites + [trna_, rrna_, mrna_, premrna_, other_rna_, protein_, orphan_protein_]
for bm in biomass_metabolites:
    reaction_ = BiomassReaction('_'.join(bm.id.split('_')[1:]) + '_biomass_to_biomass')
    reaction_.add_metabolites({bm: -1, biomass_: 1})
    biomass_reactions.append(reaction_)

# # protein biomass with unmodeled protein
upb_reaction = biomass_reactions.pop(len(biomass_reactions) - 1)
# pb_reaction = cobra.Reaction('protein_biomass_to_biomass')
# pb_reaction.add_metabolites({protein_: -1, biomass_: 1})


# The following reactions convert the biomass components which are a constant proportion from the metabolic model formulation to the ME model 
# formulation. Briefly, the coefficients of the precursor reactions must be scaled by their molecular weight, and the product must be equal 
# to the constant proportion of that class of biomass, bounded by growth (flux through reaction = growth rate).

# constant biomass reactions
def create_constant_component_formation(model_metabolites, 
                            # biomass_reactions: List[BiomassReaction] = biomass_reactions, 
                            mass_fraction: Dict[str,float] = biomass_parameters.mass_fraction, 
                            biomass_coefficients: Dict = biomass_parameters.coefficients):
    """Generations formation of biomass components for components that have constant mass (i.e., DNA, lipids, carbohydrates, and other, but not RNA/protein).

    Parameters
    ----------
    model_metabolites : _type_
        _description_
    mass_fraction : Dict[str,float], optional
        the mass fraction of the constant biomass components, by default see utils.parameters.biomass_parameters.mass_fraction
    biomass_coefficients : Dict, optional
        the  biomass component formation metabolite coefficients from the metabolic model that will be scaled by the mass fraction, by default see utils.parameters.biomass_parameters.coefficients

    Returns
    -------
    biomass_reactions : List[BiomassReaction]
        all biomass reactions
    """
    biomass_component_formations = list()
    for biomass_metabolite in constant_biomass_metabolites:
        biomass_type = biomass_metabolite.id.split('_')[1]
        biomass_component_formation = BiomassReaction(biomass_type + '_biomass_formation')

        rxn = {biomass_metabolite: mass_fraction[biomass_type]}
        if biomass_type != 'other':
            for metabolite_id, coef in biomass_coefficients[biomass_type].items():
                rxn[model_metabolites.__dict__[metabolite_id]] = coef*mass_fraction[biomass_type]
        biomass_component_formation.add_metabolites(rxn)
        biomass_component_formation._lower_bound, biomass_component_formation._upper_bound = params.mu, params.mu

        biomass_component_formations.append(biomass_component_formation)

    return biomass_component_formations

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
    biomass_change_ = dict()
    # sorting needed for precision (order of adding masses effects final sum)
    md_map = {m.id: m for m in reaction.metabolites}
    md = OrderedDict({md_map[m_id]: reaction.metabolites[md_map[m_id]] for m_id in sorted(md_map)})

    for metabolite, type in reaction.coupled_metabolites.items():
        md[metabolite] -= metabolite.coupling_coefficient[type]  # coupling not part of mass balance

    for m, count in md.items():
        if m.compartment != 'e': # secreted proteins are not contributing to biomass
            if isinstance(m, Macromolecule):
                if (m.type not in  ['complex', 'protein']) or (m.type == 'protein' and not m.dummy):  # includes ribosomes and non-dummy proteins
                    if m.type in biomass_change_:
                        biomass_change_[m.type] += (count * m.formula_weight / 1000)
                    else:
                        biomass_change_[m.type] = (count * m.formula_weight / 1000)
                elif m.type == 'protein': # dummy proteins
                    if m.dummy_type in biomass_change_:
                        biomass_change_[m.dummy_type] += (count * m.formula_weight / 1000)
                    else:
                        biomass_change_[m.dummy_type] = (count * m.formula_weight / 1000)
                else:  # complexes
                    for type_, mass_ in m.get_complex_biomass().items():
                        if type_ in biomass_change_:
                            biomass_change_[type_] += (count * mass_)
                        else:
                            biomass_change_[type_] = (count * mass_)

    # exclude trna charging/uncharging from change in trna biomass
    # this removes tradeoffs between generating protein biomass and maintaining trna biomass
    if (hasattr(reaction, 'trna_charging') and reaction.trna_charging) or (hasattr(reaction, 'translation') and reaction.translation):
        del biomass_change_['trna']

    # proxy metabolites do not contribute to bimoass
    if 'proxy' in biomass_change_:
        del biomass_change_['proxy']

    biomass_change_ = {biomass_mapper[k]: v for k, v in biomass_change_.items()}

    if inplace:
        reaction.add_metabolites(biomass_change_, combine=False)
    else:
        return biomass_change_

def check_m_biomass(m_model: cobra.Model):
    """Performs sanity checks on biomass objective of metabolic model.
    If using Recon2.2 biomass formulation, this will give warnings. If warnings are given, we recommend using our correct_m_biomass function.

    Parameters
    ----------
    m_model : cobra.Model
        the metabolic model to check. Expect biomass formatting to be consistent with Recon2.2 (especially with regards to reaction and metabolite IDs)
    """
    total_mass = abs(sum([coef for coef in m_model.reactions.biomass_reaction.metabolites.values() if coef < 1]))
    if total_mass != 1:
        warnings.warn('The total mass fraction does not sum to 1')

    # do the substrats for component formation add up to 1g?
    wrn = False
    tol = 10
    for biomass_type in ['lipid', 'DNA', 'carbohydrate']:
        tot_mass = abs(sum([coef*metabolite_.formula_weight for metabolite_, coef in \
            m_model.reactions.get_by_id('biomass_' + biomass_type ).metabolites.items() if not metabolite_.id.startswith('biomass_')]))
        if (tot_mass > 1000 + tol) or (tot_mass < 1000 - tol):
            wrn = True
    if wrn:
        warnings.warn('Some of the biomass component formation reactions are not properly mass balances')

def correct_m_biomass(m_model: cobra.Model):
    """Sets biomass objective to the default one that is used by ME Model. 
    Will correct any mass balance issues in the biomass objective (these mass balance issues are present in Recon2.2).

    Parameters
    ----------
    m_model : cobra.Model
        the metabolic model to correct

    Returns
    -------
    corrected_model : cobra.Model
        the corrected metabolic model
    """
    corrected_model = m_model.copy()


    corrected_model.reactions.EX_biomass_c.lower_bound = 0 # won't effect things, but technically more correct
    for biomass_metabolite_id in ['biomass_DNA_c', 'biomass_lipid_c']:
        biomass_metabolite = corrected_model.metabolites.get_by_id(biomass_metabolite_id)
        biomass_type = biomass_metabolite.id.split('_')[1]

        rxn = {biomass_metabolite: 1}
        for metabolite_id, coef in biomass_parameters.coefficients[biomass_type].items():
            rxn[corrected_model.metabolites.get_by_id(metabolite_id)] = coef

        corrected_model.reactions.get_by_id('biomass_' + biomass_type).add_metabolites(rxn, combine = False)
    return corrected_model


def check_me_biomass(me_model) -> Dict[str, float]:
    """Sanity chack that the ME Model biomass component formation reactions for DNA, lipid, and carbohydrate are formulated correctly.
    This is done by checking mass balance -- that the metabolites sum(stoichiometric coefficient * metabolic weight) = mass fraction

    Parameters
    ----------
    me_model
        the generated ME Model

    Returns
    -------
    Dict[str, float]
        expected mass fractions based on the implemented biomass formation reactions in the ME Model
    """
    expected_mass_fraction = dict()
    for reaction in me_model.reactions:
        if reaction.id.endswith('_biomass_formation') and not reaction.id.startswith('other_'): # filter for biomass component formation reactions except other
            biomass_type = reaction.id.split('_biomass_formation')[0]
            mass_fraction =  abs(sum({coef*metabolite_.formula_weight for metabolite_, coef in reaction.metabolites.items() if not metabolite_.id.startswith('biomass_')})/1e3)
            expected_mass_fraction[biomass_type] = mass_fraction
    return expected_mass_fraction