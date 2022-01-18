#!/usr/bin/env python
# coding: utf-8

import random
from typing import List, Dict, Optional, Union, Tuple

import cobra
from Bio.Seq import Seq
import numpy as np
import pandas as pd

import human_me.utils.machinery as mach
import human_me.utils.metabolites as metab
import human_me.utils.parameters as params
from human_me.utils.load_environmental_variables import build_files_path

HiddenPrints = params.HiddenPrints


def flatten_list(list1: List[list]) -> list:
    """Create a single list froma  list of lists

    Parameters
    ----------
    list1 : List[list]
        a list of lists

    Returns
    -------
    list
        a single list
    """
    # https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
    return [item for sublist in list1 for item in sublist]


def create_gene_reaction_map(reactions: List[cobra.core.reaction.Reaction]) -> Dict[str, List[cobra.core.reaction.Reaction]]:
    """Maps each gene to all associated reaction.

    Parameters
    ----------
    reactions : List[cobra.core.reaction.Reaction]
        list of cobra reactions

    Returns
    -------
    gene_reaction_map: Dict[str, List[cobra.core.reaction.Reaction]]
        keys are HGNC ID strings, values are a list of reactions from the reactions input in which the HGNC ID
        helps catalyze that reaction
    """
    gene_reaction_map = dict()  # type: Dict[str, List[cobra.core.reaction.Reaction]]
    for r in reactions:
        for g in r.genes:
            if g.id in gene_reaction_map:
                gene_reaction_map[g.id] += [r]
            else:
                gene_reaction_map[g.id] = [r]

    return gene_reaction_map


def convert_gi(gi, non_machinery: List[str]):
    """[summary]

    Parameters
    ----------
    gi : human_me.core.expression.gene_expression.GeneInformation
        [description]
    non_machinery : List[str]
        list of HGNC IDs

    Returns
    -------
    gi : human_me.core.expression.gene_expression.GeneInformation
        updated gene information object
    """
    gi.machinery = True
    gi.all_locations, gi.machinery_locations = gi.nonmachinery_locations.copy(), gi.nonmachinery_locations.copy()
    gi.nonmachinery_locations = dict()
    if gi.hgnc_id in non_machinery:  # create_gene_reaction_map
        gi.non_machinery_locations = gi.format_final_locations(
            final_locations=list(set(non_machinery[gi.hgnc_id]).difference(gi.machinery_locations.keys())),
            sp=True, hgnc_id=gi.hgnc_id)

    return gi


def get_reaction_compartment(reaction: cobra.core.reaction.Reaction, stochastic: bool = False, seed: Optional[int] = None) -> str:
    """Map reactions to a particular compartment according to some rules, informing the
    compartment the enzyme catalyzing the reaction should be transported to.

    Parameters
    ----------
    reaction : cobra.core.reaction.Reaction
        [description]
    stochastic : bool, optional
        In the presence of multiple compartments for a single reaction, whether one should randomly be chosen, by default False
    seed : Optional[int], optional
        A seed for if stochastic is set to True, by default None

    Returns
    -------
    compartment: str
        a singular compartment representing the location of the enzyme catalyzing the reaction
    """
    # only include metabolites with assigned compartments that are not a coupling metabolite
    compartments_ = list(set([m.compartment for m in reaction.metabolites.keys() if m.compartment is not None and not (
        hasattr(m, 'coupling_coefficient') and m.coupling_coefficient is not None)]))
    # sorted to choose the first one in alphabetical order given a tie
    if len(compartments_) > 1:  # for reactions that occur in more than one compartment
        if 'e' in compartments_:
            compartments_ = ['pm']
        else:
            if 'c' in compartments_:  # remove cytoplasmic compartment as a choice in multi-machinery
                compartments_.remove('c')
            if len(set(compartments_)) > 1:
                if not stochastic:
                    seed = 888
                #                     compartments_ = [max(sorted(compartments_), key = compartments_.count)]
                #                 else:
                max_ = max([compartments_.count(i) for i in compartments_])
                compartments_ = list(set([i for i in compartments_ if compartments_.count(i) == max_]))
                random.seed(seed)
                compartments_ = [random.choice(compartments_)]

    #     compartments_ = sorted(set(compartments_))
    if len(compartments_) != 1:
        raise ValueError('Failed to map reaction to a singular compartment')
    if compartments_[0] not in params.compartments:
        raise ValueError('Mapped reaction to a compartment that is not allowed in ME model')
    return compartments_[0]


def hydrolyze_atp(metabolites_to_add: Dict[cobra.core.metabolite.Metabolite, Union[float, int]],
                  n_atp: Union[float, int],
                  compartment: str) -> Dict[cobra.core.metabolite.Metabolite, Union[float, int]]:
    """

    Parameters
    ----------
    metabolites_to_add : Dict[cobra.core.metabolite.Metabolite, Union[float, int]]
        dictionary with metabolite objects or metabolite identifiers as
        keys and reaction coefficients as values
    n_atp : Union[float, int]
        number of atps to hydrolyze
    compartment : str
        location of hydrolysis

    Returns
    -------
    metabolites_to_add : Dict[cobra.core.metabolite.Metabolite, Union[float, int]]
        updated metabolites_to_add dictionary with n_atps hydrolyzed
    """
    n_atp = round(n_atp)

    if metab.atp_compartments[compartment] in metabolites_to_add.keys():
        metabolites_to_add[metab.atp_compartments[compartment]] -= n_atp
    else:
        metabolites_to_add[metab.atp_compartments[compartment]] = -n_atp

    if metab.h2o_compartments[compartment] in metabolites_to_add.keys():
        metabolites_to_add[metab.h2o_compartments[compartment]] -= n_atp
    else:
        metabolites_to_add[metab.h2o_compartments[compartment]] = -n_atp

    if metab.adp_compartments[compartment] in metabolites_to_add.keys():
        metabolites_to_add[metab.adp_compartments[compartment]] += n_atp
    else:
        metabolites_to_add[metab.adp_compartments[compartment]] = n_atp

    if metab.pi_compartments[compartment] in metabolites_to_add.keys():
        metabolites_to_add[metab.pi_compartments[compartment]] += n_atp
    else:
        metabolites_to_add[metab.pi_compartments[compartment]] = n_atp

    if metab.h_compartments[compartment] in metabolites_to_add.keys():
        metabolites_to_add[metab.h_compartments[compartment]] += n_atp
    else:
        metabolites_to_add[metab.h_compartments[compartment]] = n_atp

    return metabolites_to_add


def get_base_counts_and_elements(seq: Union[Seq, str], triphosphate: bool = True) -> Tuple[Dict[str, int], Dict[str, int]]:
    """Count the number of bases and atom elements in an RNA sequence

    Parameters
    ----------
    seq : Union[Seq, str]
        a gene's mRNA sequence
    triphosphate : bool, optional
        whether the mRNA has a 5' triphosphate, by default True

    Returns
    -------
    base_counts: Dict[str, int]
        the number of each RNA base in the sequence
    elements: Dict[str, int]
        the number of atom elements in the sequence
    """

    base_counts = dict()
    for base_letter in metab.seq_element_map:
        base_counts[base_letter] = seq.count(base_letter)

    elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0}
    for base_letter in metab.seq_element_map:
        for element in elements:
            elements[element] += base_counts[base_letter] * metab.seq_element_map[base_letter][element]

            # 3' OH end
    elements['H'] += 1
    elements['O'] += 1

    # 5' end
    if triphosphate:
        elements['P'] += 2
        elements['O'] += 6
    else:
        elements['H'] += 1

    return base_counts, elements


def parse_me_reaction_id(x: str) -> str:
    """Get HGNC ID associated with an expression module reaction

    Parameters
    ----------
    x : str
        cobra.core.reaction.Reaction.id

    Returns
    -------
    x : str
        HGNC ID
    """
    # used to store complex information in build_me
    # since duplicate complexes are checked for, perfect parsing is not necessary, except
    # in the case of protein degradation, since it must match parsing of complex degradation reactions

    if 'HGNC' in x.split('_')[0]:
        pdr = False
        for exception in ['_POLYUBIQUITINATION', '_DEUBIQUITINATION', '_PROTEASOMAL_DEGRADATION']:
            if exception in x:
                pdr = True
        if not pdr:
            return '_'.join(x.split('_')[1:])
        if not x.endswith('pm'):
            return '_'.join(x[x.index('_' + x[-1] + '_') + 2:].split('_')[1:])
        return '_'.join(x[x.index('_' + x[-2:] + '_') + 2:].split('_')[1:])
    return x


def SASA(mw: float) -> float:
    """Estimate the protein solvent-accessible surface area from the molecular weight

    Parameters
    ----------
    mw : float
        protein molecular weight (in kDa)

    Returns
    -------
    float
        approximate protein solvent accesible surface area
    """
    return mw ** 0.75


def average_protein_features(psim_me: pd.DataFrame, context_specific: bool = False) -> pd.DataFrame:
    """Function to get the average protein features from the proteins used in a specific ME model being generated.
    *Note, we filter for metabolic enzymes only, because most orphan reactions come from the metabolic sector.

    Parameters
    ----------
    psim_me : pd.DataFrame
        protein specific information matrix, same as corrected input file (see preprocessing output)
    context_specific : bool, optional
        whether the representative dummy protein is calculated for only genes in the user-provided context specific model from
        the user provided PSIM (True) or for all recon2.2 machinery proteins in the gold-standard PSIM , by default False

    Returns
    -------
    pd.DataFrame
        [description]
    """
    if context_specific:
        psim = psim_me.copy()
        psim = psim[psim.HGNC_ID.isin(mach.metabolic_machinery)]  # filter for metabolic machinery only
    else:
        psim = pd.read_csv(build_files_path + 'recon2_2_only_psim.csv')  # load recon2.2 metabolic machinery gold standard PSIM

    res = pd.DataFrame()
    res['premrna_counts'] = psim.PREMRNA_SEQ.dropna().apply(lambda x: {ntp: x.count(ntp) for ntp in set(x)})
    res['premrna_length'] = psim.PREMRNA_SEQ.dropna().apply(lambda x: len(x))
    res['premrna_prop'] = res.apply(lambda x: {k: v / x.premrna_length for k, v in x.premrna_counts.items()}, axis=1)

    premrna_L = res['premrna_length'].median()
    premrna_avg_prop = {ntp: res['premrna_prop'].apply(lambda x: x[ntp]).median() for ntp in ['A', 'U', 'C', 'G']}

    premrna_seq = ''
    for ntp in ['A', 'U', 'C', 'G']:
        premrna_seq += ntp * int(round(premrna_avg_prop[ntp] * premrna_L))

    res = pd.DataFrame()
    res['mrna_counts'] = psim.MRNA_SEQ.dropna().apply(lambda x: {ntp: x.count(ntp) for ntp in set(x)})
    res['mrna_length'] = psim.MRNA_SEQ.dropna().apply(lambda x: len(x))
    res['mrna_prop'] = res.apply(lambda x: {k: v / x.mrna_length for k, v in x.mrna_counts.items()}, axis=1)

    mrna_L = res['mrna_length'].median()
    mrna_avg_prop = {ntp: res['mrna_prop'].apply(lambda x: x[ntp]).median() for ntp in ['A', 'U', 'C', 'G']}

    mrna_seq = ''
    for ntp in ['A', 'U', 'C', 'G']:
        mrna_seq += ntp * int(round(mrna_avg_prop[ntp] * mrna_L))

    res = pd.DataFrame()
    res['protein_counts'] = psim.PROTEIN_SEQ.dropna().apply(lambda x: {ntp: x.count(ntp) for ntp in set(x)})
    res['protein_length'] = psim.PROTEIN_SEQ.dropna().apply(lambda x: len(x))
    res['protein_prop'] = res.apply(lambda x: {k: v / x.protein_length for k, v in x.protein_counts.items()}, axis=1)

    protein_L = int(round(mrna_L) / 3)

    def get_prop(x, aa):
        if aa in x.keys():
            return x[aa]
        return 0

    protein_avg_prop = {aa: res['protein_prop'].apply(lambda x: get_prop(x, aa)).median() for aa in params.amino_acids}

    protein_seq = ''
    for aa in params.amino_acids:
        protein_seq += aa * int(round(protein_avg_prop[aa] * protein_L))

    dummy_psim = pd.DataFrame(columns=psim.columns)
    dummy_psim.loc[0, :] = float('nan')
    dummy_psim.loc[0, ['HGNC_ID', 'PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ']] = ['HGNC:DUMMY', premrna_seq, mrna_seq,
                                                                                protein_seq]
    dummy_psim.LOCATION = ['[c]']

    # secretory args will be disregarded anyways for now
    median_vals = ['POLYA_LENGTH', 'N_EXONS', 'ALPHA_M', 'ALPHA_P', 'TMD', 'DSB', 'OG', 'NG']
    for col in median_vals:
        dummy_psim[col] = psim[col].median()

    argmax_vals = ['SP', 'GPI']

    # deal with PTR column
    if psim.PTR.dropna().convert_dtypes().dtype is np.dtype('float64'):
        dummy_psim.PTR = psim.PTR.median()
    elif isinstance(psim.PTR.dropna().convert_dtypes().dtype, pd.StringDtype):
        argmax_vals += ['PTR']

    for col in argmax_vals:
        if psim[col].dropna().shape[0] > 0:
            val = psim[col].dropna().value_counts().idxmax()
        else:
            val = float('nan')
        dummy_psim[col] = val

    return dummy_psim


def determine_transport(r: cobra.core.reaction.Reaction) -> List[str]:
    """Checks which metabolites are transported in a cobra reaction

    Parameters
    ----------
    r : cobra.core.reaction.Reaction

    Returns
    -------
    actual_transport_m : List[str]
        a list of metabolites that are actually transported across compartments each element is a string of the metabolite id without the compartment ('_compartment' ending)
    """

    sm_reactants = dict()
    sm_prod = dict()
    for m in r.reactants:
        m_id = m.id.split('_')[:-1][0]
        if m_id not in sm_reactants:
            sm_reactants[m_id] = [m.id.split('_')[-1]]
        else:
            sm_reactants[m_id] += [m.id.split('_')[-1]]
    for m in r.products:
        m_id = m.id.split('_')[:-1][0]
        if m_id not in sm_prod:
            sm_prod[m_id] = [m.id.split('_')[-1]]
        else:
            sm_prod[m_id] += [m.id.split('_')[-1]]
    potential_transport_m = set(sm_prod).intersection(sm_reactants)
    actual_transport_m = [m for m in potential_transport_m if len(set(sm_reactants[m]).intersection(sm_prod[m])) == 0
                          and (sm_reactants[m] + sm_prod[m] != ['e', 'b'])]  # this accounts for LParen_EParen reactions
    return actual_transport_m
