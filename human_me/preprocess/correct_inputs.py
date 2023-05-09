#!/usr/bin/env python
# coding: utf-8
import itertools
import json
import logging
import os
import warnings
import urllib.request
from typing import Optional, Tuple, Union, Dict, List

import cobra
import cobra.manipulation.delete as c_del
import pandas as pd
import numpy as np

from human_me.preprocess import parse_complex
from human_me.data.file_paths import build_files_url, build_local_path, input_local_path
from human_me.io import load_metabolic_model, load_psim

logging.basicConfig()
logger = logging.getLogger(cobra.__name__)
logger.setLevel(logging.ERROR)

compartments_me = {'c': 'cytosol', 'l': 'lysosome', 'm': 'mitochondria', 'r': 'endoplasmic reticulum',
                   'e': 'extracellular space', 'x': 'peroxisome/glyoxysome', 'n': 'nucleus', 'g': 'golgi apparatus',
                   'i': 'inner mitochondrial compartment', 'pm': 'plasma membrane', 'b': 'boundary'}


def bool_metabolite(m_id: str, compartment: str, m_model: cobra.Model) -> Tuple[bool, cobra.Metabolite]:
    """Checks if a metabolite is in the model.

    Parameters
    ----------
    m_id : str
        the metabolite id, without the "_compartment" at the end
    compartment : str
        one letter compartment code
    m_model : cobra.Model
        a Cobra Metabolic Model

    Returns
    -------
    bool
        whether the metabolite in that compartment is present in the model
    cobra.core.metabolite
        cobra metabolite in the m_model, if it is present
    """
    try:
        m_id = '_'.join(m_id.split('_')[:-1])
        m = m_model.metabolites.get_by_id(m_id + '_' + compartment)
        return True, m
    except:
        return False, None

def correct_model(model_file: Union[cobra.Model, str] = input_local_path + 'recon2_2.xml', biomass_reaction_id: str = 'biomass_reaction') -> Tuple[cobra.Model]:
    """Makes necessary changes to cobrapy model, largely based on issues encountered with Recon2.2.    
    
    Note that because the ME Model will create a new objective function, the input model's biomass objective function
    is removed. The returned models cm_1 and cm_2 allow the user to compare corrected metabolic models with intact
    biomass objective functions with the output ME Model from human_me.build.build_me_model.build_me. We recommend
    using cm_1 for comparisons.

    Parameters
    ----------
    model_file : Union[cobra.Model, str], optional
        string must be 'full/path/to/input_model.xml', by default input_local_path + 'recon2_2.xml' (full Recon2.2 model)
    biomass_reaction : str, optional
        the id of the input metabolic model's biomass reaction

    Returns
    -------
    cm_1 : cobra.core.Model
        Cobra model with GPRs corrected (according to observed problems with Recon2.2 GPRs) and biomass objective intact
    cm_2 : cobra.core.Model
        Cobra model with GPRs corrected, metabolite transport reactions needed for ME Model incorporated,
        and biomass objective intact
    me_input_model : cobra.core.Model
        Cobra model with GPRs corrected, metabolite transport reactions needed for ME Model incorporated,
        and removed biomass objective. Used as input to building the ME Model.
    """
    with urllib.request.urlopen(build_files_url + "required_metabolic_model_metabolites.json") as url:
        required_metabolites = json.loads(url.read().decode())
    rmd = pd.read_csv(build_files_url + 'required_metabolic_model_metabolites.csv', index_col=0)

    m_model = load_metabolic_model(model_file)

    # check for correct compartments
    different_compartments = list(set(m_model.compartments.keys()).difference(compartments_me.keys()))
    if len(different_compartments) > 0 and different_compartments != ['']:
        err = 'The input metabolic model contains compartments not considered by the ME model. '
        err += 'Please remove the following compartments from your model: ' + ', '.join(different_compartments)
        raise ValueError(err)

    for r in m_model.reactions:
        # incase GPR has redundant complexes (recon2.2 had atleast one instance of this - id = OIVD1m)
        if 'and' in r.gene_reaction_rule and 'or' in r.gene_reaction_rule: 
            machinery_final = parse_complex.eval_complex(r.gene_reaction_rule)

            idx = list(itertools.combinations(range(len(machinery_final)), 2))

            rm = []
            for i in idx:
                if machinery_final[i[0]] == machinery_final[i[1]]:
                    rm.append(i[1])

            if len(rm) > 0:
                msg = r.id + ' contains redundant complexes according to GPR, editing GPR'
                warnings.warn(msg)
                machinery_final = [machinery_final[i] for i in range(len(machinery_final)) if i not in rm]
                new_gpr = ''
                for complex_ in machinery_final:
                    if type(complex_) == list:
                        new_gpr += '(' + ' and '.join(complex_) + ')' + ' or '
                    else:
                        new_gpr += complex_ + ' or '
                new_gpr = new_gpr[:-4]

                m_model.reactions.get_by_id(r.id).gene_reaction_rule = new_gpr
    #ensure reversibility is respected and exchanges do not have a GPR
    # exchange reactions are just feeding into the model from boundary, and boundary to ECM, so should not be enzyme catalyzed; otherwise, coupling will be wrong
    for r in m_model.reactions:
        if r.reversibility: 
            lb, ub = r.bounds
            if not lb < 0:
                warnings.warn(r.id + ' reaction is supposed to be reversible, but bounds do not agree, changing lower bound to -1000')
                lb = -1000
                r._bounds = (lb, ub)
            if not ub > 0:
                warnings.warn(r.id + ' reaction is supposed to be reversible, but bounds do not agree, changing upper bound to 1000')
                r._bounds = (lb, 1000)
        elif not r.reversibility:
            lb, ub = r.bounds
            if lb < 0 and ub <= 0:
                if r.gene_reaction_rule != '':
                    msg = 'One-directional reactions going in the reverse reaction should only be for exchanges, which should not be enzyme catalyzed'
                    msg = 'Please change your reaction, ' + r.id + ', from a "product <-- substrate" to "substrate --> product" format'
                    raise ValueError(msg)
        if r.id.startswith('EX_'):
            if r.gene_reaction_rule != '':
                warnings.warn('Exchange reaction ' + r.id + ' has a catalyzing enzyme, but exchange reactions should be diffusion processes, removing')
                r.gene_reaction_rule = ''


    # remove psuedogene w/ no sequences
    g = [g for g in m_model.genes if g.id == 'HGNC:4686']
    if len(g) > 0:
        c_del.remove_genes(cobra_model=m_model, gene_list=g, remove_reactions=True)

    print('Check for the recon2.2 HGNC:HGNC error')
    # correct genes with HGNC:HGNC:, recon2.2 has this
    genes_to_format = [g.id for g in m_model.genes if 'HGNC:HGNC:' in g.id]

    if len(genes_to_format) > 0:
        warnings.warn('Your metabolic model contains genes with HGNC:HGNC:####, changing to HGNC:####')
        correct_format = ['HGNC:' + g.split('HGNC:')[-1] for g in genes_to_format]
        formatting_dict = dict(zip(genes_to_format, correct_format))

        for g in genes_to_format:
            reactions = list(m_model.genes.get_by_id(g).reactions)
            for r in reactions:
                r.gene_reaction_rule = r.gene_reaction_rule.replace(g, formatting_dict[g])

    print('Remove genes not participating in reactions')
    genes_to_remove = [g for g in m_model.genes if len(g.reactions) == 0]
    if len(genes_to_remove) > 0:
        if sorted(genes_to_format) != sorted([g.id for g in genes_to_remove]):
            warnings.warn('Your metabolic model contains genes not involved in reactions, removing these genes')

        c_del.remove_genes(cobra_model=m_model, gene_list=genes_to_remove, remove_reactions=True)
    # check for minimum required metabolites
    all_required_metabolites = [item for sublist in list(required_metabolites.values()) for item in sublist]
    all_metabolites = [m.id for m in m_model.metabolites]
    missing_metabolites = sorted(set(all_required_metabolites).difference(all_metabolites))

    cm_1 = m_model.copy()  # metabolic model with incorrect GPRs corrected

    comp_ = ['c', 'n', 'r', 'g', 'm', 'l', 'x', 'i', 'e', 'b', 'pm']

    for m_id in missing_metabolites:
        nc = m_id.split('_')[-1]

        counter_ = 0
        found = False
        while not found and counter_ < len(comp_):
            cp_ = comp_[counter_]
            found, m_t = bool_metabolite(m_id, cp_, m_model)
            counter_ += 1

        if found:  # metabolite was found in another compartment of the model, add a transport
            new_metab = m_t.copy()
            new_metab.compartment = nc
            new_metab.id = m_id

            transport_rxn = cobra.Reaction('_'.join(m_t.id.split('_')[:-1]) + 't' + nc)
            transport_rxn.add_metabolites({m_t: -1, new_metab: 1})
            transport_rxn.lower_bound = -1000

            m_model.add_reactions([transport_rxn])
        else:
            core_id = '_'.join(m_id.split('_')[:-1])

            wrn = core_id + ' does not exist in model. Adding to compartment ' + nc + ' via sink '
            wrn += 'This allows ' + core_id + ' to be in the model at no cost.'
            warnings.warn(wrn)

            new_metab = cobra.Metabolite(m_id)
            info = rmd.loc[core_id, ]

            new_metab.name = info['name']
            new_metab.compartment = nc
            new_metab.charge = int(info['charge'])
            new_metab.elements = eval(info['elements'])
            new_metab.formula = info['formula']

            if new_metab.compartment != 'e':
                m_model.add_boundary(metabolite=new_metab, type='sink')
            else:
                m_model.add_boundary(metabolite=new_metab, type='exchange')
 
    # add reactions required for me modeling that are not present in recon2 but should be
    reactions_to_add = []

    # 1) nuclear proton transport (diffusion)
    nt_present = [r for r in m_model.metabolites.get_by_id('h_n').reactions if len(r.compartments) > 1] # check whether transport exists
    if len(nt_present) == 0:
        nuclear_proton_transport = cobra.Reaction(id = 'Htn', name = 'H transporter, nucleus', lower_bound = -1000, upper_bound = 1000)
        nuclear_proton_transport.add_metabolites({m_model.metabolites.get_by_id('h_c'): -1, 
                                                m_model.metabolites.get_by_id('h_n'): 1})
        nuclear_proton_transport.bounds = (-1000, 1000)
        reactions_to_add.append(nuclear_proton_transport)
    m_model.add_reactions(reactions_to_add)

    cm_2 = m_model.copy()  # metabolic model with incorrect GPRs corrected and ME-Model required metabolite transport added
    # remove biomass
    metabolites_1 = [m.id for m in m_model.metabolites]
    try:
        biomass_reaction = m_model.reactions.get_by_id(biomass_reaction_id)
        try:
            biomass_m = list(biomass_reaction.metabolites.keys())

            m_model.remove_reactions([biomass_reaction_id], remove_orphans=True)

            for bm in biomass_m:
                reactions = list(bm.reactions)
                if len(reactions) != 1:
                    raise ValueError('Unexpected  reactions associated with biomass constituent')
                bmr = reactions[0]
                m_model.remove_reactions([bmr], remove_orphans=True)
        except:
            raise ValueError('Biomass reactions formatted differently than expected, please emulate RECON2.2')
    except:
        wrn_ = 'No biomass reaction identified in input model. Assuming that biomass reactions and metabolites '
        wrn_ = " are not present. If present, please make sure the biomass reaction id is 'biomass_reaction'"
        warnings.warn(wrn_)

    if len([m.id for m in m_model.metabolites if 'biomass' in m.id]) != 0:
        err_ = 'Extraneous biomass metabolites not associated with the biomass reaction are present,'
        err_ += ' please remove from input cobrapy model'
        raise ValueError(err_)

    rm = sorted(set(metabolites_1).difference([m.id for m in m_model.metabolites]))
    if len([i for i in rm if 'biomass' not in i]) != 0:
        warnings.warn('Non biomass metabolites removed as orphan metabolites from biomass reactions')
    me_input_model = m_model.copy()
    del m_model

    return cm_1, cm_2, me_input_model


def check_non_machinery(non_machinery: Optional[Dict[str, List[str]]] = None) -> Dict[str, List[str]]:
    """Runs checks on non-machinery input list. Non-machinery are categorized as any proteins to express in the ME-Model that are not catalyzed in the reaction.

    Parameters
    ----------
    non_machinery : Dict[str, List[str]], optional
        keys are HGNC IDs, values are a list of strings, each element of which represents a compartment within the 
        metabolic model for the gene to be expressed, by default None
        We define machinery as proteins that are utilized in the reaction GPRs. In its current format, it is not possible for a protein to both be machinery and non-machinery

    Returns
    -------
    non_machinery : Dict[str, List[str]]
        removes any non-machinery that would not have worked with model
    """
    if non_machinery is None:
        return dict()

    for hgnc_id, compartments in non_machinery.items():
        if not hgnc_id.startswith('HGNC:'):
            raise ValueError('All non-machinery must be in HGNC ID format')
        non_machinery[hgnc_id] = list(set(compartments).intersection(compartments_me))

    return non_machinery


def get_status(psim_me: pd.DataFrame) -> pd.DataFrame:
    """Checks sequence columns for validity."""
    psim = psim_me.copy()
    psim['Status'] = 1

    # annotate invalid entries
    for col in ['PROTEIN_SEQ', 'MRNA_SEQ', 'PREMRNA_SEQ', 'HGNC_ID']:
        psim.loc[psim[psim[col].isna()].index, 'Status'] = 0

    premrna_l = psim.PREMRNA_SEQ.apply(lambda x: len(x) if type(x) == str else float('nan'))
    mrna_l = psim.MRNA_SEQ.apply(lambda x: len(x) if type(x) == str else float('nan'))
    protein_l = psim.PROTEIN_SEQ.apply(lambda x: len(x) if type(x) == str else float('nan'))
    psim.loc[psim[(premrna_l < mrna_l) | (mrna_l < (3 * protein_l))].index, 'Status'] = 0

    temp = psim[(premrna_l == mrna_l)]
    psim.loc[temp[temp.MRNA_SEQ != temp.PREMRNA_SEQ].index, 'Status'] = 0

    for col in ['MRNA_SEQ', 'PREMRNA_SEQ']:
        idx = psim[psim[col].apply(
            lambda x: len(set(list(x)).difference(['A', 'U', 'G', 'C', 'N'])) > 0 if type(x) == str else True)].index
        psim.loc[idx, 'Status'] = 0

    amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    idx = psim[psim['PROTEIN_SEQ'].apply(
        lambda x: len(set(list(x)).difference(amino_acids + ['X', 'U'])) > 0 if type(x) == str else True)].index
    psim.loc[idx, 'Status'] = 0

    return psim


def correct_psim(me_input_model: Union[cobra.Model, str],
                psim_df: Union[pd.DataFrame, str] = build_local_path + 'psim_gold.h5',
                 fill_na: str = 'default',
                 non_machinery: Optional[Dict[str, List[str]]] = None):
    """Makes sure PSIM has all necessary correct information to build ME Model.

    Note, the default psim_file, build/psim_gold.h5, is a PSIM generated from MANE/RefSeq Select isoforms.
    We refer to this as the gold standard "psim_gold".

    Parameters
    ----------
    me_input_model : Union[cobra.Model, str]
        Cobra model with GPRs corrected, metabolite transport reactions needed for ME Model incorporated,
        and removed biomass objective. Used as input to building the ME Model. Output of correct_model function.
        Can be 'full/path/to/corrected_model.xml'
    psim_df : Union[pd.DataFrame, str], optional
        See PSIM_README for details on format of psim, by default build_local_path +'psim_gold.h5' (gold-standard PSIM)
    fill_na : str, optional
        options ['default', 'select'], by default 'default'
        if default: will fill incomplete values with default values (see PSIM_README for details)
        if select: will fill incomplete values with the gold standard PSIM when available, otherwise with default

        Note: this will not deal with user-provided incorrect values for optional columns, those will revert to default in the GeneInformation class

        Exceptions:
            for required columns, if incorrect in input psim, will fill with  the gold standard PSIM. required columns include *_SEQ and LOCATION for non-machinery
            if PTR is present in input PSIM and a tissue is specified most often in that column (ignoring NaN), all nan values in that column will default to that tissue
            if non-machinery do not specify an appropriate LOCATION, will fill with  the gold standard PSIM
    non_machinery : dict, optional
        keys are HGNC IDs, values are a list of strings, each element of which represents a compartment
        within the metabolic model for the gene to be expressed, by default None
        We define machinery as proteins that are utilized in the reaction GPRs. In its current format, it is not possible for a protein to both be machinery and non-machinery

    Returns
    -------
    psim_me : pd.DataFrame
        corrected version of input PSIM
    non_machinery : Dict[str, List[str]]
        removes any non-machinery that would not have worked with model
    revised_genes : Dict[str, List[str]]
        'added' is a list of genes missing (relative to me_input_model and non-machinery genes) in PSIM that were added
        'sequences' key is a superset of added, includes all genes with adjusted sequences
        'non-machinery locations' is for genes in non-machinery that did not have an appropriately specified location
    """
    me_input_model = load_metabolic_model(me_input_model)
    # run basic non-machinery check
    non_machinery = check_non_machinery(non_machinery=non_machinery)

    expression_machinery = sorted(pd.read_csv(build_files_url + 'machinery/expression_machinery.txt', header = None)[0].tolist())
    metabolic_machinery = [g.id for g in me_input_model.genes]

    # define the required/optional columns-----------------------------------------------------------------
    user_provided = ['HGNC_ID']  # must be fully provided by user
    essential_cols = ['PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ']  # required but can be filled by build/psim_me
    optional_cols = ['POLYA_LENGTH', 'N_EXONS', 'TMD', 'SP', 'DSB', 'GPI', 'OG', 'NG',
                     'ALPHA_M', 'ALPHA_P', 'PTR']
    nm_cols = ['LOCATION']
    all_columns = user_provided + essential_cols + optional_cols + nm_cols

    # load the MANE/RefSEQ Select PSIM---------------------------------------------------------------------
    psim_gold = pd.read_hdf(build_local_path + 'psim_gold.h5')
    psim_gold = psim_gold[psim_gold.Status != 0]  # drop genes that won't work with model
    psim_gold = psim_gold[all_columns]

    psim_me = load_psim(psim_df)

    psim_me.reset_index(inplace=True, drop=True)
    # check that required cols are all present and appropriately formatted------------------------------------
    err = False
    for col in user_provided:
        if col not in psim_me.columns.tolist():
            err = True
            break
        if psim_me[col].isna().any():
            err = True
            break
        if psim_me[col].unique().shape[0] != psim_me.shape[0]:
            err = True
            break
    if err:
        raise ValueError(
            'The following columns must be present, unique, and completely filled in: ' + ', '.join(user_provided))

    for col in essential_cols:
        if col not in psim_me.columns:
            psim_me[col] = float('nan')
    missing_cols = sorted(set(all_columns).difference(psim_me.columns))

    # check genes are present--------------------------------------------------------------------------
    proteins = list(set(metabolic_machinery + list(non_machinery)))
    psim_me_genes = psim_me.HGNC_ID.tolist()
    psim_gold_genes = psim_gold.HGNC_ID.tolist()

    missing_genes = sorted(set(proteins).difference(psim_me_genes + psim_gold_genes))
    if len(missing_genes) > 0:  # check all but expression machinery
        raise ValueError(
            'The following specified genes are not present in the provided psim or build/psim_gold.h5: ' + ', '.join(
                missing_genes))
    # add expression machinery to lists
    proteins = sorted(set(proteins + expression_machinery))
    revised_genes = {'missing': missing_genes}
    missing_genes = sorted(set(proteins).difference(psim_me_genes))

    # initialize missing columns and genes----------------------------------------------------------------
    for col in missing_cols:
        psim_me[col] = float('nan')
    psim_me = psim_me[all_columns]  # filter out excess columns and order

    # make sure RPL40 and RPS27A are ubiquitin fusions
    rbp_ubs = ['HGNC:10417', 'HGNC:12458']
    ub_psim = psim_gold[psim_gold.HGNC_ID.isin(rbp_ubs)]
    SINGLE_UB_SEQ = ub_psim['PROTEIN_SEQ'].iloc[0, ][:76]

    for ru in rbp_ubs:
        if ru not in missing_genes:
            rui = psim_me[psim_me.HGNC_ID == ru]
            if rui.PROTEIN_SEQ.iloc[0, ][:len(SINGLE_UB_SEQ)] != SINGLE_UB_SEQ:
                msg = 'Provided sequence for ribosomal protein ' + ru + ' must contain specific ubiquitin fusion.'
                msg += 'Replacing with correct sequences.'
                warnings.warn(msg)
                ui = ub_psim[ub_psim.HGNC_ID == ru]
                psim_me.loc[rui.index, ['PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ']] = \
                    ub_psim.loc[ui.index, ['PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ']]

    # initialize missing genes
    for gene in missing_genes:
        psim_me.loc[psim_me.shape[0], :] = [gene] + ([float('nan')] * (psim_me.shape[1] - 1))

    # fill in required columns independent of fill_na arg--------------------------------------------------
    psim_me = psim_me[psim_me.HGNC_ID.isin(proteins)]
    psim_gold = psim_gold[psim_gold.HGNC_ID.isin(proteins)]

    if psim_me.HGNC_ID.unique().shape[0] != psim_me.shape[0]:
        raise ValueError('Duplicate HGNC IDs')
    psim = get_status(psim_me)
    fix = psim[psim.Status == 0].HGNC_ID.tolist()

    psim_gold.index = psim_gold.HGNC_ID
    if len(set(fix).difference(psim_gold_genes)) == 0:
        psim_me.index = psim_me.HGNC_ID
        psim_me.loc[fix, essential_cols] = psim_gold.loc[fix, essential_cols]
        psim_me.reset_index(drop=True, inplace=True)
        revised_genes['sequences'] = fix
    else:
        mssg = ' The following required genes have incorrect sequences in the PSIM: ' + ', '.join(
            list(set(fix).difference(psim_gold_genes)))
        mssg += ' mRNA sequence lengths must be <= premrna sequence lengths and >= 3*protein sequence length'
        mssg += ' Furthermore, they must have only the allowed letters'
        raise ValueError(mssg)

    # ptr specific changes-----------------------------------------------------------------------------
    fill_ptr = True
    if 'PTR' not in missing_cols:
        max_val = psim_me['PTR'].dropna().value_counts().index.tolist()
        if len(max_val) > 0 and type(max_val[0]) == str:
            ptr = pd.read_csv(build_files_url + 'PTR_Gagneur_processed.tsv', sep='\t', index_col=0)
            ptr.drop(columns=['ENSG_ID'], inplace=True)
            ptr.columns = pd.Series(ptr.columns).apply(lambda x: x.split('_')[0] if '_PTR' in x else x).tolist()
            if max_val[0] in ptr.columns.tolist():
                psim_me['PTR'] = max_val[0]
                fill_ptr = False

    # fill_na values-----------------------------------------------------------------------------------
    psim_me.index = psim_me.HGNC_ID
    psim_gold.index = psim_gold.HGNC_ID

    if fill_na == 'select':
        for col in optional_cols:
            if col != 'PTR' or fill_ptr:
                temp = psim_me[psim_me[col].isna()]
                temp = psim_gold[psim_gold.HGNC_ID.isin(temp.HGNC_ID)]
                temp = temp[temp[col].notna()]
                psim_me.loc[temp.index, col] = temp.loc[:, col]
    elif fill_na != 'default':
        raise ValueError('Invalid value passed for fill_na')

    # deal with non_machinery location---------------------------------------------------------------
    if len(non_machinery) > 0 and 'LOCATION' not in psim_me.columns.tolist():
        psim_me['LOCATION'] = float('nan')
    # compartments = ['[c]', '[e]', '[l]', '[m]', '[r]', '[n]', '[g]', '[x]', '[b]', '[i]', '[pm]']
    #     temp = psim_me[psim_me.HGNC_ID.isin(non_machinery) & ((psim_me.LOCATION.isna()) | psim_me.LOCATION.apply(lambda x: x not in compartments))]
    #     err = False
    #     if len(set(temp.HGNC_ID).difference(psim_gold.HGNC_ID)) > 0:
    #         raise ValueError('The final location for the following non-machinery must be specified: ' + ', '.join(list(set(temp.HGNC_ID).difference(psim_gold.HGNC_ID))))
    #     temp = psim_gold.loc[temp.index, :]
    #     if temp.LOCATION.dropna().shape[0] < temp.shape[0]:
    #         raise ValueError('The final location for the following non-machinery must be specified: ' + ', '.join(temp[temp.LOCATION.isna()].HGNC_ID.tolist()))

    #     revised_genes['non-machinery locations'] = []
    #     if temp.shape[0] > 0:
    #         psim_me.loc[temp.index, 'LOCATION'] = temp.LOCATION.tolist()
    #         revised_genes['non-machinery locations'] = temp.HGNC_ID.tolist()

    # filter out unecessary gene entries
    psim_me = psim_me.loc[sorted(set(expression_machinery + metabolic_machinery + list(non_machinery))), :]
    psim_me.reset_index(inplace=True, drop=True)

    del psim_gold

    return psim_me, non_machinery, revised_genes

def format_exchanges(m_model):
    """Formats typical exchange reactions (transport across ECM) as those in recon2.2 (LPAREN_e_RPAREN_, transport across boundary)
    for the input metabolic model. 

    Parameters
    ----------
    m_model : cobra Model
        model to format
    """
    rm = []
    new_exchanges = []
    for er_e in m_model.exchanges:
        if len(er_e.metabolites) != 1:
            raise ValueError('Unexpected metabolites')
        em_e = list(er_e.metabolites)[0]

        # make the exchange reaction with the boundary rather than ECM
        em_b = em_e.copy()
        em_b.compartment = 'b'
        em_b.id = '_'.join(em_b.id.split('_')[:-1]) + '_b'
        m_model.add_metabolites([em_b])

        er_b_metabolites = dict()
        for metab, coef in er_e.metabolites.items():
            if metab.id != em_e.id:
                er_b_metabolites[metab] = coef
            else:
                er_b_metabolites[em_b] = coef

        er_b = cobra.Reaction(id = '_'.join(er_e.id.split('_')[:-1]) + '_b', 
                             name = '_'.join(er_e.id.split('_')[:-1]) + '_b', 
                             lower_bound = er_e.lower_bound, upper_bound = er_e.upper_bound)
        m_model.add_reactions([er_b])
        er_b.add_metabolites(er_b_metabolites)
        new_exchanges.append(er_b)

        er_e_2 = cobra.Reaction(id = '_'.join(er_b.id.split('_')[:-1]) + '_LPAREN_e_RPAREN_', 
                               name = er_e.name, 
                               lower_bound = -1000, upper_bound = 1000)
        m_model.add_reactions([er_e_2])
        er_e_2.add_metabolites({em_e: -1, em_b: 1})
        new_exchanges.append(er_e_2)

        rm.append(er_e)
    m_model.remove_reactions(rm, remove_orphans=True)
    m_model.formatted_exchanges = new_exchanges # .exchanges is an attribute and can't be overwritten
