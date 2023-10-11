#!/usr/bin/env python
# coding: utf-8

import os

from tqdm import tqdm
from typing import List, Dict, Union
import warnings

import pandas as pd
import numpy as np
import cobra

from human_me.core.macromolecules.macromolecule import Macromolecule, Proxy
from human_me.core.biomass import Biomass
from human_me import io
from human_me.io import HiddenPrints
from human_me.preprocess.parse_complex import eval_complex
from human_me.utils.functions import flatten_list
from human_me.utils.machinery import rbps
from human_me.core.macromolecules.protein import Protein
from human_me.core.macromolecules.complex import Complex

e_aa = ['his_L', 'ile_L', 'leu_L', 'lys_L', 'met_L', 'phe_L', 'thr_L', 'trp_L', 'val_L']
ne_aa = ['ala_L', 'arg_L', 'asn_L', 'asp_L_c', 'cys_L', 'glu_L', 'gln_L', 'gly', 'pro_L', 'ser_L', 'tyr_L']

ordered_aa = e_aa + ne_aa # essential aas first

sugar_carbon_source_order = ['g1p', 'gal', 'fuc_L', 'fru', 'sucr', 'man', 'rib_D', 
                             'arab_L', 'inost', 'sbt_D', 'xyl_D', 'tre', 'pyr']
recon2_sugar_carbons = ['CE2838', 'CE2839', 'HC00229', 'HC00822', 'HC01440', 'HC01441', 'HC01446', 'adprbp', 'adprib', 
                        'arab_L', 'drib', 'fru', 'fuc_L', 'g1p', 'gal', 'glc_D', 'hom_L', 'lcts', 'malt', 'malthp', 'malthx', 'maltpt', 
                        'malttr', 'maltttr', 'man', 'nrvnc', 'rib_D', 'sucr', 'tagat_D', 'tre', 'udpg', 'udpgal', 'xyl_D']

other_carbons = sorted(set(recon2_sugar_carbons).difference(sugar_carbon_source_order))
ordered_sugar = sugar_carbon_source_order + other_carbons

def _order_metabolites(metabolites_list: List[str], me_model, 
                       order_aa: bool = True, order_sugar: bool = True, order_excreted: bool = True, 
                       compartment_order: List[str] = None, ):
    """Order the metabolites to test in troubleshoot_me by compartment, with the first compartment being most likely 
    source of infeasability. 

    Parameters
    ----------
    metabolites_list : List[str]
        full list of metabolite IDs from the model
    order_aa : bool, optional
        whether to sort amino acids in each compartment to be first in that list (most likely to influence 
        feasability)and within amino acids in the compartment, whether to put those that are essential 
        before those that are non-essential (more likely source of infeasability), by default True
    order_sugar: bool, optional
        whether to sort sugar carbon source by which are most important, by default True
    order_excreted: bool, optional
        whether to deprioritize exchanged metabolites that can only be excreted, by default True 
    compartment_order : List[str], optional
        ME Model compartments orderd by which are likely to cause feasibility issues, by default List[str]

    Returns
    -------
    ordered_metabolites : List[str]
        the input metabolites but sorted
    """
    if compartment_order is None:
        compartment_order = ['b', 'e', 'c', 'n', 'm', 'i', 'g', 'r', 'l', 'x', 'pm']

    ordered_model_metabolite_ids = {}
    for compartment in compartment_order:
        ordered_metabolites_compartment = sorted([m_id for m_id in metabolites_list if m_id.endswith('_' + compartment)])

        if order_aa:
            ordered_aa_compartment = [aa + '_' + compartment for aa in ordered_aa]
            aa_to_add = []
            for aa in ordered_aa_compartment:
                if aa in ordered_metabolites_compartment:
                    aa_to_add.append(aa)
            ordered_metabolites_compartment = [m_id for m_id in ordered_metabolites_compartment if m_id not in aa_to_add]
            ordered_metabolites_compartment = aa_to_add + ordered_metabolites_compartment

        if order_sugar:
            ordered_sugar_compartment = [sugar + '_' + compartment for sugar in ordered_sugar]
            sugar_to_add = []
            for sugar in ordered_sugar_compartment:
                if sugar in ordered_metabolites_compartment:
                    sugar_to_add.append(sugar)
            ordered_metabolites_compartment = [m_id for m_id in ordered_metabolites_compartment if m_id not in sugar_to_add]
            ordered_metabolites_compartment = sugar_to_add + ordered_metabolites_compartment

        ordered_model_metabolite_ids[compartment] = ordered_metabolites_compartment
    
    if order_excreted: # deprioritize adding sinks of metabolites that can only be excreted
        model_reactions = [r.cobra_id for r in me_model.reactions if hasattr(r, 'cobra_id')]
        secreted_metabs = []
        ordered_metabolites_b = ordered_model_metabolite_ids['b'].copy()

        for m_id in ordered_model_metabolite_ids['b']:
            if 'EX_' + m_id in model_reactions: #if the exchange reaction exists
                if me_model.reactions.get_by_id('EX_' + m_id).lower_bound >= 0:
                    secreted_metabs.append(m_id)
                    ordered_metabolites_b.remove(m_id)
        
        ordered_model_metabolite_ids['b'] = ordered_metabolites_b
 
    ordered_metabolites = []
    for compartment in compartment_order:
        ordered_metabolites += ordered_model_metabolite_ids[compartment]
    
    if order_excreted:
        ordered_metabolites += secreted_metabs
    
    return ordered_metabolites


def troubleshoot_me(me_model_file: str, mu_val: float = 1e-9, model_metabolite_ids: List[str] = None,
                    out_path: str = None, *args, **kwargs):
    """Helps troubleshoot ME Model if solver finds infeasible solution at growth rate of "mu_val". Basic principle is to 
    iterate through the metabolites from the input metabolic model, add them as sinks, and see which are required as sinks to 
    make the model feasible.

    *Note, this results in identified metabolites required for feasability being dependent on the order they were iterated 
    through. The metabolite list is run through _order_metabolites() to prioritize metabolites of certain types and
    in certain compartments.

    Parameters
    ----------
    me_model : str
        full path to pickled me model
    mu_val : float, optional
        growth rate to solve lp at, by default 1e-9
    metab_ids : List[str], optional
        a list of metabolite IDs (subset of all) to test whether adding the sinks makes the model feasible, by default all metabolites
    out_path : str, optional
        where to save the output files, by default current working directory
    *args : 
        into "_order_metabolites" function
    *kwargs : 
        into "_order_metabolites" function

    Returns
    -------
    sink_sln : pd.DataFrame
        the sink reactions and simulated fluxes required to make the ME Model feasible at mu
    metab_include_df : pd.DataFrame
        the metabolites that need sink solutions to make the ME Model feasible at mu
    """
    me_model = io.read_pickled_me_model(me_model_file)

    if out_path is None:
        out_path = os.getcwd()

    # get all the metabolic model metabolites
    if model_metabolite_ids is None:
        model_metabolites = []
        for metab in me_model.metabolites:
            if not isinstance(metab, Macromolecule):
                if not isinstance(metab, Proxy):
                    if not isinstance(metab, Biomass):
                        model_metabolites.append(metab)
        model_metabolite_ids = [m.id for m in model_metabolites]
    else:
        model_metabolites = [me_model.metabolites.get_by_id(m_id) for m_id in model_metabolite_ids]

    for metab in tqdm(model_metabolites):
        me_model.add_boundary(metab, type = 'sink')

    with HiddenPrints():
        sln, stat, _ = me_model.solve_lp(mu_val = mu_val)
    if counter == 1 and stat != 0:
        raise ValueError('The ME Model is not feasible at {:.2f} growth even when all metabolite IDs are added as sinks'.format(mu_val))

    print('Initial flux blocking')

    n_excluded_metabolites_0 = -1
    n_excluded_metabolites = 0
    stat = 0

    counter = 1

    metab_exclude_0, metab_exclude = [], []

    while (n_excluded_metabolites > n_excluded_metabolites_0) and stat == 0:
        print('Iteration: {} for initial flux solving'.format(counter))

        n_excluded_metabolites_0 = n_excluded_metabolites
        metab_exclude_0 = metab_exclude

        formatted_sln = me_model.format_solution(sln)

        sink_sln = formatted_sln[formatted_sln.reaction_id.apply(lambda x: x.startswith('SK_'))]  # added before
        sink_sln = sink_sln[sink_sln.flux == 0]

        metab_exclude = sink_sln.reaction_id.apply(lambda x: x.split('SK_')[1]).tolist()
        metab_exclude = list(set(metab_exclude).intersection([m.id for m in model_metabolites]))
        n_excluded_metabolites = len(metab_exclude)
        print('{} of {} sink reactions can already be blocked'.format(n_excluded_metabolites, len(model_metabolites)))

        for metabolite_id in metab_exclude:
            rxn = me_model.reactions.get_by_id('SK_' + metabolite_id)
            rxn._lower_bound, rxn._upper_bound = 0, 0

        with HiddenPrints():
            sln, stat, _ = me_model.solve_lp(mu_val = mu_val)

        counter += 1

    with open(os.path.join(out_path, 'metab_exclude_0.txt'), 'w') as f:
        for line in metab_exclude_0:
            f.write(f"{line}\n")

    print('Iterate through individual metabolite sinks')
    me_model = io.read_pickled_me_model(me_model_file)
    metab_exclude_df = pd.DataFrame(data = {'metabolite_id': metab_exclude_0})
    metab_exclude_df['iteration'] = 0

    metab_include = list(set(model_metabolite_ids).difference(metab_exclude_df.metabolite_id))

    # sort the metabolites since this algorithm depends on order you iterate through
    # setting those that are simplest/expected to be most essential to be tested last
    # that way, if earlier metabolites show up, we know it's not an issue of transport or exchange
    metab_include = _order_metabolites(metab_include, me_model, *args, **kwargs)[::-1] # backward to test expected most important last

    for metab_id in tqdm(metab_include):
        me_model.add_boundary(me_model.metabolites.get_by_id(metab_id), type = 'sink')
        
    metab_include_df = pd.DataFrame(columns = ['metabolite_id', 'iteration'])

    metab_exclude = set()
    counter = 1

    for metab_id in metab_include:
        print('{} of {} metabolites tested'.format(counter, len(metab_include)))
        
        if metab_id not in metab_exclude:
            rxn = me_model.reactions.get_by_id('SK_' + metab_id)
            lb, ub = rxn.bounds

            rxn._lower_bound, rxn._upper_bound = 0, 0

            with HiddenPrints():
                sln, stat, _ = me_model.solve_lp(mu_val = mu_val)

            if stat == 0: # if feasible, doesn't need this reaction
                metab_exclude_df.loc[metab_exclude_df.shape[0], :] = metab_id, counter

                # check for any additional 0 flux sinks from this iteration -- hopefully speeds up iterations
                formatted_sln = me_model.format_solution(sln)

                sink_sln = formatted_sln[formatted_sln.reaction_id.apply(lambda x: x.startswith('SK_'))] # added before
                sink_sln = sink_sln[sink_sln.flux == 0]
                metab_exclude_additional = sink_sln.reaction_id.apply(lambda x: x.split('SK_')[1]).tolist()
                metab_exclude_additional = set(metab_exclude_additional).intersection(model_metabolite_ids)
                metab_exclude_additional = metab_exclude_additional.difference(metab_exclude_df.metabolite_id.tolist())

                for metab_id_additional in metab_exclude_additional:
                    rxn = me_model.reactions.get_by_id('SK_' + metab_id_additional)
                    rxn._lower_bound, rxn._upper_bound = 0, 0

                if len(metab_exclude_additional) >  0:
                    metab_exclude_df = pd.concat([metab_exclude_df, 
                                                pd.DataFrame(data = {'metabolite_id': sorted(metab_exclude_additional), 'iteration': counter})], 
                                                ignore_index=True)

                metab_exclude = metab_exclude.union(metab_exclude_additional) # will avoid iterating through these 

            else:  # reset reaction bounds to original
                rxn._lower_bound, rxn._upper_bound = lb, ub
                metab_include_df.loc[metab_include_df.shape[0], :] = metab_id, counter

            metab_exclude_df.to_csv(os.path.join(out_path, 'excluded_metabolites.csv'))
            metab_include_df.to_csv(os.path.join(out_path, 'excluded_metabolites.csv'))

            counter += 1
    
    me_model = io.read_pickled_me_model(me_model_file)

    for metab_id in metab_include_df.metabolite_id:
        metab = me_model.metabolites.get_by_id(metab_id)
        me_model.add_boundary(metab, type = 'sink')

    sln, stat, _ = me_model.solve_lp(mu_val = mu_val)
    formatted_sln = me_model.format_solution(sln)
    sink_sln = formatted_sln[formatted_sln.reaction_id.apply(lambda x: x.startswith('SK_'))]

    return sink_sln, metab_include_df


def get_limiting_nutrients(me_model_file: str, max_feasible_mu: float, uptake_only: bool = False):
    """Identifies a list of limiting extracellular nutrients that prevent the ME Model from growing at 2.5% higher growth rate

    Parameters
    ----------
    me_model_file : str
        full path to pickled me model
    max_feasible_mu : float
        the maximum growth rate at which the ME Model is still feasible (see ME_Model.maximize_growth method)
    uptake_only : bool, optional
        whether to only consider extracellular nutrients that can be taken up (e.g., exclude those that can
        only be excreted according to exchange reaction bounds), by default False

    Returns
    -------
    see "troubleshoot_me" function

    """
    mu_val = max_feasible_mu * 1.025
    me_model = io.read_pickled_me_model(me_model_file)

    boundary_metabolite_ids = {m.id for m in me_model.metabolites if m.compartment == 'b'}
    exch_metab_ids = {'_'.join(r.id.split('EX_')[1:]) for r in me_model.reactions if r.id.startswith('EX_') and \
                        not r.id.endswith('LPAREN_e_RPAREN_')}

    if len(exch_metab_ids.difference(boundary_metabolite_ids)) > 0 or len(boundary_metabolite_ids.difference(exch_metab_ids)) > 0:
        warnings.warn('Unexpected discordance between boundary metabolites and exchanged metbaolites')

    metab_ids = exch_metab_ids.intersection(boundary_metabolite_ids)

    if uptake_only:
        exclude_metabs = []
        for metab_id in metab_ids:
            if me_model.reactions.get_by_id('EX_' + metab_id).lower_bound >= 0:
                exclude_metabs.append(metab_id)
        metab_ids = metab_ids.difference(exclude_metabs)
    metab_ids = sorted(metab_ids)

    sink_sln, metab_include_df = troubleshoot_me(model_file = me_model_file, mu_val = mu_val, model_metabolite_ids = metab_ids)

    return sink_sln, metab_include_df

def format_reaction_as_metabolic(reaction, expression_module: bool = False) -> Dict[str, Union[str, int]]:
    """Returns a dictionary mapping a MEReaction to its respective metabolic model reaction, if it exists

    Parameters
    ----------
    reaction : cobra.core.reaction.ME_Reaction
        the ME Model reaction to format
    expression_module : bool
        whether to also parse ExpressionReactions, which are not necessary for comparison with the metabolic model, by default False

    Returns
    -------
    metabolic_format : Dict[str]
        a dictionary mapping the MEReaction to its respective metabolic reaction
            "metabolic_reaction_id": the metabolic model reaction ID
            "reaction_no": if the reaction has an OR GPR, the ME Model building will have created one reaction for each GPR and separted them out by a number
            "reaction_direction": if the reaction was reversible and not spontaneous, the ME Model will have separated out the forward and reverse directions
    """
    if not expression_module and not hasattr(reaction, 'cobra_id'):
        metabolic_format = {'metabolic_reaction_id': np.nan, 'reaction_number': np.nan, 
                           'reaction_direction': np.nan}
    else:
        reaction_id = reaction.id
        if hasattr(reaction, 'cobra_id') and reaction_id == reaction.cobra_id:
            cobra_id = reaction_id
            reaction_no = 0
            reaction_direction = 'F'
        else:
            reaction_id_split = reaction_id.split('_')

            if reaction_id_split[-1].isdigit(): # when minimal_proteome = False in building
                reaction_no = int(reaction_id_split[-1])
                reaction_id_split = reaction_id_split[:-1]
            else:
                reaction_no = 0

            if reaction_id_split[-1] in ['F', 'R']:
                reaction_direction = reaction_id_split[-1]
                reaction_id_split = reaction_id_split[:-1]
            else:
                reaction_direction = 'F' # works for reversible reactions that weren't split too (doesn't have GPRs)

            cobra_id = '_'.join(reaction_id_split)

        metabolic_format = {'metabolic_reaction_id': cobra_id, 'reaction_number': reaction_no, 
                           'reaction_direction': reaction_direction}
    return metabolic_format

def _format_reaction_id_df(reactions: List, expression_module: bool = False):
    """Create a DataFrame mapping the ME Model reaction IDs to their respective metabolic model reaction IDs

    Parameters
    ----------
    reaction_ids :List[MEReaction]
        a list of ME Model reactions 
    expression_module : bool
        whether to also parse ExpressionReactions, which are not necessary for comparison with the metabolic model, by default False

     Returns
    -------
    m_reaction_df : pd.DataFrame
        dataframe with the following columns:
            "reaction_id": the ME Model reaction ID
            "cobra_id": the metabolic model reaction ID
            "reaction_no": if the reaction has an OR GPR, the ME Model building will have created one reaction for each GPR and separted them out by a number
            "reaction_direction": if the reaction was reversible and not spontaneous, the ME Model will have separated out the forward and reverse directions
    """

    m_reaction_df = pd.DataFrame(index = range(len(reactions)),
                             columns = ['cobra_id', 'reaction_no', 'reaction_direction'])
    m_reaction_df.insert(0, 'reaction_id', [r.id for r in reactions])
    m_reaction_df['cobra_id'], m_reaction_df['reaction_no'], \
    m_reaction_df['reaction_direction'] = zip(*pd.Series(reactions).apply(lambda x: \
                                                                          format_reaction_as_metabolic(x, expression_module = expression_module).values()))

    m_reaction_df = m_reaction_df[m_reaction_df.cobra_id.notna()]

    rev_map = dict(m_reaction_df.groupby('cobra_id')['reaction_direction'].apply(lambda x: x.eq('R').any()))
    m_reaction_df['reversible'] = m_reaction_df.cobra_id.map(rev_map).tolist()
    m_reaction_df = m_reaction_df[m_reaction_df.cobra_id.notna()]

    return m_reaction_df

def format_sln_as_metabolic(me_model, me_sln: pd.DataFrame):
    """For metabolic reactions, formats and aggregates (e.g. across reversibility or multiple reactions due to "OR" GPRs) 
    the ME Model solution to be comparable to the metabolic model solution. This should then be readily compared to the 
    output of m_model.optimize().to_frame().
    
    We recommend using the cm_2 output from preprocess.correct_inputs.correct_model as the comparable metabolic model.
    Though not necessary, a more fair comparison may be to set the growth rates (mu) of the metabolic and me model to be
    the same. In the cobra.Model metabolic model this can be done by setting the growth reaction bounds to be mu. 
    In the me model, this can be specified with the "mu_val" parameter in the ME_Model.solve_lp method.

    Parameters
    ----------
    me_model : ME_Model
        the constructed ME Model
    me_sln : pd.DataFrame
        the solution to the ME Model LP (output of ME_Model.format_solution method)

    Returns
    -------
    me_metab_sln : pd.DataFrame
        dataframe containing the me model (column name "me_flux") flux solutions with formatted reaction IDs
    """

    me_metab_sln = me_sln.copy()
    m_reaction_df = _format_reaction_id_df(reactions = [me_model.reactions.get_by_id(r_id) for r_id in me_sln.reaction_id])
    me_metab_sln = pd.concat([me_metab_sln, m_reaction_df[['cobra_id', 'reaction_no', 'reaction_direction']]], axis = 1)
    me_metab_sln['me_flux'] = me_metab_sln.apply(lambda x: x.flux if x.reaction_direction == 'F' else -x.flux, axis = 1).tolist()

    me_metab_sln = me_metab_sln.groupby(['cobra_id', 'reaction_direction', 'reaction_no'])['me_flux'].sum()
    me_metab_sln = pd.DataFrame(me_metab_sln.groupby(['cobra_id']).sum())

    return(me_metab_sln)

def _get_expression_flux(me_model, me_sln: pd.DataFrame, hgnc_id: str, 
                         molecule_type: str, group_by: str = 'sum', 
                        consider_degradation: bool = True) -> float:
    """Calculate the flux through gene expression for a given gene.

    Parameters
    ----------
    me_model : 
        the ME Model
    me_sln : pd.DataFrame
        the solution to the ME Model LP (output of ME_Model.format_solution method)
    hgnc_id : str
        the hgnc ID of a gene expressed in the input me_model
    molecule_type : str
        one of 'mrna' or 'protein' to get transcriptional or translational fluxes
    group_by : str, optional
        if multiple reactions, aggregate fluxes by group_by as input to the "func" argument of  by default 'sum'
    consider_degradation : bool, optional
        whether to subtract the degradation fluxes from the synthesis fluxes or only consider synthesis fluxes, 
        by default True

    Returns
    -------
    tot_flux : float
        the net flux for expression of a given gene
    """

    _syn_key = 'synthesis' if molecule_type == 'mrna' else 'translation'

    synthesis_reactions = me_model.expressed_genes[hgnc_id].reactions['ExpressionReactions'][molecule_type][_syn_key]
    degradation_reactions = None
    if consider_degradation: 
        degradation_reactions = me_model.expressed_genes[hgnc_id].reactions['ExpressionReactions'][molecule_type]['sink']

    if molecule_type == 'mrna':
        synthesis_reactions, degradation_reactions = [synthesis_reactions], [degradation_reactions]

    tot_flux = me_sln.loc[me_sln.reaction_id.isin(synthesis_reactions)]['flux'].aggregate(func = group_by)
    if consider_degradation:
        tot_flux -= me_sln.loc[me_sln.reaction_id.isin(degradation_reactions)]['flux'].aggregate(func = group_by)
        
    return tot_flux

def _drop_or_expression_fluxes(me_model, flux_df):
    """If multiple reactions were created due to OR GPR, between those multiple reactions, 
    will only retain the reaction whose catalyzing enzymes' genes have maximal expression flux across the 
    multiple reactions.
    """
    # get the genes that will definitely be retained (participate in a AND or monomeric reaction)
    not_or_genes = []
    for reaction in tqdm(me_model.reactions):
        if hasattr(reaction, 'cobra_gpr'): # MetabolicReaction
            gpr = reaction.cobra_gpr
        else: # Expressionreaction
            gpr = reaction.gene_reaction_rule

        if gpr != '':
            parsed_gpr = flatten_list([g if (type(g) == list) else [g] for g in eval_complex(gpr)])
            if 'or' not in gpr:
                not_or_genes += parsed_gpr
    #             all_genes += parsed_gpr
    #         else: # complexes
    #             all_genes += parsed_gpr

    not_or_genes += rbps # ribosomal proteins aren't in gprs ("ribosome") 
    not_or_genes.remove('ribosome')
    not_or_genes += ['HGNC:12468', 'HGNC:12463']   # ubiquitin genes
    not_or_genes = list(set(not_or_genes))


    # get the total expression flux for genes participating in OR reactions
    # and only retain those genes that have the maximum expression flux between all reactions 
    m_reaction_df = _format_reaction_id_df(reactions = me_model.reactions, 
                                          expression_module = True)

    # filter for reactions with OR gprs
    or_reactions = m_reaction_df[m_reaction_df.reaction_no > 0].cobra_id
    m_reaction_df = m_reaction_df[m_reaction_df.cobra_id.isin(or_reactions)]

    # add the catalyzing genes for each reaction
    respective_genes = []
    for reaction_id in m_reaction_df.reaction_id:
        reaction = me_model.reactions.get_by_id(reaction_id)
        cm = [k for k,v in reaction.coupled_metabolites.items() if v == 'catalysis']
        if len(cm)>1:
            print(reaction_id)
            raise ValueError('Unexpecetd multiple catalyzing enzymes')
        if isinstance(cm[0], Protein):
            gene = [cm[0].id.split('_')[0]]
        elif isinstance(cm[0], Complex):
            gene = [g.id.split('_')[0] for g in cm[0].components]
        respective_genes.append(gene)
    m_reaction_df['catalyzing_genes'] = respective_genes

    # get total expression flux through all genes (sum for complexes) for a reaction
    flux_df.set_index('HGNC_ID', inplace = True)
    m_reaction_df['tot_expression_flux'] = m_reaction_df['catalyzing_genes'].apply(lambda x: flux_df.loc[x, flux_df.columns[0]].sum())

    # for each reaction derived due to an OR gpr, identify the one with the maximum flux going through it
    max_fluxes = m_reaction_df.groupby(['cobra_id']).tot_expression_flux.idxmax()
    m_reaction_df['max_flux'] = False
    m_reaction_df.loc[max_fluxes, :] = True
    ors_to_drop = flatten_list(m_reaction_df[~m_reaction_df.max_flux].catalyzing_genes.tolist())

    # retain those that participate as and only reactions in other reactions
    ors_to_drop = list(set(ors_to_drop).difference(not_or_genes))
    print('Dropping {} genes from the analysis that participate in OR reactions'.format(len(ors_to_drop)))
    flux_df = flux_df.loc[~flux_df.index.isin(ors_to_drop), ]
    flux_df.reset_index(inplace = True)
    return flux_df

def get_expression_fluxes(me_model, me_sln: pd.DataFrame,  
                         molecule_type: str, group_by: str = 'sum', 
                        consider_degradation: bool = True, 
                        drop_ors: bool = False) -> float:
    """Calculate the flux through gene expression for each gene in ME Model.

    Parameters
    ----------
    me_model : ME_Model
        the constructed ME Model
    me_sln : pd.DataFrame
        the solution to the ME Model LP (output of ME_Model.format_solution method)
    molecule_type : str
        one of 'mrna' or 'protein' to get transcriptional or translational fluxes
    group_by : str, optional
        if multiple reactions, aggregate fluxes by group_by as input to the "func" argument of  by default 'sum'
    consider_degradation : bool, optional
        whether to subtract the degradation fluxes from the synthesis fluxes or only consider synthesis fluxes, 
        by default True
    drop_ors : bool, optional
        If multiple reactions were created due to OR GPR, between those multiple reactions, will only retain the 
        reaction whose catalyzing enzymes' genes have maximal expression flux across the multiple reactions.

    Returns
    -------
    flux_df : pd.DataFrame
        the net flux for expression of each gene
    """

    expression_fluxes = dict()
    for hgnc_id in tqdm(me_model.expressed_genes):
        expression_fluxes[hgnc_id] = _get_expression_flux(me_model = me_model, me_sln = me_sln, hgnc_id = hgnc_id, 
                                 molecule_type = molecule_type, group_by = group_by, 
                                consider_degradation = consider_degradation)

    flux_df = pd.DataFrame(data = {'HGNC_ID': expression_fluxes.keys(), 
                        molecule_type + '_flux': expression_fluxes.values()})
    if drop_ors:
        flux_df = _drop_or_expression_fluxes(me_model = me_model, flux_df = flux_df)
    return flux_df

def _fill_by_bounds(fva_df, sln):
    """If any other reactions in the solution are at their boundary, add these to avoid iterating through them."""
    fva_df = fva_df.copy()
    # skip some iterations if possible 
    reactions_min = fva_df[fva_df.minimum.isna()].index
    reactions_min = reactions_min[np.where(sln.loc[reactions_min, 'flux'] == fva_df.loc[reactions_min, 'min_bound'])[0]]
    fva_df.loc[reactions_min, 'minimum'] = fva_df.loc[reactions_min, 'min_bound']

    reactions_max = fva_df[fva_df.maximum.isna()].index
    reactions_max = reactions_max[np.where(sln.loc[reactions_max, 'flux'] == fva_df.loc[reactions_max, 'max_bound'])[0]]
    fva_df.loc[reactions_max, 'maximum'] = fva_df.loc[reactions_max, 'max_bound'] 
    
    return fva_df

def flux_variability_analysis(me_model, mu_val: float, reactions: List[str], n_cores: int = 0) -> pd.DataFrame:
    """Runs FVA on reactions of interest at a given growth rate. 

    **Note: this currently runs FVA only for the primary objective of growth.

    Parameters
    ----------
    me_model : 
        the me model
    mu_val : float
        the growth rate at which to run the LP (should be <= maximum feasible growth rate)
    reactions : List[str]
        a list of reaction IDs from the ME Model. Solving the ME Model takes much longer than a standard
        metabolic model, so we highly recommend limiting this list to a few reactions. 
    n_cores : int, optional
        number of cores to parallelize with, by default no parallelization. This may take longer than iterating
        through each one at a time (see use of _fill_by_bounds function when n_cores <= 1). 
        *Currently not implemented.

    Returns
    -------
    fva_df : pd.DataFrame
        first two columns represent the minimum and maximum possible fluxes at that growth rate
    """
    fva_df = pd.DataFrame(index = reactions, columns = ['minimum', 'maximum', 'solution_stats'])

    max_bounds, min_bounds = dict(), dict()
    for r_id in reactions:
        reaction_bounds = me_model.reactions.get_by_id(r_id).bounds
        max_bounds[r_id] = reaction_bounds[1]
        min_bounds[r_id] = reaction_bounds[0]
    fva_df['min_bound'] = fva_df.index.map(min_bounds)
    fva_df['max_bound'] = fva_df.index.map(max_bounds)

    if n_cores <= 1:
        for reaction in tqdm(reactions):
            if np.isnan(fva_df.loc[reaction, 'minimum']):
                sln, stat, _ = me_model.solve_lp(mu_val = mu_val, objective = {reaction: -1})
                sln = me_model.format_solution(sln)
                sln.set_index('reaction_id', inplace = True)
                fva_df.loc[reaction, ['minimum', 'solution_stats']] = [sln.loc[reaction, 'flux'], stat[()]]

                # skip some iterations by getting other reactions that are at their boundary
                fva_df = _fill_by_bounds(fva_df, sln)
            if np.isnan(fva_df.loc[reaction, 'maximum']):
                sln, stat, _ = me_model.solve_lp(mu_val = mu_val, objective = {reaction: 1})
                sln = me_model.format_solution(sln)
                sln.set_index('reaction_id', inplace = True)
                fva_df.loc[reaction, ['maximum', 'solution_stats']] = [sln.loc[reaction, 'flux'], stat[()]]

                # skip some iterations by getting other reactions that are at their boundary
                fva_df = _fill_by_bounds(fva_df, sln)
    else: #TODO
        raise ValueError('Internal: need to implement parallel solving; for now, set 0 <= n_cores <= 1')
    if (fva_df.solution_stats == 1).any():
        warnings.warn('An unexpected infeasible solution occured')

    return fva_df

def _format_reaction_bounds(x, type = 'minimum'):
    if type == 'minimum':
        if x.reversible:
            if x.reaction_direction == 'F':
                return(np.nan)
            else:
                return(-x.maximum)
        else:
            return(x.minimum)
    elif type == 'maximum':
        if x.reaction_direction == 'F':
            return(x.maximum)
        else: 
            return(np.nan)

def format_fva_as_metabolic(me_model, me_fva_df: pd.DataFrame, concat_type: str = 'lenient'):
    """Format the FVA to be comparable to the metabolic model reactions

    Parameters
    ----------
    me_model : ME_Model
        the constructed ME Model
    me_fva_df : pd.DataFrame
        the output of the "flux_variability_analysis" function
    concat_type : str, optional
        in the presence of OR GPRs, how to aggregate reaction bounds across the multiple reactions generated, by default 'lenient'
        options include:
            "lenient": will take the maximum upper bound and minimum lower bound across the multiple reactions
            "stringent": will take the minimum upper bound and maximum lower bound across the multiple reactions (note that this may cause situations in which min > max)
            func: input to pd.DataFrame.aggregate's "func" argument
    Returns
    -------
    comp_fva_df : pd.DataFrame
        A dataframe that should be directly comparable to the output of running 
        cobra.flux_analysis.variability.flux_variability_analysis on the m_model
    """
    comp_fva_df = me_fva_df.copy()
    m_reaction_df = _format_reaction_id_df(reactions = [me_model.reactions.get_by_id(r_id) for r_id in comp_fva_df.index])

    comp_fva_df.reset_index(drop = True, inplace = True)
    comp_fva_df = pd.concat([comp_fva_df[['minimum', 'maximum']], 
              m_reaction_df], axis = 1)

    comp_fva_df = comp_fva_df[['reaction_id', 'minimum', 'maximum', 'cobra_id', 'reaction_no', 
                               'reaction_direction', 'reversible']]

    comp_fva_df['minimum'] = comp_fva_df.apply(lambda x: _format_reaction_bounds(x, type = 'minimum'), axis = 1).tolist()
    comp_fva_df['maximum'] = comp_fva_df.apply(lambda x: _format_reaction_bounds(x, type = 'maximum'), axis = 1).tolist()

    if concat_type == 'lenient':
        comp_fva_df = pd.concat([comp_fva_df.groupby('cobra_id')['minimum'].min(), 
                   comp_fva_df.groupby('cobra_id')['maximum'].max()], axis = 1)
    elif concat_type == 'stringent':
        comp_fva_df = pd.concat([comp_fva_df.groupby('cobra_id')['minimum'].max(), 
               comp_fva_df.groupby('cobra_id')['maximum'].max()], axis = 1).min()
    else:
        comp_fva_df.groupby('cobra_id')['minimum'].aggregate(func = concat_type)

    return comp_fva_df