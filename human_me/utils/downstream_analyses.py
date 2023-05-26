#!/usr/bin/env python
# coding: utf-8

import os

from tqdm import tqdm
from typing import List
import warnings

import pandas as pd

from human_me.core.macromolecules.macromolecule import Macromolecule, Proxy
from human_me.core.biomass import Biomass
from human_me import io
from human_me.io import HiddenPrints

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
