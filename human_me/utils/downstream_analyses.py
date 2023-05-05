#!/usr/bin/env python
# coding: utf-8

import os

from tqdm import tqdm
from typing import List

import pandas as pd

from human_me.core.macromolecules.macromolecule import Macromolecule, Proxy
from human_me.core.biomass import Biomass
from human_me import io
from human_me.io import HiddenPrints

def _order_metabolites(metabolites_list: List[str], compartment_order: List[str] = None):
    """Order the metabolites to test in troubleshoot_me by compartment, with the first compartment being most likely 
    source of infeasability. 

    Parameters
    ----------
    metabolites_list : List[str]
        full list of metabolite IDs from the model
    compartment_order : List[str], optional
        ME Model compartments orderd by which are likely to cause feasibility issues, by default List[str]

    Returns
    -------
    _type_
        _description_
    """
    if compartment_order is None:
        compartment_order = ['b', 'e', 'c', 'n', 'm', 'i', 
                        'g', 'r', 'l', 'x', 'pm']
    
    ordered_model_metabolite_ids = {}
    for compartment in compartment_order:
        ordered_model_metabolite_ids[compartment] = sorted([m_id for m_id in metabolites_list if m_id.endswith('_' + compartment)])
    
    ordered_metabolites = []
    for compartment in compartment_order:
        ordered_metabolites += ordered_model_metabolite_ids[compartment]
    
    return ordered_metabolites

def troubleshoot_me(me_model_file: str, mu_val: float = 1e-9, out_path: str = None):
    """_summary_

    Parameters
    ----------
    me_model : str
        full path to pickled me model
    mu_val : float, optional
        growth rate to solve lp at, by default 1e-9
    out_path : str, optional
        where to save the output files, by default current working directory

    Returns
    -------
    sink_sln : pd.DataFrame
        the sink reactions and simulated fluxes required to make the ME Model feasible at mu
    metab_include_df : pd.DataFrame
        the metabolites that need sink solutions to make the ME Model feasible at mu
    """

    """Troubleshooting of feasibility. Iterates through all input metabolic model metabolites to identify necessariy sinks 
    to make model feasible at growth rate of mu_val. 

    Parameters
    ----------

    """

    me_model = io.read_pickled_me_model(me_model_file)

    if out_path is None:
        out_path = os.getcwd()

    # get all the metabolic model metabolites
    model_metabolites = []
    for metab in me_model.metabolites:
        if not isinstance(metab, Macromolecule):
            if not isinstance(metab, Proxy):
                if not isinstance(metab, Biomass):
                    model_metabolites.append(metab)
    model_metabolite_ids = [m.id for m in model_metabolites]

    for metab in tqdm(model_metabolites):
        me_model.add_boundary(metab, type = 'sink')

    with HiddenPrints():
        sln, stat, _ = me_model.solve_lp(mu_val = mu_val)
    if stat != 0:
        raise ValueError('The ME Model is not feasible at {:.2f} growth even when all metabolites are added as sinks'.format(mu_val))

    print('Initial flux blocking')

    n_excluded_metabolites_0 = -1
    n_excluded_metabolites = 0

    counter = 1

    metab_exclude_0, metab_exclude = [], []

    while (n_excluded_metabolites > n_excluded_metabolites_0) and stat == 0:
        print('Iteration: {} for initial flux solving'.format(counter))

        n_excluded_metabolites_0 = n_excluded_metabolites
        metab_exclude_0 = metab_exclude

        formatted_sln = me_model.format_solution(sln)

        sink_sln = formatted_sln[formatted_sln.reaction_id.apply(lambda x: x.startswith('SK_'))] # added before
        sink_sln = sink_sln[sink_sln.flux == 0]

        metab_exclude = sink_sln.reaction_id.apply(lambda x: x.split('SK_')[1]).tolist()
        metab_exclude = list(set(metab_exclude).intersection([m.id for m in model_metabolites]))
        n_excluded_metabolites = len(metab_exclude)
        print('{} of {} sink reactions can already be blocked'.format(n_excluded_metabolites, len(model_metabolites)))

        for metabolite_id in metab_exclude:
            rxn = me_model.reactions.get_by_id('SK_' + metabolite_id)
            rxn._lower_bound, rxn._upper_bound = 0,0

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
    metab_include = _order_metabolites(metab_include)[::-1] # backward to test expected most essential last

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

            rxn._lower_bound, rxn._upper_bound = 0,0

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

                if len(metab_exclude_additional) >  0:
                    metab_exclude_df = pd.concat([metab_exclude_df, 
                                                pd.DataFrame(data = {'metabolite_id': sorted(metab_exclude_additional), 'iteration': counter})], 
                                                ignore_index=True)

                metab_exclude = metab_exclude.union(metab_exclude_additional) # will avoid iterating through these 

            else: # reset reaction bounds to original
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