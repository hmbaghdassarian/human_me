#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import warnings
from tqdm import tqdm
import ast
import copy

import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np

import sys
sys.path.insert(1, '../../scripts/')
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from utils import parameters as params
    from utils import machinery as mach
    from utils import functions as func  
    
    from preprocess import parse_complex
    
    import expression.build_mrna_expression_reactions as build_mrna
    from expression import gene_information
    from expression.protein_expression import build_protein_expression_reactions as build_protein


# In[ ]:


def get_all_expression_reactions(hgnc_id, reactions, ub_args, psim = params.psim_me, machinery_list = mach.metabolic_machinery, 
                             compress_mrna = False, nonmachinery_locations = list(), stochastic = False, seed = None):
    '''Generates all the expression reactions for a given protein from the HGNC ID and the PSIM'''
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with func.HiddenPrints():
            gene_info = gene_information.generate_from_psim(hgnc_id, psim, machinery_list, reactions = reactions, nonmachinery_locations = nonmachinery_locations, 
                                                           stochastic = stochastic, seed = seed)
            mrna_reactions, mrna_transcript_c, mrna_deg_proxy  = build_mrna.get_mrna_expression_reactions(gene_info, compress_mrna = compress_mrna)
            protein_reactions, protein_metabolites = build_protein.get_protein_expression_reactions(gene_info, mrna_transcript_c, mrna_deg_proxy, 
                                                                                                    ub_args = ub_args)
            
    return mrna_reactions + protein_reactions, protein_metabolites

def emm_par(hgnc_id, gene_reaction_map, ub_args, compress_mrna, non_machinery, stochastic, seed):
    # None bc will add later for expression model specific to this
    nml = list()
    if hgnc_id in non_machinery:
        nml = non_machinery[hgnc_id]
    expr_reactions, protein_metabolites = get_all_expression_reactions(hgnc_id, 
                                        reactions = gene_reaction_map[hgnc_id], compress_mrna = compress_mrna, 
                                            ub_args = ub_args, nonmachinery_locations = nml, stochastic = stochastic, seed = seed)
    id_protein_map = {p.compartment: p for p in protein_metabolites} # store compartments and metabolite objects for each gene
    return id_protein_map, expr_reactions

def get_expression_machinery(reactions):
    gene_reaction_map = func.create_gene_reaction_map(reactions)
    expression_machinery_me = list(gene_reaction_map)
    if 'ribosome' in expression_machinery_me:
        expression_machinery_me.remove('ribosome')
    return gene_reaction_map, expression_machinery_me

def parse_complex_degradation_reaction_id(r_id):
    '''Generates universal reaction ID analagous to func.parse_me_reaction_id specifically for Complex_Degradation_Reaction'''
    if r_id.count('_COMPLEX_') == 1:
        return r_id[r_id.index('_COMPLEX_') + len('_COMPLEX_'):]
    else:
        return '_'.join(r_id.split('_')[[i for i in range(len(r_id.split('_'))) if r_id.split('_')[i] == 'COMPLEX'][-1]+1:])

def get_ko(mach, knock_out):
    '''Determine whether a machinery list intersects with a knock_out genes list'''
    if len(knock_out) == 0 or len(set(mach).intersection(knock_out))==0:
        return False
    else:
        return True

def get_complex_df(reactions, knock_out, stochastic):
    '''Generate the complex df
    
    Paramaters
    ----------
    reactions: list
        each element is a cobra.core.Reaction object
    knock_out: list
        each element is a string representing a gene expressed in the model which should be knocked out
    '''
    
    complex_df = pd.DataFrame(columns = ['reaction_id', 'compartment', 'machinery', 'is_complex', 'creates_multiple_reactions', 
                                            'knock_out'])

    for r in tqdm(reactions):
        compartment_ = r.enzyme_compartment

        ko = True
        #YOU ARE HERE
        if len(r.genes) == 1:
            ko = get_ko(mach = [list(r.genes)[0].id], knock_out = knock_out)
            complex_df.loc[complex_df.shape[0], :] = [r.id, compartment_, list(r.genes)[0].id, False, False, ko]
        elif 'and' in r.gene_reaction_rule and 'or' in r.gene_reaction_rule: 
            machinery_final = parse_complex.eval_complex(r.gene_reaction_rule)
            for m in machinery_final:
                if type(m) == list:
                    ko = get_ko(mach = m, knock_out = knock_out)
                    complex_df.loc[complex_df.shape[0], :] = [r.id, compartment_, ';'.join(m), True, True, ko]
                else:
                    ko = get_ko(mach = [m], knock_out = knock_out)
                    complex_df.loc[complex_df.shape[0], :] = [r.id, compartment_, m, False, True, ko]
        elif 'or' in r.gene_reaction_rule:
            machinery = [g.id for g in list(r.genes)]
            for m in machinery:
                ko = get_ko(mach = [m], knock_out = knock_out)
                complex_df.loc[complex_df.shape[0], :] = [r.id, compartment_, m, False, True, ko]
        elif 'and' in r.gene_reaction_rule:
            m = sorted([g.id for g in r.genes])
            ko = get_ko(mach = m, knock_out = knock_out)
            complex_df.loc[complex_df.shape[0], :] = [r.id, compartment_, ';'.join(m), True, False, ko]
        else: # no genes
            pass

    return complex_df


# In[ ]:


# functions for minimize proteome method
def explode(df_, explode_cols, sep = None, fill_value=float('nan'), preserve_index=False):
    '''Split entries with multiple values separated by a separator into multiple rows (or entries in lists)
    
    https://stackoverflow.com/questions/12680754/split-explode-pandas-dataframe-string-entry-to-separate-rows
    
    df_: pd.DataFrame
    explode_cols: list
        columns to split
    sep: separator
        if None, already in list format
    fill_value: 
        if empty entry, what to fill it with
    preserve_index: bool
        keep original index values for each new row
    
    '''
    res = None
    for col in explode_cols:
        if res is None:
            df = df_.copy()
        else:
            df = res.copy()
        
        df[col]=df[col].str.split(sep)

        idx_cols = df.columns.difference([col])
        # calculate lengths of lists
        lens = df[col].apply(lambda x: len(x))

        idx = np.repeat(df.index.values, lens)
        # create "exploded" DF
        res = (pd.DataFrame({
                    col_:np.repeat(df[col_].values, lens)
                    for col_ in idx_cols},
                    index=idx)
                 .assign(**{col_:np.concatenate(df.loc[lens>0, col_].values)
                                for col_ in [col]}))
        # append those rows that have empty lists
        if (lens == 0).any():
            # at least one list in cells is empty
            res = (res.append(df.loc[lens==0, idx_cols], sort=False)
                      .fillna(fill_value))
        # revert the original index order
        res = res.sort_index()
        # reset index if requested
        if not preserve_index:        
            res.reset_index(drop=True, inplace = True)
    return res

def map_machinery_compartment(df):
    map_ = df.groupby('machinery')['compartment'].apply(lambda x: list(set(x))).reset_index()
    map_ = dict(zip(map_.machinery, map_.compartment))
    return map_

def map_complex_machinery_compartment(df):
    df = df[df.is_complex & ~df.knock_out.astype(bool)][['compartment', 'machinery']]
    df = explode(df, explode_cols = ['machinery'], sep = ';')
    return map_machinery_compartment(df)

def merge_maps(map1, map2):
    merged_map = map1.copy()
    for k,v in map2.items():
        if k not in map1:
            merged_map[k] = v
        else:
            merged_map[k] = list(set(map1[k] + v))
    return merged_map

