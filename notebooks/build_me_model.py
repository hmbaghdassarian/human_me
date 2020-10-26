#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
from tqdm import tqdm
import warnings
import ast
import os
import pandas as pd
import copy

import sys
sys.path.insert(1, '../scripts/') # comment out in python script
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
#     from load_environmental_variables import *
#     from utils import * 
    from utils import parameters as params
    from utils import machinery as mach
    
    from utils import functions as func  
    
    with func.HiddenPrints():
        from utils.functions import COMPLEX
        from utils import utils_2

        import build_mrna_expression_reactions as build_mrna
        import build_protein_expression_reactions as build_protein

        from uniform_processes.build_ribosome_biogenesis_reactions import ribosomal_reactions, ribosome_complex_c
        from uniform_processes.build_trna_expression_reactions import trna_biogenesis_reactions
        from uniform_processes.biomass import biomass_reactions


# In[2]:


import time
start_time = time.time()


# In[3]:


# get pre-generated reactions
me_reactions = ribosomal_reactions + trna_biogenesis_reactions + build_protein.ub_reactions


# # Generate Protein Expression Reactions for All Machinery

# In[4]:


# non_machinery = list(set(params.psim_me.loc[params.psim_me.loc[:,'LOCATION'].dropna().index, 'HGNC_ID'].dropna().tolist()).difference(mach.metabolic_machinery + mach.expression_machinery))[:10]
# non_machinery += ['HGNC:3765']
non_machinery = []


# In[5]:


def get_all_expression_reactions(hgnc_id, psim = params.psim_me, machinery_list = mach.metabolic_machinery, 
                             metabolic_model = params.human_model):
    '''Generates all the expression reactions for a given protein from the HGNC ID and the PSIM'''
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with func.HiddenPrints():
            gene_info = utils_2.generate_geneinfo_object(hgnc_id, psim, machinery_list, metabolic_model)
            mrna_reactions, mrna_transcript_c, mrna_deg_proxy  = build_mrna.get_mrna_expression_reactions(gene_info)
            protein_reactions, protein_metabolites = build_protein.get_protein_expression_reactions(gene_info, mrna_transcript_c, mrna_deg_proxy)

    return mrna_reactions + protein_reactions, protein_metabolites


# In[6]:


# get protein expression for all metabolic reactions
print('Generate protein expression reactions for metabolic enzymes and non-machinery')
id_protein_map = dict() # map HGNC ID to a dictionary of compartments and cobra.Metabolite proteins

loop_machinery = mach.metabolic_machinery + non_machinery

for hgnc_id in tqdm(loop_machinery):
    # None bc will add later for expression model specific to this
    expr_reactions, protein_metabolites = get_all_expression_reactions(hgnc_id)
    id_protein_map[hgnc_id] = {p.compartment: p for p in protein_metabolites} # store compartments and metabolite objects for each gene

    me_reactions += expr_reactions


# In[7]:


# get protein expression reactions for all expression module reactions
print('Generate protein expression reactions for expression module enzymes')

# to generate the cobra.Model, must remove the mu values; need the cobra.Model to get reaction compartments from gene_info
# this code can be modified in the future to be faster by modifying gene info to take a reaction list rather than a cobra.Model
# or by using ME_Model class one it is implemented


# 1) initialize model
expression_module = cobra.Model('expression_module')

# 2) replace mu values
me_reactions_copy = copy.deepcopy(me_reactions) # this step is slow
idx = [i for i in range(len(me_reactions_copy)) if isinstance(me_reactions_copy[i], func.ME_Reaction)]
for i in idx:
     me_reactions_copy[i].replace_coefficient_mu(1)

# 3) create model and get machinery genes
expression_module.add_reactions(me_reactions_copy) # don't use utils.expression_machinery since not all machinery may be included
expression_machinery_me = [g.id for g in expression_module.genes]
if 'ribosome' in expression_machinery_me:
    expression_machinery_me.remove('ribosome')
# for internal use
if len(set(expression_machinery_me).difference(mach.expression_machinery)) > 0:
       raise ValueError('The expression module model contains unexpected machinery')
        

# expression_reactions = list()
for hgnc_id in tqdm(list(set(expression_machinery_me))):
    expr_reactions, protein_metabolites = get_all_expression_reactions(hgnc_id, machinery_list = expression_machinery_me,
                                          metabolic_model = expression_module)
    
    
    if hgnc_id not in set(expression_machinery_me).intersection(mach.metabolic_machinery):
        if hgnc_id in id_protein_map.keys():
            raise ValueError('Some genes not accounted for when generating metabolic machinery expression reactions')
        else:
            id_protein_map[hgnc_id] = {p.compartment:p for p in protein_metabolites}
    
    # when there is machinery overlap between metabolic and expression module, deal with compartment overlap 
    else:
        ids_to_keep = list(set([r.id for r in expr_reactions]).difference([r.id for r in me_reactions]))
        expr_reactions = [r for r in expr_reactions if r.id in ids_to_keep]
        
        temp_map = {p.compartment:p for p in protein_metabolites}
        for comp, met in temp_map.items():
            if comp not in id_protein_map[hgnc_id].keys(): 
                id_protein_map[hgnc_id][comp] = met

    me_reactions += expr_reactions
#     expression_reactions += expr_reactions


# The following cell continues to add any expression module machinery that may have arised from adding expression reactions for expression machinery. The simpler solution to above cell and below is to just use expression_machinery from utils.machinery rather than expression_machinery_me, but doing this would possibly create expression reactions for unused expression machinery (since expression_machinery is a list of ALL expression machinery in all possible reactions). This takes no iterations for recon2.2

# In[8]:


# get protein expression reactions for all expression module reactions
expression_module = cobra.Model('expression_module')

me_reactions_copy = copy.deepcopy(me_reactions) # this step is slow
idx = [i for i in range(len(me_reactions_copy)) if isinstance(me_reactions_copy[i], func.ME_Reaction)]
for i in idx:
     me_reactions_copy[i].replace_coefficient_mu(1)

expression_module.add_reactions(me_reactions_copy)
expression_machinery_me_2 = [g.id for g in expression_module.genes]
if 'ribosome' in expression_machinery_me_2:
    expression_machinery_me_2.remove('ribosome')

new_expression_machinery = list(set(expression_machinery_me_2).difference(expression_machinery_me + mach.metabolic_machinery))

counter = 1
while len(new_expression_machinery)>0:  # this condition leaves possibility that an existing machinery but with a new compartment is added and not accounted for
    print(counter)
    for hgnc_id in tqdm(list(set(expression_machinery_me_2))):
        expr_reactions, protein_metabolites = get_all_expression_reactions(hgnc_id, machinery_list = expression_machinery_me_2,
                                              metabolic_model = expression_module)


        if hgnc_id not in set(expression_machinery_me_2).intersection(expression_machinery_me + mach.metabolic_machinery):
            if hgnc_id in id_protein_map.keys():
                raise ValueError('Some genes not accounted for when generating metabolic machinery expression reactions')
            else:
                id_protein_map[hgnc_id] = {p.compartment:p for p in protein_metabolites}

        # when there is machinery overlap between metabolic and expression module, deal with compartment overlap 
        else:
            ids_to_keep = list(set([r.id for r in expr_reactions]).difference([r.id for r in me_reactions]))
            expr_reactions = [r for r in expr_reactions if r.id in ids_to_keep]

            temp_map = {p.compartment:p for p in protein_metabolites}
            for comp, met in temp_map.items():
                if comp not in id_protein_map[hgnc_id].keys(): 
                    id_protein_map[hgnc_id][comp] = met

        me_reactions += expr_reactions
    
    # get protein expression reactions for all expression module reactions
    expression_machinery_me = copy.deepcopy(expression_machinery_me_2)
    
    expression_module = cobra.Model('expression_module')
    me_reactions_copy = copy.deepcopy(me_reactions) # this step is slow
    idx = [i for i in range(len(me_reactions_copy)) if isinstance(me_reactions_copy[i], func.ME_Reaction)]
    for i in idx:
         me_reactions_copy[i].replace_coefficient_mu(1)
    expression_module.add_reactions(me_reactions_copy) # don't use utils.expression_machinery_2 since not all machinery may be included
    expression_machinery_me_2 = [g.id for g in expression_module.genes]
    if 'ribosome' in expression_machinery_me_2:
        expression_machinery_me_2.remove('ribosome')
    
    new_expression_machinery = list(set(expression_machinery_me_2).difference(expression_machinery_me + mach.metabolic_machinery))
    counter += 1
    
del expression_module


# In[9]:


expression_machinery_me = expression_machinery_me_2
# list of hgnc ids of machinery that overlap with metabolic reactions and expression reactions 
# but that are not used in the expression reactions for this specific model
excess_reactions = sorted(set(mach.metabolic_machinery).intersection(mach.expression_machinery).difference(set(mach.metabolic_machinery).intersection(expression_machinery_me)))
# filter for excess reactions with different protein localization for the expression reactions they catalyze
# than the metabolic reaction they catalyze
excess_reactions = [hgnc_id for hgnc_id in excess_reactions if len(id_protein_map[hgnc_id])>1]
if len(excess_reactions) > 0:
    raise ValueError('There are excess expression machinery reactions')


# # Complex Formation

# In[10]:


print('Get metabolic model complex information')
complex_df = pd.DataFrame(columns = ['reaction_id', 'compartment', 'machinery', 'is_complex', 'creates_multiple_reactions'])

for r in tqdm(params.human_model.reactions):
    compartment_ = func.get_reaction_compartment(r)
    if len(r.genes) == 1: 
        complex_df.loc[complex_df.shape[0], :] = [r.id, compartment_, list(r.genes)[0].id, False, False]
    elif 'and' in r.gene_reaction_rule and 'or' in r.gene_reaction_rule: 
        machinery_final = func.eval_complex(r.gene_reaction_rule)
        for m in machinery_final:
            if type(m) == list:
                complex_df.loc[complex_df.shape[0], :] = [r.id, compartment_, ';'.join(m), True, True]
            else:
                complex_df.loc[complex_df.shape[0], :] = [r.id, compartment_, m, False, True]
    elif 'or' in r.gene_reaction_rule:
        machinery = [g.id for g in list(r.genes)]
        for m in machinery:
            complex_df.loc[complex_df.shape[0], :] = [r.id, compartment_, m, False, True]
    elif 'and' in r.gene_reaction_rule:
        m = sorted([g.id for g in r.genes])
        complex_df.loc[complex_df.shape[0], :] = [r.id, compartment_, ';'.join(m), True, False]
    else:
        pass

complex_df['category'] = 'metabolic_reaction'


######------------------------------------------------
print('Get me reaction complex information')
me_complex_df = pd.DataFrame(columns = ['reaction_id', 'compartment', 'machinery', 'is_complex', 'creates_multiple_reactions'])

for r in tqdm(me_reactions):
    compartment_ = func.get_reaction_compartment(r)
    if len(r.genes) == 1: 
        me_complex_df.loc[me_complex_df.shape[0], :] = [r.id, compartment_, list(r.genes)[0].id, False, False]
    elif 'and' in r.gene_reaction_rule and 'or' in r.gene_reaction_rule: 
        machinery_final = func.eval_complex(r.gene_reaction_rule)
        for m in machinery_final:
            if type(m) == list:
                me_complex_df.loc[me_complex_df.shape[0], :] = [r.id, compartment_, ';'.join(m), True, True]
            else:
                me_complex_df.loc[me_complex_df.shape[0], :] = [r.id, compartment_, m, False, True]
    elif 'or' in r.gene_reaction_rule:
        machinery = [g.id for g in list(r.genes)]
        for m in machinery:
            me_complex_df.loc[me_complex_df.shape[0], :] = [r.id, compartment_, m, False, True]
    elif 'and' in r.gene_reaction_rule:
        m = sorted([g.id for g in r.genes])
        me_complex_df.loc[me_complex_df.shape[0], :] = [r.id, compartment_, ';'.join(m), True, False]
    else:
        pass

me_complex_df.reaction_id = me_complex_df.reaction_id.apply(lambda x: func. parse_me_reaction_id(x))
me_complex_df.drop_duplicates(keep = 'first', inplace = True) # if all but reaction id HGNC were the same
me_complex_df.reset_index(inplace = True, drop = True)
me_complex_df['category'] = 'expression_reaction'

######------------------------------------------------
# merge bc will deal with duplicate complexes, in case there is duplicates b/w metabolic and expression module
complex_df = pd.concat([complex_df, me_complex_df], axis = 0)
complex_df.reset_index(inplace = True, drop = True)


# In[11]:


# assign complex ids for reactions that have complexes in them
complex_df['complex_id'] = float('nan')
complex_df.loc[complex_df[complex_df.is_complex].index, 'complex_id'] = complex_df.loc[complex_df[complex_df.is_complex].index, 'reaction_id']

# creates one singular complex id for complexes shared across multiple reactions (same compartment, same machinery)
dup_complexes = complex_df[complex_df.is_complex].drop_duplicates(subset = ['compartment', 'machinery'], keep = 'first')
dup_complexes.reset_index(inplace = True, drop = True)
for i in dup_complexes.index:
    dups = complex_df[(complex_df.compartment == dup_complexes.loc[i,'compartment']) & (complex_df.machinery == dup_complexes.loc[i, 'machinery'])]
    complex_df.loc[dups.index,'complex_id'] = '_'.join(dups.reaction_id)


# In[12]:


# create a mapping of the unique complex_df ids to the actual complex metabolite
unique_complexes = complex_df[complex_df.is_complex]
unique_complexes = unique_complexes.drop_duplicates(subset = 'complex_id', keep = 'first')
unique_complexes.reset_index(inplace = True, drop = True)

complex_formation_reactions = list() # store all complex formation reactions
complex_id_metabolite_map = dict() # map complex id to the complex cobra.Metabolite
new_complex_ids = dict()

counter = 0
for i in unique_complexes.index:
    complex_id = unique_complexes.loc[i, 'complex_id']
    compartment = unique_complexes.loc[i, 'compartment']
    machinery = unique_complexes.loc[i, 'machinery'].split(';')
    
    if len(complex_id) > 256-8-4-(7*len(machinery)): # ids that are too long
        new_complex_ids[complex_id] = str(counter)
        complex_id = str(counter)
        counter += 1
    
    
    machinery_metabolites = list()
    metabolite_types = list()
    for m in machinery:
        if m != 'ribosome':
            machinery_metabolites.append(id_protein_map[m][compartment])
            metabolite_types.append('protein')
        else:
            machinery_metabolites.append(ribosome_complex_c)
            metabolite_types.append('complex')

    complex_info = {'METABOLITES': machinery_metabolites, 'IDS': [m.id for m in machinery_metabolites], 
                   'METABOLITE_TYPES': metabolite_types}
    complex_metabolite = COMPLEX(reaction_id = complex_id, complex_id = complex_id, **complex_info)
    complex_reaction = complex_metabolite.form_complex()

    complex_formation_reactions.append(complex_reaction)
    complex_id_metabolite_map[complex_id] = complex_metabolite

# ids that were too long
for k,v in new_complex_ids.items():
    complex_df.loc[complex_df[complex_df.complex_id == k].index, 'complex_id'] = v


# # Insert Machinery in Reactions

# In[13]:


# for reactions that show up more than once
reactions_to_track = params.human_model.reactions + me_reactions
reaction_counter = dict(zip(sorted(set([r.id for r in reactions_to_track])), [0]*len(reactions_to_track))) 


# In[14]:


# get SASA and keff values for coupling
print('Calculate enzyme k_effs')
# deal with metabolic reactions first
complex_df['MW_kDa'] = float('nan')
for i in tqdm(complex_df.index):
    if not complex_df.loc[i, 'is_complex']:
        enzyme_to_couple = id_protein_map[complex_df.loc[i, 'machinery']][complex_df.loc[i, 'compartment']]
    else:
        enzyme_to_couple = complex_id_metabolite_map[complex_df.loc[i, 'complex_id']]
    
    complex_df.loc[i, 'MW_kDa'] = func.get_metabolite_mw(enzyme_to_couple)
    
complex_df['SASA'] = complex_df.MW_kDa.apply(lambda x: func.SASA(x))
median_SASA = complex_df.SASA.median()
complex_df['keff'] = complex_df['SASA'].apply(lambda x: x*(params.keff_median/median_SASA))

# # for analysis
# from utils.load_environmental_variables import local_data_path
# complex_df.to_csv(local_data_path + 'interim/SASA_recon2_2.csv')


# In[15]:


print('Add machinery to metabolic module reactions')
metabolic_reactions = [r.id for r in params.human_model.reactions]
final_reactions = complex_formation_reactions # ALL reactions list, to create the ME model

# deal with metabolic reactions first
for i in tqdm(complex_df[complex_df.category == 'metabolic_reaction'].index):
    reaction_id = complex_df.loc[i, 'reaction_id'] # original reaction id
    r_ = params.human_model.reactions.get_by_id(reaction_id).copy() # original reaction

    r = func.ME_Reaction(type_ = ['catalysis'], 
                    id = r_.id, name = r_.name, subsystem = r_.subsystem, lower_bound = r_.lower_bound, 
                    upper_bound = r_.upper_bound, 
                        cobra_id = r_.id)
    r.add_metabolites(r_.metabolites)
    r.gene_reaction_rule = r_.gene_reaction_rule
    
    metabolites = r.metabolites.copy() # original reaction metabolites

    if not complex_df.loc[i, 'is_complex']:
        enzyme_to_couple = id_protein_map[complex_df.loc[i, 'machinery']][complex_df.loc[i, 'compartment']]
    else:
        enzyme_to_couple = complex_id_metabolite_map[complex_df.loc[i, 'complex_id']]

    # add machinery to substrate side
    c3 = (params.mu + params.alpha_p)/complex_df.loc[i, 'keff']
    r.add_metabolites({enzyme_to_couple: -c3}, combine = True)
    
    
    if not r_.reversibility:
        reactions = [r]
    else: # add a forward and reverse reaction for reversible reactions
        r_f, r_r = r.copy(), r.copy()
        r_f.lower_bound, r_r.upper_bound = 0, 0
        r_f.id, r_r.id = r_f.id + '_F', r_r.id + '_R'
        reactions = [r_f, r_r]
    
    # if multiple of the same reaction with different machinery due to OR GPR, add a different id for each
    if complex_df.loc[i, 'creates_multiple_reactions']:
        for j in range(len(reactions)):
            r_ = reactions[j]
            r_.id = r_.id + '_' + str(reaction_counter[reaction_id])
            reactions[j] = r_
        if reaction_counter[reaction_id] == 0: # tracking that all metabolic reactions are added
            metabolic_reactions.remove(reaction_id)
        reaction_counter[reaction_id] += 1
    else:
        metabolic_reactions.remove(reaction_id) # tracking that all metabolic reactions are added
    final_reactions += reactions

# only reactions without machinery should be left
if sorted(metabolic_reactions) != sorted([r.id for r in params.human_model.reactions if len(r.genes) == 0]):
    raise ValueError('Not all metabolic reactions that require machinery have been accounted for')
final_reactions += [r.copy() for r in params.human_model.reactions if len(r.genes) == 0]


# In[16]:


# filter out metabolic reactions
complex_df = complex_df[complex_df.category == 'expression_reaction']
complex_df.reset_index(inplace = True, drop = True)

print('Add machinery to expression module reactions')
for rxn in tqdm([r__ for r__ in me_reactions if len(r__.genes) > 0]):
    reaction_id_short = func.parse_me_reaction_id(rxn.id) # abbreviated version
    reaction_id = rxn.id # original reaction id
    idx = complex_df[complex_df.reaction_id == reaction_id_short].index.tolist()
    
    if not isinstance(rxn, func.ME_Reaction):
        rxn_me = func.ME_Reaction(type_ = ['catalysis'], 
                        id = rxn.id, name = rxn.name, subsystem = rxn.subsystem, lower_bound = rxn.lower_bound, 
                        upper_bound = rxn.upper_bound)
        rxn_me.add_metabolites(rxn.metabolites)
        rxn_me.gene_reaction_rule = rxn.gene_reaction_rule
    else: # translation reactions
        rxn_me = func.ME_Reaction(type_ = rxn.type + ['catalysis'], 
                        id = rxn.id, name = rxn.name, subsystem = rxn.subsystem, lower_bound = rxn.lower_bound, 
                        upper_bound = rxn.upper_bound)
        rxn_me.add_metabolites(rxn.metabolites)
        rxn_me.gene_reaction_rule = rxn.gene_reaction_rule
    
    for i in idx:
        r = copy.deepcopy(rxn_me)
        metabolites = r.metabolites.copy() # original reaction metabolites

        if not complex_df.loc[i, 'is_complex']:
            metabolite_to_add = id_protein_map[complex_df.loc[i, 'machinery']][complex_df.loc[i, 'compartment']]
        else:
            metabolite_to_add = complex_id_metabolite_map[complex_df.loc[i, 'complex_id']]

        # add machinery to substrate side
        c3 = (params.mu + params.alpha_p)/complex_df.loc[i, 'keff']
        r.add_metabolites({metabolite_to_add: -c3}, combine = True) # combine true in case machinery and substrate are same (e.g., ribosome translating its own proetin - but only for non-complex proteins)


        if not rxn.reversibility:
            reactions = [r]
        else: # add a forward and reverse reaction for reversible reactions
            r_f, r_r = r.copy(), r.copy()
            r_f.lower_bound, r_r.upper_bound = 0, 0
            r_f.id, r_r.id = r_f.id + '_F', r_r.id + '_R'
            reactions = [r_f, r_r]

        # if multiple of the same reaction with different machinery due to OR GPR, add a different id for each
        if complex_df.loc[i, 'creates_multiple_reactions']:
            for j in range(len(reactions)):
                r_ = reactions[j]
                r_.id = r_.id + '_' + str(reaction_counter[reaction_id])
                reactions[j] = r_
            reaction_counter[reaction_id] += 1
        final_reactions += reactions
        
final_reactions += [r__ for r__ in me_reactions if len(r__.genes) == 0]


# In[17]:


final_reactions += biomass_reactions


# In[18]:


end_time = time.time()
print('Time to build: {} (hrs)'.format((end_time - start_time)/3600))


# In[19]:


# metabolic_reactions = [r.id for r in params.human_model.reactions]
# A,B,C,D = 0,0,0,0

# failed_A, failed_B = list(), list()
# err = False
# for r in tqdm(final_reactions):
#     if isinstance(r, func.ME_Reaction):
#         if r.cobra_id is None and len(r.check_mass_balance())>0 and r.type != ['biomass']:
#             err = True
#             A += 1
#             failed_A.append(r)
#         elif r.cobra_id is not None:
#             ogr = params.human_model.reactions.get_by_id(r.cobra_id).copy()
#             if (len([k for k in ogr.metabolites.keys() if k.elements is None]) == 0) and (r.check_mass_balance() != ogr.check_mass_balance()):
#                 err = True
#                 B += 1
#                 failed_B.append(r)
#     else:
#         if r.id in metabolic_reactions:
#             ogr = params.human_model.reactions.get_by_id(r.id).copy()
#             if (len([k for k in ogr.metabolites.keys() if k.elements is None]) == 0) and (r.check_mass_balance() != ogr.check_mass_balance()):
#                 err = True
#                 C += 1
#         elif len(r.check_mass_balance())>0:
#             err = True
#             D += 1
# if err:
#     raise ValueError('Not all expression module reactions are mass balanced') 


# In[20]:


# test if can add reactions to model, and try to solve
final_reactions_og = copy.deepcopy(final_reactions)
for r in tqdm(final_reactions):
    if isinstance(r, func.ME_Reaction):
        if 'biomass' not in r.type:
            r.replace_coefficient_mu(mu_val = 0.03)
        else:
            r.replace_bound_mu(mu_val = 0.03, inplace = True)


# In[21]:


print('Generate and save ME-Model')
me_model = cobra.Model('HUMAN_ME_MODEL')
me_model.add_reactions(final_reactions)
me_model.objective = {me_model.reactions.biomass_dilution: -1}
# cobra.io.json.save_json_model(me_model, local_data_path + 'processed/human_me_model.json')


# In[ ]:


print('Optimize')
print('Growth rate = 0.03 hr^-1')
start_time = time.time()
test = me_model.optimize()
end_time = time.time()
print('Time to optimize: {} (hrs)'.format((end_time - start_time)/3600))


# In[ ]:


# print('Slim optimize')
# start_time = time.time()
# me_model.slim_optimize()
# end_time = time.time()
# print('Time to slim optimize: {} (hrs)'.format((end_time - start_time)/3600))


# In[ ]:


# print('Optimize')
# start_time = time.time()
# test = me_model.optimize()
# end_time = time.time()
# print('Time optimize: {} (hrs)'.format((end_time - start_time)/3600))

