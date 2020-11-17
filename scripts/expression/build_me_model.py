#!/usr/bin/env python
# coding: utf-8

# In[76]:


import cobra
from tqdm import tqdm
import warnings
import ast
import os
import pandas as pd
import numpy as np
import copy
import time

import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from utils.load_environmental_variables import build_files_path
    from utils import parameters as params
    from utils import machinery as mach
    
    from utils import functions as func  
    
    with func.HiddenPrints():
        from macromolecules.complex import Complex
        from utils import utils_2

        import expression.build_mrna_expression_reactions as build_mrna
        import expression.build_protein_expression_reactions as build_protein

        from uniform_processes.build_ribosome_biogenesis_reactions import ribosomal_reactions, ribosome_complex_c
        from uniform_processes.build_trna_expression_reactions import trna_biogenesis_reactions
        from uniform_processes import biomass


# # Generate Protein Expression Reactions for All Machinery

# In[2]:


def get_all_expression_reactions(hgnc_id, psim = params.psim_me, machinery_list = mach.metabolic_machinery, 
                             metabolic_model = params.human_model, compress_mrna = False):
    '''Generates all the expression reactions for a given protein from the HGNC ID and the PSIM'''
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with func.HiddenPrints():
            gene_info = utils_2.generate_geneinfo_object(hgnc_id, psim, machinery_list, metabolic_model)
            mrna_reactions, mrna_transcript_c, mrna_deg_proxy  = build_mrna.get_mrna_expression_reactions(gene_info, compress_mrna = compress_mrna)
            protein_reactions, protein_metabolites = build_protein.get_protein_expression_reactions(gene_info, mrna_transcript_c, mrna_deg_proxy)

    return mrna_reactions + protein_reactions, protein_metabolites

def generate_expression_module(me_reactions):
    # 1) initialize model
    expression_module = cobra.Model('expression_module')

    # 2) replace mu values
    me_reactions_copy = copy.deepcopy(me_reactions) # this step is slow - can be faster with an inplace argument in .replace_coefficient_mu
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
    return expression_machinery_me, expression_module


# In[3]:


class me_builder():
    def __init__(self, non_machinery = [], psim_me = params.psim_me, human_model = params.human_model, 
                 compress_mrna = False):
        self.non_machinery = non_machinery
        self.psim_me = psim_me
        self.human_model = human_model
        # get pre-generated reactions
        self.me_reactions = ribosomal_reactions + trna_biogenesis_reactions + build_protein.ub_reactions
        # map HGNC ID to a dictionary of compartments and cobra.Metabolite proteins
        self.id_protein_map = dict() 
        self.complex_id_metabolite_map = dict() # map complex id to the complex cobra.Metabolite
        
        self.id_reactions_map = dict()
        self.complex_reactions_map = dict()
        
        self.compress_mrna = compress_mrna
    
    def express_metabolic_enzymes(self):
        # get protein expression for all metabolic reactions
        print('Generate protein expression reactions for metabolic enzymes and non-machinery')

        loop_machinery = mach.metabolic_machinery + self.non_machinery

        for hgnc_id in tqdm(loop_machinery):
            # None bc will add later for expression model specific to this
            expr_reactions, protein_metabolites = get_all_expression_reactions(hgnc_id, compress_mrna = self.compress_mrna)
            self.id_protein_map[hgnc_id] = {p.compartment: p for p in protein_metabolites} # store compartments and metabolite objects for each gene
            
            self.id_reactions_map[hgnc_id] = expr_reactions
            
            
            self.me_reactions += expr_reactions
    def express_expression_enzymes(self):
        
        #This method continues to add any expression module machinery that may have arisen from adding expression 
        #reactions for expression machinery. The simpler solution to above cell and below is to just use 
        #expression_machinery from utils.machinery rather than expression_machinery_me, but doing this would 
        #possibly create expression reactions for unused expression machinery (since expression_machinery is a 
        #list of ALL expression machinery in all possible reactions). This takes 0 iterations for recon2.2
        
        # get protein expression reactions for all expression module reactions
        print('Generate protein expression reactions for expression module enzymes')

        # to generate the cobra.Model, must remove the mu values; need the cobra.Model to get reaction compartments from gene_info
        # this code can be modified in the future to be faster by modifying gene info to take a reaction list rather than a cobra.Model
        # or by using ME_Model class one it is implemented

        expression_machinery_me, expression_module = generate_expression_module(self.me_reactions)
        
        for hgnc_id in tqdm(list(set(expression_machinery_me))):
            expr_reactions, protein_metabolites = get_all_expression_reactions(hgnc_id, machinery_list = expression_machinery_me,
                                                  metabolic_model = expression_module, compress_mrna = self.compress_mrna)


            if hgnc_id not in set(expression_machinery_me).intersection(mach.metabolic_machinery):
                if hgnc_id in self.id_protein_map.keys():
                    raise ValueError('Some genes not accounted for when generating metabolic machinery expression reactions')
                else:
                    self.id_protein_map[hgnc_id] = {p.compartment:p for p in protein_metabolites}
            # when there is machinery overlap between metabolic and expression module, deal with compartment overlap 
            else:
                ids_to_keep = list(set([r.id for r in expr_reactions]).difference([r.id for r in self.me_reactions]))
                expr_reactions = [r for r in expr_reactions if r.id in ids_to_keep]

                temp_map = {p.compartment:p for p in protein_metabolites}
                for comp, met in temp_map.items():
                    if comp not in self.id_protein_map[hgnc_id].keys(): 
                        self.id_protein_map[hgnc_id][comp] = met
            
            if hgnc_id not in self.id_reactions_map.keys():
                self.id_reactions_map[hgnc_id] = expr_reactions
            else:
                self.id_reactions_map[hgnc_id] += expr_reactions
                
            self.me_reactions += expr_reactions
            
        expression_machinery_me_2, expression_module = generate_expression_module(self.me_reactions)
        new_expression_machinery = list(set(expression_machinery_me_2).difference(expression_machinery_me + mach.metabolic_machinery))

        counter = 1
        while len(new_expression_machinery)>0:  # this condition leaves possibility that an existing machinery but with a new compartment is added and not accounted for
            print('No. iterations for new expression machinery: {}'.format(counter))
            for hgnc_id in tqdm(list(set(expression_machinery_me_2))):
                expr_reactions, protein_metabolites = get_all_expression_reactions(hgnc_id, machinery_list = expression_machinery_me_2,
                                                      metabolic_model = expression_module, compress_mrna = self.compress_mrna)


                if hgnc_id not in set(expression_machinery_me_2).intersection(expression_machinery_me + mach.metabolic_machinery):
                    if hgnc_id in self.id_protein_map.keys():
                        raise ValueError('Some genes not accounted for when generating metabolic machinery expression reactions')
                    else:
                        self.id_protein_map[hgnc_id] = {p.compartment:p for p in protein_metabolites}

                # when there is machinery overlap between metabolic and expression module, deal with compartment overlap 
                else:
                    ids_to_keep = list(set([r.id for r in expr_reactions]).difference([r.id for r in self.me_reactions]))
                    expr_reactions = [r for r in expr_reactions if r.id in ids_to_keep]

                    temp_map = {p.compartment:p for p in protein_metabolites}
                    for comp, met in temp_map.items():
                        if comp not in self.id_protein_map[hgnc_id].keys(): 
                            self.id_protein_map[hgnc_id][comp] = met
                
                if hgnc_id not in self.id_reactions_map.keys():
                    self.id_reactions_map[hgnc_id] = expr_reactions
                else:
                    self.id_reactions_map[hgnc_id] += expr_reactions

                self.me_reactions += expr_reactions

            # get protein expression reactions for all expression module reactions
            expression_machinery_me = copy.deepcopy(expression_machinery_me_2)
            expression_machinery_me_2, expression_module = generate_expression_module(self.me_reactions)
            new_expression_machinery = list(set(expression_machinery_me_2).difference(expression_machinery_me + mach.metabolic_machinery))
            counter += 1

        del expression_module
        
        expression_machinery_me = expression_machinery_me_2
        # list of hgnc ids of machinery that overlap with metabolic reactions and expression reactions 
        # but that are not used in the expression reactions for this specific model
        excess_reactions = sorted(set(mach.metabolic_machinery).intersection(mach.expression_machinery).difference(set(mach.metabolic_machinery).intersection(expression_machinery_me)))
        # filter for excess reactions with different protein localization for the expression reactions they catalyze
        # than the metabolic reaction they catalyze
        excess_reactions = [hgnc_id for hgnc_id in excess_reactions if len(self.id_protein_map[hgnc_id])>1]
        if len(excess_reactions) > 0:
            raise ValueError('There are excess expression machinery reactions')
    
    def express_dummy_protein(self):
        ups = pd.read_csv(build_files_path + 'dummy_protein_features.tab', sep = '\t', header = None, index_col = 0)
        ups_ = pd.DataFrame(columns = params.psim_me.columns)
        ups_.loc[0, 'PREMRNA_SEQ'], ups_.loc[0,'MRNA_SEQ'], ups_.loc[0,'PROTEIN_SEQ'] = ups.loc['premrna_seq', 1], ups.loc['mrna_seq', 1], ups.loc['protein_seq', 1]
        dummy_id = 'HGNC:DUMMY'
        ups_['HGNC_ID'], ups_['LOCATION'] = dummy_id, '[c]'
        dummy_reactions, dm = get_all_expression_reactions(hgnc_id = dummy_id, psim = ups_, machinery_list = [], 
                                                            metabolic_model = cobra.Model(''), compress_mrna = self.compress_mrna) 

        for r in dummy_reactions:
            if biomass.protein_ in r.metabolites.keys():
                rxn = r.metabolites.copy()
                rxn[biomass.unmodeled_protein_] = rxn[biomass.protein_]
                rxn[biomass.protein_] = 0
                r.add_metabolites(rxn, combine = False)
        self.me_reactions += dummy_reactions
        self.dummy_protein = dm[0]
        
    def get_complex_info(self):
        print('Get metabolic model complex information')
        complex_df = pd.DataFrame(columns = ['reaction_id', 'compartment', 'machinery', 'is_complex', 'creates_multiple_reactions'])

        for r in tqdm(self.human_model.reactions):
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

        for r in tqdm(self.me_reactions):
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
        
        print('Assign unique complex ids for unique machinery-compartment sets across all reactions')
        # assign complex ids for reactions that have complexes in them
        complex_df['complex_id'] = float('nan')
        complex_df.loc[complex_df[complex_df.is_complex].index, 'complex_id'] = complex_df.loc[complex_df[complex_df.is_complex].index, 'reaction_id']
        
        # if a reaction generates multiple complexes, make sure each complex has a unique ID
        crm_ = complex_df[(complex_df.creates_multiple_reactions) & (complex_df.is_complex)].reaction_id.unique()
        for crm in crm_:
            df = complex_df[(complex_df.reaction_id == crm) & (complex_df.is_complex)]
            if df.shape[0]>1: # reaction creates multiple complexes
                counter = 0
                for i in df.index:
                    complex_df.loc[i, 'complex_id'] = complex_df.loc[i, 'complex_id'] + '_' + str(counter)
                    counter += 1
        
        dup_complexes = complex_df[complex_df.is_complex].duplicated(subset = ['compartment', 'machinery'], keep = 'first')
        dup_complexes = complex_df.loc[dup_complexes.index[np.where(dup_complexes)]]
        for i in dup_complexes.index:
            dups = complex_df[(complex_df.compartment == dup_complexes.loc[i,'compartment']) & (complex_df.machinery == dup_complexes.loc[i, 'machinery'])]
            complex_df.loc[dups.index,'complex_id'] = '_'.join(dups.reaction_id)
        
        self.complex_df = complex_df    
        
    def generate_complex_reactions(self):
        # create a mapping of the unique self.complex_df ids to the actual complex metabolite
        unique_complexes = self.complex_df[self.complex_df.is_complex]
        unique_complexes = unique_complexes.drop_duplicates(subset = 'complex_id', keep = 'first')
        unique_complexes.reset_index(inplace = True, drop = True)

        complex_formation_reactions = list() # store all complex formation reactions
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
                    machinery_metabolites.append(self.id_protein_map[m][compartment])
                    metabolite_types.append('protein')
                else:
                    machinery_metabolites.append(ribosome_complex_c)
                    metabolite_types.append('complex')

            complex_info = {'METABOLITES': machinery_metabolites, 'IDS': [m.id for m in machinery_metabolites], 
                           'METABOLITE_TYPES': metabolite_types}
            complex_metabolite = Complex(reaction_id = complex_id, complex_id = complex_id, **complex_info)
            complex_reaction = complex_metabolite.form_complex()

            complex_formation_reactions.append(complex_reaction)
            self.complex_id_metabolite_map[complex_id] = complex_metabolite
            self.complex_reactions_map[complex_id] = complex_reaction

        # ids that were too long
        for k,v in new_complex_ids.items():
            self.complex_df.loc[self.complex_df[self.complex_df.complex_id == k].index, 'complex_id'] = v
        
        self.complex_formation_reactions = complex_formation_reactions
    def get_keff(self):
        # for reactions that show up more than once
        reactions_to_track = self.human_model.reactions + self.me_reactions
        self.reaction_counter = dict(zip(sorted(set([r.id for r in reactions_to_track])), [0]*len(reactions_to_track))) 

        # get SASA and keff values for coupling
        print('Calculate enzyme k_effs')
        # deal with metabolic reactions first
        self.complex_df['MW_kDa'] = float('nan')
        for i in tqdm(self.complex_df.index):
            if not self.complex_df.loc[i, 'is_complex']:
                enzyme_to_couple = self.id_protein_map[self.complex_df.loc[i, 'machinery']][self.complex_df.loc[i, 'compartment']]
            else:
                enzyme_to_couple = self.complex_id_metabolite_map[self.complex_df.loc[i, 'complex_id']]

            self.complex_df.loc[i, 'MW_kDa'] = enzyme_to_couple.formula_weight/1000 

        self.complex_df['SASA'] = self.complex_df.MW_kDa.apply(lambda x: func.SASA(x))
        median_SASA = self.complex_df.SASA.median()
        self.complex_df['keff'] = self.complex_df['SASA'].apply(lambda x: x*(params.keff_median/median_SASA))
    def minimize_proteome(self):
        c_og = self.complex_df.copy()
        n_reactions_og = len(self.me_reactions) + len(self.complex_formation_reactions)

        drop_index = list()
        reaction_multiple = self.complex_df[self.complex_df.creates_multiple_reactions].reaction_id.unique().tolist()

        for rm in reaction_multiple:
            df = self.complex_df[self.complex_df.reaction_id == rm]
            # don't directly drop machinery in case they are used in multiple reactions and are minimal in 
            # another one of those reactions
            to_drop = df[df.MW_kDa != df.MW_kDa.min()].index.tolist() 
            if df.shape[0] - len(to_drop) == 1:
                drop_index += to_drop
            else:
                raise ValueError('Something went wrong in selecting a complex by lowerst molecular weight')

        self.complex_df.drop(index = drop_index, inplace = True)

        # backtrack and remove all protein expression and complex formation reactions of dropped enzymes

        # get rid of redundant complexes
        complexes_to_drop = sorted(set(c_og[c_og.is_complex].complex_id).difference(self.complex_df.complex_id))#sorted(set(self.complex_id_metabolite_map.keys()).difference(self.complex_df.complex_id))
        complexes_to_drop_id = [self.complex_reactions_map[c_id].id for c_id in complexes_to_drop]
        self.complex_formation_reactions = [r for r in self.complex_formation_reactions if r.id not in complexes_to_drop_id]
        for c_id in complexes_to_drop:
            del self.complex_reactions_map[c_id]
            del self.complex_id_metabolite_map[c_id]

        # individual proteins is more complcated because we have to check by compartment and if they are 
        # used in complexes
        prot_to_drop = set(c_og[c_og.is_complex == False].machinery).difference(self.complex_df[self.complex_df.is_complex == False].machinery)
        reactions_to_remove = []
        id_protein_map = self.id_protein_map.copy()
        for hgnc_id in prot_to_drop:
            complex_machinery = [i.split(';') for i in self.complex_df[(self.complex_df.is_complex)].machinery.tolist()]
            complex_machinery = sorted(set([item for sublist in complex_machinery for item in sublist]))
            if hgnc_id not in complex_machinery:
                reactions_to_remove += [r.id for r in self.id_reactions_map[hgnc_id]]
                del self.id_reactions_map[hgnc_id] 
                del self.id_protein_map[hgnc_id]
            else:
                for comp in id_protein_map[hgnc_id].keys():
                    complex_machinery = [i.split(';') for i in self.complex_df[(self.complex_df.is_complex) & (self.complex_df.compartment == comp)].machinery.tolist()]
                    complex_machinery = sorted(set([item for sublist in complex_machinery for item in sublist]))

                    if hgnc_id not in complex_machinery:
                        rr = [r.id for r in self.id_reactions_map[hgnc_id] if len(r.compartments.intersection([comp])) > 0]
                        reactions_to_remove += rr
                        self.id_reactions_map[hgnc_id] = [r for r in self.id_reactions_map[hgnc_id] if r.id not in rr]
                        self.id_protein_map[hgnc_id] = {k:v for k,v in self.id_protein_map[hgnc_id].items() if k != comp}
        self.me_reactions = [r for r in self.me_reactions if r.id not in reactions_to_remove]
        n_reactions = len(self.me_reactions) + len(self.complex_formation_reactions)

        print('A total of {} reactions were dropped when forming a minimal proteome'.format(n_reactions_og - n_reactions))
    def add_metabolic_machinery(self):
        print('Add machinery to metabolic module reactions')
        metabolic_reactions = [r.id for r in params.human_model.reactions]
        final_reactions = self.complex_formation_reactions # ALL reactions list, to create the ME model

        # deal with metabolic reactions first
        for i in tqdm(self.complex_df[self.complex_df.category == 'metabolic_reaction'].index):
            reaction_id = self.complex_df.loc[i, 'reaction_id'] # original reaction id
            r_ = params.human_model.reactions.get_by_id(reaction_id).copy() # original reaction

            r = func.ME_Reaction(type_ = ['catalysis'], 
                            id = r_.id, name = r_.name, subsystem = r_.subsystem, lower_bound = r_.lower_bound, 
                            upper_bound = r_.upper_bound, 
                                cobra_id = r_.id)
            r.add_metabolites(r_.metabolites)
            r.gene_reaction_rule = r_.gene_reaction_rule

            metabolites = r.metabolites.copy() # original reaction metabolites

            if not self.complex_df.loc[i, 'is_complex']:
                enzyme_to_couple = self.id_protein_map[self.complex_df.loc[i, 'machinery']][self.complex_df.loc[i, 'compartment']]
            else:
                enzyme_to_couple = self.complex_id_metabolite_map[self.complex_df.loc[i, 'complex_id']]

            # add machinery to substrate side
            c3 = (params.mu + params.alpha_p)/self.complex_df.loc[i, 'keff']

            if not r_.reversibility:
                r.add_metabolites({enzyme_to_couple: -c3}, combine = True) # combine true in case machinery and substrate are same (e.g., ribosome translating its own proetin - but only for non-complex proteins)
                reactions = [r]
            else: # add a forward and reverse reaction for reversible reactions
                r_f,r_r = r.copy(), r.copy()
                r_f.lower_bound, r_r.lower_bound, r_r.upper_bound = 0,0, abs(r.lower_bound)
                r_r.add_metabolites({metab: -coeff for metab, coeff in r_r.metabolites.items()}, combine = False)

                r_f.add_metabolites({enzyme_to_couple: -c3}, combine = True) # combine true in case machinery and substrate are same (e.g., ribosome translating its own proetin - but only for non-complex proteins)
                r_r.add_metabolites({enzyme_to_couple: -c3}, combine = True)

                r_f.id, r_r.id = r_f.id + '_F', r_r.id + '_R'
                reactions = [r_f, r_r]

            # if multiple of the same reaction with different machinery due to OR GPR, add a different id for each
            if self.complex_df.loc[i, 'creates_multiple_reactions']:
                for j in range(len(reactions)):
                    r_ = reactions[j]
                    r_.id = r_.id + '_' + str(self.reaction_counter[reaction_id])
                    reactions[j] = r_
                if self.reaction_counter[reaction_id] == 0: # tracking that all metabolic reactions are added
                    metabolic_reactions.remove(reaction_id)
                self.reaction_counter[reaction_id] += 1
            else:
                metabolic_reactions.remove(reaction_id) # tracking that all metabolic reactions are added
            final_reactions += reactions

        # only reactions without machinery should be left
        if sorted(metabolic_reactions) != sorted([r.id for r in params.human_model.reactions if len(r.genes) == 0]):
            raise ValueError('Not all metabolic reactions that require machinery have been accounted for')
        final_reactions += [r.copy() for r in params.human_model.reactions if len(r.genes) == 0]

        self.final_reactions = final_reactions
    def add_expression_machinery(self):
        # filter out metabolic reactions
        self.complex_df = self.complex_df[self.complex_df.category == 'expression_reaction']
        self.complex_df.reset_index(inplace = True, drop = True)

        print('Add machinery to expression module reactions')
        for rxn in tqdm([r__ for r__ in self.me_reactions if len(r__.genes) > 0]):
            reaction_id_short = func.parse_me_reaction_id(rxn.id) # abbreviated version
            reaction_id = rxn.id # original reaction id
            idx = self.complex_df[self.complex_df.reaction_id == reaction_id_short].index.tolist()

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

                if not self.complex_df.loc[i, 'is_complex']:
                    metabolite_to_add = self.id_protein_map[self.complex_df.loc[i, 'machinery']][self.complex_df.loc[i, 'compartment']]
                else:
                    metabolite_to_add = self.complex_id_metabolite_map[self.complex_df.loc[i, 'complex_id']]

                # add machinery to substrate side
                c3 = (params.mu + params.alpha_p)/self.complex_df.loc[i, 'keff']

                if not rxn.reversibility:
                    r.add_metabolites({metabolite_to_add: -c3}, combine = True) # combine true in case machinery and substrate are same (e.g., ribosome translating its own proetin - but only for non-complex proteins)
                    reactions = [r]
                else: # add a forward and reverse reaction for reversible reactions
                    r_f,r_r = r.copy(), r.copy()
                    r_f.lower_bound, r_r.lower_bound, r_r.upper_bound = 0,0, abs(r.lower_bound)
                    r_r.add_metabolites({metab: -coeff for metab, coeff in r_r.metabolites.items()}, combine = False)

                    r_f.add_metabolites({metabolite_to_add: -c3}, combine = True) # combine true in case machinery and substrate are same (e.g., ribosome translating its own proetin - but only for non-complex proteins)
                    r_r.add_metabolites({metabolite_to_add: -c3}, combine = True)

                    r_f.id, r_r.id = r_f.id + '_F', r_r.id + '_R'
                    reactions = [r_f, r_r]

                # if multiple of the same reaction with different machinery due to OR GPR, add a different id for each
                if self.complex_df.loc[i, 'creates_multiple_reactions']:
                    for j in range(len(reactions)):
                        r_ = reactions[j]
                        r_.id = r_.id + '_' + str(self.reaction_counter[reaction_id])
                        reactions[j] = r_
                    self.reaction_counter[reaction_id] += 1
                self.final_reactions += reactions

        self.final_reactions += [r__ for r__ in self.me_reactions if len(r__.genes) == 0]
#     def check_me_mass_balance(self):
#         metabolic_reactions = [r.id for r in self.human_model.reactions]
#         err = False
#         for r in tqdm(self.final_reactions):
#             if isinstance(r, func.ME_Reaction):
#                 if r.cobra_id is None and len(r.check_mass_balance())>0 and r.type != ['biomass']:
#                     err = True
#                     break
#                 elif r.cobra_id is not None:
#                     ogr = self.human_model.reactions.get_by_id(r.cobra_id).copy()
#                     if (len([k for k in ogr.metabolites.keys() if k.elements is None]) == 0) and (r.check_mass_balance() != ogr.check_mass_balance()):
#                         err = True
#                         break
#             else:
#                 if r.id in metabolic_reactions:
#                     ogr = self.human_model.reactions.get_by_id(r.id).copy()
#                     if (len([k for k in ogr.metabolites.keys() if k.elements is None]) == 0) and (r.check_mass_balance() != ogr.check_mass_balance()):
#                         err = True
#                 elif len(r.check_mass_balance())>0:
#                     err = True
#         if err:
#             raise ValueError('Not all expression module reactions are mass balanced') 
    def build_me_model(self, model_id = 'HUMAN_ME_MODEL'):
        self.final_reactions += biomass.biomass_reactions
#         self.check_me_mass_balance()

        print('Generate ME-Model')
        me_model = func.ME_Model(model_id)
        me_model.add_reactions(self.final_reactions)

        return me_model


# In[ ]:





# In[4]:


def build_me(non_machinery = [], minimal_proteome = False, model_id = 'HUMAN_ME_MODEL', compress_mrna = False,
                  psim_me = params.psim_me, human_model = params.human_model):
    '''
    Returns a human ME_model. 
    
    Inputs:
        non_machinery: a list of HGNC_IDs of non_machinery proteins
        minimal_proteome: bool; For reactions with OR in the GPR, the builder by default (False) generates a 
        separate reaction for each protein complex (False). If True, builder instead will create one reaction, 
        choosing the protein complex with the lowest molecular weight to catalyze the reaction.
        model_id: string; id for the me model
        compress_mrna: boolean, whether to merge the 3 transcription, mrna processing, and mrna export to cytosol 
        reactions into one single reaction
    
    '''
    start = time.time()
    
    builder = me_builder(non_machinery = non_machinery, psim_me = psim_me, human_model = human_model, 
                        compress_mrna = compress_mrna)
    builder.express_metabolic_enzymes()
    builder.express_expression_enzymes()
    builder.express_dummy_protein()
    builder.get_complex_info()
    builder.generate_complex_reactions()
    builder.get_keff()
    if minimal_proteome:
        builder.minimize_proteome()
    builder.add_metabolic_machinery()
    builder.add_expression_machinery()
    me_model = builder.build_me_model(model_id = model_id)

    end = time.time()
    print('Time to build: {} minutes'.format((end-start)/60))


    return me_model, builder


# In[6]:


# non_machinery = []
# minimal_proteome = True
# model_id = 'HUMAN_ME_MODEL'
# compress_mrna = False,
# psim_me = params.psim_me
# human_model = params.human_model


# In[12]:


# builder = me_builder(non_machinery = non_machinery, psim_me = psim_me, human_model = human_model, 
#                     compress_mrna = compress_mrna)
# builder.express_metabolic_enzymes()
# builder.express_expression_enzymes()
# builder.get_complex_info()


# In[10]:


# builder_ = copy.copy(builder)

