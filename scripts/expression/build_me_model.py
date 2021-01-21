#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra

import os
import gc
import warnings
from tqdm import tqdm
import ast
import copy
import time

import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np


import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from utils.load_environmental_variables import build_files_path, processed_data_path
    from utils import parameters as params
    from utils import machinery as mach
    
    from utils import functions as func  
    from preprocess import parse_complex
    
    from core.reaction import ME_Reaction
    from core.model import ME_Model
    
    with func.HiddenPrints():
        from macromolecules.protein import Protein
        from macromolecules.complex import Complex
        
        import expression.build_mrna_expression_reactions as build_mrna
        from expression import gene_information
        from expression.protein_expression import ubiquitin
        from expression.protein_expression import build_protein_expression_reactions as build_protein
        
        from uniform_processes.build_ribosome_biogenesis_reactions import build_ribosome
        from uniform_processes.build_trna_expression_reactions import trna_biogenesis_reactions
        from uniform_processes import biomass


# # Generate Protein Expression Reactions for All Machinery

# In[2]:


def get_all_expression_reactions(hgnc_id, ub_args, psim = params.psim_me, machinery_list = mach.metabolic_machinery, 
                             metabolic_model = params.human_model, compress_mrna = False):
    '''Generates all the expression reactions for a given protein from the HGNC ID and the PSIM'''
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with func.HiddenPrints():
            gene_info = gene_information.generate(hgnc_id, psim, machinery_list, metabolic_model)
            mrna_reactions, mrna_transcript_c, mrna_deg_proxy  = build_mrna.get_mrna_expression_reactions(gene_info, compress_mrna = compress_mrna)
            protein_reactions, protein_metabolites = build_protein.get_protein_expression_reactions(gene_info, mrna_transcript_c, mrna_deg_proxy, 
                                                                                                    ub_args = ub_args)
            

    return mrna_reactions + protein_reactions, protein_metabolites

def generate_expression_module(me_reactions):
    # 1) initialize model
    expression_module = cobra.Model('expression_module')

    # 2) replace mu values
    me_reactions_copy = copy.deepcopy(me_reactions) # this step is slow - can be faster with an inplace argument in .replace_coefficient_mu
    idx = [i for i in range(len(me_reactions_copy)) if isinstance(me_reactions_copy[i], ME_Reaction)]
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
    def __init__(self, compress_mrna = False, dummy_protein = True, 
                psim_me = params.psim_me, m_model = params.human_model):
        
        self.non_machinery = open(processed_data_path + 'corrected_non_machinery.txt').read().splitlines()
        self.psim_me = psim_me
        self.m_model = m_model
        
        # get pre-generated reactions - the compress_mrna arg requires that they be run with that input
        self.compress_mrna = compress_mrna
        print('Generate ubiquitin reactions for proteasomal degradation')
        self.ub_args = ubiquitin.express_ubiquitin(compress_mrna = self.compress_mrna)
        print('Generate ribosome')
        ribosomal_reactions, self.ribosome_complex_c = build_ribosome(self.ub_args, self.compress_mrna )
        
        self.dummy_protein = dummy_protein
            
        
        self.deorphaned = None
        self.orphan = None
        
        self.me_reactions = [copy.deepcopy(r) for r in trna_biogenesis_reactions] + ribosomal_reactions + self.ub_args['ub_reactions']
        # map HGNC ID to a dictionary of compartments and cobra.Metabolite proteins
        self.id_protein_map = dict() 
        self.complex_id_metabolite_map = dict() # map complex id to the complex cobra.Metabolite
        
        self.id_reactions_map = dict()
        self.complex_reactions_map = dict()
    
    def express_metabolic_enzymes(self):
        # get protein expression for all metabolic reactions
        print('Generate protein expression reactions for metabolic enzymes and non-machinery')

        loop_machinery = mach.metabolic_machinery + self.non_machinery

        for hgnc_id in tqdm(loop_machinery):
            # None bc will add later for expression model specific to this
            expr_reactions, protein_metabolites = get_all_expression_reactions(hgnc_id, compress_mrna = self.compress_mrna, 
                                                                              ub_args = self.ub_args)
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
        print('Generate protein expression reactions for expression module enzymes, this step may take a few minutes')

        # to generate the cobra.Model, must remove the mu values; need the cobra.Model to get reaction compartments from gene_info
        # this code can be modified in the future to be faster by modifying gene info to take a reaction list rather than a cobra.Model
        # or by using ME_Model class one it is implemented

        expression_machinery_me, expression_module = generate_expression_module(self.me_reactions)
        
        for hgnc_id in tqdm(list(set(expression_machinery_me))):
            expr_reactions, protein_metabolites = get_all_expression_reactions(hgnc_id, machinery_list = expression_machinery_me,
                                                  metabolic_model = expression_module, compress_mrna = self.compress_mrna, 
                                                  ub_args = self.ub_args)


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
                                                      metabolic_model = expression_module, compress_mrna = self.compress_mrna, 
                                                      ub_args = self.ub_args)


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
        if self.dummy_protein:
            print('Express dummy protein')
            dummy_psim = func.average_protein_features(psim_me = self.psim_me, 
                                                      protein_ids = sorted(self.id_protein_map.keys()), 
                                                     context_specific = True)
            dummy_reactions, dm = get_all_expression_reactions(hgnc_id = 'HGNC:DUMMY', psim = dummy_psim, machinery_list = [], 
                                                                metabolic_model = cobra.Model(''), compress_mrna = self.compress_mrna, 
                                                              ub_args = self.ub_args) 
            for r in dummy_reactions:
                for m in r.metabolites:
                    if isinstance(m, Protein) and 'HGNC:DUMMY' in m.id: # str requirement to avoid converting ub proteins
                        m.type = 'dummy_protein'
                        
            self.dummy_protein = {'protein_metabolite': dm[0], 'dummy_expression_reactions': dummy_reactions}
            self.me_reactions += self.dummy_protein['dummy_expression_reactions']
            
        else:
            self.dummy_protein = None
            
    def get_complex_info(self):
        print('Get metabolic module complex information')
        complex_df = pd.DataFrame(columns = ['reaction_id', 'compartment', 'machinery', 'is_complex', 'creates_multiple_reactions'])

        for r in tqdm(self.m_model.reactions):
            compartment_ = func.get_reaction_compartment(r)
            if len(r.genes) == 1: 
                complex_df.loc[complex_df.shape[0], :] = [r.id, compartment_, list(r.genes)[0].id, False, False]
            elif 'and' in r.gene_reaction_rule and 'or' in r.gene_reaction_rule: 
                machinery_final = parse_complex.eval_complex(r.gene_reaction_rule)
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
        print('Get expression module complex information')
        me_complex_df = pd.DataFrame(columns = ['reaction_id', 'compartment', 'machinery', 'is_complex', 'creates_multiple_reactions'])

        for r in tqdm(self.me_reactions):
            compartment_ = func.get_reaction_compartment(r)
            if len(r.genes) == 1: 
                me_complex_df.loc[me_complex_df.shape[0], :] = [r.id, compartment_, list(r.genes)[0].id, False, False]
            elif 'and' in r.gene_reaction_rule and 'or' in r.gene_reaction_rule: 
                machinery_final = parse_complex.eval_complex(r.gene_reaction_rule)
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
        
        # if complexes are duplicated across different reactions assigned to the same compartment, 
        # generate a singular unique id
        dup_complexes = complex_df[complex_df.is_complex].duplicated(subset = ['compartment', 'machinery'], keep = 'first')
        dup_complexes = complex_df.loc[dup_complexes.index[np.where(dup_complexes)]]
        track_dups = dict()
        for i in dup_complexes.index:
            dups = complex_df[(complex_df.compartment == dup_complexes.loc[i,'compartment']) & (complex_df.machinery == dup_complexes.loc[i, 'machinery'])]
            new_id = '_'.join(dups.reaction_id)
            if new_id in track_dups.keys():
                track_dups[new_id] += 1
            else:
                track_dups[new_id] = 0

            new_id = new_id + '_' + str(track_dups[new_id])


            complex_df.loc[dups.index,'complex_id'] = new_id

        # those that didn't need _0
        simplify = [k for k,v in track_dups.items() if v == 0]
        for cid in simplify:
            complex_df.loc[complex_df[complex_df.complex_id == cid + '_0'].index, 'complex_id'] = complex_df[complex_df.complex_id == cid + '_0'].complex_id.apply(lambda x: x[:-2]).tolist()
        
        self.complex_df = complex_df    
        
    def generate_complex_reactions(self):
        # create a mapping of the unique self.complex_df ids to the actual complex metabolite
        unique_complexes = self.complex_df[self.complex_df.is_complex]
        unique_complexes = unique_complexes.drop_duplicates(subset = 'complex_id', keep = 'first')
        unique_complexes.reset_index(inplace = True, drop = True)

        complex_formation_reactions = list() # store all complex formation reactions

        counter = 0
        for i in unique_complexes.index:
            complex_id = unique_complexes.loc[i, 'complex_id']
            compartment = unique_complexes.loc[i, 'compartment']
            machinery = unique_complexes.loc[i, 'machinery'].split(';')

            machinery_metabolites = list()
            for m in machinery:
                if m != 'ribosome':
                    machinery_metabolites.append(self.id_protein_map[m][compartment])
                else:
                    machinery_metabolites.append(self.ribosome_complex_c)

            complex_metabolite = Complex(metabolites = machinery_metabolites, complex_id = complex_id)
            if len(complex_id) > 247: # ids that are too long
                complex_metabolite.update_id(new_id = str(counter)) # complex_metabolite.udate_id()
                counter += 1

                new_id = complex_metabolite.id
                self.complex_df.complex_id.replace(to_replace = complex_id, value = complex_metabolite.temp_id, 
                                                   inplace = True)

            complex_reaction = complex_metabolite.form_complex()

            complex_formation_reactions.append(complex_reaction)
            self.complex_id_metabolite_map[complex_metabolite.temp_id] = complex_metabolite
            self.complex_reactions_map[complex_metabolite.temp_id] = complex_reaction.id

        self.complex_formation_reactions = complex_formation_reactions
        
    def get_keff(self):
        # for reactions that show up more than once
        reactions_to_track = self.m_model.reactions + self.me_reactions
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
        
        if self.dummy_protein is not None:
            self.dummy_protein['keff'] = func.SASA(self.dummy_protein['protein_metabolite'].formula_weight/1000)*(params.keff_median/median_SASA)
    
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
            elif df.machinery.unique().shape[0] == df.shape[0]: 
                # rare case where two different complexes have the same MW
                # instead of randomly choosing, to have consistent results, just choose the first option that appears
                drop_index += df.index.tolist()[1:]
            else:
                raise ValueError('Something went wrong in selecting a complex by lowest molecular weight')

        self.complex_df.drop(index = drop_index, inplace = True)

        # backtrack and remove all protein expression and complex formation reactions of dropped enzymes

        # get rid of redundant complexes
        complexes_to_drop = sorted(set(c_og[c_og.is_complex].complex_id).difference(self.complex_df.complex_id))#sorted(set(self.complex_id_metabolite_map.keys()).difference(self.complex_df.complex_id))
        complexes_to_drop_id = [self.complex_reactions_map[c_id] for c_id in complexes_to_drop]
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
        metabolic_reactions = [r.id for r in self.m_model.reactions]
        final_reactions = []

        # deal with metabolic reactions first
        for i in tqdm(self.complex_df[self.complex_df.category == 'metabolic_reaction'].index):
            reaction_id = self.complex_df.loc[i, 'reaction_id'] # original reaction id
            r_ = self.m_model.reactions.get_by_id(reaction_id).copy() # original reaction

            r = ME_Reaction(type_ = ['catalysis'], 
                            id = r_.id, name = r_.name, subsystem = r_.subsystem, lower_bound = r_.lower_bound, 
                            upper_bound = r_.upper_bound, 
                                cobra_id = r_.id)
            r.add_metabolites(r_.metabolites, combine = False)
            r.gene_reaction_rule = r_.gene_reaction_rule

            metabolites = r.metabolites.copy() # original reaction metabolites

            if not self.complex_df.loc[i, 'is_complex']:
                enzyme_to_couple = self.id_protein_map[self.complex_df.loc[i, 'machinery']][self.complex_df.loc[i, 'compartment']]
                alpha_p = enzyme_to_couple.alpha_p
            else:
                enzyme_to_couple = self.complex_id_metabolite_map[self.complex_df.loc[i, 'complex_id']]
                alpha_p = np.median([p.alpha_p for p in enzyme_to_couple.decompose_complex() if isinstance(p, Protein)])


            # add machinery to substrate side
            c3 = (params.mu + alpha_p)/self.complex_df.loc[i, 'keff']
            if c3.subs(params.mu, 1) <= 0:
                raise ValueError('The catalysis coupling constraint is negative for ' + enzyme_to_couple.id)
            enzyme_to_couple.couple(type = 'catalysis', value = -c3)

            if not r_.reversibility:
                r.couple(metabolites = enzyme_to_couple, types = 'catalysis')
                reactions = [r]
            else: # add a forward and reverse reaction for reversible reactions
                r_f,r_r = r.copy(), r.copy()
                r_f.lower_bound, r_r.lower_bound, r_r.upper_bound = 0,0, abs(r.lower_bound)
                r_r.add_metabolites({metab: -coeff for metab, coeff in r_r.metabolites.items()}, combine = False)
                
                r_f.couple(metabolites = enzyme_to_couple, types = 'catalysis')
                r_r.couple(metabolites = enzyme_to_couple, types = 'catalysis')
                

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

        # dummy protein for orphan reactions (see deorphan)
        if sorted(metabolic_reactions) != sorted([r.id for r in self.m_model.reactions if len(r.genes) == 0]):
            raise ValueError('Not all metabolic reactions that require machinery have been accounted for')
        if self.dummy_protein is None:
            final_reactions += [r.copy() for r in self.m_model.reactions if len(r.genes) == 0]
        self.final_reactions = final_reactions

        
    def add_expression_machinery(self):
        # filter out metabolic reactions
        backup = self.complex_df.copy()
        self.complex_df = self.complex_df[self.complex_df.category == 'expression_reaction']
        self.complex_df.reset_index(inplace = True, drop = True)

        print('Add machinery to expression module reactions')
        for rxn in tqdm([r__ for r__ in self.me_reactions if len(r__.genes) > 0]):
            reaction_id_short = func.parse_me_reaction_id(rxn.id) # abbreviated version
            reaction_id = rxn.id # original reaction id
            idx = self.complex_df[self.complex_df.reaction_id == reaction_id_short].index.tolist()

            if not isinstance(rxn, ME_Reaction):
                rxn_me = ME_Reaction(type_ = ['catalysis'], 
                                id = rxn.id, name = rxn.name, subsystem = rxn.subsystem, lower_bound = rxn.lower_bound, 
                                upper_bound = rxn.upper_bound)
                rxn_me.add_metabolites(rxn.metabolites, combine = False)
                rxn_me.gene_reaction_rule = rxn.gene_reaction_rule
            else: # translation reactions
                rxn_me = ME_Reaction(type_ = rxn.type + ['catalysis'], 
                                id = rxn.id, name = rxn.name, subsystem = rxn.subsystem, lower_bound = rxn.lower_bound, 
                                upper_bound = rxn.upper_bound)
                rxn_me.add_metabolites(rxn.metabolites, combine = False)
                rxn_me.coupled_metabolites = rxn.coupled_metabolites
                rxn_me.gene_reaction_rule = rxn.gene_reaction_rule

            for i in idx:
                # necessary to keep same metabolite type
                r = ME_Reaction(type_ = rxn_me.type,id = rxn_me.id, name = rxn_me.name, 
                                subsystem = rxn_me.subsystem, lower_bound = rxn_me.lower_bound, 
                                upper_bound = rxn_me.upper_bound)
                r.add_metabolites(rxn_me.metabolites, combine = False)
                r.coupled_metabolites = rxn_me.coupled_metabolites
                r.gene_reaction_rule = rxn_me.gene_reaction_rule

                if not self.complex_df.loc[i, 'is_complex']:
                    enzyme_to_couple = self.id_protein_map[self.complex_df.loc[i, 'machinery']][self.complex_df.loc[i, 'compartment']]
                    alpha_p = enzyme_to_couple.alpha_p
                else:
                    enzyme_to_couple = self.complex_id_metabolite_map[self.complex_df.loc[i, 'complex_id']]
                    alpha_p = np.median([p.alpha_p for p in enzyme_to_couple.decompose_complex() if isinstance(p, Protein)])
                # add machinery to substrate side
                c3 = (params.mu + alpha_p)/self.complex_df.loc[i, 'keff']
                if c3.subs(params.mu, 1) <= 0:
                    raise ValueError('The catalysis coupling constraint is negative for ' + enzyme_to_couple.id)
                enzyme_to_couple.couple(type = 'catalysis', value = -c3)

                if not rxn.reversibility:
                    r.couple(metabolites = enzyme_to_couple, types = 'catalysis')
                    reactions = [r]
                else: # add a forward and reverse reaction for reversible reactions
                    r_f,r_r = r.copy(), r.copy()
                    r_f.lower_bound, r_r.lower_bound, r_r.upper_bound = 0,0, abs(r.lower_bound)
                    r_r.add_metabolites({metab: -coeff for metab, coeff in r_r.metabolites.items()}, combine = False)

                    r_f.couple(metabolites = enzyme_to_couple, types = 'catalysis')
                    r_r.couple(metabolites = enzyme_to_couple, types = 'catalysis')
                    
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
        
        if self.dummy_protein is None:
            self.final_reactions += [r for r in self.me_reactions if len(r.genes) == 0]
            self.final_reactions += self.complex_formation_reactions 
            
        self.complex_df = backup.copy()
        del backup
    
    def deorphan(self, exclude = None):
        '''Couples dummy protein to reactions that don't have specified genes ("de-orphaning")

        Parameters
        ----------
        exclude: list, default None
            a list of M_model reaction ids to exclude from coupling to dummy protein (if None, defaults to 
            non-exchange/demand metabolic model reactions and transport reactions; recommended default)

        Returns
        ----------
        deorphaned: list
            a list of ME_Model reaction IDs for reactions that were de-orphaned
        self.orphan: list
            a list of ME_Model reaction IDs for reactions there were not de-orphaned despite having 0 specified genes 

        '''

        if self.dummy_protein is not None:
            print('Deorphan enzymeless reactions')
            enzymeless_reactions = [r for r in self.m_model.reactions if len(r.genes) == 0]
            enzymeless_reactions += [r for r in self.me_reactions if len(r.genes) == 0]
            enzymeless_reactions += self.complex_formation_reactions

            if len(set([r.id for r in enzymeless_reactions]).intersection(self.final_reactions)) > 0:
                raise ValueError('Incorrect parsing of reaction lists for dummy protein')
            if exclude is None:
                # metabolic module enzymes to exclude
                self.orphan = [r for r in self.m_model.exchanges + self.m_model.demands if r in enzymeless_reactions]
                orphan = [] 
                for r in self.orphan: # secondary exchange reactions
                    if len(r.metabolites) > 1 or list(r.metabolites)[0].compartment != 'b':
                        raise ValueError('Incorrectly formatted exchange reaction: ' + r.id + '. Must follow Recon2.2 format.')

                    assoc_rxn = list(list(r.metabolites)[0].reactions)
                    assoc_rxn.remove(r)
                    if len(assoc_rxn) > 0:
                        for r_ in assoc_rxn: # id the second exchange reaction (Recon2.2 format)
                            cond1 = (sorted(r_.compartments) == ['b', 'e'])
                            cond2 = (len(set(['_'.join(m.id.split('_')[:-1]) for m in list(r_.metabolites)])) == 1)
                            cond3 = (len(r_.genes) == 0)
                            if cond1 and cond2 and cond3 :
                                orphan.append(r_)
                self.orphan += orphan
                del orphan
            else:
                m_ids = [r.id for r in self.m_model.reactions]
                for r_id in exclude:
                    if r_id not in m_ids:
                        raise ValueError('The list of metabolic reactions to exclude from dummy catalysis must be in the metabolic model reaction list')
                    if len(self.m_model.reactions.get_by_id(r_id).genes)>0:
                        raise ValueError('The list of metabolic reactions to exclude from dummy catalysis must not have an associated GPR')

                self.orphan = exclude 
            
            # expression module enzymes to exclude
            expression_rids = ['CYTOSOLIC_PROTEIN_FOLDING', 'IMPORTtn', 
                  'RIBOSOME_COMPLEX_DISSOCIATIONc', 'UNFOLDr',
                  'POLYUBIQUITIN_MOIETY_EXPORTtn', 'COMPLEX_FORMATION'] 
            for r in enzymeless_reactions:
                for expr_rid in expression_rids:
                    if r.id.__contains__(expr_rid):
                        self.orphan.append(r)
                        break

            # do not deorphan transport reactions for small molecules (can passively diffuse)
            transport = [r for r in enzymeless_reactions if len(r.compartments) > 1 and r not in self.orphan]

            remove_idx = list()
            for i in range(len(transport)):
                r = transport[i]
                tm = dict()
                mc = dict()
                counter = 0
                for m in r.metabolites:
                    if m.formula_weight <= params.membrane_diffusion_limit: # all metabolites within diffusion limit
                        counter += 1
                    m_id = '_'.join(m.id.split('_')[:-1]) # atleast one metabolite is transported across compartments
                    if m_id not in tm:
                        tm[m_id] = 1
                    else:
                        tm[m_id] += 1
                        mc[m_id] = m.charge

                uncharged = True
                for m_id,v in tm.items():
                    if v >= 2 and mc[m_id] != 0:
                        uncharged = False
                        break
                # uncharged, all metabolites that are transported are < 504 Da, and atleast one metabolite is transported 
                if max(list(tm.values())) >= 2 and counter == len(r.metabolites) and uncharged:
                    self.orphan.append(r)            

            deorphan = [r for r in enzymeless_reactions if r not in self.orphan]
            self.deorphaned = list()

            if len(deorphan) > 0:
                c3 = (params.mu + self.dummy_protein['protein_metabolite'].alpha_p)/self.dummy_protein['keff']
                self.dummy_protein['protein_metabolite'].couple(type = 'catalysis', value = -c3)

                metabolic_ids = [r.id for r in self.m_model.reactions]
                for r_ in deorphan:
                    if r_.id in metabolic_ids:
                        cid = r_.id
                    else:
                        cid = None
                    r = ME_Reaction(type_ = ['catalysis'], 
                                    id = r_.id, name = r_.name, subsystem = r_.subsystem, lower_bound = r_.lower_bound, 
                                    upper_bound = r_.upper_bound, 
                                        cobra_id = cid)
                    r.add_metabolites(r_.metabolites, combine = False)
                    r.gene_reaction_rule = r_.gene_reaction_rule
                    metabolites = r.metabolites.copy() # original reaction metabolites

                    if not r_.reversibility:
                        r.couple(metabolites = self.dummy_protein['protein_metabolite'], types = 'catalysis')
                        reactions = [r]
                    else: # add a forward and reverse reaction for reversible reactions
                        r_f,r_r = r.copy(), r.copy()
                        r_f.lower_bound, r_r.lower_bound, r_r.upper_bound = 0,0, abs(r.lower_bound)
                        r_r.add_metabolites({metab: -coeff for metab, coeff in r_r.metabolites.items()}, combine = False)

                        r_f.couple(metabolites = self.dummy_protein['protein_metabolite'], types = 'catalysis')
                        r_r.couple(metabolites = self.dummy_protein['protein_metabolite'], types = 'catalysis')
                        r_f.id, r_r.id = r_f.id + '_F', r_r.id + '_R'
                        reactions = [r_f, r_r]

                    self.deorphaned += [r.id for r in reactions]
                    self.final_reactions += reactions


            self.final_reactions += self.orphan
            self.orphan = [r.id for r in self.orphan + biomass.biomass_reactions + [biomass.upb_reaction]]

    def build_me_model(self, model_id = 'HUMAN_ME_MODEL'):
        print('Add biomass component to reactions')
        for r in self.final_reactions:
            biomass.add_biomass_change(r)

        br = [copy.deepcopy(r) for r in biomass.biomass_reactions]
#         br.append(self.pb_reaction) 
        if self.dummy_protein is not None:
            br.append(copy.deepcopy(biomass.upb_reaction)) 
        
        self.final_reactions += br

        print('Generate ME-Model')
        me_model = ME_Model(m_model = self.m_model, id_or_model = model_id)
        me_model.add_reactions(self.final_reactions)
        me_model.orphan = self.orphan
        me_model.deorphaned = self.deorphaned
        me_model.check()
        
#         del self.pb_reaction
        del self.ub_args
        del self.me_reactions
        del self.final_reactions
        del self.complex_formation_reactions
        del self.m_model
        del self.orphan
        del self.deorphaned

        return me_model


# In[4]:


def build_me(minimal_proteome = False, compress_mrna = False,
             dummy_protein = True, model_id = 'HUMAN_ME_MODEL'):
    '''Generates a human ME_model. 
    
    Parameters
    ----------
    minimal_proteome: bool
        For reactions with OR in the GPR, the builder by default (False) generates a 
        separate reaction for each protein complex (False). If True, builder instead will create one reaction, 
        choosing the protein complex with the lowest molecular weight to catalyze the reaction. If a reaction
        has multiple enzyme options with the same molecular weight, will randomly choose one. 
    compress_mrna: bool
        If true, will merge the 3 linear mrna reactions--transcription, processing, and nuclear export--for each
        gene into a single reaction
    dummy_protein: bool [True]
        whether to add a representative dummy protein to catalyze orphan reactions 
    model_id: string; id for the me model

    
    '''
    start = time.time()
    
    builder = me_builder(compress_mrna = compress_mrna, 
                         dummy_protein = dummy_protein, psim_me = params.psim_me, 
                         m_model = params.human_model)
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
    builder.deorphan()
    me_model = builder.build_me_model(model_id = model_id)

    end = time.time()
    print('Time to build: {} minutes'.format((end-start)/60))


    return me_model, builder


# In[5]:


gc.collect()


# In[6]:


# minimal_proteome = True
# compress_mrna = True 
# dummy_protein = True 
# model_id = 'HUMAN_ME_MODEL'


# In[7]:


# builder = me_builder(compress_mrna = compress_mrna, 
#                      dummy_protein = dummy_protein, psim_me = params.psim_me, 
#                      m_model = params.human_model)
# builder.express_metabolic_enzymes()
# builder.express_expression_enzymes()
# builder.express_dummy_protein()
# builder.get_complex_info()
# builder.generate_complex_reactions()
# builder.get_keff()
# if minimal_proteome:
#     builder.minimize_proteome()

# builder.add_metabolic_machinery()
# builder.add_expression_machinery()
# builder.deorphan()
# me_model = builder.build_me_model(model_id = model_id)
# print('solve')
# sln, stat, _ = me_model.solve_lp(mu_val = 1e-9)


# In[ ]:




