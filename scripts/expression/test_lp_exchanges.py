#!/usr/bin/env python
# coding: utf-8

# In[163]:


import pickle
import sys
import copy
import time

import cobra
import sympy
import pandas as pd
import numpy as np

from tqdm import tqdm

sys.path.insert(1, '/home/hratch/Projects/human_me/scripts/')
from utils import functions as func
from utils import parameters as params


# In[175]:


lp_path = '/data2/hratch/human_me/test_lp/'
with open(lp_path + 'toy_me_model.pickle', 'rb') as handle:
    me_model = pickle.load(handle)


# In[73]:


scale_proxy = 1/1e2
scale_mrna = None
with open(lp_path + 'toy_me_model_2.pickle', 'rb') as handle:
    me_model = pickle.load(handle)

new_reactions=[]

for r in me_model.reactions:
    if ('TRANSLATION_ELONGATION' in r.id or '_co_TRANSLOC_IMPORTtr' in r.id) and r.id != 'TRANSLATION_ELONGATIONc_COMPLEX_FORMATIONc':
        rxn = {k:v for k,v in r.metabolites.items() if '_mrna[c]' not in k.id and 'mrna_deg_proxy' not in k.id}

        to_add = set(r.metabolites.keys()).difference(rxn.keys())
        for m in to_add:
            if scale_mrna is not None:
                if '_mrna[c]' in m.id:
                    rxn[m] = r.metabolites[m]*scale_mrna
            if 'mrna_deg_proxy' in m.id:
                rxn[m] = r.metabolites[m]*scale_proxy

        r_ = func.ME_Reaction(r.id, type_ = r.type)
        r_.lower_bound = r.lower_bound
        r_.upper_bound = r.upper_bound
        r_.gene_reaction_rule = r.gene_reaction_rule
        r_.add_metabolites(rxn)

        new_reactions.append(r_)
    else: 
        new_reactions.append(r.copy())

me_model = func.ME_Model('')
me_model.add_reactions(new_reactions)
sln, status, _ = me_model.solve_lp(mu_val = 1e-9) 


# In[119]:


exchanges = [r for r in me_model.reactions if 'EX_' in r.id and sorted(r.compartments) == ['b', 'e']]
ef = [r for r in full_model.reactions if 'EX_' in r.id and sorted(r.compartments) == ['b', 'e']]


# In[120]:


lb = []
ub = []
for r in ef:
    lb.append(r.lower_bound)
    ub.append(r.upper_bound)
print(set(lb))
print(set(ub))


# In[110]:





# In[121]:


list(set([r.id for r in ef]).intersection([r.id for r in exchanges]))[:10]


# In[122]:


ex_id = 'EX_HC00342_LPAREN_e_RPAREN_'
full_model.reactions.get_by_id(ex_id).bounds


# In[123]:


me_model.reactions.get_by_id(ex_id).bounds


# In[17]:


# res = pd.DataFrame(columns = ['mrna_coef', 'proxy_coef', 'status'])
# counter = 0

# replace_mrna = 0
# for replace_proxy in tqdm([1,10]):

#     with open(lp_path + 'toy_me_model.pickle', 'rb') as handle:
#         me_model = pickle.load(handle)

#     new_reactions=[]

#     for r in me_model.reactions:
#         if ('TRANSLATION_ELONGATION' in r.id or '_co_TRANSLOC_IMPORTtr' in r.id) and r.id != 'TRANSLATION_ELONGATIONc_COMPLEX_FORMATIONc':
#             rxn = {k:v for k,v in r.metabolites.items() if '_mrna[c]' not in k.id and 'mrna_deg_proxy' not in k.id}

#             to_add = set(r.metabolites.keys()).difference(rxn.keys())
#             for m in to_add:
#                 if replace_mrna is not None:
#                     if '_mrna[c]' in m.id:
#                         rxn[m] = -replace_mrna
#                     else:
#                         rxn[m] = -replace_proxy

#             r_ = func.ME_Reaction(r.id, type_ = r.type)
#             r_.lower_bound = r.lower_bound
#             r_.upper_bound = r.upper_bound
#             r_.gene_reaction_rule = r.gene_reaction_rule
#             r_.add_metabolites(rxn)

#             new_reactions.append(r_)
#         else: 
#             new_reactions.append(r.copy())

#     me_model = func.ME_Model('')
#     me_model.add_reactions(new_reactions)
#     sln, status, _ = me_model.solve_lp(mu_val = 1e-9) 
    
#     res.loc[counter, :] = [replace_mrna, replace_proxy, status.max()]
#     counter +=1


# In[131]:


PTR = 875
x = 4e-3*PTR
((0.02*x)-0.069)/(1-x)


# In[132]:


0.02*3.5


# In[128]:


0.08-0.069


# In[129]:


1-4


# In[ ]:





# In[ ]:





# In[59]:


# res = pd.DataFrame(columns = ['mrna_coef', 'status'])

# replace_mrna = 0.004
# replace_proxy = 0#0.60



# with open(lp_path + 'toy_me_model.pickle', 'rb') as handle:
#     me_model = pickle.load(handle)

# new_reactions=[]

# for r in me_model.reactions:
#     if ('TRANSLATION_ELONGATION' in r.id or '_co_TRANSLOC_IMPORTtr' in r.id) and r.id != 'TRANSLATION_ELONGATIONc_COMPLEX_FORMATIONc':
#         rxn = {k:v for k,v in r.metabolites.items() if '_mrna[c]' not in k.id and 'mrna_deg_proxy' not in k.id}

#         to_add = set(r.metabolites.keys()).difference(rxn.keys())
#         for m in to_add:
#             if '_mrna[c]' in m.id:
#                 rxn[m] = -replace_mrna
#             else:
#                 rxn[m] = -replace_proxy

#         r_ = func.ME_Reaction(r.id, type_ = r.type)
#         r_.lower_bound = r.lower_bound
#         r_.upper_bound = r.upper_bound
#         r_.gene_reaction_rule = r.gene_reaction_rule
#         r_.add_metabolites(rxn)

#         new_reactions.append(r_)
#     else: 
#         new_reactions.append(r.copy())

# me_model = func.ME_Model('')
# me_model.add_reactions(new_reactions)
# sln, status, _ = me_model.solve_lp(mu_val = 1e-9) 


# In[73]:


m = [m for m in me_model.metabolites if 'HGNC:DUMMY_unfolded_protein[c]' == m.id][0]


# In[87]:


list(me_model.metabolites.get_by_id('biomass_unmodeled_protein').reactions)[1].reaction


# In[ ]:





# In[ ]:





# In[ ]:





# In[27]:


# lp_path = '/data2/hratch/human_me/test_lp/'

# # #For some reason, this protein complex used to catalyze the formation of the pre40s complex causes infeasiblity. 
# # #Including any of the 2 of the complex subcomponents in the precursors_fail list makes it work. Or their earlier 
# # #versions (up to unfolded_protein_c).

# # error_metabolites = ['pre40s_rrna_protein_COMPLEX_FORMATIONn_protein_complex[n]']
# # precurors_work = ['HGNC:21173_folded_protein[n]', 'HGNC:32790_folded_protein[n]']
# # precursors_fail = ['HGNC:25542_folded_protein[n]', 'HGNC:29100_folded_protein[n]']

# # # folded_protein[n] <-- folded_protein[c] <-- unfolded_protein[c] <-- 
# # # adding unfolded_protein[c] of precursors fail works too, don't need both of the precursors fail, just 1...
# # error_metabolites = precursors_fail.copy()

# def remove_metabolite(test_metabolites = [], mu_val = 0.01):
    
#     with open(lp_path + 'toy_me_model.pickle', 'rb') as handle:
#         me_model = pickle.load(handle)
    
# #     me_metabolites = [m.id for m in me_model.metabolites if 'deg_proxy' in m.id or ('mrna[n]' in m.id and 'premrna' not in m.id and 'lariats' not in m.id)]
# #     me_metabolites = []
#     me_metabolites += error_metabolites

#     for tm in test_metabolites:
#         me_metabolites.remove(tm)
    
#     ra = []
#     for mm_id in me_metabolites: #me_metabolites:
#         try:
#             mm_obj = me_model.metabolites.get_by_id(mm_id)
#         except:
#             mm_obj = params.human_model.metabolites.get_by_id(mm_id)
#         r = cobra.Reaction('TEST_' + mm_obj.id)
#         r.add_metabolites({mm_obj: 1}, reversibly = True)
#         ra.append(r)
#     if len(ra) > 0:
#         me_model.add_reactions(ra)
#     sln, status, _ = me_model.solve_lp(mu_val = mu_val)
#     return sln, status

# #max 1e-4
# sln,status = remove_metabolite(mu_val = 0)

