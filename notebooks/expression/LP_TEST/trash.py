#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle
import sys
import copy
import time
import os

import cobra
import sympy
import pandas as pd
import numpy as np

import multiprocessing
import gc

from tqdm import tqdm

sys.path.insert(1, '/home/hratch/Projects/human_me/scripts/')
from utils import functions as func
from utils import parameters as params
from macromolecules.macromolecule import Macromolecule
from expression.build_me_model import flatten_list
from utils.parameters import human_model as m_model
from core.model import load_pickled_model


# In[2]:


mu_val = 1e-9
n_cores = 5
counter = 0

lp_path = '/data2/hratch/human_me/other/test_lp/'
tme0 = load_pickled_model(lp_path + 'working_version_' + str(counter) + '.pickle')

def par_sinks(m_id):
    '''Add m_id and solve'''
    
    tme1 = tme0.copy()
    
    model_reaction_ids = [r.id for r in tme1.reactions]
    
    sinks = list()
    r = cobra.Reaction('SK_' + m_id)
    r.add_metabolites({tme1.metabolites.get_by_id(m_id): -1})
    r._lower_bound = -1000
    r._upper_bound = 1000
    if r.id not in model_reaction_ids:
        sinks.append(r)

    tme1.add_reactions(sinks)
    sln, stat, _ = tme1.solve_lp(mu_val = mu_val)
    
    fn = '/data2/hratch/human_me/other/test_lp/test_sinks.tab'
    if not os.path.isfile(fn):
        with open(fn, 'a+') as f:
            f.write('metabolite_id' + '\t' + 'status' + '\n')
        
        
    with open(fn, 'a+') as f:
        f.write(m_id + '\t' + str(stat.max()) + '\n')

m_ids = list(set([m.id for m in tme0.metabolites if (not hasattr(m, 'type')) and ('biomass' not in m.id)]))
 

import multiprocessing
import gc
n_cores = 18

print('Start parallelization')
pool = multiprocessing.Pool(processes = n_cores)
try:
    res = pool.map(par_sinks, m_ids)
    pool.close()
    pool.join()
    gc.collect()
except:
    pool.close()
    pool.join()
    gc.collect()
    raise ValueError('Parallelization failed')                       