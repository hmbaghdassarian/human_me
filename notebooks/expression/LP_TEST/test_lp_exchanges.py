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

from tqdm import tqdm

sys.path.insert(1, '/home/hratch/Projects/human_me/scripts/')
from utils import functions as func
from utils import parameters as params
from utils.parameters import human_model as m_model
from expression.build_me_model import flatten_list


mu_val = 1e-9
n_cores = 2
base = 0
counter = 3

lp_path = '/data2/hratch/human_me/other/test_lp/'

test_metab = list()
for r in m_model.exchanges:
    m = list(r.metabolites.keys())[0].copy()
    m.compartment = 'c'
    m.id = m.id.replace('_b', '_c')
    test_metab.append(m)

def add_sink(m, tme_new = None):
    '''m is a cobra.metabolite or metabolite ID'''
    
    if tme_new is None:
        with open(lp_path + 'working_version_' + str(counter) + '.pickle', 'rb') as handle:
            tme_new = pickle.load(handle)
    
    if isinstance(m, cobra.Metabolite): # object
        m_id = m.id
    else: # string
        m_id = m
    
    tme_new.add_boundary(tme_new.metabolites.get_by_id(m_id), type ='sink')
    sln, stat, _ = tme_new.solve_lp(mu_val = mu_val)
    
    if not os.path.isfile('test_sinks.tab'):
        with open('test_sinks.tab', 'a+') as f:
            f.write('metabolite_id' + '\t' + 'status' + '\n')
        
        
    with open('test_sinks.tab', 'a+') as f:
        f.write(m_id + '\t' + str(stat.max()) + '\n')
        

with open(lp_path + 'working_version_' + str(counter) + '.pickle', 'rb') as handle:
    tme_new = pickle.load(handle)

metabs_ = ['lys_L_e', 'gln_L_e', 'glu_L_l', 'asp_L_c', 'glu_L_c', 'lys_L_c', 'arg_L_c', 'asp_L_l', 'adp_c', 
           'pi_c', 'glc_D_e', 'glc_D_c', 'adp_l', 'gsn_l', 'arg_L_l', 'atp_c', 'atp_l', 'gln_L_c', 'pi_l', 
           'gsn_c'] 
adp_c = tme_new.metabolites.get_by_id('adp_c')
test = flatten_list([[m.id for m in r.products] for r in tme_new.reactions if adp_c in r.reactants])
metabs_ = sorted(set(test).difference(metabs_))

dadp_c = tme_new.metabolites.get_by_id('dadp_c')
test = flatten_list([[m.id for m in r.products] for r in tme_new.reactions if adp_c in r.reactants])
test = sorted(set(test).difference(metabs_))

import multiprocessing
import gc

pool = multiprocessing.Pool(processes = n_cores)
try:
    res = pool.map(add_sink, test)
    pool.close()
    pool.join()
    gc.collect()
except:
    pool.close()
    pool.join()
    gc.collect()
    raise ValueError('Parallelization failed')