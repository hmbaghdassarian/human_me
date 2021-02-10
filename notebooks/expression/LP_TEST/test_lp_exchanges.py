#!/usr/bin/env python
# coding: utf-8
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
from macromolecules.macromolecule import Macromolecule
from expression.build_me_model import flatten_list
from utils.parameters import human_model as m_model

mu_val = 1e-9
n_cores = 10
counter = 7

lp_path = '/data2/hratch/human_me/other/test_lp/'

def add_boundary(m_id, type, tme_new = None):
    '''m is a cobra.metabolite or metabolite ID'''
    
    if tme_new is None:
        with open(lp_path + 'working_version_' + str(counter) + '.pickle', 'rb') as handle:
            tme_new = pickle.load(handle)

    tme_new.add_boundary(tme_new.metabolites.get_by_id(m_id), type =type)
    sln, stat, _ = tme_new.solve_lp(mu_val = mu_val)
    
    fn = '/data2/hratch/human_me/other/test_lp/test_' + type + '.tab'
    if not os.path.isfile(fn):
        with open(fn, 'a+') as f:
            f.write('metabolite_id' + '\t' + 'status' + '\n')
        
        
    with open(fn, 'a+') as f:
        f.write(m_id + '\t' + str(stat.max()) + '\n')


#######---------------------------------------------------------------------------------------------------
import multiprocessing
import gc
metabs_ = [m.id for m in m_model.metabolites]

# demands
print('Start demands')
pool = multiprocessing.Pool(processes = n_cores)
try:
    res = pool.starmap(add_boundary, zip(metabs_, ['demand']*len(metabs_)))
    pool.close()
    pool.join()
    gc.collect()
except:
    pool.close()
    pool.join()
    gc.collect()
    raise ValueError('Parallelization failed')

#sinks 
# demands
print('Start sinks')
pool = multiprocessing.Pool(processes = n_cores)
try:
    res = pool.starmap(add_boundary, zip(metabs_, ['sink']*len(metabs_)))
    pool.close()
    pool.join()
    gc.collect()
except:
    pool.close()
    pool.join()
    gc.collect()
    raise ValueError('Parallelization failed')                       

