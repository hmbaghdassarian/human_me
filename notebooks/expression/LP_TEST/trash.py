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

import pickle
from utils.parameters import human_model as m_model
mu_val = 1e-9

lp_path = '/data2/hratch/human_me/other/test_lp/'

test_metab = list()
for r in m_model.exchanges:
    m = list(r.metabolites.keys())[0].copy()
    m.compartment = 'c'
    m.id = m.id.replace('_b', '_c')
    test_metab.append(m)

def add_sink(m):
    with open(lp_path + 'dummy_protein_infeasible.pickle', 'rb') as handle:
        tme3 = pickle.load(handle)
    
    tme3.add_boundary(m, type ='sink')
    sln3, stat, _ = tme3.solve_lp(mu_val = mu_val)
    if stat.max() == 0:
        with open('trash.txt', 'a+') as f:
            f.write(m.id + '\n')

import multiprocessing
import gc

n_cores = 25
pool = multiprocessing.Pool(processes = n_cores)
try:
    res = pool.map(add_sink, test_metab)
    pool.close()
    pool.join()
    gc.collect()
except:
    pool.close()
    pool.join()
    gc.collect()
    raise ValueError('Parallelization failed')
