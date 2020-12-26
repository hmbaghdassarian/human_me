#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra

import pandas as pd

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

import multiprocessing
from multiprocessing import Pool
from tqdm import tqdm

import scipy.stats as st
import urllib
# import seaborn as sns
# import matplotlib.pyplot as plt
import numpy as np
import warnings
import gc
from tqdm import tqdm

import requests, sys, json, re
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils.load_environmental_variables import *
prebuild = '/data2/hratch/human_me/prebuild/'
# human_model = cobra.io.load_json_model(local_data_path + 'processed/corrected_model.json')


# # PREMRNA Sequence

# In[650]:


from pandarallel import pandarallel
pandarallel.initialize(nb_workers = 20, verbose = 0)

psim_me = pd.read_csv('/data2/hratch/human_me/temp_psim.csv', index_col = 0)

def get_premrna_seq(ensg_id):
    try:
        hyperlink = 'https://rest.ensembl.org/sequence/id/' + ensg_id + '?' 
        x = requests.get(hyperlink, headers={ "Content-Type" : "text/plain"}).text # introns and UTRs
        return str(Seq(x).transcribe())
    except:
        return float('nan')


# In[706]:


fail_mssg = 'You have exceeded uhe limiu of 15 requesus per second; please reduce your concurrenu connecuions'
psim_me['PREMRNA_SEQ'] = fail_mssg
fail_idx = psim_me.index
print('Begin getting sequence for {} genes'.format(psim_me.shape[0]))


fail = True
counter = 1


while fail and (counter < 11):
    print('Iteration: {}'.format(counter))
    
    if counter > 7:
        pandarallel.initialize(nb_workers = 5, verbose = 0)
    elif counter > 4:
        pandarallel.initialize(nb_workers = 10, verbose = 0)

    psim_temp = psim_me.loc[fail_idx, :]
    psim_temp.ENSG_ID = psim_temp.ENSG_ID.apply(lambda x: x.split('.')[0])
    psim_temp['PREMRNA_SEQ'] = psim_temp.ENSG_ID.parallel_apply(lambda x: get_premrna_seq(x))
    psim_me.loc[fail_idx, 'PREMRNA_SEQ'] = psim_temp.PREMRNA_SEQ.tolist()
    
    fail_idx = psim_me[psim_me['PREMRNA_SEQ'] == fail_mssg].index
    fail = (fail_idx.shape[0] > 0)
    print('Missing sequences: {}'.format(fail_idx.shape[0]))
    print(fail)

    counter += 1
    print('---------------')
psim_me.to_csv('/data2/hratch/human_me/temp_psim2.csv')