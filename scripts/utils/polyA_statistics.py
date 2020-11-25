#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import statsmodels.api as sm
import scipy.stats as st
import pandas as pd

import requests, sys, json, re
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils.load_environmental_variables import build_files_path


# In[1]:


# polyA polyA_params
polyA = pd.read_csv(build_files_path + 'polyA_length.csv', index_col = 0)
polyA_params = st.johnsonsu.fit(polyA.MEAN)
idx = sorted(set(polyA.SD.dropna().index.tolist()).intersection(polyA.MEAN.dropna().index.tolist()))
reg_data = polyA.loc[idx, ['SD', 'MEAN']]
X = sm.add_constant(reg_data.MEAN, prepend = False)
polyA_mod = sm.OLS(reg_data.SD, X).fit()
min_polyA_mean = -polyA_mod.params['const']/polyA_mod.params['MEAN']


def calculate_polyA_length(polyA_length = None):
    '''Calculates expected length of polyA tail based on input float and data distribution'''

    if polyA_length == None or pd.isna(polyA_length):
         polyA_length = round(st.johnsonsu.rvs(loc=polyA_params[-2], scale=polyA_params[-1], *polyA_params[:-2]))
    else:
        if  polyA_length > min_polyA_mean:
             polyA_length = round(polyA_mod.predict((polyA_length, 1))[0])
        else:
             polyA_length = round(polyA_length)
    return polyA_length

