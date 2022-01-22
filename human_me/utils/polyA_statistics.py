#!/usr/bin/env python
# coding: utf-8
from typing import SupportsRound, Optional, Union

import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm

from human_me.data.data import build_files_url

# polyA polyA_params
polyA = pd.read_csv(build_files_url + 'polyA_length.csv', index_col=0)
polyA_params = st.johnsonsu.fit(polyA.MEAN)
idx = sorted(set(polyA.SD.dropna().index.tolist()).intersection(polyA.MEAN.dropna().index.tolist()))
reg_data = polyA.loc[idx, ['SD', 'MEAN']]
X = sm.add_constant(reg_data.MEAN, prepend=False)
polyA_mod = sm.OLS(reg_data.SD, X).fit()
min_polyA_mean = -polyA_mod.params['const'] / polyA_mod.params['MEAN']


def calculate_polyA_length(polyA_length: Optional[Union[int, float]] = None, stochastic: bool = False, seed: Optional[int] = None) -> int:
    """Calculates expected length of polyA tail based on input float and data distribution

    Parameters
    ----------
    polyA_length : Union[int, float], optional
        nt length of the polyA tail, by default None 
    stochastic : bool, optional
        whether to estimate the polyA tail length from a regression model based on length input if length is provided
        or draw from a distribution if length is not provided, by default False
    seed : int, optional
        ensures the same polyA length on separate runs if polyA_length is not provided and stochastic is True, by default None

    Returns
    -------
    int
        nt length of the polyA tail
    """
    if not stochastic:
        if polyA_length is None or pd.isna(polyA_length):
            polyA_length = polyA_params[-2]
    else:
        if polyA_length is None or pd.isna(polyA_length):
            np.random.seed(seed)
            polyA_length = st.johnsonsu.rvs(loc=polyA_params[-2], scale=polyA_params[-1], *polyA_params[:-2])
        else:
            if polyA_length > min_polyA_mean:
                polyA_length = polyA_mod.predict((polyA_length, 1))[0]
    return int(round(polyA_length))
