import sys
sys.path.insert(1, '../../../scripts/')
from preprocess import preprocess

from preprocess import correct_inputs as ci
from utils.load_environmental_variables import *

import pandas as pd
from expression import build_me_model

dummy_protein = [True, False]
minimal_proteome = [True, False]
compress_mrna = [True, False]

counter = 0
res = pd.DataFrame(columns = ['dp', 'mp', 'cm', 'status'])
for dp in dummy_protein:
    for mp in minimal_proteome:
        for cm in compress_mrna:
            try:
                toy_me_model, builder = build_me_model.build_me(minimal_proteome = mp, compress_mrna = cm, 
                                                            dummy_protein = dp)
                sln, stat, _ = toy_me_model.solve_lp(mu_val =  1e-9)
                res.loc[counter, : ] = [dp, mp, cm, stat.max()]
            except:
                res.loc[counter, : ] = [dp, mp, cm, float('nan')]
            counter += 1
res.to_csv('trash2.csv')