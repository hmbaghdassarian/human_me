#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import copy
import warnings
import time
from tqdm import tqdm

from pathos.multiprocessing import ProcessingPool as Pool
import gc

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')

import pandas as pd
from math import log
import numpy as np
import scipy
from collections import OrderedDict

from qminospy.solver import QMINOS # need solveME installed and working

import sys
sys.path.insert(1, '../../scripts/')
from core.reaction import Biomass_Reaction
from utils import functions as func


# In[ ]:


class qminos_solver():
    def __init__(self, precision = 'quad'):
        '''Initializer for solving with qMINOS
        
        Parameters
        ----------  
        precision: string, default "quad"
            The precision for the qminos solver (options ['double', 'quad', 'dq', 'dqq'])
        
        '''
        
        self.precision = 'quad'
    
    def solve_lp(self, me_model, mu_val, objective = {'biomass_dilution': 1}, tolerance = 0,
                 close_biomass_dilution = True):   
        '''Solves the linear program for a specified objective at a specified growth rate

        Parameters
        ----------
        me_model: human_me.core.model.ME_Model
            ME Model to solve
        mu_val: float
            The growth value for which to solve the linear program
        objective: dict, default {'bimoass_dilution': 1}
            The objective function to optimize. Dictionary represent a linear combination of reactions to optimize,
            with reaction ids as keys and the coefficient of the lin. comb. as the values. 
            Values must either be 1 for maximization or -1 for minimization.

            Example: simplest case, to maximize reaction with id 'A', objective = {'A': 1}
        tolerance: float; default 0
            Threshold below which expected sensitivity of solver is too low to detect infeasibility
        close_biomass_dilution: bool, default True
            Internal use only, whether to constrain the biomass_dilution reaction bounds by mu

        Returns
        ----------
        sln: 1D np.array
            the vector of fluxes in the optimal solution
        stat: 
            the solver status 
                0     Optimal solution found.
                1     The problem is infeasible.
                2     The problem is unbounded (or badly scaled).
                3     Too many iterations.
                4     Apparent stall.  The solution has not changed
                      for a large number of iterations (e.g. 1000).
        hsq: 
            optimal basis (see qminospy.solver.QMINOS)
        '''

        # normalize objective to be 1
        if not all([np.sign(v) == 1 for v in objective.values()]) or all([np.sign(v) == -1 for v in objective.values()]):
            raise ValueError('Current version can only maximize or minimize combinations of objectives, not do both')
        tot = abs(sum(objective.values()))
        objective = {k: v/tot for k,v in objective.items()}
        
        # get stoichiometric matrix at mu_val
        S = me_model.create_stoichiometric_matrix(mu_val = mu_val, array_type = 'numpy', inplace = False)

        # get equality and inequality matrices to format (Ev=b; lb <= Iv <= ub; A = concat[E,I])
        zero_tol=1e-6
        inequality_index = []
        equality_index = []
        b_equal = []
        b_less = []
        b_greater = []

        counter = 0
        for metab in me_model.metabolites:
            lb = -np.inf if metab.constraint.lb is None else metab.constraint.lb
            ub = np.inf if metab.constraint.ub is None else metab.constraint.ub
            equality = (ub - lb) < zero_tol
            if equality:
                b_equal.append(lb if abs(lb) > zero_tol else 0.0)
                equality_index.append(counter)
            else:
                b_less.append(ub)
                b_greater.append(lb)
                inequality_index.append(counter)
            counter += 1

        E = S[equality_index,]
        I = S[inequality_index,]

        A = scipy.sparse.dok_matrix(np.concatenate([E,I,I])) # E, L, G
        b = np.array(b_equal + b_less + b_greater)
        csense = np.array(['E']*len(b_equal) + ['L']*len(b_less) + ['G']*len(b_greater))
        del b_less
        del b_greater
        del b_equal
        del S
        
        if tolerance is None:
            tolerance = 0
        if tolerance < 0:
            tolerance = abs(tolerance)

        # reaction bounds at mu
        xl, xu = np.zeros(len(me_model.reactions)), np.zeros(len(me_model.reactions))
        xl[:], xu[:] = np.nan, np.nan
        counter = 0
        for r in me_model.reactions:
            if not(isinstance(r, Biomass_Reaction)):
                xl[counter] = copy.copy(r.lower_bound) - tolerance
                xu[counter] = copy.copy(r.upper_bound) + tolerance
            else:
                if (r.id != 'biomass_dilution') or close_biomass_dilution:
                    bounds = r.replace_bound_mu(mu_val = mu_val, inplace = False)
                    xl[counter] = bounds[0] - tolerance
                    xu[counter] = bounds[1] + tolerance
                else:
                    xl[counter] = 0 - tolerance
                    xu[counter] = 1000 + tolerance
            counter +=1
        
        # objective vector (max c.T*v)
        c = np.zeros(len(me_model.reactions))
        for r_id, coeff in objective.items():
            try:
                r_index = me_model.reactions.index(r_id)
            except:
                raise ValueError('Specified reaction id(s) not in model')
            c[r_index] = coeff


        qs = QMINOS()

        sln, stat, hsq = qs.solvelp(A,b,c,xl,xu,csense,precision=self.precision) 

        # remove unwanted output files
        abspath = os.path.abspath(os.getcwd())
        for fn in [os.path.join(abspath, 'fort.9'), os.path.join(abspath, 'fort.11')]:
            if os.path.isfile(fn):
                os.remove(fn)
                
        return sln, stat, hsq   
    def maximize_growth(self, me_model, min_mu=0, max_mu=0.05, mu_accuracy=1e-10, increment = 0.02,
                        tolerance = 0, verbose=True):
        '''Binary search to find the maximum feasible growth rate

        Parameters
        ----------
        me_model: human_me.core.model.ME_Model
            ME Model to solve
        min_mu: float, default 0
            Expected minimum feasible growth rate (~0)
        max_mu: float, default 0.05
            Expected minimum infeasible growth rate (i.e., just above expected maximum feasible growth rate)
        mu_accuracy: float, default 1e-4
            The maximum error in mu after the binary search
        increment: float, default 1
            The amount to increase growth by when searching for maximum infeasible growth rate from max_mu
        tolerance: float; default 0
            Threshold below which expected sensitivity of solver is too low to detect infeasibility
        verbose: bool, default True
            Prints information about each linear program iteration

        Returns
        ----------
        mu_max: int
            the maximum feasible growth value (in hours)
        res: dict
            keys are all attempted growth values, values are dictionaries with keys as output from self.solve_lp
        '''

        objective = {'biomass_dilution': 1} # maximizing for growth
        feasible_mu = []
        infeasible_mu = []
        res = dict()


        def try_mu(mu_val):
            with func.HiddenPrints():
                sln,stat,hsq = self.solve_lp(me_model, mu_val, objective = objective)
            if stat.max() == 1 and mu_val < 1e-9:
                warnings.warn('Model is infeasible at mu = 0. Trying mu = 1e-9 instead')
                mu_val = 1e-9
                with func.HiddenPrints():
                    sln,stat,hsq = self.solve_lp(me_model, mu_val, objective = objective)
                if stat.max() == 1:
                    raise ValueError('Provided minimum mu is infeasible')
                    
            res[mu_val] = {'solution': sln, 'status': stat.max(), 'basis': hsq}

            if stat.max() == 0:#"optimal":
                if verbose:
                    print('The problem has an optimal solution at mu = {} (hrs)'.format(mu_val))
                feasible_mu.append(mu_val)
                return True, sln, stat, hsq 
            elif stat.max() == 1:
                infeasible_mu.append(mu_val)
                if verbose:
                    print('The problem is infeasible at mu = {} (hrs)'.format(mu_val))
                return False, None, None, None
            else:
                raise valueError('The problem returned with stat: {}'.format(stat.max()))

        start = time.time()
        
        
        if verbose:
            print('Trying mu: {}'.format(min_mu))
        bool_, sln,stat, hsq = try_mu(min_mu) # start with minimal
        while try_mu(max_mu)[0]:  # If max_mu was feasible, keep increasing
            max_mu += increment
            if verbose:
                print('Trying mu: {}'.format(max_mu))
        while (infeasible_mu[-1] - feasible_mu[-1]) > mu_accuracy:
            if verbose:
                print('Trying mu: {}'.format((infeasible_mu[-1] - feasible_mu[-1])*0.5))
            bool_, sln,stat, hsq = try_mu((infeasible_mu[-1] + feasible_mu[-1]) * 0.5)

        if verbose:
            tot = ((time.time() - start)/3600)
            print("completed in {:.2f} hours and {} iterations".format(tot, len(feasible_mu+ infeasible_mu)))
      
        mu_max = np.max(feasible_mu)
        res_ = OrderedDict({k: res[k] for k in sorted(list(res.keys()))})
        return mu_max, res_
    def optimize(self, me_model, objective, mu_max, n_points = 10, 
                 tolerance = 0, n_cores = None, graph = True, fig_name = None):
        '''General optimization of any non-growth objective
        
        Parameters
        ----------
        me_model: human_me.core.model.ME_Model
            ME Model to solve
        objective: dict
            The objective function to optimize. Dictionary represent a linear combination of reactions to optimize,
            with reaction ids as keys and the coefficient of the lin. comb. as the values. 
            Values can only be all 1 for maximization, or all -1 for minimization. 

            Example: simplest case, to maximize reaction with id 'A', objective = {'A': 1}
        mu_max: float
            the maximum growth value at which the model is feasible; use .maximize_growth() method to identify (should be <= mu_max output of .maxmimize_growth() method)
        tolerance: float; default 0
            Threshold below which expected sensitivity of solver is too low to detect infeasibility
        n_cores: int, default None
            the number of workers to use for parallelization
        graph: bool; default True
            plot the relationship between growth and the objective function of interest
        fig_name: str; default None
            save the plotted figure to 
            
        Returns
        ----------
        sln: tuple 
            first element is the growth value at which the non-growth objective is optimized
            second element is the optimized non-growth objective value
        predicted: pd.DataFrame
            1000 growth values between 0 and mu_max, with corresponding interpolated objective values
        interp_fit: output of scipy.interpolate.interp1d
            a function to interpolate objective values from growth values, used to generated predicted
        optimal_vals: collections.OrderedDict
            keys are n_points growth values between 0 and mu_max, values are the objective value optimized at 
            the corresponding growth value
        res: dict
            keys are n_points growth values between 0 and mu_max, values are the output of .solve_lp at 
            corresponding growth values with the objective set to the non-growth objective input
        '''
        
        obj_keys = list(objective.keys())
        if len(obj_keys) == 1 and obj_keys[0] == 'biomass_dilution':
            raise ValueError('To optimize for growth, use the .maximize_growth() method')
        elif len(set(list(objective.values()))) != 1:
            raise ValueError('Currently, NLP objectives must be either only maximization or minimization (only all 1s or only all -1s in objective dictionary values)')
             
        
        growth_vals = np.arange(0,mu_max + mu_max/n_points, mu_max/(n_points-1))
        if (n_cores <= 1) or (n_cores is None):
            res = list()
            for mu_val in tqdm(growth_vals):
                sln, stat, hsq = self.solve_lp(me_model = me_model, mu_val = mu_val, objective = objective, tolerance = tolerance,
                             close_biomass_dilution = True)
                res.append([sln, stat, hsq])
        else:
            # msg: Currently, parallelization errors out at /data2/hratch/Software/qminos_solver/solvemepy/qminospy/solver.py
            # at line 241-243. for some reason, parallelization doesn't recognize the self.precision = 'quad' as == 'quad'. 
            # This does not occur in the serial loop. hard-coding the commented lines 244-250 (same as line 193-199)
            # fixes this issue, but is only a temporary solution. 
#             raise ValueError('me_solver/solve_me .optimize() method, see message above for error')
            n_cores = min([n_cores, n_points])
    #         args_ = zip([me_model]*n_points, list(growth_vals), [objective]*n_points, [tolerance]*n_points, 
    #                    [True]*n_points)
            pool = Pool(n_cores)
            try:
                res = pool.map(self.solve_lp, [me_model]*n_points, list(growth_vals), [objective]*n_points, [tolerance]*n_points, 
                       [True]*n_points)
                pool.close()
                pool.join()
                pool.restart()
                gc.collect()
            except:
                pool.close()
                pool.join()
                pool.restart()
                gc.collect()
                raise ValueError('Parallelization failed')
            res = [list(r) for r in res]

        reaction_indeces = [me_model.reactions.index(j) for j in sorted(objective.keys())]
        res = OrderedDict({i[0]: dict(zip(['val', 'sln', 'stat', 'hsq'], 
                    [(dict(zip(sorted(objective.keys()), i[1][0][reaction_indeces])))] + i[1])) \
             for i in zip(growth_vals, res)})
        optimal_vals = OrderedDict({k: sum(v['val'].values()) for k,v in res.items()})

        # estimate objective values across growth using interpolation
        interp_fit = scipy.interpolate.interp1d(x = list(optimal_vals.keys()),y = list(optimal_vals.values()),
                                 bounds_error=False)
        obj_label = 'predicted_' + '_'.join(list(objective.keys()))
        if len(obj_label) > 30:
            obj_label = 'predicted_Objective'
        predicted = pd.DataFrame(data = {'growth': np.arange(0,mu_max + mu_max/1000, mu_max/(1000-1))})
        predicted[obj_label] = predicted.growth.apply(lambda x: interp_fit(x).item()).values.tolist()
        
        if list(objective.values())[0] == 1:
            optimal_val = predicted[obj_label].max()
        elif list(objective.values())[0] == -1:
            optimal_val = predicted[obj_label].min()
          
        optimal_val_growth = predicted[predicted[obj_label] == optimal_val].growth.values.tolist()[0]
        sln = (optimal_val_growth, optimal_val)

        if graph:
            fig, ax = plt.subplots(figsize = (5,5))
            sns.lineplot(x = 'growth', y = obj_label, data = predicted, ax = ax)
            sns.scatterplot(x = list(optimal_vals.keys()), y = list(optimal_vals.values()), color = 'black', ax = ax)
            plt.plot([optimal_val_growth], [optimal_val], marker='o', markersize=3, color="red")
            ax.legend(handles=[mpatches.Patch(color='black', label='Solved'), 
                              mpatches.Patch(color=sns.color_palette('tab10')[0], label='Interpolated')], 
                     fancybox = True, fontsize = 12)
            ax.set_xlabel(r'$\mu$ $[hr^{-1}]$', fontsize = 15, labelpad = 5)
            ax.set_ylabel(ax.get_ylabel().split('predicted_')[1] + r'$\;\frac{mmol}{gDw \;hr}$', fontsize = 15, labelpad = 5)
            ax.tick_params(axis='both', labelsize=12)

            ax.vlines(x = optimal_val_growth, ymin = ax.get_ylim()[0], ymax = optimal_val, linestyles = '--', 
                  color = sns.color_palette('pastel')[3])
            ax.hlines(y = optimal_val, xmin = ax.get_xlim()[0], xmax = optimal_val_growth, linestyles = '--', 
                      color = sns.color_palette('pastel')[3])

            fig.tight_layout();
            if fig_name is not None:
                fig_name = os.path.splitext(fig_name)[0]
                for ext in ['.png', '.pdf', '.svg']:
                    plt.savefig(fig_name + ext, bbox_inches = 'tight')

        return sln, predicted, interp_fit, optimal_vals, res


# In[ ]:


# counter = 0
# lp_path = '/data2/hratch/human_me/other/test_lp/'
# import pickle
# import time

# with open(lp_path + 'working_version_' + str(counter) + '.pickle', 'rb') as handle:
#     me_model = pickle.load(handle)

# # inputs
# mu_val = 1e-9
# mu_max = 0.03
# n_points = 10
# objective = {'ATPS4m_0': 1, 'DTMPK': 1}
# tolerance = 1e-20
# n_cores = 10
# fig_name = None #path/to/fig_nam (no extension)


# In[ ]:


# solver = qminos_solver()
# sln, stat, _ = solver.solve_lp(me_model, mu_val = mu_val, objective = objective, 
#                                  tolerance = tolerance)

# solver = qminos_solver()
# sln, predicted, interp_fit, optimal_vals, res = solver.optimize(me_model = me_model, objective = objective, 
#                                                                 mu_max = mu_max, 
#                n_points = n_points, tolerance = tolerance, n_cores = n_cores)


# In[ ]:


# import pickle
# import sys
# import copy
# import time

# import cobra

# sys.path.insert(1, '/home/hratch/Projects/human_me/scripts/')
# from utils import functions as func


# In[ ]:


# # test LP
# test_model = func.ME_Model(cobra.io.load_json_model('/data2/hratch/Software/qminos_solver/solvemepy/examples/models/iJO1366.json'))
# # test_model.constraints[0].lb = -5
# # test_model.constraints[0].ub = 2

# # test_model.constraints[20].ub = 50
# # test_model.constraints[20].lb = 30
# xq,statq,hsq = test_model.solve_lp(mu_val = 0.03, objective = {'BIOMASS_Ec_iJO1366_core_53p95M': 1}, 
#                                   precision = 'quad')


# In[ ]:


# import cobra
# from utils.functions import ME_Model
# from utils import parameters as params

# test_reaction = ME_Reaction('test', type_ = ['catalysis'])


# A, B, C = cobra.Metabolite('mA'), cobra.Metabolite('mB'), cobra.Metabolite('mC')
# rxn_A, rxn_B = ME_Reaction('rA', type_ = ['translation']), ME_Reaction('rB', type_ = ['translation'])
# rxn_C = ME_Reaction('rC', type_ = ['biomass'])
# rxn_A.add_metabolites({A: -2*params.mu, B: 3*params.mu, C: 2})
# rxn_B.add_metabolites({C: -5*params.mu, B: 10})
# rxn_C._lower_bound, rxn_C._upper_bound = params.mu, params.mu

# test_model = ME_Model('test')
# test_model.add_reactions([rxn_A, rxn_B, rxn_C])

# xq,statq,hsq = test_model.solve_lp(mu_val = 0.03, objective = {'rA': 1})


# In[ ]:


# test_model.solve_lp(mu_val = 1e10, objective = {'rA': 1})


# In[18]:


objective = {1: -2, 3: -4}


# In[19]:





# In[ ]:




