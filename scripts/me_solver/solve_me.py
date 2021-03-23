#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import copy
import warnings
import time

from math import log
import numpy as np
import scipy
from collections import OrderedDict

from qminospy.solver import QMINOS # need solveME installed and working

import sys
sys.path.insert(1, '../../scripts/')
from core.reaction import ME_Reaction
from utils import functions as func


# In[3]:


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
            if not(isinstance(r, ME_Reaction) and r.type == ['biomass']):
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


        qminos_solver = QMINOS()

        sln, stat, hsq = qminos_solver.solvelp(A,b,c,xl,xu,csense,precision=self.precision) 

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


# In[4]:


# import pickle
# import sys
# import copy
# import time

# import cobra

# sys.path.insert(1, '/home/hratch/Projects/human_me/scripts/')
# from utils import functions as func


# In[5]:


# # test LP
# test_model = func.ME_Model(cobra.io.load_json_model('/data2/hratch/Software/qminos_solver/solvemepy/examples/models/iJO1366.json'))
# # test_model.constraints[0].lb = -5
# # test_model.constraints[0].ub = 2

# # test_model.constraints[20].ub = 50
# # test_model.constraints[20].lb = 30
# xq,statq,hsq = test_model.solve_lp(mu_val = 0.03, objective = {'BIOMASS_Ec_iJO1366_core_53p95M': 1}, 
#                                   precision = 'quad')


# In[6]:


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


# In[7]:


# test_model.solve_lp(mu_val = 1e10, objective = {'rA': 1})


# In[ ]:




