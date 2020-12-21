#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import scipy
from qminospy.solver import QMINOS # need solveME installed and working

import sys
sys.path.insert(1, '../../scripts/')
from core.reaction import ME_Reaction


# In[8]:


def solve_lp(me_model, mu_val, objective = {'biomass_dilution': 1}, close_biomass_dilution = True,
             solver_type = 'qminos', precision = 'quad'):
        
        '''
        
        mu_val is the growth value at which to optimize. 
        objective is a dictionary with keys as reaction ids to maximize as some linear combination and values as the coefficient for the linear objective
        close_biomass_dilution is a boolean indicating whether to bound biomass_dilution by mu (True) or by [0,1000] (False)
        solver_type is a string, options of [qminos] - must have solveME and qMINOS installed
        precision options for solver_type as in solveME
        
        Returns same outputs as qminospy.solver.solvelp:
        x: optimal solution
        stat: status
        hs: optimal basis
        
        
        stat:
        0     Optimal solution found.
        1     The problem is infeasible.
        2     The problem is unbounded (or badly scaled).
        3     Too many iterations.
        4     Apparent stall.  The solution has not changed
              for a large number of iterations (e.g. 1000).
        
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
        
        # reaction bounds at mu
        xl, xu = np.zeros(len(me_model.reactions)), np.zeros(len(me_model.reactions))
        xl[:], xu[:] = np.nan, np.nan
        counter = 0
        for r in me_model.reactions:
            if not(isinstance(r, ME_Reaction) and r.type == ['biomass']):
                xl[counter] = r.lower_bound
                xu[counter] = r.upper_bound
            else:
                if (r.id != 'biomass_dilution') or close_biomass_dilution:
                    bounds = r.replace_bound_mu(mu_val = mu_val, inplace = False)
                    xl[counter] = bounds[0]
                    xu[counter] = bounds[1]
                else:
                    xl[counter] = 0
                    xu[counter] = 1000
            counter +=1
        
        # objective vector (max c.T*v)
        c = np.zeros(len(me_model.reactions))
        for r_id, coeff in objective.items():
            try:
                r_index = me_model.reactions.index(r_id)
            except:
                raise ValueError('Specified reaction id(s) not in model')
            c[r_index] = coeff

        if solver_type == 'qminos':
            qminos_solver = QMINOS()
            
            sln, stat, hsq = qminos_solver.solvelp(A,b,c,xl,xu,csense,precision=precision) 
            
            # remove unwanted output files
            abspath = os.path.abspath(os.getcwd())
            for fn in [os.path.join(abspath, 'fort.9'), os.path.join(abspath, 'fort.11')]:
                if os.path.isfile(fn):
                    os.remove(fn)
            return sln, stat, hsq   
        else:
            raise ValueError('Only qminos solver is implemented for now')


# In[9]:


# adapted from: https://github.com/SBRG/cobrame/blob/master/cobrame/solve/algorithms.py
from math import log
# from tempfile import mkdtemp
# from os.path import join
import warnings
import time

def binary_search(me_model, min_mu=0, max_mu=0.05, mu_accuracy=1e-4, increment = 0.02,
                  solver_type='qminos', precision = 'quad', objective = {'biomass_dilution': 1}, # solver args
                  verbose=True):
    """Computes maximum feasible growth rate (mu) through a binary search
    The objective function of the model should be set to a dummy
    reaction which forces translation of a dummy protein.
    :param float max_mu: A guess for a growth rate which will be infeasible
    :param float min_mu: A guess for a growth rate which will be feasible
    :param float mu_accuracy: The final error in mu after the binary search
    :param boolean verbose: will print out each mu in the binary search
    """
    if solve not in ['qminos']:
        raise ValueError('Only qMINOS solver is available')
    
    feasible_mu = []
    infeasible_mu = []

    # String formatting for display
    str_places = int(abs(round(log(mu_accuracy)/log(10)))) + 1
    num_format = "%." + str(str_places) + "f"
    mu_str = "mu".ljust(str_places + 2)

    
    def try_mu(mu_val):
        if mu_val == 0:
            warnings.warn('model is infeasible at mu = 0. Using mu = 1e-9 instead.')
            mu = 1e-9
        
        xq,status,hsq = me_model.solve_lp(mu_val, objective = objective, solver_type = solver_type, precision = precision)
       
        if status.max() == 0:#"optimal":
            if verbose:
                print('The problem has an optimal solution at mu = ' + num_format.format(mu_val) + ' (hrs)')
            feasible_mu.append(mu)
            return True, xq, status, hsq 
        elif status.max() == 1:
            infeasible_mu.append(mu)
            if verbose:
                print('The problem is infeasible at mu = ' + num_format.format(mu_val) + ' (hrs)')
            return False
        else:
            raise valueError('The problem returned with status: {}'.format(status.var()))

    start = time.time()
    # find highest possible value of mu try the edges of binary search
    if not try_mu(min_mu)[0]:
        # Try 0 if min_mu failed
        if min_mu <= 1e-9 or not try_mu(0):
            raise ValueError("0 needs to be feasible")
    while try_mu(max_mu)[0]:  # If max_mu was feasible, keep increasing
        max_mu += increment
    while (infeasible_mu[-1] - feasible_mu[-1]) > mu_accuracy:
        bool_, xq,status, hsq = try_mu((infeasible_mu[-1] + feasible_mu[-1]) * 0.5)

    if verbose:
        tot = ((time.time() - start)/3600)
        print("completed in {:.2f} seconds and {} iterations".format(tot, len(feasible_mu) + len(infeasible_mu)))
              
    return xq, status, hsq, feasible_mu, infeasible_mu

            


# In[ ]:


# import pickle
# import sys
# import copy
# import time

# import cobra

# sys.path.insert(1, '/home/hratch/Projects/human_me/scripts/')
# from utils import functions as func


# In[49]:


# # test LP
# test_model = func.ME_Model(cobra.io.load_json_model('/data2/hratch/Software/qminos_solver/solvemepy/examples/models/iJO1366.json'))
# # test_model.constraints[0].lb = -5
# # test_model.constraints[0].ub = 2

# # test_model.constraints[20].ub = 50
# # test_model.constraints[20].lb = 30
# xq,statq,hsq = test_model.solve_lp(mu_val = 0.03, objective = {'BIOMASS_Ec_iJO1366_core_53p95M': 1}, 
#                                   precision = 'quad')


# In[1]:


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


# In[2]:


# test_model.solve_lp(mu_val = 1e10, objective = {'rA': 1})


# In[ ]:




