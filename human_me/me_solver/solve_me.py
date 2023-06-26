#!/usr/bin/env python
# coding: utf-8

import copy
import gc
import os
import time
import warnings
from collections import OrderedDict
from typing import Any, Dict, List, Optional, SupportsFloat
from itertools import repeat

import cobra
from cobra.util.array import create_stoichiometric_matrix
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from multiprocessing import Pool
# from pathos.multiprocessing import ProcessingPool as Pool
from qminospy.solver import QMINOS  # need solveME (https://github.com/SBRG/solvemepy) installed and working
from tqdm import tqdm

from human_me.core.reaction import BiomassReaction
from human_me.io import HiddenPrints

def _solve_lp_par(me_model, mu_val: SupportsFloat, objective: Optional[Dict[str, int]] = None,
             tolerance: SupportsFloat = 0, close_biomass_dilution: bool = True):
    """Solves the linear program for a specified objective at a specified growth rate.
    Same as "solve_lp" method, but for calling parallel processes in "optimize" method, for some reason does 
    not recogize the precision = self.precision assignment and errors out. This runs the solver
    with precision quad by default rather than passing the precision argument.
    """
    if objective is None:
        objective = {'biomass_dilution': 1}
    # normalize objective to be 1
    if not all([np.sign(v) == 1 for v in objective.values()]) and not all(
            [np.sign(v) == -1 for v in objective.values()]):
        raise ValueError('Current version can only maximize or minimize combinations of objectives, not do both')
    tot = abs(sum(objective.values()))
    objective = {k: v / tot for k, v in objective.items()}

    # get stoichiometric matrix at mu_val
    if type(me_model) != cobra.Model: # ME Model
        S = me_model.create_stoichiometric_matrix(mu_val=mu_val, array_type='numpy', inplace=False)
    else: # a regular metabolic model
        S = create_stoichiometric_matrix(me_model)

    # get equality and inequality matrices to format (Ev=b; lb <= Iv <= ub; A = concat[E,I])
    zero_tol = 1e-6
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

    E = S[equality_index, ]
    I = S[inequality_index, ]

    A = scipy.sparse.dok_matrix(np.concatenate([E, I, I]))  # E, L, G
    b = np.array(b_equal + b_less + b_greater)
    csense = np.array(['E'] * len(b_equal) + ['L'] * len(b_less) + ['G'] * len(b_greater))
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
        if not isinstance(r, BiomassReaction):
            xl[counter] = copy.copy(r.lower_bound) - tolerance
            xu[counter] = copy.copy(r.upper_bound) + tolerance
        else:
            if (r.id != 'biomass_dilution') or close_biomass_dilution:
                bounds = r.replace_bound_mu(mu_val=mu_val, inplace=False)
                xl[counter] = bounds[0] - tolerance
                xu[counter] = bounds[1] + tolerance
            else:
                xl[counter] = 0 - tolerance
                xu[counter] = 1000 + tolerance
        counter += 1

    # objective vector (max c.T*v)
    c = np.zeros(len(me_model.reactions))
    for r_id, coeff in objective.items():
        try:
            r_index = me_model.reactions.index(r_id)
        except:
            raise ValueError('Specified reaction id(s) not in model')
        c[r_index] = coeff

    qs = QMINOS()

    sln, stat, hsq = qs.solvelp(A, b, c, xl, xu, csense)

    # remove unwanted output files
    abspath = os.path.abspath(os.getcwd())
    for fn in [os.path.join(abspath, 'fort.9'), os.path.join(abspath, 'fort.11')]:
        if os.path.isfile(fn):
            os.remove(fn)

    return sln, stat, hsq

class qminosSolver:
    """Solving human ME Model with qMINOS."""
    def __init__(self, precision: str = 'quad') -> None:
        """"Initializer for solving with qMINOS

        Parameters
        ----------
        precision : str, optional
            The precision for the qminos solver (options ['double', 'quad', 'dq', 'dqq']), by default 'quad'
        """
        self.precision = precision

    def solve_lp(self, me_model, mu_val: SupportsFloat, objective: Optional[Dict[str, int]] = None,
                 tolerance: SupportsFloat = 0, close_biomass_dilution: bool = True):
        """Solves the linear program for a specified objective at a specified growth rate.

        Parameters
        ----------
        me_model : human_me.core.model.ME_Model
            ME Model to solve 
            (can also input a metabolic model with some of the parameter no longer being relevant)
        mu_val : SupportsFloat
            The growth value for which to solve the linear program [hr^-1]
        objective : Dict[str, int], optional
            The objective function to optimize. Dictionary represent a linear combination of reactions to optimize,
            with reaction ids as keys and the coefficient of the linear combination as the values.
            Values must either be 1 for maximization or -1 for minimization, by default {'biomass_dilution': 1}
        tolerance : SupportsFloat, optional
            Threshold below which expected sensitivity of solver is too low to detect infeasibility, by default 0
        close_biomass_dilution : bool, optional
            Internal use only, whether to constrain the biomass_dilution reaction bounds by mu, by default True

        Returns
        -------
        sln: 1D np.array
            the vector of fluxes in the optimal solution
        stat: int
            the solver status
                0     Optimal solution found.
                1     The problem is infeasible.
                2     The problem is unbounded (or badly scaled).
                3     Too many iterations.
                4     Apparent stall.  The solution has not changed
                      for a large number of iterations (e.g. 1000).
        hsq:
            optimal basis (see qminospy.solver.QMINOS)
        """
        if objective is None:
            objective = {'biomass_dilution': 1}
        # normalize objective to be 1
        if not all([np.sign(v) == 1 for v in objective.values()]) and not all(
                [np.sign(v) == -1 for v in objective.values()]):
            raise ValueError('Current version can only maximize or minimize combinations of objectives, not do both')
        tot = abs(sum(objective.values()))
        objective = {k: v / tot for k, v in objective.items()}

        # get stoichiometric matrix at mu_val
        if type(me_model) != cobra.Model: # ME Model
            S = me_model.create_stoichiometric_matrix(mu_val=mu_val, array_type='numpy', inplace=False)
        else: # a regular metabolic model
            S = create_stoichiometric_matrix(me_model)

        # get equality and inequality matrices to format (Ev=b; lb <= Iv <= ub; A = concat[E,I])
        zero_tol = 1e-6
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

        E = S[equality_index, ]
        I = S[inequality_index, ]

        A = scipy.sparse.dok_matrix(np.concatenate([E, I, I]))  # E, L, G
        b = np.array(b_equal + b_less + b_greater)
        csense = np.array(['E'] * len(b_equal) + ['L'] * len(b_less) + ['G'] * len(b_greater))
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
            if not isinstance(r, BiomassReaction):
                xl[counter] = copy.copy(r.lower_bound) - tolerance
                xu[counter] = copy.copy(r.upper_bound) + tolerance
            else:
                if (r.id != 'biomass_dilution') or close_biomass_dilution:
                    bounds = r.replace_bound_mu(mu_val=mu_val, inplace=False)
                    xl[counter] = bounds[0] - tolerance
                    xu[counter] = bounds[1] + tolerance
                else:
                    xl[counter] = 0 - tolerance
                    xu[counter] = 1000 + tolerance
            counter += 1

        # objective vector (max c.T*v)
        c = np.zeros(len(me_model.reactions))
        for r_id, coeff in objective.items():
            try:
                r_index = me_model.reactions.index(r_id)
            except:
                raise ValueError('Specified reaction id(s) not in model')
            c[r_index] = coeff

        qs = QMINOS()

        sln, stat, hsq = qs.solvelp(A, b, c, xl, xu, csense, precision=self.precision)

        # remove unwanted output files
        abspath = os.path.abspath(os.getcwd())
        for fn in [os.path.join(abspath, 'fort.9'), os.path.join(abspath, 'fort.11')]:
            if os.path.isfile(fn):
                os.remove(fn)

        return sln, stat, hsq

    def _try_mu(self, me_model, mu_val, objective, tolerance, 
                res: Dict[SupportsFloat, Dict[str, Any]], feasible_mu: List[SupportsFloat], infeasible_mu: List[SupportsFloat], verbose: bool = True):
        """To be used with maximize_growth method"""
        with HiddenPrints():
            sln, stat, hsq = self.solve_lp(me_model, mu_val, objective=objective, tolerance = tolerance)
        if stat.max() == 1 and mu_val < 1e-9:
            warnings.warn('Model is infeasible at mu = 0. Trying mu = 1e-9 instead')
            mu_val = 1e-9
            with HiddenPrints():
                sln, stat, hsq = self.solve_lp(me_model, mu_val, objective=objective, tolerance = tolerance)
            if stat.max() == 1:
                raise ValueError('Provided minimum mu is infeasible')

        res[mu_val] = {'solution': sln, 'status': stat.max(), 'basis': hsq}

        if stat.max() == 0:  # "optimal":
            if verbose:
                print('The problem has an optimal solution at mu = {} (hrs)'.format(mu_val))
            feasible_mu.append(mu_val)
            return True, sln, stat, hsq, res, feasible_mu, infeasible_mu
        if stat.max() == 1:
            infeasible_mu.append(mu_val)
            if verbose:
                print('The problem is infeasible at mu = {} (hrs)'.format(mu_val))
            return False, None, None, None, res, feasible_mu, infeasible_mu
        raise ValueError('The problem returned with stat: {}'.format(stat.max()))

    def maximize_growth(self, me_model, min_mu: SupportsFloat = 0, max_mu: SupportsFloat = 0.05,
                        mu_accuracy: SupportsFloat = 1e-10, increment: SupportsFloat = 0.02, tolerance: SupportsFloat = 0,
                        verbose: bool = True):
        """Binary search to find the maximum feasible growth rate.

        Parameters
        ----------
        me_model : human_me.core.model.ME_Model
            ME Model to solve
        min_mu : SupportsFloat, optional
            Expected minimum feasible growth rate, by default 0
        max_mu : SupportsFloat, optional
            Expected minimum infeasible growth rate (i.e., just above expected maximum feasible growth rate), by default 0.05
        mu_accuracy : SupportsFloat, optional
            The maximum error in mu after the binary search, by default 1e-10
        increment : SupportsFloat, optional
            The amount to increase growth by when searching for maximum infeasible growth rate from max_mu, by default 0.02
        tolerance : SupportsFloat, optional
            Threshold below which expected sensitivity of solver is too low to detect infeasibility, by default 0
        verbose : bool, optional
            Prints information about each linear program iteration, by default True

        Returns
        -------
        mu_max: int
            the maximum feasible growth value (in hours)
        res: Dict[float, Tuple[np.array, int]]]
            keys are all attempted growth values, values are dictionaries with keys as output from self.solve_lp
        """
        objective = {'biomass_dilution': 1}  # maximizing for growth
        feasible_mu = []
        infeasible_mu = []
        res = dict()

        start = time.time()

        if verbose:
            print('Trying mu: {}'.format(min_mu))
        bool_, sln, stat, hsq, res, feasible_mu, infeasible_mu = self._try_mu(me_model, min_mu, objective = objective, tolerance = tolerance, res = res, 
                                                                            feasible_mu = feasible_mu, infeasible_mu = infeasible_mu, verbose = verbose)  # start with minimal


        bool_max, sln, stat, hsq, res, feasible_mu, infeasible_mu = self._try_mu(me_model, max_mu, objective = objective, tolerance = tolerance, res = res, 
                                                                            feasible_mu = feasible_mu, infeasible_mu = infeasible_mu, verbose = verbose)
        while bool_max:  # If max_mu was feasible, keep increasing
            max_mu += increment
            if verbose:
                print('Trying mu: {}'.format(max_mu))
            bool_max, sln, stat, hsq, res, feasible_mu, infeasible_mu = self._try_mu(me_model, max_mu, objective = objective, tolerance = tolerance, 
                                                                                    res = res, feasible_mu = feasible_mu, infeasible_mu = infeasible_mu, verbose = verbose)
        while (infeasible_mu[-1] - feasible_mu[-1]) > mu_accuracy:
            if verbose:
                print('Trying mu: {}'.format((infeasible_mu[-1] - feasible_mu[-1]) * 0.5))
            bool_, sln, stat, hsq, res, feasible_mu, infeasible_mu = self._try_mu(me_model, (infeasible_mu[-1] + feasible_mu[-1]) * 0.5,  
                                                                                objective = objective, tolerance = tolerance, res = res, 
                                                                                feasible_mu = feasible_mu, infeasible_mu = infeasible_mu, verbose = verbose)
        if verbose:
            tot = ((time.time() - start) / 3600)
            print("completed in {:.2f} hours and {} iterations".format(tot, len(feasible_mu + infeasible_mu)))

        mu_max = np.max(feasible_mu)
        res_ = OrderedDict({k: res[k] for k in sorted(list(res.keys()))})

        return mu_max, res_

    def optimize(self, me_model, objective: Dict[str, int], mu_max: SupportsFloat,
                 n_points: int = 10, tolerance: SupportsFloat = 0, n_cores: Optional[int] = None,
                 visualize: bool = True, fig_name: Optional[str] = None):
        """General optimization of any non-growth objective.

        Parameters
        ----------
        me_model : human_me.core.model.ME_Model
            ME Model to solve
        objective : Dict[str, int]
            The objective function to optimize. Dictionary represent a linear combination of reactions to optimize,
            with reaction ids as keys and the coefficient of the linear combination as the values.
            Values must either be 1 for maximization or -1 for minimization.
        mu_max : SupportsFloat
            the maximum growth value at which the model is feasible [hr^-1]; use .maximize_growth() method to identify (should be <= mu_max output)
            if using an experimental value, make sure it is feasible using the .solve_lp() method 
        n_points : int, optional
            # of growth values to consider between 0 and mu_max, by default 10
        tolerance : SupportsFloat, optional
            Threshold below which expected sensitivity of solver is too low to detect infeasibility, by default 0
        n_cores : Optional[int], optional
            the number of workers to use for parallelization, by default None
        visualize : bool, optional
            plot the relationship between growth and the objective function of interest, by default True
        fig_name : Optional[str], optional
            save the plotted figure to 'path/to/filename.ext', by default None

        Returns
        -------
        sln: Tuple[float]
            first element is the growth value at which the non-growth objective is optimized
            second element is the optimized non-growth objective value
        predicted: pandas.DataFrame
            1000 growth values between 0 and mu_max, with corresponding interpolated objective values
        interp_fit: scipy.interpolate.interp1d
            a function to interpolate objective values from growth values, used to generated predicted
        optimal_vals: collections.OrderedDict
            keys are n_points growth values between 0 and mu_max, values are the objective value optimized at
            the corresponding growth value
        res: Dict[float, Tuple[np.array, int]]
            keys are n_points growth values between 0 and mu_max, values are the output of .solve_lp at
            corresponding growth values with the objective set to the non-growth objective input
        """
        obj_keys = list(objective.keys())
        if len(obj_keys) == 1 and obj_keys[0] == 'biomass_dilution':
            raise ValueError('To optimize for growth, use the .maximize_growth() method')
        if len(set(list(objective.values()))) != 1:
            raise ValueError(
                'Currently, NLP objectives must be either only maximization or minimization (only all 1s or only all -1s in objective dictionary values)')

        growth_vals = np.arange(0, mu_max + mu_max / n_points, mu_max / (n_points - 1))
        if (n_cores is None) or (n_cores <= 1):
            res = list()
            for mu_val in tqdm(growth_vals):
                sln, stat, hsq = self.solve_lp(me_model=me_model, mu_val=mu_val, objective=objective,
                                               tolerance=tolerance,
                                               close_biomass_dilution=True)
                res.append([sln, stat, hsq])
        else:
            n_cores = min([n_cores, n_points])
            pool = Pool(n_cores)
            try:
                res = pool.starmap(_solve_lp_par, 
                                    zip(repeat(me_model), list(growth_vals), repeat(objective),
                                        repeat(tolerance), repeat(True)))
                # res = pool.map(self.solve_lp, [me_model] * n_points, list(growth_vals), [objective] * n_points,
                #                [tolerance] * n_points,
                #                [True] * n_points)
                pool.close()
                pool.join()
                # pool.restart()
                gc.collect()
            except:
                pool.close()
                pool.join()
                # pool.restart()
                gc.collect()
                raise ValueError('Parallelization failed')
            res = [list(r) for r in res]

        reaction_indeces = [me_model.reactions.index(j) for j in sorted(objective.keys())]
        res = OrderedDict({i[0]: dict(zip(['val', 'sln', 'stat', 'hsq'],
                                          [(dict(zip(sorted(objective.keys()), i[1][0][reaction_indeces])))] + i[1]))
                           for i in zip(growth_vals, res)})
        optimal_vals = OrderedDict({k: sum(v['val'].values()) for k, v in res.items()})

        # estimate objective values across growth using interpolation
        interp_fit = scipy.interpolate.interp1d(x=list(optimal_vals.keys()), y=list(optimal_vals.values()),
                                                bounds_error=False)
        obj_label = 'predicted_' + '_'.join(list(objective.keys()))
        if len(obj_label) > 30:
            obj_label = 'predicted_Objective'
        predicted = pd.DataFrame(data={'growth': np.arange(0, mu_max + mu_max / 1000, mu_max / (1000 - 1))})
        predicted[obj_label] = predicted.growth.apply(lambda x: interp_fit(x).item()).values.tolist()

        if np.sign(list(objective.values())[0]) == 1:
            optimal_val = predicted[obj_label].max()
        elif np.sign(list(objective.values())[0]) == -1:
            optimal_val = predicted[obj_label].min()

        optimal_val_growth = predicted[predicted[obj_label] == optimal_val].growth.values.tolist()[0]
        sln = (optimal_val_growth, optimal_val)

        if visualize:
            fig, ax = plt.subplots(figsize=(8, 5))
            sns.lineplot(x='growth', y=obj_label, data=predicted, ax=ax)
            sns.scatterplot(x=list(optimal_vals.keys()), y=list(optimal_vals.values()), color='black', ax=ax)
            plt.plot([optimal_val_growth], [optimal_val], marker='o', markersize=3, color="red")
            ax.legend(handles=[mpatches.Patch(color='black', label='Solved'),
                               mpatches.Patch(color=sns.color_palette('tab10')[0], label='Interpolated')],
                      fancybox=True, fontsize=12, bbox_to_anchor = [1,1])
            ax.set_xlabel(r'$\mu$ $[hr^{-1}]$', fontsize=15, labelpad=5)
            ax.set_ylabel(ax.get_ylabel().split('predicted_')[1] + r'$\;\frac{mmol}{gDw \;hr}$', fontsize=15,
                          labelpad=5)
            ax.tick_params(axis='both', labelsize=12)

            ax.vlines(x=optimal_val_growth, ymin=ax.get_ylim()[0], ymax=optimal_val, linestyles='--',
                      color=sns.color_palette('pastel')[3])
            ax.hlines(y=optimal_val, xmin=ax.get_xlim()[0], xmax=optimal_val_growth, linestyles='--',
                      color=sns.color_palette('pastel')[3])

            fig.tight_layout()
            if fig_name is not None:
                fig_name = os.path.splitext(fig_name)[0]
                for ext in ['.png', '.pdf', '.svg']:
                    plt.savefig(fig_name + ext, bbox_inches='tight')

        return sln, predicted, interp_fit, optimal_vals, res
