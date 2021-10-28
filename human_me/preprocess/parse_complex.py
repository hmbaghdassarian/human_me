#!/usr/bin/env python
# coding: utf-8

import cobra
from cobra.core.gene import parse_gpr
import ast
import itertools

def eval_complex_recur_full(expr):
    """
    
    Recursive parsing of gprs into lists of complexes. Input expr is a cobra.parse_gpr(gpr_string), 
    output is a list of lists, each entry of which is a complex joined by 'AND'; netsted lists are joined 
    with each other by 'OR'. 
    
    Inspired by corda source code, should cite them.
    
    """
    
    # corda: https://github.com/resendislab/corda/blob/master/corda/util.py
    if isinstance(expr, ast.Expression):
        return eval_complex_recur_full(expr.body) # Here!
    elif isinstance(expr, ast.Name):
        return cobra.core.gene.ast2str(expr)
    elif isinstance(expr, ast.BoolOp):
        op = expr.op
        if isinstance(op, ast.Or):
            return [eval_complex_recur_full(i) for i in expr.values] # Here!
        elif isinstance(op, ast.And):
            or_op = False
            names = []
            bool_ors = []
            bool_ands = []
            for e in expr.values:
                if isinstance(e, ast.BoolOp):
                    if isinstance(e.op, ast.Or):
                        or_op = True
                        bool_ors.append(e)
                    else:
                        bool_ands.append(e)
                elif isinstance(e, ast.Name):
                    names.append(eval_complex_recur_full(e))
            
            if or_op:
                product = []
                if len(bool_ors) > 1:
                    product = list(itertools.product(*[eval_complex_recur_full(i) for i in bool_ors]))
                    
                    
                ba_lists = []
                for ba in bool_ands:
                    ba_list = [eval_complex_recur_full(j) for i in bool_ands for j in i.values]
                    ba_lists.append(ba_list)
                
                result = []    
                if len(ba_lists) > 0:
                    bal_results = [eval_complex_recur_full(bal) for bal in ba_lists]
                    for br in bal_results:
                        if len(product) == 0:
                            result += [[eval_complex_recur_full(j)] + names + br for i in bool_ors for j in i.values] # Here!
                        else:
                            result += [list(p) + names + br for p in product]
                else:
                    if len(product) == 0:
                        result += [[eval_complex_recur_full(j)] + names for i in bool_ors for j in i.values] # Here!
                    else:
                        result += [list(p) + names for p in product]
                return result
                    
            else:
                return [eval_complex_recur_full(i) for i in expr.values] # Here!

def unnested_list(nested_list):
    lists = []
    non_lists = []
    for l in nested_list:
        if isinstance(l, list):
            lists.append(l)
        else:
            non_lists.append(l)
    size = len(lists)
    if size == len(nested_list):
        new_list = []
        for l in lists:
            new_list += unnested_list(l)
        return new_list
    elif size == 0:
        return [nested_list]
    else:
        return [item for sublist in [unnested_list(l) for l in lists] for item in sublist] + non_lists
    
def eval_complex(expr):
    """Input is a gene reaction rule, output is a list. If the length of the list is longer than 1, these are
    due to the presence of OR. If the entry itself is a list longer than one, this entry is a complex."""
    
    complexes = unnested_list(eval_complex_recur_full(parse_gpr(expr)[0]))
    
    # make sure machinery is always output in the same order
    for i in range(len(complexes)):
        if isinstance(complexes[i],list):
            complexes[i] = sorted(complexes[i])
    return complexes

