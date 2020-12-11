#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import itertools

import os
import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils import metabolites as metab
from utils import parameters as params
from utils import machinery as mach


# # Functions

# In[2]:


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


# In[3]:


def get_reaction_compartment(reaction):
    '''Input is a cobra.Reaction, output is a singular compartment. This function maps reactions to a particular 
    compartment according to some rules'''
    
    compartments_ = sorted([m.compartment for m in reaction.metabolites.keys() if m.compartment is not None])
    # sorted to choose the first one in alphabetical order given a tie
    if len(set(compartments_)) > 1: # for reactions that occur in more than one compartment
        if 'c' in compartments_: # remove cytoplasmic compartment as a choice in multi-machinery
            compartments_ = sorted([c for c in compartments_ if c != 'c'])
        if len(set(compartments_)) > 1:
            compartments_ = [max(compartments_, key = compartments_.count)]
    
    compartments_ = sorted(set(compartments_))
    if len(compartments_) != 1:
        raise ValueError('Failed to map reaction to a singular compartment')
    elif compartments_[0] not in params.compartments.keys():
        raise ValueError('Mapped reaction to a compartment that is not allowed in ME model')
    else:
        return compartments_[0]


# In[4]:


def hydrolyze_atp(rxn, n_atp, compartment):
    '''
    Rxn is a dict for the cobra.Reaction.add_metabolite function.
    n_atp is the # of atp to hydrolyze
    compartment is the compartment for hydrolysis
    
    '''
    n_atp = round(n_atp)
    
    if metab.atp_compartments[compartment] in rxn.keys():
        rxn[metab.atp_compartments[compartment]] -= n_atp 
    else:
        rxn[metab.atp_compartments[compartment]] = -n_atp 

    if metab.h2o_compartments[compartment] in rxn.keys():
        rxn[metab.h2o_compartments[compartment]] -= n_atp 
    else:
        rxn[metab.h2o_compartments[compartment]] = -n_atp 

    if metab.adp_compartments[compartment] in rxn.keys():
        rxn[metab.adp_compartments[compartment]] += n_atp 
    else:
        rxn[metab.adp_compartments[compartment]] = n_atp

    if metab.pi_compartments[compartment] in rxn.keys():
        rxn[metab.pi_compartments[compartment]] += n_atp 
    else:
        rxn[metab.pi_compartments[compartment]] = n_atp

    if metab.h_compartments[compartment] in rxn.keys():
        rxn[metab.h_compartments[compartment]] += n_atp 
    else:
        rxn[metab.h_compartments[compartment]] = n_atp
    
    return rxn


# In[ ]:


def get_base_counts_and_elements(seq, triphosphate = True):
    '''
    
    Inputs:
    1) Seq is a Bio.Seq object or a string representing an RNA sequence. 
    2) triphosphate is a boolean. If true, assume the 5' end has a triphosphate, otherwise assume monophosphate. 
   
   Outputs:
    1) base_counts is a dictionary with keys as one-letter RNA base characters and values as the # of 
    occurences of that base in the RNA sequence
    2) elements is a dictionary emulating cobra.Metabolite.elements
   
   '''
    base_counts = dict()
    for base_letter in metab.seq_element_map.keys():
        base_counts[base_letter] = seq.count(base_letter)
        
    elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0}
    for base_letter in metab.seq_element_map.keys():
        for element in elements.keys():
            elements[element] += base_counts[base_letter]* metab.seq_element_map[base_letter][element]   
    
    #3' OH end
    elements['H'] += 1 
    elements['O'] += 1
    
    # 5' end
    if triphosphate:
        elements['P'] += 2
        elements['O'] += 6
    else:
        elements['H'] += 1
      
        
    return base_counts, elements


# In[12]:


def parse_me_reaction_id(x):
    if 'HGNC' in x.split('_')[0]:
        return '_'.join(x.split('_')[1:])
    else:
        return x


# In[14]:


def SASA(mw):
    '''Estimate the protein solvent-accessible surface area from the molecular weight'''
    return mw**(0.75)

