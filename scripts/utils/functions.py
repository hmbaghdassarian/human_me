#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
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


# In[5]:


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


# In[6]:


def parse_me_reaction_id(x):
    '''Get HGNC ID associated with an expression module reaction
    
    Parameters
    ----------
    x: str
        cobra.Reaction.id
    
    Returns
    ----------
    parsed x: str
        HGNC ID
    '''
    
    
    if 'HGNC' in x.split('_')[0]:
        return '_'.join(x.split('_')[1:])
    else:
        return x


# In[7]:


def SASA(mw):
    '''Estimate the protein solvent-accessible surface area from the molecular weight
    
    Parameters
    ----------
    mw: float
        protein molecular weight (in kDa)
    
    Returns
    ----------
    SASA: float
        approximate protein solvent accesible surface area
    '''
    return mw**(0.75)


# In[9]:


def average_protein_features(psim_me, protein_ids, context_specific = True):
    '''Function to get the average protein features from the proteins used in a specific ME model being generated.
    This is explicitly written to help generate the dummy protein. 
    
    Parameters
    ----------
    psim_me: pd.DataFrame
        protein specific information matrix, same as corrected input file (see preprocessing output)
    protein_ids: list
        each entry is a string protein HGNC ID, the list should include all proteins included in the ME Model
    context_specific: bool, default True
        whether to use the use provided input M-model and PSIM, or get an average dummy representative of the whole proteome


    Returns
    ----------
    dummy_psim: pd.DataFrame
        same as PSIM but with one row, representing the average features of all proteins 
    '''
    if context_specific:
        psim = psim_me.copy()
    else:
        psim = pd.read_csv(build_files_path + 'psim_recon2_2.csv')
    
    res = pd.DataFrame()
    res['premrna_counts'] = psim.PREMRNA_SEQ.dropna().apply(lambda x: {ntp: x.count(ntp)for ntp in set(x)})
    res['premrna_length'] = psim.PREMRNA_SEQ.dropna().apply(lambda x: len(x))
    res['premrna_prop'] = res.apply(lambda x: {k: v/x.premrna_length for k,v in x.premrna_counts.items()}, axis = 1)

    premrna_L = res['premrna_length'].median()
    premrna_avg_prop = {ntp: res['premrna_prop'].apply(lambda x: x[ntp]).median() for ntp in ['A', 'U', 'C', 'G']}

    premrna_seq = ''
    for ntp in ['A', 'U', 'C', 'G']:
        premrna_seq += ntp*int(round(premrna_avg_prop[ntp]*premrna_L))


    res = pd.DataFrame()
    res['mrna_counts'] = psim.MRNA_SEQ.dropna().apply(lambda x: {ntp: x.count(ntp)for ntp in set(x)})
    res['mrna_length'] = psim.MRNA_SEQ.dropna().apply(lambda x: len(x))
    res['mrna_prop'] = res.apply(lambda x: {k: v/x.mrna_length for k,v in x.mrna_counts.items()}, axis = 1)

    mrna_L = res['mrna_length'].median()
    mrna_avg_prop = {ntp: res['mrna_prop'].apply(lambda x: x[ntp]).median() for ntp in ['A', 'U', 'C', 'G']}

    mrna_seq = ''
    for ntp in ['A', 'U', 'C', 'G']:
        mrna_seq += ntp*int(round(mrna_avg_prop[ntp]*mrna_L))

    res = pd.DataFrame()
    res['protein_counts'] = psim.PROTEIN_SEQ.dropna().apply(lambda x: {ntp: x.count(ntp)for ntp in set(x)})
    res['protein_length'] = psim.PROTEIN_SEQ.dropna().apply(lambda x: len(x))
    res['protein_prop'] = res.apply(lambda x: {k: v/x.protein_length for k,v in x.protein_counts.items()}, axis = 1)

    protein_L = int(round(mrna_L)/3)

    def get_prop(x, aa):
        if aa in x.keys():
            return x[aa]
        else:
            return 0

    protein_avg_prop = {aa: res['protein_prop'].apply(lambda x: get_prop(x, aa)).median() for aa in params.amino_acids}

    protein_seq = ''
    for aa in params.amino_acids:
        protein_seq += aa*int(round(protein_avg_prop[aa]*protein_L))

    dummy_psim = pd.DataFrame(columns = psim.columns)
    dummy_psim.loc[0,:] = float('nan')
    dummy_psim.loc[0,['HGNC_ID', 'PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ']] = ['HGNC:DUMMY', premrna_seq, mrna_seq, protein_seq] 
    dummy_psim.LOCATION = ['[c]']
    
    # secretory args will be disregarded anyways for now
    median_vals = ['POLYA_LENGTH', 'N_EXONS', 'ALPHA_M', 'ALPHA_P', 'TMD', 'DSB', 'OG', 'NG']
    for col in median_vals:
        dummy_psim[col] = psim[col].median()
        
    
    argmax_vals = ['SP', 'GPI']
    
    # deal with PTR column
    if psim.PTR.dropna().convert_dtypes().dtype is np.dtype('float64'):
        dummy_psim.PTR = psim.PTR.median()
    elif isinstance(psim.PTR.dropna().convert_dtypes().dtype, pd.StringDtype):
        argmax_vals += ['PTR']
    
    for col in argmax_vals:
        if psim[col].dropna().shape[0] > 0:
            val = psim[col].dropna().value_counts().idxmax()
        else:
            val = float('nan')
        dummy_psim[col] = val    
        
    return dummy_psim

