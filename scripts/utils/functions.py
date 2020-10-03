#!/usr/bin/env python
# coding: utf-8

# In[10]:


import cobra
import pandas as pd
import os

import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils.load_environmental_variables import *
from utils import metabolites as metab
from utils import parameters as params


# # Functions

# In[ ]:


# def blockPrint():
#     sys.stdout = open(os.devnull, 'w')
# def enablePrint():
#     sys.stdout = sys.__stdout__

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


# In[ ]:


def get_reaction_compartment(reaction):
    '''Input is a cobra.Reaction, output is a singular compartment. This function maps reactions to a particular 
    compartment according to some rules'''
    
    compartments_ = list(reaction.compartments.copy())
    if len(compartments_) > 1: # for reactions that occur in more than one compartment
        if 'c' in compartments_ and len(compartments_) == 2: # remove cytoplasmic compartment between the two for machinery
            compartments_.remove('c')
        else: # choose most common compartment 
            compartments_ = [max(compartments_, key = compartments_.count)]
    if len(compartments_) > 1:
        raise ValueError('Failed to map reaction to a singular compartment')
    elif compartments_[0] not in params.compartments.keys():
        raise ValueError('Mapped reaction to a compartment that is not allowed in ME model')
    else:
        return compartments_[0]


# In[ ]:


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


# In[ ]:


def make_rna_metabolite(metabolite_name, seq, molecule_type, compartment = 'n', triphosphate = True):
    
    '''
    Inputs:
    1) metabolite_name is the name of the RNA molecule (unique ID)
    2) seq is a string representing the one-letter sequence of the RNA molecule.
    3) molecule type = ['mrna', 'trna', 'rrna']
    4) compartment is the one-letter string representing the location of the RNA molecule (usually 'n' or 'c')
    5) triphosphate is a boolean. If true, assume the 5' end has a triphosphate, otherwise assume monophosphate.
    
    Outputs:
    1) rna_metabolite is an object of cobra.Metabolite representing teh RNA molecule
    2) base_counts is a dictionary with keys as one-letter RNA base characters and values as the # of 
    occurences of that base in the RNA sequence 
    
    '''
    
    if molecule_type not in ['mrna', 'trna', 'rrna']:
        raise ValueError('molecule_type must be mrna, trna, or rrna')

    rna_metabolite = cobra.Metabolite(metabolite_name + '_' + molecule_type + '[' + compartment + ']')
    rna_metabolite.compartment = compartment
    base_counts, elements = get_base_counts_and_elements(seq, triphosphate = triphosphate) # utils function

    rna_metabolite.elements = elements
    rna_metabolite.charge = -len(seq)
    
    if triphosphate:
        rna_metabolite.charge -= 3
    
    return rna_metabolite, base_counts


# In[ ]:


def rna_exonucleolytic_degradation(rna_metabolite, rna_base_counts, rna_sequence, reaction_name, 
                                   triphosphate = True, nucleus = True):
    ''' 
    
    Generates a reaction for exonucleolytic cleavage of an RNA molecule (RNA-->NMPs).
    Inputs:
    1) rna_metabolite is a cobra.Metabolite object representing the rna molecule to be degraded.
    2) rna_base_counts is a dictionary with keys as one-letter RNA base characters and values as the # of 
    occurences of that base in the RNA sequence.
    3) rna_sequence is the ordered (5'-->3') sequence of the RNA molecule, as a string of one-letter bases
    4) reaction_name is a string representing the name you want to give the reaction
    5) triphosphate is a boolean. If true, assume the 5' end has a triphosphate, otherwise assume monophosphate. 
    6) nucleus is a boolean. If true, the degradation reaction is taking place in the nucleus. Otherwise, assume
    it takes place in the cytoplasm (current iteration of model only degrades RNA in these two compartments).
    
    Output: a degradation reaction of type cobra.Reaction
    no GPRs or subsystems added to reaction
    
    '''
    # exonucleolytic cleavage of RNA reaction

    if nucleus: 
        rna_degradation = cobra.Reaction(reaction_name + '_DEGRADATIONn')
        rxn = dict()
        rxn[metab.h2o_n] = -sum(rna_base_counts.values())+1
        rxn[rna_metabolite] = -1
        for k,v in metab.nmp_map_n.items():
            rxn[v] = rna_base_counts[k]

        # triphosphate on 5' end
        if triphosphate:
            rxn[metab.nmp_map_n[rna_sequence[0]]] -= 1
            rxn[metab.ntp_map_n[rna_sequence[0]]]  = 1  
            rxn[metab.h_n] = sum(rna_base_counts.values())-1
        else:
            rxn[metab.h_n] = sum(rna_base_counts.values()) # extra H on 5' end <--unsure about this

        rna_degradation.add_metabolites(rxn)

        
    else:
        rna_degradation = cobra.Reaction(reaction_name + '_DEGRADATIONc')
        rxn = dict()
        rxn[metab.h2o_c] = -sum(rna_base_counts.values())+1
        rxn[rna_metabolite] = -1
        for k,v in metab.nmp_map_c.items():
            rxn[v] = rna_base_counts[k]

        # triphosphate on 5' end
        if triphosphate:
            rxn[metab.nmp_map_c[rna_sequence[0]]] -= 1
            rxn[metab.ntp_map_c[rna_sequence[0]]]  = 1  
            rxn[metab.h_c] = sum(rna_base_counts.values())-1
        else:
            rxn[metab.h_c] = sum(rna_base_counts.values()) # extra H on 5' end <--unsure about this

        rna_degradation.add_metabolites(rxn)
        
    return rna_degradation 


# In[ ]:


def make_protein_metabolite(id_, amino_acid_counts, L_protein, compartment):
    '''
    ID is a string to name the protein metabolite. 
    Amino acid counts is a dictionary with keys as the aa one letter code and counts as the number of occurences of that amino acid in the protein sequence
    L_protein is the length of the protein
    Compartment is the location of the protein (one letter string, corresponds to Recon2.2s compartments)
    
    Will return a cobra.Metabolite object with relevant charge and elements.
    
    '''
    
    protein_metabolite = cobra.Metabolite(id_ + '_protein[' + compartment + ']')
    protein_metabolite.compartment = compartment
    
    elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}
    if compartment in metab.seq_amino_acid_map_compartments.keys():
        for aa_code, aa_count in amino_acid_counts.items():
            aa_elements = metab.seq_amino_acid_map_compartments[compartment][aa_code].elements
            for element in aa_elements:
                elements[element] += aa_count*aa_elements[element]
    else:
        raise ValueError('Must add this compartment to make_protein_metabolite function')

    # peptide bond formation
    elements['H'] -= 2*(L_protein-1)
    elements['O'] -= 1*(L_protein-1)

    protein_metabolite.elements = elements
    # assume charge of amino acid is the ssame regardless of metabolite
    protein_metabolite.charge = sum([metab.seq_amino_acid_map_compartments[compartment][aa_code].charge*aa_count for aa_code, aa_count in amino_acid_counts.items()])
    return protein_metabolite


# In[1]:


def make_complex_metabolite(complex_id = None, **complex_info):# metabolites, *ids, *metabolite_types):
    '''
    
    Inputs:
    Complex info is a dictionary with three keys ['METABOLITES', 'IDS', 'METABOLITE_TYPES']
    Each value is a list:
        metabolites is a list of cobra.Metabolite objects
        IDs is a list of string identifiers corresponding to each metabolite object
        Metabolite_types is a list of strings; possible values are ['protein', 'rrna', 'trna', 'mrna',  'metabolite']
        This means complexes can form between any of these species, including other complexes; metabolite is a M-model metabolite
    complex_id is a string for the id of the complex metabolite, otherwise will form one from metabolite ids
    Output:
    A cobra.Metabolite object representing the complex formed between metabolites
    
    '''
    if sorted(set(complex_info.keys())) != ['IDS', 'METABOLITES', 'METABOLITE_TYPES']:
        raise ValueError('Invalid complex information keys or insufficient complex information keys')
    
    metabolites_, ids, metabolite_types = complex_info['METABOLITES'], complex_info['IDS'], complex_info['METABOLITE_TYPES']
    
    if len(set(metabolite_types).difference(['protein', 'rrna', 'trna', 'mrna',  'metabolite', 'complex']))>1:
        raise ValueError('At least one of the metabolite types is not considered in complex formation currently')
    
    
    compartments = list(set([m.compartment for m in metabolites_]))
    if len(compartments) == 1:
        compartment = compartments[0]
    # exception of ribosome complex
    elif (len(compartments) == 2) and ('c' in compartments) and ('mature_ribosome_complex_complex[c]' in [m.id for m in metabolites_]):
        compartment = 'c'
    else:
        raise ValueError('metabolites are not in the same compartment')
    
    mt_type = '_'.join(list(set(metabolite_types)))
    
    ids_ = '_'.join(ids)
    
    if complex_id == None:
        id_ = ids_ + '_' + mt_type
    else: 
        id_ = complex_id + '_' + mt_type
    
    complex_id = id_ + '_complex' + '[' + compartment + ']'
    if len(complex_id)>(256-8-4-len(mt_type)): #-8 and -4 for _complex and compartment appended to end
        err_msg = 'Cobrapy requires metabolite ids to be less than 256 characters, please specify a '
        err_msg += 'shorter user-defined complex id'
        raise ValueError(err_msg)
    
        
    complex_metabolite = cobra.Metabolite(complex_id)
    complex_metabolite.compartment = compartment
    complex_metabolite.charge = sum([m.charge for m in metabolites_])
    
    elements = dict()
    for m in metabolites_:
        for k,v in m.elements.items():
            if k in elements.keys():
                elements[k] += v
            else:
                elements[k] = v
    complex_metabolite.elements = elements
    
    return complex_metabolite, id_

def form_complex(reaction_id = None, complex_id = None, **complex_info):
    
    '''
    
    Inputs:
    Complex info is a dictionary with three keys ['METABOLITES', 'IDS', 'METABOLITE_TYPES']
    Each value is a list:
        metabolites is a list of cobra.Metabolite objects
        IDs is a list of string identifiers corresponding to each metabolite object
        Metabolite_types is a list of strings; possible values are ['protein', 'rrna', 'trna', 'mrna',  'metabolite', 'complex']
  
    Output:
    A cobra.Reaction object representing the complex formation between metabolites
    
    '''
    
    complex_metabolite, id_ = make_complex_metabolite(complex_id, **complex_info)
    metabolites_ = complex_info['METABOLITES']
    compartment = list(set([m.compartment for m in metabolites_]))[0]

    if reaction_id == None:
        reaction_id = id_ + '_COMPLEX_FORMATION' + compartment
    else:
        reaction_id = reaction_id + '_COMPLEX_FORMATION' + compartment
    complex_formation = cobra.Reaction(reaction_id)
    
    rxn = {m: -1 for m in metabolites_}
    rxn[complex_metabolite] = 1
    complex_formation.add_metabolites(rxn)
    complex_formation.lower_bound = -1000 # reversible
    
    return complex_formation, complex_metabolite


# In[1]:


def get_metabolite_mw(metabolite, no_copies = 1, metabolite_elements = None, 
                      element_mw = {'C': 0.0120107, 'H': 0.00100784, 'N': 0.0140067, 'O': 0.015999, 
                                    'P': 0.030973762, 'S': 0.032065}):
    '''Input is a cobra.Metabolite object. 
    Alternatively, a dicitonary (metabolite_elements) with elements as keys and element counts as values can be provided. 
    The cobra.Metabolite object takes precedent over the dictionary if both are provided.
    no_copies is the number of molecules of that metabolite (i.e., stoichiometric coefficient in the reaction)
    output is the molecular weight of that metabolite in kDa'''
    
    if metabolite != None:
        return no_copies*sum([element_mw[element]*count for element, count in metabolite.elements.items()])

    else:
        if metabolite_elements != None:
            mw = no_copies*sum([element_mw[element]*count for element, count in metabolite_elements.items()])
        else:
            raise ValueError('Must provide a cobra.Metabolite object or dictionary of elements')
    

