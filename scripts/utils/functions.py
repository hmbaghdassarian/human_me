#!/usr/bin/env python
# coding: utf-8

# In[20]:


import cobra
from cobra.core.gene import parse_gpr

import pandas as pd
import itertools
import ast


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
    


# In[8]:


class COMPLEX(cobra.Metabolite):
    def __init__(self, complex_id = None, reaction_id = None, **complex_info):
        '''
        Inputs:
        Complex info is a dictionary with three keys ['METABOLITES', 'IDS', 'METABOLITE_TYPES']
        Each value is a list:
            metabolites is a list of cobra.Metabolite objects
            IDs is a list of string identifiers corresponding to each metabolite object
            component_types is a list of strings; possible values are ['protein', 'rrna', 'trna', 'mrna', 'metabolite']
            This means complexes can form between any of these species, including other complexes; metabolite is a M-model metabolite
        complex_id is a string for the id of the complex metabolite, otherwise will form one from metabolite ids

        Output:
        A cobra.Metabolite object representing the complex formed between metabolites stored in self.complex_metabolite
        
        '''
        # checks
        if sorted(set(complex_info.keys())) != ['IDS', 'METABOLITES', 'METABOLITE_TYPES']:
            raise ValueError('Invalid complex information keys or insufficient complex information keys')

        self.subcomponents, ids, self.component_types = complex_info['METABOLITES'], complex_info['IDS'], complex_info['METABOLITE_TYPES']
        
        if len(set(self.component_types).difference(['protein', 'rrna', 'trna', 'mrna',  'metabolite', 'complex']))>1:
            raise ValueError('At least one of the metabolite types is not considered in complex formation currently')
        
        if len(self.subcomponents) != len(self.component_types):
            raise ValueError('Each metabolite must have a corresponding metabolite type')
        # parse compartment    
        compartments = list(set([m.compartment for m in self.subcomponents]))
        if len(compartments) == 1:
            compartment = compartments[0]
        # exception of ribosome complex
        elif (len(compartments) == 2) and ('c' in compartments) and ('mature_ribosome_complex_complex[c]' in [m.id for m in self.subcomponents]):
            compartment = 'c'
        else:
            raise ValueError('metabolites are not in the same compartment')

        # parse metabolite id
        mt_type = '_'.join(list(set(self.component_types)))
        ids_ = '_'.join(ids)
        if complex_id == None:
            id_ = ids_ + '_' + mt_type
        else: 
            id_ = complex_id + '_' + mt_type
            
        if reaction_id == None:
            self.reaction_id = id_ + '_COMPLEX_FORMATION' + compartment
        else:
            self.reaction_id = reaction_id + '_COMPLEX_FORMATION' + compartment
            
        
        complex_id = id_ + '_complex' + '[' + compartment + ']'
        if len(complex_id)>(256-8-4-len(mt_type)): #-8 and -4 for _complex and compartment appended to end
            err_msg = 'Cobrapy requires metabolite ids to be less than 256 characters, please specify a '
            err_msg += 'shorter user-defined complex id'
            raise ValueError(err_msg)
        
        # make the metabolite
        cobra.Metabolite.__init__(self, id = complex_id, 
                                  compartment = compartment, charge = sum([m.charge for m in self.subcomponents]))
        elements = dict()
        for m in self.subcomponents:
            for k,v in m.elements.items():
                if k in elements.keys():
                    elements[k] += v
                else:
                    elements[k] = v
        self.elements = elements
        
    def form_complex(self):

        '''
        Output: A cobra.Reaction object representing the complex formation between metabolites stored in self.complex_formation

        '''
        
        # expected no biomass change
        complex_formation = cobra.Reaction(self.reaction_id)
        rxn = {m: -1 for m in self.subcomponents}
        rxn[self] = 1
        complex_formation.add_metabolites(rxn)
        complex_formation.lower_bound = -1000 # reversible
        
        return complex_formation
#         self.complex_formation = complex_formation
        
        
    def get_complex_biomass(self, decomposed_complex = None):
        '''Recursive method to get the complex biomass by its individual components'''
    
        if decomposed_complex == None:
            decomposed_complex = COMPLEX(complex_id = 'foo', reaction_id = 'foo',
                                         **{'METABOLITES': self.subcomponents, 
                                           'IDS': [str(i) for i in list(range(len(self.subcomponents)))],
                                           'METABOLITE_TYPES': self.component_types})
            
            
        if 'complex' not in decomposed_complex.component_types:
            total_components = len(decomposed_complex.subcomponents)

            biomass_by_type = dict(zip(decomposed_complex.component_types, [0]* total_components))
            for i in range(total_components): # in case repeated metabolites in complex (homodimers)
                m_ = decomposed_complex.subcomponents[i]
                mt_ = decomposed_complex.component_types[i]
                biomass_by_type[mt_] += get_metabolite_mw(m_)
            
#             if sum(biomass_by_type.values()) != get_metabolite_mw(self): # sanity check
#                 raise ValueError('Something went wrong in decomposing complex biomass')
            
            return biomass_by_type
        else: # unpack nested complexes
            complex_idx = [i for i in range(len(decomposed_complex.component_types)) if decomposed_complex.component_types[i] == 'complex']
            non_complex_idx = sorted(set(list(range(len(decomposed_complex.component_types)))).difference(complex_idx))


            metabolites_ = [decomposed_complex.subcomponents[i] for i in non_complex_idx] + [item for sublist in [decomposed_complex.subcomponents[i].subcomponents for i in complex_idx] for item in sublist]
            ids = [str(i) for i in list(range(len(metabolites_)))]
            metabolite_types = [decomposed_complex.component_types[i] for i in non_complex_idx] + [item for sublist in [decomposed_complex.subcomponents[i].component_types for i in complex_idx] for item in sublist]
            new_complex_info = {'METABOLITES': metabolites_, 'IDS': ids, 'METABOLITE_TYPES': metabolite_types}

            return self.get_complex_biomass(decomposed_complex = COMPLEX(complex_id = 'foo', reaction_id = 'foo', 
                                                                         **new_complex_info))  


# In[ ]:


def get_complex_biomass_change(complex_products, complex_reactants):
    '''Input is two lists of type COMPLEX, one representing those on the product side, one representing those on the reactant side
    output is a dictionary of biomass change for each respective biomass type.'''
    
    product_biomass = dict()
    for cp in complex_products:
        if type(cp)!= COMPLEX:
            raise TypeError('All complex products must be a COMPLEX object')
        for bt, mw in cp.get_complex_biomass().items():
            if bt in product_biomass.keys():
                product_biomass[bt] += mw
            else:
                product_biomass[bt] = mw
    
    reactant_biomass = dict()
    for cr in complex_reactants:
        if type(cr)!= COMPLEX:
            raise TypeError('All complex reactants must be a COMPLEX object')
        for bt, mw in cr.get_complex_biomass().items():
            if bt in reactant_biomass.keys():
                reactant_biomass[bt] += mw
            else:
                reactant_biomass[bt] = mw
    
    for bt in set(product_biomass.keys()).difference(reactant_biomass.keys()):
        reactant_biomass[bt] = 0
    for bt in set(reactant_biomass.keys()).difference(product_biomass.keys()):
        product_biomass[bt] = 0    
    
    return {bt: product_biomass[bt] - reactant_biomass[bt] for bt in product_biomass.keys() if product_biomass[bt] - reactant_biomass[bt] != 0}


# In[ ]:


def parse_me_reaction_id(x):
    if 'HGNC' in x.split('_')[0]:
        return '_'.join(x.split('_')[1:])
    else:
        return x


# In[47]:


# def eval_complex(expr):
#     '''
    
#     Recursive parsing of gprs into lists of complexes. Input expr is a cobra.parse_gpr(gpr_string), 
#     output is a list of netsted lists, each entry of which is a complex joined by 'AND'; netsted lists are joined 
#     with each other by 'OR'
    
#     Inspired by corda source code, should cite them.
    
#     '''
    
#     # corda: https://github.com/resendislab/corda/blob/master/corda/util.py
#     if isinstance(expr, ast.Expression):
#         return eval_complex_recur(expr.body)
#     elif isinstance(expr, ast.Name):
#         return cobra.core.gene.ast2str(expr)
#     elif isinstance(expr, ast.BoolOp):
#         op = expr.op
#         if isinstance(op, ast.Or):
#             return [eval_complex_recur(i) for i in expr.values]
#         elif isinstance(op, ast.And):
#             return [eval_complex_recur(i) for i in expr.values]



def eval_complex_recur_full(expr):
    '''
    
    Recursive parsing of gprs into lists of complexes. Input expr is a cobra.parse_gpr(gpr_string), 
    output is a list of lists, each entry of which is a complex joined by 'AND'; netsted lists are joined 
    with each other by 'OR'. 
    
    Inspired by corda source code, should cite them.
    
    '''
    
    
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
    '''Input is a gene reaction rule, output is a list. If the length of the list is longer than 1, these are
    due to the presence of OR. If the entry itself is a list longer than one, this entry is a complex.'''
    
    complexes = unnested_list(eval_complex_recur_full(parse_gpr(expr)[0]))
    
    # make sure machinery is always output in the same order
    for i in range(len(complexes)):
        if isinstance(complexes[i],list):
            complexes[i] = sorted(complexes[i])
    return complexes


# In[ ]:


def SASA(mw):
    return mw**(0.75)


# In[ ]:


def check_me_mass_balance(r, metabolic_model = params.human_model):
    '''r is a cobra.Reaction object'''
    if len(r.genes) == 0:
        return r.check_mass_balance()
    else:
    
        metabolic_reaction_names = [r.name for r in metabolic_model.reactions if len(r.genes)>0]
        if r.name in metabolic_reaction_names: # metabolic reactions
            # remove coupling constraint
            rxn = {m:c for m,c in r.metabolites.items()  if ('protein' not in m.id) and ('complex' not in m.id)}
        else: # expression reactions
            raise ValueError('Do not currently have code base to get mass balance of expression reactions')
#             if 'TRANSLATION' in r.id:
#                 rxn = {m:c for m,c in r[0].metabolites.items()  if ('protein' not in m.id) and ('complex' not in m.id)}
#                 rxn = {m:c for m,c in rxn.items()  if ('mrna' not in m.id) or ('proxy' in m.id)}
    r_ = cobra.Reaction(' ')
    r_.add_metabolites(rxn) 
    
    return r_.check_mass_balance()

