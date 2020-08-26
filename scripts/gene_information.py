#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import random

import cobra
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_rna
from Bio.SeqUtils import molecular_weight as calculate_molecular_weight

import requests, sys, json, re, warnings

import sys
sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *
from utils import *


# In[8]:


class gene_information():
    '''This class compiles all the necessary information for a given protein to be expressed in the 
    ME model. 
    
    Notes: As of right now, machinery PTMs are not considered. Proteins processed through the secretory pathway includes 
    1) anything assigned to the Secreted Protein module, or 2) machinery with a final location
    in the following compartments: ['l', 'r', 'e', 'x', 'g', 'pm']. See get_final_locations() method
    for details.'''
    
    
    def __init__(self, hgnc_id, premrna_seq, mrna_seq, protein_seq,
                 metabolic_machinery = metabolic_machinery,
                 ptms = {}, tmd = 0, sp = False, keff = None, polyA_length = None, n_introns = None):
        '''
        
        1) Metabolic model is a cobrapy model - required. 
        
        2) HGNC ID is a string in the format HGNC:#### - required.
        
        3-5) Relevant string representing sequence - required
        
        6) Metabolic machinery is a list of HGNC IDs of metabolic enzymes in the metabolic model. 
        The defaul is a list generated from the input CobraPy model in utils.
        
        7) PTMs is a dictionary with keys as a string representing the ptm and values as an integer
        representing the number of that ptms of that kind for that gene. The exception here is gpi, which is binary 
        with 0 for no GPI Anchor and 1 indicating GPI Anchor presence. The keys of the dictionary allowed_ptms 
        show all possible key values here. PTMs are not currently considered for machinery.  - optional
        
        8) TMD is an integer indicating the number of transmembrane domains the protein has. This is only relevant
        for proteins processed into secretory pathway. - optional
        
        9) SP is a boolean indicating whether a protein has a signal peptide. 
        Not used in current format - unimplemented
        
        10) keff is a float representing the kinetic constant the enzyme in [units]. - optional
        
        11) polyA_length is an floating point representing the length of the polyA tail. This information will be 
        estimated by a statistical model if not provided. - optional
        
        12) n_introns is an integer representing the length of the polyA tail. This information will be estimated
        if not provided. Should be specific to the transcript isoform. - optional
        
        '''
        
        self.hgnc_id = hgnc_id
        
        # current structure assumes that a protein is either machinery (catalyzing a reaction) or
        # a secreted protein (processed through secretory pathway, does not catalyze reaction) but not both
        
#         if expression_model != None:
#             expression_machinery = [g.id for g in expression_model.genes]
#             if 'ribosome' in expression_machinery:
#                 expression_machinery.remove('ribosome')
#         else:
#             expression_machinery = list()
            
        self.module = list()
        if hgnc_id in metabolic_machinery or hgnc_id in expression_machinery:
            if hgnc_id in metabolic_machinery:
                self.module += ['Metabolic Machinery']
            else:
                self.module += ['Expression Machinery']
        else:
            self.module += ['Non-Machinery']
        
        # sequence check
        if premrna_seq == None or mrna_seq == None or protein_seq == None:
            raise ValueError(self.hgnc_id + ': All of the sequence types (premrna, mrna, protein) must be provided')
        if 'N' in mrna_seq:
            warnings.warn(self.hgnc_id + ': The letter N is in the mrna sequence. Replacing with a random nucleotide')
            mrna_seq = mrna_seq.replace('N', random.choice(['A', 'U', 'G', 'C']))
        if 'N' in premrna_seq:
            warnings.warn(self.hgnc_id + ': The letter N is in the premrna sequence. Replacing with a random nucleotide')
            premrna_seq = premrna_seq.replace('N', random.choice(['A', 'U', 'G', 'C']))
        if len(set(premrna_seq).difference(['A', 'U', 'G', 'C'])) > 0:
            raise ValueError(self.hgnc_id + ': The premrna sequence contains bases which are not allowed')
        if len(set(mrna_seq).difference(['A', 'U', 'G', 'C'])) > 0:
            raise ValueError(self.hgnc_id + ': The mrna sequence contains bases which are not allowed')
        
        if 'X' in protein_seq:
            warnings.warn(self.hgnc_id + ': The letter X is in the protein sequence. Replacing with a random amino acid')
            protein_seq = protein_seq.replace('X', random.choice(amino_acids))
        if 'U' in protein_seq:
            warnings.warn(self.hgnc_id + ': Selenocysteine not currently considered by model, replacing with cysteine')
            protein_seq = protein_seq.replace('U', 'C')
        
        if len(set(protein_seq).difference(amino_acids)) > 0:
            raise ValueError(self.hgnc_id + ': The protein sequence contains amino acids which are not allowed')
        
        if len(premrna_seq) < len(mrna_seq):
            raise ValueError(self.hgnc_id + ': The premrna sequence provided is shorter than the mrna sequence provided')
        elif len(premrna_seq) == len(mrna_seq):
            if premrna_seq != mrna_seq:
                raise ValueError(self.hgnc_id + ': Premrna and mrna sequences are the same length, but not the same sequence')
            if n_introns != None and n_introns > 0:
                warning_ = 'Premrna and mrna sequences are the same length, but you have indicated this is' 
                warning_ += 'not an intronless gene. Setting n_introns to None'
                warnings.warn(warning_)
                n_introns = None
#         else:
#             premrna_base_counts, mrna_base_counts = dict(), dict()
#             for base_letter in seq_element_map.keys():
#                 premrna_base_counts[base_letter] = premrna_seq.count(base_letter)
#                 mrna_base_counts[base_letter] = mrna_seq.count(base_letter)
#             for k,v in mrna_base_counts.items():
#                 if v > premrna_base_counts[k]:
#                     raise ValueError(self.hgnc_id + ': Number of ' + k + ' bases in premrna sequence less than that of mrna sequence')
            
#             self.premrna_base_counts = premrna_base_counts
#             self.mrna_base_counts = mrna_base_counts
            
            
        if len(mrna_seq) < len(protein_seq)*3:
            warnings.warn(self.hgnc_id + ': The mrna and protein sequence lengths are inconsistent')

        self.premrna_seq = Seq(premrna_seq, generic_rna)
        self.mrna_seq = Seq(mrna_seq, generic_rna)
        self.protein_seq = protein_seq
        
        self.protein_mass = calculate_molecular_weight(seq=self.protein_seq, seq_type='protein')
        self.L_protein = len(self.protein_seq)
        self.amino_acid_counts = {k: self.protein_seq.count(k) for k in amino_acids}
        
        remove_ptms = list()
        for k,v in ptms.items():
            if v == None or pd.isna(v) or v == 0:
                remove_ptms.append(k)
        for k in remove_ptms:
            del ptms[k]
        self.ptms = ptms
        
        
        if pd.isna(tmd) or tmd == None or round(tmd) == 0:
            self.tmd = 0
        else:
            self.tmd = round(tmd) # must be an integer
        
        if pd.isna(sp) or sp == None:
            self.sp = False 
        elif type(sp) != bool:
            raise ValueError(self.hgnc_id + ': SP must be a boolean')
        else:
            self.sp = sp
        
        if polyA_length == None or pd.isna(polyA_length):
            self.polyA_length = None
        elif polyA_length >= 0:
            self.polyA_length = polyA_length # can be floating point, rounded in polyA_statistics script
        else:
            raise ValueError(self.hgnc_id + ': polyA_length must either be an integer >= 0 or None/nan')
        
        if n_introns == None or pd.isna(n_introns): # or round(n_introns) == n_introns):
            self.n_introns = None
        elif n_introns >= 0:
            self.n_introns = round(n_introns) # must be an integer
        else:
            raise ValueError(self.hgnc_id + ': n_introns must either be an integer >= 0 or None/nan')
        
#         self.keff = keff
#         if self.keff == None and self.module == 'Machinery':
#             warnings.warn(self.hgnc_id + ': No keff specified for this enzyme, will assume a value in model building')
#         elif self.keff != None and self.module == 'Non-Machinery':
#             warnings.warn(self.hgnc_id + ': keff specified for a non-machinery protein, will not be used')
            
            
        self.final_locations = None
       
    def get_final_locations(self, metabolic_model = human_model, final_locations = None):
        '''Assigns a final compartment for proteins. For machinery, extracts this from the model. 
        For secreted proteins, final_locations should be specified by a list of strings
        within the allowable compartments. This method helps define necessary transport reactions.
        
        The final output will be a dictionary with keys as the final locations and values as the method of 
        synthesis (Cytosolic Transport, Mitochondrial Expression, Canonical Secretion, Non-Canonical Secretion) 
        depending on Boolean rules. Traditional are those that don't go through the secretory pathway.'''
        
        if self.module != 'Non-Machinery':
            if final_locations != None:
                warnings.warn(self.hgnc_id + ': Final location extacted from cobrapy model, will disregard user input.')
  
            rxns = list()
            if 'Metabolic Machinery' in self.module:
                rxns += list(metabolic_model.genes.get_by_id(self.hgnc_id).reactions)
            if 'Expression Machinery' in self.module:
                rxns += list(expression_model_2.genes.get_by_id(self.hgnc_id).reactions)
#             fl = final_locations.copy()
            final_locations = []
            
            for r in rxns:
            # proteins can be associated with multiple locations due to multiple reactions, but for each reaction
            # we want that protein to be associated with one compartment
#                 compartments_ = r.compartments.copy()
#                 if len(compartments_) == 1: # not needed but more efficient
#                     final_locations += list(compartments_)
#                     pass
#                 elif len(compartments_) == 2: # for reactions that occur in more than one compartment
#                     if 'c' in compartments_: # remove cytoplasmic compartment between the two for machinery
#                         compartments_.remove('c')
#                     else: # choose compartment on reactant side if no cytoplasmic compartment
#                         reactant_compartments_ = set([m.compartment for m in r.reactants])
#                         if len(reactant_compartments_) == 1:
#                             compartments_ = reactant_compartments_
#                         else:
#                             compartments_ = max(reactant_compartments_, key = list(reactant_compartments_).count)
#                 elif len(compartments_) > 2: # hardcoded for ASPGLUm reaction
#                     compartments_ = {'i'}
                final_locations += [get_reaction_compartment(r)]
            final_locations = sorted(set(final_locations))# + fl)) # redundancy from multiple reactions

                 
        if 'Non-Machinery' in self.module:
            if final_locations == None:
                raise ValueError(self.hgnc_id + ': For non-machinery, must specify the final locations')
            if type(final_locations) != list:
                raise ValueError(self.hgnc_id + ': Final locations must be a list of string')
            if len(set(final_locations).difference(compartments.keys())) > 0:
                error = 'At least one of the locations specified is not allowed in this model.'
                raise ValueError(error + ' Allowable comparments include: ' + ', '.join(list(compartments.keys())))

   
        # transport rules
        # assume location dictates transport pathway ind of sp;
        # assume all genes are transported to mitochondria 
        # thus, two modes of transport:
        # 1) cytosolic transport: cytosolic translation-->import to final compartment
        # 2) canonical secretion: transport/translation via secretory pathway to final compartment
        
        # can expand on these based on signal peptide and transmembrane domain logic in the future
        self.final_locations = {}    
        for loc in final_locations: # no signal peptide consideration
            if loc in ['n', 'c', 'x', 'm', 'i']: 
                # mitochondrial expression not considered
                self.final_locations[loc] = 'Cytosolic Tranport'
                if self.sp: 
                    warnings.warn(self.hgnc_id + ': Signal peptides not considered for these compartments')
            else:
                if not self.sp:
                    # add non-canonical in future
                    
                    # current structure assumes signal peptide presence for multi-localizing proteins with atleast
                    # one compartment in secretory pathway. in the future, presence of signal peptide could be 
                    # conditional for each location, somewhat analogous to transcript isoforms
                    
                    warning_ = 'Final location is part of secretory pathway, but no signal peptide indicated.'
                    warning_ += 'Non canonical secretion is not considered currently. Changing sp to True'
                    warnings.warn(warning_)
                    self.sp = True
                if self.sp:
                    self.final_locations[loc] = 'Canonical Secretion'
                else:
                    self.final_locations[loc] = 'Non-Canonical Secretion'

    def check_gene_information(self):
        if self.final_locations == None:
            raise ValueError(self.hgnc_id + ': Must specify a final location for the gene. Use the get_final_locations() method')
        if len(self.ptms) > 0:
            if 'Metabolic Machinery' in self.module or 'Expression Machinery' in self.module:
                # change in the future
                warnings.warn(self.hgnc_id + ': PTMs are not considered for machinery proteins currently')
                self.ptms = {}
            elif len(set(self.ptms.keys()).difference(allowed_ptms.keys())) > 0:
                warnings.warn(self.hgnc_id + ': Atleast one of the PTMs provided will not be considered in this model')
                self.ptms = {k:v for k in self.ptms.keys() if k in allowed_ptms.keys()}
            elif 'gpi' in self.ptms.keys() and self.ptms['gpi'] > 1:
                warnings.warn(self.hgnc_id + ': GPI is binary, 1 for presence or 0 for absence. Changing to 1')
                self.ptms['gpi'] = 1


        print('No errors raised')

