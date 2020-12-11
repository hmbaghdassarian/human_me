#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import random

import cobra
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_rna

import requests, sys, json, re, warnings

import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils import machinery as mach
from utils import parameters as params
from utils import functions as func


# In[30]:


class gene_information():
    '''This class compiles all the necessary information for a given protein to be expressed in the 
    ME model. 
    
    Notes: 
    
    1) As of right now, machinery PTMs are not considered. Proteins processed through the secretory pathway includes 
    those with a final location in the following compartments: ['l', 'r', 'e', 'x', 'g', 'pm']. See get_final_locations() method
    for details.
    
    2) Expression machinery is checked from the utils expression model; there is no input for this. Further, 
    the expression model may include expression machinery that the final specific ME model does not include.'''
    
    
    def __init__(self, hgnc_id, premrna_seq, mrna_seq, protein_seq,
                 machinery_list = mach.metabolic_machinery,
                 ptms = {}, tmd = 0, sp = False, polyA_length = None, n_introns = None, 
                coupling_params = None):
        '''
        
        1) HGNC ID is a string in the format HGNC:#### - required.
        
        2-4) Relevant string representing sequence - required
        
        5) Metabolic machinery is a list of HGNC IDs of metabolic enzymes in the metabolic model. 
        The default is a list generated from the input CobraPy model in utils.
        
        6) PTMs is a dictionary with keys as a string representing the ptm and values as an integer
        representing the number of that ptms of that kind for that gene. The exception here is gpi, which is binary 
        with 0 for no GPI Anchor and 1 indicating GPI Anchor presence. The keys of the dictionary utils.allowed_ptms 
        show all possible key values here. PTMs are not currently considered for machinery.  - optional
        
        7) TMD is an integer indicating the number of transmembrane domains the protein has. This is only relevant
        for proteins processed into secretory pathway. - optional
        
        8) SP is a boolean indicating whether a protein has a signal peptide. 
        Not used in current format (automatically defaults to True for secretory pathway proteins) - unimplemented
                
        9) polyA_length is an floating point representing the length of the polyA tail. This information will be 
        estimated by a statistical model if not provided. - optional
        
        10) n_introns is an integer representing the length of the polyA tail. This information will be estimated
        if not provided. Should be specific to the transcript isoform. - optional
        
        11) coupling_params is a dictionary with required parameters for coupling constraints. The key-value pairs are as follows:
            a) 'mrna_half_life': The half life for the mrna in units of hours. If not provided, defaults to 10.
            b) 'alpha_p': The protein first-order degradation constant in units of hours^-1. If not provided, defaults to 0.02. 
            c) 'ptr': A value for the protein-to-RNA ratio
            d) 'ptr_tissue': A tissue from which to get or estimate the PTR, if a ptr value is not provided. Options include:
            ['Median', Adrenal', 'Appendices', 'Brain', 'Colon', 'Duodenum', 'Endometrium', 'Esophagus', 'Fallopiantube', 'Fat', 'Gallbladder', 'Heart', 
            'Kidney', 'Liver', 'Lung', 'Lymphnode', 'Ovary', 'Pancreas', 'Placenta', 'Prostate', 'Rectum', 'Salivarygland', 'Smallintestine', 'Smoothmuscle', 
            'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Tonsil', 'Urinarybladder']
            e) 'constant_ptr': bool - whether to estimate the same PTR for all genes (True). Not recommended. Will override provided ptr/ptr_tissue
        '''
        
        self.hgnc_id = hgnc_id
        
        # current structure assumes that a protein is either machinery (catalyzing a reaction) or
        # a secreted protein (processed through secretory pathway, does not catalyze reaction) but not both
        
            
        #self.module = list()
        if hgnc_id in machinery_list: #or hgnc_id in expression_machinery:
            self.module = 'Machinery'
#             if hgnc_id in machinery_list:
#                 self.module += ['Metabolic Machinery']
#             else:
#                 self.module += ['Expression Machinery']
        else:
            self.module = 'Non-Machinery'
        
        # sequence check
        if premrna_seq == None or mrna_seq == None or protein_seq == None:
            raise ValueError(self.hgnc_id + ': All of the sequence types (premrna, mrna, protein) must be provided')
        if 'N' in mrna_seq:
            warnings.warn(self.hgnc_id + ': The letter N is in the mrna sequence. Replacing with a random nucleotide')
            mrna_seq = mrna_seq.replace('N', 'U')#mrna_seq.replace('N', random.choice(['A', 'U', 'G', 'C']))
        if 'N' in premrna_seq:
            warnings.warn(self.hgnc_id + ': The letter N is in the premrna sequence. Replacing with a random nucleotide')
            premrna_seq = premrna_seq.replace('N', 'C')#('N', random.choice(['A', 'U', 'G', 'C']))
        if len(set(premrna_seq).difference(['A', 'U', 'G', 'C'])) > 0:
            raise ValueError(self.hgnc_id + ': The premrna sequence contains bases which are not allowed')
        if len(set(mrna_seq).difference(['A', 'U', 'G', 'C'])) > 0:
            raise ValueError(self.hgnc_id + ': The mrna sequence contains bases which are not allowed')
        
        if 'X' in protein_seq:
            warnings.warn(self.hgnc_id + ': The letter X is in the protein sequence. Replacing with a random amino acid')
            protein_seq = protein_seq.replace('X', 'A')#('X', random.choice(params.amino_acids))
        if 'U' in protein_seq:
            warnings.warn(self.hgnc_id + ': Selenocysteine not currently considered by model, replacing with cysteine')
            protein_seq = protein_seq.replace('U', 'C')
        
        if len(set(protein_seq).difference(params.amino_acids)) > 0:
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
#         self.premrna_mass = calculate_molecular_weight(seq = self.premrna_seq)/1000 #kDa
        self.mrna_seq = Seq(mrna_seq, generic_rna)
        self.protein_seq = protein_seq
        
#         self.protein_mass = calculate_molecular_weight(seq=self.protein_seq, seq_type='protein')/1000 #kDa
        self.L_protein = len(self.protein_seq)
        self.amino_acid_counts = {k: self.protein_seq.count(k) for k in params.amino_acids}
        
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
            raise TypeError(self.hgnc_id + ': SP must be a boolean')
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
            raise TypeError(self.hgnc_id + ': n_introns must either be an integer >= 0 or None/nan')
        
#         self.keff = keff
#         if self.keff == None and self.module == 'Machinery':
#             warnings.warn(self.hgnc_id + ': No keff specified for this enzyme, will assume a value in model building')
#         elif self.keff != None and self.module == 'Non-Machinery':
#             warnings.warn(self.hgnc_id + ': keff specified for a non-machinery protein, will not be used')
            
            
        self.final_locations = None
        
        
        # coupling parameters
        ptr_tissue_orig = None
        if coupling_params == None or pd.isna(coupling_params):
            coupling_params = params.coupling_params
        else:
            if 'ptr_tissue' in coupling_params.keys():
                ptr_tissue_orig = coupling_params['ptr_tissue']
            for k in params.coupling_params.keys():
                if (k not in coupling_params.keys()) or (coupling_params[k] == None) or (pd.isna(coupling_params[k])):
                    coupling_params[k] = params.coupling_params[k]
                    
        # get the PTR
        if not coupling_params['constant_ptr']:
            if not (coupling_params['ptr'] is None or pd.isna(coupling_params['ptr'])):
                self.ptr = coupling_params['ptr']
                if not (ptr_tissue_orig is None or pd.isna(ptr_tissue_orig)):
                    warnings.warn('You have indicated using a specific PTR, ignoring user input PTR tissue')
                
            else:
                if coupling_params['ptr_tissue'] not in params.ptr.columns:
                    raise ValueError('Specified tissue type for PTR is not available')
                else: # tissue estimation
                    if self.hgnc_id in params.ptr[coupling_params['ptr_tissue']].dropna().index:
                        self.ptr = params.ptr.loc[self.hgnc_id, coupling_params['ptr_tissue']]
                    else:
                        self.ptr = params.ptr[coupling_params['ptr_tissue']].median()
        else:
            self.ptr = params.constant_ptr
            if not (coupling_params['ptr'] is None or pd.isna(coupling_params['ptr'])):
                warnings.warn('You have indicated using a constant PTR, ignoring user input PTR')
            if not (ptr_tissue_orig is None or pd.isna(ptr_tissue_orig)):
                warnings.warn('You have indicated using a constant PTR, ignoring user input PTR tissue')
        
        self.coupling = dict()
        self.coupling['c2'] = (np.log(2)/coupling_params['mrna_half_life'])/((coupling_params['alpha_p'] + params.mu)*self.ptr)
        # c1c
#         self.coupling['c1'] = params.mu/((coupling_params['alpha_p'] + params.mu)*self.ptr)
        # c1b
        self.coupling['c1'] = ((np.log(2)/coupling_params['mrna_half_life']) + params.mu)/((coupling_params['alpha_p'] + params.mu)*self.ptr)

       
    def get_final_locations(self, metabolic_model = params.human_model, final_locations = None):
        '''Assigns a set of final compartments for the protein. For machinery, extracts this from the inputer
        cobrapy model. For non-machinery, final_locations should be specified by a list of strings
        within the allowable compartments. This method helps define necessary transport reactions.
        
        The final output will be a dictionary with keys as the final locations and values as the method of 
        synthesis (Cytosolic Transport, Mitochondrial Expression - unimplemented, Canonical Secretion, Non-Canonical Secretion) 
        depending on Boolean rules.'''
        
        if self.module == 'Machinery':
            if final_locations != None:
                warnings.warn(self.hgnc_id + ': Final location extacted from cobrapy model, will disregard user input.')
              
            rxns = list(metabolic_model.genes.get_by_id(self.hgnc_id).reactions)
            final_locations = sorted(set([func.get_reaction_compartment(r) for r in rxns])) # redundancy from multiple reactions

        elif self.module == 'Non-Machinery':
    
            if final_locations == None:
                raise ValueError(self.hgnc_id + ': For non-machinery, must specify the final locations')
            if type(final_locations) != list:
                raise ValueError(self.hgnc_id + ': Final locations must be a list of string')
            if len(set(final_locations).difference(params.compartments.keys())) > 0:
                error = 'At least one of the locations specified is not allowed in this model.'
                raise ValueError(error + ' Allowable comparments include: ' + ', '.join(list(params.compartments.keys())))
        else:
            raise ValueError('Model does not currently deal with both non-machinery and machinery')
   
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
                    warnings.warn(self.hgnc_id + ': Signal peptides not considered for ' + params.compartments[loc])
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
        
        # OLD
        # in the case that protein synthesis flux spread across multiple reactions due to multi-localization
        if len(set(self.final_locations.values())) > 1:
            if len(set(self.final_locations.values())) == 2:
                self.coupling['c2'] = 0.5*self.coupling['c2']
                self.coupling['c1'] = 0.5*self.coupling['c1']
            else:
                raise ValueError('Have not yet accounted for Non-Canonical Secretion or other synthesis forms in coupling of mrna degradataion to protein synthesis')

    def check_gene_information(self):
        if self.final_locations == None:
            raise ValueError(self.hgnc_id + ': Must specify a final location for the gene. Use the get_final_locations() method')
        if len(self.ptms) > 0:
            if self.module == 'Machinery':
                # change in the future
                warnings.warn(self.hgnc_id + ': PTMs are not considered for machinery proteins currently')
                self.ptms = {}
            if len(set(self.ptms.keys()).difference(params.allowed_ptms.keys())) > 0:
                warnings.warn(self.hgnc_id + ': Atleast one of the PTMs provided will not be considered in this model')
                self.ptms = {k:v for k in self.ptms.keys() if k in params.allowed_ptms.keys()}
            if 'gpi' in self.ptms.keys() and self.ptms['gpi'] > 1:
                warnings.warn(self.hgnc_id + ': GPI is binary, 1 for presence or 0 for absence. Changing to 1')
                self.ptms['gpi'] = 1


        print('No errors raised')


# In[ ]:


ptm_cols = ['DSB', 'GPI', 'NG', 'OG']
ptm_keys = list(params.allowed_ptms.keys())
cp_keys = ['mrna_half_life', 'alpha_p', 'ptr', 'ptr_tissue', 'constant_ptr']

def generate(hgnc_id, psim = params.psim_me, machinery_list = mach.metabolic_machinery, 
                             metabolic_model = params.human_model):
    '''Generates gene information object from PSIM'''
    
    idx = psim[psim.HGNC_ID == hgnc_id].index.tolist()
    if len(idx) == 0:
        raise ValueError(hgnc_id + ' is not in the PSIM')
    if len(idx) > 1:
        warnings.warn('More than one entry of this gene by HGNC ID in PSIM, taking the first')

    entries = psim.loc[idx[0],:]
    if type(entries['LOCATION']) == str:
        entries['LOCATION'] = list(entries['LOCATION'].split(']')[0].split('[')[1].split(','))
    
    cp_values = entries['MRNA_HALF_LIFE'], entries['ALPHA_P'], entries['PTR'], entries['PTR_TISSUE'], entries['CONSTANT_PTR']

    gene_info = gene_information(hgnc_id = entries['HGNC_ID'], 
                    premrna_seq = entries['PREMRNA_SEQ'], mrna_seq = entries['MRNA_SEQ'], 
                    protein_seq = entries['PROTEIN_SEQ'], 
                    machinery_list = machinery_list,
                    ptms = dict(zip(['dsb', 'og', 'gpi'],[entries['DSB'], entries['OG'], entries['GPI']])),
                    tmd = entries['TMD'], sp = entries['SP'], polyA_length = entries['POLYA_LENGTH'], 
                    n_introns = entries['N_INTRONS'], 
                    coupling_params = dict(zip(cp_keys, cp_values)))
    gene_info.get_final_locations(metabolic_model = metabolic_model, 
                                  final_locations = entries['LOCATION'])
    gene_info.check_gene_information()
    return gene_info

