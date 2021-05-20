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


# In[5]:


def format_final_locations(final_locations, sp, hgnc_id):
    '''
    Parameters
    ----------
    final_locations: list 
        each element is a string representing a compartment in the model
    sp: bool
        whether a signal peptide is present
    hgnc_id: str
        the hgnc id of the associated gene
    
    Returns
    -------
    final_locations_dict: dict
        a dictionary with keys as elements in final_locations and values as means of translation
    '''
    
    # transport rules
    # assume location dictates transport pathway ind of sp;
    # assume all genes are transported to mitochondria 
    # thus, two modes of transport:
    # 1) cytosolic transport: cytosolic translation-->import to final compartment
    # 2) canonical secretion: transport/translation via secretory pathway to final compartment


    # can expand on these based on signal peptide and transmembrane domain logic in the future
    final_locations_dict = {}    
    for loc in final_locations: # no signal peptide consideration
        if loc in ['n', 'c', 'x', 'm', 'i']: 
            # mitochondrial expression not considered
            final_locations_dict[loc] = 'Cytosolic Tranport'
            if sp: 
                warnings.warn(hgnc_id + ': Signal peptides not considered for ' + params.compartments[loc])
        else:
            if not sp:
                # add non-canonical in future

                # current structure assumes signal peptide presence for multi-localizing proteins with atleast
                # one compartment in secretory pathway. in the future, presence of signal peptide could be 
                # conditional for each location, somewhat analogous to transcript isoforms

                warning_ = 'Final location is part of secretory pathway, but no signal peptide indicated.'
                warning_ += 'Non canonical secretion is not considered currently. Changing sp to True'
                warnings.warn(warning_)
                sp = True
            if sp:
                final_locations_dict[loc] = 'Canonical Secretion'
            else:
                final_locations_dict[loc] = 'Non-Canonical Secretion'
    return final_locations_dict


# In[6]:


class gene_information():
    '''This class compiles all the necessary information for a given transcript/protein to be expressed in the 
    ME model. 
    
    Notes: 
    
    1) As of right now, machinery PTMs are not considered. Only non-machinery proteins processed via the 
    secretory pathway can have PTMs. 
    '''
    
    
    def __init__(self, hgnc_id, premrna_seq, mrna_seq, protein_seq,
                 machinery_list = mach.metabolic_machinery,
                 ptms = {}, tmd = 0, sp = False, polyA_length = None, n_exons = None, 
                coupling_params = None):
        '''Initialize object
        
        Parameters
        ----------
        hgnc_id: str
            gene HGNC ID in the format HGNC:#### 
        premrna_seq: str
            the premrna sequence
        mrna_seq: str
            the mrna sequence (length must be <= premrna_seq)
        protein_seq: str
            the protein sequence (length must be <= mrna_seq/3)
        machinery_list: list, default to input M_model genes
            each entry is the HGNC ID of a protein that should be considered as catalyzing an M_model reaction
        ptms: dict, default {}, optional
            keys (str) representing the ptm and values (int) representing the number of that ptms of that kind for 
            that gene. The exception here is gpi, which is binary with 0 for no GPI Anchor and 1 indicating GPI 
            Anchor presence. Key options include ['dsb', 'gpi', 'og'] for 
            ['disulfide bond formation', 'GPI Anchor','O-linked glycosylation'] respectively. PTMs are not currently considered for machinery.  - optional
        tmd: int, default 0, optional
            the number of transmembrane domains the protein has. This is only relevant for proteins processed into 
            the secretory pathway.
        sp: bool, default False, optional
            whether a protein has a signal peptide for the secretory pathway.
            Current format disregards this input, automatically defaulting to True for secretory pathway proteins. 
            Will be used in future for non-canonical secretion. 
        polyA_length: float, default estimated from statistical model
            length of mrna polyA tail
        n_exons: int, default estimated from premrna_seq length
            the number of exons in the mrna isoform. Used to estimate the number of introns (n_exons -1) in the premrna
        coupling_params: dict
            keys (str) specify the parameter, values the value for the parameter. Used to calculate the coupling 
            coefficients. The key-value pairs are as follows:
                a) 'alpha_m': The mrna first-order degradation constant (hrs^-1). If not provided, defaults to a gene-specific value from build/Gregersen_mrna_turnover_processed.tsv if hgnc_id is in the dataframe, otherwise to 0.061
                b) 'alpha_p': The protein first-order degradation constant (hrs^-1). If not provided, defaults to a gene-specific value from build/protein_turnover.csv if hgnc_id is in the dataframe, otherwise to 0.018
                c) 'ptr': float or str
                        the protein-to-RNA ratio. default values drawn from: build/PTR_Gagneur_processed.tsv
                        if float, assumes it is a specified ptr value
                        if str, assumes it is a tissue (must be one of tissues in the .tsv column, or "Median"); 
                        will take the gene-specific tissue-median if the HGNC ID is present in the dataframe, otherwise the median across all genes in the tissue
                        if None/nan, will take the gene-specific median across all tissues if HGNC ID is present in the dataframe, otherwise the whole dataframe median
                        
                        str tissue options: ['Median', Adrenal', 'Appendices', 'Brain', 'Colon', 'Duodenum', 'Endometrium', 'Esophagus', 'Fallopiantube', 'Fat', 'Gallbladder', 'Heart', 
                                        'Kidney', 'Liver', 'Lung', 'Lymphnode', 'Ovary', 'Pancreas', 'Placenta', 'Prostate', 'Rectum', 'Salivarygland', 'Smallintestine', 'Smoothmuscle', 
                                        'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Tonsil', 'Urinarybladder']
              
        '''
        
        self.hgnc_id = hgnc_id
        
        # current structure assumes that a protein is either machinery (catalyzing a reaction) or
        # a secreted protein (processed through secretory pathway, does not catalyze reaction) but not both
        
        if hgnc_id in machinery_list: 
            self.machinery = True
        else:
            self.machinery = False
        
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
            if n_exons != None and n_exons > 1:
                warning_ = 'Premrna and mrna sequences are the same length, but you have indicated this gene' 
                warning_ += 'has atleast 1 intron (n_exons > 1, n_introns = n_exons - 1).' 
                warning_ += 'Setting n_exons to 1 for premrna and mrna of the same length.'
                warnings.warn(warning_)
            n_exons = 1  # n_introns is 0 when the sequence lengths are equal         
            
        if len(mrna_seq) < len(protein_seq)*3:
            raise ValueError(self.hgnc_id + ': The mrna and protein sequence lengths are inconsistent')
        
        # complete checks, assign attributes
        self.premrna_seq = Seq(premrna_seq, generic_rna)
        self.mrna_seq = Seq(mrna_seq, generic_rna)
        self.protein_seq = protein_seq
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
        else: 
            self.sp = bool(sp)
        
        if polyA_length == None or pd.isna(polyA_length):
            self.polyA_length = None
        elif polyA_length >= 0:
            self.polyA_length = polyA_length # can be floating point, rounded in polyA_statistics script
        else:
            raise ValueError(self.hgnc_id + ': polyA_length must either be an integer >= 0 or None/nan')
        
        if n_exons == None or pd.isna(n_exons): 
            self.n_introns = len(self.premrna_seq) * params.rate_intron
        elif n_exons >= 1:
            self.n_introns = n_exons - 1
        else:
            raise TypeError(self.hgnc_id + ': n_exons must either be an integer >= 1 or None/nan')

        # check for mismatch in premrna and mrna seq for sequences with very similar lengths
        valid_seq = True
        for nt in ['A', 'U', 'G', 'C']:
            if self.premrna_seq.count(nt) - self.mrna_seq.count(nt) < 0:
                warnings.warn('Mismatch between premrna and mrna seq, settting both to mrna seq')
                valid_seq = False
                break
        if not valid_seq:
            if self.n_introns > 1:
                raise ValueError('Internal: Unexpected behavior')
            self.premrna_seq = self.mrna_seq
            self.n_introns = 0


        self.all_locations = None
        
        
        # coupling parameters
        if type(coupling_params) != dict:
            if not(coupling_params == None or pd.isna(coupling_params)):
                raise ValueError('Provided coupling parameters are not formatted correctly')
            self.coupling_params = {'alpha_m': None, 'alpha_p': None, 'ptr': None}  
        else: 
            self.coupling_params = coupling_params
        self.get_coupling()
    
    def get_coupling(self):
        # ptr-----    
        if 'ptr' not in self.coupling_params:
            self.coupling_params['ptr'] = None    
        if type(self.coupling_params['ptr']) == float or type(self.coupling_params['ptr']) == int: # specified ptr
            if not (self.coupling_params['ptr'] > 0): #float('nan') will return True as desired
                warnings.warn('Invalid ptr value provided, will use default instead')
                self.coupling_params['ptr'] = None
        elif type(self.coupling_params['ptr']) == str: # unspecified ptr, but tissue specified
            if self.coupling_params['ptr'] not in params.ptr.columns.tolist():
                mssg = 'The specified tissue, type' + self.coupling_params['ptr'] + ' is not available'
                mssg += 'Will use default values instead'
                warnings.warn(mssg)
                self.coupling_params['ptr'] = None # will forward to next statement
            else: # tissue 
                val = float('nan')
                if self.hgnc_id in params.ptr.HGNC_ID.tolist():
                    val = params.ptr.groupby(params.ptr.HGNC_ID).median().loc[self.hgnc_id, self.coupling_params['ptr']]

                if not pd.isna(val): # the if statement above can still product a nan
                    self.coupling_params['ptr'] = val
                else:
                    self.coupling_params['ptr'] = params.ptr.groupby(params.ptr.HGNC_ID).median().loc[:, self.coupling_params['ptr']].median()
        if self.coupling_params['ptr'] is None or pd.isna(self.coupling_params['ptr']): # totally unspecified ptr
            if self.hgnc_id in params.ptr.HGNC_ID.tolist():
                self.coupling_params['ptr'] = params.ptr.groupby(params.ptr.HGNC_ID).median().loc[self.hgnc_id, :].median()
            else:
                self.coupling_params['ptr'] = params.ptr_median
        #alpha_p and # alpha_m-----------
        for tp in ['alpha_p', 'alpha_m']: 
            if tp not in self.coupling_params:
                self.coupling_params[tp] = None
            if type(self.coupling_params[tp]) == float or type(self.coupling_params[tp]) == int: # specified params.turnover[tp]
                if not (self.coupling_params[tp] > 0): #float('nan') will return True as desired
                    warnings.warn('Invalid ' + tp + ' value provided, will use default instead')
                    self.coupling_params[tp] = None
            if self.coupling_params[tp] is None or pd.isna(self.coupling_params[tp]): # totally unspecified ptr
                if self.hgnc_id in params.turnover[tp].index.tolist():
                    self.coupling_params[tp] = params.turnover[tp][self.hgnc_id]
                else:
                    self.coupling_params[tp] = params.turnover[tp + '_median']

        self.coupling = dict()
        denom = (self.coupling_params['alpha_p'] + params.mu)*self.coupling_params['ptr']
        self.coupling['mrna_degradation'] = self.coupling_params['alpha_m']/denom
        self.coupling['mrna_formation'] = ((self.coupling_params['alpha_m']) + params.mu)/denom
#         self.coupling['mrna_dilution'] = params.mu/denom

        for k,v in self.coupling.items():
            if v.subs(params.mu, 1) <= 0:
                raise ValueError('The coupling constraint "' + k + '" must be positive for gene ' + self.hgnc_id)

    def get_final_locations(self, reactions = None, nonmachinery_locations = []):
        '''Assigns a set of final compartments for the protein. For machinery, extracts this from the input
        cobrapy model. This method helps define necessary transport reactions.
        
        Parameters
        ----------
        reactions: list
            each element is a cobra.core.reactions.Reaction associated with the gene
        nonmachinery_locations: list
            a list of strings of final compartments for non-machinery
        
        Returns
        -------
        self.machinery_locations: dictionary
            keys as the final compartments for enzymes and values as the method of 
            synthesis (Cytosolic Transport, Mitochondrial Expression - unimplemented, Canonical Secretion, Non-Canonical Secretion) 
            depending on Boolean rules.
        self.nonmachinery_locations: dictionary
            same as self.machinery, but for non-enzyme proteins
        self.all_locations: dictionary
            combined machinery and nonmachinery locations
        
        '''
        self.machinery_locations = list()
        if self.machinery:
            if reactions is None:
                raise ValueError('For machinery, need associated reactions')
#             rxns = list(metabolic_model.genes.get_by_id(self.hgnc_id).reactions)
            self.machinery_locations = sorted(set([func.get_reaction_compartment(r) for r in reactions])) 
   
        if not self.machinery and len(nonmachinery_locations) == 0:
                raise ValueError(self.hgnc_id + ': For non-machinery, must specify the final compartments')
        
        # no overlap in compartments of machinery and non-machinery
        self.nonmachinery_locations = sorted(set(nonmachinery_locations).difference(self.machinery_locations))

        if len(set(self.nonmachinery_locations).difference(params.compartments.keys())) > 0:
            error = 'At least one of the locations specified is not allowed in this model.'
            raise ValueError(error + ' Allowable comparments include: ' + ', '.join(list(params.compartments.keys())))
   
        self.machinery_locations = format_final_locations(final_locations = self.machinery_locations, sp = self.sp, hgnc_id = self.hgnc_id)
        self.nonmachinery_locations = format_final_locations(final_locations = self.nonmachinery_locations, sp = self.sp, hgnc_id = self.hgnc_id)
        
        # scale coupling
        self.all_locations = self.machinery_locations.copy()
        for k,v in self.nonmachinery_locations.items():
            self.all_locations[k] = v
        # in the case that protein synthesis flux spread across multiple reactions due to multi-localization
        if len(set(self.all_locations.values())) > 1:
            if len(set(self.all_locations.values())) == 2:
                self.coupling['mrna_degradation'] = 0.5*self.coupling['mrna_degradation']
                self.coupling['mrna_formation'] = 0.5*self.coupling['mrna_formation']
            else:
                raise ValueError('Have not yet accounted for Non-Canonical Secretion or other synthesis forms in coupling of mrna degradataion to protein synthesis')

    def check_gene_information(self):
        if self.all_locations == None:
            raise ValueError(self.hgnc_id + ': Must specify a final location for the gene. Use the get_final_locations() method')
        if len(self.ptms) > 0:
            if self.machinery:
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


# In[10]:


ptm_cols = ['DSB', 'GPI', 'NG', 'OG']
ptm_keys = list(params.allowed_ptms.keys())
cp_keys = ['alpha_m', 'alpha_p', 'ptr']

def generate(hgnc_id, psim = params.psim_me, machinery_list = mach.metabolic_machinery, 
                             reactions = None, nonmachinery_locations = []):
    '''Generates gene information object from PSIM. Assumes the gene information object being 
    generated is not for a non-machinery protein.'''
    
    idx = psim[psim.HGNC_ID == hgnc_id].index.tolist()

    entries = psim.loc[idx[0],:]
    if type(entries['LOCATION']) == str:
        entries['LOCATION'] = list(entries['LOCATION'].split(']')[0].split('[')[1].split(','))
    
    cp_values = entries['ALPHA_M'], entries['ALPHA_P'], entries['PTR']

    gene_info = gene_information(hgnc_id = entries['HGNC_ID'], 
                    premrna_seq = entries['PREMRNA_SEQ'], mrna_seq = entries['MRNA_SEQ'], 
                    protein_seq = entries['PROTEIN_SEQ'], 
                    machinery_list = machinery_list,
                    ptms = dict(zip(['dsb', 'og', 'gpi'],[entries['DSB'], entries['OG'], entries['GPI']])),
                    tmd = entries['TMD'], sp = entries['SP'], polyA_length = entries['POLYA_LENGTH'], 
                    n_exons = entries['N_EXONS'], 
                    coupling_params = dict(zip(cp_keys, cp_values)))
    gene_info.get_final_locations(reactions = reactions, 
                                  nonmachinery_locations = nonmachinery_locations)
    gene_info.check_gene_information()
    return gene_info

