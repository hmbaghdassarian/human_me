#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
import pandas as pd

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import molecular_weight as calculate_molecular_weight

import requests, sys, json, re, warnings

import sys
sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *


# In[2]:


# define necessary variables 
human_model_2 = cobra.io.load_json_model(local_data_path + 'raw/RECON3D.json')
compartments = human_model_2.compartments
compartments['pm'] = 'plasma membrane'
allowed_ptms = {'dsb': 'disulfide bond formation', 'gpi': 'GPI Anchor', 
               'ng': 'N-linked glycosylation', 'og': 'O-linked glycosylation'}

amino_acids = ['A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 
              'P', 'S', 'T', 'W', 'Y', 'V']



# In[3]:


def get_premrna_seq(ensg_id):
    try:
        hyperlink = 'https://rest.ensembl.org/sequence/id/' + ensg_id + '?' 
        gene_sequence = requests.get(hyperlink, headers={ "Content-Type" : "text/plain"}).text
        return str(Seq(gene_sequence).transcribe()) # introns and UTRs
    except:
        raise ValueError('This ensg id cannot be queried for a premrna sequence')

def get_mrna_seq(ensg_id):
    try:
        hyperlink = 'https://rest.ensembl.org/sequence/id/' + ensg_id + '?' 
        cdna = requests.get(hyperlink+'type=cdna;multiple_sequences=2', 
                            headers={ "Content-Type" : "text/plain"}).text.splitlines()[0]
        return str(Seq(cdna, generic_dna).transcribe()) # UTRs, no introns, no polyA tail right now 
    except:
        raise ValueError('This ensg id cannot be queried for a mrna sequence')

def get_protein_seq(ensg_id):
    try:
        hyperlink = 'https://rest.ensembl.org/sequence/id/' + ensg_id + '?' 
        return requests.get(hyperlink+'type=protein;multiple_sequences=2', 
                            headers={ "Content-Type" : "text/plain"}).text.splitlines()[0]
    except:
        raise ValueError('This ensg id cannot be queried for a protein sequence')



class gene_information():
    '''This class compiles all the necessary information for a given protein to be expressed in the 
    ME model. 
    
    Notes: As of right now, PTMs are not considered. Proteins processed through the secretory pathway includes 
    1) anything assigned to the Secreted Protein module, 2) machinery with a signal peptide, or 3) machinery with a final location
    in the following compartments: ['l', 'r', 'e', 'x', 'g', 'pm']. See get_final_locations() method
    for details.'''
    
    
    def __init__(self, metabolic_model, hgnc_id, ptms = {}, tmd = 0, sp = False, keff = None):
        '''
        1) Metabolic model is a cobrapy model. 
        2) HGNC ID is a string in the format HGNC:####.
        
        3) PTMs is a dictionary with keys as the ptm and values as the number of that ptm for that gene.
        PTMs are not currently considered . 
               
        4) TMD is an integer indicating the number of transmembrane domains the protein has. This is only relevant
        for proteins processed into secretory pathway. 
        
        5) SP is a boolean indicating whether a protein has a signal peptide.
        
        6) keff is a float representing the kinetic constant the enzyme in [units]'''
        
        self.hgnc_id = hgnc_id
        
        # current structure assumes that a protein is either machinery (catalyzing a reaction) or
        # a secreted protein (processed through secretory pathway, does not catalyze reaction) but not both
        
        machinery = [g.id for g in metabolic_model.genes] # not super efficient to do this each time
        if hgnc_id in machinery:
            self.module = 'Machinery'
        else:
            self.module = 'Non-Machinery'
        
        self.ptms = ptms
        self.tmd = tmd
        self.sp = sp
        self.keff = keff
        
        if self.keff == None and self.module == 'Machinery':
            warnings.warn('No keff specified for this enzyme, will assume a value in model building')
        
        self.final_locations = None
        self.premrna_seq = None
        self.mrna_seq = None
        self.protein_seq = None
        self.protein_mass = None
       
    def get_final_locations(self, metabolic_model, final_locations = None):
        '''Assigns a final compartment for proteins. For machinery, extracts this from the model. 
        For secreted proteins, final_locations should be specified by a list of strings
        within the allowable compartments. This method helps define necessary transport reactions.
        
        The final output will be a dictionary with keys as the final locations and values as the method of 
        synthesis (Traditional Expression, Mitochondrial Expression, Canonical Secretion, Non-Canonical Secretion) 
        depending on Boolean rules. Traditional are those that don't go through the secretory pathway.'''
        
        if self.module == 'Machinery':
            if final_locations != None:
                warnings.warn('Final location extacted from cobrapy model, will disregard user input.')

            rxns = list(metabolic_model.genes.get_by_id(self.hgnc_id).reactions)
            final_locations = []
            for r in rxns:
                final_locations += list(r.compartments)
            final_locations = sorted(set(final_locations))

                 
        if self.module == 'Non-Machinery':
            if final_locations == None:
                raise ValueError('For non-machinery, must specify the final locations')
            if type(final_locations) != list:
                raise ValueError('Final locations must be a list of string')
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
            else:
                self.final_locations[loc] = 'Canonical Secretion'
                if not self.sp:
                    # add non-canonical in future
                    
                    # current structure assumes signal peptide presence for multi-localizing proteins with atleast
                    # one compartment in secretory pathway. in the future, presence of signal peptide could be 
                    # conditional for each location, somewhat analogous to transcript isoforms
                    
                    warning_ = 'Final location is part of secretory pathway, but no signal peptide indicated.'
                    warning_ += 'Non canonical secretion is not considered currently. Changing sp to True'
                    warnings.warn(warning_)
                    self.sp = True

    
    # could potentially add translation/transcription of sequences in the future
    
    def get_sequences(self, ensg_id = None, premrna_seq = None, mrna_seq = None, protein_seq = None):
        '''For user provided information, must input strings for each of the sequence types 
        (premrna, mrna, and protein). Otherwise, must provide the ensembl gene id and this method will use the 
        Ensembl REST API to get the necessary information.'''
        
        # premrna is expected to include UTRs and introns, mrna to include UTRs
        # since no transcription and translation reactions are being encoded, this is not strict
        
        if premrna_seq == None and mrna_seq == None and protein_seq == None:
            if ensg_id == None:
                raise ValueError('Must provide either gene sequences or ensg id')
            else:
                self.premrna_seq, self.mrna_seq, self.protein_seq = get_premrna_seq(ensg_id), get_mrna_seq(ensg_id), get_protein_seq(ensg_id)
        else:
            if premrna_seq == None or mrna_seq == None or protein_seq == None:
                raise ValueError('All of the sequence types (premrna, mrna, protein) must be provided')
            elif len(set(premrna_seq).difference(['A', 'U', 'G', 'C'])) > 0:
                raise ValueError('The premrna sequence contains bases which are not allowed')
            elif len(set(mrna_seq).difference(['A', 'U', 'G', 'C'])) > 0:
                raise ValueError('The mrna sequence contains bases which are not allowed')
            elif len(set(protein_seq).difference(amino_acids)) > 0:
                raise ValueError('The protein sequence contains amino acids which are not allowed')
            else:
                self.premrna_seq = premrna_seq
                self.mrna_seq = mrna_seq
                self.protein_seq = protein_seq
        
        self.protein_mass = calculate_molecular_weight(seq=self.protein_seq, seq_type='protein')

    def check_gene_information(self):
        if self.final_locations == None:
            raise ValueError('Must specify a final location for the gene. Use the get_final_locations() method')
        if self.premrna_seq == None or self.mrna_seq == None or self.protein_seq == None or self.protein_mass == None:
            raise ValueError('Must specify transcriptional and translational sequences. Use the get_sequences() method')
        
        
        if len(self.ptms) > 0:
            if self.module == 'Machinery':
                # change in the future
                warnings.warn('PTMs are not considered for machinery proteins currently')
            elif len(set(self.ptms.keys()).difference(allowed_ptms.keys())) > 0:
                warnings.warn('Atleast one of the PTMs provided will not be considered in this model')
            

        print('No errors raised')


# # Usage

# # In[4]:


# psim_me = pd.read_csv(local_data_path + 'processed/psim_me.csv', index_col = 0)
# human_model = cobra.io.load_json_model(local_data_path + 'processed/corrected_recon2_2.json')
# sp_dict = {1: True, 0: False}
# ptm_cols = ['DSB', 'GPI', 'NG', 'OG']
# ptm_keys = list(allowed_ptms.keys())


# # In[5]:


# psim_me.head()


# # Usage example 

# # In[6]:


# # gene catalyzing metabolic reaction, not processed via secretory pathway
# gene1_id = human_model.genes[1].id
# # get information from the ME psim that was built (can do user provided information instead)
# idx  = psim_me[psim_me['HGNC ID'] == gene1_id].index
# ptms_ = dict(zip(ptm_keys, psim_me.loc[idx, ptm_cols].iloc[0,:].tolist()))
# ptms_ = {k:v for k,v in ptms_.items() if v != 0}
# ptms_['Phosphorylation'] = 3 # example of one that will not be considered
# psim_me.Location = psim_me.Location.replace(float('nan'), '0')
# fl = [i.replace('[', '') for i in psim_me.loc[idx, 'Location'].tolist()]
# fl = [i.replace(']', '') for i in fl]
# ensg_ = psim_me.loc[idx, 'Ensembl gene ID'].tolist()[0]
# pm,m,p = psim_me.loc[idx, 'premrna_seq'].tolist()[0], psim_me.loc[idx, 'mrna_seq'].tolist()[0], psim_me.loc[idx, 'protein_seq'].tolist()[0]


# # In[7]:


# # initialize the gene class
# gene1 = gene_information(human_model, hgnc_id = gene1_id, ptms = ptms_, 
#                          tmd = psim_me.loc[idx,'TMD'].tolist()[0], sp = sp_dict[psim_me.loc[idx,'SP'].tolist()[0]])
# print(gene1.module)
# print(gene1.hgnc_id)
# print(gene1.sp)
# print(gene1.ptms)
# print(gene1.tmd)


# # In[8]:


# # get the gene's final locations, final_locations list does not need to be specified for machinery
# gene1.get_final_locations(metabolic_model = human_model, final_locations=fl)
# print(gene1.final_locations)


# # In[9]:


# # method 1 - user provided
# gene1.get_sequences(premrna_seq=pm, mrna_seq=m, protein_seq=p)
# print(gene1.premrna_seq[1:20])
# # method 2 for gene sequences - REST API
# gene1.get_sequences(ensg_id = ensg_)
# print(gene1.mrna_seq[1:20])


# # In[10]:


# gene1.protein_mass


# # In[11]:


# gene1.check_gene_information()

