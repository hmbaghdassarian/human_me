#!/usr/bin/env python
# coding: utf-8

# In[68]:


from Bio.SeqUtils import molecular_weight as calculate_molecular_weight
import random
import itertools
import numpy as np
import pandas as pd
import gc
import multiprocessing
import os
import warnings
import time
n_cores = multiprocessing.cpu_count() # number of cores to use in parallelization

import sys
sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *
from utils import *
from utils_2 import *
import build_mrna_expression_reactions as build_mrna
import build_protein_expression_reactions as build_protein
from polyA_statistics import min_polyA_mean


# Trying to get all combinations of universal variables/edge cases to test expression module. Two such are the nuclear_diffusion_limit and the ptt_length variables. As we can see, there will be no case such that a protein has a length to undergo post-transcriptional translation but not undergo reversible diffusion as the means of nuclear transport. This is because even the heaviest by molecular weight amino acid, at the length for post-transcriptional translation, is not above the nuclear diffusion limit.  

# In[2]:


# order the amino acids by molecular weight
mw_map = dict()
for k in seq_amino_acid_map_c.keys():
    k_ = calculate_molecular_weight(k, seq_type = 'protein')
    if k_ not in mw_map.keys():
        mw_map[k_] = [k]
    else:
        mw_map[k_] += [k]
mw_map_2 = dict()    
for m in sorted(set(mw_map.keys())):
    aa_code = mw_map[m]
    if len(mw_map[m]) == 1:
        mw_map_2[aa_code[0]] = m
    else:
        for aa_code_ in aa_code:
            mw_map_2[aa_code_] = m
mw_map = mw_map_2
mw_map


# In[3]:


protein_lengths = [80, ptt_length, ptt_length + 1, 300,900] 
# first three list element lengths will always be under nuclear diffusion limit
# fourth element can be either over or under, depending on amino acids included
# fifth will always be over nuclear diffusion limit
# ptt_length is those that undergo post or co-translational translocation
# we want to create two random protein sequences per situation

max_aa_map = dict()
for l in protein_lengths:
    max_amino_acid = None
    for k,v in mw_map.items():
        if v*l < nuclear_diffusion_limit:
            max_amino_acid = k
        else:
            break
    max_aa_map[l] = max_amino_acid
max_aa_map


# In[4]:


protein_sequences = list()
n_sequences = 1#2
sorted_amino_acids = list(mw_map.keys())
for l in protein_lengths:
    # this provides all the possible amino acids to choose to keep under or overthe nuclear diffusion limit
    # at a given length
    if max_aa_map[l] != None:
        aa_under = sorted_amino_acids[:sorted_amino_acids.index(max_aa_map[l])+1]
    else:
        aa_under = []
    aa_over = sorted(set(sorted_amino_acids).difference(aa_under))

    if len(aa_over) + len(aa_under) != len(amino_acids):
        print(len(aa_over))
        print(len(aa_under))
        raise ValueError('Not all amino acids are considered')
    
    # create two protein sequences of each protein length for both those under and over nuclear diffusion limit
    if len(aa_under) > 0:
        for i in range(n_sequences):
            protein_sequences.append(''.join(random.choices(aa_under, k = l)))
    if len(aa_over) > 0:
        for i in range(n_sequences):
            protein_sequences.append(''.join(random.choices(aa_over, k = l)))
            
        
if len(protein_sequences) != n_sequences*(len(protein_lengths)+1):
    raise ValueError('Did not get all sequences')


# Since there is no check that the mrna sequence corresponds to the protein sequence EXCEPT that it is atleast 3x the length, we will generate random nucleotide sequences for each protein sequence of that minimum length, as well as an mrna sequence that is 200 nucleutides longer and 500 nts longer than the minimum length. 
# 
# Additionally, we will create 3 premrna sequences per mrna sequence, one with 0 introns, with 1 intron (under the rate), 1 intron (at the rate), 2 introns (just over the rate), and with 6 introns. 
# 
# Thus we should have many possible combinations of: post vs co transcriptional translation, nuclear diffusion vs active nuclear transport, intron and no intron. 
# 
# *Note had to exclude some of these to limit number of total combinations

# In[5]:


seq_df = pd.DataFrame(columns = ['PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ'])
counter = 0
for i in range(len(protein_sequences)):
    ps = protein_sequences[i]
    min_mrna_l = len(ps)*3

    mrna_seq_min = ''.join(random.choices(['A', 'U', 'C', 'G'], k = min_mrna_l)) 
#     mrna_seq_larger = mrna_seq_min + ''.join(random.choices(['A', 'U', 'C', 'G'], k = 200))
#     mrna_seq_largest = mrna_seq_larger + ''.join(random.choices(['A', 'U', 'C', 'G'], k = 300))
    mrna_seq_largest = mrna_seq_min + ''.join(random.choices(['A', 'U', 'C', 'G'], k = 500))

    for mrna_seq in [mrna_seq_min, mrna_seq_largest]:#,mrna_seq_larger] :
        L_mrna = len(mrna_seq)
        pre_mrna_sequences = list()
        premrna_seq_no_intron = mrna_seq
        if (L_mrna + 1) * rate_intron < 1: # those with introns less than minimal length
            premrna_seq_one_intron_A = mrna_seq + ''.join(random.choices(['A', 'U', 'C', 'G'], k = 1))
            pre_mrna_sequences.append(premrna_seq_one_intron_A)
        else: 
            premrna_seq_one_intron_A = None
        premrna_seq_one_intron_B = mrna_seq + ''.join(random.choices(['A', 'U', 'C', 'G'], k = round(1/rate_intron) - L_mrna))
#         premrna_seq_two_introns = mrna_seq + ''.join(random.choices(['A', 'U', 'C', 'G'], k = round(2.001/rate_intron) - L_mrna))
        premrna_seq_six_introns = mrna_seq + ''.join(random.choices(['A', 'U', 'C', 'G'], k = round(6.001/rate_intron) - L_mrna))

        pre_mrna_sequences += [premrna_seq_one_intron_B, premrna_seq_six_introns]#, premrna_seq_two_introns]
        for pms in pre_mrna_sequences:
            seq_df.loc[counter, :] = [pms, mrna_seq, ps]
            counter += 1

if seq_df.shape[0] <  len(protein_sequences)*3*2:
    raise ValueError('Not all combinations of sequences accounted for')
    
locations = list(compartments.keys())
location_combos = list()
for l in range(1,len(locations)):
    location_combos += [list(item) for item in list(itertools.combinations(locations, l))]
location_combos.append(locations)



# In[6]:


sp = [None, True, False]
polyA_length = [None, 0, min_polyA_mean - 10, min_polyA_mean + 100]
gpi = [None, 1]
dsb = [None, 3]
og = dsb.copy()
n_introns = [None, 0, 3]
tmd = dsb.copy()
seq_combos = [seq_df.loc[i,:].tolist() for i in seq_df.index]

all_combos = [seq_combos, sp, polyA_length, dsb, og, tmd, location_combos, gpi, n_introns]
all_combos = list(itertools.product(*all_combos))

psim_test = pd.DataFrame(all_combos)
del all_combos
psim_test.columns = ['Seq', 'SP', 'POLYA_LENGTH', 'DSB', 'OG', 'TMD', 'LOCATION', 'GPI', 'N_INTRONS']

df2 = pd.DataFrame(psim_test['Seq'].to_list(), columns=seq_df.columns.tolist())
psim_test = pd.concat([df2, psim_test[psim_test.columns.tolist()[1:]]], axis = 1)

del df2
del seq_df
del seq_combos
del location_combos

# psim_test['HGNC_ID'] = ['HGNC_' + str(i) for i in psim_test.index] # underscore instead of colon to not be in machinery


# In[7]:


if not os.path.isdir(local_data_path + 'processed/test_modules'):
    os.mkdir(local_data_path + 'processed/test_modules/')

n_iter = psim_test.shape[0]
psim_test_cols = psim_test.columns.tolist()
psim_test = np.array(psim_test)
n_splits = 10000 # number of files to create, limits ram on any given run
split_length = round(n_iter/n_splits)
idxes = sorted(set(list(range(0,n_iter,split_length)) + [n_iter]))
gc.collect()



def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__
warnings.filterwarnings("ignore")
def test_expression_module(i):
    if (i+1)%1000 == 0:
        print('{} of {}. {}% complete'.format(i+1, split_length, (i+1)/split_length * 100))
    error = list()
    entries = list(psim_test[i])#psim_test.loc[i,]
    
    
    try:
        blockPrint()
        gene_info = gene_information(metabolic_model=human_model, hgnc_id = 'HGNC_' + str(i), 
                    premrna_seq = entries[0], mrna_seq = entries[1], 
                    protein_seq = entries[2], 
                    ptms = dict(zip(['dsb', 'og', 'gpi'],[entries[5], entries[6], entries[9]])),
                    tmd = entries[7], sp = entries[3], polyA_length = entries[4], 
                    n_introns = entries[10])
        gene_info.get_final_locations(metabolic_model = human_model, final_locations = entries[8])
        gene_info.check_gene_information()
        enablePrint()
        try:
            mrna_reactions = build_mrna.mrna_expression(gene_info)
            if len([r.id for r in mrna_reactions if len(r.check_mass_balance()) > 0]) != 0:
                error += ['mrna reaction mass balance']
            del mrna_reactions
            try:
                protein_expression_reactions, protein_metabolites = build_protein.get_protein_expression_reactions(gene_info)
                if len([r.id for r in protein_expression_reactions if len(r.check_mass_balance()) > 0]) != 0:
                    error += ['protein reaction mass balance']

                if sorted([p.compartment for p in protein_metabolites]) != sorted(entries[8]):
                    error += ['protein metabolite incorrect compartment']
                del protein_expression_reactions
                del protein_metabolites
                del gene_info
            except:
                error += ['could not create protein reactions']


        except:
            error += ['could not create mrna reactions']

    except: 
        error += ['could not create gene information object']
    gc.collect()
    if len(error)>0:
        return (i, ';'.join(error))
    else:
        del error
        return (float('nan'),float('nan'))


# In[14]:


# for i in list(range(n_iter)):#range(10):
#     res.append(test_expression_module(i))
start_time = time.time()
if __name__ == '__main__':
    for i in range(len(idxes)-1):
        print('Starting data frame number: {} of {}'.format(i+1, n_splits))
        print('------------------------------------------------------------')
        idx = list(range(idxes[i], idxes[i+1]+1))
        
        pool = multiprocessing.Pool(processes=round(n_cores/2))
        res = pool.map(test_expression_module, idx)#range(20))
        pool.close()

        res = [i for i in res if not pd.isna(i[0])]
        if len(res) > 0:
            fail_index, error_message = list(zip(*res))
            psim_fail = psim_test[fail_index].copy()
            psim_fail = pd.DataFrame(psim_fail, columns = psim_test_cols)
            psim_fail['ERROR'] = error_message
        else:
            psim_fail = pd.DataFrame(columns = psim_test_cols + ['ERROR'])
        
        gc.collect()
        file_name = local_data_path + 'processed/test_modules/test_expression_module' + str(i) + '.csv'
        psim_fail.to_csv(file_name, compression = 'gzip')
        del psim_fail
        gc.collect()
        
        print("--- Time passed so far: {:.3f} hours ---".format((time.time() - start_time)/3600))
        
    del psim_test
    print('COMPLETE')
    gc.collect()

