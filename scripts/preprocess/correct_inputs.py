#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
import cobra.manipulation.delete as c_del

import pandas as pd
import warnings
import os
import json
import itertools
import more_itertools as mit
import copy

import logging
logging.basicConfig()
logger = logging.getLogger(cobra.__name__)
logger.setLevel(logging.ERROR)

# make sure the creat_environment function from the preprocesss script is run before thise
import sys
sys.path.insert(1, '../../scripts/')
from utils.load_environmental_variables import *
from preprocess import parse_complex 

compartments_me = {'c': 'cytosol',  'l': 'lysosome', 'm': 'mitochondria', 'r': 'endoplasmic reticulum', 
                'e': 'extracellular space', 'x': 'peroxisome/glyoxysome', 'n': 'nucleus', 'g': 'golgi apparatus',
                'i': 'inner mitochondrial compartment', 'pm': 'plasma membrane', 'b': 'boundary'}


# In[2]:


def bool_metabolite(m_id, compartment, m_model):
    try: 
        m_id = '_'.join(m_id.split('_')[:-1])
        m = m_model.metabolites.get_by_id(m_id + '_' + compartment )
        return True, m 
    except:
        return False, None

def correct_model(model = input_data_path + 'recon2_2.xml'):
    '''Makes necessary changes to cobrapy model, largely based on issues encountered with Recon2.2
    
    Parameters
    ----------
    model_file: str or cobra.Model
        if str: 'full/path/to/input_model.xml'
    
    Writes corrected model to outdir/processed_recon2_2.xml
    '''  
    
    required_metabolites = json.load(open(build_files_path + "required_metabolic_model_metabolites.json"))
    rmd = pd.read_csv(build_files_path + 'required_metabolic_model_metabolites.csv', index_col = 0)
    
    if type(model) == str:
        if not os.path.isfile(model):
            raise ValueError('Specified file does not exist')
        elif os.path.splitext(model)[1] != '.xml':
            raise ValueError('Specified file must be an sbml model with extentsion ".xml"')
        m_model = cobra.io.read_sbml_model(model)
    elif not isinstance(model, cobra.Model):
        raise ValueError('Model arg must either by a cobrapy model or specify a path to a sbml file of a cobrapy model')
    
    # check for correct compartments
    different_compartments = list(set(m_model.compartments.keys()).difference(compartments_me.keys()))
    if len(different_compartments)>0 and different_compartments != ['']:
        err = 'The input metabolic model contains compartments not considered by the ME model. '
        err += 'Please remove the following compartments from your model: ' + ', '.join(different_compartments)
        raise ValueError(err)

    # incase GPR has redundant complexes (recon2.2 had atleast one instance of this - id = OIVD1m)
    for r in m_model.reactions:
        if 'and' in r.gene_reaction_rule and 'or' in r.gene_reaction_rule:
            machinery_final = parse_complex.eval_complex(r.gene_reaction_rule)

            idx = list(itertools.combinations(range(len(machinery_final)), 2))

            rm = list()
            for i in idx:
                if machinery_final[i[0]] == machinery_final[i[1]]:
                    rm.append(i[1])

            if len(rm) > 0:
                msg = r.id + ' contains redundant complexes according to GPR, editing GPR'
                warnings.warn(msg)
                machinery_final = [machinery_final[i] for i in range(len(machinery_final)) if i not in rm]
                new_gpr = ''
                for complex_ in machinery_final:
                    if type(complex_) == list:
                        new_gpr += '(' + ' and '.join(complex_) + ')' + ' or '
                    else:
                        new_gpr += complex_ + ' or '
                new_gpr = new_gpr[:-4]

                m_model.reactions.get_by_id(r.id).gene_reaction_rule = new_gpr
    
    # remove psuedogene w/ no sequences
    g = [g for g in m_model.genes if g.id == 'HGNC:4686'] 
    if len(g) > 0:
        c_del.remove_genes(cobra_model = m_model, gene_list = g, remove_reactions = True)
            
    # check for minimum required metabolites
    all_required_metabolites = [item for sublist in list(required_metabolites.values()) for item in sublist]
    all_metabolites = [m.id for m in m_model.metabolites]
    missing_metabolites = sorted(set(all_required_metabolites).difference(all_metabolites))
    
    comp_ = ['c', 'n', 'r', 'g', 'm', 'l', 'x', 'i', 'e', 'b', 'pm']

    for m_id in missing_metabolites:
        nc = m_id.split('_')[-1]

        counter_ = 0
        found = False
        while not found and counter_ < len(comp_):
            cp_ = comp_[counter_]
            found, m_t = bool_metabolite(m_id, cp_, m_model) 
            counter_ += 1

        if found: # metabolite was found in another compartment of the model, add a transport
            new_metab = m_t.copy()
            new_metab.compartment = nc
            new_metab.id = m_id

            transport_rxn = cobra.Reaction('_'.join(m_t.id.split('_')[:-1]) + 't' + nc)
            transport_rxn.add_metabolites({m_t: -1, new_metab: 1})
            transport_rxn.lower_bound = -1000

            m_model.add_reactions([transport_rxn])
        else: 
            core_id = '_'.join(m_id.split('_')[:-1])

            wrn = core_id + ' does not exist in model. Adding to compartment ' + nc + ' via sink '
            wrn += 'This allows ' + core_id + ' to be in the model at no cost.'
            warnings.warn(wrn)

            new_metab = cobra.Metabolite(m_id)
            info = rmd.loc[core_id,]

            new_metab.name = info['name']
            new_metab.compartment = nc
            new_metab.charge = int(info['charge'])
            new_metab.elements = eval(info['elements'])
            new_metab.formula = info['formula']

            if new_metab.compartment != 'e':
                m_model.add_boundary(metabolite = new_metab, type = 'sink')
            else:
                m_model.add_boundary(metabolite = new_metab, type = 'exchange')

    print('Check for the recon2.2 HGNC:HGNC error')
    # correct genes with HGNC:HGNC:, recon2.2 has this
    genes_to_format = [g.id for g in m_model.genes if 'HGNC:HGNC:' in g.id]
    
    if len(genes_to_format) > 0:
        warnings.warn('Your metabolic model contains genes with HGNC:HGNC:####, changing to HGNC:####')
        correct_format = ['HGNC:' + g.split('HGNC:')[-1] for g in genes_to_format]
        formatting_dict = dict(zip(genes_to_format, correct_format))

        for g in genes_to_format:
            reactions = list(m_model.genes.get_by_id(g).reactions)
            for r in reactions:
                r.gene_reaction_rule = r.gene_reaction_rule.replace(g, formatting_dict[g])

    print('Remove genes not participating in reactions')
    genes_to_remove = [g for g in m_model.genes if len(g.reactions)==0]
    if len(genes_to_remove) > 0:
        if sorted(genes_to_format) != sorted([g.id for g in genes_to_remove]):
            warnings.warn('Your metabolic model contains genes not involved in reactions, removing these genes')
        
        c_del.remove_genes(cobra_model = m_model, gene_list = genes_to_remove, remove_reactions = True)
    
    # remove biomass        
    metabolites_1 = [m.id for m in m_model.metabolites]
    try:
        biomass_reaction = m_model.reactions.get_by_id('biomass_reaction')
        try:
            biomass_m = list(biomass_reaction.metabolites.keys())

            m_model.remove_reactions(['biomass_reaction'], remove_orphans = True)

            for bm in biomass_m:
                reactions = list(bm.reactions)
                if len(reactions) != 1:
                    raise ValueError('Unexpected  reactions associated with biomass constituent')
                bmr = reactions[0]
                m_model.remove_reactions([bmr], remove_orphans = True)
        except:
            raise ValueError('Biomass reactions formatted differently than expected, please emulate RECON2.2')
    except:
        wrn_ = 'No biomass reaction identified in input model. Assuming that biomass reactions and metabolites '
        wrn_ = " are not present. If present, please make sure the biomass reaction id is 'biomass_reaction'"
        warnings.warn(wrn_)

    if len([m.id for m in m_model.metabolites if 'biomass' in m.id]) != 0:
        err_ = 'Extraneous biomass metabolites not associated with the biomass reaction are present,'
        err_ += ' please remove from input cobrapy model'
        raise ValueError(err_)

    rm = sorted(set(metabolites_1).difference([m.id for m in m_model.metabolites]))
    if len([i for i in rm if 'biomass' not in i]) != 0:
        warnings.warn('Non biomass metabolites removed as orphan metabolites from biomass reactions')
    cobra.io.write_sbml_model(cobra_model = m_model, filename = processed_data_path + 'corrected_model.xml')
    del m_model
    


# In[6]:


def check_non_machinery(non_machinery = input_data_path + 'non_machinery.txt'):
    '''Runs checks on non-machinery input list. Non-machinery are categorized as any proteins
    to express in the ME-Model that are not listed in the GPR. 
    
    Parameters
    ----------
    non_machinery: str, list, or None
        if list: each entry is a string representing a gene id (HGNC ID format) to express in ME-Model
        if str: "path/to/non_machinery.txt" is a text file containing the same gene ids as described in the list, with separator = '\n'
            
    '''
    
    # load files------------------------------------
    if type(non_machinery) == str:
        if os.path.isfile(non_machinery):
            non_machinery = open(non_machinery).read().splitlines()
        else:
            raise ValueError('Non-machinery file does not exist')
    elif non_machinery is None:
        non_machinery = []
    elif type(non_machinery) != list:
        raise ValueError('The passed non_machinery argument is invalid')
    
    expression_machinery = list(open(build_files_path + 'expression_machinery.txt').read().splitlines())
    if os.path.isfile(processed_data_path + 'corrected_model.xml'):
        m_model = cobra.io.read_sbml_model(processed_data_path + 'corrected_model.xml')
    else:
        raise ValueError('Please run preprocess.correct_inputs.correct_model first')
    metabolic_machinery = [g.id for g in m_model.genes]
    #------------------------------------
    if len([i for i in non_machinery if not i.startswith('HGNC:')]) > 0:
        raise ValueError('All non-machinery must be in HGNC ID format')
    
    
    if len(list(set(non_machinery).intersection(expression_machinery + metabolic_machinery))) > 0:
        msg = 'You have specificied some non-machinery genes which overlap with genes in either the metabolic model '
        msg += 'or the expression module of the ME model. The current format of the model only allows non GPR '
        msg += 'genes to be categorizes as non-machinery. Removing these from the non-machinery list'
        warnings.warn(msg)

        non_machinery = list(set(non_machinery).difference(expression_machinery + metabolic_machinery))
    
    del m_model
    return non_machinery, expression_machinery, metabolic_machinery


# In[20]:


def get_status(psim_me):
    '''Checks sequence columns for validity'''
    psim = psim_me.copy()
    psim['Status'] = 1

    # annotate invalid entries
    for col in ['PROTEIN_SEQ', 'MRNA_SEQ', 'PREMRNA_SEQ', 'HGNC_ID']:
        psim.loc[psim[psim[col].isna()].index,'Status'] = 0
    
    premrna_l = psim.PREMRNA_SEQ.apply(lambda x: len(x) if type(x) == str else float('nan'))
    mrna_l = psim.MRNA_SEQ.apply(lambda x: len(x) if type(x) == str else float('nan'))
    protein_l = psim.PROTEIN_SEQ.apply(lambda x: len(x) if type(x) == str else float('nan'))
    psim.loc[psim[(premrna_l < mrna_l) | (mrna_l < (3*protein_l))].index, 'Status'] = 0

    temp = psim[(premrna_l == mrna_l)]
    psim.loc[temp[temp.MRNA_SEQ != temp.PREMRNA_SEQ].index, 'Status'] = 0
    
    for col in ['MRNA_SEQ', 'PREMRNA_SEQ']:
        idx = psim[psim[col].apply(lambda x: len(set(list(x)).difference(['A', 'U', 'G', 'C', 'N'])) > 0 if type(x) == str else True)].index
        psim.loc[idx, 'Status'] = 0

    amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    idx = psim[psim['PROTEIN_SEQ'].apply(lambda x: len(set(list(x)).difference(amino_acids + ['X', 'U'])) > 0 if type(x) == str else True)].index
    psim.loc[idx, 'Status'] = 0
    
    
    return psim

def correct_psim(psim_df = input_data_path + 'psim_me.h5', fill_na = 'default', 
                non_machinery = None):
    '''Makes sure PSIM has all necessary correct information to build ME Model
    
    *Note, the default psim_file, build/psim_me.h5, is a PSIM generated from MANE/RefSeq Select isoforms.
    We refer to this as the gold standard. 
    
    Parameters
    ----------
    psim_df: pd.DataFrame or str, defaults to 'input_path/psim_me.h5' 
        See PSIM_README for details on format of psim
        str: full/path/to/psim.csv or psim.h5
    fill_na: str
        options ['default', 'select']
        if default: will fill incomplete values with default values (see PSIM_README for details)
        if select: will fill incomplete values with the gold standard PSIM when available, otherwise with default
        
        Note: this will not deal with user-provided incorrected values for optional columns, those will revert to default in the gene_information class
        
        Exceptions:
            for required columns, if incorrect in input psim, will fill with  the gold standard PSIM. required columns include *_SEQ and LOCATION for non-machinery
            if PTR is present in input PSIM and a tissue is specified most often in that column (ignoring NaN), all nan values in that column will default to that tissue
            if non-machinery do not specify an appropriate LOCATION, will fill with  the gold standard PSIM
    non_machinery: str, list, or None
        if list: each entry is a string representing a gene id (HGNC ID format) to express in ME-Model
        if str: "path/to/non_machinery.txt" is a text file containing the same gene ids as described in the list, with separator = '\n'
    
    Returns
    ----------
    revised_genes: dict
        'added' is a list of genes missing (relative to m_model and non-machinery genes) in PSIM that were added 
        'sequences' key is a superset of added, includes all genes with adjusted sequences
        'non-machinery locations' is for genes in non-machinery that did not have an appropriately specified location
        
    Also writes corrected PSIM to outdir/corrected_psim_me.h5 (specified in preprocess.create_environment)
    '''
    # run basic non-machinery check
    non_machinery, expression_machinery, metabolic_machinery = check_non_machinery(non_machinery = non_machinery)

    # define the required/optional columns-----------------------------------------------------------------
    user_provided = ['HGNC_ID'] # must be fully provided by user
    essential_cols = ['PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ'] # required but can be filled by build/psim_me
    optional_cols = ['POLYA_LENGTH', 'N_EXONS', 'TMD', 'SP', 'DSB', 'GPI', 'OG', 'NG', 
                    'ALPHA_M', 'ALPHA_P', 'PTR']
    nm_cols = ['LOCATION']
    all_columns = user_provided + essential_cols + optional_cols + nm_cols
    
    # load the MANE/RefSEQ Select PSIM---------------------------------------------------------------------
    psim_gold = pd.read_hdf(build_files_path + 'psim_me.h5')
    psim_gold = psim_gold[psim_gold.Status != 0] # drop genes that won't work with model
    psim_gold = psim_gold[all_columns]
    
    # load input PSIM---------------------------------------------------------------------------------------
    if type(psim_df) == str:
        if not os.path.isfile(psim_df):
            raise ValueError('The specified PSIM file does not exist')
        filename, file_extension = os.path.splitext(psim_df)
        if file_extension == '.h5':
            psim_me = pd.read_hdf(psim_df)
        elif file_extension == '.csv':
            psim_me = pd.read_csv(psim_df, index_col = 0)
        else:
            raise ValueError('PSIM must be in .csv or .h5 format')
    elif not isinstance(psim_df, pd.DataFrame):
        raise ValueError('The specified psim_df arg is invalid')
    
    psim_me.reset_index(inplace = True, drop = True) 
    # check that required cols are all present and appropriately formatted------------------------------------
    err = False
    for col in user_provided:
        if col not in psim_me.columns.tolist():
            err = True
            break
        if psim_me[col].isna().any():
            err = True
            break
        if psim_me[col].isna().any():
            err = True
            break
        if psim_me[col].unique().shape[0] != psim_me.shape[0]:
            err = True
            break
    if err:
        raise ValueError('The following columns must be present, unique, and completely filled in: ' + ', '.join(user_provided))

    for col in essential_cols:
        if col not in psim_me.columns:
            psim_me[col] = float('nan')
    missing_cols = sorted(set(all_columns).difference(psim_me.columns))
   
    # check genes are present--------------------------------------------------------------------------
    proteins = metabolic_machinery + non_machinery
    psim_me_genes = psim_me.HGNC_ID.tolist()
    psim_gold_genes = psim_gold.HGNC_ID.tolist()

    missing_genes = sorted(set(proteins).difference(psim_me_genes + psim_gold_genes))
    if len(missing_genes)>0: # check all but expression machinery
        raise ValueError('The following specified genes are not present in the provided psim or build/psim: ' + ', '.join(missing_genes))
    else: # add expression machinery to lists
        proteins = sorted(set(proteins + expression_machinery)) 
        revised_genes = {'missing': missing_genes}
        missing_genes = sorted(set(proteins).difference(psim_me_genes))
        
    
    # initialize missing columns and genes----------------------------------------------------------------
    for col in missing_cols:
        psim_me[col] = float('nan')
    psim_me = psim_me[all_columns] # filter out excess columns and order

    for gene in missing_genes:
        psim_me.loc[psim_me.shape[0], :] = [gene] + ([float('nan')]*(psim_me.shape[1]-1))
        
    #fill in required columns independent of fill_na arg--------------------------------------------------
    psim_me = psim_me[psim_me.HGNC_ID.isin(proteins)]
    psim_gold = psim_gold[psim_gold.HGNC_ID.isin(proteins)]
    
    if psim_me.HGNC_ID.unique().shape[0] != psim_me.shape[0]:
        raise ValueError('Duplicate HGNC IDs')
    psim = get_status(psim_me)
    fix = psim[psim.Status == 0].HGNC_ID.tolist()
    
    
    psim_gold.index = psim_gold.HGNC_ID
    if len(set(fix).difference(psim_gold_genes)) == 0:
        psim_me.index = psim_me.HGNC_ID
        psim_me.loc[fix, essential_cols] = psim_gold.loc[fix, essential_cols] 
        psim_me.reset_index(drop = True, inplace = True)
        revised_genes['sequences'] = fix   
    else:
        mssg = ' The following required genes have incorrect sequences in the PSIM: ' + ', '.join(list(set(fix).difference(psim_gold_genes)))
        mssg += ' mRNA sequence lengths must be <= premrna sequence lengths and >= 3*protein sequence length'
        mssg += ' Furthermore, they must have only the allowed letters'
        raise ValueError(mssg)
    
    # ptr specific changes-----------------------------------------------------------------------------
    fill_ptr = True
    if 'PTR' not in missing_cols:
        max_val = psim_me['PTR'].dropna().value_counts().index.tolist()
        if len(max_val) > 0 and type(max_val[0]) == str:
            ptr = pd.read_csv(build_files_path + 'PTR_Gagneur_processed.tsv', sep = '\t', index_col = 0)
            ptr.drop(columns = ['ENSG_ID'], inplace = True)
            ptr.columns = pd.Series(ptr.columns).apply(lambda x: x.split('_')[0] if '_PTR' in x else x).tolist()
            if max_val[0] in ptr.columns.tolist():
                psim_me['PTR'] = max_val[0]
                fill_ptr = False
    
    # fill_na values-----------------------------------------------------------------------------------
    psim_me.index = psim_me.HGNC_ID
    psim_gold.index = psim_gold.HGNC_ID
    
    if fill_na == 'select':
        for col in optional_cols:
            if col != 'PTR' or fill_ptr:
                temp = psim_me[psim_me[col].isna()]
                temp = psim_gold[psim_gold.HGNC_ID.isin(temp.HGNC_ID)]
                temp = temp[temp[col].notna()]
                psim_me.loc[temp.index, col] = temp.loc[:, col]
    elif fill_na != 'default':
        raise ValueError('Invalid value passed for fill_na')
    
    # deal with non_machinery location---------------------------------------------------------------
    if len(non_machinery) > 0 and 'LOCATION' not in psim_me.columns.tolist():
        psim_me['LOCATION'] = float('nan')
    compartments = ['[c]', '[e]', '[l]', '[m]', '[r]', '[n]', '[g]', '[x]', '[b]', '[i]', '[pm]']
    temp = psim_me[psim_me.HGNC_ID.isin(non_machinery) & ((psim_me.LOCATION.isna()) | psim_me.LOCATION.apply(lambda x: x not in compartments))]
    err = False
    if len(set(temp.HGNC_ID).difference(psim_gold.HGNC_ID)) > 0:
        raise ValueError('The final location for the following non-machinery must be specified: ' + ', '.join(list(set(temp.HGNC_ID).difference(psim_gold.HGNC_ID))))
    temp = psim_gold.loc[temp.index, :]
    if temp.LOCATION.dropna().shape[0] < temp.shape[0]:
        raise ValueError('The final location for the following non-machinery must be specified: ' + ', '.join(temp[temp.LOCATION.isna()].HGNC_ID.tolist()))
    
    if temp.shape[0] > 0:
        psim_me.loc[temp.index, 'LOCATION'] = temp.LOCATION.tolist()
        revised_genes['non-machinery locations'] = temp.HGNC_ID.tolist()
    psim_me.reset_index(inplace = True, drop = True)
    
    del psim_gold
    
    psim_me.to_hdf(processed_data_path + 'corrected_psim.h5', key = 'corrected')
    with open(processed_data_path + 'corrected_non_machinery.txt', 'w'):
        for i in non_machinery:
            f.write(i + '\n')
    return revised_genes

