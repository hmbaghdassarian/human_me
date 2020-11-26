#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra

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
from utils.load_environmental_variables import root_path, build_files_path, processed_data_path
from preprocess import parse_complex 

compartments_me = {'c': 'cytosol',  'l': 'lysosome', 'm': 'mitochondria', 'r': 'endoplasmic reticulum', 
                'e': 'extracellular space', 'x': 'peroxisome/glyoxysome', 'n': 'nucleus', 'g': 'golgi apparatus',
                'i': 'inner mitochondrial compartment', 'pm': 'plasma membrane', 'b': 'boundary'}


# In[2]:


required_metabolites = json.load(open(build_files_path + "required_metabolic_model_metabolites.json"))
rmd = pd.read_csv(build_files_path + 'required_metabolic_model_metabolites.csv', index_col = 0)


# In[3]:


def bool_metabolite(m_id, compartment, m_model):
    try: 
        m_id = '_'.join(m_id.split('_')[:-1])
        m = m_model.metabolites.get_by_id(m_id + '_' + compartment )
        return True, m 
    except:
        return False, None

def correct_model(model_file = root_path + 'recon2_2.xml', 
                 psim_file = root_path + 'psim_recon2_2.csv'):
    '''
    Makes some necessary changes to cobrapy model, largely based on issues encountered with Recon2.2.
    model_file is path/to/cobra_smbl_model
    psim_file is path/to/psim_csv
    
    '''
    
    human_model = cobra.io.read_sbml_model(model_file)
    
    psim_me = pd.read_csv(psim_file)
    psim_me.reset_index(inplace = True, drop = True)
    if not ('HGNC_ID' in psim_me.columns):
        raise ValueError('PSIM_ME must have "HGNC_ID" as a column')
    psim_me_genes = psim_me.HGNC_ID.tolist()

    
    # check for correct compartments
    different_compartments = list(set(human_model.compartments.keys()).difference(compartments_me.keys()))
    if len(different_compartments)>0 and different_compartments != ['']:
        err = 'The input metabolic model contains compartments not considered by the ME model. '
        err += 'Please remove the following compartments from your model: ' + ', '.join(different_compartments)
        raise ValueError(err)

    # incase GPR has redundant complexes (recon2.2 had atleast one instance of this - id = OIVD1m)
    for r in human_model.reactions:
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

                human_model.reactions.get_by_id(r.id).gene_reaction_rule = new_gpr
            
    # check for minimum required metabolites
    all_required_metabolites = [item for sublist in list(required_metabolites.values()) for item in sublist]
    all_metabolites = [m.id for m in human_model.metabolites]
    missing_metabolites = sorted(set(all_required_metabolites).difference(all_metabolites))
    
    comp_ = ['c', 'n', 'r', 'g', 'm', 'l', 'x', 'i', 'pm']
#     reactions_to_add = list()

    for m_id in missing_metabolites:
        nc = m_id.split('_')[-1]

        counter_ = 0
        found = False
        while not found and counter_ < len(comp_):
            cp_ = comp_[counter_]
            found, m_t = bool_metabolite(m_id, cp_, human_model) 
            counter_ += 1

        if m_t is not None: # metabolite was found in another compartment of the model, add a transport
            new_metab = m_t.copy()
            new_metab.compartment = nc
            new_metab.id = m_id

            transport_rxn = cobra.Reaction('_'.join(m_t.id.split('_')[:-1]) + 't' + nc)
            transport_rxn.add_metabolites({m_t: -1, new_metab: 1})
            transport_rxn.lower_bound = -1000
            
            human_model.add_reactions([transport_rxn])#reactions_to_add.append(transport_rxn)
        else: # metabolite does not exist anywhere in model, add an exchange and transport from [e] to compartment
            core_id = '_'.join(m_id.split('_')[:-1])

            wrn = core_id + ' does not exist in model. Adding to compartment ' + nc + ' via exchange and '
            wrn += ' transport reactions. This allows ' + core_id + ' to be in the model at no cost.'
            warnings.warn(wrn)

            new_metab = cobra.Metabolite(m_id)
            info = rmd.loc[core_id,]

            new_metab.name = info['name']
            new_metab.compartment = nc
            new_metab.charge = int(info['charge'])
            new_metab.elements = eval(info['elements'])
            new_metab.formula = info['formula']

            exchange_rxn_1 = cobra.Reaction('EX_' + core_id + '_b')
            exchange_rxn_1.name = exchange_rxn_1.id
            exchange_rxn_1.lower_bound = -1000

            if new_metab.compartment != 'b':
                m_t_1, m_t_2 = new_metab.copy(), new_metab.copy()
                m_t_1.compartment, m_t_2.compartment = 'b', 'e'
                m_t_1.id, m_t_2.id = core_id + '_b', core_id + '_e'

                # boundary exchange
                exchange_rxn_1.add_metabolites({m_t_1: -1})

                # extracellular exchange
                exchange_rxn_2 = cobra.Reaction('EX_' + core_id + '_LPAREN_e_RPAREN_')
                exchange_rxn_2.lower_bound,exchange_rxn_2.upper_bound  = -float('inf'), float('inf')
                exchange_rxn_2.name = 'exchange reaction for ' + core_id
                exchange_rxn_2.add_metabolites({m_t_2: -1, m_t_1: 1})

                # exchange to compartment
                transport_rxn = cobra.Reaction(core_id + 't' + nc)
                transport_rxn.add_metabolites({m_t_2: -1, new_metab: 1})
                transport_rxn.lower_bound = -1000
                
                human_model.add_reactions([exchange_rxn_1, exchange_rxn_2, transport_rxn])#reactions_to_add += [exchange_rxn_1, exchange_rxn_2, transport_rxn]

            else:
                
                exchange_rxn_1.add_metabolites({new_metab: -1})
                human_model.add_reactions([exchange_rxn_1])#reactions_to_add.append(exchange_rxn_1)

#     human_model.add_reactions(reactions_to_add)

    print('Check for the recon2.2 HGNC:HGNC error')
    # correct genes with HGNC:HGNC:, recon2.2 has this
    genes_to_format = [g.id for g in human_model.genes if 'HGNC:HGNC:' in g.id]
    # make sure this is not a user formatting thing
    # assuming corrected of just one HGNC version is in PSIM, will check later (correct_format[i] in psim.HGNC_ID)
    genes_to_format = [g for g in genes_to_format if g not in psim_me_genes] 

    if len(genes_to_format) > 0:
        warnings.warn('Your metabolic model contains genes with HGNC:HGNC:####, changing to HGNC:####')
        correct_format = ['HGNC:' + g.split('HGNC:')[-1] for g in genes_to_format]
        formatting_dict = dict(zip(genes_to_format, correct_format))

        for g in genes_to_format:
            reactions = list(human_model.genes.get_by_id(g).reactions)
            for r in reactions:
                r.gene_reaction_rule = r.gene_reaction_rule.replace(g, formatting_dict[g])

    print('Remove genes not participating in reactions')
    genes_to_remove = [g for g in human_model.genes if len(g.reactions)==0]
    if len(genes_to_remove) > 0:
        if sorted(genes_to_format) != sorted([g.id for g in genes_to_remove]):
            warnings.warn('Your metabolic model contains genes not involved in reactions, removing these genes')

        for g in genes_to_remove:
            human_model.genes.remove(g)
    
    # remove biomass        
    metabolites_1 = [m.id for m in human_model.metabolites]
    try:
        biomass_reaction = human_model.reactions.get_by_id('biomass_reaction')
        try:
            biomass_m = list(biomass_reaction.metabolites.keys())

            human_model.remove_reactions(['biomass_reaction'], remove_orphans = True)

            for bm in biomass_m:
                reactions = list(bm.reactions)
                if len(reactions) != 1:
                    raise ValueError('Unexpected  reactions associated with biomass constituent')
                bmr = reactions[0]
                human_model.remove_reactions([bmr], remove_orphans = True)
        except:
            raise ValueError('Biomass reactions formatted differently than expected, please emulate RECON2.2')
    except:
        wrn_ = 'No biomass reaction identified in input model. Assuming that biomass reactions and metabolites '
        wrn_ = " are not present. If present, please make sure the biomass reaction id is 'biomass_reaction'"
        warnings.warn(wrn_)

    if len([m.id for m in human_model.metabolites if 'biomass' in m.id]) != 0:
        err_ = 'Extraneous biomass metabolites not associated with the biomass reaction are present,'
        err_ += ' please remove from input cobrapy model'
        raise ValueError(err_)

    rm = sorted(set(metabolites_1).difference([m.id for m in human_model.metabolites]))
    if len([i for i in rm if 'biomass' not in i]) != 0:
        warnings.warn('Non biomass metabolites removed as orphan metabolites from biomass reactions')
    cobra.io.write_sbml_model(cobra_model = human_model, filename = processed_data_path + 'corrected_model.xml')
    del human_model
    


# In[4]:


def correct_psim(psim_file = root_path + 'psim_recon2_2.csv'):
    '''
    Runs checks on psim csv.
    psim_file is path/to/psim_csv
    
    '''
    # load files------------------------------------
    psim_me = pd.read_csv(psim_file, index_col = 0)
    psim_me.reset_index(inplace = True, drop = True)
    psim_me_genes = psim_me.HGNC_ID.tolist()
    
    essential_cols = ['HGNC_ID', 'PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ']
    optional_cols = ['POLYA_LENGTH', 'TMD', 'SP', 'N_INTRONS', 'DSB', 'GPI', 'OG', 
                    'MRNA_HALF_LIFE', 'ALPHA_P', 'PTR', 'PTR_TISSUE', 'CONSTANT_PTR'] # NG
    other_cols = ['LOCATION']
    all_columns = essential_cols + optional_cols + other_cols

    expression_psim = pd.read_csv(build_files_path + 'expression_module_psim.csv', index_col = 0)
    expression_machinery = sorted(open(build_files_path + 'expression_machinery.txt').read().splitlines())
    expression_psim = expression_psim[all_columns] 
    
    if not os.path.isfile(processed_data_path + 'corrected_model.xml'):
        raise ValueError('Run correct_inputs.correct_model befor running this function')
    human_model = cobra.io.read_sbml_model(processed_data_path + 'corrected_model.xml')
    metabolic_machinery = sorted([g.id for g in human_model.genes])
    #------------------------------------
    
    # check that essential cols are all present and appropriately formatted
    if len(set(essential_cols).difference(psim_me.columns.tolist()))>0:
        raise ValueError('Must specify the following columns in PSIM: ' + ', '.join(essential_cols))

    # add optional cols that aren't present - FIX THIS
    em_0_idx = psim_me[psim_me.HGNC_ID.isin(expression_machinery)].index
    em_temp = expression_psim.copy()
    em_temp.index = em_temp.HGNC_ID

    for col in optional_cols + other_cols:
        if col not in psim_me.columns.tolist():
            warnings.warn(col + ' not in provided PSIM, adding as default values to gene_information')
            psim_me[col] = float('nan')
            # make sure default values from expression matrix are includes
            psim_me.loc[em_0_idx, col] = em_temp.loc[psim_me.loc[em_0_idx, 'HGNC_ID'].tolist(), col]

    # get rid of excess psim columns and put columns in order
    if len(set(psim_me).difference(all_columns)) > 1:
        msg = 'Some of the columns in the provided PSIM are not used in ME model building and are being removed'
        warnings.warn(msg)
    # do regardless of above if statement to put columns in order 
    psim_me = psim_me[all_columns]
    # similarly put expression psim columns in order; this shouldn't be necessary but just in case
    
    
    # add missing expression module details to psim if it is not already there - CHECK THIS
    em_overlap = set(expression_machinery).intersection(psim_me.HGNC_ID.tolist())
    if len(em_overlap) > 0:
        for e in em_overlap:
            i_p = psim_me[psim_me.HGNC_ID == e].index
            i_e = expression_psim[expression_psim.HGNC_ID == e].index

            vals = psim_me.loc[i_p,:].T.iloc[:,0]

            for col_ in vals[vals.isna()].index:
                psim_me.loc[i_p, col_] = expression_psim.loc[i_e, col_]
    
    # add expression module machinery to psim if it is not already there
    if len(set(expression_machinery).difference(psim_me.HGNC_ID.tolist())) > 0:
        msg = 'Some machinery from expression module of ME model not in PSIM, adding using expression_psim.csv'
        print(msg)

        for i in expression_psim.index:
            if expression_psim.loc[i, 'HGNC_ID'] not in psim_me_genes:
                psim_me.loc[psim_me.shape[0], :] = expression_psim.loc[i,:]

    # check all genes in metabolic model are in the psim
    if len(set(metabolic_machinery).difference(psim_me.HGNC_ID.tolist())) > 0:
        msg = 'Not all genes in provided metabolic model are in the the PSIM '
        msg += 'Please add their information to the PSIM'
        raise ValueError(msg)
    
    psim_me.to_csv(processed_data_path + 'corrected_psim_me.csv')
    del psim_me
    del human_model


# In[5]:


def check_non_machinery(nonmachinery_file = root_path + 'non_machinery.txt'):
    '''
    Runs checks on non-machinery input list
    nonmachinery_file is path/to/nonmachinery_txt_file
    
    '''
    
    # load files------------------------------------
    non_machinery = open(nonmachinery_file).read().splitlines()
    
    expression_psim = pd.read_csv(build_files_path + 'expression_module_psim.csv', index_col = 0)
    expression_machinery = sorted(open(build_files_path + 'expression_machinery.txt').read().splitlines())
    
    if not os.path.isfile(processed_data_path + 'corrected_model.json'):
        raise ValueError('Run correc_inputs.correct_model befor running this function')
    human_model = cobra.io.read_sbml_model(processed_data_path + 'corrected_model.json')
    metabolic_machinery = sorted([g.id for g in human_model.genes])
    
    if not os.path.isfile(processed_data_path + 'corrected_psim_me.csv'):
        raise ValueError('Run correc_inputs.correct_model befor running this function')
    psim_me = pd.read_csv(processed_data_path + 'corrected_psim_me.csv', index_col = 0)
    #------------------------------------
    
    
    # location must be specified for non-machinery
    if len(list(set(non_machinery).intersection(expression_machinery + metabolic_machinery))) > 0:
        msg = 'You have specificied some non-machinery genes which overlap with genes in either the metabolic model '
        msg += 'or the expression module of the ME model. The current format of the model only allows non GPR '
        msg += 'genes to be categorizes as non-machinery. Removing these from the non-machinery list'
        warnings.warn(msg)

        non_machinery = list(set(non_machinery).difference(expression_machinery + metabolic_machinery))

    if len(set(non_machinery).difference(psim_me.HGNC_ID.tolist())) > 0:
        msg = 'Not all genes in provided non-machinery list are in the the PSIM '
        msg += 'Please either remove these from the non-machinery list or add their information to the PSIM'
        raise ValueError(msg)

    loc_vals = psim_me[psim_me.HGNC_ID.isin(non_machinery)].LOCATION.unique().tolist()
    for loc in loc_vals:
        if pd.isna(loc):
            msg = 'All locations for non-machinery must be specificied in the PSIM LOCATION column. '
            msg += 'Each location entry must be formatted as a list of compartments corresponding to the metabolic model compartments'
            raise ValueError(msg)

        msg = 'All locations for non-machinery must be specificied in the PSIM LOCATION column '
        msg += 'Each location entry must be formatted as a list of compartments corresponding to the metabolic model compartments'
        if type(loc) == str:
            if loc[0] != '[' or loc[-1] != ']':
                raise ValueError(msg)
            else:
                loc = loc[1:-1].split(',')
        if type(loc) != list:
            raise ValueError(msg)


        for c in loc:
            if c not in compartments_me.keys():
                raise ValueError('The non-machinery final locations specified are not a compartment considered in the ME model')
    
    del psim_me
    del human_model

