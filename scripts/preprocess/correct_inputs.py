#!/usr/bin/env python
# coding: utf-8

# In[61]:


import cobra
import pandas as pd
import warnings
import os

# make sure the creat_environment function from the preprocesss script is run before thise
import sys
sys.path.insert(1, '../../scripts/')
from utils.load_environmental_variables import root_path, build_files_path, processed_data_path

compartments_ = {'c': 'cytosolic',  'l': 'lysosomal', 'm': 'mitochondrial', 'r': 'endoplasmic reticulum', 
                'e': 'extracellular space', 'x': 'peroxisomal', 'n': 'nuclear', 'g': 'golgi apparatus',
                'i': 'inner mitochondrial compartment', 'pm': 'plasma membrane'}

compartments_me = {'c': 'cytosol',  'l': 'lysosome', 'm': 'mitochondria', 'r': 'endoplasmic reticulum', 
                'e': 'extracellular space', 'x': 'peroxisome/glyoxysome', 'n': 'nucleus', 'g': 'golgi apparatus',
                'i': 'inner mitochondrial compartment', 'pm': 'plasma membrane'}


# In[82]:


def add_required_metabolites(metabolite_id, compartment, human_model, all_metabolites):
    '''Assures metabolite with id metabolite_id exists, and that an exchange reaction with cytoplasm exists.'''
    
    # metabolite exists
    if metabolite_id not in all_metabolites: # assume in 'C'
        met_ = human_model.metabolites.get_by_id(metabolite_id.replace('[' + compartment + ']', '[c]')).copy()
        met_.compartment = compartment
        met_.id = met_.id.replace('[c]', '[' + compartment + ']')
        human_model.add_metabolites([met_])

    # transport from cytosol exists
    met_1 = human_model.metabolites.get_by_id(metabolite_id)
    met_2 = human_model.metabolites.get_by_id(metabolite_id.replace('[' + compartment + ']', '[c]'))
    
    if met_2.id not in all_metabolites:
        err = met_1.id + ' must be in the metabolic model. Attempted to add, but for this, must have the'
        err += 'cytoplasmic version of this metabolite. Please manually add ' + met_2.id + ' to model'
        raise ValueError(err)
    
    all_reactions = list(human_model.reactions)
    r_ = [r for r in all_reactions if met_1 in r.metabolites.keys() and met_2 in r.metabolites.keys()]
    if len(r_) == 0:
        # add transport
        m_name = met_2.id.replace('[c]', '').upper()
        transport = cobra.Reaction(m_name + 't' + compartment)
        transport.name = m_name + ' ' + compartments_[compartment] + ' transport'
        transport.add_metabolites({met_2: -1, met_1: 1})
        transport.lower_bound = -1000
        human_model.add_reactions([transport])
    elif len(r_) == 1:
        r_ = r_[0]
        if not r_.reversibility: # make all exchanges reversible
            r_.lower_bound = -1000
    else:
        reactants, products = list(), list()
        for r in r_:
            reactants += r.reactants
            products += r.products
        if not((met_1 in reactants and met_2 in products) or (met_2 in reactants and met_1 in products)):
            raise ValueError('No reversibility')



def correct_model(model_file = root_path + 'recon2_2.json', 
                 psim_file = root_path + 'psim_recon2_2.csv'):
    '''
    Makes some necessary changes to cobrapy model, largely based on issues encountered with Recon2.2.
    model_file is path/to/cobra_json_model
    psim_file is path/to/psim_csv
    
    '''
    
    human_model = cobra.io.json.load_json_model(model_file)
    
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


    # check for minimum required metabolites
    required_metabolites = pd.read_csv(build_files_path + 'required_metabolic_model_metabolites.csv', index_col = 0)
    all_metabolites = [m.id for m in human_model.metabolites]

    missing_metabolites = sorted(set(required_metabolites.index.tolist()).difference(all_metabolites))
    mm_obj = []
    if len(missing_metabolites) > 0:
        err = 'The input CobraPy model is missing the minimal required metabolites. Manually adding them to the model.'
        warnings.warn(err)
        for mm in missing_metabolites:
            met = cobra.Metabolite(mm)
            info = required_metabolites.loc[mm]

            met.name = info['name']
            met.compartment = info['compartment']
            met.charge = float(info['charge'])
            met.elements = eval(info['elements'])
            met.formula = info['formula']

            mm_obj.append(met)
            all_metabolites.append(mm)
    human_model.add_metabolites(mm_obj) 
        
    print('Adding necessary metabolites and reactions')
    # h2o and pi
    for met_name_ in ['h2o', 'pi']:
        for compartment in ['l', 'm', 'n', 'r', 'x']:
            metabolite_id = met_name_ + '[' + compartment + ']'
            add_required_metabolites(metabolite_id, compartment, human_model, all_metabolites)

    # h
    for compartment in ['l', 'g', 'n', 'r', 'x']: # i and m must be in model for proton gradient appropriate rxns
        metabolite_id = 'h' + '[' + compartment + ']'
        add_required_metabolites(metabolite_id, compartment, human_model, all_metabolites)

    #ppi
    add_required_metabolites('ppi[n]', 'n', human_model, all_metabolites)

    # nucleotides
    compartments1, compartments2 = ['n'], ['l', 'm', 'r', 'x']
    for nucleotide in ['a', 'c', 'u', 'g']:
        for phosphate in ['tp', 'dp', 'mp']:
            if nucleotide == 'u' and phosphate == 'dp':
                compartments = ['g', 'r', 'l']
            elif nucleotide != 'a' or phosphate == 'mp' or nucleotide == 'g':
                compartments = compartments1
            else:
                compartments = compartments1 + compartments2
            for compartment in compartments:
                metabolite_id = nucleotide + phosphate + '[' + compartment + ']'
                add_required_metabolites(metabolite_id, compartment, human_model, all_metabolites)

    compartments = ['n', 'm', 'x', 'r', 'l']
    amino_acids_ = ['ala_L', 'arg_L', 'asn_L', 'asp_L', 'cys_L', 'glu_L', 'gln_L', 'gly', 'his_L', 'ile_L', 'leu_L', 
                  'lys_L', 'met_L', 'phe_L', 'pro_L', 'ser_L', 'thr_L', 'trp_L', 'tyr_L', 'val_L']
    for amino_acid in amino_acids_:
            for compartment in compartments:
                metabolite_id = amino_acid + '[' + compartment + ']'
                add_required_metabolites(metabolite_id, compartment, human_model, all_metabolites)


    print('Check for the recon2.2 HGNC:HGNC error')
    # correct genes with HGNC:HGNC:, recon2.2 has this
    genes_to_format = [g.id for g in human_model.genes if 'HGNC:HGNC:' in g.id]
    # make sure this is not a user formatting thing
    # assuming corrected of just one HGNC version is in PSIM, will check later (correct_format[i] in psim.HGNC_ID)
    genes_to_format = [g for g in genes_to_format if g not in psim_me_genes] 

    if len(genes_to_format) > 0:
        warnings.warn('Your metabolic model contains genes with HGNC:HGNC:####, changing to HGNC:####')
        correct_format = [g[5:] for g in genes_to_format]
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
        raise ValueErorr()

    if len([m.id for m in human_model.metabolites if 'biomass' in m.id]) != 0:
        err_ = 'Extraneous biomass metabolites not associated with the biomass reaction are present,'
        err_ += ' please remove from input cobrapy model'
        raise ValueError(err_)

    rm = sorted(set(metabolites_1).difference([m.id for m in human_model.metabolites]))
    if len([i for i in rm if 'biomass' not in i]) != 0:
        warnings.warn('Non biomass metabolites removed as orphan metabolites from biomass reactions')
    
    cobra.io.save_json_model(human_model, processed_data_path + 'corrected_model.json')
    del human_model
    


# In[58]:


def correct_psim(psim_file = root_path + 'psim_recon2_2.csv'):
    '''
    Runs checks on psim csv.
    psim_file is path/to/psim_csv
    
    '''
    # load files------------------------------------
    psim_me = pd.read_csv(psim_file)
    psim_me.reset_index(inplace = True, drop = True)
    psim_me_genes = psim_me.HGNC_ID.tolist()
    
    essential_cols = ['HGNC_ID', 'PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ']
    optional_cols = ['POLYA_LENGTH', 'TMD', 'SP', 'N_INTRONS', 'DSB', 'GPI', 'OG'] # NG
    other_cols = ['LOCATION']
    all_columns = essential_cols + optional_cols + other_cols

    expression_psim = pd.read_csv(build_files_path + 'expression_module_psim.csv', index_col = 0)
    expression_machinery = sorted(open(build_files_path + 'expression_machinery.txt').read().splitlines())
    expression_psim = expression_psim[all_columns] 
    
    if not os.path.isfile(processed_data_path + 'corrected_model.json'):
        raise ValueError('Run correc_inputs.correct_model befor running this function')
    human_model = cobra.io.load_json_model(processed_data_path + 'corrected_model.json')
    metabolic_machinery = sorted([g.id for g in human_model.genes])
    #------------------------------------
    
    # check that essential cols are all present and appropriately formatted
    if len(set(essential_cols).difference(psim_me.columns.tolist()))>0:
        raise ValueError('Must specify the following columns in PSIM: ' + ', '.join(essential_cols))

    # add optional cols that aren't present
    for col in optional_cols + other_cols:
        if col not in psim_me.columns.tolist():
            warnings.warn(col + ' not in provided PSIM, adding as default values to gene_information')
            psim_me[col] = float('nan')

    # get rid of excess psim columns and put columns in order
    if len(set(psim_me).difference(all_columns)) > 1:
        msg = 'Some of the columns in the provided PSIM are not used in ME model building and are being removed'
        warnings.warn(msg)
    # do regardless of above if statement to put columns in order 
    psim_me = psim_me[all_columns]
    # similarly put expression psim columns in order; this shouldn't be necessary but just in case
    
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


# In[59]:


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
    human_model = cobra.io.load_json_model(processed_data_path + 'corrected_model.json')
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

