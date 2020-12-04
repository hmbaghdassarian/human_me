#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra

import sys
sys.path.insert(1, '../../../scripts/') # comment out in python script
from utils import machinery as mach
from utils import parameters as params
from utils import metabolites as metab
from utils import functions as func
from utils import utils_2

from macromolecules.protein import Protein


# In[8]:


def translate_protein_cytosolic(gene_info, mrna_transcript_c, mrna_deg_proxy):    

    # peptide bond formation: https://d1j63owfs0b5j3.cloudfront.net/pop-quiz/answerImage/Amino-Acid-1-popquiz.png
    # tRNA amino acide release: https://rnajournal.cshlp.org/content/14/8/1526/F1.expansion.html

    rxn = dict()
    for aa_code, aa_count in gene_info.amino_acid_counts.items():
        rxn[utils_2.charged_trna_map[aa_code]] = -aa_count # tRNA consumption
    
    rxn[utils_2.modified_trna_transcript_c]  = gene_info.L_protein 
    
    rxn[metab.h2o_c] = -gene_info.L_protein # release of peptide from tRNA, addition of -OH to uncharged tRNA
    rxn[metab.h_c] = gene_info.L_protein # release of peptide from tRNA, addition of -OH to uncharged tRNA
    rxn[metab.h2o_c] += gene_info.L_protein - 1 # peptide bond formation (hydrolysis)
    
    # gtp hydrolysis per aa added
    rxn[metab.ntp_map_c['G']] = -gene_info.L_protein 
    rxn[metab.h2o_c] -= gene_info.L_protein
    rxn[metab.ndp_map_c['G']] = gene_info.L_protein
    rxn[metab.pi_c] = gene_info.L_protein
    rxn[metab.h_c] += gene_info.L_protein
    

    rxn_c = rxn.copy()
    unfolded_protein_c = Protein(compartment = 'c', id_ = 'unfolded', gene_info = gene_info)
    rxn_c[unfolded_protein_c] = 1
    
    # coupling
    rxn_c[mrna_deg_proxy] = -gene_info.coupling['c2']# couple mrna degradation to protein synthesis 
    rxn_c[mrna_transcript_c] = -gene_info.coupling['c1'] # couple mrna dilution to protein synthesis
    
    # biomas    
    translation_elongation = func.ME_Reaction(gene_info.hgnc_id + '_TRANSLATION_ELONGATIONc', 
                                             type_ = ['translation'])
    translation_elongation.subsytem = 'Protein_Expression'
    translation_elongation.add_metabolites(rxn_c)

    translation_elongation.gene_reaction_rule = ' and '.join(mach.translation_efs + ['ribosome']) # GPRs

    return translation_elongation, unfolded_protein_c


# In[9]:


# import random
# import cobra
# import pandas as pd
# from utils import parameters as params
# from utils import functions as func


# psim_toy = pd.DataFrame(columns = ['HGNC_ID', 'PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ', 'POLYA_LENGTH', 'TMD', 
#                                'SP', 'N_INTRONS', 'DSB', 'GPI', 'OG', 'LOCATION'])

# hgnc_id, premrna_seq = 'HGNC:TOY', ''.join(random.choices(['U', 'C', 'G', 'A'], k = 100))
# mrna_seq = premrna_seq[25:75]
# # note that there is no check that the protein_sequence corresponds to the mrna_sequence beyond checking for the length
# protein_seq = ''.join(random.choices(params.amino_acids, k = int(len(mrna_seq)/3)))
# polyA_length, tmd, sp, n_introns, dsb, gpi, og  = None, 1, True, 0, 2, 2, 2
# location = ['c'] # cytoplasm and golgi

# psim_toy.loc[0,:] = [hgnc_id, premrna_seq, mrna_seq, protein_seq, polyA_length, tmd, sp, n_introns, dsb, gpi, og, location]
# from expression.gene_information import gene_information
# gene_info = gene_information(hgnc_id, premrna_seq, mrna_seq, protein_seq,
#                  ptms = {'dsb': dsb, 'og': og, 'gpi': gpi}, tmd = tmd, sp = sp, polyA_length = polyA_length, 
#                  n_introns = n_introns)
# gene_info.get_final_locations(metabolic_model = cobra.Model(''), final_locations = location)
# import expression.build_mrna_expression_reactions as build_mrna
# mrna_reactions, mrna_transcript_c, mrna_deg_proxy = build_mrna.get_mrna_expression_reactions(gene_info)

