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

from uniform_processes import biomass


# In[2]:


def translate_protein_cytosolic(gene_info, mrna_transcript_c, mrna_deg_proxy):    

    # peptide bond formation: https://d1j63owfs0b5j3.cloudfront.net/pop-quiz/answerImage/Amino-Acid-1-popquiz.png
    # tRNA amino acide release: https://rnajournal.cshlp.org/content/14/8/1526/F1.expansion.html

    tb_reactant_bm = 0
    rxn = dict()
    for aa_code, aa_count in gene_info.amino_acid_counts.items():
        rxn[utils_2.charged_trna_map[aa_code]] = -aa_count # tRNA consumption
        tb_reactant_bm += utils_2.charged_trna_map_mw[aa_code]*aa_count # trna biomass consumed kDa
    
    rxn[utils_2.modified_trna_transcript_c]  = gene_info.L_protein 
    biomass_change = (gene_info.L_protein*utils_2.modified_trna_transcript_c_mw)-tb_reactant_bm
    rxn[biomass.trna_] = biomass_change
    
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
    unfolded_protein_c = func.make_protein_metabolite(id_ = gene_info.hgnc_id + '_unfolded', 
                amino_acid_counts = gene_info.amino_acid_counts, L_protein = gene_info.L_protein,
                compartment = 'c')
    rxn_c[unfolded_protein_c] = 1
    
    # coupling
    rxn_c[mrna_deg_proxy] = -gene_info.coupling['c2']# couple mrna degradation to protein synthesis 
    rxn_c[mrna_transcript_c] = -gene_info.coupling['c1'] # couple mrna dilution to protein synthesis
    
    # biomas
    rxn_c[biomass.protein_] = gene_info.protein_mass
    
    translation_elongation = func.ME_Reaction(gene_info.hgnc_id + '_TRANSLATION_ELONGATIONc', 
                                             type_ = ['translation'])
    translation_elongation.subsytem = 'Protein_Expression'
    translation_elongation.add_metabolites(rxn_c)

    translation_elongation.gene_reaction_rule = ' and '.join(mach.translation_efs + ['ribosome']) # GPRs

    return translation_elongation, unfolded_protein_c

def fold_protein_cytosolic(gene_info, unfolded_protein_c):
    # extending proteostasis network in the future would be good
    # will need to make sure inputs to each compartment-specific reactions are at the correct folding stage
    # e.g., mitochondria currently takes unfolded protein, and in future we may want it to take a partially folded
    
    folded_protein_c = unfolded_protein_c.copy()
    folded_protein_c.id = folded_protein_c.id.replace('unfolded', 'folded')
    rxn = {unfolded_protein_c: -1, folded_protein_c: 1}
    protein_folding = cobra.Reaction(gene_info.hgnc_id + '_CYTOSOLIC_PROTEIN_FOLDING')
    protein_folding.subsytem = 'Protein_Expression'
    
    if gene_info.L_protein > 100: #chaperone assisted for larger proteins - https://www.nature.com/articles/nature10317
        rxn = func.hydrolyze_atp(rxn, n_atp = gene_info.L_protein*params.proteolysis_translocation_atp_cost, compartment = 'c')
        protein_folding.gene_reaction_rule = ' and '.join(mach.HSP40_c + mach.HSP70_c) # GPRs
    
    
    protein_folding.add_metabolites(rxn)

    
    
    return protein_folding, folded_protein_c

