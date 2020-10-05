#!/usr/bin/env python
# coding: utf-8

# In[22]:


import cobra
from sympy.parsing.sympy_parser import parse_expr

import sys
sys.path.insert(1, '../../scripts/') # comment out in python script
from utils import parameters as params
from utils import metabolites as metab
from utils import functions as func


# In[23]:


# make the metabolites
biomass_ = cobra.Metabolite('biomass')

# constant
dna_ = cobra.Metabolite('biomass_DNA')
carb_ = cobra.Metabolite('biomass_carbohydrate')
lipid_ = cobra.Metabolite('biomass_lipid')
other_ = cobra.Metabolite('biomass_other')

# variable
protein_, unmodeled_protein_ = cobra.Metabolite('biomass_protein'),  cobra.Metabolite('biomass_unmodeled_protein')
trna_ = cobra.Metabolite('biomass_tRNA')
rrna_ = cobra.Metabolite('biomass_rRNA')
mrna_ = cobra.Metabolite('biomass_mRNA')
premrna_ = cobra.Metabolite('biomass_premRNA')
other_rna_ = cobra.Metabolite('biomass_other_RNA')

type_to_object = {'rrna': rrna_, 'protein': protein_, 'mrna': mrna_, 'trna': trna_}


# In[24]:


# biomass formation reactions
biomass_reactions = list()
biomass_metabolites = [dna_, carb_, lipid_, other_, trna_, rrna_, mrna_, premrna_, other_rna_]
for bm in biomass_metabolites:
    reaction_ = cobra.Reaction(bm.id.split('_')[1] + '_biomass_to_biomass')
    reaction_.add_metabolites({bm: -1, biomass_: 1})
    biomass_reactions.append(reaction_)

# protein biomass with unmodeled protein 
upc=(params.unmodeled_protein_frac)/(1-params.unmodeled_protein_frac)
reaction_ = cobra.Reaction('protein_biomass_to_biomass')
reaction_.add_metabolites({protein_: -1, 
                         unmodeled_protein_: upc, 
                         biomass_: 1 + upc})
biomass_reactions.append(reaction_)


# The following reactions convert the biomass components which are a constant proportion from the metabolic model formulation to the ME model formulation. Briefly, the coefficients of the precursor reactions must be scaled by their molecular weight, and the product must be equal to the constant proportion of that class of biomass, bounded by growth (flux through reaction = growth rate). 

# In[25]:


# constant biomass reactions


#DNA------------------------------------------------------
dna_reaction = cobra.Reaction('DNA_biomass_formation')

# coefs from original RECON2.2
datp_coef = 0.941642857142857
dctp_coef = 0.674428571428572
dgtp_coef = 0.707
dttp_coef = 0.935071428571429

# original coefficient from DNA biomass formation reaction*metabolite molecular weight
rxn = {metab.datp_n: -datp_coef*func.get_metabolite_mw(metab.datp_n),
      metab.dctp_n: -dctp_coef*func.get_metabolite_mw(metab.dctp_n),
      metab.dgtp_n: -dgtp_coef*func.get_metabolite_mw(metab.dgtp_n),
      metab.dttp_n: -dttp_coef*func.get_metabolite_mw(metab.dttp_n),
      dna_: params.dna_frac}
dna_reaction.add_metabolites(rxn)

# CARBOHYDRATE------------------------------------------------------
g6p_coef = 3.87591549295775
carbohydrate_reaction = cobra.Reaction('carbohydrate_biomass_formation')
rxn = {metab.g6p_c: -g6p_coef*func.get_metabolite_mw(metab.g6p_c), 
      carb_: params.carb_frac}
carbohydrate_reaction.add_metabolites(rxn)


# LIPID------------------------------------------------------
chsterol_coef = 0.210319587628866
clpn_hs_coef = 0.120185567010309
pail_hs_coef = 0.240360824742268
pchol_hs_coef = 1.59237113402062
pe_hs_coef = 0.570865979381443
pglyc_hs_coef = 0.0300412371134021
ps_hs_coef = 0.0600927835051546
sphmyln_hs_coef = 0.180268041237113

# lipid_reaction = cobra.Reaction('carbohydrate_biomass_formation')
# rxn = {metab.chsterol_c: -chsterol_coef*func.get_metabolite_mw(metab.chsterol_c),
#        metab.clpn_hs_c: -clpn_hs_coef*1.466059,
#        metab.pail_hs_c: -pail_hs_coef*0.47021,
#        metab.pchol_hs_c: -pchol_hs_coef*func.get_metabolite_mw(metab.pchol_hs_c),
#        metab.pe_hs_c: -pe_hs_coef*func.get_metabolite_mw(metab.pe_hs_c),
#        metab.pglyc_hs_c: -pglyc_hs_coef*func.get_metabolite_mw(metab.pglyc_hs_c),
#        metab.ps_hs_c: -ps_hs_coef*func.get_metabolite_mw(metab.ps_hs_c),
#        metab.sphmyln_hs_c: -sphmyln_hs_coef*func.get_metabolite_mw(metab.sphmyln_hs_c),
#       lipid_: params.lipid_frac}
# lipid_reaction.add_metabolites(rxn)

biomass_reactions += [dna_reaction, carbohydrate_reaction]#, lipid_reaction]

