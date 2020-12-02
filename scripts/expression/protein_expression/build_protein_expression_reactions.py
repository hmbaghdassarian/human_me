#!/usr/bin/env python
# coding: utf-8

# In[208]:


import cobra

import sys
sys.path.insert(1, '../../../scripts/') # comment out in python script
from utils import machinery as mach
from utils import parameters as params
from utils import metabolites as metab
from utils import functions as func
from utils import utils_2

from macromolecules.protein import Protein

from uniform_processes import biomass
from expression.protein_expression import cytosolic_translation as c_trln


# In[209]:


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


# # Degradation (Ubiquitin-Proteasome)

# In[210]:


def protein_polyubiquitination(gene_info, protein_metabolite, compartment, ub_args):
    
    if compartment == 'c' or compartment == 'n':
        polyu_protein_aa_counts = gene_info.amino_acid_counts.copy()
        for aa_code,aa_counts in ub_args['monoub_aa_counts'].items():
            if aa_code in polyu_protein_aa_counts:
                polyu_protein_aa_counts[aa_code] += aa_counts*params.n_ub
            else: 
                polyu_protein_aa_counts[aa_code] = aa_counts*params.n_ub
    
    if compartment == 'c':
        if protein_metabolite.compartment != 'c':
            raise ValueError('Compartment mismatch for polyubiquitination')
    
        polyubiquitinate_protein = cobra.Reaction(protein_metabolite.id + '_POLYUBIQUITINATIONc')
        polyubiquitinate_protein.subsytem = 'Protein_Expression'

        polyub_protein_c = Protein(id_ = protein_metabolite.id + '_polyub', compartment = 'c',
                           amino_acid_counts = polyu_protein_aa_counts) 
        
        
        biomass_products = polyub_protein_c.mass
        # NOTE: ORIGINAL CODE SUBTRACTED UB_C FROM PROTEIN METABOLITE....I think this was an error
        biomass_substrates = protein_metabolite.mass + (params.n_ub*ub_args['ub_c'].mass) 
        biomass_change = biomass_products - biomass_substrates
        
        
        rxn = {protein_metabolite: -1, ub_args['ub_c']: -params.n_ub, polyub_protein_c:1, metab.h2o_c: params.n_ub,
              biomass.protein_: biomass_change}
        # 1 ATP hydrolysis per ubiquitin monomer added (https://link.springer.com/article/10.1007/s10637-020-00894-6)
        rxn = func.hydrolyze_atp(rxn, n_atp = params.n_ub, compartment = 'c')

        polyubiquitinate_protein.add_metabolites(rxn)
        polyubiquitinate_protein.gene_reaction_rule = ' and '.join(mach.UB_ligases_c)
        return polyubiquitinate_protein, polyub_protein_c
    elif compartment == 'n':
        if protein_metabolite.compartment != 'n':
            raise ValueError('Compartment mismatch for polyubiquitination')

        polyubiquitinate_protein = cobra.Reaction(protein_metabolite.id + '_POLYUBIQUITINATIONn')
        polyubiquitinate_protein.subsytem = 'Protein_Expression'

        polyub_protein_n = Protein(id_ = protein_metabolite.id + '_polyub', compartment = 'n',
                           amino_acid_counts = polyu_protein_aa_counts) 
        
        biomass_products = polyub_protein_n.mass 
        # NOTE ORIGINAL CODE substracted ub_n weight, i think this was an error
        biomass_substrates = (protein_metabolite.mass + (params.n_ub*ub_args['ub_n'].mass))        
        biomass_change = biomass_products - biomass_substrates
        
    
        rxn = {protein_metabolite: -1, ub_args['ub_n']: -params.n_ub, polyub_protein_n: 1, metab.h2o_n: params.n_ub, 
              biomass.protein_: biomass_change}
        # 1 ATP hydrolysis per ubiquitin monomer added (https://link.springer.com/article/10.1007/s10637-020-00894-6)
        rxn = func.hydrolyze_atp(rxn, n_atp = params.n_ub, compartment = 'n')


        polyubiquitinate_protein.add_metabolites(rxn)
        polyubiquitinate_protein.gene_reaction_rule = ' and '.join(mach.UB_ligases_n)
        
        return polyubiquitinate_protein, polyub_protein_n
    elif compartment == 'pm':
        if protein_metabolite.compartment != 'pm':
            raise ValueError('Compartment mismatch for polyubiquitination')
    
        polyubiquitinate_protein = cobra.Reaction(protein_metabolite.id + '_POLYUBIQUITINATIONpm')
        polyubiquitinate_protein.subsytem = 'Protein_Expression'
        
        polyub_protein_pm = protein_metabolite.change_compartment('pm')
        
        elements = polyub_protein_pm.elements.copy()
        for aa_code, aa_count in ub_args['monoub_aa_counts'].items():
            aa_elements = metab.seq_amino_acid_map_c[aa_code].elements
            for element in aa_elements:
                elements[element] += aa_count*aa_elements[element]*params.n_ub
            polyub_protein_pm.charge += metab.seq_amino_acid_map_c[aa_code].charge*aa_count
        # peptide bond formation
        elements['H'] -= 2*(ub_args['L_monoub']*params.n_ub) # no -1 bc already accounted for in copying elements
        elements['O'] -= 1*(ub_args['L_monoub']*params.n_ub)
        polyub_protein_pm.elements = elements
        
        biomass_products = polyub_protein_pm.mass 
        biomass_substrates = protein_metabolite.mass + (params.n_ub*ub_args['ub_c'].mass)
        biomass_change = (biomass_products - biomass_substrates)

        rxn = {protein_metabolite: -1, ub_args['ub_c']: -params.n_ub, polyub_protein_pm: 1, metab.h2o_c: params.n_ub, 
              biomass.protein_: biomass_change}
        # 1 ATP hydrolysis per ubiquitin monomer added (https://link.springer.com/article/10.1007/s10637-020-00894-6)
        rxn = func.hydrolyze_atp(rxn, n_atp = params.n_ub, compartment = 'c')

        polyubiquitinate_protein.add_metabolites(rxn)
        polyubiquitinate_protein.gene_reaction_rule = ' and '.join(mach.UB_ligases_c + mach.HSP70_c + mach.HSP90AB1)
        
        return polyubiquitinate_protein, polyub_protein_pm
    
    else:
        raise ValueError('Current compartment does not have polyubiquitination')

def proteasomal_degradation(gene_info, protein_metabolite, polyub_protein_metabolite, compartment, 
                            ub_args):
    if compartment == 'c':
        if (protein_metabolite.compartment != 'c') or (polyub_protein_metabolite.compartment != 'c'):
            raise ValueError('Compartment mismatch for cytoplasmic proteasomal degradation')
        
        
        deubiquitination = cobra.Reaction(protein_metabolite.id + '_DEUBIQUITINATIONc')
        deubiquitination.subsytem = 'Protein_Expression'
        
        biomass_products = protein_metabolite.mass + ub_args['polyub_c'].mass
        biomass_substrates = polyub_protein_metabolite.mass 
        biomass_change = (biomass_products - biomass_substrates)
        
        deubiquitination.add_metabolites({polyub_protein_metabolite: -1, metab.h2o_c: -1, 
                                          protein_metabolite: 1, ub_args['polyub_c']: 1, biomass.protein_: biomass_change})
        deubiquitination.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)

        protein_degradation = cobra.Reaction(protein_metabolite.id + '_PROTEASOMAL_DEGRADATIONc')
        protein_degradation.subsytem = 'Protein_Expression'
        rxn = {metab.seq_amino_acid_map_c[aa_code]: aa_counts for aa_code, aa_counts in gene_info.amino_acid_counts.items()}
        rxn[polyub_protein_metabolite], rxn[metab.h2o_c], rxn[ub_args['polyub_c']] = -1, -gene_info.L_protein, 1
        rxn[biomass.protein_] = (ub_args['polyub_c'].mass - polyub_protein_metabolite.mass)
        
        # atp hydrolysis for translocation/unfolding  - known 1 ATP per 2 residues - https://www.nature.com/articles/s41586-018-0736-4
        # L_polub_protein = (gene_info.L_protein + (ub_args['L_monoub']*params.n_ub)) 
        rxn = func.hydrolyze_atp(rxn, n_atp = gene_info.L_protein/2, compartment = 'c')


        protein_degradation.add_metabolites(rxn)
        protein_degradation.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)

        protein_degradation_reactions = [deubiquitination, protein_degradation]

        return protein_degradation_reactions
    elif compartment == 'n':
        if (protein_metabolite.compartment != 'n') or (polyub_protein_metabolite.compartment != 'n'):
            raise ValueError('Compartment mismatch for nuclear proteasomal degradation')


        deubiquitination = cobra.Reaction(protein_metabolite.id + '_DEUBIQUITINATIONn')
        deubiquitination.subsytem = 'Protein_Expression'
        
        biomass_products = protein_metabolite.mass + ub_args['polyub_n'].mass
        biomass_substrates = polyub_protein_metabolite.mass
        biomass_change = (biomass_products - biomass_substrates)
        
        
        deubiquitination.add_metabolites({polyub_protein_metabolite: -1, metab.h2o_n: -1, protein_metabolite: 1, ub_args['polyub_n']: 1, 
                                         biomass.protein_: biomass_change})
        deubiquitination.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)

        protein_degradation = cobra.Reaction(protein_metabolite.id + '_PROTEASOMAL_DEGRADATIONn')
        protein_degradation.subsytem = 'Protein_Expression'
        rxn = {metab.seq_amino_acid_map_n[aa_code]: aa_counts for aa_code, aa_counts in gene_info.amino_acid_counts.items()}
        rxn[polyub_protein_metabolite], rxn[metab.h2o_n], rxn[ub_args['polyub_n']] = -1, -gene_info.L_protein, 1
        rxn[biomass.protein_] = (ub_args['polyub_n'].mass - polyub_protein_metabolite.mass)

        # atp hydrolysis for translocation/unfolding  - known 1 ATP per 2 residues - https://www.nature.com/articles/s41586-018-0736-4
#         L_polub_protein = (gene_info.L_protein + (ub_args['L_monoub']*params.n_ub)) 
        rxn = func.hydrolyze_atp(rxn, n_atp = gene_info.L_protein/2, compartment = 'n')


        protein_degradation.add_metabolites(rxn)
        protein_degradation.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)

        protein_degradation_reactions = [deubiquitination, protein_degradation]

        return protein_degradation_reactions
    else:
        raise ValueError('Current compartment does not have proteasomal degradation')


# # Cytosolic Degradation

# In[211]:


def degrade_cytosolic_protein(gene_info, folded_protein_c, ub_args):
    polyubiquitinate_folded_protein_c, polyub_protein_c = protein_polyubiquitination(gene_info, 
                                                          protein_metabolite = folded_protein_c, compartment = 'c', 
                                                                                    ub_args = ub_args) 
    cytosolic_proteasomal_degradation_reactions = proteasomal_degradation(gene_info, 
                                                  protein_metabolite = folded_protein_c, 
                                                  polyub_protein_metabolite = polyub_protein_c, compartment = 'c', 
                                                  ub_args = ub_args)
    
    return [polyubiquitinate_folded_protein_c] + cytosolic_proteasomal_degradation_reactions


# # Nuclear Reactions

# In[212]:


def transport_nuclear_protein(gene_info, folded_protein_c):

    folded_protein_n = folded_protein_c.change_compartment('n')

    nuclear_import = cobra.Reaction(gene_info.hgnc_id + '_IMPORTtn')
    nuclear_import.subsytem = 'Protein_Expression'
#     nuclear_export = nuclear_import.copy()
#     nuclear_export.id = nuclear_export.id.replace('IMPORT', 'EXPORT')

    import_rxn = {folded_protein_c: -1, folded_protein_n: 1}
#     export_rxn = {folded_protein_c: 1, folded_protein_n: -1}

    if folded_protein_c.mass > params.nuclear_diffusion_limit:
        # gtp hydrolysis per import
        import_rxn[metab.gtp_n], import_rxn[metab.h2o_n], import_rxn[metab.gdp_n], import_rxn[metab.pi_n], import_rxn[metab.h_n]  = -1, -1, 1, 1, 1
        nuclear_import.add_metabolites(import_rxn)
        nuclear_import.gene_reaction_rule = ' and '.join(mach.importins + mach.RAN)
#         export_rxn[gtp_c], export_rxn[h2o_c], export_rxn[gdp_c], export_rxn[pi_c], export_rxn[h_c]  = -1, -1, 1, 1, 1
#         nuclear_export.add_metabolites(export_rxn)
#         nuclear_export.gene_reaction_rule = ' and '.join(XPO1 + RAN)

    else: # diffusion
        nuclear_import.add_metabolites(import_rxn)
        nuclear_import.lower_bound = -1000 # reversible
        
    return nuclear_import, folded_protein_n

def degrade_nuclear_protein(gene_info, folded_protein_n, ub_args):
    polyubiquitinate_folded_protein_n, polyub_protein_n = protein_polyubiquitination(gene_info, 
                                                          protein_metabolite = folded_protein_n, compartment = 'n', 
                                                          ub_args = ub_args)
    nuclear_proteasomal_degradation_reactions = proteasomal_degradation(gene_info, 
                                                protein_metabolite = folded_protein_n, 
                                                polyub_protein_metabolite = polyub_protein_n, compartment = 'n', 
                                                ub_args = ub_args)
    return [polyubiquitinate_folded_protein_n] + nuclear_proteasomal_degradation_reactions

def get_nuclear_reactions(gene_info, folded_protein_c, ub_args):
    nuclear_import, folded_protein_n = transport_nuclear_protein(gene_info, folded_protein_c)
    nuclear_degradation_reactions = degrade_nuclear_protein(gene_info, folded_protein_n, 
                                    ub_args = ub_args)
    
    return [nuclear_import] + nuclear_degradation_reactions, folded_protein_n


# # Mitochondrial Reactions
# 

# In[213]:


# i is intermembrane space, but called inner in compartments BIGG
# stick to notation and use inner instead of inter in reaction naming

def transport_mitochondrial_matrix(gene_info, unfolded_protein_c):
    # transport and folding
    if unfolded_protein_c.compartment != 'c':
        raise ValueError('Only cytoplasmic proteins can be transported to mitochondrial matrix')
    
    mitochondrial_matrix_transport = cobra.Reaction(gene_info.hgnc_id + '_IMPORTtm')
    mitochondrial_matrix_transport.subsytem = 'Protein_Expression'
    pre_protein_m = unfolded_protein_c.change_compartment('m')
    pre_protein_m.id = pre_protein_m.id.replace('unfolded', 'folded_pre')
    
    rxn = {unfolded_protein_c: -1, pre_protein_m: 1}
    # ATP hydrolysis for transport, assums 1 ATP consumed per 2 residues
    rxn = func.hydrolyze_atp(rxn, n_atp = gene_info.L_protein*params.transport_translocation_atp_cost, compartment = 'm')
    
    
    mitochondrial_matrix_transport.add_metabolites(rxn)
    mitochondrial_matrix_transportgene_reaction_rule = ' and '.join(mach.TOM + mach.TIM23_PAM + mach.HSP70_m)
    
    return mitochondrial_matrix_transport, pre_protein_m

def mitochondrial_matrix_protein_processing(gene_info, pre_protein_m):
    # implement this in the future: cleavage of MTS (and degradation of MTS)
    
    # add biomass changes in the future
    processed_protein_m, aa_counts_processed_m, L_processed_protein_m = pre_protein_m, gene_info.amino_acid_counts.copy(), gene_info.L_protein
    process_mitochondrial_matrix_protein = None
    return process_mitochondrial_matrix_protein, processed_protein_m, aa_counts_processed_m, L_processed_protein_m


def degrade_mitochondrial_protein(gene_info, protein_metabolite, compartment, L_protein, amino_acid_counts):
    rxn = {metab.seq_amino_acid_map_m[aa_code]: aa_counts for aa_code, aa_counts in amino_acid_counts.items()}
    rxn[protein_metabolite], rxn[metab.h2o_m] = -1, -(L_protein-1)
    rxn[biomass.protein_] = -protein_metabolite.mass
    
    if compartment == 'm':
        mitochondrial_degradation = cobra.Reaction(gene_info.hgnc_id + '_DEGRADATIONm')
        mitochondrial_degradation.gene_reaction_rule = mach.mLON[0]
        
        # ATP hydrolysis by LON: 2 ATP per residue - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518814/
        rxn = func.hydrolyze_atp(rxn, n_atp = L_protein*2, compartment = 'm')
        

    elif compartment == 'i':
        mitochondrial_degradation = cobra.Reaction(gene_info.hgnc_id + '_DEGRADATIONi')
        mitochondrial_degradation.gene_reaction_rule = mach.iAAA[0]#' and '.join(mAAA + iAAA)
        # in the future, may want to add ubqituin-proteasome: 
        
        # ATP hydrolysis by m/i-AAA: 1 ATP per 2 residues -- no source, assumes same as 26S proteasome
        rxn = func.hydrolyze_atp(rxn, n_atp = L_protein*params.proteolysis_translocation_atp_cost, compartment = 'i')
  
    mitochondrial_degradation.subsytem = 'Protein_Expression'
    mitochondrial_degradation.add_metabolites(rxn)
    
    return mitochondrial_degradation

def transport_mitochondrial_inter(gene_info, processed_protein_m):
    # upper left Fig 12-29 https://www.ncbi.nlm.nih.gov/books/NBK26828/ 
    # import to matrix then re-export to inter membrane space
    
    if processed_protein_m.compartment != 'm':
        raise ValueError('Only the mechanism of mitochondrial matrix import and re-export to inter membrane is considered')
    
    mitochondrial_inter_transport = cobra.Reaction(gene_info.hgnc_id + '_IMPORTti')
    mitochondrial_inter_transport.subsytem = 'Protein_Expression'
    pre_protein_i = processed_protein_m.change_compartment('i')
    
    rxn = {processed_protein_m: -1, pre_protein_i: 1}
    
    mitochondrial_inter_transport.add_metabolites(rxn)
    mitochondrial_inter_transportgene_reaction_rule = mach.OXA[0]
    
    return mitochondrial_inter_transport, pre_protein_i

def mitochondrial_inter_protein_processing(gene_info, pre_protein_i):
    # implement this in the future: cleavage of secondary sequence (and degradation)
    
    # add biomass reactions if actually processing
    processed_protein_i, aa_counts_processed_i, L_processed_protein_i = pre_protein_i, gene_info.amino_acid_counts.copy(), gene_info.L_protein
    process_mitochondrial_matrix_protein = None
    return process_mitochondrial_matrix_protein, processed_protein_i, aa_counts_processed_i, L_processed_protein_i

def get_mitochondrial_reactions(gene_info, unfolded_protein_c, compartments):
    mitochondrial_matrix_transport, pre_protein_m = transport_mitochondrial_matrix(gene_info, unfolded_protein_c)
    process_mitochondrial_matrix_protein, processed_protein_m, aa_counts_processed_m, L_processed_protein_m = mitochondrial_matrix_protein_processing(gene_info, pre_protein_m)
    
    mitochondrial_reactions = [mitochondrial_matrix_transport]
    mitochondrial_protein_metabolites = list()
    if process_mitochondrial_matrix_protein != None:
        mitochondrial_reactions += [process_mitochondrial_matrix_protein]
    
    if 'm' in compartments:
        mitochondrial_matrix_degradation = degrade_mitochondrial_protein(gene_info, protein_metabolite = processed_protein_m, compartment = 'm', L_protein = L_processed_protein_m, amino_acid_counts = aa_counts_processed_m)
        mitochondrial_reactions += [mitochondrial_matrix_degradation]
        mitochondrial_protein_metabolites += [processed_protein_m]
    if 'i' in compartments:
        mitochondrial_inter_transport, pre_protein_i = transport_mitochondrial_inter(gene_info, processed_protein_m)
        process_mitochondrial_inter_protein, processed_protein_i, aa_counts_processed_i, L_processed_protein_i = mitochondrial_inter_protein_processing(gene_info, pre_protein_i)
        if process_mitochondrial_matrix_protein != None:
                mitochondrial_reactions += [process_mitochondrial_inter_protein]        
        mitochondrial_inter_degradation = degrade_mitochondrial_protein(gene_info, protein_metabolite = processed_protein_i, compartment = 'i', L_protein = L_processed_protein_i, amino_acid_counts = aa_counts_processed_i)
        mitochondrial_reactions += [mitochondrial_inter_transport, mitochondrial_inter_degradation]
        mitochondrial_protein_metabolites += [processed_protein_i]

    return mitochondrial_reactions, mitochondrial_protein_metabolites


# # Peroxisomal

# In[214]:


def transport_peroxisome(gene_info, folded_protein_c):
    if folded_protein_c.compartment != 'c':
        raise ValueError('Only cytoplasmic proteins can be transported to mitochondrial matrix')
    
    peroxisomal_transport = cobra.Reaction(gene_info.hgnc_id + '_IMPORTtx')
    peroxisomal_transport.subsytem = 'Protein_Expression'
    folded_protein_x = folded_protein_c.change_compartment('x')
    
    rxn = {folded_protein_c: -1, folded_protein_x: 1}
    # ATP hydrolysis for transport--translocation of protein, export of PEX5S receptor
    rxn = func.hydrolyze_atp(rxn, n_atp = (gene_info.L_protein+mach.L_PEX5)*params.transport_translocation_atp_cost, 
                        compartment = 'x')
    
    
    peroxisomal_transport.add_metabolites(rxn)
    peroxisomal_transport.gene_reaction_rule = ' and '.join(mach.peroxins + mach.AWP1)
    
    return peroxisomal_transport, folded_protein_x

def degrade_peroxisomal_protein(gene_info, folded_protein_x):
    
    peroxisomal_degradation = cobra.Reaction(gene_info.hgnc_id + '_DEGRADATIONx')
    peroxisomal_degradation.subsytem = 'Protein_Expression'
    peroxisomal_degradation.gene_reaction_rule = mach.LONP2[0]

    rxn = {metab.seq_amino_acid_map_x[aa_code]: aa_counts for aa_code, aa_counts in gene_info.amino_acid_counts.items()}
    rxn[folded_protein_x], rxn[metab.h2o_x] = -1, -(gene_info.L_protein-1)
    rxn[biomass.protein_] = -folded_protein_x.mass
    # ATP hydrolysis by LON: 2 ATP per residue - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518814/
    rxn = func.hydrolyze_atp(rxn, n_atp = gene_info.L_protein*2, compartment = 'x')
    peroxisomal_degradation.add_metabolites(rxn)
    
    return peroxisomal_degradation

def get_peroxisomal_reactions(gene_info, folded_protein_c):
    peroxisomal_transport, folded_protein_x = transport_peroxisome(gene_info, folded_protein_c)
    peroxisomal_degradation = degrade_peroxisomal_protein(gene_info, folded_protein_x)
    
    return [peroxisomal_transport, peroxisomal_degradation], folded_protein_x


# # Secretory Pathway
# 
# adapted from Jahir's Recon2_2s

# # ER transport

# In[215]:


def post_translational_translocation(gene_info, unfolded_protein_c):
    if gene_info.L_protein > params.ptt_length:
        raise ValueError('This protein is too long for post-translational translocation')
    if unfolded_protein_c.compartment != 'c':
        raise ValueError('Protein metabolite is not in cytosolic compartment')
    
    ptt_reactions = list()
    
    folded_protein_r = unfolded_protein_c.change_compartment('r')
    folded_protein_r.id = folded_protein_r.id.replace('unfolded', 'folded')
    folded_protein_r.compartment = 'r'
    rxn = {unfolded_protein_c: -1, folded_protein_r: 1}
    rxn = func.hydrolyze_atp(rxn, n_atp = 1, compartment = 'c')
    
    if gene_info.tmd > 0 or 'pm' in gene_info.final_locations.keys(): # membrane secreted protein
        post_translational_translocation_r = cobra.Reaction(gene_info.hgnc_id + '_post_TRANSLOC_3A_IMPORTtr')
        post_translational_translocation_r.subsytem = 'Protein_Expression'
        post_translational_translocation_r.gene_reaction_rule = ' and '.join(mach.ASNA1 + mach.WRB + mach.translation_efs + ['ribosome'])
        
        # complex cleavage from jahir's (+h2o_c, +h_c, + pi_c) not included 
        post_translational_translocation_r.add_metabolites(rxn)
        ptt_reactions += [post_translational_translocation_r]
     
    else: #non membrane secreted protein
        number_BiP = gene_info.L_protein/40
        post_translational_translocation_r = cobra.Reaction(gene_info.hgnc_id + '_post_TRANSLOC_3B_IMPORTtr')
        post_translational_translocation_r.subsytem = 'Protein_Expression'
        post_translational_translocation_r.gene_reaction_rule = ' and '.join(mach.ptnm+mach.translation_efs + ['ribosome'])
        
        rxn = func.hydrolyze_atp(rxn, n_atp = number_BiP, compartment = 'r')
        post_translational_translocation_r.add_metabolites(rxn)
        ptt_reactions += [post_translational_translocation_r]

    return ptt_reactions, folded_protein_r

def co_translational_translocation(gene_info, mrna_transcript_c, mrna_deg_proxy):
    if gene_info.L_protein <= params.ptt_length:
        raise ValueError('This protein is too short for co-translational translocation')    
    
    ctt_reactions = list()
    
    # reaction metabolites------------------------------------------------------------------------------------
    number_BiP = gene_info.L_protein/40
    
#     rxn = {utils_2.charged_trna_map[aa_code]: -aa_count for aa_code, aa_count in gene_info.amino_acid_counts.items()} # tRNA consumption
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
    unprocessed_protein_r = Protein(compartment = 'r', id_ = 'unprocessed_folded', gene_info = gene_info)
    
    mw_upr = unprocessed_protein_r.mass
    rxn[unprocessed_protein_r], rxn[biomass.protein_] = 1, mw_upr
    rxn = func.hydrolyze_atp(rxn, n_atp = number_BiP, compartment = 'r')
    
    # coupling
    rxn[mrna_deg_proxy] = -gene_info.coupling['c2']# couple mrna degradation to protein synthesis 
    rxn[mrna_transcript_c] = -gene_info.coupling['c1'] # couple mrna dilution to protein synthesis


    #------------------------------------------------------------------------------------

    co_translational_translocation_r = func.ME_Reaction(gene_info.hgnc_id + '_co_TRANSLOC_IMPORTtr', 
                                                       type_ = ['translation'])
    co_translational_translocation_r.subsytem = 'Protein_Expression'
    co_translational_translocation_r.add_metabolites(rxn)
    co_translational_translocation_r.gene_reaction_rule = ' and '.join(mach.ctnm + mach.translation_efs + ['ribosome'])
    ctt_reactions += [co_translational_translocation_r]
    
    # sp degradation
    sp_seq = gene_info.protein_seq[:params.L_sp]
    sp_aa_counts = {k: sp_seq.count(k) for k in params.amino_acids}
    gene_info.protein_seq = gene_info.protein_seq[params.L_sp:]
    gene_info.amino_acid_counts = {k: gene_info.protein_seq.count(k) for k in params.amino_acids}
    gene_info.L_protein = len(gene_info.protein_seq)

    folded_protein_r = Protein(compartment = 'r', id_ = 'folded', gene_info = gene_info)

    rxn = {metab.seq_amino_acid_map_r[aa]: count for aa, count in sp_aa_counts.items()}
    rxn[metab.h2o_r] = -params.L_sp
    rxn[unprocessed_protein_r],rxn[folded_protein_r] = -1, 1
    rxn[biomass.protein_] = (folded_protein_r.mass) - mw_upr
    
    sp_degradation = cobra.Reaction(gene_info.hgnc_id + '_SP_degradationr')
    sp_degradation.subsystem = 'Protein Expression'
    sp_degradation.add_metabolites(rxn)
    sp_degradation.gene_reaction_rule = mach.sp_rule
    ctt_reactions += [sp_degradation]

    return ctt_reactions, folded_protein_r, gene_info


# # ER Modifications

# In[239]:


def form_disulfide_bond(gene_info, folded_protein_r):
    number_DSB = gene_info.ptms['dsb']
    disulfide_bond_formation = cobra.Reaction(gene_info.hgnc_id + '_DSBr')
    disulfide_bond_formation.subsystem = 'Protein Expression'
    modified_protein_dsb_r = folded_protein_r.copy()
    modified_protein_dsb_r.id = modified_protein_dsb_r.id.replace('folded', 'folded_DSB')
    elements = folded_protein_r.elements.copy()
    elements['H'] -= 2*number_DSB
    modified_protein_dsb_r.elements = elements
    # diagram https://www.google.com/url?sa=i&url=https%3A%2F%2Fen.wikipedia.org%2Fwiki%2FProtein_disulfide-isomerase&psig=AOvVaw0bGpff4XX1eYEF61H1RJKw&ust=1597273135069000&source=images&cd=vfe&ved=0CAIQjRxqFwoTCJi6l6GglOsCFQAAAAAdAAAAABAJ
    # incorporate exchange with PDI in future versions
    rxn = {folded_protein_r: -1, modified_protein_dsb_r: 1, metab.o2_r: -number_DSB, metab.h2o2_r: number_DSB, 
          biomass.protein_: (modified_protein_dsb_r.mass - folded_protein_r.mass)}
    disulfide_bond_formation.add_metabolites(rxn)
    disulfide_bond_formation.gene_reaction_rule = mach.P4HB[0]
    
    return disulfide_bond_formation, modified_protein_dsb_r

def form_gpi(gene_info, modified_protein_r):
    gpi_formation = cobra.Reaction(gene_info.hgnc_id + '_GPIr')
    gpi_formation.subsystem = 'Protein Expression'
    modified_protein_gpi_r = modified_protein_r.copy()
    modified_protein_gpi_r.id = modified_protein_gpi_r.id.replace('folded', 'folded_GPI')
    
#     # need to figure this out correctly!!
#     elements = modified_protein_r.elements.copy()
#     for e,c in metab.balanced_gpi.items():
#         if e in elements.keys():
#             elements[e] += c
#         else:
#             elements[e] = c
            
#     modified_protein_gpi_r.elements = elements

#     rxn = metab.M4ATAer.copy() # need these additional metabolties to get mass balance with gpi_sig[r]
    rxn = dict()
    rxn[metab.hdca_r], rxn[metab.gpi_hs_r], rxn[metab.h_r], rxn[metab.h2o_r] = 1,-1,1,-1
    rxn[modified_protein_r], rxn[modified_protein_gpi_r]= -1, 1
#     rxn[biomass.protein_] = func.get_metabolite_mw(modified_protein_gpi_r) - func.get_metabolite_mw(modified_protein_r)

    gpi_formation.add_metabolites(rxn)
    gpi_formation.gene_reaction_rule = ' and '.join(mach.gpi_machinery)
    
    raise ValueError('This reaction is not balanced, must fix')
    
    return gpi_formation, modified_protein_gpi_r

def glycosylate_n_linked(gene_info, modified_protein_r):
    raise ValueError('N-glycosylation not yet incorporated')
#     n_glycosylation = cobra.Reaction(gene_info.hgnc_id + 'NGLYCOr')
#     n_glycosylation.subsystem = 'Protein Expression'
#     modified_protein_ng_r = modified_protein_r.copy()
#     modified_protein_ng_r.id = modified_protein_ng_r.id.replace('folded', 'folded_NG')

    # # add metabolites and GPRS
    # return n_glycosylation, modified_protein_ng_r
    


def modify_protein_er(gene_info, folded_protein_r):
    if folded_protein_r.compartment != 'r':
        raise ValueError('Only er compartment proteins can get disulfide bonds, GPI anchors, or n glycosylation')
    
    
    modification_reactions = list()
    
    if 'dsb' in gene_info.ptms.keys() and gene_info.ptms['dsb'] > 0:
        disulfide_bond_formation, modified_protein_r = form_disulfide_bond(gene_info, folded_protein_r)
        modification_reactions += [disulfide_bond_formation]
    else:
        modified_protein_r = folded_protein_r # these lines update the protein metabolite to be appropriate inputs to proceeding functions

    if 'gpi' in gene_info.ptms.keys() and gene_info.ptms['gpi'] > 0: # == 1 but doesn't matter bc of gene_info checks
        
        gpi_formation, modified_protein_r = form_gpi(gene_info, modified_protein_r)
        modification_reactions += [gpi_formation]
    else:
        modified_protein_r = modified_protein_r
    
    # N GLYCOSYLATION here must be updated   
    if 'ng' in gene_info.ptms.keys() and gene_info.ptms['ng'] > 0: 
        n_glycosylation, modified_protein_r = glycosylate_n_linked(gene_info, modified_protein_r)
        modification_reactions += [n_glycosylation]
    else:
        modified_protein_r = modified_protein_r
    
    
    return modification_reactions, modified_protein_r


# # Golgi Reactions

# In[217]:


def import_golgi(gene_info, modified_protein_r):
    V = modified_protein_r.mass * 1.21 / 1000.0 # Protein Volume in nm^3
    copii_coeff = int(round(268082.35 * params.Kv / V))
    
    protein_g = modified_protein_r.change_compartment('g')
    
    rxn = {modified_protein_r: -copii_coeff, protein_g: copii_coeff}
    # gtp hydrolysis
    rxn[metab.ntp_map_c['G']], rxn[metab.h2o_c], rxn[metab.ndp_map_c['G']], rxn[metab.pi_c], rxn[metab.h_c]  = -94, -94, 94, 94, 94

    golgi_import = cobra.Reaction(gene_info.hgnc_id + '_COPII_IMPORTtg')
    golgi_import.subsystem = 'Protein Expression'
    golgi_import.add_metabolites(rxn)
    
    
    if 'gpi' in gene_info.ptms.keys() and 'ng' not in gene_info.ptms.keys(): # this if statement is analogous to Recon2.2S's connector statements in copii reactions
        golgi_import.gene_reaction_rule = ' and '.join(mach.copii_gpi_m)
    else:
        golgi_import.gene_reaction_rule = ' and '.join(mach.copii_r_m)
        
    
    return golgi_import, protein_g

def glycosylate_o_linked(gene_info, protein_g):
    number_Oglycans = gene_info.ptms['og']
    o_glycosylation = cobra.Reaction(gene_info.hgnc_id + '_OGg')
    o_glycosylation.subsystem = 'Protein Expression'

    # metabolites
    modified_protein_og_g = protein_g.copy()
    modified_protein_og_g.id = modified_protein_og_g.id.replace('folded', 'folded_OG')

    balance_og = {'C': (8 + 6 + 8)*number_Oglycans, # each 1/3 entry is for each 1/3 reactions in Jahir's model, in case want to separate in the future 
                  'H': (13 + 10 + 13)*number_Oglycans, 
                  'N': (1 + 0 + 1)*number_Oglycans, 
                  'O': (5 + 5 + 5)*number_Oglycans}
    elements = modified_protein_og_g.elements.copy()
    for e,c in balance_og.items():
        if e in elements.keys():
            elements[e] += c
        else:
            elements[e] = c
    modified_protein_og_g.elements = elements

    rxn = {protein_g: -1, modified_protein_og_g: 1, metab.udpacgal_g: -number_Oglycans, 
           metab.udpgal_g: -number_Oglycans, metab.uacgam_g: -number_Oglycans, metab.h_g: 3*number_Oglycans, 
           metab.udp_g: 3* number_Oglycans,
           biomass.protein_: (protein_g.mass - modified_protein_og_g.mass)}
    
    
    o_glycosylation.add_metabolites(rxn)
    
    o_glycosylation.gene_reaction_rule = mach.og_rule
    
    return o_glycosylation, modified_protein_og_g
    

def modify_protein_golgi(gene_info, protein_g):
    if protein_g.compartment != 'g':
        raise ValueError('Only golgi compartment proteins can be O-glycosylated')
    
    # this set up allows for incorporation of other Golgi PTMs in the future, similar to modify_protein_er fct
    modification_reactions = list()
    if 'og' in gene_info.ptms.keys() and gene_info.ptms['og'] > 0:
        o_glycosylation, modified_protein_g = glycosylate_o_linked(gene_info, protein_g)
        modification_reactions += [o_glycosylation]
    else:
        modified_protein_g = protein_g

    return modification_reactions, modified_protein_g


def retrograde_er(gene_info, modified_protein_g):
    V = modified_protein_g.mass * 1.21 / 1000.0 # Protein Volume in nm^3
    copi_coeff = int(round(143793.19 * params.Kv / V))

    retro_protein_r = modified_protein_g.change_compartment('r')

    rxn = {modified_protein_g: -copi_coeff, retro_protein_r: copi_coeff}
    # gtp hydrolysis
    rxn[metab.ntp_map_c['G']], rxn[metab.h2o_c], rxn[metab.ndp_map_c['G']], rxn[metab.pi_c], rxn[metab.h_c]  = -127, -127, 127, 127, 127

    retrograde_transport = cobra.Reaction(gene_info.hgnc_id + '_COPI_RETROtr')
    retrograde_transport.subsystem = 'Protein Expression'
    retrograde_transport.add_metabolites(rxn)
    retrograde_transport.gene_reaction_rule = ' and '.join(mach.copi_m)
    
    return retrograde_transport, retro_protein_r


# # Lysosomal, Extracellular, and Plasma Membrane Transport

# In[218]:


def secrete_protein(gene_info, modified_protein_g):
    V = modified_protein_g.mass * 1.21 / 1000.0 # Protein Volume in nm^3
    clathrin_coeff = int(round(29880.01 * params.Kv / V)) # Number of proteins per clathrin vesicle  

    secreted_proteins = list()
    for nc in ['e', 'pm', 'l']:
        if nc in gene_info.final_locations.keys():
            secreted_protein = modified_protein_g.change_compartment(nc)
            secreted_proteins += [secreted_protein]
    
#     statement1 = 'pm' in gene_info.final_locations.keys()
#     lysosomal_degradation_ptm_condition = 'gpi' in gene_info.ptms.keys() and len(gene_info.ptms.keys()) == 1
#     statment2 = statement2 and ('r' in gene_info.final_locations.keys() or 'g' in gene_info.final_locations.keys())
#     if statement1 or lysosomal_degradation_ptm_condition:
   

    secreted_protein_reactions = list()
    for secreted_protein in secreted_proteins:

        rxn = {modified_protein_g: -clathrin_coeff, secreted_protein: clathrin_coeff}
        # gtph hydrolysis
        rxn[metab.ntp_map_c['G']], rxn[metab.h2o_c], rxn[metab.ndp_map_c['G']], rxn[metab.pi_c], rxn[metab.h_c]  = -44, -44, 44, 44, 44

        secrete_protein = cobra.Reaction(gene_info.hgnc_id + '_Clathrin_IMPORTt' + secreted_protein.compartment)
        secrete_protein.subsystem = 'Protein Expression'
        if secreted_protein.compartment == 'e':
            rxn[biomass.protein_] = -clathrin_coeff*secreted_protein.mass
        secrete_protein.add_metabolites(rxn)
        secrete_protein.gene_reaction_rule = ' and '.join(mach.clathrin_m)
        secreted_protein_reactions += [secrete_protein]

    return secreted_protein_reactions, secreted_proteins


# # Secretory Pathway Protein Degradation

# In[219]:


def unfold_secretory_protein(gene_info, protein_metabolite):
    '''Remove PTMs and unfold proteins for lysosomal and secretory compartments. For lysosomal degradation, 
    only remove PTMs, there is no unfolding/misfolding as in ERAD.'''
    
    if not (protein_metabolite.compartment == 'r' or protein_metabolite.compartment == 'l'):
        raise ValueError('Protein metabolite does not have correct compartment')
    
    
    unfold_protein = cobra.Reaction(gene_info.hgnc_id + '_UNFOLD' + protein_metabolite.compartment)
    unfold_protein.subsystem = 'Protein Expression'
    unfolded_protein = protein_metabolite.copy()
    if protein_metabolite.compartment == 'r': # not "unfolding/misfolding" for lysosomal degradation
        unfolded_protein.id = unfolded_protein.id.replace('folded', 'unfolded')
    
    rxn = dict()
    rxn[protein_metabolite] = -1

    unfold_mach = list()
    elements = unfolded_protein.elements.copy()
#     lysosomal_degradation_ptm_condition = 'gpi' in gene_info.ptms.keys() and len(gene_info.ptms.keys()) == 1

    # PTM removals HERE #YOU ARE HERE 
    if 'ng' in gene_info.ptms.keys():
        raise ValueError('N-glycosylation not yet incorporated')
#     if lysosomal_degradation_ptm_condition:
#         raise ValueError('GPI-anchored proteins with no other ptms should be degraded via lysosomal pathway')
    if 'gpi' in gene_info.ptms.keys():
        unfolded_protein.id = unfolded_protein.id.replace('_GPI', '')
#         for e,c in metab.balanced_gpi.items():
#             elements[e] -= c
        
#         for m,s in metab.M4ATAer.copy().items(): # no lysosomal compartment metabolites
#             rxn[m] = -s
        rxn[metab.gpi_hs_r] = 1
        if protein_metabolite.compartment == 'r':
            rxn[metab.hdca_r], rxn[metab.h_r], rxn[metab.h2o_r] = -1,-1,1
        elif protein_metabolite.compartment == 'l':
            rxn[metab.hdca_l], rxn[metab.h_l], rxn[metab.h2o_l] = -1,-1,1
        

    if 'dsb' in gene_info.ptms.keys():
        unfolded_protein.id = unfolded_protein.id.replace('_DSB', '')
        number_DSB = gene_info.ptms['dsb']
        elements['H'] += 2*number_DSB
        # incorporate exchange with reductase in future versions
        if protein_metabolite.compartment == 'r':
            rxn[metab.o2_r], rxn[metab.h2o2_r] = number_DSB, -number_DSB
        elif protein_metabolite.compartment == 'l':
            rxn[metab.o2_l], rxn[metab.h2o2_l] = number_DSB, -number_DSB
        unfold_mach += mach.ERDJ5
    if 'og' in gene_info.ptms.keys():
        unfolded_protein.id = unfolded_protein.id.replace('_OG', '')
        number_Oglycans = gene_info.ptms['og']
        balance_og = {'C': (8 + 6 + 8)*number_Oglycans, # each 1/3 entry is for each 1/3 reactions in Jahir's model, in case want to separate in the future 
                      'H': (13 + 10 + 13)*number_Oglycans, 
                      'N': (1 + 0 + 1)*number_Oglycans, 
                      'O': (5 + 5 + 5)*number_Oglycans}
        for e,c in balance_og.items():
            elements[e] -= c
        
        if protein_metabolite.compartment == 'r':
            rxn[metab.udpacgal_r], rxn[metab.udpgal_r], rxn[metab.uacgam_r] = number_Oglycans, number_Oglycans, number_Oglycans
            rxn[metab.udp_r] = -3*number_Oglycans
            if metab.h_r in rxn.keys():
                rxn[metab.h_r] -= 3*number_Oglycans
            else:
                rxn[metab.h_r] = -3*number_Oglycans
        elif protein_metabolite.compartment == 'l':
            rxn[metab.udpacgal_l], rxn[metab.udpgal_g], rxn[metab.uacgam_g] = number_Oglycans, number_Oglycans, number_Oglycans
            rxn[metab.udp_l] = -3*number_Oglycans
            if metab.h_l in rxn.keys():
                rxn[metab.h_l] -= 3*number_Oglycans
            else:
                rxn[metab.h_l] = -3*number_Oglycans

    ###########
    
    unfolded_protein.elements = elements
    rxn[unfolded_protein] = 1
    
    biomass_change = (unfolded_protein.mass - protein_metabolite.mass)
    if biomass_change != 0:
        rxn[biomass.protein_] = biomass_change
        
    unfold_protein.add_metabolites(rxn)
    if len(unfold_mach) > 1:
        unfold_protein.gene_reaction_rule = ' and '.join(unfold_mach)
    elif len(unfold_mach) == 1:
        unfold_protein.gene_reaction_rule = unfold_mach[0]
    
    # PTMs - merged with unfolding reaction for now  <--same structure as adding the PTMs
    #     if 'ng' in gene_info.ptms.keys():
    #         raise ValueError('N-glycosylation not yet incorporated')
    #     if 'dsb' in gene_info.ptms.keys():
    #     else:
    #         unmodified_protein_r = unfolded_protein
    #     if 'gpi' in gene_info.ptms.keys():
    #     else:
    #         unmodified_protein_r = unmodified_protein_r
    #     if 'og' in gene_info.ptms.keys()
    #     else:
    #         unmodified_protein_r = unmodified_protein_r
    #     if len(gene_info.ptms) == 0:
    #         unmodified_protein_r = unfolded_protein

    unmodified_protein = unfolded_protein # for adding PTMs as separate reactions in future, if want to
    
    return unfold_protein, unmodified_protein

def build_erad_reactions(gene_info, retro_protein_r, unfolded_protein_c = None):
    if retro_protein_r.compartment != 'r':
        raise ValueError('ERAD can only occur with proteins in ER compartment')
    unfold_er_protein, unmodified_protein_r = unfold_secretory_protein(gene_info, retro_protein_r)
    

    # Retro-translocation
    retrotranslocate_protein = cobra.Reaction(gene_info.hgnc_id + '_RETROTRANSLOCATION')
    retrotranslocate_protein.subsystem = 'Protein Expression'
    if unfolded_protein_c == None: # those that underwent co-translational rather than post-translational
        unfolded_protein_c = unmodified_protein_r.change_compartment('c')
        # replace unfolded bc, if it underwent co-translational translocation
        # the signal peptide degradation step makes it such that it is not the same as any 
        # cytosolically translated unfolded proteins (multi-localization), i.e., 'i' destined proteins
        # and their degradation
        unfolded_protein_c.id = unfolded_protein_c.id.replace('unfolded', 'retrotranslocated_unfolded')
    
    rxn = {unmodified_protein_r: -1, unfolded_protein_c: 1}
    
    # FUTURE: separate nonglycosylated and glycosylated ERAD
    # if 'ng' in gene_info.ptms.keys() or 'og' in gene_info.ptms.keys():
    
     # glyco based ERAD from Jahir
    rxn = func.hydrolyze_atp(rxn, n_atp = 6, compartment = 'c') # from Jahir retro_TRANSLOC_2
    retrotranslocate_protein.add_metabolites(rxn)
    retrotranslocate_protein.gene_reaction_rule = ' and '.join(mach.retro_mach_glyco)
    

    erad_reactions = [unfold_er_protein, retrotranslocate_protein]

    return erad_reactions, unfolded_protein_c


# In[220]:


def build_endocytosis_reactions(gene_info, protein_pm, ub_args, protein_l = None):
    ##polyubiquitination for lysosomal targetting-------------------------------------
    polyubiquitinate_protein, polyub_protein_pm = protein_polyubiquitination(gene_info, protein_pm, compartment = 'pm', 
                                                                            ub_args = ub_args)
    
    polyub_protein_l = polyub_protein_pm.change_compartment('l')
    
    
    ##endocytosis--------------------------------------------------------------------------
    # combine dequbiquitination with endocytosis https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3987138/
    if protein_l == None:
        protein_l = protein_pm.change_compartment('l')
        

    endocytosis = cobra.Reaction(gene_info.hgnc_id + '_CLATHRIN_ENDOCYTOSIS')
    endocytosis.subsystem = 'Protein_Expression'
    
    
    rxn = {polyub_protein_pm: -1, protein_l: 1, metab.h2o_c: -1, ub_args['polyub_c']: 1}
    
    biomass_change = (protein_l.mass - polyub_protein_pm.mass)
    if biomass_change != 0:
        rxn[biomass.protein_] = biomass_change
    
    # gtp hydrolysis for vesicle scission
    rxn[metab.ntp_map_c['G']] = -round(gene_info.L_protein) * params.transport_translocation_atp_cost
    rxn[metab.h2o_c] -= round(gene_info.L_protein) * params.transport_translocation_atp_cost
    rxn[metab.ndp_map_c['G']] = round(gene_info.L_protein) * params.transport_translocation_atp_cost
    rxn[metab.pi_c] = round(gene_info.L_protein) * params.transport_translocation_atp_cost
    rxn[metab.h_c] = round(gene_info.L_protein) * params.transport_translocation_atp_cost

    
    endocytosis.add_metabolites(rxn)
    endocytosis.gene_reaction_rule = ' and '.join(mach.endocytic_machinery)
    
    return [polyubiquitinate_protein, endocytosis], protein_l


# In[221]:


def lysosomal_degradation(gene_info, protein_l):
    if protein_l.compartment != 'l':
        raise ValueError('This reaction only occurs in the lysosome')
    lysosomal_degradation_reactions = list()
    
    if len(gene_info.ptms.keys())>0:
        # not unfolding lysosomal protein is just removing the PTMs, not an actual unfolding reaction
        unfold_lysosomal_protein, unmodified_protein_l = unfold_secretory_protein(gene_info, protein_l)
        lysosomal_degradation_reactions += [unfold_lysosomal_protein]
    else:
        unmodified_protein_l = protein_l
    
    
    degrade_lysosomal_protein = cobra.Reaction(gene_info.hgnc_id + '_LYSOSOMAL_DEGRADATION')
    degrade_lysosomal_protein.subsystem = 'Protein_Expression'
    
    rxn = {metab.seq_amino_acid_map_l[aa_code]: aa_counts for aa_code, aa_counts in gene_info.amino_acid_counts.items()}
    rxn[unmodified_protein_l], rxn[metab.h2o_l] = -1, -(gene_info.L_protein-1)
    rxn[biomass.protein_] = -unmodified_protein_l.mass
    rxn = func.hydrolyze_atp(rxn, n_atp=gene_info.L_protein*params.proteolysis_translocation_atp_cost, compartment = 'l')
    
    degrade_lysosomal_protein.add_metabolites(rxn)
    degrade_lysosomal_protein.gene_reaction_rule = ' and '.join(mach.cathepsins)
    lysosomal_degradation_reactions += [degrade_lysosomal_protein]

    return lysosomal_degradation_reactions


# In[222]:


# # Jahir's NCBI GPRs to HGNC GPRs
# import re
# ehm = pd.read_csv(local_data_path + 'raw/identifiers.txt', sep = '\t')
# ehm = ehm.loc[ehm['NCBI gene ID'].dropna().index,:]
# ehm['NCBI gene ID'] = ehm['NCBI gene ID'].astype('int64').astype(str)


# test = ['(6400) and (84447) and (55666) and (7353) and (7415) and (79139) and (55829) and (91319) and (10134)']
# test = [re.findall(r'\d+', i) for i in test]
# test = sorted(set([item for sublist in test for item in sublist]))
# L_test = len(test)
# ehm = ehm[ehm['NCBI gene ID'].isin(test)]

# if len(ehm['NCBI gene ID'].unique()) != L_test:
#     print(set(test).difference(ehm['NCBI gene ID'].tolist()))
# if len(ehm['NCBI gene ID'].unique()) != ehm.shape[0]:
#     print('Redundant genes')
    
# mapper = dict(zip(ehm['NCBI gene ID'], ehm['HGNC ID']))


# # Protein Expression All

# In[223]:


def get_protein_expression_reactions(gene_info, mrna_transcript_c, mrna_deg_proxy, ub_args):
    # after transport, expand these to secretory pathways
    protein_expression_reactions, protein_metabolites = list(), list()
    
    # cytosolic transport: c, n, m, i, x and post-translational translocation
    if 'Cytosolic Tranport' in gene_info.final_locations.values() or gene_info.L_protein <= params.ptt_length: 
        translation_elongation_c, unfolded_protein_c = c_trln.translate_protein_cytosolic(gene_info, mrna_transcript_c, mrna_deg_proxy)
        protein_expression_reactions.append(translation_elongation_c)

        if 'Cytosolic Tranport' in gene_info.final_locations.values():
            if 'c' in gene_info.final_locations.keys() or 'x' in gene_info.final_locations.keys() or 'n' in gene_info.final_locations.keys():
                protein_folding_cytosolic, folded_protein_c = fold_protein_cytosolic(gene_info, unfolded_protein_c)
                protein_expression_reactions += [protein_folding_cytosolic]


                if 'c' in gene_info.final_locations.keys() or 'x' in gene_info.final_locations.keys() or ('n' in gene_info.final_locations.keys() and folded_protein_c.mass <= params.nuclear_diffusion_limit):
                   # cytoplasmic degradation of folded proteins: cytoplasmic proteins, peroxisomal proteins, or nuclear proteins undergoing passive diffusion
                    protein_expression_reactions += degrade_cytosolic_protein(gene_info, folded_protein_c, 
                                                                              ub_args = ub_args)

                    if 'c' in gene_info.final_locations.keys():
                        protein_metabolites += [folded_protein_c]

                    if 'x' in gene_info.final_locations.keys():
                        peroxisomal_reactions, folded_protein_x = get_peroxisomal_reactions(gene_info, folded_protein_c)
                        protein_expression_reactions += peroxisomal_reactions
                        protein_metabolites += [folded_protein_x]

                if 'n' in gene_info.final_locations.keys():
                    nuclear_reactions, folded_protein_n = get_nuclear_reactions(gene_info, folded_protein_c,
                                                                                ub_args = ub_args)
                    protein_expression_reactions += nuclear_reactions
                    protein_metabolites += [folded_protein_n]


            if 'i' in gene_info.final_locations.keys(): # no folding for i, but cytoplasmic degradation
                protein_expression_reactions += degrade_cytosolic_protein(gene_info, unfolded_protein_c, 
                                                                             ub_args = ub_args)
            # mitochondrial transport and degradation ('i' and 'm')
            if ('m' in gene_info.final_locations.keys()) or ('i' in gene_info.final_locations.keys()):
                if ('m' in gene_info.final_locations.keys()) and ('i' in gene_info.final_locations.keys()):
                    mitochondrial_reactions, mitochondrial_protein_metabolites = get_mitochondrial_reactions(gene_info, unfolded_protein_c, compartments = ['m','i'])
                elif 'm' in gene_info.final_locations.keys():
                    mitochondrial_reactions, mitochondrial_protein_metabolites = get_mitochondrial_reactions(gene_info, unfolded_protein_c, compartments = ['m'])
                elif 'i' in gene_info.final_locations.keys():
                    mitochondrial_reactions, mitochondrial_protein_metabolites = get_mitochondrial_reactions(gene_info, unfolded_protein_c, compartments = ['i'])
                protein_expression_reactions += mitochondrial_reactions
                protein_metabolites += mitochondrial_protein_metabolites
    
    # SECRETORY PATHWAY: r, g, l, e, pm proteins             
    if 'Canonical Secretion' in gene_info.final_locations.values():
        if gene_info.L_protein > params.ptt_length:
            ptt_ = False
        else:
            ptt_ = True
        # ptt_ variable is for ERAD reactions, on the off chance that the protein length was <= 160+22 
        # residues before signal peptide degradation occured since we update gene_info object during cotranslation
        # to be the new protein length after degradation of 22 residues

        if not ptt_: # co translational translocation
            ctt_reactions, folded_protein_r, gene_info = co_translational_translocation(gene_info, mrna_transcript_c, mrna_deg_proxy)
            protein_expression_reactions += ctt_reactions
        else: # post translational translocation
            ptt_reactions, folded_protein_r = post_translational_translocation(gene_info, unfolded_protein_c)
            protein_expression_reactions += ptt_reactions
            
        
        # er ptms
        if 'dsb' in gene_info.ptms.keys() or 'gpi' in gene_info.ptms.keys() or 'ng' in gene_info.ptms.keys():
            modification_er_reactions, modified_protein_r = modify_protein_er(gene_info, folded_protein_r)
            protein_expression_reactions += modification_er_reactions
        else:
            modified_protein_r = folded_protein_r
        
        # golgi and beyond transport; og ER resident proteins are retro-translocated;
        # ER/Golgi resident proteins with only a GPI anchor PTM will undergo lysosomal degradation rather than ERAD
        # lysosomal degradation only imported as part of its degradation pathway

#         lysosomal_degradation_ptm_condition = 'gpi' in gene_info.ptms.keys() and len(gene_info.ptms.keys()) == 1
        if len(set(['g', 'pm', 'e', 'l']).intersection(gene_info.final_locations.keys())) > 0 or 'og' in gene_info.ptms.keys():# or lysosomal_degradation_ptm_condition:
            golgi_import, protein_g = import_golgi(gene_info, modified_protein_r)
            protein_expression_reactions += [golgi_import]
            
            # golgi ptms
            if 'og' in gene_info.ptms.keys():
                modification_golgi_reactions, modified_protein_g = modify_protein_golgi(gene_info, protein_g)
                protein_expression_reactions += modification_golgi_reactions
            else: 
                modified_protein_g = protein_g
                
            # transport to plasma membrane, ECM, and lysosome 
            if len(set(['pm', 'e', 'l']).intersection(gene_info.final_locations.keys())) > 0:# or lysosomal_degradation_ptm_condition:
                secreted_protein_reactions, secreted_proteins = secrete_protein(gene_info, modified_protein_g)
                protein_expression_reactions += secreted_protein_reactions
                protein_metabolites += secreted_proteins
                            
            
            # retrograde transport
            if 'r' in gene_info.final_locations.keys() or 'g' in gene_info.final_locations.keys():# and not lysosomal_degradation_ptm_condition:
                # golgi retrograde transport for degradation 
                retrograde_transport, retro_protein_r = retrograde_er(gene_info, modified_protein_g)
                protein_expression_reactions += [retrograde_transport]
                if 'g' in gene_info.final_locations.keys():
                    protein_metabolites += [modified_protein_g]

        else:
            retro_protein_r = modified_protein_r # for ER resident proteins with no O-glycosylation, they are not transported to Golgi and retrograde transported
        
            
        # ERAD: ER and Golgi-resident proteins 
        if ('r' in gene_info.final_locations.keys() or 'g' in gene_info.final_locations.keys()):# and not lysosomal_degradation_ptm_condition:
            if 'r' in gene_info.final_locations.keys():
                protein_metabolites += [retro_protein_r]
            if ptt_:
                erad_reactions, unfolded_protein_c = build_erad_reactions(gene_info, retro_protein_r, 
                                                                          unfolded_protein_c = unfolded_protein_c)
                protein_expression_reactions += erad_reactions
                if 'i' not in gene_info.final_locations.keys(): # this reaction doesn't already exist
                    protein_expression_reactions += degrade_cytosolic_protein(gene_info, unfolded_protein_c, ub_args)
            else:
                erad_reactions, unfolded_protein_c = build_erad_reactions(gene_info, retro_protein_r, 
                                                                          unfolded_protein_c = None)
                protein_expression_reactions += erad_reactions
                # since metabolite id is different for unfolded_protein_c (see erad) and proteasomal degradation
                # reactions use metabolite id rather than gene_info.hgnc_id, 
                # don't need to worry about overlap with 'i' compartment degradation reactions
                # in the case of multi-localization 
                #(this unfolded protein is different than cytosolically translated ones bc of the)
                # signal peptide degradation reaction
                protein_expression_reactions += degrade_cytosolic_protein(gene_info, unfolded_protein_c, ub_args)
        # PM/L degradation needed
        # endocytosis of plasma membrane proteins
        if 'l' in gene_info.final_locations.keys() or 'pm' in gene_info.final_locations.keys():
            if 'l' in gene_info.final_locations.keys():
                protein_l = [p for p in secreted_proteins if p.compartment == 'l'][0]
            else: 
                protein_l = None
            
            # endocytosis
            if 'pm' in gene_info.final_locations.keys():
                protein_pm = [p for p in secreted_proteins if p.compartment == 'pm'][0]
                endocytosis_reactions, protein_l = build_endocytosis_reactions(gene_info, protein_pm, 
                                                                               protein_l = protein_l,
                                                                               ub_args = ub_args)                
                protein_expression_reactions += endocytosis_reactions
            
            # lysosomal degradation
            protein_expression_reactions += lysosomal_degradation(gene_info, protein_l)

    elif 'Non-Canonical Secretion' in gene_info.final_locations.values():
        raise ValueError('Model does not currently account for non-canonical secretion')
    
        
    return protein_expression_reactions, protein_metabolites


# In[224]:


# import random
# import pandas as pd
# from expression.gene_information import gene_information
# import expression.build_mrna_expression_reactions as build_mrna
# from expression.protein_expression import ubiquitin

# psim_toy = pd.DataFrame(columns = ['HGNC_ID', 'PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ', 'POLYA_LENGTH', 'TMD', 
#                                'SP', 'N_INTRONS', 'DSB', 'GPI', 'OG', 'LOCATION'])

# hgnc_id, premrna_seq = 'HGNC:TOY', ''.join(random.choices(['U', 'C', 'G', 'A'], k = 100))
# mrna_seq = premrna_seq[25:75]
# # note that there is no check that the protein_sequence corresponds to the mrna_sequence beyond checking for the length
# protein_seq = ''.join(random.choices(params.amino_acids, k = int(len(mrna_seq)/3)))
# polyA_length, tmd, sp, n_introns, dsb, gpi, og  = None, 1, True, 0, 2, 2, 2
# ub_args = ubiquitin.express_ubiquitin(compress_mrna = False)

# import itertools
# reactions = list()
# for l in list(itertools.combinations(params.compartments.keys(),2)):
#     location = list(l)
#     psim_toy.loc[0,:] = [hgnc_id, premrna_seq, mrna_seq, protein_seq, polyA_length, tmd, sp, n_introns, dsb, gpi, og, location]
#     gene_info = gene_information(hgnc_id, premrna_seq, mrna_seq, protein_seq,
#                      ptms = {}, tmd = tmd, sp = sp, polyA_length = polyA_length, 
#                      n_introns = n_introns) 
#     gene_info.get_final_locations(metabolic_model = cobra.Model(''), final_locations = location)

#     transcription_reactions, mrna_transcript_c, mrna_deg_proxy = build_mrna.get_mrna_expression_reactions(gene_info)
#     protein_expression_reactions, protein_metabolites = get_protein_expression_reactions(gene_info, 
#                                                      mrna_transcript_c, mrna_deg_proxy, 
#                                                     ub_args = ub_args)
#     reactions += protein_expression_reactions


