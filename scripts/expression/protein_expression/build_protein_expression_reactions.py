#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.insert(1, '../../../scripts/') # comment out in python script
from utils import machinery as mach
from utils import parameters as params
from utils import metabolites as metab
from utils import functions as func
from core.reaction import Protein_Expression_Reaction

from macromolecules.protein import Protein
from expression.protein_expression import cytosolic_translation as c_trln
from expression.protein_expression import degradation

from uniform_processes.build_trna_expression_reactions import modified_trna_transcript_c, charged_trna_map


# In[2]:


def fold_protein_cytosolic(gene_info, unfolded_protein_c):
    # extending proteostasis network in the future would be good
    # will need to make sure inputs to each compartment-specific reactions are at the correct folding stage
    # e.g., mitochondria currently takes unfolded protein, and in future we may want it to take a partially folded
    
    folded_protein_c = unfolded_protein_c.copy()
    folded_protein_c.id = folded_protein_c.id.replace('unfolded', 'folded')
    rxn = {unfolded_protein_c: -1, folded_protein_c: 1}
    protein_folding = Protein_Expression_Reaction(gene_info.hgnc_id + '_CYTOSOLIC_PROTEIN_FOLDING', 
                                         hgnc_id = gene_info.hgnc_id)
    
    if gene_info.L_protein > 100: #chaperone assisted for larger proteins - https://www.nature.com/articles/nature10317
        rxn = func.hydrolyze_atp(rxn, n_atp = gene_info.L_protein*params.proteolysis_translocation_atp_cost, compartment = 'c')
        protein_folding.gene_reaction_rule = ' and '.join(mach.HSP40_c + mach.HSP70_c) # GPRs
    
    protein_folding.add_metabolites(rxn)

    return protein_folding, folded_protein_c


# # Nuclear Reactions

# In[3]:


def transport_nuclear_protein(gene_info, folded_protein_c):

    folded_protein_n = folded_protein_c.change_compartment('n')

    nuclear_import = Protein_Expression_Reaction(gene_info.hgnc_id + '_IMPORTtn', hgnc_id = gene_info.hgnc_id)
    nuclear_import.subsytem = 'Protein_Expression'
    
    import_rxn = {folded_protein_c: -1, folded_protein_n: 1}

    if folded_protein_c.formula_weight/1000 > params.nuclear_diffusion_limit:
        # gtp hydrolysis per import
        import_rxn[metab.gtp_n], import_rxn[metab.h2o_n], import_rxn[metab.gdp_n], import_rxn[metab.pi_n], import_rxn[metab.h_n]  = -1, -1, 1, 1, 1
        nuclear_import.add_metabolites(import_rxn)
        nuclear_import.gene_reaction_rule = ' and '.join(mach.importins + mach.RAN)

    else: # diffusion
        nuclear_import.add_metabolites(import_rxn)
        nuclear_import.lower_bound = -1000 # reversible
        
    return nuclear_import, folded_protein_n

def get_nuclear_reactions(gene_info, folded_protein_c, ub_args):
    nuclear_import, folded_protein_n = transport_nuclear_protein(gene_info, folded_protein_c)
    nuclear_degradation_reactions = degradation.degrade(macromolecule = folded_protein_n, **{'ub_args': ub_args})
    
    return [nuclear_import] + nuclear_degradation_reactions, folded_protein_n


# # Mitochondrial Reactions
# 

# In[4]:


# i is intermembrane space, but called inner in compartments BIGG
# stick to notation and use inner instead of inter in reaction naming

def transport_mitochondrial_matrix(gene_info, unfolded_protein_c):
    # transport and folding
    if unfolded_protein_c.compartment != 'c':
        raise ValueError('Only cytoplasmic proteins can be transported to mitochondrial matrix')
    
    mitochondrial_matrix_transport = Protein_Expression_Reaction(gene_info.hgnc_id + '_IMPORTtm', hgnc_id = gene_info.hgnc_id)
    mitochondrial_matrix_transport.subsytem = 'Protein_Expression'
    pre_protein_m = unfolded_protein_c.change_compartment('m')
    pre_protein_m.id = pre_protein_m.id.replace('unfolded', 'folded_pre')
    
    rxn = {unfolded_protein_c: -1, pre_protein_m: 1}
    # ATP hydrolysis for transport, assums 1 ATP consumed per 2 residues
    rxn = func.hydrolyze_atp(rxn, n_atp = gene_info.L_protein*params.transport_translocation_atp_cost, compartment = 'm')
    
    
    mitochondrial_matrix_transport.add_metabolites(rxn)
    mitochondrial_matrix_transport.gene_reaction_rule = ' and '.join(mach.TOM + mach.TIM23_PAM + mach.HSP70_m)
    
    return mitochondrial_matrix_transport, pre_protein_m

def mitochondrial_matrix_protein_processing(gene_info, pre_protein_m):
    # implement this in the future: cleavage of MTS (and degradation of MTS)
    
    processed_protein_m, aa_counts_processed_m, L_processed_protein_m = pre_protein_m, gene_info.amino_acid_counts.copy(), gene_info.L_protein
    process_mitochondrial_matrix_protein = None
    return process_mitochondrial_matrix_protein, processed_protein_m, aa_counts_processed_m, L_processed_protein_m


def transport_mitochondrial_inter(gene_info, processed_protein_m):
    # upper left Fig 12-29 https://www.ncbi.nlm.nih.gov/books/NBK26828/ 
    # import to matrix then re-export to inter membrane space
    
    if processed_protein_m.compartment != 'm':
        raise ValueError('Only the mechanism of mitochondrial matrix import and re-export to inter membrane is considered')
    
    mitochondrial_inter_transport = Protein_Expression_Reaction(gene_info.hgnc_id + '_IMPORTti', hgnc_id = gene_info.hgnc_id)
    mitochondrial_inter_transport.subsytem = 'Protein_Expression'
    pre_protein_i = processed_protein_m.change_compartment('i')
    
    rxn = {processed_protein_m: -1, pre_protein_i: 1}
    
    mitochondrial_inter_transport.add_metabolites(rxn)
    mitochondrial_inter_transportgene_reaction_rule = mach.OXA[0]
    
    return mitochondrial_inter_transport, pre_protein_i

def mitochondrial_inter_protein_processing(gene_info, pre_protein_i):
    # implement this in the future: cleavage of secondary sequence (and degradation)
    
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
        if aa_counts_processed_m != processed_protein_m._amino_acid_counts:
            raise ValueError('Internal: Unaccounted for change in amino acid counts for: ' + processed_protein_m.id)
        mitochondrial_reactions += degradation.degrade(macromolecule = processed_protein_m)
        for r in mitochondrial_reactions:
            r._final_compartments.append('m')
        mitochondrial_protein_metabolites += [processed_protein_m]
    if 'i' in compartments:
        mitochondrial_inter_transport, pre_protein_i = transport_mitochondrial_inter(gene_info, processed_protein_m)
        process_mitochondrial_inter_protein, processed_protein_i, aa_counts_processed_i, L_processed_protein_i = mitochondrial_inter_protein_processing(gene_info, pre_protein_i)
        if process_mitochondrial_matrix_protein != None:
            mitochondrial_reactions += [process_mitochondrial_inter_protein]  
        if aa_counts_processed_i != processed_protein_i._amino_acid_counts:
            raise ValueError('Internal: Unaccounted for change in amino acid counts for: ' + processed_protein_i.id)
        
        mitochondrial_reactions += [mitochondrial_inter_transport] + degradation.degrade(macromolecule = processed_protein_i)
        for r in mitochondrial_reactions:
            r._final_compartments.append('i')
        mitochondrial_protein_metabolites += [processed_protein_i]

    return mitochondrial_reactions, mitochondrial_protein_metabolites


# # Peroxisomal

# In[5]:


def transport_peroxisome(gene_info, folded_protein_c):
    if folded_protein_c.compartment != 'c':
        raise ValueError('Only cytoplasmic proteins can be transported to mitochondrial matrix')
    
    peroxisomal_transport = Protein_Expression_Reaction(gene_info.hgnc_id + '_IMPORTtx', hgnc_id = gene_info.hgnc_id)
    peroxisomal_transport.subsytem = 'Protein_Expression'
    folded_protein_x = folded_protein_c.change_compartment('x')
    
    rxn = {folded_protein_c: -1, folded_protein_x: 1}
    # ATP hydrolysis for transport--translocation of protein, export of PEX5S receptor
    rxn = func.hydrolyze_atp(rxn, n_atp = (gene_info.L_protein+mach.L_PEX5)*params.transport_translocation_atp_cost, 
                        compartment = 'x')
    
    
    peroxisomal_transport.add_metabolites(rxn)
    peroxisomal_transport.gene_reaction_rule = ' and '.join(mach.peroxins + mach.AWP1)
    
    return peroxisomal_transport, folded_protein_x

def get_peroxisomal_reactions(gene_info, folded_protein_c):
    pr = list()
    peroxisomal_transport, folded_protein_x = transport_peroxisome(gene_info, folded_protein_c)
    pr = [peroxisomal_transport] + degradation.degrade(macromolecule = folded_protein_x)
    
    return pr, folded_protein_x


# # Secretory Pathway
# 
# adapted from Jahir's Recon2_2s

# # ER transport

# In[6]:


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
    
    if gene_info.tmd > 0 or 'pm' in gene_info.all_locations.keys(): # membrane secreted protein
        post_translational_translocation_r = Protein_Expression_Reaction(gene_info.hgnc_id + '_post_TRANSLOC_3A_IMPORTtr', 
                                                                 hgnc_id = gene_info.hgnc_id)
        post_translational_translocation_r.subsytem = 'Protein_Expression'
        post_translational_translocation_r.gene_reaction_rule = ' and '.join(mach.ASNA1 + mach.WRB)
        
        # complex cleavage from jahir's (+h2o_c, +h_c, + pi_c) not included 
        post_translational_translocation_r.add_metabolites(rxn)
        ptt_reactions += [post_translational_translocation_r]
     
    else: #non membrane secreted protein
        number_BiP = gene_info.L_protein/40
        post_translational_translocation_r = Protein_Expression_Reaction(gene_info.hgnc_id + '_post_TRANSLOC_3B_IMPORTtr', 
                                                                 hgnc_id = gene_info.hgnc_id)
        post_translational_translocation_r.subsytem = 'Protein_Expression'
        post_translational_translocation_r.gene_reaction_rule = ' and '.join(mach.ptnm)
        
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
    
    rxn = dict()
    for aa_code, aa_count in gene_info.amino_acid_counts.items():
        rxn[charged_trna_map[aa_code]] = -aa_count # tRNA consumption
    
    rxn[modified_trna_transcript_c]  = gene_info.L_protein 
    
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
    
    rxn[unprocessed_protein_r] = 1
    rxn = func.hydrolyze_atp(rxn, n_atp = number_BiP, compartment = 'r')
    
    #------------------------------------------------------------------------------------

    co_translational_translocation_r = Protein_Expression_Reaction(gene_info.hgnc_id + '_co_TRANSLOC_IMPORTtr', 
                                                           translation = True, hgnc_id = gene_info.hgnc_id)
    co_translational_translocation_r.gene_reaction_rule = ' and '.join(mach.ctnm + mach.translation_efs + ['ribosome'])
    
    co_translational_translocation_r.add_metabolites(rxn)
    #coupling
    mrna_deg_proxy.couple(value = -gene_info.coupling['mrna_degradation'])
    mrna_transcript_c.couple(type = 'mrna_formation', value = -gene_info.coupling['mrna_formation'])
    co_translational_translocation_r.couple(metabolites = [mrna_deg_proxy, mrna_transcript_c], 
                                 types = ['mrna_degradation', 'mrna_formation'])
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
    
    sp_degradation = Protein_Expression_Reaction(gene_info.hgnc_id + '_SP_degradationr', hgnc_id = gene_info.hgnc_id)
    sp_degradation.add_metabolites(rxn)
    sp_degradation.gene_reaction_rule = mach.sp_rule
    ctt_reactions += [sp_degradation]

    return ctt_reactions, folded_protein_r, gene_info


# # ER Modifications

# In[7]:


def form_disulfide_bond(gene_info, folded_protein_r):
    number_DSB = gene_info.ptms['dsb']
    disulfide_bond_formation = Protein_Expression_Reaction(gene_info.hgnc_id + '_DSBr', hgnc_id = gene_info.hgnc_id)
    modified_protein_dsb_r = folded_protein_r.copy()
    modified_protein_dsb_r.id = modified_protein_dsb_r.id.replace('folded', 'folded_DSB')
    elements = folded_protein_r.elements.copy()
    elements['H'] -= 2*number_DSB
    modified_protein_dsb_r.elements = elements
    # diagram https://www.google.com/url?sa=i&url=https%3A%2F%2Fen.wikipedia.org%2Fwiki%2FProtein_disulfide-isomerase&psig=AOvVaw0bGpff4XX1eYEF61H1RJKw&ust=1597273135069000&source=images&cd=vfe&ved=0CAIQjRxqFwoTCJi6l6GglOsCFQAAAAAdAAAAABAJ
    # incorporate exchange with PDI in future versions
    rxn = {folded_protein_r: -1, modified_protein_dsb_r: 1, metab.o2_r: -number_DSB, metab.h2o2_r: number_DSB}
    disulfide_bond_formation.add_metabolites(rxn)
    disulfide_bond_formation.gene_reaction_rule = mach.P4HB[0]
    
    return disulfide_bond_formation, modified_protein_dsb_r

def form_gpi(gene_info, modified_protein_r):
    gpi_formation = Protein_Expression_Reaction(gene_info.hgnc_id + '_GPIr', hgnc_id = gene_info.hgnc_id)
    modified_protein_gpi_r = modified_protein_r.copy() 
    modified_protein_gpi_r.id = modified_protein_gpi_r.id.replace('folded', 'folded_GPI')

    # need to figure this out correctly!!
    elements = modified_protein_r.elements.copy()
    for e,c in metab.balanced_gpi.items():
        if e in elements.keys():
            elements[e] += c
        else:
            elements[e] = c
    modified_protein_gpi_r.elements = elements

    rxn = dict()
#     rxn[metab.hdca_r], rxn[metab.gpi_hs_r], rxn[metab.h_r], rxn[metab.h2o_r], rxn[metab.gpi_sig_r] = 1,-1,1,-1,1
    rxn[metab.hdca_r], rxn[metab.gpi_hs_r], rxn[metab.h_r], rxn[metab.h2o_r] = 1,-1,1,-1
    rxn[modified_protein_r], rxn[modified_protein_gpi_r] = -1, 1

    gpi_formation.add_metabolites(rxn)
    gpi_formation.gene_reaction_rule = ' and '.join(mach.gpi_machinery)
        
    return gpi_formation, modified_protein_gpi_r

def glycosylate_n_linked(gene_info, modified_protein_r):
    raise ValueError('N-glycosylation not yet incorporated')
#     n_glycosylation = Protein_Expression_Reaction(gene_info.hgnc_id + 'NGLYCOr')
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

# In[8]:


def import_golgi(gene_info, modified_protein_r):
    V = modified_protein_r.formula_weight/1000 * 1.21 / 1000.0 # Protein Volume in nm^3
    copii_coeff = int(round(268082.35 * params.Kv / V))
    
    protein_g = modified_protein_r.change_compartment('g')
    
    rxn = {modified_protein_r: -copii_coeff, protein_g: copii_coeff}
    # gtp hydrolysis
    rxn[metab.ntp_map_c['G']], rxn[metab.h2o_c], rxn[metab.ndp_map_c['G']], rxn[metab.pi_c], rxn[metab.h_c]  = -94, -94, 94, 94, 94

    golgi_import = Protein_Expression_Reaction(gene_info.hgnc_id + '_COPII_IMPORTtg', hgnc_id = gene_info.hgnc_id)
    golgi_import.add_metabolites(rxn)
    
    
    if 'gpi' in gene_info.ptms.keys() and 'ng' not in gene_info.ptms.keys(): # this if statement is analogous to Recon2.2S's connector statements in copii reactions
        golgi_import.gene_reaction_rule = ' and '.join(mach.copii_gpi_m)
    else:
        golgi_import.gene_reaction_rule = ' and '.join(mach.copii_r_m)
        
    
    return golgi_import, protein_g

def glycosylate_o_linked(gene_info, protein_g):
    number_Oglycans = gene_info.ptms['og']
    o_glycosylation = Protein_Expression_Reaction(gene_info.hgnc_id + '_OGg', hgnc_id = gene_info.hgnc_id)

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
           metab.udp_g: 3* number_Oglycans}
    
    
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


# # Lysosomal, Extracellular, and Plasma Membrane Transport

# In[9]:


def secrete_protein(gene_info, modified_protein_g):
    V = modified_protein_g.formula_weight/1000 * 1.21 / 1000.0 # Protein Volume in nm^3
    clathrin_coeff = int(round(29880.01 * params.Kv / V)) # Number of proteins per clathrin vesicle  

    secreted_proteins = list()
    for nc in ['e', 'pm', 'l']:
        if nc in gene_info.all_locations.keys():
            secreted_protein = modified_protein_g.change_compartment(nc)
            secreted_proteins += [secreted_protein]
    
#     statement1 = 'pm' in gene_info.all_locations.keys()
#     lysosomal_degradation_ptm_condition = 'gpi' in gene_info.ptms.keys() and len(gene_info.ptms.keys()) == 1
#     statment2 = statement2 and ('r' in gene_info.all_locations.keys() or 'g' in gene_info.all_locations.keys())
#     if statement1 or lysosomal_degradation_ptm_condition:
   

    secreted_protein_reactions = list()
    for secreted_protein in secreted_proteins:

        rxn = {modified_protein_g: -clathrin_coeff, secreted_protein: clathrin_coeff}
        # gtph hydrolysis
        rxn[metab.ntp_map_c['G']], rxn[metab.h2o_c], rxn[metab.ndp_map_c['G']], rxn[metab.pi_c], rxn[metab.h_c]  = -44, -44, 44, 44, 44

        secrete_protein = Protein_Expression_Reaction(gene_info.hgnc_id + '_Clathrin_IMPORTt' + secreted_protein.compartment, 
                                             hgnc_id = gene_info.hgnc_id)
        secrete_protein.add_metabolites(rxn)
        secrete_protein.gene_reaction_rule = ' and '.join(mach.clathrin_m)
        secreted_protein_reactions += [secrete_protein]

    return secreted_protein_reactions, secreted_proteins


# # Secretory Pathway Protein Degradation

# In[10]:


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

# In[11]:


def get_protein_expression_reactions(gene_info, mrna_transcript_c, mrna_deg_proxy, ub_args):
    # after transport, expand these to secretory pathways
    protein_expression_reactions, protein_metabolites = list(), list()
    
    # cytosolic transport: c, n, m, i, x and post-translational translocation
    if 'Cytosolic Tranport' in gene_info.all_locations.values() or gene_info.L_protein <= params.ptt_length: 
        translation_elongation_c, unfolded_protein_c = c_trln.translate_protein_cytosolic(gene_info, mrna_transcript_c, mrna_deg_proxy)
        translation_elongation_c._final_compartments += [comp for comp, v in gene_info.all_locations.items() if v == 'Cytosolic Tranport']
        protein_expression_reactions.append(translation_elongation_c)

        if 'Cytosolic Tranport' in gene_info.all_locations.values():
            if 'c' in gene_info.all_locations.keys() or 'x' in gene_info.all_locations.keys() or 'n' in gene_info.all_locations.keys():
                protein_folding_cytosolic, folded_protein_c = fold_protein_cytosolic(gene_info, unfolded_protein_c)
                protein_folding_cytosolic._final_compartments += list(set(gene_info.all_locations).intersection(['c', 'x', 'n']))
                protein_expression_reactions += [protein_folding_cytosolic]


                if 'c' in gene_info.all_locations or 'x' in gene_info.all_locations or ('n' in gene_info.all_locations and folded_protein_c.formula_weight/1000 <= params.nuclear_diffusion_limit):
                   # cytoplasmic degradation of folded proteins: cytoplasmic proteins, peroxisomal proteins, or nuclear proteins undergoing passive diffusion
                    dr = degradation.degrade(folded_protein_c, **{'ub_args':ub_args})
                    fc = list(set(gene_info.all_locations).intersection(['c', 'x', 'n']))
                    for r in dr:
                        r._final_compartments += fc
                    protein_expression_reactions += dr

                    if 'c' in gene_info.all_locations.keys():
                        protein_metabolites += [folded_protein_c]

                    if 'x' in gene_info.all_locations.keys():
                        peroxisomal_reactions, folded_protein_x = get_peroxisomal_reactions(gene_info, folded_protein_c)
                        for r in peroxisomal_reactions:
                            r._final_compartments.append('x')
                        protein_expression_reactions += peroxisomal_reactions
                        protein_metabolites += [folded_protein_x]

                if 'n' in gene_info.all_locations.keys():
                    nuclear_reactions, folded_protein_n = get_nuclear_reactions(gene_info, folded_protein_c,
                                                                                ub_args = ub_args)
                    for r in nuclear_reactions:
                        r._final_compartments.append('n')
                    protein_expression_reactions += nuclear_reactions
                    protein_metabolites += [folded_protein_n]


            if 'i' in gene_info.all_locations.keys(): # no folding for i, but cytoplasmic degradation
                dr = degradation.degrade(macromolecule = unfolded_protein_c, **{'ub_args': ub_args})
                for r in dr:
                    r._final_compartments.append('i')
                protein_expression_reactions += dr
            # mitochondrial transport and degradation ('i' and 'm')
            if ('m' in gene_info.all_locations.keys()) or ('i' in gene_info.all_locations.keys()):
                if ('m' in gene_info.all_locations.keys()) and ('i' in gene_info.all_locations.keys()):
                    mitochondrial_reactions, mitochondrial_protein_metabolites = get_mitochondrial_reactions(gene_info, unfolded_protein_c, compartments = ['m','i'])
                elif 'm' in gene_info.all_locations.keys():
                    mitochondrial_reactions, mitochondrial_protein_metabolites = get_mitochondrial_reactions(gene_info, unfolded_protein_c, compartments = ['m'])
                elif 'i' in gene_info.all_locations.keys():
                    mitochondrial_reactions, mitochondrial_protein_metabolites = get_mitochondrial_reactions(gene_info, unfolded_protein_c, compartments = ['i'])
                protein_expression_reactions += mitochondrial_reactions
                protein_metabolites += mitochondrial_protein_metabolites
    
    # SECRETORY PATHWAY: r, g, l, e, pm proteins             
    if 'Canonical Secretion' in gene_info.all_locations.values():
        if gene_info.L_protein > params.ptt_length:
            ptt_ = False
        else:
            ptt_ = True
        # ptt_ variable is for ERAD reactions, on the off chance that the protein length was <= 160+22 
        # residues before signal peptide degradation occured since we update gene_info object during cotranslation
        # to be the new protein length after degradation of 22 residues
        fc = [comp for comp, v in gene_info.all_locations.items() if v == 'Canonical Secretion']
        if not ptt_: # co translational translocation
            ctt_reactions, folded_protein_r, gene_info = co_translational_translocation(gene_info, mrna_transcript_c, mrna_deg_proxy)
            for r in ctt_reactions:
                r._final_compartments += fc
            protein_expression_reactions += ctt_reactions
        else: # post translational translocation
            ptt_reactions, folded_protein_r = post_translational_translocation(gene_info, unfolded_protein_c)
            for r in ptt_reactions:
                r._final_compartments += fc
            protein_expression_reactions += ptt_reactions
            
        
        # er ptms
        if 'dsb' in gene_info.ptms.keys() or 'gpi' in gene_info.ptms.keys() or 'ng' in gene_info.ptms.keys():
            modification_er_reactions, modified_protein_r = modify_protein_er(gene_info, folded_protein_r)
            for r in modification_er_reactions:
                r._final_compartments += fc
            protein_expression_reactions += modification_er_reactions
        else:
            modified_protein_r = folded_protein_r
        
        # golgi and beyond transport; og ER resident proteins are retro-translocated;
        # ER/Golgi resident proteins with only a GPI anchor PTM will undergo lysosomal degradation rather than ERAD
        # lysosomal degradation only imported as part of its degradation pathway

#         lysosomal_degradation_ptm_condition = 'gpi' in gene_info.ptms.keys() and len(gene_info.ptms.keys()) == 1
        if len(set(['g', 'pm', 'e', 'l']).intersection(gene_info.all_locations.keys())) > 0 or 'og' in gene_info.ptms.keys():# or lysosomal_degradation_ptm_condition:
            golgi_import, protein_g = import_golgi(gene_info, modified_protein_r)
            golgi_import._final_compartments += fc
            protein_expression_reactions += [golgi_import]
            
            # golgi ptms
            if 'og' in gene_info.ptms.keys():
                modification_golgi_reactions, modified_protein_g = modify_protein_golgi(gene_info, protein_g)
                for r in modification_golgi_reactions:
                    r._final_compartments += fc
                protein_expression_reactions += modification_golgi_reactions
            else: 
                modified_protein_g = protein_g
                
            # transport to plasma membrane, ECM, and lysosome 
            if len(set(['pm', 'e', 'l']).intersection(gene_info.all_locations.keys())) > 0:# or lysosomal_degradation_ptm_condition:
                secreted_protein_reactions, secreted_proteins = secrete_protein(gene_info, modified_protein_g)
                fc = list(set(gene_info.all_locations).intersection(['pm', 'e', 'l']))
                for r in secreted_protein_reactions:
                    r._final_compartments += fc
                protein_expression_reactions += secreted_protein_reactions
                protein_metabolites += secreted_proteins
                            
            
            # retrograde transport
            if 'r' in gene_info.all_locations or 'g' in gene_info.all_locations:# and not lysosomal_degradation_ptm_condition:
                # golgi retrograde transport for degradation 
                retrograde_transport, retro_protein_r = degradation.retrograde_er(modified_protein_g)
                retrograde_transport._final_compartments += list(set(gene_info.all_locations).intersection(['r', 'g']))
                protein_expression_reactions += [retrograde_transport]
                if 'g' in gene_info.all_locations.keys():
                    protein_metabolites += [modified_protein_g]

        else:
            retro_protein_r = modified_protein_r # for ER resident proteins with no O-glycosylation, they are not transported to Golgi and retrograde transported
        
            
        # ERAD: ER and Golgi-resident proteins 
        if ('r' in gene_info.all_locations or 'g' in gene_info.all_locations):# and not lysosomal_degradation_ptm_condition: 
            if 'r' in gene_info.all_locations.keys():
                protein_metabolites += [retro_protein_r]
            rpdr = list()
            fc = list(set(gene_info.all_locations).intersection(['r', 'g']))
            if ptt_:
                erad_reactions, unfolded_protein_c = degradation.degrade(macromolecule = retro_protein_r,
                                                     **{'unfolded_protein_c': unfolded_protein_c})
                for r in erad_reactions:
                    r._final_compartments += fc
                rpdr += erad_reactions
                if 'i' not in gene_info.all_locations.keys(): # this reaction doesn't already exist
                    dr = degradation.degrade(unfolded_protein_c, **{'ub_args': ub_args})
                    for r in dr:
                        r._final_compartments += fc
                    rpdr += dr
            else:
                erad_reactions, unfolded_protein_c = degradation.degrade(macromolecule = retro_protein_r,
                                                                          **{'unfolded_protein_c': None})
                for r in erad_reactions:
                    r._final_compartments += fc
                rpdr += erad_reactions
                # since metabolite id is different for unfolded_protein_c (see erad) and proteasomal degradation
                # reactions use metabolite id rather than gene_info.hgnc_id, 
                # don't need to worry about overlap with 'i' compartment degradation reactions
                # in the case of multi-localization 
                #(this unfolded protein is different than cytosolically translated ones bc of the)
                # signal peptide degradation reaction
                dr = degradation.degrade(macromolecule = unfolded_protein_c, **{'ub_args': ub_args})
                for r in dr:
                    r._final_compartments += fc
                rpdr += dr
                for r in rpdr:
                    if 'g' in gene_info.all_locations.keys():
                        r._update_tracking(modified_protein_g) # not explicitly accounted for for protein monomers
                    r._update_tracking([retro_protein_r, unfolded_protein_c])
                    r._consolidate_macromolecules()    
            
            protein_expression_reactions += rpdr
            
            
            
        # PM/L degradation needed
        # endocytosis of plasma membrane proteins
        if 'l' in gene_info.all_locations.keys() or 'pm' in gene_info.all_locations.keys():
            if 'l' in gene_info.all_locations.keys():
                protein_l = [p for p in secreted_proteins if p.compartment == 'l'][0]
            else: 
                protein_l = None
            
            # endocytosis
            ppdr = list()
            if 'pm' in gene_info.all_locations.keys():
                protein_pm = [p for p in secreted_proteins if p.compartment == 'pm'][0]
                endocytosis_reactions, protein_l = degradation.build_endocytosis_reactions(                                                   macromolecule_pm = protein_pm, 
                                                   **{'ub_args': ub_args, 'macromolecule_l': protein_l})
                for r in endocytosis_reactions:
                    r._final_compartments.append('pm')
                        
                ppdr += endocytosis_reactions
            
            # lysosomal degradation
            fc = list(set(gene_info.all_locations).intersection(['l', 'pm']))
            dr = degradation.degrade(macromolecule = protein_l)
            for r in dr:
                r._final_compartments += fc
            ppdr += dr
            
            if 'pm' in gene_info.all_locations.keys():
                for r in ppdr:
                    r._update_tracking(protein_pm) # not explicitly accounted for for protein monomers
                    r._consolidate_macromolecules()
            
            protein_expression_reactions += ppdr

    elif 'Non-Canonical Secretion' in gene_info.all_locations.values():
        raise ValueError('Model does not currently account for non-canonical secretion')
    
    for r in protein_expression_reactions:
        if r.subsystem != 'Protein_Degradation':
            r.subsystem = 'Protein_Expression'
        r._final_compartments = list(set(r._final_compartments))

    for m in protein_metabolites:
        if m.compartment not in gene_info.machinery_locations.keys():
            m.non_machinery = True
    
    return protein_expression_reactions, protein_metabolites


# In[ ]:


# import random
# import cobra
# import pandas as pd
# from expression.gene_information import gene_information
# import expression.build_mrna_expression_reactions as build_mrna
# from expression.protein_expression import ubiquitin

# psim_toy = pd.DataFrame(columns = ['HGNC_ID', 'PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ', 'POLYA_LENGTH', 'TMD', 
#                                'SP', 'n_exons', 'DSB', 'GPI', 'OG', 'LOCATION'])

# hgnc_id, premrna_seq = 'HGNC:TOY', ''.join(random.choices(['U', 'C', 'G', 'A'], k = 100))
# mrna_seq = premrna_seq[25:75]
# # note that there is no check that the protein_sequence corresponds to the mrna_sequence beyond checking for the length
# protein_seq = ''.join(random.choices(params.amino_acids, k = int(len(mrna_seq)/3)))
# polyA_length, tmd, sp, n_exons, dsb, gpi, og  = None, 1, True, None, 2, 2, 2
# ub_args = ubiquitin.express_ubiquitin(compress_mrna = False)

# import itertools
# reactions = list()
# for l in list(itertools.combinations(params.compartments.keys(),2)):
#     location = list(l)
#     psim_toy.loc[0,:] = [hgnc_id, premrna_seq, mrna_seq, protein_seq, polyA_length, tmd, sp, n_exons, dsb, gpi, og, location]
#     gene_info = gene_information(hgnc_id, premrna_seq, mrna_seq, protein_seq,
#                      ptms = {}, tmd = tmd, sp = sp, polyA_length = polyA_length, 
#                      n_exons = n_exons) 
#     gene_info.get_final_locations(reactions = None, nonmachinery_locations = location)

#     transcription_reactions, mrna_transcript_c, mrna_deg_proxy = build_mrna.get_mrna_expression_reactions(gene_info)
#     protein_expression_reactions, protein_metabolites = get_protein_expression_reactions(gene_info, 
#                                                      mrna_transcript_c, mrna_deg_proxy, 
#                                                     ub_args = ub_args)
#     reactions += protein_expression_reactions


