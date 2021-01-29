#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra

import sys
sys.path.insert(1, '../../../scripts/')
from macromolecules.protein import Protein
from macromolecules.complex import Complex

from utils import machinery as mach
from utils import parameters as params
from utils import metabolites as metab
from utils import functions as func


# # Cytosol and Nucleus

# In[4]:


def protein_polyubiquitination(macromolecule, **kwargs):
    ub_args = kwargs['ub_args']

    polyubiquitinate_protein = cobra.Reaction(macromolecule.id + '_POLYUBIQUITINATION' + macromolecule.compartment)
    polyubiquitinate_protein.subsytem = 'Protein_Expression'
    
    if macromolecule.compartment in ['c', 'n']:
        polyu_protein_aa_counts = macromolecule.amino_acid_counts.copy()
        for aa_code,aa_counts in ub_args['monoub_aa_counts'].items():
            if aa_code in polyu_protein_aa_counts:
                polyu_protein_aa_counts[aa_code] += aa_counts*params.n_ub
            else: 
                polyu_protein_aa_counts[aa_code] = aa_counts*params.n_ub
        
        
        polyub_macromolecule = Protein(id_ = macromolecule.id + '_polyub', compartment = macromolecule.compartment,
                           amino_acid_counts = polyu_protein_aa_counts) 
        
        rxn = {macromolecule: -1, ub_args['ub_' + macromolecule.compartment]: -params.n_ub, 
               polyub_macromolecule:1, metab.h2o_compartments[macromolecule.compartment]: params.n_ub}
            # 1 ATP hydrolysis per ubiquitin monomer added (https://link.springer.com/article/10.1007/s10637-020-00894-6)
        rxn = func.hydrolyze_atp(rxn, n_atp = params.n_ub, compartment = macromolecule.compartment)
        polyubiquitinate_protein.add_metabolites(rxn)
        
        polyubiquitinate_protein.gene_reaction_rule = ' and '.join(mach.UB_ligases[macromolecule.compartment])
        return polyubiquitinate_protein, polyub_macromolecule
    
    elif macromolecule.compartment == 'pm':
        polyub_macromolecule_pm = macromolecule.copy() #macromolecule.change_compartment('pm')
        elements = polyub_macromolecule_pm.elements.copy()
        for aa_code, aa_count in ub_args['monoub_aa_counts'].items():
            aa_elements = metab.seq_amino_acid_map_c[aa_code].elements
            for element in aa_elements:
                elements[element] += aa_count*aa_elements[element]*params.n_ub
            polyub_macromolecule_pm.charge += metab.seq_amino_acid_map_c[aa_code].charge*aa_count
        # peptide bond formation
        elements['H'] -= 2*(ub_args['L_monoub']*params.n_ub) # no -1 bc already accounted for in copying elements
        elements['O'] -= 1*(ub_args['L_monoub']*params.n_ub)
        polyub_macromolecule_pm.elements = elements

        rxn = {macromolecule: -1, ub_args['ub_c']: -params.n_ub, polyub_macromolecule_pm: 1, metab.h2o_c: params.n_ub}
        # 1 ATP hydrolysis per ubiquitin monomer added (https://link.springer.com/article/10.1007/s10637-020-00894-6)
        rxn = func.hydrolyze_atp(rxn, n_atp = params.n_ub, compartment = 'c')

        polyubiquitinate_protein.add_metabolites(rxn)
        polyubiquitinate_protein.gene_reaction_rule = ' and '.join(mach.UB_ligases_c + mach.HSP70_c + mach.HSP90AB1)
        
        return polyubiquitinate_protein, polyub_macromolecule_pm
    
    else:
        raise ValueError(macromolecule.id + ': Current compartment, ' + macromolecule.compartment + ' does not have polyubiquitination reactions available')

def proteasomal_degradation(macromolecule, **kwargs):
    polyub_macromolecule = kwargs['polyub_macromolecule']
    ub_args = kwargs['ub_args']
    
    if macromolecule.compartment != polyub_macromolecule.compartment:
        raise ValueError('Polyubiquitinated and deubiquitinated protein compartment does not match')
    if macromolecule.compartment not in ['c', 'n']:
        raise ValueError(macromolecule.id + ': Only proteins/complexes in nucleus and cytosol are considered for proteasomal degradation')
    #------------------------------deubiquitination------------------------------------  
    deubiquitination = cobra.Reaction(macromolecule.id + '_DEUBIQUITINATION' + macromolecule.compartment)
    deubiquitination.subsytem = 'Protein_Expression'
    deubiquitination.add_metabolites({polyub_macromolecule: -1, 
                                      metab.h2o_compartments[macromolecule.compartment]: -1, 
                                        macromolecule: 1, ub_args['polyub_' + macromolecule.compartment]: 1})
    
    deubiquitination.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)
    
    #------------------------------degradation------------------------------------
    protein_degradation = cobra.Reaction(macromolecule.id + '_PROTEASOMAL_DEGRADATION' + macromolecule.compartment)
    protein_degradation.subsytem = 'Protein_Expression'
    
    rxn = {metab.seq_amino_acid_map_compartments[macromolecule.compartment][aa_code]: aa_counts for aa_code, aa_counts in macromolecule.amino_acid_counts.items()}
    rxn[polyub_macromolecule] = -1
    rxn[metab.h2o_compartments[macromolecule.compartment]] = -macromolecule.L_protein
    rxn[ub_args['polyub_' + macromolecule.compartment]] =  1
    # atp hydrolysis for translocation/unfolding  - known 1 ATP per 2 residues - https://www.nature.com/articles/s41586-018-0736-4
    rxn = func.hydrolyze_atp(rxn, n_atp = macromolecule.L_protein/2, 
                             compartment = macromolecule.compartment)

    protein_degradation.add_metabolites(rxn)
    protein_degradation.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)

    return [deubiquitination, protein_degradation]

def degrade_cytosolic_nuclear_protein(macromolecule, **kwargs):
    '''Degradation reactions for cytosolic or nuclear proteins'''
    ub_args = kwargs['ub_args']
    polyubiquitinate_macromolecule, polyub_macromolecule = protein_polyubiquitination(macromolecule = macromolecule, 
                                                                                  ub_args = ub_args) 
    proteasomal_degradation_reactions = proteasomal_degradation(macromolecule = macromolecule, 
                                                  polyub_macromolecule = polyub_macromolecule, 
                                                  ub_args = ub_args)
    
    return [polyubiquitinate_macromolecule] + proteasomal_degradation_reactions


# # Mitochondria and Intermembrane Space

# In[ ]:


def degrade_mitochondrial_protein(macromolecule):
    rxn = {metab.seq_amino_acid_map_m[aa_code]: aa_counts for aa_code, aa_counts in macromolecule.amino_acid_counts.items()}
    rxn[macromolecule], rxn[metab.h2o_m] = -1, -(macromolecule.L_protein-1)
    
    if macromolecule.compartment == 'm':
        mitochondrial_degradation = cobra.Reaction(macromolecule.hgnc_id + '_DEGRADATIONm')
        mitochondrial_degradation.gene_reaction_rule = mach.mLON[0]
        
        # ATP hydrolysis by LON: 2 ATP per residue - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518814/
        rxn = func.hydrolyze_atp(rxn, n_atp = macromolecule.L_protein*2, compartment = 'm')
        

    elif macromolecule.compartment == 'i':
        mitochondrial_degradation = cobra.Reaction(macromolecule.hgnc_id + '_DEGRADATIONi')
        mitochondrial_degradation.gene_reaction_rule = mach.iAAA[0]#' and '.join(mAAA + iAAA)
        # in the future, may want to add ubqituin-proteasome: 
        
        # ATP hydrolysis by m/i-AAA: 1 ATP per 2 residues -- no source, assumes same as 26S proteasome
        rxn = func.hydrolyze_atp(rxn, n_atp = macromolecule.L_protein*params.proteolysis_translocation_atp_cost, compartment = 'i')
  
    mitochondrial_degradation.subsytem = 'Protein_Expression'
    mitochondrial_degradation.add_metabolites(rxn)
    
    return mitochondrial_degradation


def degrade_peroxisomal_protein(macromolecule):
    
    peroxisomal_degradation = cobra.Reaction(macromolecule.hgnc_id + '_DEGRADATIONx')
    peroxisomal_degradation.subsytem = 'Protein_Expression'
    peroxisomal_degradation.gene_reaction_rule = mach.LONP2[0]

    rxn = {metab.seq_amino_acid_map_x[aa_code]: aa_counts for aa_code, aa_counts in macromolecule.amino_acid_counts.items()}
    rxn[macromolecule], rxn[metab.h2o_x] = -1, -(macromolecule.L_protein-1)
    # ATP hydrolysis by LON: 2 ATP per residue - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518814/
    rxn = func.hydrolyze_atp(rxn, n_atp = macromolecule.L_protein*2, compartment = 'x')
    peroxisomal_degradation.add_metabolites(rxn)
    
    return peroxisomal_degradation


# # Secretory Pathway Degradation

# In[ ]:


def unfold_secretory_protein(macromolecule):
    '''Remove PTMs and unfold proteins for lysosomal and secretory compartments. For lysosomal degradation, 
    only remove PTMs, there is no unfolding/misfolding as in ERAD.'''
    
    if not (macromolecule.compartment == 'r' or macromolecule.compartment == 'l'):
        raise ValueError('Protein metabolite does not have correct compartment')
    
    
    unfold_protein = cobra.Reaction(macromolecule.hgnc_id + '_UNFOLD' + macromolecule.compartment)
    unfolded_protein = macromolecule.copy()
    if macromolecule.compartment == 'r': # not "unfolding/misfolding" for lysosomal degradation
        unfolded_protein.id = unfolded_protein.id.replace('folded', 'unfolded')
    
    rxn = dict()
    rxn[macromolecule] = -1

    unfold_mach = list()
    elements = unfolded_protein.elements.copy()
#     lysosomal_degradation_ptm_condition = 'gpi' in macromolecule.ptms.keys() and len(macromolecule.ptms.keys()) == 1

    # PTM removals HERE #YOU ARE HERE 
    if 'ng' in macromolecule.ptms.keys():
        raise ValueError('N-glycosylation not yet incorporated')
#     if lysosomal_degradation_ptm_condition:
#         raise ValueError('GPI-anchored proteins with no other ptms should be degraded via lysosomal pathway')
    if 'gpi' in macromolecule.ptms.keys():
        unfolded_protein.id = unfolded_protein.id.replace('_GPI', '')
#         for e,c in metab.balanced_gpi.items():
#             elements[e] -= c
        
#         for m,s in metab.M4ATAer.copy().items(): # no lysosomal compartment metabolites
#             rxn[m] = -s
        rxn[metab.gpi_hs_r] = 1
        if macromolecule.compartment == 'r':
            rxn[metab.hdca_r], rxn[metab.h_r], rxn[metab.h2o_r] = -1,-1,1
        elif macromolecule.compartment == 'l':
            rxn[metab.hdca_l], rxn[metab.h_l], rxn[metab.h2o_l] = -1,-1,1
        

    if 'dsb' in macromolecule.ptms.keys():
        unfolded_protein.id = unfolded_protein.id.replace('_DSB', '')
        number_DSB = macromolecule.ptms['dsb']
        elements['H'] += 2*number_DSB
        # incorporate exchange with reductase in future versions
        if macromolecule.compartment == 'r':
            rxn[metab.o2_r], rxn[metab.h2o2_r] = number_DSB, -number_DSB
        elif macromolecule.compartment == 'l':
            rxn[metab.o2_l], rxn[metab.h2o2_l] = number_DSB, -number_DSB
        unfold_mach += mach.ERDJ5
    if 'og' in macromolecule.ptms.keys():
        unfolded_protein.id = unfolded_protein.id.replace('_OG', '')
        number_Oglycans = macromolecule.ptms['og']
        balance_og = {'C': (8 + 6 + 8)*number_Oglycans, # each 1/3 entry is for each 1/3 reactions in Jahir's model, in case want to separate in the future 
                      'H': (13 + 10 + 13)*number_Oglycans, 
                      'N': (1 + 0 + 1)*number_Oglycans, 
                      'O': (5 + 5 + 5)*number_Oglycans}
        for e,c in balance_og.items():
            elements[e] -= c
        
        if macromolecule.compartment == 'r':
            rxn[metab.udpacgal_r], rxn[metab.udpgal_r], rxn[metab.uacgam_r] = number_Oglycans, number_Oglycans, number_Oglycans
            rxn[metab.udp_r] = -3*number_Oglycans
            if metab.h_r in rxn.keys():
                rxn[metab.h_r] -= 3*number_Oglycans
            else:
                rxn[metab.h_r] = -3*number_Oglycans
        elif macromolecule.compartment == 'l':
            rxn[metab.udpacgal_l], rxn[metab.udpgal_g], rxn[metab.uacgam_g] = number_Oglycans, number_Oglycans, number_Oglycans
            rxn[metab.udp_l] = -3*number_Oglycans
            if metab.h_l in rxn.keys():
                rxn[metab.h_l] -= 3*number_Oglycans
            else:
                rxn[metab.h_l] = -3*number_Oglycans

    ###########
    
    unfolded_protein.elements = elements
    rxn[unfolded_protein] = 1
        
    unfold_protein.add_metabolites(rxn)
    if len(unfold_mach) > 1:
        unfold_protein.gene_reaction_rule = ' and '.join(unfold_mach)
    elif len(unfold_mach) == 1:
        unfold_protein.gene_reaction_rule = unfold_mach[0]
    
    # PTMs - merged with unfolding reaction for now  <--same structure as adding the PTMs
    #     if 'ng' in macromolecule.ptms.keys():
    #         raise ValueError('N-glycosylation not yet incorporated')
    #     if 'dsb' in macromolecule.ptms.keys():
    #     else:
    #         unmodified_protein_r = unfolded_protein
    #     if 'gpi' in macromolecule.ptms.keys():
    #     else:
    #         unmodified_protein_r = unmodified_protein_r
    #     if 'og' in macromolecule.ptms.keys()
    #     else:
    #         unmodified_protein_r = unmodified_protein_r
    #     if len(macromolecule.ptms) == 0:
    #         unmodified_protein_r = unfolded_protein

    unmodified_protein = unfolded_protein # for adding PTMs as separate reactions in future, if want to
    
    return unfold_protein, unmodified_protein

def retrograde_er(macromolecule):
    if macromolecule.compartment != 'g':
        raise ValueError('ER retrograde transport can only occur for Golgi macromolecules')
        
    V = macromolecule.formula_weight/1000 * 1.21 / 1000.0 # Protein Volume in nm^3
    copi_coeff = int(round(143793.19 * params.Kv / V))

    retro_protein_r = macromolecule.change_compartment('r')

    rxn = {macromolecule: -copi_coeff, retro_protein_r: copi_coeff}
    # gtp hydrolysis
    rxn[metab.ntp_map_c['G']], rxn[metab.h2o_c], rxn[metab.ndp_map_c['G']], rxn[metab.pi_c], rxn[metab.h_c]  = -127, -127, 127, 127, 127

    retrograde_transport = cobra.Reaction(macromolecule.hgnc_id + '_COPI_RETROtr')
    retrograde_transport.add_metabolites(rxn)
    retrograde_transport.gene_reaction_rule = ' and '.join(mach.copi_m)
    
    return retrograde_transport, retro_protein_r

def build_erad_reactions(macromolecule, **kwargs):
    if 'unfolded_protein_c' in kwargs.keys():
        unfolded_protein_c = kwargs['unfolded_protein_c']
    else:
        unfolded_protein_c = None
        
    if macromolecule.compartment == 'g':
        raise ValueError('This situation is only for complexes, see commented out code below')
        #retrograde_transport, macromolecule = retrograde_er(macromolecule)
    else:
        retrograde_transport = None
    
    if macromolecule.compartment != 'r':
        raise ValueError('ERAD can only occur with proteins in ER compartment')
    unfold_er_protein, unmodified_protein_r = unfold_secretory_protein(macromolecule)
    

    # Retro-translocation
    retrotranslocate_protein = cobra.Reaction(macromolecule.hgnc_id + '_RETROTRANSLOCATION')
    if unfolded_protein_c == None: # those that underwent co-translational rather than post-translational
        unfolded_protein_c = unmodified_protein_r.change_compartment('c')
        # replace unfolded bc, if it underwent co-translational translocation
        # the signal peptide degradation step makes it such that it is not the same as any 
        # cytosolically translated unfolded proteins (multi-localization), i.e., 'i' destined proteins
        # and their degradation
        unfolded_protein_c.id = unfolded_protein_c.id.replace('unfolded', 'retrotranslocated_unfolded')
    
    rxn = {unmodified_protein_r: -1, unfolded_protein_c: 1}
    
    # FUTURE: separate nonglycosylated and glycosylated ERAD
    # if 'ng' in macromolecule.ptms.keys() or 'og' in macromolecule.ptms.keys():
    
     # glyco based ERAD from Jahir
    rxn = func.hydrolyze_atp(rxn, n_atp = 6, compartment = 'c') # from Jahir retro_TRANSLOC_2
    retrotranslocate_protein.add_metabolites(rxn)
    retrotranslocate_protein.gene_reaction_rule = ' and '.join(mach.retro_mach_glyco)
    

    erad_reactions = [unfold_er_protein, retrotranslocate_protein]
    if retrograde_transport is not None:
        erad_reactions += [retrograde_transport]

    return erad_reactions, unfolded_protein_c

# def degrade_er_golgi(macromolecule, unfolded_protein_c = None, ptt = None):
    
#     if ptt is None:
#         if macromolecule.L_protein > params.ptt_length:
#             ptt_ = False
#         else:
#             ptt_ = True
    
#     if ptt_:
#         if unfolded_protein_c is None:
#             print(macromolecule.hgnc_id)
#             print(macromolecule.L_protein)
#             raise ValueError('Need an unfolded_protein_c metabolite')
#         erad_reactions, unfolded_protein_c = build_erad_reactions(macromolecule, 
#                                                                   unfolded_protein_c = unfolded_protein_c)
#     else:
#         erad_reactions, unfolded_protein_c = build_erad_reactions(macromolecule, 
#                                                                           unfolded_protein_c = None)
#     return erad_reactions, unfolded_protein_c


# In[ ]:


def build_endocytosis_reactions(macromolecule_pm, **kwargs):
    ub_args = kwargs['ub_args']
    if 'macromolecule_l' in kwargs.keys():
        macromolecule_l = kwargs['macromolecule_l']
    else:
        macromolecule_l = None    
    
    ##polyubiquitination for lysosomal targetting-------------------------------------
    polyubiquitinate_protein, polyub_macromolecule_pm = protein_polyubiquitination(macromolecule = macromolecule_pm, 
                                                                                         ub_args = ub_args)
    
    polyub_macromolecule_l = polyub_macromolecule_pm.change_compartment('l')
    
    
    ##endocytosis--------------------------------------------------------------------------
    # combine dequbiquitination with endocytosis https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3987138/
    if macromolecule_l == None:
        macromolecule_l = macromolecule_pm.change_compartment('l')
    else:
        if macromolecule_l.L_protein != macromolecule_pm.L_protein:
            raise ValueError('Endocytosis: lysosomal and plasma membrane protein lengths disagree')
        

    endocytosis = cobra.Reaction(macromolecule_pm.hgnc_id + '_CLATHRIN_ENDOCYTOSIS')    
    
    rxn = {polyub_macromolecule_pm: -1, macromolecule_l: 1, metab.h2o_c: -1, ub_args['polyub_c']: 1}
    
    # gtp hydrolysis for vesicle scission
    rxn[metab.ntp_map_c['G']] = -round(macromolecule_pm.L_protein) * params.transport_translocation_atp_cost
    rxn[metab.h2o_c] -= round(macromolecule_pm.L_protein) * params.transport_translocation_atp_cost
    rxn[metab.ndp_map_c['G']] = round(macromolecule_pm.L_protein) * params.transport_translocation_atp_cost
    rxn[metab.pi_c] = round(macromolecule_pm.L_protein) * params.transport_translocation_atp_cost
    rxn[metab.h_c] = round(macromolecule_pm.L_protein) * params.transport_translocation_atp_cost

    
    endocytosis.add_metabolites(rxn)
    endocytosis.gene_reaction_rule = ' and '.join(mach.endocytic_machinery)
    
    return [polyubiquitinate_protein, endocytosis], macromolecule_l

def lysosomal_degradation(macromolecule):
    if macromolecule.compartment != 'l':
        raise ValueError('This reaction only occurs in the lysosome')
    lysosomal_degradation_reactions = list()
    
    if len(macromolecule.ptms.keys())>0:
        # not unfolding lysosomal protein is just removing the PTMs, not an actual unfolding reaction
        unfold_lysosomal_protein, unmodified_macromolecule = unfold_secretory_protein(macromolecule)
        lysosomal_degradation_reactions += [unfold_lysosomal_protein]
    else:
        unmodified_macromolecule = macromolecule
    
    
    degrade_lysosomal_protein = cobra.Reaction(macromolecule.hgnc_id + '_LYSOSOMAL_DEGRADATION')
    
    rxn = {metab.seq_amino_acid_map_l[aa_code]: aa_counts for aa_code, aa_counts in macromolecule.amino_acid_counts.items()}
    rxn[unmodified_macromolecule], rxn[metab.h2o_l] = -1, -(macromolecule.L_protein-1)
    rxn = func.hydrolyze_atp(rxn, n_atp=macromolecule.L_protein*params.proteolysis_translocation_atp_cost, compartment = 'l')
    
    degrade_lysosomal_protein.add_metabolites(rxn)
    degrade_lysosomal_protein.gene_reaction_rule = ' and '.join(mach.cathepsins)
    lysosomal_degradation_reactions += [degrade_lysosomal_protein]

    return lysosomal_degradation_reactions

def degrade_lysosomal_pm_protein(macromolecule, **kwargs):
    if 'ub_args' in kwargs.keys():
        ub_args = kwargs['ub_args']
    else:
        ub_args = None  
    if macromolecule.compartment == 'pm':
        raise ValueError('This is for complexes in the future (see commented code below)')
        # endocytosis_reactions, macromlecule = build_endocytosis_reactions(macromolecule_pm = protein_pm, 
                                                                               #protein_l = protein_l,
                                                                               # ub_args = ub_args)   
    else:
        endocytosis_reactions = None
    if macromolecule.compartment != 'l':
        raise ValueError('Must be a lysosomal protein/complex being degraded')
    
    deg_reactions = lysosomal_degradation(macromolecule)
    if endocytosis_reactions is not None:
        deg_reactions += endocytosis_reactions
    return deg_reactions
    
    


# In[3]:


degrade_map = {'c': degrade_cytosolic_nuclear_protein, 'n': degrade_cytosolic_nuclear_protein, 
              'm': degrade_mitochondrial_protein, 'i': degrade_mitochondrial_protein, 
              'x': degrade_peroxisomal_protein, 
               'r': build_erad_reactions, 'g': build_erad_reactions, 
              'l': degrade_lysosomal_pm_protein, 'pm': degrade_lysosomal_pm_protein}

def degrade(macromolecule, **kwargs):
    '''Compartment-specific degradation reactions for proteins or protein-protein complexes'''
    
    if not type(macromolecule) == Protein or type(macromolecule) == Complex:
        raise ValueError('Macromolecule to degrade must be protein or complex')
    
    deg_reactions = degrade_map[macromolecule.compartment](macromolecule, **kwargs)
    return deg_reactions

