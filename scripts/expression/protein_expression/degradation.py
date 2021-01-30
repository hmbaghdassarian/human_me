#!/usr/bin/env python
# coding: utf-8

# In[48]:


import cobra

import sys
sys.path.insert(1, '../../../scripts/')
from utils import machinery as mach
from utils import parameters as params
from utils import metabolites as metab
from utils import functions as func

from macromolecules.protein import Protein
from macromolecules.complex import Complex
from macromolecules.complex import add_complex_metabolites


# In[49]:


class Protein_Degradation_Reaction(cobra.Reaction):
    def __init__(self, id=None, name='', subsystem='', lower_bound=0.0, upper_bound=None):
        cobra.Reaction.__init__(self, id=id, name=name, subsystem=subsystem, lower_bound=lower_bound, 
                                upper_bound=upper_bound)
        self._macromolecules = [] # list of macromolecule ids associated with this degradation reaction
        self._enzymes = None # list of enzyme ids associated with this degradation reaction
        self.final_reaction = False # whether the reaction is the "final" (to amino acids) degradation reaction
        self.subsystem = 'Protein_Degradation'
        
    def _update_tracking(self, macromolecules):
        '''Mutual tracking of degradation reactions associated with a macromolecule and vice-versa'''
        if type(macromolecules) != list:
            macromolecules._degradation_reactions.append(self.id)
            self._macromolecules.append(macromolecules)
        else:
            for macromolecule in macromolecules:
                macromolecule._degradation_reactions.append(self.id)
                self._macromolecules.append(macromolecule)
    def _consolidate_macromolecules(self):
        '''Remove redundant macromolecules'''
        for m in self._macromolecules:
            m._consolidate_degradation_rxns()
        self._macromolecules = list(set(self._macromolecules))
        
    def _update_enzymes(self):
        '''Update enzymes list to include macromolecules that are classified as enzymes'''
        self._enzymes = [m for m in self._macromolecules if m.enzyme]
        for m in self.enzymes:
            if self.id not in m._degradation_reactions:
                raise ValueErorr('Improper tracking of degradation reactions and associated macromolecules')

class Complex_Degradation_Reaction(cobra.Reaction):
    def __init__(self, id=None, name='', subsystem='', lower_bound=0.0, upper_bound=None):
        cobra.Reaction.__init__(self, id=id, name=name, subsystem=subsystem, lower_bound=lower_bound, 
                                upper_bound=upper_bound)
        self._macromolecules = [] # list of macromolecule ids associated with this degradation reaction
        self._enzymes = None # list of enzyme ids associated with this degradation reaction
        self.final_reaction = False # whether the reaction is the "final" (to amino acids) degradation reaction
        self.subsystem = 'Complex_Degradation'
        
    def _update_tracking(self, macromolecules):
        '''Mutual tracking of degradation reactions associated with a macromolecule and vice-versa'''
        if type(macromolecules) != list:
            macromolecules._degradation_reactions.append(self.id)
            self._macromolecules.append(macromolecules)
        else:
            for macromolecule in macromolecules:
                macromolecule._degradation_reactions.append(self.id)
                self._macromolecules.append(macromolecule)
    def _consolidate_macromolecules(self):
        '''Remove redundant macromolecules'''
        for m in self._macromolecules:
            m._consolidate_degradation_rxns()
        self._macromolecules = list(set(self._macromolecules))
        
    def _update_enzymes(self):
        '''Update enzymes list to include macromolecules that are classified as enzymes'''
        self._enzymes = [m for m in self._macromolecules if m.enzyme]
        for m in self.enzymes:
            if self.id not in m._degradation_reactions:
                raise ValueErorr('Improper tracking of degradation reactions and associated macromolecules')

deg_reaction_map = {'protein': Protein_Degradation_Reaction, 'complex': Complex_Degradation_Reaction}


# # Cytosol and Nucleus

# In[50]:


def protein_polyubiquitination(macromolecule, **kwargs):
    ub_args = kwargs['ub_args']

    polyubiquitinate_protein = deg_reaction_map[macromolecule.type](macromolecule.id + '_POLYUBIQUITINATION' + macromolecule.compartment)
    cmap = {'n': 'n', 'c': 'c', 'pm': 'c'}
    if macromolecule.compartment not in cmap.keys():
        raise ValueError(macromolecule.id + ': Current compartment, ' + macromolecule.compartment + ' does not have polyubiquitination reactions available')
    
    if macromolecule.type == 'protein':
        polyu_protein_aa_counts = macromolecule._amino_acid_counts.copy()
        for aa_code,aa_counts in ub_args['monoub_aa_counts'].items():
            if aa_code in polyu_protein_aa_counts:
                polyu_protein_aa_counts[aa_code] += aa_counts*params.n_ub
            else: 
                polyu_protein_aa_counts[aa_code] = aa_counts*params.n_ub
        polyub_macromolecule = Protein(id_ = '_'.join(macromolecule.id.split('_')[:-1]) + '_polyub', 
                                       compartment = cmap[macromolecule.compartment],
                           amino_acid_counts = polyu_protein_aa_counts) 
        if macromolecule.compartment == 'pm':
            polyub_macromolecule.change_compartment('pm')
    else:
        fused_poly_ub = ub_args['polyub_' + cmap[macromolecule.compartment]].copy()
        fused_poly_ub.id = fused_poly_ub.id.replace('cleaved_', 'fused_')
        # peptide bond formation
        elements = fused_poly_ub.elements.copy()
        elements['H'] -= 2
        elements['O'] -= 1
        fused_poly_ub.elements = elements
        fused_poly_ub.compartment = macromolecule.compartment
        
        polyub_macromolecule = add_complex_metabolites(cplx = macromolecule, 
                               met_to_add = {fused_poly_ub: 1}, 
                                complex_id = macromolecule.temp_id + '_polyub')

    rxn = {macromolecule: -1, ub_args['ub_' + cmap[macromolecule.compartment]]: -params.n_ub, 
           polyub_macromolecule:1, metab.h2o_compartments[cmap[macromolecule.compartment]]: params.n_ub}
    # 1 ATP hydrolysis per ubiquitin monomer added (https://link.springer.com/article/10.1007/s10637-020-00894-6)
    rxn = func.hydrolyze_atp(rxn, n_atp = params.n_ub, compartment = cmap[macromolecule.compartment])
    polyubiquitinate_protein.add_metabolites(rxn)

    if macromolecule.compartment in ['c', 'n']:
        polyubiquitinate_protein.gene_reaction_rule = ' and '.join(mach.UB_ligases[macromolecule.compartment])
    else:
        polyubiquitinate_protein.gene_reaction_rule = ' and '.join(mach.UB_ligases_c + mach.HSP70_c + mach.HSP90AB1)

    polyubiquitinate_protein._update_tracking(macromolecules = [macromolecule, polyub_macromolecule])
    return polyubiquitinate_protein, polyub_macromolecule

def proteasomal_degradation(macromolecule, **kwargs):
    polyub_macromolecule = kwargs['polyub_macromolecule']
    ub_args = kwargs['ub_args']
    
    if macromolecule.compartment != polyub_macromolecule.compartment:
        raise ValueError('Polyubiquitinated and deubiquitinated protein compartment does not match')
    if macromolecule.compartment not in ['c', 'n']:
        raise ValueError(macromolecule.id + ': Only proteins/complexes in nucleus and cytosol are considered for proteasomal degradation')
    #------------------------------deubiquitination------------------------------------  
    deubiquitination = deg_reaction_map[macromolecule.type](macromolecule.id + '_DEUBIQUITINATION' + macromolecule.compartment)
    deubiquitination.add_metabolites({polyub_macromolecule: -1, 
                                      metab.h2o_compartments[macromolecule.compartment]: -1, 
                                        macromolecule: 1, ub_args['polyub_' + macromolecule.compartment]: 1})
    
    deubiquitination.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)
    
    #------------------------------degradation------------------------------------
    protein_degradation = deg_reaction_map[macromolecule.type](macromolecule.id + '_PROTEASOMAL_DEGRADATION' + macromolecule.compartment)
    
    h2o_length = macromolecule._L_protein
    if macromolecule.type == 'protein':
        aac = macromolecule._amino_acid_counts.copy()
    else:
        aac = polyub_macromolecule._amino_acid_counts.copy()
        aac.subtract(ub_args['polyub_' + macromolecule.compartment]._amino_acid_counts)
        h2o_length -= sum(macromolecule.decompose_complex().values()) # non-covalent bonds don't require h2o, 
        h2o_length += 1 # 1 required for fused polyub
        
    
    rxn = {metab.seq_amino_acid_map_compartments[macromolecule.compartment][aa_code]: aa_counts for aa_code, aa_counts in aac.items()}
    rxn[polyub_macromolecule] = -1
    print(macromolecule._L_protein)
    rxn[metab.h2o_compartments[macromolecule.compartment]] = -h2o_length
    rxn[ub_args['polyub_' + macromolecule.compartment]] =  1
    # atp hydrolysis for translocation/unfolding  - known 1 ATP per 2 residues - https://www.nature.com/articles/s41586-018-0736-4
    rxn = func.hydrolyze_atp(rxn, n_atp = macromolecule._L_protein/2, 
                             compartment = macromolecule.compartment)

    protein_degradation.add_metabolites(rxn)
    protein_degradation.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)
    
    # tracking
    protein_degradation.final_reaction = True
    for r in [deubiquitination, protein_degradation]:
        r._update_tracking(macromolecules = [macromolecule, polyub_macromolecule])
    
    return [deubiquitination, protein_degradation]

def degrade_cytosolic_nuclear_protein(macromolecule, **kwargs):
    '''Degradation reactions for cytosolic or nuclear proteins'''
    ub_args = kwargs['ub_args']
    polyubiquitinate_macromolecule, polyub_macromolecule = protein_polyubiquitination(macromolecule = macromolecule, 
                                                                                  ub_args = ub_args) 
    proteasomal_degradation_rxns = proteasomal_degradation(macromolecule = macromolecule, 
                                                  polyub_macromolecule = polyub_macromolecule, 
                                                  ub_args = ub_args)
    
    return [polyubiquitinate_macromolecule] + proteasomal_degradation_rxns


# # Mitochondria and Intermembrane Space

# In[51]:


def degrade_mitochondrial_protein(macromolecule):
    rxn = {metab.seq_amino_acid_map_m[aa_code]: aa_counts for aa_code, aa_counts in macromolecule._amino_acid_counts.items()}
    
    h2o_length = macromolecule._L_protein - 1 if macromolecule.type == 'protein' else     macromolecule._L_protein - sum(macromolecule.decompose_complex().values()) #non-covalent bonds, +1,-1 cancel out
    rxn[macromolecule], rxn[metab.h2o_m] = -1, -h2o_length
    
    if macromolecule.compartment == 'm':
        mitochondrial_degradation = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_DEGRADATIONm')
        mitochondrial_degradation.gene_reaction_rule = mach.mLON[0]
        
        # ATP hydrolysis by LON: 2 ATP per residue - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518814/
        rxn = func.hydrolyze_atp(rxn, n_atp = macromolecule._L_protein*2, compartment = 'm')
        

    elif macromolecule.compartment == 'i':
        mitochondrial_degradation = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_DEGRADATIONi')
        mitochondrial_degradation.gene_reaction_rule = mach.iAAA[0]#' and '.join(mAAA + iAAA)
        # in the future, may want to add ubqituin-proteasome: 
        
        # ATP hydrolysis by m/i-AAA: 1 ATP per 2 residues -- no source, assumes same as 26S proteasome
        rxn = func.hydrolyze_atp(rxn, n_atp = macromolecule._L_protein*params.proteolysis_translocation_atp_cost, compartment = 'i')
  
    mitochondrial_degradation.add_metabolites(rxn)
    
    mitochondrial_degradation.final_reaction = True
    mitochondrial_degradation._update_tracking(macromolecule)
    
    return [mitochondrial_degradation]


def degrade_peroxisomal_protein(macromolecule):
    
    peroxisomal_degradation = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_DEGRADATIONx')
    peroxisomal_degradation.gene_reaction_rule = mach.LONP2[0]
    
    h2o_length = macromolecule._L_protein - 1 if macromolecule.type == 'protein' else     macromolecule._L_protein - sum(macromolecule.decompose_complex().values()) #non-covalent bonds +1,-1 cancel out
    
    rxn = {metab.seq_amino_acid_map_x[aa_code]: aa_counts for aa_code, aa_counts in macromolecule._amino_acid_counts.items()}
    rxn[macromolecule], rxn[metab.h2o_x] = -1, -h2o_length
    # ATP hydrolysis by LON: 2 ATP per residue - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518814/
    rxn = func.hydrolyze_atp(rxn, n_atp = macromolecule._L_protein*2, compartment = 'x')
    peroxisomal_degradation.add_metabolites(rxn)
    
    peroxisomal_degradation.final_reaction = True
    peroxisomal_degradation._update_tracking(macromolecule)
    
    return [peroxisomal_degradation]


# # Secretory Pathway Degradation

# In[52]:


def unfold_secretory_protein(macromolecule):
    '''Remove PTMs and unfold proteins for lysosomal and secretory compartments. For lysosomal degradation, 
    only remove PTMs, there is no unfolding/misfolding as in ERAD.'''
    
    if not (macromolecule.compartment == 'r' or macromolecule.compartment == 'l'):
        raise ValueError('Protein metabolite does not have correct compartment')
    
    
    unfold_protein = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_UNFOLD' + macromolecule.compartment)
    unfolded_protein = macromolecule.copy()
    if macromolecule.compartment == 'r': # not "unfolding/misfolding" for lysosomal degradation
        unfolded_protein.id = unfolded_protein.id.replace('folded', 'unfolded')
    
    rxn = dict()
    rxn[macromolecule] = -1

    unfold_mach = list()
    elements = unfolded_protein.elements.copy()
#     lysosomal_degradation_ptm_condition = 'gpi' in macromolecule._ptms.keys() and len(macromolecule._ptms.keys()) == 1

    # PTM removals HERE #YOU ARE HERE 
    if 'ng' in macromolecule._ptms.keys():
        raise ValueError('N-glycosylation not yet incorporated')
#     if lysosomal_degradation_ptm_condition:
#         raise ValueError('GPI-anchored proteins with no other ptms should be degraded via lysosomal pathway')
    if 'gpi' in macromolecule._ptms.keys():
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
        

    if 'dsb' in macromolecule._ptms.keys():
        unfolded_protein.id = unfolded_protein.id.replace('_DSB', '')
        number_DSB = macromolecule._ptms['dsb']
        elements['H'] += 2*number_DSB
        # incorporate exchange with reductase in future versions
        if macromolecule.compartment == 'r':
            rxn[metab.o2_r], rxn[metab.h2o2_r] = number_DSB, -number_DSB
        elif macromolecule.compartment == 'l':
            rxn[metab.o2_l], rxn[metab.h2o2_l] = number_DSB, -number_DSB
        unfold_mach += mach.ERDJ5
    if 'og' in macromolecule._ptms.keys():
        unfolded_protein.id = unfolded_protein.id.replace('_OG', '')
        number_Oglycans = macromolecule._ptms['og']
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
    #     if 'ng' in macromolecule._ptms.keys():
    #         raise ValueError('N-glycosylation not yet incorporated')
    #     if 'dsb' in macromolecule._ptms.keys():
    #     else:
    #         unmodified_protein_r = unfolded_protein
    #     if 'gpi' in macromolecule._ptms.keys():
    #     else:
    #         unmodified_protein_r = unmodified_protein_r
    #     if 'og' in macromolecule._ptms.keys()
    #     else:
    #         unmodified_protein_r = unmodified_protein_r
    #     if len(macromolecule._ptms) == 0:
    #         unmodified_protein_r = unfolded_protein

    unmodified_protein = unfolded_protein # for adding PTMs as separate reactions in future, if want to
    unfold_protein._update_tracking(macromolecules = [macromolecule, unfolded_protein, unmodified_protein])
    
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

    retrograde_transport = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_COPI_RETROtr')
    retrograde_transport.add_metabolites(rxn)
    retrograde_transport.gene_reaction_rule = ' and '.join(mach.copi_m)
    retrograde_transport._update_tracking([macromolecule, retro_protein_r])
    
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
    retrotranslocate_protein = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_RETROTRANSLOCATION')
    if unfolded_protein_c == None: # those that underwent co-translational rather than post-translational
        unfolded_protein_c = unmodified_protein_r.change_compartment('c')
        # replace unfolded bc, if it underwent co-translational translocation
        # the signal peptide degradation step makes it such that it is not the same as any 
        # cytosolically translated unfolded proteins (multi-localization), i.e., 'i' destined proteins
        # and their degradation
        unfolded_protein_c.id = unfolded_protein_c.id.replace('unfolded', 'retrotranslocated_unfolded')
    
    rxn = {unmodified_protein_r: -1, unfolded_protein_c: 1}
    
    # FUTURE: separate nonglycosylated and glycosylated ERAD
    # if 'ng' in macromolecule._ptms.keys() or 'og' in macromolecule._ptms.keys():
    
     # glyco based ERAD from Jahir
    rxn = func.hydrolyze_atp(rxn, n_atp = 6, compartment = 'c') # from Jahir retro_TRANSLOC_2
    retrotranslocate_protein.add_metabolites(rxn)
    retrotranslocate_protein.gene_reaction_rule = ' and '.join(mach.retro_mach_glyco)
    

    erad_reactions = [unfold_er_protein, retrotranslocate_protein]
    if retrograde_transport is not None:
        erad_reactions.append(retrograde_transport)
    
    for r in erad_reactions:
        r._update_tracking([macromolecule, unmodified_protein_r, unfolded_protein_c])
    return erad_reactions, unfolded_protein_c


# In[53]:


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
        if macromolecule_l._L_protein != macromolecule_pm._L_protein:
            raise ValueError('Endocytosis: lysosomal and plasma membrane protein lengths disagree')
        

    endocytosis = deg_reaction_map[type(macromolecule_pm)](macromolecule_pm._deg_id + '_CLATHRIN_ENDOCYTOSIS')    
    
    rxn = {polyub_macromolecule_pm: -1, macromolecule_l: 1, metab.h2o_c: -1, ub_args['polyub_c']: 1}
    
    # gtp hydrolysis for vesicle scission
    rxn[metab.ntp_map_c['G']] = -round(macromolecule_pm._L_protein) * params.transport_translocation_atp_cost
    rxn[metab.h2o_c] -= round(macromolecule_pm._L_protein) * params.transport_translocation_atp_cost
    rxn[metab.ndp_map_c['G']] = round(macromolecule_pm._L_protein) * params.transport_translocation_atp_cost
    rxn[metab.pi_c] = round(macromolecule_pm._L_protein) * params.transport_translocation_atp_cost
    rxn[metab.h_c] = round(macromolecule_pm._L_protein) * params.transport_translocation_atp_cost

    
    endocytosis.add_metabolites(rxn)
    endocytosis.gene_reaction_rule = ' and '.join(mach.endocytic_machinery)
    
    for r in [polyubiquitinate_protein, endocytosis]:
        r._update_tracking([macromolecule_pm, macromolecule_l, polyub_macromolecule_pm])
    return [polyubiquitinate_protein, endocytosis], macromolecule_l

def lysosomal_degradation(macromolecule):
    if macromolecule.compartment != 'l':
        raise ValueError('This reaction only occurs in the lysosome')
    lysosomal_degradation_rxns = list()
    
    if len(macromolecule._ptms.keys())>0:
        # not unfolding lysosomal protein is just removing the PTMs, not an actual unfolding reaction
        unfold_lysosomal_protein, unmodified_macromolecule = unfold_secretory_protein(macromolecule)
        lysosomal_degradation_rxns += [unfold_lysosomal_protein]
    else:
        unmodified_macromolecule = macromolecule
    
    
    degrade_lysosomal_protein = deg_reaction_map[macromolecule.type](macromolecule._deg_id + '_LYSOSOMAL_DEGRADATION')
    
    rxn = {metab.seq_amino_acid_map_l[aa_code]: aa_counts for aa_code, aa_counts in macromolecule._amino_acid_counts.items()}
    rxn[unmodified_macromolecule], rxn[metab.h2o_l] = -1, -(macromolecule._L_protein-1)
    rxn = func.hydrolyze_atp(rxn, n_atp=macromolecule._L_protein*params.proteolysis_translocation_atp_cost, compartment = 'l')
    
    degrade_lysosomal_protein.add_metabolites(rxn)
    degrade_lysosomal_protein.gene_reaction_rule = ' and '.join(mach.cathepsins)
    degrade_lysosomal_protein.final_reaction = True
    
    lysosomal_degradation_rxns += [degrade_lysosomal_protein]
    
    for r in lysosomal_degradation_rxns:
        r._update_tracking([macromolecule, unmodified_macromolecule])
        
    return lysosomal_degradation_rxns

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
    
    for r in deg_reactions:
        r._update_tracking([macromolecule])
    return deg_reactions
    
    


# In[54]:


degrade_reaction_map = {'c': degrade_cytosolic_nuclear_protein, 'n': degrade_cytosolic_nuclear_protein, 
              'm': degrade_mitochondrial_protein, 'i': degrade_mitochondrial_protein, 
              'x': degrade_peroxisomal_protein, 
               'r': build_erad_reactions, 'g': build_erad_reactions, 
              'l': degrade_lysosomal_pm_protein, 'pm': degrade_lysosomal_pm_protein}

def degrade(macromolecule, **kwargs):
    '''Compartment-specific degradation reactions for proteins or protein-protein complexes'''
    
    if not (macromolecule.type == 'protein' or macromolecule.type == 'complex'):
        raise ValueError('Macromolecule to degrade must be protein or complex')
    if macromolecule.type == 'complex' and not macromolecule._deg_initialized:
        macromolecule._initialize_deg_params()
    
    deg_reactions = degrade_reaction_map[macromolecule.compartment](macromolecule, **kwargs)
    
    if macromolecule.compartment not in ['r', 'g']:
        dr = deg_reactions
    else:
        dr = deg_reactions[0]
    for r in dr:
        r._consolidate_macromolecules()    
    
    if len([r for r in deg_reactions if len(r.check_mass_balance()) > 0])>0:
        raise ValueError(macromolecule.id + ': Degradation reactions are unbalanced')
    return deg_reactions


# In[55]:


# import random
# import pandas as pd
# from expression.gene_information import gene_information
# import expression.build_mrna_expression_reactions as build_mrna
# from expression.protein_expression import ubiquitin
# from macromolecules.protein import Protein
# from macromolecules.complex import Complex

# proteins = []
# for i in range(2):
#     psim_toy = pd.DataFrame(columns = ['HGNC_ID', 'PREMRNA_SEQ', 'MRNA_SEQ', 'PROTEIN_SEQ', 'POLYA_LENGTH', 'TMD', 
#                                    'SP', 'N_INTRONS', 'DSB', 'GPI', 'OG', 'LOCATION'])

#     hgnc_id, premrna_seq = 'HGNC:TOY', ''.join(random.choices(['U', 'C', 'G', 'A'], k = 100))
#     mrna_seq = premrna_seq[25:75]
#     # note that there is no check that the protein_sequence corresponds to the mrna_sequence beyond checking for the length
#     protein_seq = ''.join(random.choices(params.amino_acids, k = int(len(mrna_seq)/3)))
#     polyA_length, tmd, sp, n_exons, dsb, gpi, og  = None, 1, True, None, 2, 2, 2
#     ub_args = ubiquitin.express_ubiquitin(compress_mrna = False)

#     import itertools
#     reactions = list()

#     location = list(params.compartments.keys())
#     psim_toy.loc[0,:] = [hgnc_id, premrna_seq, mrna_seq, protein_seq, polyA_length, tmd, sp, n_exons, dsb, gpi, og, location]
#     gene_info = gene_information(hgnc_id, premrna_seq, mrna_seq, protein_seq,
#                      ptms = {}, tmd = tmd, sp = sp, polyA_length = polyA_length, 
#                      n_exons = n_exons) 
#     gene_info.get_final_locations(metabolic_model = cobra.Model(''), final_locations = location)
#     proteins.append(Protein(id_ = 'a', compartment = 'c', gene_info = gene_info))
    
#     proteins_ = list()
#     for p in proteins:
#         proteins_.append(p.change_compartment('x'))
# #     proteins = proteins_

# cplx = Complex(metabolites = proteins_+proteins_, complex_id = 'test')
# # cplx._initialize_deg_params()


# In[56]:


# rxns = degrade(proteins_[0])#, **{'ub_args': ub_args})
# print([r.check_mass_balance() for r in rxns])
# rxns[0]


# In[57]:


# rxns = degrade(cplx)
# print([r.check_mass_balance() for r in rxns])
# rxns[0]

