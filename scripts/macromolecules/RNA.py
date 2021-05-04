#!/usr/bin/env python
# coding: utf-8

# In[2]:


import cobra
import sys
sys.path.insert(1, '../../scripts/')
from utils import metabolites as metab
from utils import machinery as mach
from utils import functions as func
from utils import parameters as params
# from uniform_processes.biomass import biomass_rna_mapper
from macromolecules.macromolecule import Macromolecule
from core.reaction import Expression_Reaction


# In[ ]:


molecule_type = ['premrna', 'rrna', 'trna']
synth_mach_map = dict(zip(molecule_type, [' and '.join(mach.ec), '', ' and '.join(mach.rnap3_transcription_machinery)]))

rule_part2 = ' and '.join(mach.lariat_machinery['Exosome'] + mach.lariat_machinery['NEXT Complex']) + ')'
rule = mach.lariat_machinery["5' Degradation"][0] + ' or (' + rule_part2
deg_mach_map = dict(zip(molecule_type, 
                       ['', {'n': rule, 'c': ' and '.join(mach.exosome['HGNC ID (gene)'].tolist())}, 
                        ' and '.join(mach.lariat_machinery['Exosome'])]))

molecule_type.append('mrna')
subsystem_map = dict(zip(molecule_type, ['mRNA_expression', 'rRNA_expression', 'tRNA_Biogenesis', 
                                    'mRNA_expression']))
rb_map = dict(zip(molecule_type, [False, True, False, False]))



# In[ ]:


class RNA(Macromolecule):
    def __init__(self, metabolite_name, seq, compartment = 'n', triphosphate = True, hgnc_id = None):
        '''
        
        Generates an RNA metabolite. Seq is a string of the mrna sequence, triphosphate indicates 
        the presence (True) or absence (False) of a triphosphate on the 5' end of the molecule. Otherwise, 
        assumes a monophosphate. 
        
        '''
        rna_id = metabolite_name + '_RNA_' + compartment 
        
        self.sequence = seq
        self.triphosphate = triphosphate
        self.length = len(self.sequence)
        self.get_base_counts_and_elements()

        Macromolecule.__init__(self, id = rna_id, compartment = compartment, charge = -self.length, elements = self.elements, 
                              hgnc_id = hgnc_id)
        if triphosphate:
            self.charge -= 3
    
    def get_base_counts_and_elements(self):
        '''
        Updates RNA metabolite to have appropriate formula according to sequence
        
        '''
        self.base_counts, self.elements = func.get_base_counts_and_elements(seq = self.sequence, 
                                                                  triphosphate = self.triphosphate)        

    def synthesize(self, id_):
        '''
        Generates a reaction for transcription of an RNA molecule (NTPs-->RNA).
        Inputs:
        1) self is an object of type RNA representing the rna molecule to be degraded.
        2) id_ is a string representing the name you want to give the reaction
        hgnc_id: str
            HGNC ID of the gene for which this reaction is being generated
        
        Output: a degradation reaction of type Expression_Reaction
        '''
        
        if self.type not in ['premrna', 'rrna', 'trna']:
            raise ValueError('Only premrna, rrna, or trna can be synthesized')  

#         hgnc_id = None
#         if self.type == 'premrna':
#             hgnc_id = self.id.split('_')[0]
        
        rna_synthesis = Expression_Reaction(id_, subsystem = subsystem_map[self.type], hgnc_id = self.hgnc_id,
                                           ribosome_biogenesis = rb_map[self.type])
        rxn = dict()
        for ntp, base_letter in metab.seq_metabolite_map.items():
            rxn[ntp] = -1*self.base_counts[base_letter]
        # pyrophosphate released per base added, -1 for 3/5' ends
        rxn[metab.ppi_n] = self.length - 1
        rxn[self] = 1
        rna_synthesis.add_metabolites(rxn)
        rna_synthesis.gene_reaction_rule = synth_mach_map[self.type]

        if len(rna_synthesis.check_mass_balance()) > 0:
            raise ValueError('RNA synthesis for ' + id_ + ' is unbalanced')
        elif list(rna_synthesis.compartments) != ['n']:
            raise ValueError('RNA synthesis must be confined to nuclear compartment')
        else:
            return rna_synthesis
            
        
    def exonucleolytic_degradation(self, reaction_name, balanced = True, update = False):
        ''' 

        Generates a reaction for exonucleolytic cleavage of an RNA molecule (RNA-->NMPs).
        Inputs:
        1) self is an object of type RNA representing the rna molecule to be degraded.
        2) reaction_name is a string representing the name you want to give the reaction
        hgnc_id: str
            HGNC ID of the gene for which this reaction is being generated
        3) balanced is a boolean indicating whether to check for mass balance
        4) update is a boolean indicating whether to update degradation reaction with subsystem and machinery

        Output: a degradation reaction of type Expression_Reaction
        no GPRs or subsystems added to reaction

        '''
        # exonucleolytic cleavage of RNA reaction
        if self.type in molecule_type:
            _type = self.type
        elif self.type == 'fragment_rna':
            if self.fragment_type == 'lariat':
                _type = 'premrna'
            elif self.fragment_type in ['5_leader', '3_trailer', 'trna_intron']:
                _type = 'trna'
            elif self.fragment_type in ['its', 'ets']:
                _type = 'rrna'
            else:
                raise ValueError('Only premrna, rrna, or trna can be degraded')  
        else:
            raise ValueError('Only premrna, rrna, or trna can be degraded')   
        
#         hgnc_id = None
#         if _type in ['premrna', 'mrna']:
#             hgnc_id = self.id.split('_')[0]
        rna_degradation = Expression_Reaction(reaction_name + '_DEGRADATION' + self.compartment, 
                                             subsystem = subsystem_map[_type], hgnc_id = self.hgnc_id, 
                                             ribosome_biogenesis = rb_map[_type])

        rxn = dict()
        rxn[metab.h2o_compartments[self.compartment]] = -self.length + 1# -sum(self.base_counts.values()) + 1
        rxn[self] = -1
        for k,v in metab.nmp_map[self.compartment].items():
            rxn[v] = self.base_counts[k]

        if self.triphosphate: # triphosphate on 5' end
            rxn[metab.nmp_map[self.compartment][self.sequence[0]]] -= 1
            rxn[metab.ntp_map[self.compartment][self.sequence[0]]]  = 1  
            rxn[metab.h_compartments[self.compartment]] = self.length - 1 #sum(self.base_counts.values())-1
        else:
            rxn[metab.h_compartments[self.compartment]] = self.length #sum(self.base_counts.values()) # extra H on 5' end <--unsure about this

        rna_degradation.add_metabolites(rxn)
        
        if update:
            if _type == 'rrna':
                rna_degradation.gene_reaction_rule = deg_mach_map[_type][self.compartment]
            elif _type == 'trna':
                rna_degradation.gene_reaction_rule = deg_mach_map[_type]
            else:
                raise ValueError('Only trna or rrna degradation reactions can be updated')
        

        if list(rna_degradation.compartments) != [self.compartment]:
            raise ValueError('RNA degradation has incorrect compartments')
            
        if balanced and len(rna_degradation.check_mass_balance())>0:
            raise ValueError('RNA degradation is not mass balanced')
        return rna_degradation 
    
    def update_metabolite(self, seq, append = False, append_to = None):
        '''Updates the RNA metabolite information according to an input sequence string. If add is True, 
        assumes string is being added, else, assumes string is being removed.'''
        
        if append:
            if append_to is None:
                raise ValueError("Must specify an end to append to append [5_primed, 3_primed]")
                
            if append_to == "5_primed":
                self.sequence = seq + self.sequence
            elif append_to == "3_primed":
                self.sequence += seq
            else:
                raise ValueError('The specified end to append to should be one of [5_primed, 3_primed]')
                
            self.length += len(seq)
            self.charge += -len(seq)
            
            base_counts = dict()
            for base_letter in metab.seq_element_map.keys():
                base_counts[base_letter] = seq.count(base_letter)
                self.base_counts[base_letter] += seq.count(base_letter)
            
            new_elements = self.elements.copy()
            for base_letter in metab.seq_element_map.keys():
                for element in new_elements.keys():
                    new_elements[element] += base_counts[base_letter]* metab.seq_element_map[base_letter][element]    
            self.elements = new_elements
            
        else:
            raise ValueError('Situation in which RNA sequence is removed or replaced has not been implemented yet')


# In[ ]:


class pre_mRNA(RNA):
    def __init__(self, gene_info, compartment = 'n', triphosphate = True):
        if compartment != 'n':
            raise ValueError("Premrna's outside of the nucleus are not currently considered")
        
        RNA.__init__(self, metabolite_name = gene_info.hgnc_id, seq = gene_info.premrna_seq, 
                     compartment = compartment, triphosphate = triphosphate, hgnc_id = gene_info.hgnc_id)
        self.type = 'premrna'
        self.id = self.id.replace('RNA', self.type)
        self.hgnc_id = gene_info.hgnc_id

class mRNA(RNA):
    def __init__(self, gene_info, compartment = 'n', triphosphate = True):
        if compartment not in ['n', 'c']:
            raise ValueError('mRNA must either be in nucleus or cytosol')
        
        
        RNA.__init__(self, metabolite_name = gene_info.hgnc_id, seq = gene_info.mrna_seq, 
                     compartment = compartment, triphosphate = triphosphate, hgnc_id = gene_info.hgnc_id)
        self.type = 'mrna'
        self.id = self.id.replace('RNA', self.type)
        self.hgnc_id = gene_info.hgnc_id
        

class tRNA(RNA):
    def __init__(self, metabolite_name, seq, compartment = 'n', triphosphate = True):
        RNA.__init__(self, metabolite_name = metabolite_name, seq = seq, compartment = compartment, 
             triphosphate = triphosphate)
        
        self.type = 'trna'
        self.id = self.id.replace('RNA', self.type)
        
    
class rRNA(RNA):
    def __init__(self, metabolite_name, seq, compartment = 'n', triphosphate = True):
        RNA.__init__(self, metabolite_name = metabolite_name, seq = seq, compartment = compartment, 
             triphosphate = triphosphate)
        
        self.type = 'rrna'
        self.id = self.id.replace('RNA', self.type)
        self.k_deg = params.rrna_degradation_constant
        
        
class RNA_fragment(RNA):
    def __init__(self, metabolite_name, seq, fragment_type, compartment = 'n', triphosphate = True, 
                hgnc_id = None):
        
        if fragment_type not in ['lariat', 'its', 'ets', '5_leader', '3_trailer', 'trna_intron']:
            raise ValueError('RNA fragment type specified is not considered')
        
        if compartment not in ['n', 'c']:
            raise ValueError('Only nuclear or cytosolic RNA fragments are incorporated for now')
        
        RNA.__init__(self, metabolite_name = metabolite_name, seq = seq, compartment = compartment, 
                     triphosphate = triphosphate, hgnc_id = hgnc_id)
        
        self.type = 'fragment_rna'
        self.fragment_type = fragment_type
        if self.fragment_type == 'lariat' and hgnc_id is None:
            raise ValueError('Must specify hgnc ID for lariats')
        self.id = self.id.replace('RNA', self.fragment_type)

