#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import cobra

import sys
sys.path.insert(1, '../../scripts/')
from utils import parameters as params


# In[ ]:


class Macromolecule(cobra.Metabolite):
    def __init__(self, id=None, formula=None, name="",charge=None, compartment=None, elements = None,
                 hgnc_id = None):
        '''Inherits from cobra.Metabolite. See help(cobra.Metabolite) for additional parameters
        
        Parameters
        ----------
        proxy: bool
            whether the object is a proxy metabolite (for coupling purposes, no associated mass/charge)
        hgnc_id: str
            the associated gene HGNC ID of the macromolecule (HGNC:####)
        
        '''
        
        if compartment not in params.compartments.keys():
            err = 'Specified compartment is not considered in the ME Model. Please input one of the following: ' 
            err += ', '.join(list(params.compartments.keys()))
            raise ValueError(err)
        
        
        cobra.Metabolite.__init__(self, id = id, charge = charge, compartment = compartment)
        
        if elements is not None:
            self.elements = elements
        if self.id.split('_')[-1] != self.compartment:
            raise ValueError('Macromolecules must syntactically have compartment as part of id')

        self.coupling_coefficient = None
        self.hgnc_id = hgnc_id
    
    def change_compartment(self, new_compartment):
        '''Returns a copy of the macromolecule metabolite, but in new compartment'''
        
        if new_compartment == self.compartment:
            raise ValueError('The macromolecule is already in this compartment')
        if new_compartment not in params.compartments.keys():
            err = 'Specified compartment is not considered in the ME Model. Please input one of the following: ' 
            err += ', '.join(list(params.compartments.keys()))
        
        new_macromolecule = self.copy()
        new_macromolecule.id = '_'.join(self.id.split('_')[:-1]) + '_' + new_compartment
        new_macromolecule.compartment = new_compartment
        
        return new_macromolecule    
    
    def couple(self, type, value):
        '''
        Parameters
        ----------
            type: str
                the type of reaction this macromolecule is coupled to
            value: float or sympy.Expr (function of parameters.mu)
                the coupling coefficient value 
        
        Returns
        ----------
        self.coupling_coefficient: dictionary 
            dictionary of length one, key is the type, value is the coupling coefficient 
        
        '''
        
        if type not in ['catalysis', 'enzyme_degradation', 'mrna_degradation', 'mrna_formation']:
            raise ValueError('The couple id must be one of catalysis, mrna_degradation, enzyme_degradation, or mrna_formation')
        
        if self.coupling_coefficient is None:
            self.coupling_coefficient = {type: value}
        else:
            if [type] != list(self.coupling_coefficient):
                raise ValueError('More than one coupling type associated with macromolecule: ' + self.id)
            if value != self.coupling_coefficient[type]:
                raise ValueError('More than one coupling coefficient value associated with macromolecule: ' + self.id)
        
        if type == 'catalysis':
            self.enzyme = True


# In[ ]:


class Proxy(Macromolecule):
    '''For c2/c4 coupling of degradation'''
    def __init__(self, associated_macromolecule):
        '''
        Parameters
        ----------
        associated_macromolecule: Macromolecule
            the c1/c3 associated macromolecule to the respective c2/c4 coupling
        '''
        if associated_macromolecule.type not in ['mrna', 'protein', 'complex']:
            raise ValueError('Unexpected associated macromolecule for proxy metabolite')
        key_mapper = {'mrna': 'mrna_degradation', 'protein': 'enzyme_degradation', 
                     'complex': 'enzyme_degradation'}
        
        id_ = associated_macromolecule.hgnc_id if associated_macromolecule.type == 'mrna' else associated_macromolecule._deg_id 
        Macromolecule.__init__(self, id = '_'.join([id_, 
                                                    key_mapper[associated_macromolecule.type], 'proxy', 
                                                   associated_macromolecule.compartment]), 
                               compartment = associated_macromolecule.compartment, 
                               hgnc_id = associated_macromolecule.hgnc_id)
        self.type = 'proxy'
        self.associated_macromolecule = associated_macromolecule.id
        self._amt = associated_macromolecule.type
        # only for complexes, used in Expressed_Gene class
        self._complex_hgnc_ids = [p.hgnc_id for p in associated_macromolecule.decompose_complex()                                 if p.hgnc_id is not None] if self._amt == 'complex' else [] 
        
    
    def couple(self, value):
        
        '''
        Parameters
        ----------
            type: str
                the type of reaction this macromolecule is coupled to
            value: float or sympy.Expr (function of parameters.mu)
                the coupling coefficient value 
        
        Returns
        ----------
        self.coupling_coefficient: dictionary 
            dictionary of length one, key is the type, value is the coupling coefficient 
        
        '''

        key_mapper = {'mrna': 'mrna_degradation', 'protein': 'enzyme_degradation', 
                     'complex': 'enzyme_degradation'}
        
        if self.coupling_coefficient is None:
            self.coupling_coefficient = {key_mapper[self._amt]: value}
        else:
            if value != self.coupling_coefficient[key_mapper[self._amt]]:
                raise ValueError('More than one coupling coefficient value associated with macromolecule: ' + self.id)
        
            
        

