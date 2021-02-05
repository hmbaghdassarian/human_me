#!/usr/bin/env python
# coding: utf-8

# In[7]:


import cobra

import sys
sys.path.insert(1, '../../scripts/')
from utils import parameters as params


# In[6]:


class Macromolecule(cobra.Metabolite):
    def __init__(self, id=None, formula=None, name="",charge=None, compartment=None, elements = None, 
                proxy = False):
        
        if not proxy and compartment not in params.compartments.keys():
            err = 'Specified compartment is not considered in the ME Model. Please input one of the following: ' 
            err += ', '.join(list(params.compartments.keys()))
            raise ValueError(err)
        
        
        cobra.Metabolite.__init__(self, id = id, charge = charge, compartment = compartment)
        if not proxy:
            self.elements = elements
            if self.id.split('_')[-1] != self.compartment:
                raise ValueError('Macromolecules must syntactically have compartment as part of id')
        else:
            self.type = 'proxy'
        self.coupling_coefficient = None
    
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
            raise ValueError('The couple id must be one of catalysis, mrna_degradation, or mrna_formation')
        else:
            if self.coupling_coefficient is None:
                self.coupling_coefficient = {type: value}
            else:
                self.coupling_coefficient[type] = value
        if type == 'catalysis':
            self.enzyme = True

