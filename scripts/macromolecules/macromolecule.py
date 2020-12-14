#!/usr/bin/env python
# coding: utf-8

# In[5]:


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
    
    def couple(self, type, value, combine = False):
        '''
        Input:
            type: a string, one of ['catalysis', 'mrna_degradation', 'mrna_formation']
            value: The coupling coefficient value (sympy.Expr or float)
            combine: boolean, whether to 
        
        Returns: self.coupling_coefficient, dictionary of possible coupling coefficients and their values for the macromolecule
        
        '''
        
        if type not in ['catalysis', 'mrna_degradation', 'mrna_formation']:
            raise ValueError('The couple id must be one of catalysis, mrna_degradation, or mrna_formation')
        else:
            if self.coupling_coefficient is None:
                self.coupling_coefficient = {type: value}
            else:
                if not combine or id not in self.coupling_coefficient.keys():
                    self.coupling_coefficient[type] = value
                else:
                    self.coupling_coefficient[type] += value
                    
                    

