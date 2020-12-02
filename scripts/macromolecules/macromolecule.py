#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import cobra

import sys
sys.path.insert(1, '../../scripts/')
from utils import parameters as params


# In[ ]:


class Macromolecule(cobra.Metabolite):
    def __init__(self, id=None, formula=None, name="",charge=None, compartment=None, elements = None):
        
        if compartment not in params.compartments.keys():
            err = 'Specified compartment is not considered in the ME Model. Please input one of the following: ' 
            err += ', '.join(list(params.compartments.keys()))
            raise ValueError(err)
        
        
        cobra.Metabolite.__init__(self, id = id, charge = charge, compartment = compartment)
        self.elements = elements
        
        if self.id.split('_')[-1] != self.compartment:
            raise ValueError('Macromolecules must syntactically have compartment as part of id')
    
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

