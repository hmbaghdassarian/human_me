#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import cobra
# from collections import namedtuple

# import sys
# sys.path.insert(1, '../../scripts/')
# from macromolecules.complex import Complex
# from macromolecules.protein import Protein
# from macromolecules.RNA import mRNA


# In[174]:


# class Coupled_Constraint(mRNA, Protein, Complex, cobra.Metabolite):
#     def __init__(self, metabolite):
#         self.name = ''
#         if isinstance(metabolite, mRNA):
#             gi = namedtuple('gene_info', [ 'hgnc_id', 'mrna_seq'])
#             mRNA.__init__(self, gi(metabolite.id.split('_')[0], metabolite.sequence), 
#                           compartment = metabolite.compartment, triphosphate = metabolite.triphosphate)
#         elif isinstance(metabolite, Protein):
#             Protein.__init__(self, compartment = metabolite.compartment, 
#                              id_ = metabolite.id.split('_protein_')[0], 
#                              amino_acid_counts = {})
#             self.charge = metabolite.charge
#             self.elements = metabolite.elements
#         elif isinstance(metabolite, Complex):
#             Complex.__init__(self, metabolites = list(metabolite.components.keys()), complex_id = None)
#         else:
#             cobra.Metabolite.__init__(self, id=metabolite.id, formula=metabolite.formula, name=metabolite.name, 
#                                       charge=metabolite.charge, compartment=metabolite.compartment)
        
#         # in case other methods had been run on the metabolite
#         self.compartment = metabolite.compartment
#         self.charge = metabolite.charge
#         self.elements = metabolite.elements
#         self.name = metabolite.name
#         self.id = metabolite.id
        
#         for r in list(metabolite.reactions):
#             if r.model is not None:
#                 raise ValueError('Coupling constraint metabolites must be generated before creating model')
#             r.add_metabolites({self: r.metabolites[metabolite], metabolite: 0}, combine = False)

