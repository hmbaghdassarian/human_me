#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import cobra


# In[ ]:


class Complex(cobra.Metabolite):
    def __init__(self, complex_id = None, reaction_id = None, **complex_info):
        '''
        Inputs:
        Complex info is a dictionary with three keys ['METABOLITES', 'IDS', 'METABOLITE_TYPES']
        Each value is a list:
            metabolites is a list of cobra.Metabolite objects
            IDs is a list of string identifiers corresponding to each metabolite object
            component_types is a list of strings; possible values are ['protein', 'rrna', 'trna', 'mrna', 'metabolite']
            This means complexes can form between any of these species, including other complexes; metabolite is a M-model metabolite
        complex_id is a string for the id of the complex metabolite, otherwise will form one from metabolite ids

        Output:
        A cobra.Metabolite object representing the complex formed between metabolites stored in self.complex_metabolite
        
        '''
        # checks
        if sorted(set(complex_info.keys())) != ['IDS', 'METABOLITES', 'METABOLITE_TYPES']:
            raise ValueError('Invalid complex information keys or insufficient complex information keys')

        self.subcomponents, ids, self.component_types = complex_info['METABOLITES'], complex_info['IDS'], complex_info['METABOLITE_TYPES']
        
        if len(set(self.component_types).difference(['protein', 'rrna', 'trna', 'mrna',  'metabolite', 'complex']))>0:
            raise ValueError('At least one of the metabolite types is not considered in complex formation currently')
        
        if len(self.subcomponents) != len(self.component_types):
            raise ValueError('Each metabolite must have a corresponding metabolite type')
        # parse compartment    
        compartments = list(set([m.compartment for m in self.subcomponents]))
        if len(compartments) == 1:
            compartment = compartments[0]
        # exception of ribosome complex
        elif (len(compartments) == 2) and ('c' in compartments) and ('mature_ribosome_complex_complex[c]' in [m.id for m in self.subcomponents]):
            compartment = 'c'
        else:
            raise ValueError('metabolites are not in the same compartment')

        # parse metabolite id
        mt_type = '_'.join(list(set(self.component_types)))
        ids_ = '_'.join(ids)
        if complex_id == None:
            id_ = ids_ + '_' + mt_type
        else: 
            id_ = complex_id + '_' + mt_type
            
        if reaction_id == None:
            self.reaction_id = id_ + '_COMPLEX_FORMATION' + compartment
        else:
            self.reaction_id = reaction_id + '_COMPLEX_FORMATION' + compartment
            
        
        complex_id = id_ + '_complex' + '[' + compartment + ']'
        if len(complex_id)>(256-8-4-len(mt_type)): #-8 and -4 for _complex and compartment appended to end
            err_msg = 'Cobrapy requires metabolite ids to be less than 256 characters, please specify a '
            err_msg += 'shorter user-defined complex id'
            raise ValueError(err_msg)
        
        # make the metabolite
        cobra.Metabolite.__init__(self, id = complex_id, 
                                  compartment = compartment, charge = sum([m.charge for m in self.subcomponents]))
        elements = dict()
        for m in self.subcomponents:
            for k,v in m.elements.items():
                if k in elements.keys():
                    elements[k] += v
                else:
                    elements[k] = v
        self.elements = elements
        
    def form_complex(self):

        '''
        Output: A cobra.Reaction object representing the complex formation between metabolites stored in self.complex_formation

        '''
        
        # expected no biomass change
        complex_formation = cobra.Reaction(self.reaction_id)
        rxn = {m: -1 for m in self.subcomponents}
        rxn[self] = 1
        complex_formation.add_metabolites(rxn)
        complex_formation.lower_bound = -1000 # reversible
        
        return complex_formation
#         self.complex_formation = complex_formation
        
        
    def get_complex_biomass(self, decomposed_complex = None):
        '''Recursive method to get the complex biomass by its individual components'''
    
        if decomposed_complex == None:
            decomposed_complex = Complex(complex_id = 'foo', reaction_id = 'foo',
                                         **{'METABOLITES': self.subcomponents, 
                                           'IDS': [str(i) for i in list(range(len(self.subcomponents)))],
                                           'METABOLITE_TYPES': self.component_types})
            
            
        if 'complex' not in decomposed_complex.component_types:
            total_components = len(decomposed_complex.subcomponents)

            biomass_by_type = dict(zip(decomposed_complex.component_types, [0]* total_components))
            for i in range(total_components): # in case repeated metabolites in complex (homodimers)
                m_ = decomposed_complex.subcomponents[i]
                mt_ = decomposed_complex.component_types[i]
                biomass_by_type[mt_] += m_.formula_weight/1000 
            
#             if sum(biomass_by_type.values()) != get_metabolite_mw(self): # sanity check
#                 raise ValueError('Something went wrong in decomposing complex biomass')
            
            return biomass_by_type
        else: # unpack nested complexes
            complex_idx = [i for i in range(len(decomposed_complex.component_types)) if decomposed_complex.component_types[i] == 'complex']
            non_complex_idx = sorted(set(list(range(len(decomposed_complex.component_types)))).difference(complex_idx))


            metabolites_ = [decomposed_complex.subcomponents[i] for i in non_complex_idx] + [item for sublist in [decomposed_complex.subcomponents[i].subcomponents for i in complex_idx] for item in sublist]
            ids = [str(i) for i in list(range(len(metabolites_)))]
            metabolite_types = [decomposed_complex.component_types[i] for i in non_complex_idx] + [item for sublist in [decomposed_complex.subcomponents[i].component_types for i in complex_idx] for item in sublist]
            new_complex_info = {'METABOLITES': metabolites_, 'IDS': ids, 'METABOLITE_TYPES': metabolite_types}

            return self.get_complex_biomass(decomposed_complex = Complex(complex_id = 'foo', reaction_id = 'foo', 
                                                                         **new_complex_info))  


# In[ ]:


def get_complex_biomass_change(complex_products, complex_reactants):
    '''Input is two lists of type COMPLEX, one representing those on the product side, one representing those on the reactant side
    output is a dictionary of biomass change for each respective biomass type.'''
    
    product_biomass = dict()
    for cp in complex_products:
        if type(cp)!= Complex:
            raise TypeError('All complex products must be a COMPLEX object')
        for bt, mw in cp.get_complex_biomass().items():
            if bt in product_biomass.keys():
                product_biomass[bt] += mw
            else:
                product_biomass[bt] = mw
    
    reactant_biomass = dict()
    for cr in complex_reactants:
        if type(cr)!= Complex:
            raise TypeError('All complex reactants must be a COMPLEX object')
        for bt, mw in cr.get_complex_biomass().items():
            if bt in reactant_biomass.keys():
                reactant_biomass[bt] += mw
            else:
                reactant_biomass[bt] = mw
    
    for bt in set(product_biomass.keys()).difference(reactant_biomass.keys()):
        reactant_biomass[bt] = 0
    for bt in set(reactant_biomass.keys()).difference(product_biomass.keys()):
        product_biomass[bt] = 0    
    
    return {bt: product_biomass[bt] - reactant_biomass[bt] for bt in product_biomass.keys() if product_biomass[bt] - reactant_biomass[bt] != 0}

