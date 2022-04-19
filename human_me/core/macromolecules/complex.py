#!/usr/bin/env python
# coding: utf-8

import collections
from collections import OrderedDict
import copy
from typing import Dict, List, Optional

import numpy as np
from faker import Faker

from human_me.core.macromolecules.macromolecule import Macromolecule, Proxy
from human_me.core.reaction import ExpressionReaction
from human_me.utils import machinery as mach
from human_me.utils import parameters as params
from human_me.utils.functions import flatten_list

cotransloc_ids = set([mid + '_folded_protein_r' for mid in mach.ctnm + mach.translation_efs])


class Complex(Macromolecule):
    """Complexes formed by non-covalent interactions between macromolecules."""

    type = 'complex' 

    def __init__(self, metabolites: List[Macromolecule], complex_id: Optional[str] = None, seed: Optional[int] = None):
        """Init method for Complex.

        Parameters
        ----------
        metabolites : List[Macromolecule]
            the macromolecules to form a complex
        complex_id : str, optional
            complex identifier, by default None
        seed : Optional[int], optional
            A seed for generating the complex ID if complex_id is None, by default None
        """
        # checks
        if type(metabolites) != list or len(metabolites) == 0:
            raise ValueError('Must provide a list of macromolecules to form complex')
        # cobra metabolite not set up, check for bc they don't have the attribute .type
        if len([m for m in metabolites if not isinstance(m, Macromolecule)]) > 0:
            raise ValueError('Generic cobra.Metabolite cannot form complexes with macromolecules currently')

        self.components = {m: metabolites.count(m) for m in metabolites}
        # parse compartment
        compartments = {m.compartment for m in self.components}
        if len(compartments) == 1:
            compartment = list(compartments)[0]
        else:
            raise ValueError('Metabolites forming a complex must all be in the same compartment')

        # parse metabolite id
        self._seed = seed
        if complex_id is None:
            Faker.seed(self._seed)
            f1 = Faker()
            self.temp_id = f1.uuid4().split('-')[0]
        else:
            self.temp_id = complex_id

        elements = dict()
        for m, count in self.components.items():
            for k, v in m.elements.items():
                if k in elements:
                    elements[k] += v * count
                else:
                    elements[k] = v * count

        # make the metabolite
        super().__init__(id=self.temp_id + '_complex_' + compartment, compartment=compartment,
                               charge=sum([m.charge * count for m, count in self.components.items()]),
                               elements=elements)

        self.reaction_id = None  # none before running form_complex(); this is used in update_id() method
        self._deg_initialized = False
        self.enzyme = False
        self.keff = None
        self.hgnc_id = None  # always None, for internal use with expression/protein_expression/degradation

    def update_id(self, new_id: Optional[str] = None):
        """In cases where complex id is too long (see build_me_model generate_complex_reactions method)."""
        if self.reaction_id is not None:
            raise ValueError(
                'Reaction and complex IDs will be consistent since you are updating the id after forming the reaction.')
        if new_id is None:
            Faker.seed(self._seed)
            f1 = Faker()
            self.temp_id = f1.uuid4().split('-')[0]
        else:
            self.temp_id = new_id
        self.id = self.temp_id + '_complex_' + self.compartment

    def form_complex(self, reaction_id: Optional[str] = None, reversible: bool = False,
                     synthesis: bool = True, synthesis_type: str = 'complex') -> ExpressionReaction:
        """The reaction to generate the Complex object.
        
        Note: assumes non-covalent complex formation (in terms of elemental balance)

        Parameters
        ----------
        reaction_id : Optional[str], optional
            reaction identifier, by default None
        reversible : bool, optional
            Whether the reaction is reversible, by default False
            Setting to True may make model more efficient (allows reuse of self.components if involved in other reactions)
        synthesis : bool, optional
            whether the reaction represents the "main" synthesis/production for the macromolecule (used for appropriate mapping of coupling)
            intended for use with genes (reactions with an associated hgnc id, and complexes), by default True
        synthesis_type : str, optional
            [description], by default 'complex'

        Returns
        -------
        complex_formation : ExpressionReaction
            the complex formation between metabolites stored in self.components
        """
        self._check_metabolite_types()
        if reaction_id is None:
            self.reaction_id = self.temp_id + '_COMPLEX_FORMATION' + self.compartment
        else:
            self.reaction_id = reaction_id + '_COMPLEX_FORMATION' + self.compartment

        complex_formation = ExpressionReaction(self.reaction_id, subsystem='Complex_Formation',
                                               synthesis=synthesis, synthesis_type=synthesis_type)
        rxn = {m: -count for m, count in self.components.items()}
        rxn[self] = 1
        complex_formation.add_metabolites(rxn)
        if reversible:
            complex_formation.lower_bound = -1000  # reversible

        return complex_formation

    def decompose_complex(self, decomposed_complex=None) -> OrderedDict:
        """Recursive method to get the complex by its individual components, including nested complexes."""
        # TODO: this method is very slow
        if decomposed_complex is None:
            all_metab = flatten_list([[m] * count for m, count in self.components.items()])
            decomposed_complex = Complex(metabolites=all_metab, complex_id='ignore')

        if 'complex' not in [m.type for m in decomposed_complex.components.keys()]:
            #             return decomposed_complex.components
            # order the metabolites to avoid precision issues, e.g. in .get_complex_biomass
            dc = decomposed_complex.components
            dc_map = {m.id: m for m in dc}
            return OrderedDict({dc_map[m_id]: dc[dc_map[m_id]] for m_id in sorted(dc_map)})

        metabolites_ = flatten_list(
            [[m] * count for m, count in decomposed_complex.components.items() if m.type != 'complex'])
        metabolites_ += flatten_list(flatten_list(
            [[[m_] * count_ for m_, count_ in m.components.items()] * count for m, count in
                decomposed_complex.components.items() if m.type == 'complex']))
        return self.decompose_complex(decomposed_complex=Complex(metabolites=metabolites_, complex_id='ignore'))

    def get_complex_biomass(self) -> Dict[str, float]:
        """Returns a dictionary of the complex biomass by its individual component types."""
        biomass_by_type = dict()
        for m, count in self.decompose_complex().items():
            if m.type in biomass_by_type:
                biomass_by_type[m.type] += count * (m.formula_weight / 1000)
            else:
                biomass_by_type[m.type] = count * (m.formula_weight / 1000)

        return biomass_by_type

    def _initialize_deg_params(self):
        """Initialize attributes for creating degradation reactions."""
        self._deg_initialized = True
        dc = self.decompose_complex()
        self._check_metabolite_types()

        self._amino_acid_counts = collections.Counter()
        self._ptms = collections.Counter()
        self.length = 0
        for p, c in dc.items():
            self.length += p.length * c
            for i in range(c):
                self._amino_acid_counts.update(p._amino_acid_counts) # TODO: unused var i, check
                if hasattr(p, '_ptms'):
                    self._ptms.update(p._ptms)

        if len(self._ptms) > 0:
            raise ValueError(
                'PTMs to Complexes is currently unaccounted for and will likely lead to imbalances in degradation reactions')

        self._deg_id = self.temp_id + '_COMPLEX'
        self._degradation_reactions = set()
        del dc

    def change_compartment(self, new_compartment: str):
        """Returns a copy of the complex metabolite, but in new compartment."""
        self._check_metabolite_types()
        return self._change_compartment_and_components(new_compartment)

    def _change_compartment_and_components(self, new_compartment):
        """Create a complex the same as self, but in a different compartment.

        Parameters
        ----------
        new_compartment : str
            one-letter code of compartment to change to 

        Returns
        -------
        Complex
            a copy of the macromolecule in the new compartment
        """
        if new_compartment == self.compartment:
            raise ValueError('The macromolecule is already in this compartment')
        if new_compartment not in params.compartments:
            err = 'Specified compartment is not considered in the ME Model. Please input one of the following: '
            err += ', '.join(list(params.compartments.keys()))

        metabolites = list()
        for m, c in self.components.items():
            if m.type != 'complex':
                metabolites += [m.change_compartment(new_compartment)] * c
            else:
                metabolites += [m._change_compartment_and_components(new_compartment)] * c

        new_complex = Complex(metabolites=metabolites, complex_id=self.temp_id)
        if self._deg_initialized:
            new_complex._initialize_deg_params()
        return new_complex

    def _check_metabolite_types(self):
        """Only protein-protein complexes can be formed for now."""
        types = list(set([m.type for m in self.decompose_complex()]))
        if len(types) > 1 or types[0] != 'protein':
            raise ValueError('ME-Model is currently designed to only handle protein-protein complexes')

    def get_k_deg(self):
        self.k_deg = np.median([m.k_deg * c for m, c in self.decompose_complex().items()])

    def make_proxy(self):
        """Make a proxy metabolite for coupling enzyme degradation to reaction catalysis."""
        return Proxy(associated_macromolecule=self)
    def copy(self, copy_metabolites = False):
        """Overwrite cobra.Species.copy.

        Parameters
        ----------
        copy_metabolites : bool, optional
            whether to copy the individual complex subunits, or retain the same pointer, by default False
            set to True if will be modifying the individual complex subunits in just this complex, but not others

        Returns
        -------
        _type_
            _description_
        """
        # copy while preserving metabolite pointers to avoid memory issues and not copying it to avoid speed issues
        components = self.components
        del self.components

        new_macromolecule = copy.deepcopy(self)
        
        self.components = components
        if not copy_metabolites:
            new_macromolecule.components = {subunit: stoich for subunit, stoich in components.items()} # retain same pointer to avoid memory issues
        else:
            new_macromolecule.components = {subunit.copy(): stoich for subunit, stoich in components.items()} # new pointer, but still faster 
        return new_macromolecule


class RibosomalComplex(Complex):
    """Complexes specifically associated with ribosome biogenesis, which has RNA-protein complexes and
    multiple compartments"""

    type = 'complex' 

    def __init__(self, metabolites: List[Macromolecule], complex_id: Optional[str] = None, ignore_compartment: bool = False, seed: Optional[int] = None):
        """Init method for RibosomalComplex.

        Parameters
        ----------
        metabolites : List[Macromolecule]
            each entry is a Macromolecule object (protein or RNA or complex, not generic metabolites)
        complex_id : Optional[str], optional
            the id of the complex metabolite; if None, will generate a random id, by default None
        ignore_compartment : bool, optional
            whether to ignore the metabolite compartments, mainly for internal use, by default False
        seed : Optional[int], optional
            A seed for generating the complex ID if it is None, by default None
        """
        # checks
        if type(metabolites) != list or len(metabolites) == 0:
            raise ValueError('Must provide a list of macromolecules to form complex')
        # cobra metabolite not set up, check for bc they don't have the attribute .type
        if len([m for m in metabolites if not isinstance(m, Macromolecule)]) > 0:
            raise ValueError('Generic cobra.Metabolite cannot form complexes with macromolecules currently')

        self.components = {m: metabolites.count(m) for m in metabolites}
        # parse compartment
        compartments = list(set([m.compartment for m in self.components]))

        # test compartment consistency - exception of ribosome complexes
        comp_ids = [m.id for m in self.components]
        cotransloc_cond = (len(cotransloc_ids.difference(comp_ids)) == 0) or ('mature_ribosome_complex_c' in comp_ids)

        if len(compartments) == 1:
            compartment = compartments[0]
        elif (sorted(compartments) == ['c', 'r']) and cotransloc_cond:
            compartment = 'c'
        else:
            raise ValueError('Metabolites forming a complex must all be in the same compartment')

        # parse metabolite id
        self._seed = seed
        if complex_id is None:
            Faker.seed(self._seed)
            f1 = Faker()
            self.temp_id = f1.uuid4().split('-')[0]
        else:
            self.temp_id = complex_id

        elements = dict()
        for m, count in self.components.items():
            for k, v in m.elements.items():
                if k in elements:
                    elements[k] += v * count
                else:
                    elements[k] = v * count

        # make the metabolite
        Macromolecule.__init__(self, id=self.temp_id + '_complex_' + compartment, compartment=compartment,
                               charge=sum([m.charge * count for m, count in self.components.items()]),
                               elements=elements)

        self.reaction_id = None  # none before running form_complex(); this is used in update_id() method
        self._deg_initialized = False
        self.enzyme = False
        self.keff = None
        self.hgnc_id = None  # always None, for internal use with expression/protein_expression/degradation

    def decompose_complex(self, decomposed_complex=None):
        """Recursive method to get the complex by its individual components, including nested complexes."""
        if decomposed_complex is None:
            all_metab = flatten_list([[m] * count for m, count in self.components.items()])
            decomposed_complex = RibosomalComplex(metabolites=all_metab, complex_id='ignore')

        if 'complex' not in [m.type for m in decomposed_complex.components.keys()]:
            #             return decomposed_complex.components
            # order the metabolites to avoid precision issues, e.g. in .get_complex_biomass
            dc = decomposed_complex.components
            dc_map = {m.id: m for m in dc}
            return OrderedDict({dc_map[m_id]: dc[dc_map[m_id]] for m_id in sorted(dc_map)})

        metabolites_ = flatten_list(
            [[m] * count for m, count in decomposed_complex.components.items() if m.type != 'complex'])
        metabolites_ += flatten_list(flatten_list(
            [[[m_] * count_ for m_, count_ in m.components.items()] * count for m, count in
                decomposed_complex.components.items() if m.type == 'complex']))
        return self.decompose_complex(decomposed_complex=RibosomalComplex(metabolites=metabolites_, complex_id='ignore'))

    def form_complex(self, reaction_id: Optional[str] = None, reversible: bool = False,
                     synthesis: bool = False, synthesis_type: Optional[str] = None) -> ExpressionReaction:
        """The reaction to generate the Complex object.
        
        Note: assumes non-covalent complex formation (in terms of elemental balance)

        Parameters
        ----------
        reaction_id : Optional[str], optional
            reaction identifier, by default None
        reversible : bool, optional
            Whether the reaction is reversible, by default False
            Setting to True may make model more efficient (allows reuse of self.components if involved in other reactions)
        synthesis : bool, optional
            whether the reaction represents the "main" synthesis/production for the macromolecule (used for appropriate mapping of coupling)
            intended for use with genes (reactions with an associated hgnc id, and complexes), by default False
        synthesis_type : str, optional
            [description], by default None

        Returns
        -------
        complex_formation : ExpressionReaction
            the complex formation between metabolites stored in self.components
        """
        self._check_metabolite_types()
        if reaction_id is None:
            self.reaction_id = self.temp_id + '_COMPLEX_FORMATION' + self.compartment
        else:
            self.reaction_id = reaction_id + '_COMPLEX_FORMATION' + self.compartment

        complex_formation = ExpressionReaction(self.reaction_id, subsystem='Complex_Formation',
                                               synthesis=synthesis, synthesis_type=synthesis_type,
                                               ribosome_biogenesis=True)
        rxn = {m: -count for m, count in self.components.items()}
        rxn[self] = 1
        complex_formation.add_metabolites(rxn)
        if reversible:
            complex_formation.lower_bound = -1000  # reversible

        return complex_formation

    def _check_metabolite_types(self):
        """Only protein and rRNA can be included in ribosomal complexes."""
        if len((set([m.type for m in self.decompose_complex()])).difference(['protein', 'rrna'])) > 0:
            raise ValueError('Ribosomal complexes can only include proteins and rRNA')

    def _initialize_deg_params(self):
        """Initialize attributes for creating degradation reactions."""

        self._deg_initialized = True
        self.keff = None
        dc = self.decompose_complex()
        self._check_metabolite_types()

        self._amino_acid_counts = collections.Counter()
        self._ptms = collections.Counter()
        self.length = {'protein': 0, 'rrna': 0}
        for p, c in dc.items():
            self.length[p.type] += p.length * c
            if p.type == 'protein':
                for i in range(c):
                    self._amino_acid_counts.update(p._amino_acid_counts)
                    if hasattr(p, '_ptms'):
                        self._ptms.update(p._ptms)

        if len(self._ptms) > 0:
            raise ValueError(
                'PTMs to Complexes is currently unaccounted for and will likely lead to imbalances in degradation reactions')

        self._deg_id = self.temp_id + '_COMPLEX'
        self._degradation_reactions = set()

        del dc

    def _change_compartment_and_components(self, new_compartment: str):
        """Returns a copy of the complex metabolite, but in new compartment. Recursive to change all components (nested complexes and their components)."""
        if new_compartment == self.compartment:
            raise ValueError('The macromolecule is already in this compartment')
        if new_compartment not in params.compartments:
            err = 'Specified compartment is not considered in the ME Model. Please input one of the following: '
            err += ', '.join(list(params.compartments.keys()))

        metabolites = list()
        for m, c in self.components.items():
            if m.type != 'complex':
                metabolites += [m.change_compartment(new_compartment)] * c
            else:
                metabolites += [m._change_compartment_and_components(new_compartment)] * c

        new_complex = RibosomalComplex(metabolites=metabolites, complex_id=self.temp_id)
        if self._deg_initialized:
            new_complex._initialize_deg_params()
        return new_complex


def add_complex_metabolites(cplx: Complex, met_to_add: Dict[Macromolecule, int], complex_id: str) -> Complex:
    """Add a metabolite to an existing complex object, returning it as a new complex.

    Parameters
    ----------
    cplx : Complex
        The Complex object to be appended
    met_to_add : Dict[Macromolecule, int]
        The metabolites to add to the complex. Keys are the macromolecules and values are the number of copies to add.
    complex_id : str
        The new id for the new complex

    Returns
    -------
    cplx2 : Complex
        The new Complex object with the added metabolites
    """
    mtblts = list()
    for m, c in cplx.components.items():
        mtblts += [m] * c
    for m, c in met_to_add.items():
        mtblts += [m] * c

    if not isinstance(cplx, RibosomalComplex):
        cplx2 = Complex(metabolites=mtblts,
                        complex_id=complex_id)
    else:
        cplx2 = RibosomalComplex(metabolites=mtblts,
                                 complex_id=complex_id)

    if cplx._deg_initialized:
        cplx2._initialize_deg_params()
        cplx2._degradation_reactions = cplx2._degradation_reactions.union(cplx._degradation_reactions)  # inherit degradation reactions
    return cplx2
