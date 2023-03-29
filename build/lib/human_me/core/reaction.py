#!/usr/bin/env python
# coding: utf-8
import copy
from collections import defaultdict
from collections.abc import Iterable
from operator import attrgetter
from typing import Dict, List, Optional, SupportsFloat, Union

import cobra
import numpy as np
import sympy
from six import iteritems

from human_me.utils import machinery as mach
from human_me.utils import parameters as params


class ME_Reaction(cobra.Reaction):
    """Allows stochiometric coefficient to be a function of mu."""

    def __init__(self, id, name='', subsystem='', lower_bound=0.0, upper_bound=None):
        """Helps distinguish between these reactions, which have mu in bounds, and coupling reactions, which
        have mu in stochiometric coefficient."""

        super().__init__(id, name, subsystem, lower_bound, upper_bound)
        self.coupled_metabolites = {}
        self._protein_deg_proxy = False
    
    def copy(self):
        """Overwrite cobra.Species copy method.

        Returns
        -------
        new_rxn : ME_Reaction
            a copy of the reaction
        """

        # avoid deeopcopying some attributes for:
        # 1) compute time
        # 2) maintaining pointers (same metabolites in the new reaction)

        # additionally, make sure metabolites add the new copied reaction to their attributes (mutual tracking)
        attr_keep = ['_metabolites', 'coupled_metabolites'] # in the future, this can be a parameter (will have to check type - list or dict)
        stored_attrs = {}
        for ak in attr_keep:
            stored_attrs[ak] = self.__dict__[ak]
            self.__dict__[ak] = dict()

        new_rxn = copy.deepcopy(self)
        metabolites = set()
        for ak, av in stored_attrs.items():
            self.__dict__[ak] = av
            new_rxn.__dict__[ak] = {k:v for k,v in av.items()} # pointer: same objects, different bin
            metabolites = metabolites.union(av.keys())
        
        # mutual tracking
        for metabolite in metabolites: 
            metabolite._reaction.add(new_rxn)

        return new_rxn

    def _couple(self, metabolite, type: str):
        """Add coupling coefficient and associated metadata to reaction for a coupled metabolite.

        Parameters
        ----------
        metabolite : Macromolecule
            a macromolecule with an associated coupling coefficient
        type : str
            the type of reactions that are being coupled (one of ['mrna_degradation', 'mrna_formation', 'catalysis',
            'enzyme_degradation'])
        """
        if metabolite.coupling_coefficient is None:
            raise ValueError(
                'Cannot add coupling metadata to reaction for a metabolite without coupling coefficient metadata')
        if type not in metabolite.coupling_coefficient:  # this also checks correct coupling types defined
            raise ValueError('Incorrect coupling coefficient type specified for this metabolite')

        self.add_metabolites({metabolite: metabolite.coupling_coefficient[type]}, combine=True)
        mmap = {m.id: m for m in self._metabolites}
        if metabolite.id in mmap:  # maintain the coupling attributes this way
            self._metabolites[metabolite] = self._metabolites.pop(mmap[metabolite.id])
            metabolite._reaction.add(self)

        # if self.coupled_metabolites == dict():
        #     self.coupled_metabolites = {metabolite: type}
        # else:
        self.coupled_metabolites[metabolite] = type

    def couple(self, metabolites, types: Union[List[str], str]):
        """Add coupling coefficient and associated metadata to reaction for a coupled metabolite.

        Parameters
        ----------
        metabolites : Union[List[Macromolecule], Macromolecule]
            macromolecules with an associated coupling coefficient
        types : Union[List[str], str]
            the type of reactions that are being coupled (one of ['mrna_degradation', 'mrna_formation', 'catalysis',
            'enzyme_degradation'])
        """

        if isinstance(metabolites, list):
            for metabolite, type in dict(zip(metabolites, types)).items():
                self._couple(metabolite, type)
        else:
            self._couple(metabolites, types)

    @staticmethod
    def _check_me_bounds(lb, ub):
        if isinstance(lb, sympy.Expr) or isinstance(ub, sympy.Expr):
            raise ValueError('Reaction bounds can only be a function of mu for reactions of type biomass')

    def build_reaction_string(self, use_metabolite_names: bool = False):
        """Generate a human readable reaction string."""

        def format(number):
            return "" if number == 1 else str(number).rstrip(".") + " "

        id_type = 'id'
        if use_metabolite_names:
            id_type = 'name'
        reactant_bits = []
        product_bits = []
        for met in sorted(self._metabolites, key=attrgetter("id")):
            coefficient = self._metabolites[met]
            name = str(getattr(met, id_type))

            if not isinstance(coefficient, sympy.Expr):
                if coefficient >= 0:
                    product_bits.append(format(coefficient) + name)
                else:
                    reactant_bits.append(format(abs(coefficient)) + name)
            else:
                if float(coefficient.subs(params.mu, 1)) >= 0:
                    product_bits.append(format(coefficient) + name)
                else:
                    reactant_bits.append(format(coefficient).replace('-', '') + name)

        reaction_string = ' + '.join(reactant_bits)
        if not self.reversibility:
            if self.lower_bound < 0 and self.upper_bound <= 0:
                reaction_string += ' <-- '
            else:
                reaction_string += ' --> '
        else:
            reaction_string += ' <=> '
        reaction_string += ' + '.join(product_bits)
        return reaction_string

    # def _map_coupled_metabolites(self):
    #     mmap = {m.id: m for m in self.metabolites}
    #     cm = dict()
    #     for md, type_ in self.coupled_metabolites.items():
    #         cm[mmap[md.id]] = type_
    #     self.coupled_metabolites = cm

    def check_mass_balance(self, tol: SupportsFloat = 0, sympy_tol: SupportsFloat = 1e-15) -> Dict[str, float]:
        """Compute mass and charge balance for the reaction.

        Parameters
        ----------
        tol : SupportsFloat, optional
            precision tolerance for categorizing an element as unbalanced, by default 0
        sympy_tol : SupportsFloat, optional
            sympy.Expr conversions may result in some error, account for this when getting rid of the 
            coupling coefficient value, by default 1e-15

        Returns
        -------
        Dict[str, float]
            dict of {element: amount} for unbalanced elements. "charge" is treated as an element in this dict. 
            This should be empty for balanced reactions.
        """
        reaction_element_dict = defaultdict(int)
        mmap = {m.id: m for m in self._metabolites}

        md = dict()
        for m, c in self.metabolites.items():
            if m.id not in md:
                md[m.id] = c
            else:
                raise ValueError('Same metabolite id, different objects in reaction')  # md[m.id] += c

        # deal with coupled metabolites (also required id mapping above)
        for metabolite, type in self.coupled_metabolites.items():
            md[metabolite.id] -= metabolite.coupling_coefficient[type]  # coupling not part of mass balance
            md[metabolite.id] = float(md[metabolite.id])
            for val in [1, 0]:
                if abs((np.sign(md[metabolite.id]) * val) - md[metabolite.id]) < sympy_tol:
                    md[metabolite.id] = np.sign(md[metabolite.id]) * val

        for m_id, coefficient in iteritems(md):
            metabolite = mmap[m_id]
            if metabolite.charge is not None:
                reaction_element_dict["charge"] += coefficient * metabolite.charge
            for element, amount in iteritems(metabolite.elements):
                reaction_element_dict[element] += coefficient * amount

        return {k: v for k, v in iteritems(reaction_element_dict) if abs(v) > tol}

    def replace_coefficient_mu(self, mu_val: SupportsFloat, inplace: bool = True):
        """Replace mu coefficients with an actual value.

        Parameters
        ----------
        mu_val : SupportsFloat
            value to replace mu with
        inplace : bool, optional
            whether to update value in place or return a copy of the reaction with updated value (False), by default True
        """
        if not mu_val > 0:
            raise ValueError('Mu must be > 0')

        if inplace:
            new_rxn = {k:v for k,v in self.metabolites.items()}
        else:
            new_rxn = self.metabolites.copy()
        
        for met, coeff in self.metabolites.items():
            if isinstance(coeff, sympy.Expr):
                new_rxn[met] = float(coeff.subs(params.mu, mu_val))

        if inplace:
            self.add_metabolites(new_rxn, combine=False)
        else:
            reaction = copy.deepcopy(self)
            reaction.add_metabolites(new_rxn, combine=False)
            return reaction

    @property
    def reactants(self):
        """Return a list of reactants for the reaction."""

        reactants_ = list()
        for k, v in iteritems(self._metabolites):
            if not isinstance(v, sympy.Expr) and v < 0:
                reactants_.append(k)
            elif isinstance(v, sympy.Expr):
                if float(v.subs(params.mu, 1)) < 0:
                    reactants_.append(k)
        return reactants_

    @property
    def products(self):
        """Return a list of products for the reaction."""
        products_ = list()
        for k, v in iteritems(self._metabolites):
            if not isinstance(v, sympy.Expr) and v >= 0:
                products_.append(k)
            elif isinstance(v, sympy.Expr):
                if float(v.subs(params.mu, 1)) >= 0:
                    products_.append(k)
        return products_

    def _add_protein_deg_proxy(self, protein_deg_proxy):
        """Proxy metabolite for protein degradation

        Parameters
        ----------
        protein_deg_proxy : Proxy
            proxy metabolite for coupling protein degradation 
        """
        if not protein_deg_proxy.type == 'proxy':
            raise ValueError('Expected proxy macromolecule')
        if self._protein_deg_proxy:
            raise ValueError('Protein degradation proxy already added')

        self.add_metabolites({protein_deg_proxy: 1})
        self._protein_deg_proxy = True
        self.protein_deg_proxy = protein_deg_proxy


class MetabolicReaction(ME_Reaction):
    """Inherited from ME_Reaction, specifies the metabolic reactions in the model."""

    def __init__(self, id, cobra_id: str, name='', subsystem='', lower_bound=0.0, upper_bound=None):
        """cobra_id specifies the original reaction name in the M-Model"""

        super().__init__(id, name, subsystem, lower_bound, upper_bound)
        self.cobra_id = cobra_id

def to_metabolic_reaction(model_metabolites, reaction: cobra.Reaction, id: Optional[str] = None) -> MetabolicReaction:
    """Convert a cobrapy Reaction to a ME_Model MetabolicReaction.

    Parameters
    ----------
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin
    reaction : cobra.Reaction
        cobra reaction to be converted
    id : str, optional
        reaction id (defaults to reaction.id), by default None

    Returns
    -------
    MetabolicReaction
        Converted reaction
    """
    if type(reaction) != cobra.Reaction:
        raise TypeError('Reaction must be a cobra.Reaction')
    if id is None:
        id = reaction.id
    new_rxn = MetabolicReaction(id=id, cobra_id=reaction.id, name=reaction.name, subsystem=reaction.subsystem,
                            lower_bound=reaction.lower_bound, upper_bound=reaction.upper_bound)
    new_rxn.add_metabolites({model_metabolites.id_object_map[m.id]: stoich for m, stoich in reaction.metabolites.items()}, 
                        combine=False)

    if hasattr(reaction, 'enzyme_compartment'):
        new_rxn.enzyme_compartment = reaction.enzyme_compartment # compartment of the enzyme catalyzing this reaction
    # for k in set(reaction.__dict__.keys()).difference(['_id', 'name', 'subsystem',
    #                                                    '_lower_bound', '_upper_bound',
    #                                                    '_model']):
    #     rxn.__dict__[k] = copy.deepcopy(reaction.__dict__[k])
    return new_rxn


class ExpressionReaction(ME_Reaction):
    """Inherited from ME_Reaction, specifies the expression reactions in the model."""

    def __init__(self, id: str, subsystem: str, name: str = '', lower_bound: SupportsFloat = 0.0, upper_bound: Optional[SupportsFloat] = None,
                 hgnc_id: Optional[str] = None,
                 synthesis: bool = False, synthesis_type: Optional[str] = None, sink: bool = False, sink_type: Optional[str] = None,
                 ubiquitin_biogenesis: bool = False, ribosome_biogenesis: bool = False, trna_charging: bool = False):
        """Initialize the Expression reaction.

        Parameters
        ----------
        id : str
            reaction id
        subsystem : str
            one of ['tRNA_Biogenesis', 'rRNA_expression', 'mRNA_expression', 'Protein_Expression', 'Protein_Degradation', 'Complex_Formation', 'Complex_Degradation']
        name : str, optional
            reaction name, by default ''
        lower_bound : SupportsFloat, optional
            lower flux bound constraint, by default 0.0
        upper_bound : Optional[SupportsFloat], optional
            upper flux bound constraint, by default None
        hgnc_id : Optional[str], optional
            HGNC ID of gene being synthesized, by default None
        synthesis : bool, optional
            whether the reaction represents the "main" synthesis/production for the macromolecule (used for appropriate mapping of coupling)
            intended for use with genes (reactions with an associated hgnc id, and complexes), by default False
        synthesis_type : Optional[str], optional
            one of ['mRNA', 'protein', 'complex'], by default None
            if synthesis is True, the type of macromolecule being synthesized should also be specified
                *for mRNA, the synthesis reaction is coupled to its respective protein translation reaction
                *for proteins and complexes, the synthesis reaction is the final reaction producing the enzyme which
                will be coupled to the metabolic catalysis reaction
        sink : bool, optional
            whether the reaction represents the "main" sink/degradation for the macromolecule, by default False
            intended for use with genes (reactions with an associated hgnc id, and complexes)

                *for mRNA, the degradation reaction will be coupled to the protein translation reaction
                *for proteins and complexes, the degradation reaction will be coupled to the respective
                metabolic catalysis reaction
                *exceptions are synthesis and sink of reactions in ubiquitin_biogenesis (True); these are
                assigned as synthesis and sink to track ubiquitin, but are not themselves coupled to anything
        sink_type : Optional[str], optional
            one of ['mRNA', 'protein', 'complex'], by default None
            if sink True, the type of macromolecule being degraded should also be specified
        ubiquitin_biogenesis : bool, optional
            whether the ExpressionReaction is part of ubiquitin_biogenesis reactions, only used to ignore hgnc_id is None, by default False
        ribosome_biogenesis : bool, optional
            whether the ExpressionReaction is part of ribosome_biogenesis reactions, only used to ignore hgnc_id is None, by default False
        trna_charging : bool, optional
            specifies that the ExpressionReaction is a trna charging reaction (True), for use with calculating biomass change, by default False
        """
        if subsystem not in ['tRNA_Biogenesis', 'rRNA_expression', 'mRNA_expression', 'Protein_Expression',
                             'Protein_Degradation', 'Complex_Formation', 'Complex_Degradation']:
            raise ValueError('Must specify an appropriate expression subsystem')

        super().__init__(id, name, subsystem, lower_bound, upper_bound)

        self.ubiquitin_biogenesis = ubiquitin_biogenesis
        if (not (self.subsystem in ['tRNA_Biogenesis', 'rRNA_expression', 'Complex_Formation',
                                    'Complex_Degradation'] or self.ubiquitin_biogenesis)) and (hgnc_id is None):
            raise ValueError('Must specify hgnc_id of the gene being expressed')

        self.hgnc_id = hgnc_id
        self.synthesis = synthesis
        if self.synthesis and synthesis_type not in ['mRNA', 'protein', 'complex']:
            raise ValueError('The synthesis type must be specified')
        self.synthesis_type = synthesis_type

        self.sink = sink
        if self.sink and sink_type not in ['mRNA', 'protein', 'complex']:
            raise ValueError('The synthesis type must be specified')
        self.sink_type = sink_type

        self.ribosome_biogenesis = ribosome_biogenesis
        self.trna_charging = trna_charging


class ProteinExpressionReaction(ExpressionReaction):
    """Inherited from ExpressionReaction, specifies the protein expression reactions in the model."""

    def __init__(self, id: str, name='', lower_bound: SupportsFloat = 0.0, upper_bound: Optional[SupportsFloat] = None,
                 hgnc_id: Optional[str] = None, translation: bool = False, synthesis: bool = False,
                 ubiquitin_biogenesis: bool = False, ribosome_biogenesis: bool = False):
        """Initialize ProteinExpressionReaction. 

        Parameters
        ----------
        id : str
            reaction id
        name : str, optional
            reaction name, by default ''
        lower_bound : SupportsFloat, optional
            lower flux bound constraint, by default 0.0
        upper_bound : Optional[SupportsFloat], optional
            upper flux bound constraint, by default None
        hgnc_id : [type], optional
            HGNC ID of protein being synthesized, by default None
        translation : bool, optional
            whether the reaction represents the "main" synthesis/production for the protein, by default False
            represents coupling of initial protein product to mRNA (mRNA-->protein coupling)
        synthesis : bool, optional
            whether the reaction represents the "main" synthesis/production for the macromolecule (used for appropriate mapping of coupling)
            intended for use with genes (reactions with an associated hgnc id, and complexes), by default False
        ubiquitin_biogenesis : bool, optional
            whether the ExpressionReaction is part of ubiquitin_biogenesis reactions, only used to ignore hgnc_id is None, by default False
        ribosome_biogenesis : bool, optional
            whether the ExpressionReaction is part of ribosome_biogenesis reactions, only used to ignore hgnc_id is None, by default False
        """
        synthesis_type = None
        if synthesis:
            synthesis_type = 'protein'
        super().__init__(id=id,
                         name=name, lower_bound=lower_bound, upper_bound=upper_bound, hgnc_id=hgnc_id,
                         synthesis=synthesis, synthesis_type=synthesis_type,
                         ubiquitin_biogenesis=ubiquitin_biogenesis, ribosome_biogenesis=ribosome_biogenesis,
                         sink=False, sink_type=None, subsystem='Protein_Expression')
        self.translation = translation

        # self is needed for expressing the associated protein metabolites in the these compartments:
        self._final_compartments = set()


class ProteinDegradationReaction(ExpressionReaction):
    def __init__(self, id: str, hgnc_id: str, sink: bool = False, sink_type: Optional[str] = None,
                 name: str = '', lower_bound: SupportsFloat = 0.0, upper_bound: Optional[SupportsFloat] = None):
        """See ExpressionReaction for parameter details."""
        super().__init__(id=id, subsystem='Protein_Degradation', sink=sink, sink_type=sink_type,
                         name=name, lower_bound=lower_bound, upper_bound=upper_bound, hgnc_id=hgnc_id)
        self._macromolecules = set()  # set of macromolecule ids associated with this degradation reaction
        self._enzymes = set()  # set of enzyme ids associated with this degradation reaction
        self._ribosomal_degradation = False  # see complex_degradation_reaction for details
        self._final_compartments = set()

    def copy(self):
        """Also couple ._macromolecules."""
         # avoid deeopcopying some attributes for:
        # 1) compute time
        # 2) maintaining pointers (same metabolites in the new reaction)

        # additionally, make sure metabolites add the new copied reaction to their attributes (mutual tracking)

        dict_attr_keep = ['_metabolites', 'coupled_metabolites'] #_macromolecules # in the future, this can be a parameter (will have to check type - list or dict)
        list_attr_keep = ['_macromolecules', '_enzymes']
        stored_dict_attrs = {}
        for ak in dict_attr_keep:
            stored_dict_attrs[ak] = self.__dict__[ak]
            self.__dict__[ak] = {}
        
        stored_list_attrs = {}
        for ak in list_attr_keep:
            stored_list_attrs[ak] = self.__dict__[ak]
            self.__dict__[ak] = []

        new_rxn = copy.deepcopy(self)

        metabolites = set()
        for ak, av in stored_dict_attrs.items():
            self.__dict__[ak] = av
            new_rxn.__dict__[ak] = {k:v for k,v in av.items()} # different pointer to bin (dictionary) but not elements in bin for updating of metabolite objecs
            metabolites = metabolites.union(av.keys())
        for ak, av in stored_list_attrs.items():
            self.__dict__[ak] = av
            new_rxn.__dict__[ak] = [item for item in av] # different pointer to bin (list) but not elements in bin
            # metabolites = metabolites.union(av)
        
        # mutual tracking
        for metabolite in metabolites: 
            metabolite._reaction.add(new_rxn)

        return new_rxn

       
        # attr_keep = ['_metabolites', 'coupled_metabolites'] # in the future, this can be a parameter (will have to check type - list or dict)
        # stored_attrs = {}
        # for ak in attr_keep:
        #     stored_attrs[ak] = self.__dict__[ak]
        #     self.__dict__[ak] = dict()

        # new_rxn = copy.deepcopy(self)
        # metabolites = set()
        # for ak, av in stored_attrs.items():
        #     self.__dict__[ak] = av
        #     new_rxn.__dict__[ak] = {k:v for k,v in av.items()} # pointer: same objects, different bin
        #     metabolites = metabolites.union(av.keys())
        
        # # mutual tracking
        # for metabolite in metabolites: 
        #     metabolite._reaction.add(new_rxn)

        return new_rxn

    def _update_tracking(self, macromolecules):
        """Mutual tracking of degradation reactions associated with a macromolecule and vice-versa.

        Parameters
        ----------
        macromolecules : Union[Macromolecule, Set[Macromolecule]]
            marcomolecule to be tracked
        """
        if not isinstance(macromolecules, Iterable):
            macromolecules = {macromolecules}

        for macromolecule in macromolecules:
            macromolecule._degradation_reactions.add(self.id)

        self._macromolecules = self._macromolecules.union(macromolecules)

    def _update_enzymes(self):
        """Update enzymes list to include macromolecules that are classified as enzymes."""
        self._enzymes = {m for m in self._macromolecules if m.enzyme}
        for m in self._enzymes:
            if self.id not in m._degradation_reactions:
                raise ValueError('Improper tracking of degradation reactions and associated macromolecules')

    def _set_proteasomal_degradation(self, **kwargs):
        """For code consistency, mainly for ComplexDegradationReaction, see that method."""
        self.gene_reaction_rule = ' and '.join(mach.proteasome_machinery)


class ComplexDegradationReaction(ExpressionReaction):
    def __init__(self, id: Optional[str] = None, sink: bool = False, sink_type: Optional[str] = None,
                 name: str = '', lower_bound: SupportsFloat = 0.0, upper_bound: Optional[SupportsFloat] = None, hgnc_id: Optional[str] = None):
        """See ExpressionReaction for parameter details."""
        super().__init__(id=id, subsystem='Complex_Degradation', sink=sink, sink_type=sink_type,
                         name=name, lower_bound=lower_bound, upper_bound=upper_bound)
        self._macromolecules = set()  # set of macromolecule ids associated with this degradation reaction
        self._enzymes = set()  # set of enzyme ids associated with this degradation reaction
        self._ribosomal_degradation = False

    def copy(self):
        """Also couple ._macromolecules."""
         # avoid deeopcopying some attributes for:
        # 1) compute time
        # 2) maintaining pointers (same metabolites in the new reaction)

        # additionally, make sure metabolites add the new copied reaction to their attributes (mutual tracking)

        dict_attr_keep = ['_metabolites', 'coupled_metabolites'] #_macromolecules # in the future, this can be a parameter (will have to check type - list or dict)
        list_attr_keep = ['_macromolecules', '_enzymes']
        stored_dict_attrs = {}
        for ak in dict_attr_keep:
            stored_dict_attrs[ak] = self.__dict__[ak]
            self.__dict__[ak] = {}
        
        stored_list_attrs = {}
        for ak in list_attr_keep:
            stored_list_attrs[ak] = self.__dict__[ak]
            self.__dict__[ak] = []

        new_rxn = copy.deepcopy(self)

        metabolites = set()
        for ak, av in stored_dict_attrs.items():
            self.__dict__[ak] = av
            new_rxn.__dict__[ak] = {k:v for k,v in av.items()} # different pointer to bin (dictionary) but not elements in bin for updating of metabolite objecs
            metabolites = metabolites.union(av.keys())
        for ak, av in stored_list_attrs.items():
            self.__dict__[ak] = av
            new_rxn.__dict__[ak] = [item for item in av] # different pointer to bin (list) but not elements in bin
            # metabolites = metabolites.union(av)
        
        # mutual tracking
        for metabolite in metabolites: 
            metabolite._reaction.add(new_rxn)

        return new_rxn

    def _update_tracking(self, macromolecules):
        """Mutual tracking of degradation reactions associated with a macromolecule and vice-versa.

        Parameters
        ----------
        macromolecules : Union[Macromolecule, Set[Macromolecule]]
            marcomolecule to be tracked
        """
        if not isinstance(macromolecules, Iterable):
            macromolecules = {macromolecules}

        for macromolecule in macromolecules:
            macromolecule._degradation_reactions.add(self.id)

        self._macromolecules = self._macromolecules.union(macromolecules)

    def _update_enzymes(self):
        """Update enzymes list to include macromolecules that are classified as enzymes."""
        self._enzymes = {m for m in self._macromolecules if m.enzyme}
        for m in self._enzymes:
            if self.id not in m._degradation_reactions:
                raise ValueError('Improper tracking of degradation reactions and associated macromolecules')

    def _set_proteasomal_degradation(self, macromolecule, ribosomal_complex: bool):
        """Quick addition of attribute for build_me script, since current format has the machinery for
        the proteosomal degradation different than standard complexes (to degrade rRNAs as well).
        Change in machinery hard-coded into degradation.degrade script and double-checked in build_me script.

        Parameters
        ----------
        macromolecule : Union[Protein, Complex]
        ribosomal_complex : bool
            whether the macromolecule is a ribosomal complex (True) or not (False)
        """

        if not ribosomal_complex:
            machinery_ = mach.proteasome_machinery
        else:
            self._ribosomal_degradation = True
            # hard-coded

            # Option 1: degrade rRNA with ribosomal degradation - see also expression.gene_expression.protein_expression.degradation.proteasomal_degradation
            machinery_ = list(set(mach.proteasome_machinery + mach.exosome['HGNC ID (gene)'].tolist()))

        #             # # Option 2: degrade proteins with ribosomal degradation, releasing rRNA as intact - see also expression.gene_expression.protein_expression.degradation.proteasomal_degradation
        #             machinery_ = mach.proteasome_machinery

        # this is more flexible, but hardcoded check in degradation.proteasomal_degradation
        # renders this unecessary (keep for future iterations)

        #             # add machinery
        #             counter = 0
        #             machinery_ = list()
        #             if len(set([m for m in mdc if m.type == 'protein'])) > 0:
        #                 machinery_ += mach.proteasome_machinery
        #                 counter += 1
        #             if len(rm) > 0:
        #                 machinery_ += mach.exosome['HGNC ID (gene)'].tolist() #rrna degradation machinery
        #                 counter += 1

        #             if counter != 2:
        #                 err = 'Internal: Only expect mature ribosome complex with rRNA and protein to be degraded.
        #                 err += 'Should work with current code, but double check'
        #                 raise ValueError(err)

        self.gene_reaction_rule = ' and '.join(machinery_)

class BiomassReaction(cobra.Reaction):
    """Specifies biomass reactions in the model, allowing reaction bounds to be a function of mu."""

    def __init__(self, id: str, name: str = '', subsystem: str = '',
                 lower_bound: Union[SupportsFloat, sympy.Expr] = 0.0, upper_bound: Optional[Union[SupportsFloat, sympy.Expr]] = None):
        super().__init__(id, name, subsystem, lower_bound, upper_bound)

    @staticmethod
    def _check_me_bounds(lb, ub):
        if isinstance(lb, sympy.Expr) or isinstance(ub, sympy.Expr):
            if params.mu not in lb.free_symbols and params.mu in ub.free_symbols:
                raise ValueError(
                    'Currently, if reaction bounds are a function of mu, they must be for both the upper and lower bound')

    def replace_bound_mu(self, mu_val: SupportsFloat = 1, values: Optional[List[sympy.Expr]] = None, inplace: bool = False, _ub: bool = True):
        """Gives numeric values to bounds. Assums growth is always > 0. 

        Parameters
        ----------
        mu_val : SupportsFloat, optional
            The value for mu to replace the bounds that contain a mu expression with, by default 1
        values : Optional[List[sympy.Expr]], optional
            Each entry is an expression containing mu to be replaced by mu_val; inplace must be False, by default None
        inplace : bool, optional
            Whether to replace the bounds inplace on the reaction object (True), or return the bounds, by default False
        _ub : bool, optional
            internal use, whether to use cobra.Reaction._upper_bound or cobra.Reaction.upper_bound, by default True
        """
        if _ub:
            lb, ub = copy.copy(self._lower_bound), copy.copy(self._upper_bound)
        else:
            lb, ub = copy.copy(self.lower_bound), copy.copy(self.upper_bound)

        self._check_me_bounds(lb, ub)

        if isinstance(lb, sympy.Expr):  # _check_me_bounds makes sure both lb and ub are symp.Expr objects
            # replace growth with input mu val (assuming growth always > 0)
            lb, ub = float(lb.subs(params.mu, mu_val)), float(ub.subs(params.mu, mu_val))
        #         else:
        #             warnings.warn('Bounds do not have a mu value')

        if values is None:
            if not inplace:
                return lb, ub
            self._lower_bound, self._upper_bound = lb, ub
        else:
            if not isinstance(values, list):
                raise TypeError('values must a list')

            for i, value in enumerate(values):
                if isinstance(value, sympy.Expr):  # assumes the sympy expression always contains params.mu
                    values[i] = float(value.subs(params.mu, mu_val))
            if not inplace:
                return lb, ub, values
            raise ValueError('Either values must be None or inplace False')

    @property
    def reversibility(self) -> bool:
        """Whether the reaction can proceed in both directions (reversible). This is computed from the current upper and lower bounds."""
        lb, ub = self.replace_bound_mu()
        return lb < 0 < ub

    def build_reaction_string(self, use_metabolite_names: bool = False):
        """Generate a human readable reaction string."""
        def format(number):
            return "" if number == 1 else str(number).rstrip(".") + " "

        id_type = 'id'
        if use_metabolite_names:
            id_type = 'name'
        reactant_bits = []
        product_bits = []
        for met in sorted(self._metabolites, key=attrgetter("id")):
            coefficient = self._metabolites[met]
            name = str(getattr(met, id_type))
            if coefficient >= 0:
                product_bits.append(format(coefficient) + name)
            else:
                reactant_bits.append(format(abs(coefficient)) + name)

        reaction_string = ' + '.join(reactant_bits)
        if not self.reversibility:
            lb, ub = self.replace_bound_mu(_ub=False)
            if lb < 0 and ub <= 0:
                reaction_string += ' <-- '
            else:
                reaction_string += ' --> '
        else:
            reaction_string += ' <=> '
        reaction_string += ' + '.join(product_bits)
        return reaction_string
