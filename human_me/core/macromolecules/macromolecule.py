#!/usr/bin/env python
# coding: utf-8
from typing import Dict, Optional, SupportsFloat, Union

import cobra
import sympy

from human_me.utils import parameters as params


class Macromolecule(cobra.Metabolite):
    def __init__(self, id: Optional[str] = None, formula: Optional[str] = None, name: str = "", charge: Optional[int] = None,
                 compartment: Optional[str] = None, elements: Optional[Dict[str, int]] = None,
                 hgnc_id: Optional[str] = None):
        """RNA, proteins, complexes, and coupling proxies
        see cobra.Metabolite for undescribed parameters
        Parameters
        ----------
        id : str, optional
            macromolecule id, by default None
        formula : Optional[str], optional
            [description], by default None
        name : str, optional
            [description], by default ""
        charge : Optional[int], optional
            [description], by default None
        compartment : str, optional
            one-letter code for macromolecule compartment, by default None
        elements : Optional[Dict[str, int]], optional
            keys are the atomic element, values are the number of occurences of that element in the macromolecule, by default None
        hgnc_id : str, optional
            the associated gene HGNC ID of the macromolecule (HGNC:####), by default None
        """

        if compartment not in params.compartments:
            err = 'Specified compartment is not considered in the ME Model. Please input one of the following: '
            err += ', '.join(list(params.compartments))
            raise ValueError(err)

        cobra.Metabolite.__init__(self, id=id, charge=charge, compartment=compartment, formula=formula, name=name)

        if elements is not None:
            self.elements = elements
        if self.id.split('_')[-1] != self.compartment:
            raise ValueError('Macromolecules must syntactically have compartment as part of id')

        self.coupling_coefficient = None
        self.hgnc_id = hgnc_id

    def change_compartment(self, new_compartment: str):
        """Create a macromolecule the same as self, but in a different compartment

        Parameters
        ----------
        new_compartment : str
            one-letter code of compartment to change to 

        Returns
        -------
        Macromolecule
            a copy of the macromolecule in the new compartment
        """

        if new_compartment == self.compartment:
            raise ValueError('The macromolecule is already in this compartment')
        if new_compartment not in params.compartments.keys():
            err = 'Specified compartment is not considered in the ME Model. Please input one of the following: '
            err += ', '.join(list(params.compartments.keys()))

        new_macromolecule = self.copy()
        new_macromolecule.id = '_'.join(self.id.split('_')[:-1]) + '_' + new_compartment
        new_macromolecule.compartment = new_compartment

        return new_macromolecule

    def couple(self, type: str, value: Union[SupportsFloat, sympy.Expr]):
        """Stores coupling information for macromolecule

        Parameters
        ----------
        type : str
            the type of reaction this macromolecule is coupled to
        value : Union[SupportsFloat, sympy.Expr]
            the coupling coefficient value

        Returns
        ----------
        self.coupling_coefficient: Dict[str, Union[SupportsFloat, sympy.Expr]]
            dictionary of length one, key is the type, value is the coupling coefficient
        """
        if hasattr(self, 'non_machinery') and self.non_machinery:
            raise ValueError('Unexpected coupling of a non_machinery protein:' + self.id)
        if type not in ['catalysis', 'enzyme_degradation', 'mrna_degradation', 'mrna_formation']:
            raise ValueError(
                'The couple id must be one of catalysis, mrna_degradation, enzyme_degradation, or mrna_formation')

        if self.coupling_coefficient is None:
            self.coupling_coefficient = {type: value}
        else:
            if [type] != list(self.coupling_coefficient):
                raise ValueError('More than one coupling type associated with macromolecule: ' + self.id)
            if value != self.coupling_coefficient[type]:
                raise ValueError('More than one coupling coefficient value associated with macromolecule: ' + self.id)

        if type == 'catalysis':
            self.enzyme = True


class Proxy(Macromolecule):
    """Proxy macromolecules for c2/c4 coupling of degradation"""

    def __init__(self, associated_macromolecule: Macromolecule):
        """Init method for Proxy.

        Parameters
        ----------
        associated_macromolecule : Macromolecule
           the c1/c3 associated macromolecule to the respective c2/c4 coupling
        """
        if associated_macromolecule.type not in ['mrna', 'protein', 'complex']:
            raise ValueError('Unexpected associated macromolecule for proxy metabolite')
        key_mapper = {'mrna': 'mrna_degradation', 'protein': 'enzyme_degradation',
                      'complex': 'enzyme_degradation'}

        id_ = associated_macromolecule.hgnc_id if associated_macromolecule.type == 'mrna' else associated_macromolecule._deg_id
        Macromolecule.__init__(self, id='_'.join([id_,
                                                  key_mapper[associated_macromolecule.type], 'proxy',
                                                  associated_macromolecule.compartment]),
                               compartment=associated_macromolecule.compartment,
                               hgnc_id=associated_macromolecule.hgnc_id)
        self.type = 'proxy'
        self.associated_macromolecule = associated_macromolecule.id
        self._amt = associated_macromolecule.type
        # only for complexes, used in ExpressedGene class
        self._complex_hgnc_ids = [p.hgnc_id for p in associated_macromolecule.decompose_complex() if
                                  p.hgnc_id is not None] if self._amt == 'complex' else []

    def couple(self, value: Union[SupportsFloat, sympy.Expr]):
        """Stores coupling information for proxy. 

        Parameters
        ----------
        value : Union[SupportsFloat, sympy.Expr]
            the coupling coefficient value

        Returns
        ----------
        self.coupling_coefficient: Dict[str, Union[SupportsFloat, sympy.Expr]]
            dictionary of length one, key is the type, value is the coupling coefficient
        """

        key_mapper = {'mrna': 'mrna_degradation', 'protein': 'enzyme_degradation',
                      'complex': 'enzyme_degradation'}

        if self.coupling_coefficient is None:
            self.coupling_coefficient = {key_mapper[self._amt]: value}
        else:
            if value != self.coupling_coefficient[key_mapper[self._amt]]:
                raise ValueError('More than one coupling coefficient value associated with macromolecule: ' + self.id)
