#!/usr/bin/env python
# coding: utf-8

from typing import Dict, Optional

from human_me.utils import metabolites as metab
from human_me.core.macromolecules.macromolecule import Macromolecule, Proxy


class Protein(Macromolecule):
    """Protein macromolecule object representation."""

    def __init__(self, compartment: str, id_: str, gene_info = None,
                amino_acid_counts: Optional[Dict[str,int]] = None, dummy: bool = False,
                 non_machinery: bool = False):
        """Generates a Macromolecule in the compartment for a protein. Generated either from gene_info or id_ and amino_acid_counts.


        Parameters
        ----------
        compartment : str
            protein subcellular location (one-letter code)
        id_ : str
            protein identifier
        gene_info : GeneInformation, optional
            gene's associated GeneInformation object, by default None
        amino_acid_counts : Dict[str, int], optional
            keys are one-letter amino acid codes, values are the number of occurences of that amino acid in the protein sequence, by default None
        dummy : bool, optional
             whether the protein is a dummy protein for the unmodeled protein fraction of the ME-Model, by default False
        non_machinery : bool, optional
            whether the protein metabolite is non_machinery, by default False
            *Note, applies to final protein product but not intermediates
        """
        if gene_info is not None and (amino_acid_counts is not None):
            raise ValueError('Please specify either gene_info only or amino_acid_counts only')
        if gene_info is None and ((id_ is None) or (amino_acid_counts is None)):
            raise ValueError('Please specify either gene_info or id_/amino_acid_counts')
        if id_ is None:
            raise ValueError('Unaccounted for condition in protein id naming')

        self.hgnc_id = None
        if gene_info is not None:
            self.hgnc_id = gene_info.hgnc_id
            # for degradation
            self.length = gene_info.L_protein
            self._amino_acid_counts = gene_info.amino_acid_counts
            self._ptms = gene_info.ptms
            self._deg_id = gene_info.hgnc_id
            # for rest of pipeline
            self.k_deg = gene_info.coupling_params['alpha_p']
            id_ = gene_info.hgnc_id + '_' + id_ + '_protein_' + compartment

        else:
            self._amino_acid_counts = amino_acid_counts
            self.length = sum(amino_acid_counts.values())
            id_ = id_ + '_protein_' + compartment

        charge = sum(
            [metab.seq_amino_acid_map_compartments[compartment][aa_code].charge * aa_count for aa_code, aa_count in
             self._amino_acid_counts.items()])

        elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}
        if compartment in metab.seq_amino_acid_map_compartments:
            for aa_code, aa_count in self._amino_acid_counts.items():
                aa_elements = metab.seq_amino_acid_map_compartments[compartment][aa_code].elements
                for element in aa_elements:
                    elements[element] += aa_count * aa_elements[element]
        else:
            raise ValueError('Internal: Must add ' + compartment + ' compartment to amino acid map in metab')

        # peptide bond formation
        elements['H'] -= 2 * (self.length - 1)
        elements['O'] -= 1 * (self.length - 1)

        Macromolecule.__init__(self, id=id_, compartment=compartment, charge=charge, elements=elements,
                               hgnc_id=self.hgnc_id)

        self.type = 'protein'
        self.dummy = dummy
        self.enzyme = False  # whether the protein is involved in catalysis of a reaction
        self.keff = None
        self._degradation_reactions = []  # associated degradation reactions for protein monomer, if any
        self.non_machinery = non_machinery

    def _consolidate_degradation_rxns(self):
        """Remove redundant IDs"""
        self._degradation_reactions = list(set(self._degradation_reactions))

    def make_proxy(self):
        """Make a proxy metabolite for coupling enzyme degradation to reaction catalysis"""
        if self.non_machinery:
            raise ValueError('Unexpected generation of coupling proxy metabolites for non-machinery')
        return Proxy(associated_macromolecule=self)
