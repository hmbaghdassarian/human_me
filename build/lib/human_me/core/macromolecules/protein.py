#!/usr/bin/env python
# coding: utf-8
import copy
from typing import Dict, Optional

from human_me.core.macromolecules.macromolecule import Macromolecule, Proxy


class Protein(Macromolecule):
    """Protein macromolecule object representation."""

    type = 'protein'

    def __init__(self, compartment: str, id_: str, model_metabolites, gene_info = None,
                amino_acid_counts: Optional[Dict[str,int]] = None, dummy: bool = False,
                 non_machinery: bool = False):
        """Generates a Macromolecule in the compartment for a protein. Generated either from gene_info or id_ and amino_acid_counts.


        Parameters
        ----------
        compartment : str
            protein subcellular location (one-letter code)
        id_ : str
            protein identifier
        model_metabolites : utils.metabolites.MetaboliteBin
            the me_input_model metabolites as specified by MetaboliteBin
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

        self.model_metabolites = model_metabolites
        charge = sum(
            [self.model_metabolites.seq_amino_acid_map_compartments[compartment][aa_code].charge * aa_count for aa_code, aa_count in
             self._amino_acid_counts.items()])

        elements = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}
        if compartment in self.model_metabolites.seq_amino_acid_map_compartments:
            for aa_code, aa_count in self._amino_acid_counts.items():
                aa_elements = self.model_metabolites.seq_amino_acid_map_compartments[compartment][aa_code].elements
                for element in aa_elements:
                    elements[element] += aa_count * aa_elements[element]
        else:
            raise ValueError('Internal: Must add ' + compartment + ' compartment to amino acid map in metab')

        # peptide bond formation
        elements['H'] -= 2 * (self.length - 1)
        elements['O'] -= 1 * (self.length - 1)

        super().__init__(id=id_, compartment=compartment, charge=charge, elements=elements,
                               hgnc_id=self.hgnc_id)

        self.dummy = dummy
        self.dummy_type = None # string to specify whether orphan or unmodeled dummy protein, assigned in me building
        self.enzyme = False  # whether the protein is involved in catalysis of a reaction
        self.keff = None
        self._degradation_reactions = set() # associated degradation reactions for protein monomer, if any
        self.non_machinery = non_machinery

    def make_proxy(self):
        """Make a proxy metabolite for coupling enzyme degradation to reaction catalysis."""
        if self.non_machinery:
            raise ValueError('Unexpected generation of coupling proxy metabolites for non-machinery')
        return Proxy(associated_macromolecule=self)
    
    def copy(self):
        """Overwrite cobra.Species.copy."""
        # copy while preserving model_metabolites pointer to avoid memory issues and not copying it to avoid speed issues
        model_metabolites = self.model_metabolites
        del self.model_metabolites

        new_macromolecule = copy.deepcopy(self)
        
        self.model_metabolites = model_metabolites
        new_macromolecule.model_metabolites = self.model_metabolites # retain same pointer to avoid memory issues

        return new_macromolecule
