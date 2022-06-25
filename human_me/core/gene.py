#!/usr/bin/env python
# coding: utf-8

from typing import Union, SupportsFloat

from human_me.core.reaction import MetabolicReaction, ExpressionReaction
from human_me.core.macromolecules.complex import Complex
from human_me.utils.functions import flatten_list

ribosomal_genes = {'HGNC:11325', 'HGNC:11740', 'HGNC:10369', 'HGNC:10414', 'HGNC:11326', 'HGNC:10448', 'HGNC:10421',
                   'HGNC:3300', 'HGNC:10346', 'HGNC:12458', 'HGNC:10304', 'HGNC:10389', 'HGNC:10311', 'HGNC:13631',
                   'HGNC:11301', 'HGNC:10410', 'HGNC:10312', 'HGNC:10364', 'HGNC:10359', 'HGNC:10759', 'HGNC:10305',
                   'HGNC:11323', 'HGNC:28962', 'HGNC:10440', 'HGNC:10385', 'HGNC:10313', 'HGNC:10404', 'HGNC:10347',
                   'HGNC:10371', 'HGNC:17976', 'HGNC:10298', 'HGNC:17094', 'HGNC:10387', 'HGNC:10327', 'HGNC:10416',
                   'HGNC:10397', 'HGNC:11324', 'HGNC:10328', 'HGNC:11302', 'HGNC:10426', 'HGNC:10384', 'HGNC:10307',
                   'HGNC:10386', 'HGNC:10360', 'HGNC:21082', 'HGNC:7670', 'HGNC:10377', 'HGNC:10396', 'HGNC:11300',
                   'HGNC:10372', 'HGNC:10344', 'HGNC:3214', 'HGNC:10442', 'HGNC:23401', 'HGNC:23400', 'HGNC:10301',
                   'HGNC:10441', 'HGNC:10348', 'HGNC:10325', 'HGNC:20090', 'HGNC:10340', 'HGNC:10331', 'HGNC:10362',
                   'HGNC:10332', 'HGNC:10345', 'HGNC:26212', 'HGNC:10317', 'HGNC:24624', 'HGNC:11307', 'HGNC:10413',
                   'HGNC:10363', 'HGNC:10388', 'HGNC:10429', 'HGNC:10419', 'HGNC:10383', 'HGNC:10302', 'HGNC:10333',
                   'HGNC:10411', 'HGNC:10418', 'HGNC:18276', 'HGNC:11846', 'HGNC:11299', 'HGNC:3208', 'HGNC:10299',
                   'HGNC:10417', 'HGNC:17718', 'HGNC:3189', 'HGNC:10354', 'HGNC:10401', 'HGNC:10336', 'HGNC:16993',
                   'HGNC:10330', 'HGNC:10368', 'HGNC:10334', 'HGNC:10329', 'HGNC:18501', 'HGNC:18277', 'HGNC:10353',
                   'HGNC:10424', 'HGNC:3597', 'HGNC:10425', 'HGNC:11303', 'HGNC:17050', 'HGNC:10316', 'HGNC:18476',
                   'HGNC:5238', 'HGNC:10315', 'HGNC:10349', 'HGNC:21370', 'HGNC:10409', 'HGNC:10405', 'HGNC:10351',
                   'HGNC:10420', 'HGNC:16931', 'HGNC:10306', 'HGNC:10402', 'HGNC:10350'}

# checks work as follows: maximum limits checked within each method, minimum limit checked in .check methods
class ExpressedGene:
    """Tracks all reactions and macromolecules associated with a ME Model gene.Designed to be used after building the full ME_Model."""

    def __init__(self, hgnc_id: str):
        """Init method for ExpressedGene.

        Parameters
        ----------
        hgnc_id:  str
            the gene HGNC ID
        """
        if not hgnc_id.startswith('HGNC:'):
            raise ValueError('Currently all genes must be in standard HGNC ID format')
        self.hgnc_id = hgnc_id
        # initialize reactions
        self.reactions = {'Catalysis_Reactions': {'Metabolic_Module': dict(), 'Expression_Module': dict()},
                          'ExpressionReactions': {'mrna': {'synthesis': None, 'sink': None, 'other': []},
                                                  'protein': {'translation': [], 'synthesis': [],
                                                              'sink': [], 'other': []},
                                                  'complex': {'synthesis': dict(), 'sink': dict(), 'other': dict()}}}
        self.ubiquitin_biogenesis = False
        self.ribosome_biogenesis = False

        if self.hgnc_id != 'HGNC:DUMMYUNMODELED':
            self._is_unmodeled_protein = False
        else:
            self._is_unmodeled_protein = True

        self.macromolecules = {'RNA': {'premrna': None, 'mrna': {'coupled': {}, 'other': None}, 'lariat': None},
                               'Protein': {'coupled': {}, 'other': [], 'non-machinery': []},
                               'Complex': {'coupled': {}, 'other': []},
                               'Proxy': {'mrna_degradation': {},
                                         'enzyme_degradation': {}}
                               }
        self._summarized = False # for __repre__ method
    # REACTIONS--------------------------------------------------------------------------------
    def add_reaction(self, r: Union[MetabolicReaction, ExpressionReaction]):
        """Organizes ME_Model reaction into self.reactions attribute.

        Parameters
        ----------
        r : Union[MetabolicReaction, ExpressionReaction]
            the reaction to add
        """
        catalysis = self.is_catalyzing(r)
        if catalysis:
            self._add_catalysis_reaction(r)

        expression = False
        not_expressing = ['tRNA_Biogenesis',
                          'rRNA_expression']  # subsystems that are not expressing any hgnc_id macromolecules, but do have catalysis
        if isinstance(r,
                      ExpressionReaction) and r.subsystem not in not_expressing:  # not elif, can be catalyzing its own expression for certain subsystems
            expression = self._add_expression_reaction(r)

        if not catalysis and not expression:
            raise ValueError('The reaction ' + r.id + ' does not appear to be associated with the gene ' + self.hgnc_id)

    def is_catalyzing(self, r) -> bool:
        """Determines whether the gene is involved in catalysis of the reaction.

        Parameters
        ----------
        r : ME_Reaction
            [description]

        Returns
        -------
        catalysis : bool
            True if the gene is involved in catalyzing the reaction, False otherwise
        """
        catalysis = True
        if not hasattr(r, '_ribosomal_degradation') or not (r._ribosomal_degradation and r.sink):
            assoc_macro = {t: m for m, t in r.coupled_metabolites.items()}
            if 'catalysis' in assoc_macro:
                if not isinstance(assoc_macro['catalysis'], Complex):  # monomers
                    for m, t in r.coupled_metabolites.items():
                        if t not in ['catalysis', 'enzyme_degradation']:
                            raise ValueError('Unexpected coupling type in catalysis reaction for ' + r.id)
                        if m.hgnc_id != self.hgnc_id:
                            catalysis = False
                else:  # complexes
                    if self.hgnc_id not in [m.hgnc_id for m in assoc_macro['catalysis'].decompose_complex()]: # TODO: .decompose_complex is very slow
                        catalysis = False
            else:
                catalysis = False
        else:
            if self.hgnc_id not in set(flatten_list([[p.hgnc_id for p in cplx.decompose_complex()] for cplx in
                                                 [m for m, t in r.coupled_metabolites.items() if t == 'catalysis']])):
                catalysis = False
        return catalysis

    def _add_catalysis_reaction(self, r: Union[MetabolicReaction, ExpressionReaction]):
        """Adds reactions that the gene catalyzes, splitting by whether it is catalyzing a metabolic orexpression module reaction.

        Hierarchy is organized as follows: {reaction_id: {catalysis: enzyme_id, deg_proxy: proxy_id}}. 

        Where catalysis: enzyme_id represents the macromolecules enzyme that is coupled to the reaction for
        protein synthesis to reaction catalysis and deg_proxy: proxy_id represents the proxy macromolecule that couples
        protein degradation to reaction catalysis.

        Parameters
        ----------
        r : Union[MetabolicReaction, ExpressionReaction]
            catalysis reaction being added
        """
        if isinstance(r, MetabolicReaction):
            self.reactions['Catalysis_Reactions']['Metabolic_Module'][r.id] = {m.id: t for m, t in
                                                                               r.coupled_metabolites.items()}
        else:
            if not hasattr(r, '_ribosomal_degradation') or not (r._ribosomal_degradation and r.sink):
                self.reactions['Catalysis_Reactions']['Expression_Module'][r.id] = {m.id: t for m, t in
                                                                                    r.coupled_metabolites.items()}
            else:
                if r.id not in self.reactions['Catalysis_Reactions']['Expression_Module']:
                    self.reactions['Catalysis_Reactions']['Expression_Module'][r.id] = dict()
                for m, t in r.coupled_metabolites.items():
                    if m.type == 'complex' and self.hgnc_id in [p.hgnc_id for p in m.decompose_complex()]:
                        self.reactions['Catalysis_Reactions']['Expression_Module'][r.id][m.id] = t
                    elif m.type == 'proxy' and self.hgnc_id in m._complex_hgnc_ids:
                        self.reactions['Catalysis_Reactions']['Expression_Module'][r.id][m.id] = t
                if len(self.reactions['Catalysis_Reactions']['Expression_Module'][r.id]) != 2:
                    raise ValueError(
                        'Unexpected number of coupled macromolecules for ribosomal degradation: ' + self.hgnc_id)

    def _assign_biogenesis(self, r):
        if r.ubiquitin_biogenesis:
            self.ubiquitin_biogenesis = True
        if r.ribosome_biogenesis:
            self.ribosome_biogenesis = True

    def _add_expression_reaction(self, r, tol: SupportsFloat = 1e-17):
        """Reactions involving expression of a gene (not catalysis, even if it is catalysis of an expression-module reaction).

        Parameters
        ----------
        tol : SupportsFloat
            tolerance threshold for determining whether a complex is self-catalysing the reaction (exceptional cases)
        """
        expression = True

        # mrna
        if r.subsystem == 'mRNA_expression':
            if r.hgnc_id != self.hgnc_id:
                expression = False

            else:
                self._assign_biogenesis(r)
                if r.synthesis:
                    if self.reactions['ExpressionReactions']['mrna']['synthesis'] is not None:
                        raise ValueError('Multiple ' + 'mrna' + '  synthesis reactions assigned to ' + self.hgnc_id)
                    self.reactions['ExpressionReactions']['mrna']['synthesis'] = r.id
                elif r.sink:
                    if self.reactions['ExpressionReactions']['mrna']['sink'] is not None:
                        raise ValueError('Multiple ' + 'mrna' + '  sink reactions assigned to ' + self.hgnc_id)
                    self.reactions['ExpressionReactions']['mrna']['sink'] = r.id
                else:
                    self.reactions['ExpressionReactions']['mrna']['other'] += [r.id]

                    # protein
        # multiple compartments allow multiple sink/synthesis reactions
        elif r.subsystem.startswith('Protein_'):
            if r.hgnc_id != self.hgnc_id:
                expression = False
            else:
                self._assign_biogenesis(r)
                if r.synthesis:
                    self.reactions['ExpressionReactions']['protein']['synthesis'] += [r.id]
                elif r.sink:
                    self.reactions['ExpressionReactions']['protein']['sink'] += [r.id]
                elif hasattr(r, 'translation') and r.translation:
                    if len(self.reactions['ExpressionReactions']['protein']['translation']) > 1:
                        raise ValueError(
                            'More than two ' + 'protein' + '  translation reactions assigned to ' + self.hgnc_id)
                    self.reactions['ExpressionReactions']['protein']['translation'] += [r.id]
                else:
                    self.reactions['ExpressionReactions']['protein']['other'] += [r.id]

        # complex
        elif r.subsystem.startswith('Complex_'):
            complexes = list()
            for cplx in r.metabolites:
                if isinstance(cplx, Complex):
                    if self.hgnc_id in [m.hgnc_id for m in cplx.decompose_complex()]:
                        if cplx in r.coupled_metabolites:  # account for self-catalysis
                            # note that the two internal abs() shouldn't need to be there, but in case of
                            # inconsistent formatting of coupling coefficient values
                            if abs(r.metabolites[cplx] - cplx.coupling_coefficient[r.coupled_metabolites[cplx]]) > tol:
                                complexes += [cplx]
                        else:
                            complexes += [cplx]
            #             complexes = [cplx for cplx in r.products if isinstance(cplx, Complex) and cplx not in r.coupled_metabolites \
            #                         and self.hgnc_id in [m.hgnc_id for m in cplx.components]]
            if len(complexes) == 0:
                expression = False
            else:
                if r.synthesis:
                    complex_product = [c for c in complexes if c in r.products]
                    if len(complex_product) > 1:
                        raise ValueError('Unexpected synthesis of multiple complex products')
                    complex_product = complex_product[0]
                    if complex_product.id in self.reactions['ExpressionReactions']['complex']['synthesis']:
                        raise ValueError(
                            'Multiple ' + 'complex' + '  synthesis reactions assigned to ' + complex_product.id)

                    self.reactions['ExpressionReactions']['complex']['synthesis'][complex_product.id] = r.id
                    for cplx in set(complexes).difference([complex_product]):
                        self._add_other_complex_expression(cplx, r)
                elif r.sink:
                    complex_reactant = [c for c in complexes if c in r.reactants]
                    if len(complex_reactant) > 1:
                        raise ValueError('Unexpected degradation of multiple complex reactants')
                    complex_reactant = complex_reactant[0]
                    if complex_reactant.id in self.reactions['ExpressionReactions']['complex']['sink']:
                        raise ValueError(
                            'Multiple ' + 'complex' + '  sink reactions assigned to ' + complex_reactant.id)

                    self.reactions['ExpressionReactions']['complex']['sink'][complex_reactant.id] = r.id
                    for cplx in set(complexes).difference([complex_reactant]):
                        self._add_other_complex_expression(cplx, r)
                else:
                    self._add_other_complex_expression(cplx, r)
        else:
            raise ValueError('Unaccounted for Expression Reaction subsystem')
        return expression

    def _add_other_complex_expression(self, cplx, r):
        if cplx.id not in self.reactions['ExpressionReactions']['complex']['other']:
            self.reactions['ExpressionReactions']['complex']['other'][cplx.id] = [r.id]
        else:
            self.reactions['ExpressionReactions']['complex']['other'][cplx.id] += [r.id]

    def _check_reactions(self):

        if len(self.reactions['Catalysis_Reactions']['Metabolic_Module']) == 0 and len(
                self.reactions['Catalysis_Reactions']['Expression_Module']) == 0 and not (
                self.ubiquitin_biogenesis or self._is_non_machinery_only or self._is_unmodeled_protein):
            raise ValueError('No catalysis reactions associated with: ' + self.hgnc_id)
        mrna = self.reactions['ExpressionReactions']['mrna']
        protein = self.reactions['ExpressionReactions']['protein']
        complex_ = self.reactions['ExpressionReactions']['complex']

        if mrna['synthesis'] is None:
            raise ValueError('No mrna synthesis reactions associated with: ' + self.hgnc_id)
        if mrna['sink'] is None:
            raise ValueError('No mrna degradation reactions associated with: ' + self.hgnc_id)
        if len(protein['translation']) == 0:
            raise ValueError('No protein translation reactions associated with: ' + self.hgnc_id)
        if len(protein['synthesis']) == 0 and len(complex_['synthesis']) == 0 and not self._is_non_machinery_only:
            raise ValueError('No enzyme synthesis reactions associated with: ' + self.hgnc_id)

        rib_deg_coupling = self.ribosome_biogenesis or self.hgnc_id in ribosomal_genes
        if len(protein['sink']) == 0 and len(complex_['sink']) == 0 and not self.ubiquitin_biogenesis \
                and self._enzyme_compartments != ['e'] and not self._is_unmodeled_protein and not rib_deg_coupling:  # TODO: TEMPORARY - no ribosomal degradation coupling
            raise ValueError('No enzyme degradation reactions associated with: ' + self.hgnc_id)

        # all proteins have atleast 1 translation reaction
        # appropriate couplings between reactions

    # MACROMOLECULES--------------------------------------------------------------------------------
    def add_macromolecule(self, m):
        """Updates self.macromolecules with m.

        Parameters
        ----------
        m: Macromolecule
            the macromolecule to add
        """
        if m.type not in ['complex', 'proxy']:
            if not hasattr(m, 'hgnc_id') or m.hgnc_id is None or m.hgnc_id != self.hgnc_id:
                raise ValueError(
                    'The macromolecule ' + m.id + ' does not appear to be associated with the gene ' + self.hgnc_id)
        elif m.type == 'complex' and self.hgnc_id not in [m.hgnc_id for m in m.decompose_complex() if
                                                          m.type == 'protein']:
            raise ValueError('The macromolecule ' + m.id + ' does not appear to be associated with the gene ' + self.hgnc_id)
        elif m.type == 'proxy':
            if (m.hgnc_id is not None and m.hgnc_id != self.hgnc_id) and (self.hgnc_id not in m._complex_hgnc_ids):
                raise ValueError(
                    'The macromolecule ' + m.id + ' does not appear to be associated with the gene ' + self.hgnc_id)

        if m.type == 'fragment_rna':
            if self.macromolecules['RNA']['lariat'] is not None:
                raise ValueError('Multiple ' + 'lariats' + ' assigned to ' + self.hgnc_id)
            self.macromolecules['RNA']['lariat'] = m.id
        elif m.type == 'premrna':
            if self.macromolecules['premrna'] is not None:
                raise ValueError('Multiple ' + 'premrna' + ' assigned to ' + self.hgnc_id)
            self.macromolecules['RNA']['premrna'] = m.id
        elif m.type == 'mrna':
            if m.coupling_coefficient is None:
                if self.macromolecules['RNA']['mrna']['other'] is not None:
                    raise ValueError('Multiple uncoupled' + 'mrna' + ' assigned to ' + self.hgnc_id)
                self.macromolecules['RNA']['mrna']['other'] = m.id
            else:
                if m.id in self.macromolecules['RNA']['mrna']['coupled']:
                    raise ValueError('Multiple reactions coupled to' + 'mrna' + ' for ' + self.hgnc_id)
                if len(self.macromolecules['RNA']['mrna']['coupled']) > 0:
                    raise ValueError('Multiple coupled' + 'mrna' + ' assigned to ' + self.hgnc_id)
                if list(m.coupling_coefficient.keys()) != ['mrna_formation']:
                    raise ValueError('Unexpected coupling type for mrna associated with ' + self.hgnc_id)

                cr = [r.id for r in m.reactions if m in r.coupled_metabolites]

                if len(cr) > 2:
                    raise ValueError(
                        'Unexpected mrna coupling to multiple reactions associated with ' + self.hgnc_id)
                self.macromolecules['RNA']['mrna']['coupled'][m.id] = cr
        elif m.type == 'protein':
            if m.coupling_coefficient is None:
                if not m.non_machinery:
                    self.macromolecules['Protein']['other'] += [m.id]
                else:
                    self.macromolecules['Protein']['non-machinery'] += [m.id]
            else:
                if list(m.coupling_coefficient.keys()) != ['catalysis']:
                    raise ValueError('Unexpected coupling type for monomer enzyme associated with ' + self.hgnc_id)
                if not m.enzyme:
                    raise ValueError('Incorrect tracking of enzymes')

                cr = [r.id for r in m.reactions if
                      m in r.coupled_metabolites and r.coupled_metabolites[m] == 'catalysis']
                self.macromolecules['Protein']['coupled'][m.id] = cr
        elif m.type == 'complex':
            if m.coupling_coefficient is None:
                self.macromolecules['Complex']['other'] += [m.id]
            else:
                if list(m.coupling_coefficient.keys()) != ['catalysis']:
                    raise ValueError('Unexpected coupling type for complex associated with ' + self.hgnc_id)
                if not m.enzyme:
                    raise ValueError('Incorrect tracking of enzymes')

                cr = [r.id for r in m.reactions if
                      m in r.coupled_metabolites and r.coupled_metabolites[m] == 'catalysis']
                self.macromolecules['Complex']['coupled'][m.id] = cr
        elif m.type == 'proxy':
            if m.coupling_coefficient is None or len(m.coupling_coefficient.keys()) != 1 or \
                    list(m.coupling_coefficient)[0] not in ['enzyme_degradation', 'mrna_degradation']:
                raise ValueError('Unexpected coupling type for proxy associated with ' + self.hgnc_id)
            proxy_type = list(m.coupling_coefficient)[0]

            if proxy_type == 'mrna_degradation':
                if len(self.macromolecules['Proxy']['mrna_degradation'].keys()) > 0:
                    raise ValueError('Multiple mrna degradation proxies associated with ' + self.hgnc_id)
                cr = [r.id for r in m.reactions if
                        m in r.coupled_metabolites and r.coupled_metabolites[m] == 'mrna_degradation']
                if len(cr) > 2:
                    raise ValueError('mRNA degradation proxy is coupled to more than 2 protein synthesis reactions')
                self.macromolecules['Proxy']['mrna_degradation'] = {m.id: cr}
            elif proxy_type == 'enzyme_degradation':
                cr = [r.id for r in m.reactions if
                      m in r.coupled_metabolites and r.coupled_metabolites[m] == 'enzyme_degradation']
                self.macromolecules['Proxy']['enzyme_degradation'][m.id] = cr
            else:
                raise ValueError('Unexpected proxy type asssociated with' + self.hgnc_id)
        else:
            raise ValueError('Unexpected metabolite type associated with ' + self.hgnc_id)

    def _check_macromolecules(self):
        if len(self.macromolecules['RNA']['mrna']['coupled']) == 0:
            raise ValueError('No coupled mrna molecule for ' + self.hgnc_id)
        if (len(self.macromolecules['Protein']['coupled']) + len(self.macromolecules['Protein']['other'])) == 0:
            raise ValueError('No protein metabolites associated with ' + self.hgnc_id)
        if (len(self.macromolecules['Protein']['coupled']) + len(
                self.macromolecules['Complex']['coupled'])) == 0 and not (
                self.ubiquitin_biogenesis or self._is_non_machinery_only or self._is_unmodeled_protein):
            raise ValueError('No coupled enzymes associated with ' + self.hgnc_id)
        if len(self.macromolecules['Proxy']['mrna_degradation']) == 0:
            raise ValueError('No mrna degradation proxy metabolite associated with ' + self.hgnc_id)

        rib_deg_coupling = self.ribosome_biogenesis or self.hgnc_id in ribosomal_genes
        if len(self.macromolecules['Proxy']['enzyme_degradation']) == 0 \
                and not (self.ubiquitin_biogenesis or self._enzyme_compartments == ['e']
                         or self._is_non_machinery_only or self._is_unmodeled_protein) and not rib_deg_coupling:  # TEMPORARY no rib deg coupling
            raise ValueError('No enzyme degradation proxy metabolite associated with ' + self.hgnc_id)

    def check(self):
        """Checks for completeness of self.reactions and self.macromolecules after adding all associated objects."""
        self._is_non_machinery_only = ((len(self.macromolecules['Protein']['coupled']) == 0) and len(
            self.macromolecules['Protein']['non-machinery']) > 0)
        self._enzyme_compartments = [e_id.split('_')[-1] for e_id in
                                     list(self.macromolecules['Protein']['coupled'].keys()) + list(
                                         self.macromolecules['Complex']['coupled'].keys())]
        self._check_reactions()
        self._check_macromolecules()
    
    def __repr__(self):

        if not self._summarized:
            n_mr = len(self.reactions['Catalysis_Reactions']['Metabolic_Module'])
            n_er = len(self.reactions['Catalysis_Reactions']['Expression_Module'])
            
            self._summ = self.hgnc_id + ' catalyzes '
            if n_mr > 0:
                self._summ += '{} metabolic reactions'.format(n_mr)
                if n_er > 0:
                    self._summ += ' and '
            if n_er > 0:
                self._summ += '{} expression reactions'.format(n_er)
        self._summarized = True

        return self._summ
