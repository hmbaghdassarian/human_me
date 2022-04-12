#!/usr/bin/env python
# coding: utf-8

import gc
import os
import random
import time
import warnings
from typing import Dict, List, Optional, Union

import cobra
import numpy as np
import pandas as pd
from faker import Faker
from tqdm import tqdm

pd.options.mode.chained_assignment = None

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from human_me.io import HiddenPrints
    from human_me import core
    from human_me.core.model import ME_Model
    from human_me.core.reaction import (ComplexDegradationReaction,
                                        ProteinDegradationReaction,
                                        to_metabolic_reaction)
    from human_me.preprocess import parse_complex
    from human_me.utils import functions as func
    from human_me.utils import machinery as mach
    from human_me.utils import parameters as params
    from human_me.utils.metabolites import MetaboliteBin
    from human_me.io import load_metabolic_model, load_psim

    with HiddenPrints():
        from human_me.build.build_utils import (
            get_all_expression_reactions, get_complex_df,
            get_expression_machinery, map_complex_machinery_compartment,
            map_machinery_compartment, merge_maps,
            parse_complex_degradation_reaction_id)
        from human_me.core import biomass
        from human_me.core.macromolecules.complex import (Complex,
                                                          RibosomalComplex)
        from human_me.core.macromolecules.protein import Protein
        from human_me.expression.build_ribosome_biogenesis_reactions import \
            build_ribosome, set_ribosomal_psim
        from human_me.expression.build_trna_expression_reactions import create_trna
        from human_me.expression.gene_expression.protein_expression import (
            degradation, ubiquitin)


class MEBuilder:
    """Builder class for ME Model."""
    def __init__(self, 
                m_model: Union[cobra.Model,str], 
                 psim_me: Union[pd.DataFrame, str],
                 model_id: str,
                 stochastic: bool = False, seed: int = 888, n_cores: int = os.cpu_count(),
                 non_machinery: Optional[Dict[str, List[str]]] = None, knock_out: Optional[List[str]] = None,
                 dummy_protein: bool = True, context_specific_dummy: bool = False,
                 minimal_proteome: bool = True, compress_mrna: bool = True,
                 check_all: bool = True,
                 deg_args: Dict[str, bool] = {'couple': True, 'reversible_complex_formation': False, 'nonenzyme_degradation': False,
                                              'complex_degradation': True}
                 ):
        """See build_me function for parameter descriptions."""

        # check deg args
        if not deg_args['complex_degradation']:
            if not deg_args['nonenzyme_degradation'] or not deg_args['reversible_complex_formation']:
                err = 'If complex degradation is not included, the formation must be reversible and the individual '
                err += 'components must be able to be degraded. Set nonenzyme_degradation and reversible_complex_formation'
                err += ' to True or complex_degradation to True'
                raise ValueError(err)
        if deg_args['couple'] and not deg_args['complex_degradation']:
            raise ValueError(
                'In order to couple metabolic catalysis to enzyme degradation, complex_degradation must be True')
        self.deg_args = deg_args

        if non_machinery is None:
            non_machinery = dict()
        for exception in ['HGNC:12468', 'HGNC:12463'] + mach.rbps:
            if exception in non_machinery:
                warnings.warn(
                    exception + ' is a ribosomal or ubiquitin-related gene, which cannot be specified as non-machinery, removing')
                del non_machinery[exception]
        self.non_machinery = non_machinery

        if knock_out is None:
            self.knock_out = list()
        else:
            self.knock_out = knock_out
        if len(set(self.knock_out).intersection(mach.expression_machinery)) > 0:
            raise ValueError(
                'Knock outs can only be applied to metabolic machinery and non-machinery, not expression machinery')
        if len(set(self.knock_out).intersection(self.non_machinery)) > 0:
            raise ValueError('Speficied knocking out of genes that are also specified to be expressed as non-machinery')

        # all parameters that use psim_me as input
        self.psim_me = load_psim(psim_me)
        psim_rib = set_ribosomal_psim(self.psim_me)

        # all parameters that use m_model as input
        self.m_model = load_metabolic_model(m_model)
        self.metabolic_machinery, self.all_machinery = mach.get_model_machinery(self.m_model)
        self.model_metabolites = MetaboliteBin(self.m_model)
        self.biomass_reactions= biomass.create_biomass_reactions(self.model_metabolites)
        self.trna_biogenesis_reactions, self.charged_trna_map, self.modified_trna_transcript_c = create_trna(self.model_metabolites)


        self.n_cores = n_cores
        if self.n_cores in [0, 1, None]:
            self._par = False
        else:
            self._par = True

        self.stochastic = stochastic
        if not self.stochastic and seed is None:
            raise ValueError('For non-stochastic outputs, please assign a seed')
        self.seed = seed
        random.seed(self.seed)
        self._seeds = random.sample(range(0, int((2 ** 32 - 1))), k=250000)  # 10x # of human protein-coding genes
        # assign a gene-specific seed, independent of the order the gene's expression reactions
        all_genes = self.all_machinery + list(self.non_machinery)
        self._gene_seed_map = dict(zip(all_genes, self._seeds[:len(all_genes)]))
        # if _gene_seed_map causes issues, it is unecessary and can be removed bc sorted() was added to
        # each gene for which expression reactions are being built. _gene_seed_map makes it independent of order,
        # but ordering alone suffices (will have to reset self._seed_idx to 0 and +=1 in each iteration)
        self._seed_idx = len(all_genes)
        # self._seed_idx = 0

        # get pre-generated reactions - the compress_mrna arg requires that they be run with that input
        self.compress_mrna = compress_mrna
        print('Generate ubiquitin reactions for proteasomal degradation')
        self.ub_args = ubiquitin.express_ubiquitin(model_metabolites=self.model_metabolites, psim_me = self.psim_me, compress_mrna=self.compress_mrna, 
                                                    charged_trna_map=self.charged_trna_map, modified_trna_transcript_c=self.modified_trna_transcript_c)
        print('Generate ribosome')

        # ribosome seeds different from seeds here for gene_info bc it also calls random.sample
        random.seed(self.seed)
        rib_seed = random.randint(0, int((2 ** 32 - 1)))
        ribosomal_reactions, self.ribosome_complex_c = build_ribosome(self.model_metabolites, psim_rib, self.modified_trna_transcript_c, self.charged_trna_map, 
                                                                        self.ub_args, self.compress_mrna,
                                                                        self.deg_args['reversible_complex_formation'],
                                                                        stochastic=self.stochastic, seed=rib_seed)

        self.dummy_protein = dummy_protein
        self.context_specific_dummy = context_specific_dummy

        self.deorphaned = None
        self.orphan = None

        self.me_reactions = self.trna_biogenesis_reactions + ribosomal_reactions + self.ub_args['ub_reactions']
        # map HGNC ID to a dictionary of compartments and cobra.Metabolite proteins
        self.id_protein_map = dict()
        self.complex_id_metabolite_map = dict()  # map complex id to the complex cobra.Metabolite

        self._ko_id_protein_map = dict()
        self._ko_complex_id_metabolite_map = dict()

        self.id_reactions_map = dict()
        self.complex_reactions_map = dict()
        self.check_all = check_all

        self.model_id = model_id
        self.minimal_proteome = minimal_proteome

    def _clean_non_machinery(self):
        ipm = self.id_protein_map.copy()
        for k, v in self._ko_id_protein_map.items():
            ipm[k] = v
        nm2 = self.non_machinery.copy()
        for hgnc_id in nm2:
            for compartment in self.non_machinery[hgnc_id]:
                if not ipm[hgnc_id][compartment].non_machinery:  # when overlaps with enzyme
                    self.non_machinery[hgnc_id].remove(compartment)
            if len(self.non_machinery[hgnc_id]) == 0:
                del self.non_machinery[hgnc_id]
        del nm2, ipm

    def express_metabolic_enzymes(self):
        """Get protein expression reactions for all metabolic enzymes and user-input non-machinery"""

        # get protein expression for all metabolic reactions
        print('Generate protein expression reactions for metabolic enzymes and non-machinery')

        self.loop_machinery = list(set(self.metabolic_machinery + list(self.non_machinery)))

        gene_reaction_map = func.create_gene_reaction_map(self.m_model.reactions)
        for hgnc_id in self.non_machinery:
            if hgnc_id not in gene_reaction_map:
                gene_reaction_map[hgnc_id] = None

        iterable = sorted(set(self.loop_machinery).difference(
            self.knock_out))  # sorted so seeds are consistent; probably unecessary given the gene seed map but keeping
        #         if not self._par:
        for hgnc_id in tqdm(iterable):
            nml = list()
            if hgnc_id in self.non_machinery:
                nml = self.non_machinery[hgnc_id]
            expr_reactions, protein_metabolites = get_all_expression_reactions(model_metabolites = self.model_metabolites, 
                                                                                hgnc_id=hgnc_id,
                                                                               reactions=gene_reaction_map[hgnc_id],
                                                                               compress_mrna=self.compress_mrna,
                                                                               ub_args=self.ub_args,
                                                                               psim=self.psim_me,
                                                                               machinery_list=self.metabolic_machinery, 
                                                                               modified_trna_transcript_c=self.modified_trna_transcript_c, 
                                                                               charged_trna_map=self.charged_trna_map,
                                                                               nonmachinery_locations=nml,
                                                                               stochastic=self.stochastic,
                                                                               seed=self._gene_seed_map[
                                                                                   hgnc_id])  # self._seeds[self._seed_idx])
            self.id_protein_map[hgnc_id] = {p.compartment: p for p in
                                            protein_metabolites}  # store compartments and metabolite objects for each gene
            self.id_reactions_map[hgnc_id] = expr_reactions
            self.me_reactions += expr_reactions
        #             self._seed_idx += 1
        #        else: # currently doesn't work bc makes a copy of objects, so does not retain assigned r.enzyme_compartment attr
        #             pool = multiprocessing.Pool(processes = self.n_cores)
        #             try:
        #                 n_iter = len(iterable)
        #                 args = zip(iterable, [self.m_model]*n_iter, [gene_reaction_map]*n_iter, [self.psim_me]*n_iter, [self.metabolic_machinery]*n_iter, 
        #                                       [self.ub_args]*n_iter, [self.modified_transcript_c]*n_iter, [self.charged_trna_map]*n_iter, [self.compress_mrna]*n_iter,
        #                            [self.non_machinery]*n_iter, [self.stochastic]*n_iter, list(range(self._seed_idx, len(iterable))))
        #                 mm = pool.starmap(emm_par, args)
        #                 self._seed_idx += len(iterable)
        #                 pool.close()
        #                 pool.join()
        #                 gc.collect()
        #             except:
        #                 pool.close()
        #                 pool.join()
        #                 gc.collect()
        #                 raise ValueError('Parallelization failed')
        #             self.id_protein_map = dict(zip(iterable, [i[0] for i in mm]))
        #             expr_reactions = [i[1] for i in mm]
        #             self.id_reactions_map = dict(zip(iterable, expr_reactions))
        #             self.me_reactions += func.flatten_list(expr_reactions)
        #             del expr_reactions

        for hgnc_id in sorted(self.knock_out):
            # None bc will add later for expression model specific to this
            expr_reactions, protein_metabolites = get_all_expression_reactions(model_metabolites = self.model_metabolites, 
                                                                                hgnc_id=hgnc_id,
                                                                               reactions=gene_reaction_map[hgnc_id],
                                                                               psim=self.psim_me,
                                                                               machinery_list=self.metabolic_machinery, 
                                                                               modified_trna_transcript_c=self.modified_trna_transcript_c, 
                                                                               charged_trna_map=self.charged_trna_map,
                                                                               compress_mrna=self.compress_mrna,
                                                                               ub_args=self.ub_args,
                                                                               stochastic=self.stochastic,
                                                                               seed=self._gene_seed_map[
                                                                                   hgnc_id])  # self._seeds[self._seed_idx])
            #             self._seed_idx += 1
            self._ko_id_protein_map[hgnc_id] = {p.compartment: p for p in
                                                protein_metabolites}  # store compartments and metabolite objects for each gene

    def express_expression_enzymes(self):
        """Get protein expression reactions for all expression machinery"""
        # This method continues to add any expression module machinery that may have arisen from adding expression
        # reactions for expression machinery.

        self._expr_rxn_cmap = dict()
        for r in self.me_reactions:
            if hasattr(r, 'enzyme_compartment'):
                generic_id = func.parse_me_reaction_id(r.id)
                if generic_id in self._expr_rxn_cmap:
                    raise ValueError('Unexpected presence of expression enzyme')
                self._expr_rxn_cmap[generic_id] = {'compartment': r.enzyme_compartment, 'seed': r._compartment_seed}

        gene_reaction_map, expression_machinery_me = get_expression_machinery(self.me_reactions)

        for hgnc_id in tqdm(
                sorted(set(expression_machinery_me).difference(self.knock_out))):  # sorted so seeds are in same order
            nml = list()
            if hgnc_id in self.non_machinery:
                nml = self.non_machinery[hgnc_id]

            # ensure same compartment for existing reactions
            not_present = dict()
            for r in gene_reaction_map[hgnc_id]:
                generic_id = func.parse_me_reaction_id(r.id)
                if generic_id in self._expr_rxn_cmap:
                    r.enzyme_compartment, r._compartment_seed = self._expr_rxn_cmap[generic_id]['compartment'], \
                        self._expr_rxn_cmap[generic_id]['seed']
                else:
                    not_present[generic_id] = r  # if shows up multiple times, will map to same compartment anyways

            expr_reactions, protein_metabolites = get_all_expression_reactions(model_metabolites = self.model_metabolites,
                                                                                hgnc_id=hgnc_id, 
                                                                                psim=self.psim_me,
                                                                               machinery_list=expression_machinery_me,
                                                                               modified_trna_transcript_c=self.modified_trna_transcript_c, 
                                                                               charged_trna_map=self.charged_trna_map,
                                                                               reactions=gene_reaction_map[hgnc_id],
                                                                               compress_mrna=self.compress_mrna,
                                                                               ub_args=self.ub_args,
                                                                               nonmachinery_locations=nml,
                                                                               stochastic=self.stochastic,
                                                                               seed=self._gene_seed_map[
                                                                                   hgnc_id])  # self._seeds[self._seed_idx])

            # track new reaction compartments
            for generic_id, r in not_present.items():
                self._expr_rxn_cmap[generic_id] = {'compartment': r.enzyme_compartment, 'seed': r._compartment_seed}

                # ensure same compartment for existing reactions
            for r in expr_reactions:
                generic_id = func.parse_me_reaction_id(r.id)
                if generic_id in self._expr_rxn_cmap:
                    rmap = self._expr_rxn_cmap[generic_id]
                    r.enzyme_compartment, r._compartment_seed = rmap['compartment'], rmap['seed']

            #             self._seed_idx += 1
            if hgnc_id not in set(expression_machinery_me).intersection(self.loop_machinery):
                if hgnc_id in self.id_protein_map:
                    raise ValueError('Some genes not accounted for when generating metabolic machinery expression reactions')
                self.id_protein_map[hgnc_id] = {p.compartment: p for p in protein_metabolites}
            # when there is machinery overlap between metabolic and expression module, deal with compartment overlap
            else:
                ids_to_keep = list(set([r.id for r in expr_reactions]).difference([r.id for r in self.me_reactions]))
                expr_reactions = [r for r in expr_reactions if r.id in ids_to_keep]

                temp_map = {p.compartment: p for p in protein_metabolites}
                for comp, met in temp_map.items():
                    if comp not in self.id_protein_map[hgnc_id]:
                        self.id_protein_map[hgnc_id][comp] = met
                    elif not met.non_machinery:  # in the case that it was specified as non-machinery during metabolic iteration
                        met = self.id_protein_map[hgnc_id][comp]
                        met.non_machinery = False

            if hgnc_id not in self.id_reactions_map:
                self.id_reactions_map[hgnc_id] = expr_reactions
            else:
                self.id_reactions_map[hgnc_id] += expr_reactions

            self.me_reactions += expr_reactions

        gene_reaction_map_2, expression_machinery_me_2 = get_expression_machinery(self.me_reactions)
        new_expression_machinery = list(
            set(expression_machinery_me_2).difference(expression_machinery_me + self.knock_out))

        counter = 1
        while len(
                new_expression_machinery) > 0:  # this condition leaves possibility that an existing machinery but with a new compartment is added and not accounted for
            print('No. iterations for new expression machinery: {}'.format(counter))
            for hgnc_id in tqdm(sorted(set(expression_machinery_me_2))):  # sorted so seeds are in same order
                nml = list()
                if hgnc_id in self.non_machinery:
                    nml = self.non_machinery[hgnc_id]

                # ensure same compartment for existing reactions
                not_present = dict()
                for r in gene_reaction_map_2[hgnc_id]:
                    generic_id = func.parse_me_reaction_id(r.id)
                    if generic_id in self._expr_rxn_cmap:
                        r.enzyme_compartment, r._compartment_seed = self._expr_rxn_cmap[generic_id]['compartment'], \
                            self._expr_rxn_cmap[generic_id]['seed']
                    else:
                        not_present[generic_id] = r  # if shows up multiple times, will map to same compartment anyways

                expr_reactions, protein_metabolites = get_all_expression_reactions( model_metabolites = self.model_metabolites, 
                                                                                hgnc_id=hgnc_id, psim=self.psim_me,
                                                                                   machinery_list=expression_machinery_me_2,
                                                                               modified_trna_transcript_c=self.modified_trna_transcript_c, 
                                                                               charged_trna_map=self.charged_trna_map,
                                                                                   reactions=gene_reaction_map_2[
                                                                                       hgnc_id],
                                                                                   compress_mrna=self.compress_mrna,
                                                                                   ub_args=self.ub_args,
                                                                                   nonmachinery_locations=nml,
                                                                                   stochastic=self.stochastic,
                                                                                   seed=self._gene_seed_map[
                                                                                       hgnc_id])  # self._seeds[self._seed_idx])

                # track new reaction compartments
                for generic_id, r in not_present.items():
                    self._expr_rxn_cmap[generic_id] = {'compartment': r.enzyme_compartment, 'seed': r._compartment_seed}

                    # ensure same compartment for existing reactions
                for r in expr_reactions:
                    generic_id = func.parse_me_reaction_id(r.id)
                    if generic_id in self._expr_rxn_cmap:
                        rmap = self._expr_rxn_cmap[generic_id]
                        r.enzyme_compartment, r._compartment_seed = rmap['compartment'], rmap['seed']

                        # self._seed_idx += 1
                if hgnc_id not in set(expression_machinery_me_2).intersection(expression_machinery_me + self.metabolic_machinery):
                    if hgnc_id in self.id_protein_map:
                        raise ValueError('Some genes not accounted for when generating metabolic machinery expression reactions')
                    self.id_protein_map[hgnc_id] = {p.compartment: p for p in protein_metabolites}
                # when there is machinery overlap between metabolic and expression module, deal with compartment overlap
                else:
                    ids_to_keep = list(
                        set([r.id for r in expr_reactions]).difference([r.id for r in self.me_reactions]))
                    expr_reactions = [r for r in expr_reactions if r.id in ids_to_keep]

                    temp_map = {p.compartment: p for p in protein_metabolites}
                    for comp, met in temp_map.items():
                        if comp not in self.id_protein_map[hgnc_id]:
                            self.id_protein_map[hgnc_id][comp] = met
                        elif not met.non_machinery:  # in the case that it was specified as non-machinery during metabolic iteration
                            met = self.id_protein_map[hgnc_id][comp]
                            met.non_machinery = False

                if hgnc_id not in self.id_reactions_map:
                    self.id_reactions_map[hgnc_id] = expr_reactions
                else:
                    self.id_reactions_map[hgnc_id] += expr_reactions

                self.me_reactions += expr_reactions

            # get protein expression reactions for all expression module reactions
            expression_machinery_me = expression_machinery_me_2
            gene_reaction_map_2, expression_machinery_me_2 = get_expression_machinery(self.me_reactions)
            new_expression_machinery = list(
                set(expression_machinery_me_2).difference(expression_machinery_me + self.knock_out))
            counter += 1

        for r in self.me_reactions:
            if len(r.genes) > 0 and not hasattr(r, 'enzyme_compartment'):
                generic_id = func.parse_me_reaction_id(r.id)
                if generic_id not in self._expr_rxn_cmap:
                    raise ValueError('Reactions with unaccounted for compartments')
                r.enzyme_compartment, r._compartment_seed = self._expr_rxn_cmap[generic_id]['compartment'], self._expr_rxn_cmap[generic_id]['seed']

        self._clean_non_machinery()

    def express_dummy_protein(self):
        """Generate the dummy protein."""
        if self.dummy_protein:
            print('Express dummy protein')
            dummy_psim = func.average_protein_features(psim_me=self.psim_me,
                                                       context_specific=self.context_specific_dummy, 
                                                       metabolic_machinery=self.metabolic_machinery)

            dummy_reactions, dm = get_all_expression_reactions( model_metabolites = self.model_metabolites, 
                                                                hgnc_id='HGNC:DUMMY', psim=dummy_psim, machinery_list=[],
                                                                modified_trna_transcript_c=self.modified_trna_transcript_c, 
                                                                charged_trna_map=self.charged_trna_map,
                                                               reactions=None, compress_mrna=self.compress_mrna,
                                                               ub_args=self.ub_args, nonmachinery_locations=['c'],
                                                               stochastic=self.stochastic,
                                                               seed=self._seeds[self._seed_idx])
            self._seed_idx += 1
            dm[0].non_machinery = False
            for r in dummy_reactions:
                for m in r.metabolites:
                    if isinstance(m, Protein) and m.id.startswith(
                            'HGNC:DUMMY'):  # str requirement to avoid converting ub proteins
                        m.dummy = True
                if len(r.genes) > 0:
                    generic_id = func.parse_me_reaction_id(r.id)
                    if generic_id not in self._expr_rxn_cmap:
                        raise ValueError('Reactions with unaccounted for compartments')
                    r.enzyme_compartment, r._compartment_seed = self._expr_rxn_cmap[generic_id]['compartment'], self._expr_rxn_cmap[generic_id]['seed']

            self.dummy_protein = {'protein_metabolite': dm[0], 'dummy_expression_reactions': dummy_reactions}

            srs = [sr for sr in list(self.dummy_protein['protein_metabolite'].reactions) if
                   self.dummy_protein['protein_metabolite'] in sr.products and not isinstance(sr, ProteinDegradationReaction)]
            if len(srs) != 1:
                raise ValueError(self.dummy_protein[
                    'protein_metabolite'].id + ' has an incorrect number of associated synthesis reactions')
            srs[0].synthesis, srs[0].synthesis_type = True, 'protein'

            self.me_reactions += self.dummy_protein['dummy_expression_reactions']

        else:
            self.dummy_protein = None

    def get_complex_info(self):
        "Parse and organize GPRs"
        # ------------Metabolic Complexes
        print('Get metabolic module complex information')
        complex_df = get_complex_df(reactions=[r for r in self.m_model.reactions if len(r.genes) > 0],
                                    knock_out=self.knock_out)
        complex_df['category'] = 'metabolic_reaction'

        # ------------Expression Complexes
        print('Get expression module complex information')
        me_complex_df = get_complex_df(reactions=[r for r in self.me_reactions if len(r.genes) > 0],
                                       knock_out=self.knock_out)
        # deal with most expression reactions having redundant machinery in a concise manner:
        me_complex_df.reaction_id = me_complex_df.reaction_id.apply(lambda x: func.parse_me_reaction_id(x))
        me_complex_df.drop_duplicates(keep='first', inplace=True)  # if all but reaction id HGNC were the same
        me_complex_df.reset_index(inplace=True, drop=True)
        me_complex_df['category'] = 'expression_reaction'

        # -------------Merge Modules
        # merge bc will deal with duplicate complexes, in case there is duplicates b/w metabolic and expression module
        complex_df = pd.concat([complex_df, me_complex_df], axis=0)
        complex_df.reset_index(inplace=True, drop=True)

        print('Assign unique complex ids for unique machinery-compartment sets across all reactions')
        # assign complex ids for reactions that have complexes in them
        complex_df['complex_id'] = float('nan')
        complex_df.loc[complex_df[complex_df.is_complex].index, 'complex_id'] = complex_df.loc[
            complex_df[complex_df.is_complex].index, 'reaction_id']

        # if a reaction generates multiple complexes, make sure each complex has a unique ID
        crm_ = complex_df[complex_df.creates_multiple_reactions & (complex_df.is_complex)].reaction_id.unique()

        for crm in crm_:
            df = complex_df[(complex_df.reaction_id == crm) & (complex_df.is_complex)]
            if df.shape[0] > 1:  # reaction creates multiple complexes
                counter = 0
                for i in df.index:
                    complex_df.loc[i, 'complex_id'] = complex_df.loc[i, 'complex_id'] + '_' + str(counter)
                    counter += 1

        # if complexes are duplicated across different reactions assigned to the same compartment,
        # generate a singular unique id
        dup_complexes = complex_df[complex_df.is_complex].duplicated(subset=['compartment', 'machinery'])
        dup_complexes = complex_df.loc[dup_complexes.index[np.where(dup_complexes)]]
        dup_idx = dup_complexes.index
        dup_complexes = dup_complexes.drop_duplicates(subset=['compartment', 'machinery'], keep='first')

        dup_complex_map = dict()

        if self.seed is not None:
            dup_seeds = self._seeds[self._seed_idx:self._seed_idx + dup_complexes.shape[0]]
            self._seed_idx += dup_complexes.shape[0]
        else:
            dup_seeds = [None] * dup_complexes.shape[0]

        dup_complexes.reset_index(inplace=True, drop=True)
        for i in dup_complexes.index:
            Faker.seed(dup_seeds[i])
            f1 = Faker()
            dup_complex_map[(dup_complexes.loc[i, 'compartment'], dup_complexes.loc[i, 'machinery'])] = f1.uuid4().split('-')[0]
            self._seed_idx += 1

        complex_df.loc[dup_idx, 'complex_id'] = complex_df.loc[dup_idx, ][['compartment', 'machinery']].apply(
            lambda x: dup_complex_map[(x[0], x[1])], axis=1).tolist()

        self.complex_df = complex_df

    def generate_complex_reactions(self):
        """Create complex formation and degradation reactions."""
        # create a mapping of the unique self.complex_df ids to the actual complex metabolite
        unique_complexes = self.complex_df[self.complex_df.is_complex]
        unique_complexes = unique_complexes.drop_duplicates(subset='complex_id', keep='first')
        unique_complexes.reset_index(inplace=True, drop=True)

        self.complex_formation_reactions = list()  # store all complex formation reactions
        complex_degradation_reactions = list()

        retain = list(set(func.flatten_list(
            [i.split(';') for i in unique_complexes[~unique_complexes.knock_out.astype(bool)].machinery.tolist()])))
        self.additional_ko = list()

        counter = 0
        for i in tqdm(unique_complexes.index):
            ko = unique_complexes.loc[i, 'knock_out']
            if ko:
                # get the  genes that are only expressed to participate as part of a complex that is knocked out
                self.additional_ko += list(
                    set(unique_complexes.loc[i, 'machinery'].split(';')).difference(retain + self.knock_out))
            complex_id = unique_complexes.loc[i, 'complex_id']
            compartment = unique_complexes.loc[i, 'compartment']
            machinery = unique_complexes.loc[i, 'machinery'].split(';')
            machinery_metabolites = list()
            if not ko:
                counter_rib = 0
                for m in machinery:
                    if m != 'ribosome':
                        machinery_metabolites.append(self.id_protein_map[m][compartment])
                    else:
                        machinery_metabolites.append(self.ribosome_complex_c)
                        counter_rib += 1

                if counter_rib == 0:
                    complex_metabolite = Complex(metabolites=machinery_metabolites, complex_id=complex_id)
                else:
                    complex_metabolite = RibosomalComplex(metabolites=machinery_metabolites, complex_id=complex_id)
                if len(complex_id) > 247:  # ids that are too long
                    complex_metabolite.update_id(new_id=str(counter))  # complex_metabolite.udate_id()
                    counter += 1

                    self.complex_df.complex_id.replace(to_replace=complex_id, value=complex_metabolite.temp_id,
                                                       inplace=True)

                complex_formation_reaction = complex_metabolite.form_complex(
                    reversible=self.deg_args['reversible_complex_formation'],
                    synthesis=True, synthesis_type='complex')
                complex_formation_reaction.synthesis = True  # for ribosomal complex formation
                if self.deg_args['complex_degradation']:
                    if complex_metabolite.compartment in ['c', 'n', 'r', 'g', 'pm']:
                        complex_degradation_reaction = degradation.degrade(complex_metabolite, model_metabolites=self.model_metabolites, 
                                                                           **{'ub_args': self.ub_args})
                    elif complex_metabolite.compartment == 'e':
                        complex_degradation_reaction = list()
                    else:
                        complex_degradation_reaction = degradation.degrade(complex_metabolite, model_metabolites=self.model_metabolites, )

                    for r in complex_degradation_reaction:
                        if len(r.genes) > 0:
                            generic_id = r.id.split(complex_metabolite._deg_id)[1][1:]
                            if generic_id not in self._expr_rxn_cmap:
                                raise ValueError('Reactions with unaccounted for compartments')
                            r.enzyme_compartment, r._compartment_seed = self._expr_rxn_cmap[generic_id]['compartment'], \
                                                                        self._expr_rxn_cmap[generic_id]['seed']

                    complex_degradation_reactions += complex_degradation_reaction
                else:
                    complex_degradation_reaction = list()

                self.complex_formation_reactions.append(complex_formation_reaction)
                self.complex_id_metabolite_map[complex_metabolite.temp_id] = complex_metabolite
                self.complex_reactions_map[complex_metabolite.temp_id] = [r.id for r in [
                    complex_formation_reaction] + complex_degradation_reaction]
            else:
                for m in machinery:
                    if m in self.id_protein_map:
                        machinery_metabolites.append(self.id_protein_map[m][compartment])
                    else:
                        machinery_metabolites.append(self._ko_id_protein_map[m][compartment])
                complex_metabolite = Complex(metabolites=machinery_metabolites, complex_id=complex_id)
                self._ko_complex_id_metabolite_map[complex_metabolite.temp_id] = complex_metabolite

        if self.deg_args['complex_degradation'] and self.check_all:
            # --check that machinery is the same for complex degradation (with exception of proteasomal degradation of ribosome)
            # as protein degradation

            # double check that only ribosomal degradation has additional machinery
            # this following code can be commented out if don't want to double check
            # keep, so in future iterations, can be used to iteratively add new machinery

            # this is all pretty hard-coded, starting from human_me.core.reaction.ComplexDegradationReaction._set_proteasomal_degration()

            reactions_to_add = []
            self.complex_degradation_reactions = complex_degradation_reactions
            for r in self.complex_degradation_reactions:
                if len(r.genes) > 0:
                    rxn_mach = parse_complex.eval_complex(r.gene_reaction_rule)
                    for rm in rxn_mach:
                        if type(rm) != list:
                            rm = [rm]
                        else:
                            rm = sorted(rm)
                        present = self.complex_df[(self.complex_df.machinery == ';'.join(rm)) & (
                            self.complex_df.compartment == r.enzyme_compartment) & (
                            self.complex_df.reaction_id == parse_complex_degradation_reaction_id(r.id))]
                        if present.shape[0] == 0:
                            reactions_to_add.append(r)
            # Option 1: degrade rRNA with ribosomal degradation
            if len(reactions_to_add) != 2 or not np.all([r._ribosomal_degradation for r in reactions_to_add]):
                err = 'Internal: Expected proteasomal degradation of ribosomal complexes to be the only difference in '
                err += 'machinery of complex degradation reactions. If other missed ones (perhaps in very small model '
                err += 'scenarios, but seems unlikely), will have to account for iteratitively adding new machinery '
                err += 'degradation reactions'
                print(err)
                raise ValueError('See internal error message above')
        #      # Option 2: degrade proteins with ribosomal degradation, releasing rRNA as intact
        #     if  len(reactions_to_add) != 0:
        #         err = 'Internal: Expected no additional machinery'
        #         print(err)
        #         raise ValueError('See internal error message above')

        # if want, in the future, can back-track and remove associated expression reactions with these
        self.additional_ko = list(set(self.additional_ko))
        # #started some code for this:
        # for hgnc_id in additional_ko:
        #     if len(self.id_protein_map[hgnc_id])>1:
        #         err = 'Internal: Have not accounted for scenario where knocked out complex has an associated gene'
        #         err = 'that was not explicitly knocked out but '
        #         raise ValueError(err)
        #     else:
        #         self._ko_protein_map[hgnc_id] = self.id_protein_map.pop(hgnc_id)
        #         # remove expression reactions...

    def get_keff(self):
        """Estimate the keff of all enzymes."""
        # calculated beforeprotein minimization to get average of all proteins
        # get SASA and keff values for coupling
        print('Calculate enzyme k_effs')
        # retain knocked-out genes in keff calculations
        # do not include non_machinery
        cplx_bool_map = {True: 'complex', False: 'monomer'}
        ko_bool_map = {True: 'ko', False: 'retain'}
        mach_col_map = {True: 'complex_id', False: 'machinery'}
        ko_map = {'complex': {'ko': self._ko_complex_id_metabolite_map, 'retain': self.complex_id_metabolite_map},
                  'monomer': {'ko': self._ko_id_protein_map, 'retain': self.id_protein_map}}

        self.complex_df['MW_kDa'] = float('nan')
        for i in tqdm(self.complex_df.index):
            cplx_bool, ko_bool, _compartment = self.complex_df.loc[i, ['is_complex', 'knock_out', 'compartment']]
            _mach = self.complex_df.loc[i, mach_col_map[cplx_bool]]

            if not cplx_bool:
                enzyme_to_couple = ko_map[cplx_bool_map[cplx_bool]][ko_bool_map[ko_bool]][_mach][_compartment]
            else:
                enzyme_to_couple = ko_map[cplx_bool_map[cplx_bool]][ko_bool_map[ko_bool]][_mach]
            self.complex_df.loc[i, 'MW_kDa'] = enzyme_to_couple.formula_weight / 1000

        self.complex_df['SASA'] = self.complex_df.MW_kDa.apply(lambda x: func.SASA(x))
        median_SASA = self.complex_df.SASA.median()
        self.complex_df['keff'] = self.complex_df['SASA'].apply(lambda x: x * (params.KEFF_MEDIAN / median_SASA))

        if self.dummy_protein is not None:
            self.dummy_protein['protein_metabolite'].keff = func.SASA(
                self.dummy_protein['protein_metabolite'].formula_weight / 1000) * (params.KEFF_MEDIAN / median_SASA)

    def minimize_proteome(self):
        """In the presence of OR GPRs, retain only one reaction, that catalyzed by the enzyme with the lowest MW"""
        if self.minimal_proteome:
            c_og = self.complex_df.copy()
            n_reactions_og = len(self.me_reactions) + len(self.complex_formation_reactions)
            if self.deg_args['complex_degradation']:
                n_reactions_og += len(self.complex_degradation_reactions)

            drop_index = list()
            reaction_multiple = self.complex_df[self.complex_df.creates_multiple_reactions].reaction_id.unique().tolist()
            for rm in reaction_multiple:
                df = self.complex_df[self.complex_df.reaction_id == rm]
                # deal w/ knockouts
                to_drop = list()

                check = False
                if df[df.knock_out].shape[0] > 0 and df[df.knock_out].shape[0] < df.shape[0]:
                    check = True
                    to_drop += df[df.knock_out.astype(bool)].index.tolist()
                    df_ = df[~df.knock_out.astype(bool)]
                else:
                    df_ = df
                # don't directly drop machinery in case they are used in multiple reactions and are minimal in
                # another one of those reactions
                to_drop += df_[df_.MW_kDa != df_.MW_kDa.min()].index.tolist()
                if df.shape[0] - len(to_drop) == 1:
                    drop_index += to_drop
                elif df.machinery.unique().shape[0] == df.shape[0]:
                    if check:
                        raise ValueError(
                            'Make sure proceeding line of code is correct for knock-out situation, currently has not been tested')
                    # rare case where two different complexes have the same MW
                    pop_ = df_.index.tolist()
                    drop_index += random.sample(population=pop_, k=len(pop_) - 1)
                else:
                    raise ValueError('Something went wrong in selecting a complex by lowest molecular weight')

            self.complex_df.drop(index=drop_index, inplace=True)

            # COMPLEXES--------------------------------------------------------------------
            # get rid of dropped complexes and reactions associated with dropped complexes
            complexes_to_drop = sorted(
                set(c_og[c_og.is_complex & ~c_og.knock_out.astype(bool)].complex_id).difference(self.complex_df.complex_id))
            complexes_to_drop_id = func.flatten_list([self.complex_reactions_map[c_id] for c_id in complexes_to_drop])
            self.complex_formation_reactions = [r for r in self.complex_formation_reactions if
                                                r.id not in complexes_to_drop_id]
            if self.deg_args['complex_degradation']:
                self.complex_degradation_reactions = [r for r in self.complex_degradation_reactions if
                                                      r.id not in complexes_to_drop_id]
            for c_id in complexes_to_drop:
                del self.complex_reactions_map[c_id]
                del self.complex_id_metabolite_map[c_id]

            # MONOMERS--------------------------------------------------------------------
            # find which active monomers were dropped, remembering to account for compartment specific expression
            # and retained active complexes that may not have active monomers
            # non-machinery is not dropped because it is not on complex_df

            # monomeric enzymes that are retained or dropped
            all_monomer_enzymes = map_machinery_compartment(
                c_og[~c_og.is_complex.astype(bool) & ~c_og.knock_out.astype(bool)])
            retain_monomer_enzymes = map_machinery_compartment(
                self.complex_df[~self.complex_df.is_complex.astype(bool) & ~self.complex_df.knock_out.astype(bool)])

            # complex monomeric components that are retained or dropped
            all_monomer_complexes = map_complex_machinery_compartment(c_og)
            retain_monomer_complexes = map_complex_machinery_compartment(self.complex_df)

            # all monomers vs the monomers that are retained, accounting for both monomeric enzymes and complexes
            all_monomers = merge_maps(all_monomer_enzymes, all_monomer_complexes)
            retain_monomers = merge_maps(retain_monomer_enzymes, retain_monomer_complexes)

            drop_monomers = dict()
            for hgnc_id, compartments in all_monomers.items():
                if hgnc_id not in retain_monomers:
                    drop_monomers[hgnc_id] = compartments
                else:
                    dropped_compartments = set(compartments).difference(retain_monomers[hgnc_id])
                    if len(dropped_compartments) > 0:
                        drop_monomers[hgnc_id] = dropped_compartments

            # drop monomers and associated reactions
            reactions_to_remove = list()
            for hgnc_id, compartments in drop_monomers.items():  # keeping compartments in mind, remove dropped monomeric proteins and reactions
                # proteins
                new_prot = {c: v for c, v in self.id_protein_map[hgnc_id].items() if c not in compartments}
                if len(new_prot) > 0:
                    self.id_protein_map[hgnc_id] = new_prot
                else:
                    del self.id_protein_map[hgnc_id]

                # reactions
                rtk = [r for r in self.id_reactions_map[hgnc_id] if
                       hasattr(r, '_final_compartments') and len(set(r._final_compartments).difference(compartments)) > 0]
                if len(rtk) > 0:
                    rtk += [r for r in self.id_reactions_map[hgnc_id] if not hasattr(r, '_final_compartments')]
                    rtk_id = [r.id for r in rtk]
                    reactions_to_remove += [r.id for r in self.id_reactions_map[hgnc_id] if r.id not in rtk_id]
                    self.id_reactions_map[hgnc_id] = rtk
                else:
                    reactions_to_remove += [r.id for r in self.id_reactions_map[hgnc_id]]
                    del self.id_reactions_map[hgnc_id]

            self.me_reactions = [r for r in self.me_reactions if r.id not in reactions_to_remove]
            n_reactions = len(self.me_reactions) + len(self.complex_formation_reactions)
            if self.deg_args['complex_degradation']:
                n_reactions += len(self.complex_degradation_reactions)

            print(
                'A total of {} reactions were dropped when forming a minimal proteome'.format(n_reactions_og - n_reactions))

    def add_metabolic_machinery(self):
        """Couple metabolic machinery."""
        # deal with metabolic reactions first
        print('Add machinery to metabolic module reactions')
        self._check_catalysis_coefficient = {}
        metabolic_reactions = [r.id for r in self.m_model.reactions]
        reaction_counter = dict(zip(sorted(set(metabolic_reactions)), [0] * len(metabolic_reactions)))
        final_reactions = []

        if self.check_all and self.complex_df[(self.complex_df.category == 'metabolic_reaction') &
                                              (self.complex_df.knock_out)].reaction_id.value_counts().unique() != np.array([1]):
            raise ValueError('Internal: Something went wrong in formatting complex df for knock out')

        for i in tqdm(self.complex_df[self.complex_df.category == 'metabolic_reaction'].index):
            reaction_id = self.complex_df.loc[i, 'reaction_id']  # original reaction id
            r = to_metabolic_reaction(model_metabolites=self.model_metabolites, 
                                    reaction=self.m_model.reactions.get_by_id(reaction_id))

            ko = self.complex_df.loc[i, 'knock_out']
            if not ko:
                if not self.complex_df.loc[i, 'is_complex']:
                    enzyme_to_couple = self.id_protein_map[self.complex_df.loc[i, 'machinery']][
                        self.complex_df.loc[i, 'compartment']]

                    # back track assign synthesis attribute to monomeric enzymes

                    # FUTURE: if not (enzyme_to_couple.compartment == 'r' and 'og' in enzyme_to_couple._ptms)
                    # if including PTMs in machinery in future, ER enzymes with OG PTM will have synthesis
                    # reaction be "COPI_RETROtr", which is a protein degradation reaction

                    sr_tracker = set()
                    for sr in list(enzyme_to_couple.reactions):
                        if enzyme_to_couple in sr.products and not isinstance(sr, ProteinDegradationReaction):
                            sr.synthesis, sr.synthesis_type = True, 'protein'
                            sr_tracker.add(sr.id)
                    if len(sr_tracker) != 1:
                        raise ValueError(enzyme_to_couple.id + ' has an incorrect number of associated synthesis reactions')
                else:
                    enzyme_to_couple = self.complex_id_metabolite_map[self.complex_df.loc[i, 'complex_id']]
                    enzyme_to_couple.get_k_deg()
                    if self.check_all and len([1 for r in list(enzyme_to_couple.reactions) if (
                            hasattr(r, 'synthesis') and r.synthesis and enzyme_to_couple in r.products)]) != 1:
                        raise ValueError(
                            enzyme_to_couple.id + ' has an incorrect number of associated synthesis reactions')
                enzyme_to_couple.keff = self.complex_df.loc[i, 'keff']

                # add machinery to substrate side
                c3 = (params.mu + enzyme_to_couple.k_deg) / enzyme_to_couple.keff

                if self.check_all:
                    if not enzyme_to_couple.enzyme:
                        if enzyme_to_couple.id in self._check_catalysis_coefficient:
                            raise ValueError('Enzyme exists but is not classified as one')
                        self._check_catalysis_coefficient[enzyme_to_couple.id] = [c3]
                    else:
                        self._check_catalysis_coefficient[enzyme_to_couple.id] += [c3]

                    if c3.subs(params.mu, 1) <= 0:
                        raise ValueError('The catalysis coupling constraint is negative for ' + enzyme_to_couple.id)
                enzyme_to_couple.couple(type='catalysis', value=-c3)

                if not r.reversibility:
                    r.couple(metabolites=enzyme_to_couple, types='catalysis')
                    reactions = [r]
                else:  # add a forward and reverse reaction for reversible reactions
                    r_f, r_r = r.copy(), r.copy()
                    r_f.lower_bound, r_r.lower_bound, r_r.upper_bound = 0, 0, abs(r.lower_bound)
                    r_r.add_metabolites({metab: -coeff for metab, coeff in r_r.metabolites.items()}, combine=False)

                    r_f.couple(metabolites=enzyme_to_couple, types='catalysis')
                    r_r.couple(metabolites=enzyme_to_couple, types='catalysis')

                    r_f.id, r_r.id = r_f.id + '_F', r_r.id + '_R'
                    reactions = [r_f, r_r]
            else:  # block metabolic reaction for knocked out genes
                r.lower_bound, r.upper_bound = 0, 0
                reactions = [r]

            # if multiple of the same reaction with different machinery due to OR GPR, add a different id for each
            if self.complex_df.loc[i, 'creates_multiple_reactions']:
                if len(reactions) > 1:
                    for j, r_ in enumerate(reactions):
                        r_.id = r_.id + '_' + str(reaction_counter[reaction_id])
                        reactions[j] = r_
                if reaction_counter[reaction_id] == 0:  # tracking that all metabolic reactions are added
                    metabolic_reactions.remove(reaction_id)
                reaction_counter[reaction_id] += 1
            else:
                metabolic_reactions.remove(reaction_id)  # tracking that all metabolic reactions are added
            final_reactions += reactions

        # dummy protein for orphan reactions (see deorphan)
        if sorted(metabolic_reactions) != sorted([r.id for r in self.m_model.reactions if len(r.genes) == 0]):
            raise ValueError('Not all metabolic reactions that require machinery have been accounted for')
        if self.dummy_protein is None:
            self.orphan = [to_metabolic_reaction(model_metabolites=self.model_metabolites, reaction=r) for r in self.m_model.reactions if len(r.genes) == 0]
            final_reactions += self.orphan
            self.deorphaned = list()
        self.final_reactions = final_reactions

    def add_expression_machinery(self):
        """Couple expression machinery."""
        # filter out metabolic reactions
        backup = self.complex_df.copy()
        self.complex_df = self.complex_df[self.complex_df.category == 'expression_reaction']
        self.complex_df.reset_index(inplace=True, drop=True)

        if not self.deg_args['complex_degradation']:
            expression_reactions = self.me_reactions
        else:
            expression_reactions = self.me_reactions + [r for r in self.complex_degradation_reactions if
                                                        not r._ribosomal_degradation]
            # filter out ribosomal_degradation reactions
        expression_reactions = [r for r in expression_reactions if len(r.genes) > 0]
        reaction_counter = dict(zip(sorted(set([r.id for r in expression_reactions])), [0] * len(expression_reactions)))

        print('Add machinery to expression module reactions')
        for rxn in tqdm(expression_reactions):
            if type(rxn) != ComplexDegradationReaction:
                reaction_id_short = func.parse_me_reaction_id(rxn.id)  # abbreviated version
            else:
                reaction_id_short = parse_complex_degradation_reaction_id(rxn.id)

            reaction_id = rxn.id  # original reaction id
            idx = self.complex_df[self.complex_df.reaction_id == reaction_id_short].index.tolist()
            for i in idx:
                r = rxn.copy()
                #                 r._metabolites = rxn.metabolites
                if not self.complex_df.loc[i, 'is_complex']:
                    enzyme_to_couple = self.id_protein_map[self.complex_df.loc[i, 'machinery']][
                        self.complex_df.loc[i, 'compartment']]
                    # back track assign synthesis attribute to monomeric enzymes

                    sr_tracker = set()
                    for sr in list(enzyme_to_couple.reactions):
                        if enzyme_to_couple in sr.products and not isinstance(sr, ProteinDegradationReaction):
                            sr.synthesis, sr.synthesis_type = True, 'protein'
                            sr_tracker.add(sr.id)
                    if len(sr_tracker) != 1:
                        raise ValueError(enzyme_to_couple.id + ' has an incorrect number of associated synthesis reactions')
                else:
                    enzyme_to_couple = self.complex_id_metabolite_map[self.complex_df.loc[i, 'complex_id']]
                    enzyme_to_couple.get_k_deg()
                    if self.check_all and len([1 for r in list(enzyme_to_couple.reactions) if (
                            hasattr(r, 'synthesis') and r.synthesis and enzyme_to_couple in r.products)]) != 1:
                        raise ValueError(
                            enzyme_to_couple.id + ' has an incorrect number of associated synthesis reactions')
                enzyme_to_couple.keff = self.complex_df.loc[i, 'keff']

                # add machinery to substrate side
                c3 = (params.mu + enzyme_to_couple.k_deg) / enzyme_to_couple.keff

                if self.check_all:
                    if not enzyme_to_couple.enzyme:
                        if enzyme_to_couple.id in self._check_catalysis_coefficient:
                            raise ValueError('Enzyme exists but is not classified as one')
                        self._check_catalysis_coefficient[enzyme_to_couple.id] = [c3]
                    else:
                        self._check_catalysis_coefficient[enzyme_to_couple.id] += [c3]

                    if c3.subs(params.mu, 1) <= 0:
                        raise ValueError('The catalysis coupling constraint is negative for ' + enzyme_to_couple.id)
                enzyme_to_couple.couple(type='catalysis', value=-c3)

                if not r.reversibility:
                    r.couple(metabolites=enzyme_to_couple, types='catalysis')
                    reactions = [r]
                else:  # add a forward and reverse reaction for reversible reactions
                    r_f, r_r = r.copy(), r.copy()
                    r_f.lower_bound, r_r.lower_bound, r_r.upper_bound = 0, 0, abs(r.lower_bound)
                    r_r.add_metabolites({metab: -coeff for metab, coeff in r_r.metabolites.items()}, combine=False)

                    r_f.couple(metabolites=enzyme_to_couple, types='catalysis')
                    r_r.couple(metabolites=enzyme_to_couple, types='catalysis')

                    r_f.id, r_r.id = r_f.id + '_F', r_r.id + '_R'
                    reactions = [r_f, r_r]

                # if multiple of the same reaction with different machinery due to OR GPR, add a different id for each
                if self.complex_df.loc[i, 'creates_multiple_reactions']:
                    if len(reactions) > 1:
                        for j,r_ in enumerate(reactions):
                            r_.id = r_.id + '_' + str(reaction_counter[reaction_id])
                            reactions[j] = r_
                    reaction_counter[reaction_id] += 1
                self.final_reactions += reactions

        if self.deg_args['complex_degradation']:
            # hard-coded for ribosomal degradation
            enzyme_to_couple = [self.complex_id_metabolite_map[self.complex_df[
                self.complex_df.reaction_id == 'PROTEASOMAL_DEGRADATIONc'].complex_id.tolist()[0]]]
            enzyme_to_couple.append(self.complex_id_metabolite_map[self.complex_df[
                self.complex_df.reaction_id == '5s_rRNA_DEGRADATIONc'].complex_id.tolist()[0]])
            ribosomal_degradation_reactions = [r for r in self.complex_degradation_reactions if
                                               r._ribosomal_degradation]

            for rxn in ribosomal_degradation_reactions:
                rxn.couple(metabolites=enzyme_to_couple, types=['catalysis', 'catalysis'])
                self.final_reactions.append(rxn)

        if self.dummy_protein is None:
            me_orphans = [r for r in self.me_reactions if len(r.genes) == 0]
            if self.deg_args['complex_degradation']:
                me_orphans += [r for r in self.complex_degradation_reactions if len(r.genes) == 0]

            me_orphans += self.complex_formation_reactions
            self.orphan += me_orphans
            self.orphan = [r.id for r in self.orphan + self.biomass_reactions]
            self.final_reactions += me_orphans
            del me_orphans

        self.complex_df = backup.copy()
        del backup

        if self.check_all:
            for k, v in self._check_catalysis_coefficient.items():
                if len(set(v)) != 1:
                    raise ValueError(k + ' received multiple coupling coefficients for different reactions')
        del self._check_catalysis_coefficient

    def deorphan(self):
        """Couples dummy protein to reactions that don't have specified genes ("de-orphaning")

        Returns
        ----------
        deorphaned: list
            a list of ME_Model reaction IDs for reactions that were de-orphaned
        self.orphan: list
            a list of ME_Model reaction IDs for reactions there were not de-orphaned despite having 0 specified genes 

        """

        if self.dummy_protein is not None:
            print('Deorphan enzymeless reactions')
            enzymeless_reactions = [to_metabolic_reaction(model_metabolites=self.model_metabolites, reaction=r) for r in self.m_model.reactions if len(r.genes) == 0]
            enzymeless_reactions_map = {r.cobra_id: r for r in enzymeless_reactions}
            enzymeless_reactions += [r for r in self.me_reactions if
                                     len(r.genes) == 0] + self.complex_formation_reactions
            if self.deg_args['complex_degradation']:
                enzymeless_reactions += [r for r in self.complex_degradation_reactions if len(r.genes) == 0]
            boundary_ids = [r.id for r in self.m_model.exchanges + self.m_model.demands]

            if len(set([r.id for r in enzymeless_reactions]).intersection([r.id for r in self.final_reactions])) > 0:
                raise ValueError('Incorrect parsing of reaction lists for dummy protein')
            # if exclude is None
            # metabolic module enzymes to exclude from deorphaning - boundary reactions
            self.orphan = [r for r in enzymeless_reactions_map.values() if
                           hasattr(r, 'cobra_id') and r.cobra_id in boundary_ids]
            _orphan = list()
            for r in self.orphan:  # secondary exchange reactions
                if len(r.metabolites) > 1 or list(r.metabolites)[0].compartment != 'b':
                    raise ValueError(
                        'Incorrectly formatted exchange reaction: ' + r.id + '. Must follow Recon2.2 format.')

                assoc_rxn = [r_.id for r_ in list(list(self.m_model.reactions.get_by_id(r.id).metabolites)[0].reactions)]
                assoc_rxn.remove(r.cobra_id)

                if len(assoc_rxn) > 0:
                    for r_id in assoc_rxn:  # id the second exchange reaction (Recon2.2 format)
                        r_ = self.m_model.reactions.get_by_id(r_id)
                        cond1 = (sorted(r_.compartments) == ['b', 'e'])
                        cond2 = (len(set(['_'.join(m.id.split('_')[:-1]) for m in list(r_.metabolites)])) == 1)
                        cond3 = (len(r_.genes) == 0)
                        if cond1 and cond2 and cond3:
                            _orphan.append(enzymeless_reactions_map[r_.id])
            self.orphan += _orphan
            del _orphan
            # # Deprecated
            # if exclude is not None:
            #     for r_id in exclude:
            #         if r_id not in m_ids:
            #             raise ValueError('The list of metabolic reactions to exclude from dummy catalysis must be in the metabolic model reaction list')
            #         if len(self.m_model.reactions.get_by_id(r_id).genes)>0:
            #             raise ValueError('The list of metabolic reactions to exclude from dummy catalysis must not have an associated GPR')
            #
            #     self.orphan = [to_metabolic_reaction(r) for r in exclude]

            # expression module enzymes to exclude
            expression_rids = ['CYTOSOLIC_PROTEIN_FOLDING', 'IMPORTtn',
                               'RIBOSOME_COMPLEX_DISSOCIATIONc', 'UNFOLDr',
                               'POLYUBIQUITIN_MOIETY_EXPORTtn', 'COMPLEX_FORMATION']
            for r in enzymeless_reactions:
                for expr_rid in expression_rids:
                    if expr_rid in r.id:
                        self.orphan.append(r)
                        break

            # do not deorphan transport reactions for small molecules (can passively diffuse)
            transport = [r for r in enzymeless_reactions if len(r.compartments) > 1 and r not in self.orphan]

            # # DEPRECATED (inaccurate groups ) -- include reactions listed in transport groups
            # m_transport = list(set(func.flatten_list([[group_rxn.id for group_rxn in group._members] for \
            #                                           group in self.m_model.groups if g.id.startswith('Transport')])))
            # transport = [r for r in enzymeless_reactions if len(r.compartments) > 1 or \
            #              (hasattr(r, 'cobra_id') and r.cobra_id in m_transport)] # not all of these are transport
            # transport = [r for r in transport if r not in self.orphan]

            # note using [g._members for g in model.groups if 'Transport' in g.name] results in some reations that are incorrectly classified as trasnport
            # remove_idx = list()
            # exclude reactions stated as "diffusion" by the reaction name
            diffusion_reactions = [r.id for r in self.m_model.reactions if 'via diffusion' in r.name]
            for r in transport:
                # filter for molecules that are actually transported (in reactants and products) and check
                # if they are over the diffusion limit or charged, or annotated as diffusion in reaction name
                diffusion_reaction = False
                if hasattr(r, 'cobra_id') and r.cobra_id in diffusion_reactions:
                    diffusion_reaction = True

                actual_transport_m = func.determine_transport(r)
                sm_transport = [m for m in r.reactants if (m.id.split('_')[:-1][0] in actual_transport_m) and (
                    (m.formula_weight is not None and m.formula_weight > params.MEMBRANE_DIFFUSION_LIMIT) or (
                        m.charge != 0))]
                if len(sm_transport) == 0:  # if all transported metabolites are under diffusion limit and uncharged,:
                    diffusion_reaction = True
                if diffusion_reaction:
                    self.orphan.append(r)

            #                 # V2
            #                 sm_prod = [m.id.split('_')[:-1][0] for m in r.products]
            #                 sm_transport = [m for m in r.reactants if (m.id.split('_')[:-1][0] in sm_prod) and \
            #                                    ((m.formula_weight is not None and m.formula_weight > membrane_diffusion_limit) \
            #                                     or (m.charge != 0))]
            #                 if len(sm_transport) == 0:
            #                     self.orphan.append(r)
            #                # V1
            #                 r = transport[i]
            #                 tm = dict()
            #                 mc = dict()
            #                 counter = 0
            #                 for m in r.metabolites:
            #                     if m.formula_weight <= params.MEMBRANE_DIFFUSION_LIMIT: # all metabolites within diffusion limit
            #                         counter += 1
            #                     m_id = '_'.join(m.id.split('_')[:-1]) # atleast one metabolite is transported across compartments
            #                     if m_id not in tm:
            #                         tm[m_id] = 1
            #                     else:
            #                         tm[m_id] += 1
            #                         mc[m_id] = m.charge

            #                 uncharged = True
            #                 for m_id,v in tm.items():
            #                     if v >= 2 and mc[m_id] != 0:
            #                         uncharged = False
            #                         break
            #                 # uncharged, all metabolites that are transported are < 504 Da, and atleast one metabolite is transported
            #                 if max(list(tm.values())) >= 2 and counter == len(r.metabolites) and uncharged:
            #                     self.orphan.append(r)

            deorphan = [r for r in enzymeless_reactions if r not in self.orphan]
            self.deorphaned = list()

            if len(deorphan) > 0:
                c3 = (params.mu + self.dummy_protein['protein_metabolite'].k_deg) / self.dummy_protein[
                    'protein_metabolite'].keff
                self.dummy_protein['protein_metabolite'].couple(type='catalysis', value=-c3)

                for r in deorphan:
                    if not r.reversibility:
                        r.couple(metabolites=self.dummy_protein['protein_metabolite'], types='catalysis')
                        reactions = [r]
                    else:  # add a forward and reverse reaction for reversible reactions
                        r_f, r_r = r.copy(), r.copy()
                        r_f.lower_bound, r_r.lower_bound, r_r.upper_bound = 0, 0, abs(r.lower_bound)
                        r_r.add_metabolites({metab: -coeff for metab, coeff in r_r.metabolites.items()}, combine=False)

                        r_f.couple(metabolites=self.dummy_protein['protein_metabolite'], types='catalysis')
                        r_r.couple(metabolites=self.dummy_protein['protein_metabolite'], types='catalysis')
                        r_f.id, r_r.id = r_f.id + '_F', r_r.id + '_R'
                        reactions = [r_f, r_r]
                    self.deorphaned += reactions
            self.final_reactions += self.orphan + self.deorphaned
            self.orphan = [r.id for r in self.orphan + self.biomass_reactions + [biomass.upb_reaction]]
            for r in self.deorphaned:
                r.enzyme_compartment = 'c'
            self.deorphaned = [r.id for r in self.deorphaned]

    def incorporate_protein_degradation(self):
        """Removes degradation reactions of inactive monomers and couples protein degradation to catalysis, depending on 
        deg_args input."""

        if self.check_all and self.deg_args['complex_degradation']:
            for r in self.complex_degradation_reactions:
                r._update_enzymes()  # updates rxn ._enzymes attribute to include all macromolecules involved in reaction catalysis
                if not len(r._enzymes) > 0:
                    raise ValueError(r.id + ': this ComplexDegradationReaction is not associated with an active enzyme')

        if not self.deg_args['nonenzyme_degradation']:  # remove degradation reactions of nonenzymes (degraded in complex)
            reactions_to_remove = list()

            nonmachinery_exceptions = list()
            for hgnc_id, compartments in self.non_machinery.items():
                compartments = set(compartments)
                nonmachinery_exceptions += [r.id for r in self.id_reactions_map[hgnc_id] if
                                            isinstance(r, ProteinDegradationReaction) and len(
                                                compartments.intersection(r._final_compartments)) > 0]

            pdr = [r for r in self.me_reactions if
                   isinstance(r, ProteinDegradationReaction) and r.id not in nonmachinery_exceptions]
            for r in pdr:
                r._update_enzymes()
                if not len(r._enzymes) > 0:
                    reactions_to_remove.append(r.id)
            print(
                '{} of {} protein degradation reactions will be removed because they are not associated with an active enzyme'.format(
                    len(reactions_to_remove), len(pdr)))

            if self.check_all and len(set(reactions_to_remove).difference([r.id for r in self.final_reactions])) > 0:
                raise ValueError('Untracked protein degradation reactions (not in final reactions list)')

            self.final_reactions = [r for r in self.final_reactions if r.id not in reactions_to_remove]
            self.orphan = [r_id for r_id in self.orphan if r_id not in reactions_to_remove]
            self.deorphaned = [r_id for r_id in self.deorphaned if r_id not in reactions_to_remove]

        if self.deg_args['couple']:
            print('Couple enzyme degradation to catalysis')

            pdr = [r for r in self.me_reactions if isinstance(r, ProteinDegradationReaction)]
            dr_map = dict()
            for r in pdr + self.complex_degradation_reactions:
                r._update_enzymes()
                dr_map[r.id] = r
            catalysis_reactions = [r for r in self.final_reactions if
                                   r.coupled_metabolites != dict() and 'catalysis' in r.coupled_metabolites.values() and r.enzyme_compartment != 'e' and not (
                                       r.lower_bound == r.upper_bound == 0)]

            # TODO: couple ribosomal degradation
            # TEMPORARY: don't couple ribosomal degradation for now - change in ExpressedGene._check_macromolecules
            catalysis_reactions = [r for r in catalysis_reactions if 'mrna_formation' not in r.coupled_metabolites.values()]

            for r in tqdm(catalysis_reactions):
                enzymes = [m for m, t in r.coupled_metabolites.items() if t == 'catalysis']

                deg_reactions = list()
                deg_proxies = list()
                for e in enzymes:  # in case multiple catalysis proteins (ribosomal degradatio)
                    deg_reactions_ = [r_id for r_id in e._degradation_reactions if (dr_map[r_id].sink) and (e in dr_map[r_id]._enzymes)]
                    if len(deg_reactions_) == 0:
                        raise ValueError('No degradation reactions associated with catalyzing enzyme')
                    if len(deg_reactions_) > 1:
                        raise ValueError('More than 1 degradation reaction associated with catalyzing enzyme')
                    deg_reactions += deg_reactions_
                    dp = e.make_proxy()
                    dp.couple(value=-e.k_deg / e.keff)
                    deg_proxies.append(dp)

                if len(deg_reactions) > 1 and not (r.id in dr_map or dr_map[r.id]._ribosomal_degradation):
                    raise ValueError('More than 1 degradation reaction associated with the catalysis reaction')

                deg_reactions = [r_ for r_ in self.final_reactions if r_.id in deg_reactions]
                for dr, dp in list(zip(deg_reactions, deg_proxies)):
                    # keeping track of whether the enzyme degradation reaction already has a proxy metabolite for
                    # coupling (occurs in scenarios where an enzyme catalyzes multiple reactions)
                    if not dr._protein_deg_proxy:
                        dr._add_protein_deg_proxy(dp)
                    else:
                        dp = dr.protein_deg_proxy

                    if (self.check_all) and ('enzyme_degradation' in r.coupled_metabolites.values()) and (not dr_map[r.id]._ribosomal_degradation):
                        raise ValueError('This reaction already is coupled to degradation')
                    r.couple(metabolites=dp, types='enzyme_degradation')
                    # .couple works in scenarios where r == dr because .couple uses .add_metabolites(combine = True)

    def build_me_model(self):
        """Generate the final ME_Model object."""
        print('Add biomass component to reactions')

        for r in self.final_reactions:
            biomass.add_biomass_change(r)

        #         br.append(self.pb_reaction)
        if self.dummy_protein is not None:
            self.biomass_reactions.append(biomass.upb_reaction)

        if len([r for r in self.final_reactions if not isinstance(r, core.reaction.ME_Reaction)]) > 0:
            raise ValueError('Internal: Reactions not of type ME_Reaction are included in the model')
        self.final_reactions += self.biomass_reactions

        print('Generate ME-Model')
        me_model = ME_Model(m_model=self.m_model, id_or_model=self.model_id, n_cores=self.n_cores,
                            knock_out=self.knock_out, non_machinery=self.non_machinery,
                            additional_ko=self.additional_ko)

        # note, at end of .add_reactions() method, we reassign .coupled_metabolites attribute
        # running .add_metabolic_reactions() code outside of ME_Builder object doesn't create disagreement
        # between r.metabolites and r.coupled_metabolites, but running the method on the object does
        me_model.add_reactions(self.final_reactions)
        # TODO: incorporate the following two lines into ME_Model class instead
        me_model.reaction_types['orphan'] = self.orphan
        me_model.reaction_types['deorphaned'] = self.deorphaned

        me_model.check()
        me_model._generate_expressed_genes()

        #         del self.pb_reaction
        #         del self.ub_args
        del self.me_reactions
        del self.final_reactions
        del self.complex_formation_reactions
        del self.complex_degradation_reactions
        del self.m_model
        del self.orphan
        del self.deorphaned

        return me_model


def build_me(me_input_model: Union[cobra.Model,str],
             psim_me: Union[pd.DataFrame, str],
             model_id: str = 'HUMAN_ME_MODEL',
             stochastic: bool = False, seed: int = 888, n_cores: int = os.cpu_count(),
             non_machinery: Optional[Dict[str, List[str]]] = None, knock_out: Optional[List[str]] = None,
             dummy_protein: bool = True, context_specific_dummy: bool = False,
             minimal_proteome: bool = True, compress_mrna: bool = True,
             check_all: bool = True,
             deg_args: Dict[str, bool] = {'couple': True, 'reversible_complex_formation': False, 'nonenzyme_degradation': False,
                                          'complex_degradation': True}
             ):
    """Build a human ME model according to input M-Model (as provided in preprocess.correct_inputs.correct_model), PSIM (as provided in preprocess.correct_inputs.correct_psim), 
    and the parameters below. 

    Parameters
    ----------
    me_input_model : cobra.Model
        the corrected input metabolic model (as provided in preprocess.correct_inputs.correct_model) or 'full/path/to/corrected_model.xml'
    psim_me : pd.DataFrame
        the corrected psim matrix (as provided in preprocess.correct_inputs.correct_psim) or 'full/path/to/corrected_psim.csv or .h5'
    model_id : str, optional
        model identifier, by default 'HUMAN_ME_MODEL'
    stochastic : bool, optional
        Whether a potentially stochastic output should be stochastic, or choose a default behavior instead, by default False
    seed : int, optional
        A seed for if stochastic is set to True, by default 888
    n_cores : int, optional
        # of cores to parallelize on, by default os.cpu_count()
    non_machinery : Dict[str, List[str]], optional
        keys are HGNC IDs, values are a list wherein element represents a compartment within the model for the gene to be expressed, by default None
        We define machinery as proteins that are utilized in the reaction GPRs. In its current format, it is not possible for a protein to both be machinery and non-machinery
    knock_out : List[str], optional
        each element is the HGNC ID of a gene expressed in the model which should be knocked out, by default None
        *Note: you may want to knock-out during building if setting minimal_proteome = True and knocking out a 
        gene that participates in a OR GPR rule (in case it is the one that is selected by minimal proteome); otherwise ME_Model.knock_out() method should suffice
    dummy_protein : bool, optional
        whether to add a representative dummy protein to catalyze orphan reactions, by default True
    context_specific_dummy : bool, optional
        whether the representative dummy protein is calculated for only genes in the user-provided context specific model from
        the user provided PSIM (True) or for all recon2.2 machinery proteins in the gold-standard PSIM (False), by default False
    minimal_proteome : bool, optional
        For reactions with OR in the GPR, whether the builder generates a separate reaction for each protein complex (False) 
        or just one reaction, choosing the protein complex with the lowest molecular weight to catalyze the reaction (True). 
        If a reaction has multiple enzyme options with the same molecular weight, will randomly choose one. 
        Will not consider a complex that contains a knocked out gene. 
        by default True
    compress_mrna : bool, optional
         whether to condense elongation, processing, and nuclear export reactions into a single reaction, by default True
    check_all : bool, optional
         Whether to check that building is proceeding correctly. Increases run time, by default True
    deg_args : Dict[str, bool], optional
        A number of options related to protein and complex degradation. Becomes important in slow growth conditions.
        Note the default values focus on coupling fluxes and degrading the specific enzymes associated with 
        each reaction. 
        By default {'couple': True, 'reversible_complex_formation': False, 'nonenzyme_degradation': False, 'complex_degradation': True}

        Key value pairs:
            "couple": bool
                Whether to explicitly couple enzyme degradation reactions to metabolic catalysis. Becomes 
                particularly important in slow growth conditions.
            "reversible_complex_formation": bool
                Whether reactions to form complexes are reversible (<->) or not (-->). Setting to True may make
                model more efficient (reuse of proteins involved in catalysis of multiple reactions in same compartment)
            "nonenzyme_degradation": bool
                Whether to retain degradation reactions (associated with the build_protein_expression script) for
                proteins that form complexes rather than become monomeric enzymes; i.e., all individual complex 
                subunits have their own protein degradation reaction. Note that even if set to False, 
                protein intermediates associated with the monomeric enzyme that had degradation rections are retained.
                Regardless of this parameter, only the specific enzymatic degradation reaction associated with the 
                catalysis reaction will be coupled. Independent of complex_degradation and 
                reversible_complex_formation arguments.
            "complex_degration": bool
                Whether to generate degradation reactions for whole complexes in addition to individual monomers
                (required for coupling)

    Returns
    -------
    me_model : ME_Model
        the fully constructed ME-Model
    builder : MEBuilder
        instance of MEBuilder used to generate the me_model
    """
    start = time.time()
    builder = MEBuilder(m_model=me_input_model, psim_me=psim_me, model_id=model_id, 
                        stochastic=stochastic, seed=seed, n_cores=n_cores,
                        non_machinery=non_machinery, knock_out=knock_out,
                        dummy_protein=dummy_protein, context_specific_dummy=context_specific_dummy,
                        minimal_proteome=minimal_proteome, compress_mrna=compress_mrna,
                        check_all=check_all, deg_args=deg_args)
    builder.express_metabolic_enzymes()
    builder.express_expression_enzymes()
    builder.express_dummy_protein()
    builder.get_complex_info()
    builder.generate_complex_reactions()
    builder.get_keff()
    builder.minimize_proteome()
    builder.add_metabolic_machinery()
    builder.add_expression_machinery()
    builder.deorphan()
    builder.incorporate_protein_degradation()
    me_model = builder.build_me_model()

    end = time.time()
    print('Time to build: {:.2f} minutes'.format((end - start) / 60))

    return me_model, builder


gc.collect()
