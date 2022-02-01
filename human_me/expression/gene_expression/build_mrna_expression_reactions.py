#!/usr/bin/env python
# coding: utf-8

from human_me.utils import machinery as mach
from human_me.utils import functions as func
from human_me.utils.polyA_statistics import calculate_polyA_length

from human_me.core.reaction import ExpressionReaction
from human_me.core.macromolecules.RNA import RNA_fragment, pre_mRNA, mRNA


class ExpressMrna:
    '''Gene-specific mRNA expression.'''

    def __init__(self, gene_info, model_metabolites) -> None:
        """Init method for ExpressMrna

        Parameters
        ----------
        gene_info : gene_information.GeneInformation
            GeneInformation object of gene to be expressed
        model_metabolites : utils.metabolites.MetaboliteBin
            the me_input_model metabolites as specified by MetaboliteBin
        """
        self.reactions = []
        self.lariat = None
        self.gene_info = gene_info
        self.model_metabolites = model_metabolites

    def transcribe_premrna(self) -> None:
        """Elongation reaction."""
        # elongation reaction
        # https://www.google.com/search?q=rna+polymerization+reaction&source=lnms&tbm=isch&sa=X&ved=2ahUKEwiN_73Vk7rqAhXOsJ4KHW5lB4UQ_AUoAXoECA4QAw&biw=1920&bih=1001#imgrc=w7XH4mHmJglCuM

        self.premrna = pre_mRNA(model_metabolites=self.model_metabolites, gene_info=self.gene_info)
        self.transcript_elongation = self.premrna.synthesize(id_=self.gene_info.hgnc_id + '_TRANSCRIPTION_ELONGATION')
        self.reactions.append(self.transcript_elongation)

    def process_mrna(self) -> None:
        """Processing includes capping, splicing, and polyA tail."""
        # combine in to one to not create too many reactions (capping itself is 4 reactions)
        # make mrna_n metabolite
        self.mrna_n = mRNA(model_metabolites=self.model_metabolites, gene_info=self.gene_info, compartment='n')

        self.polyA_length = calculate_polyA_length(self.gene_info.polyA_length, self.gene_info.stochastic,
                                                   self.gene_info.seed)
        self.mrna_n.update_metabolite(seq=''.join(['A'] * self.polyA_length),
                                      append=True, append_to='3_primed')

        # +2 for cap
        self.mrna_n.charge += 2  # (-self.polyA_length + 2)

        transcript_processing = ExpressionReaction(self.gene_info.hgnc_id + '_TRANSCRIPTION_PROCESSING',
                                                   subsystem='mRNA_expression', hgnc_id=self.gene_info.hgnc_id)
        rxn = dict()

        rxn[self.model_metabolites.atp_n], rxn[self.model_metabolites.ppi_n] = -self.polyA_length, self.polyA_length  # polyA tail

        # 5' cap: https://sites.google.com/site/learnorganicchem/organic-molecules/biomolecules/rna/rna-processing?tmpl=%2Fsystem%2Fapp%2Ftemplates%2Fprint%2F&showPrintDialog=1
        rxn[self.model_metabolites.h2o_n], rxn[self.model_metabolites.pi_n] = -1, 1  # rtpase
        rxn[self.model_metabolites.gtp_n] = -1  # gp transfer
        rxn[self.model_metabolites.ppi_n] += 1  # gp transfer
        rxn[self.model_metabolites.amet_n], rxn[self.model_metabolites.ahcys_n] = -2, 2  # methyltransferase - cap0 and cap1 structure
        rxn[self.model_metabolites.h_n] = 1  # methyltransferase cap1

        # 4 ATP consumed per capping reaction
        rxn = func.hydrolyze_atp(rxn, n_atp=4, compartment='n', model_metabolites=self.model_metabolites)

        processed_elements = self.mrna_n.elements.copy()
        for element in processed_elements.keys():
            #             processed_elements[element] += (self.polyA_length* self.model_metabolites.seq_element_map['A'][element]) # polyA tail
            processed_elements[element] += self.model_metabolites.gp[element]  # 5' cap rxn2 - addition of Gp

        # 5' cap
        processed_elements['P'] -= 1  # rxn 1: lost of third triphosphate by RTPase
        processed_elements['O'] -= 4  # rxn 1: loss of third triophosphate by RTPase
        processed_elements['C'] += 2  # rxn 3-4: methyltransferase - cap0 and cap1 structure
        processed_elements['H'] += 5  # methyltransferase - cap0 and cap1 structure
        self.mrna_n.elements = processed_elements

        rxn[self.premrna] = -1
        rxn[self.mrna_n] = 1

        # splicing
        #         if self.premrna.length > self.mrna_n.length - self.polyA_length:
        if self.gene_info.n_introns > 0:
            lariat_seq = ''
            for nt in ['A', 'U', 'G', 'C']:

                if nt != 'A':
                    diff = self.premrna.sequence.count(nt) - self.mrna_n.sequence.count(nt)
                    lariat_seq += ''.join([nt] * diff)
                else:
                    diff = self.premrna.sequence.count(nt) - (self.mrna_n.sequence.count(nt) - self.polyA_length)
                    lariat_seq += ''.join([nt] * diff)

            self.lariat = RNA_fragment(model_metabolites=self.model_metabolites, metabolite_name=self.gene_info.hgnc_id, fragment_type='lariat',
                                       seq=lariat_seq, triphosphate=False, hgnc_id=self.gene_info.hgnc_id)

            rxn[self.lariat] = 1
            rxn[self.model_metabolites.h2o_n] -= 1  # endonucleolytic cleavage
            # 10 ATP consumed per intron during splicing
            rxn = func.hydrolyze_atp(rxn, n_atp=10 * self.gene_info.n_introns, compartment='n', model_metabolites=self.model_metabolites)
            # lariat degradation - no linearization reaction (just one triphosphate consumption)
            lariat_degradation = self.lariat.exonucleolytic_degradation(
                reaction_name=self.gene_info.hgnc_id + '_lariats')
            lariat_degradation.gene_reaction_rule = mach.lm_rule
            if list(lariat_degradation.compartments) != ['n']:
                raise ValueError('Lariat degradation must be confined to nuclear compartment')
            self.reactions.append(lariat_degradation)
        else:
            lariat_degradation = None

        transcript_processing.add_metabolites(rxn)
        transcript_processing.gene_reaction_rule = ' and '.join(mach.polyA + mach.capping + mach.spliceosome)

        if len(transcript_processing.check_mass_balance()) > 0:
            raise ValueError('Transcript processing for ' + self.gene_info.hgnc_id + ' is unbalanced')
        if list(transcript_processing.compartments) != ['n']:
            raise ValueError('Transcript processing must be confined to nuclear compartment')

        self.transcript_processing = transcript_processing
        self.reactions.append(transcript_processing)

    def export_mrna(self) -> None:
        """Nuclear export of mRNA."""
        # make the cytosolic mrna metabolite
        self.mrna_c = self.mrna_n.change_compartment('c')

        # make the transport reaction
        mrna_export = ExpressionReaction(self.gene_info.hgnc_id + '_mRNA_EXPORTtn', subsystem='mRNA_expression',
                                         hgnc_id=self.gene_info.hgnc_id, synthesis=True, synthesis_type='mRNA')
        mrna_export.name = 'mRNA nuclear export'
        rxn = dict()
        rxn[self.mrna_n], rxn[self.mrna_c] = -1, 1

        # 10 ATP consumer per transcript exported
        rxn = func.hydrolyze_atp(rxn, n_atp=10, compartment='n', model_metabolites=self.model_metabolites)

        mrna_export.add_metabolites(rxn)
        # can change this GPR as an if statement in future based on following source:
        # https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.0010144.g005
        mrna_export.gene_reaction_rule = ' and '.join(mach.trex)

        if len(mrna_export.check_mass_balance()) > 0:
            raise ValueError('mRNA export for ' + self.gene_info.hgnc_id + ' is unbalanced')

        self.mrna_export = mrna_export
        self.reactions.append(mrna_export)

    def degrade_mrna(self, decapping: bool = True, three_to_five: bool = False) -> None:
        """Currently, only one of the two degradation pathways is included. We assume the 5' to 3' pathway is present.
        This is simply to limit the total number of reactions


        Parameters
        ----------
        decapping : bool, optional
            whether mRNA needs to have a decapping reaction prior to degradation, by default True
        three_to_five : bool, optional
            whether the 5' to 3' (False) or 3' to 5' (True) mRNA degradation pathway is used, by default False
        """

        rxn = self.mrna_c.exonucleolytic_degradation(reaction_name='', balanced=False)
        rxn = rxn.metabolites.copy()
        del rxn[[m for m in rxn.keys() if m.id == self.model_metabolites.ntp_map_c[self.gene_info.mrna_seq[0]].id][0]]

        # no m7g metabolite in recon2.2, so just reverse the methylation instead
        rxn[self.model_metabolites.amet_c], rxn[self.model_metabolites.ahcys_c] = 2, -2  # reverse methyltransferase - cap0 and cap1 structure

        #         proxy metabolite for coupling mRNA degradation to protein synthesis flux
        self.mrna_deg_proxy = self.mrna_c.make_proxy()
        rxn[self.mrna_deg_proxy] = 1

        h2o_c = [m for m in rxn.keys() if m.id == 'h2o_c'][0]  # won't load directly from metab for some reason
        h_c = [m for m in rxn.keys() if m.id == 'h_c'][0]

        if three_to_five:
            transcript_degradation_1 = ExpressionReaction(self.gene_info.hgnc_id + "_3'to5'_mRNA_DEGRADATIONc",
                                                          subsystem='mRNA_expression', hgnc_id=self.gene_info.hgnc_id,
                                                          sink=True, sink_type='mRNA')
            rxn_1 = rxn.copy()

            rxn_1[self.model_metabolites.ndp_map_c[self.gene_info.mrna_seq[0]]] = 1

            gmp_c = [m for m in rxn.keys() if m.id == 'gmp_c'][0]
            rxn_1[h2o_c] -= 1
            rxn_1[h_c] += 1
            rxn_1[gmp_c] += 1
            transcript_degradation_1.add_metabolites(rxn_1)
            transcript_degradation_1.gene_reaction_rule = mach.degradation_rule1

            if len(transcript_degradation_1.check_mass_balance()) > 0:
                raise ValueError('3 primed to 5 primed degradation for ' + self.gene_info.hgnc_id + ' is unbalanced')
            if list(transcript_degradation_1.compartments) != ['c']:
                raise ValueError('Transcript degradation must be confined to cytosolic compartment')

            self.reactions.append(transcript_degradation_1)
        if decapping:
            transcript_degradation_2_decapping = ExpressionReaction(
                self.gene_info.hgnc_id + "_DECAPPING_mRNA_DEGRADATIONc",
                subsystem='mRNA_expression', hgnc_id=self.gene_info.hgnc_id,
                sink=True, sink_type='mRNA')

            rxn_2 = rxn.copy()
            rxn_2[[m for m in rxn_2.keys() if m.id == self.model_metabolites.nmp_map_c[self.gene_info.mrna_seq[0]].id][0]] += 1
            # 5' cap - from 5'-->3' direction (DCP1/DCP2 - NUDIX mechanism)
            rxn_2[h2o_c] -= 1
            rxn_2[h_c] += 1
            rxn_2[self.model_metabolites.ndp_map_c['G']] = 1

            transcript_degradation_2_decapping.add_metabolites(rxn_2)
            transcript_degradation_2_decapping.gene_reaction_rule = mach.decapping_rule
            if len(transcript_degradation_2_decapping.check_mass_balance()) > 0:
                raise ValueError('Decapping degradation for ' + self.gene_info.hgnc_id + ' is unbalanced')
            if list(transcript_degradation_2_decapping.compartments) != ['c']:
                raise ValueError('Transcript degradation must be confined to cytosolic compartment')
            self.reactions.append(transcript_degradation_2_decapping)

    def compress_mrna_module(self):
        """Condense elongation, processing, and nuclear export reactions into a single reaction."""
        rxns_to_remove = [self.transcript_elongation, self.transcript_processing, self.mrna_export]
        rxn = dict()
        rxn_map = dict()
        for r in rxns_to_remove:
            for met, coeff in r.metabolites.items():
                if met.id in rxn:
                    rxn[met.id] += coeff
                else:
                    rxn[met.id] = coeff
                    rxn_map[met.id] = met
        rxn = {rxn_map[k]: v for k, v in rxn.items()}

        transcription = ExpressionReaction(self.gene_info.hgnc_id + "_TRANSCRIPTION",
                                           subsystem='mRNA_expression', hgnc_id=self.gene_info.hgnc_id,
                                           synthesis=True, synthesis_type='mRNA')
        transcription.add_metabolites(rxn)

        transcription.gene_reaction_rule = ' and '.join(
            sorted(set([item.id for sublist in [list(r.genes) for r in rxns_to_remove] for item in sublist])))

        if len(transcription.check_mass_balance()) > 0:
            raise ValueError('Condensed transcription reaction for ' + self.gene_info.hgnc_id + ' is unbalanced')

        for r in rxns_to_remove:
            self.reactions.remove(r)
        self.reactions.append(transcription)


def get_mrna_expression_reactions(gene_info, model_metabolites, compress_mrna: bool = False):
    """Generates reactions and macromolecules associated with transcription of a gene.

    Parameters
    ----------
    gene_info : GeneInformation
        representation of gene to be expressed
    model_metabolites : utils.metabolites.MetaboliteBin
        the me_input_model metabolites as specified by MetaboliteBin
    compress_mrna : bool, optional
        whether to condense elongation, processing, and nuclear export reactions into a single reaction, by default False

    Returns
    -------
    em.reactions : List[ExpressionReaction]
        the reactions to express em.mrna_c
    em.mrna_c : mRNA
        the final, cytosolic mRNA transcript
    em.mrna_deg_proxy : core.macromolecules.macromolecule.Proxy
        proxy metabolite generated in mRNA degradation reaction for coupling
    """
    em = ExpressMrna(gene_info=gene_info, model_metabolites=model_metabolites)
    em.transcribe_premrna()
    em.process_mrna()
    em.export_mrna()
    em.degrade_mrna()
    if compress_mrna:
        em.compress_mrna_module()

    return em.reactions, em.mrna_c, em.mrna_deg_proxy
