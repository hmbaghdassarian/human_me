#!/usr/bin/env python
# coding: utf-8

class MetaboliteBin:
    '''Stores all metabolites used in model building'''
    def __init__(self, me_input_model):
        """Init method for MetaboliteBin

        Parameters
        ----------
        me_input_model : cobra.Model
            the corrected input metabolic model (as provided in preprocess.correct_inputs.correct_model)
        """

        self.human_model = me_input_model

        self.atp_n = self.human_model.metabolites.get_by_id('atp_n')
        self.gtp_n = self.human_model.metabolites.get_by_id('gtp_n')
        self.ppi_n = self.human_model.metabolites.get_by_id('ppi_n')
        self.ppi_c = self.human_model.metabolites.get_by_id('ppi_c')

        # mrna expression
        # processing variables------------------------------------------------------------

        self.pi_n = self.human_model.metabolites.get_by_id('pi_n')
        self.h_n = self.human_model.metabolites.get_by_id('h_n')
        self.h2o_n = self.human_model.metabolites.get_by_id('h2o_n')
        gp = gtp_n.elements
        gp['O'] -= 6
        gp['P'] -= 2
        self.gp = gp
        self.amet_n = self.human_model.metabolites.get_by_id('amet_n')
        self.ahcys_n = self.human_model.metabolites.get_by_id('ahcys_n')
        self.adp_n = self.human_model.metabolites.get_by_id('adp_n')

        # mrna degradation------------------------------------------------------------
        self.h_c = self.human_model.metabolites.get_by_id('h_c')
        self.h2o_c = self.human_model.metabolites.get_by_id('h2o_c')
        self.pi_c = self.human_model.metabolites.get_by_id('pi_c')
        self.amet_c = self.human_model.metabolites.get_by_id('amet_c')
        self.ahcys_c = self.human_model.metabolites.get_by_id('ahcys_c')
        self.amp_c = self.human_model.metabolites.get_by_id('amp_c')

        # rrna expression
        self.atp_c = self.human_model.metabolites.get_by_id('atp_c')

        # protein expression

        # nucleus------------------------------------------------------------
        self.gdp_n = self.human_model.metabolites.get_by_id('gdp_n')

        # secretory pathway------------------------------------------------------------
        self.o2_r = self.human_model.metabolites.get_by_id('o2_r')
        self.h2o2_r = self.human_model.metabolites.get_by_id('h2o2_r')

        self.hdca_r = self.human_model.metabolites.get_by_id('hdca_r')
        self.gpi_hs_r = self.human_model.metabolites.get_by_id('gpi_hs_r')
        self.balanced_gpi = {'C': 71, 'H': 138, 'N': 4, 'O': 41, 'P': 4}  # == dgpi_prot_hs_r.elements

        self.udpacgal_g = self.human_model.metabolites.get_by_id('udpacgal_g')
        self.udpgal_g = self.human_model.metabolites.get_by_id('udpgal_g')
        self.uacgam_g = self.human_model.metabolites.get_by_id('uacgam_g')
        self.h_g = self.human_model.metabolites.get_by_id('h_g')
        self.udp_g = self.human_model.metabolites.get_by_id('udp_g')

        self.udpgal_r = self.human_model.metabolites.get_by_id('udpgal_r')
        self.uacgam_r = self.human_model.metabolites.get_by_id('uacgam_r')
        self.udpacgal_r = self.human_model.metabolites.get_by_id('udpacgal_r')
        self.udp_r = self.human_model.metabolites.get_by_id('udp_r')

        self.hdca_l = self.human_model.metabolites.get_by_id('hdca_l')
        self.h_l = self.human_model.metabolites.get_by_id('h_l')
        self.h2o_l = self.human_model.metabolites.get_by_id('h2o_l')
        self.o2_l = self.human_model.metabolites.get_by_id('o2_l')
        self.h2o2_l = self.human_model.metabolites.get_by_id('h2o2_l')

        self.udpacgal_l = self.human_model.metabolites.get_by_id('udpacgal_l')
        self.udp_l = self.human_model.metabolites.get_by_id('udp_l')

        # mrna expression
        self.seq_metabolite_map = {self.human_model.metabolites.get_by_id('utp_n'): 'U',
                            gtp_n: 'G',
                            self.human_model.metabolites.get_by_id('ctp_n'): 'C',
                            atp_n: 'A'}

        # RNA backbone elements
        seq_element_map = dict()
        for k, v in seq_metabolite_map.items():
            elements = k.elements
            elements['O'] = elements['O'] - 7  # lost from incoming ntp
            elements['P'] = elements['P'] - 2  # lost from incoming ntp
            elements['H'] = elements['H'] - 1  # lost from 3' end of growing strand
            seq_element_map[v] = elements
        self.seq_element_map = seq_element_map
        # lariat degradataion------------------------------------------------------------
        self.nmp_map_n = {'C': self.human_model.metabolites.get_by_id('cmp_n'),
                    'U': self.human_model.metabolites.get_by_id('ump_n'),
                    'G': self.human_model.metabolites.get_by_id('gmp_n'),
                    'A': self.human_model.metabolites.get_by_id('amp_n')}
        self.ntp_map_n = {v: k for k, v in seq_metabolite_map.items()}

        # mrna degradation------------------------------------------------------------
        self.nmp_map_c = {'C': self.human_model.metabolites.get_by_id('cmp_c'),
                    'U': self.human_model.metabolites.get_by_id('ump_c'),
                    'G': self.human_model.metabolites.get_by_id('gmp_c'),
                    'A': amp_c}
        self.ndp_map_c = {'C': self.human_model.metabolites.get_by_id('cdp_c'),
                    'U': self.human_model.metabolites.get_by_id('udp_c'),
                    'G': self.human_model.metabolites.get_by_id('gdp_c'),
                    'A': self.human_model.metabolites.get_by_id('adp_c')}

        # rrna expression
        self.ntp_map_c = {'C': self.human_model.metabolites.get_by_id('ctp_c'),
                    'U': self.human_model.metabolites.get_by_id('utp_c'),
                    'G': self.human_model.metabolites.get_by_id('gtp_c'),
                    'A': atp_c}

        # trna expression
        self.seq_amino_acid_map_c = {
            'A': self.human_model.metabolites.get_by_id('ala_L_c'),
            'R': self.human_model.metabolites.get_by_id('arg_L_c'),
            'N': self.human_model.metabolites.get_by_id('asn_L_c'),
            'D': self.human_model.metabolites.get_by_id('asp_L_c'),
            'C': self.human_model.metabolites.get_by_id('cys_L_c'),
            'E': self.human_model.metabolites.get_by_id('glu_L_c'),
            'Q': self.human_model.metabolites.get_by_id('gln_L_c'),
            'G': self.human_model.metabolites.get_by_id('gly_c'),
            'H': self.human_model.metabolites.get_by_id('his_L_c'),
            'I': self.human_model.metabolites.get_by_id('ile_L_c'),
            'L': self.human_model.metabolites.get_by_id('leu_L_c'),
            'K': self.human_model.metabolites.get_by_id('lys_L_c'),
            'M': self.human_model.metabolites.get_by_id('met_L_c'),
            'F': self.human_model.metabolites.get_by_id('phe_L_c'),
            'P': self.human_model.metabolites.get_by_id('pro_L_c'),
            'S': self.human_model.metabolites.get_by_id('ser_L_c'),
            'T': self.human_model.metabolites.get_by_id('thr_L_c'),
            'W': self.human_model.metabolites.get_by_id('trp_L_c'),
            'Y': self.human_model.metabolites.get_by_id('tyr_L_c'),
            'V': self.human_model.metabolites.get_by_id('val_L_c'),
        }

        self.seq_amino_acid_map_m = {aa_code: self.human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_m')
                                for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_l = {aa_code: self.human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_l')
                                for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_x = {aa_code: self.human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_x')
                                for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_n = {aa_code: self.human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_n')
                                for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_r = {aa_code: self.human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_r')
                                for aa_code, aa_metabolite in seq_amino_acid_map_c.items()}

        self.seq_amino_acid_map_compartments = {'c': seq_amino_acid_map_c, 'x': seq_amino_acid_map_x, 'r': seq_amino_acid_map_r,
                                        'm': seq_amino_acid_map_m, 'n': seq_amino_acid_map_n, 'l': seq_amino_acid_map_l}

        self.adp_c = ndp_map_c['A']

        self.atp_m = self.human_model.metabolites.get_by_id('atp_m')
        self.adp_m = self.human_model.metabolites.get_by_id('adp_m')
        self.h_m = self.human_model.metabolites.get_by_id('h_m')
        self.pi_m = self.human_model.metabolites.get_by_id('pi_m')
        self.h2o_m = self.human_model.metabolites.get_by_id('h2o_m')
        self.h_i = self.human_model.metabolites.get_by_id('h_i')

        self.h_x = self.human_model.metabolites.get_by_id('h_x')
        self.h2o_x = self.human_model.metabolites.get_by_id('h2o_x')
        self.pi_x = self.human_model.metabolites.get_by_id('pi_x')
        self.atp_x = self.human_model.metabolites.get_by_id('atp_x')
        self.adp_x = self.human_model.metabolites.get_by_id('adp_x')

        self.h_r = self.human_model.metabolites.get_by_id('h_r')
        self.h2o_r = self.human_model.metabolites.get_by_id('h2o_r')
        self.pi_r = self.human_model.metabolites.get_by_id('pi_r')
        self.atp_r = self.human_model.metabolites.get_by_id('atp_r')
        self.adp_r = self.human_model.metabolites.get_by_id('adp_r')

        self.pi_l = self.human_model.metabolites.get_by_id('pi_l')
        self.atp_l = self.human_model.metabolites.get_by_id('atp_l')
        self.adp_l = self.human_model.metabolites.get_by_id('adp_l')

        self.atp_compartments = {'c': atp_c, 'm': atp_m, 'i': atp_m, 'x': atp_x, 'n': atp_n, 'r': atp_r, 'l': atp_l}
        self.adp_compartments = {'c': adp_c, 'm': adp_m, 'i': adp_m, 'x': adp_x, 'n': adp_n, 'r': adp_r, 'l': adp_l}
        self.h2o_compartments = {'c': h2o_c, 'm': h2o_m, 'i': h2o_m, 'x': h2o_x, 'n': h2o_n, 'r': h2o_r, 'l': h2o_l}
        self.pi_compartments = {'c': pi_c, 'm': pi_m, 'i': pi_m, 'x': pi_x, 'n': pi_n, 'r': pi_r, 'l': pi_l}
        self.h_compartments = {'c': h_c, 'm': h_m, 'i': h_i, 'x': h_x, 'n': h_n, 'r': h_r, 'l': h_l}

        self.datp_n = self.human_model.metabolites.get_by_id('datp_n')
        self.dctp_n = self.human_model.metabolites.get_by_id('dctp_n')
        self.dgtp_n = self.human_model.metabolites.get_by_id('dgtp_n')
        self.dttp_n = self.human_model.metabolites.get_by_id('dttp_n')

        # carbohydrate
        self.g6p_c = self.human_model.metabolites.get_by_id('g6p_c')

        # lipid
        self.chsterol_c = self.human_model.metabolites.get_by_id('chsterol_c')
        self.clpn_hs_c = self.human_model.metabolites.get_by_id('clpn_hs_c')
        self.pail_hs_c = self.human_model.metabolites.get_by_id('pail_hs_c')
        self.pchol_hs_c = self.human_model.metabolites.get_by_id('pchol_hs_c')
        self.pe_hs_c = self.human_model.metabolites.get_by_id('pe_hs_c')
        self.pglyc_hs_c = self.human_model.metabolites.get_by_id('pglyc_hs_c')
        self.ps_hs_c = self.human_model.metabolites.get_by_id('ps_hs_c')
        self.sphmyln_hs_c = self.human_model.metabolites.get_by_id('sphmyln_hs_c')

        # h2o = {'r': h2o_r, 'c': h2o_c, 'l': h2o_l, 'm': h2o_m, 'n': h2o_n, 'x': h2o_x}
        self.nmp_map = {'n': nmp_map_n, 'c': nmp_map_c}
        self.ntp_map = {'n': ntp_map_n, 'c': ntp_map_c}
