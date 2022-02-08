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
        # copying all metabolites makes them unassociated with the input model
        self.atp_n = self.human_model.metabolites.get_by_id('atp_n').copy()
        self.gtp_n = self.human_model.metabolites.get_by_id('gtp_n').copy()
        self.ppi_n = self.human_model.metabolites.get_by_id('ppi_n').copy()
        self.ppi_c = self.human_model.metabolites.get_by_id('ppi_c').copy()

        # mrna expression
        # processing variables------------------------------------------------------------

        self.pi_n = self.human_model.metabolites.get_by_id('pi_n').copy()
        self.h_n = self.human_model.metabolites.get_by_id('h_n').copy()
        self.h2o_n = self.human_model.metabolites.get_by_id('h2o_n').copy()
        gp = self.gtp_n.elements.copy()
        gp['O'] -= 6
        gp['P'] -= 2
        self.gp = gp
        self.amet_n = self.human_model.metabolites.get_by_id('amet_n').copy()
        self.ahcys_n = self.human_model.metabolites.get_by_id('ahcys_n').copy()
        self.adp_n = self.human_model.metabolites.get_by_id('adp_n').copy()

        # mrna degradation------------------------------------------------------------
        self.h_c = self.human_model.metabolites.get_by_id('h_c').copy()
        self.h2o_c = self.human_model.metabolites.get_by_id('h2o_c').copy()
        self.pi_c = self.human_model.metabolites.get_by_id('pi_c').copy()
        self.amet_c = self.human_model.metabolites.get_by_id('amet_c').copy()
        self.ahcys_c = self.human_model.metabolites.get_by_id('ahcys_c').copy()
        self.amp_c = self.human_model.metabolites.get_by_id('amp_c').copy()

        # rrna expression
        self.atp_c = self.human_model.metabolites.get_by_id('atp_c').copy()

        # protein expression

        # nucleus------------------------------------------------------------
        self.gdp_n = self.human_model.metabolites.get_by_id('gdp_n').copy()

        # secretory pathway------------------------------------------------------------
        self.o2_r = self.human_model.metabolites.get_by_id('o2_r').copy()
        self.h2o2_r = self.human_model.metabolites.get_by_id('h2o2_r').copy()

        self.hdca_r = self.human_model.metabolites.get_by_id('hdca_r').copy()
        self.gpi_hs_r = self.human_model.metabolites.get_by_id('gpi_hs_r').copy()
        self.balanced_gpi = {'C': 71, 'H': 138, 'N': 4, 'O': 41, 'P': 4}  # == dgpi_prot_hs_r.elements

        self.udpacgal_g = self.human_model.metabolites.get_by_id('udpacgal_g').copy()
        self.udpgal_g = self.human_model.metabolites.get_by_id('udpgal_g').copy()
        self.uacgam_g = self.human_model.metabolites.get_by_id('uacgam_g').copy()
        self.h_g = self.human_model.metabolites.get_by_id('h_g').copy()
        self.udp_g = self.human_model.metabolites.get_by_id('udp_g').copy()

        self.udpgal_r = self.human_model.metabolites.get_by_id('udpgal_r').copy()
        self.uacgam_r = self.human_model.metabolites.get_by_id('uacgam_r').copy()
        self.udpacgal_r = self.human_model.metabolites.get_by_id('udpacgal_r').copy()
        self.udp_r = self.human_model.metabolites.get_by_id('udp_r').copy()

        self.hdca_l = self.human_model.metabolites.get_by_id('hdca_l').copy()
        self.h_l = self.human_model.metabolites.get_by_id('h_l').copy()
        self.h2o_l = self.human_model.metabolites.get_by_id('h2o_l').copy()
        self.o2_l = self.human_model.metabolites.get_by_id('o2_l').copy()
        self.h2o2_l = self.human_model.metabolites.get_by_id('h2o2_l').copy()

        self.udpacgal_l = self.human_model.metabolites.get_by_id('udpacgal_l').copy()
        self.udp_l = self.human_model.metabolites.get_by_id('udp_l').copy()

        # mrna expression
        self.seq_metabolite_map = {self.human_model.metabolites.get_by_id('utp_n').copy(): 'U',
                            self.gtp_n: 'G',
                            self.human_model.metabolites.get_by_id('ctp_n').copy(): 'C',
                            self.atp_n: 'A'}

        # RNA backbone elements
        self.seq_element_map = dict()
        for k, v in self.seq_metabolite_map.items():
            elements = k.elements.copy()
            elements['O'] = elements['O'] - 7  # lost from incoming ntp
            elements['P'] = elements['P'] - 2  # lost from incoming ntp
            elements['H'] = elements['H'] - 1  # lost from 3' end of growing strand
            self.seq_element_map[v] = elements
        # lariat degradataion------------------------------------------------------------
        self.nmp_map_n = {'C': self.human_model.metabolites.get_by_id('cmp_n').copy(),
                    'U': self.human_model.metabolites.get_by_id('ump_n').copy(),
                    'G': self.human_model.metabolites.get_by_id('gmp_n').copy(),
                    'A': self.human_model.metabolites.get_by_id('amp_n').copy()}
        self.ntp_map_n = {v: k for k, v in self.seq_metabolite_map.items()}

        # mrna degradation------------------------------------------------------------
        self.gmp_c = self.human_model.metabolites.get_by_id('gmp_c').copy()
        self.nmp_map_c = {'C': self.human_model.metabolites.get_by_id('cmp_c').copy(),
                    'U': self.human_model.metabolites.get_by_id('ump_c').copy(),
                    'G': self.gmp_c,
                    'A': self.amp_c}
        self.gdp_c = self.human_model.metabolites.get_by_id('gdp_c').copy()
        self.ndp_map_c = {'C': self.human_model.metabolites.get_by_id('cdp_c').copy(),
                    'U': self.human_model.metabolites.get_by_id('udp_c').copy(),
                    'G': self.gdp_c,
                    'A': self.human_model.metabolites.get_by_id('adp_c').copy()}

        # rrna expression
        self.ntp_map_c = {'C': self.human_model.metabolites.get_by_id('ctp_c').copy(),
                    'U': self.human_model.metabolites.get_by_id('utp_c').copy(),
                    'G': self.human_model.metabolites.get_by_id('gtp_c').copy(),
                    'A': self.atp_c}

        # trna expression
        self.seq_amino_acid_map_c = {
            'A': self.human_model.metabolites.get_by_id('ala_L_c').copy(),
            'R': self.human_model.metabolites.get_by_id('arg_L_c').copy(),
            'N': self.human_model.metabolites.get_by_id('asn_L_c').copy(),
            'D': self.human_model.metabolites.get_by_id('asp_L_c').copy(),
            'C': self.human_model.metabolites.get_by_id('cys_L_c').copy(),
            'E': self.human_model.metabolites.get_by_id('glu_L_c').copy(),
            'Q': self.human_model.metabolites.get_by_id('gln_L_c').copy(),
            'G': self.human_model.metabolites.get_by_id('gly_c').copy(),
            'H': self.human_model.metabolites.get_by_id('his_L_c').copy(),
            'I': self.human_model.metabolites.get_by_id('ile_L_c').copy(),
            'L': self.human_model.metabolites.get_by_id('leu_L_c').copy(),
            'K': self.human_model.metabolites.get_by_id('lys_L_c').copy(),
            'M': self.human_model.metabolites.get_by_id('met_L_c').copy(),
            'F': self.human_model.metabolites.get_by_id('phe_L_c').copy(),
            'P': self.human_model.metabolites.get_by_id('pro_L_c').copy(),
            'S': self.human_model.metabolites.get_by_id('ser_L_c').copy(),
            'T': self.human_model.metabolites.get_by_id('thr_L_c').copy(),
            'W': self.human_model.metabolites.get_by_id('trp_L_c').copy(),
            'Y': self.human_model.metabolites.get_by_id('tyr_L_c').copy(),
            'V': self.human_model.metabolites.get_by_id('val_L_c').copy(),
        }

        self.seq_amino_acid_map_m = {aa_code: self.human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_m').copy()
                                for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_l = {aa_code: self.human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_l').copy()
                                for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_x = {aa_code: self.human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_x').copy()
                                for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_n = {aa_code: self.human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_n').copy()
                                for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_r = {aa_code: self.human_model.metabolites.get_by_id('_'.join(aa_metabolite.id.split('_')[:-1]) + '_r').copy()
                                for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}

        self.seq_amino_acid_map_compartments = {'c': self.seq_amino_acid_map_c, 'x': self.seq_amino_acid_map_x, 'r': self.seq_amino_acid_map_r,
                                        'm': self.seq_amino_acid_map_m, 'n': self.seq_amino_acid_map_n, 'l': self.seq_amino_acid_map_l}

        self.adp_c = self.ndp_map_c['A']

        self.atp_m = self.human_model.metabolites.get_by_id('atp_m').copy()
        self.adp_m = self.human_model.metabolites.get_by_id('adp_m').copy()
        self.h_m = self.human_model.metabolites.get_by_id('h_m').copy()
        self.pi_m = self.human_model.metabolites.get_by_id('pi_m').copy()
        self.h2o_m = self.human_model.metabolites.get_by_id('h2o_m').copy()
        self.h_i = self.human_model.metabolites.get_by_id('h_i').copy()

        self.h_x = self.human_model.metabolites.get_by_id('h_x').copy()
        self.h2o_x = self.human_model.metabolites.get_by_id('h2o_x').copy()
        self.pi_x = self.human_model.metabolites.get_by_id('pi_x').copy()
        self.atp_x = self.human_model.metabolites.get_by_id('atp_x').copy()
        self.adp_x = self.human_model.metabolites.get_by_id('adp_x').copy()

        self.h_r = self.human_model.metabolites.get_by_id('h_r').copy()
        self.h2o_r = self.human_model.metabolites.get_by_id('h2o_r').copy()
        self.pi_r = self.human_model.metabolites.get_by_id('pi_r').copy()
        self.atp_r = self.human_model.metabolites.get_by_id('atp_r').copy()
        self.adp_r = self.human_model.metabolites.get_by_id('adp_r').copy()

        self.pi_l = self.human_model.metabolites.get_by_id('pi_l').copy()
        self.atp_l = self.human_model.metabolites.get_by_id('atp_l').copy()
        self.adp_l = self.human_model.metabolites.get_by_id('adp_l').copy()

        self.atp_compartments = {'c': self.atp_c, 'm': self.atp_m, 'i': self.atp_m, 'x': self.atp_x, 'n': self.atp_n, 'r': self.atp_r, 'l': self.atp_l}
        self.adp_compartments = {'c': self.adp_c, 'm': self.adp_m, 'i': self.adp_m, 'x': self.adp_x, 'n': self.adp_n, 'r': self.adp_r, 'l': self.adp_l}
        self.h2o_compartments = {'c': self.h2o_c, 'm': self.h2o_m, 'i': self.h2o_m, 'x': self.h2o_x, 'n': self.h2o_n, 'r': self.h2o_r, 'l': self.h2o_l}
        self.pi_compartments = {'c': self.pi_c, 'm': self.pi_m, 'i': self.pi_m, 'x': self.pi_x, 'n': self.pi_n, 'r': self.pi_r, 'l': self.pi_l}
        self.h_compartments = {'c': self.h_c, 'm': self.h_m, 'i': self.h_i, 'x': self.h_x, 'n': self.h_n, 'r': self.h_r, 'l': self.h_l}

        self.datp_n = self.human_model.metabolites.get_by_id('datp_n').copy()
        self.dctp_n = self.human_model.metabolites.get_by_id('dctp_n').copy()
        self.dgtp_n = self.human_model.metabolites.get_by_id('dgtp_n').copy()
        self.dttp_n = self.human_model.metabolites.get_by_id('dttp_n').copy()

        # carbohydrate
        self.g6p_c = self.human_model.metabolites.get_by_id('g6p_c').copy()

        # lipid
        self.chsterol_c = self.human_model.metabolites.get_by_id('chsterol_c').copy()
        self.clpn_hs_c = self.human_model.metabolites.get_by_id('clpn_hs_c').copy()
        self.pail_hs_c = self.human_model.metabolites.get_by_id('pail_hs_c').copy()
        self.pchol_hs_c = self.human_model.metabolites.get_by_id('pchol_hs_c').copy()
        self.pe_hs_c = self.human_model.metabolites.get_by_id('pe_hs_c').copy()
        self.pglyc_hs_c = self.human_model.metabolites.get_by_id('pglyc_hs_c').copy()
        self.ps_hs_c = self.human_model.metabolites.get_by_id('ps_hs_c').copy()
        self.sphmyln_hs_c = self.human_model.metabolites.get_by_id('sphmyln_hs_c').copy()

        # h2o = {'r': h2o_r, 'c': h2o_c, 'l': h2o_l, 'm': h2o_m, 'n': h2o_n, 'x': h2o_x}
        self.nmp_map = {'n': self.nmp_map_n, 'c': self.nmp_map_c}
        self.ntp_map = {'n': self.ntp_map_n, 'c': self.ntp_map_c}
