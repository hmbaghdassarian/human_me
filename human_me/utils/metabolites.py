#!/usr/bin/env python
# coding: utf-8

class MetaboliteBin:
    '''Stores all metabolites used in model building'''
    def __init__(self, me_input_model):
        """Init method for MetaboliteBin.

        Parameters
        ----------
        me_input_model : cobra.Model
            the corrected input metabolic model (as provided in preprocess.correct_inputs.correct_model)
        """
        self.human_model = me_input_model
        self.id_object_map = {m.id: m.copy() for m in self.human_model.metabolites} # copying dissociates metabolite from original cobrapy reactions/model

        # copying all metabolites makes them unassociated with the input model
        self.atp_n = self.id_object_map['atp_n']
        self.gtp_n = self.id_object_map['gtp_n']
        self.ppi_n = self.id_object_map['ppi_n']
        self.ppi_c = self.id_object_map['ppi_c']

        # mrna expression
        # processing variables------------------------------------------------------------

        self.pi_n = self.id_object_map['pi_n']
        self.h_n = self.id_object_map['h_n']
        self.h2o_n = self.id_object_map['h2o_n']
        gp = self.gtp_n.elements.copy()
        gp['O'] -= 6
        gp['P'] -= 2
        self.gp = gp
        self.amet_n = self.id_object_map['amet_n']
        self.ahcys_n = self.id_object_map['ahcys_n']
        self.adp_n = self.id_object_map['adp_n']

        # mrna degradation------------------------------------------------------------
        self.h_c = self.id_object_map['h_c']
        self.h2o_c = self.id_object_map['h2o_c']
        self.pi_c = self.id_object_map['pi_c']
        self.amet_c = self.id_object_map['amet_c']
        self.ahcys_c = self.id_object_map['ahcys_c']
        self.amp_c = self.id_object_map['amp_c']

        # rrna expression
        self.atp_c = self.id_object_map['atp_c']

        # protein expression

        # nucleus------------------------------------------------------------
        self.gdp_n = self.id_object_map['gdp_n']

        # secretory pathway------------------------------------------------------------
        self.o2_r = self.id_object_map['o2_r']
        self.h2o2_r = self.id_object_map['h2o2_r']

        self.hdca_r = self.id_object_map['hdca_r']
        self.gpi_hs_r = self.id_object_map['gpi_hs_r']
        self.balanced_gpi = {'C': 71, 'H': 138, 'N': 4, 'O': 41, 'P': 4}  # == dgpi_prot_hs_r.elements

        self.udpacgal_g = self.id_object_map['udpacgal_g']
        self.udpgal_g = self.id_object_map['udpgal_g']
        self.uacgam_g = self.id_object_map['uacgam_g']
        self.h_g = self.id_object_map['h_g']
        self.udp_g = self.id_object_map['udp_g']

        self.udpgal_r = self.id_object_map['udpgal_r']
        self.uacgam_r = self.id_object_map['uacgam_r']
        self.udpacgal_r = self.id_object_map['udpacgal_r']
        self.udp_r = self.id_object_map['udp_r']

        self.hdca_l = self.id_object_map['hdca_l']
        self.h_l = self.id_object_map['h_l']
        self.h2o_l = self.id_object_map['h2o_l']
        self.o2_l = self.id_object_map['o2_l']
        self.h2o2_l = self.id_object_map['h2o2_l']

        self.udpacgal_l = self.id_object_map['udpacgal_l']
        self.udp_l = self.id_object_map['udp_l']

        # mrna expression
        self.seq_metabolite_map = {self.id_object_map['utp_n']: 'U',
                            self.gtp_n: 'G',
                            self.id_object_map['ctp_n']: 'C',
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
        self.nmp_map_n = {'C': self.id_object_map['cmp_n'],
                    'U': self.id_object_map['ump_n'],
                    'G': self.id_object_map['gmp_n'],
                    'A': self.id_object_map['amp_n']}
        self.ntp_map_n = {v: k for k, v in self.seq_metabolite_map.items()}

        # mrna degradation------------------------------------------------------------
        self.gmp_c = self.id_object_map['gmp_c']
        self.nmp_map_c = {'C': self.id_object_map['cmp_c'],
                    'U': self.id_object_map['ump_c'],
                    'G': self.gmp_c,
                    'A': self.amp_c}
        self.gdp_c = self.id_object_map['gdp_c']
        self.ndp_map_c = {'C': self.id_object_map['cdp_c'],
                    'U': self.id_object_map['udp_c'],
                    'G': self.gdp_c,
                    'A': self.id_object_map['adp_c']}

        # rrna expression
        self.ntp_map_c = {'C': self.id_object_map['ctp_c'],
                    'U': self.id_object_map['utp_c'],
                    'G': self.id_object_map['gtp_c'],
                    'A': self.atp_c}

        # trna expression
        self.seq_amino_acid_map_c = {
            'A': self.id_object_map['ala_L_c'],
            'R': self.id_object_map['arg_L_c'],
            'N': self.id_object_map['asn_L_c'],
            'D': self.id_object_map['asp_L_c'],
            'C': self.id_object_map['cys_L_c'],
            'E': self.id_object_map['glu_L_c'],
            'Q': self.id_object_map['gln_L_c'],
            'G': self.id_object_map['gly_c'],
            'H': self.id_object_map['his_L_c'],
            'I': self.id_object_map['ile_L_c'],
            'L': self.id_object_map['leu_L_c'],
            'K': self.id_object_map['lys_L_c'],
            'M': self.id_object_map['met_L_c'],
            'F': self.id_object_map['phe_L_c'],
            'P': self.id_object_map['pro_L_c'],
            'S': self.id_object_map['ser_L_c'],
            'T': self.id_object_map['thr_L_c'],
            'W': self.id_object_map['trp_L_c'],
            'Y': self.id_object_map['tyr_L_c'],
            'V': self.id_object_map['val_L_c'],
        }

        self.seq_amino_acid_map_m = {aa_code: self.id_object_map['_'.join(aa_metabolite.id.split('_')[:-1]) + '_m']
                                for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_l = {aa_code: self.id_object_map['_'.join(aa_metabolite.id.split('_')[:-1]) + '_l']
                                for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_x = {aa_code: self.id_object_map['_'.join(aa_metabolite.id.split('_')[:-1]) + '_x']
                                for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_n = {aa_code: self.id_object_map['_'.join(aa_metabolite.id.split('_')[:-1]) + '_n']
                                for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
        self.seq_amino_acid_map_r = {aa_code: self.id_object_map['_'.join(aa_metabolite.id.split('_')[:-1]) + '_r']
                                for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}

        self.seq_amino_acid_map_compartments = {'c': self.seq_amino_acid_map_c, 'x': self.seq_amino_acid_map_x, 'r': self.seq_amino_acid_map_r,
                                        'm': self.seq_amino_acid_map_m, 'n': self.seq_amino_acid_map_n, 'l': self.seq_amino_acid_map_l}

        self.adp_c = self.ndp_map_c['A']

        self.atp_m = self.id_object_map['atp_m']
        self.adp_m = self.id_object_map['adp_m']
        self.h_m = self.id_object_map['h_m']
        self.pi_m = self.id_object_map['pi_m']
        self.h2o_m = self.id_object_map['h2o_m']
        self.h_i = self.id_object_map['h_i']

        self.h_x = self.id_object_map['h_x']
        self.h2o_x = self.id_object_map['h2o_x']
        self.pi_x = self.id_object_map['pi_x']
        self.atp_x = self.id_object_map['atp_x']
        self.adp_x = self.id_object_map['adp_x']

        self.h_r = self.id_object_map['h_r']
        self.h2o_r = self.id_object_map['h2o_r']
        self.pi_r = self.id_object_map['pi_r']
        self.atp_r = self.id_object_map['atp_r']
        self.adp_r = self.id_object_map['adp_r']

        self.pi_l = self.id_object_map['pi_l']
        self.atp_l = self.id_object_map['atp_l']
        self.adp_l = self.id_object_map['adp_l']

        self.atp_compartments = {'c': self.atp_c, 'm': self.atp_m, 'i': self.atp_m, 'x': self.atp_x, 'n': self.atp_n, 'r': self.atp_r, 'l': self.atp_l}
        self.adp_compartments = {'c': self.adp_c, 'm': self.adp_m, 'i': self.adp_m, 'x': self.adp_x, 'n': self.adp_n, 'r': self.adp_r, 'l': self.adp_l}
        self.h2o_compartments = {'c': self.h2o_c, 'm': self.h2o_m, 'i': self.h2o_m, 'x': self.h2o_x, 'n': self.h2o_n, 'r': self.h2o_r, 'l': self.h2o_l}
        self.pi_compartments = {'c': self.pi_c, 'm': self.pi_m, 'i': self.pi_m, 'x': self.pi_x, 'n': self.pi_n, 'r': self.pi_r, 'l': self.pi_l}
        self.h_compartments = {'c': self.h_c, 'm': self.h_m, 'i': self.h_i, 'x': self.h_x, 'n': self.h_n, 'r': self.h_r, 'l': self.h_l}

        self.datp_n = self.id_object_map['datp_n']
        self.dctp_n = self.id_object_map['dctp_n']
        self.dgtp_n = self.id_object_map['dgtp_n']
        self.dttp_n = self.id_object_map['dttp_n']

        # carbohydrate
        self.g6p_c = self.id_object_map['g6p_c']

        # lipid
        self.chsterol_c = self.id_object_map['chsterol_c']
        self.clpn_hs_c = self.id_object_map['clpn_hs_c']
        self.pail_hs_c = self.id_object_map['pail_hs_c']
        self.pchol_hs_c = self.id_object_map['pchol_hs_c']
        self.pe_hs_c = self.id_object_map['pe_hs_c']
        self.pglyc_hs_c = self.id_object_map['pglyc_hs_c']
        self.ps_hs_c = self.id_object_map['ps_hs_c']
        self.sphmyln_hs_c = self.id_object_map['sphmyln_hs_c']

        # h2o = {'r': h2o_r, 'c': h2o_c, 'l': h2o_l, 'm': h2o_m, 'n': h2o_n, 'x': h2o_x}
        self.nmp_map = {'n': self.nmp_map_n, 'c': self.nmp_map_c}
        self.ntp_map = {'n': self.ntp_map_n, 'c': self.ntp_map_c}


# class MetaboliteBin:
#     '''Stores all metabolites used in model building'''
#     def __init__(self, me_input_model):
#         """Init method for MetaboliteBin

#         Parameters
#         ----------
#         me_input_model : cobra.Model
#             the corrected input metabolic model (as provided in preprocess.correct_inputs.correct_model)
#         """

#         self.human_model = me_input_model
#         # copying all metabolites makes them unassociated with the input model
#         self.atp_n = self.id_object_map['atp_n']
#         self.gtp_n = self.id_object_map['gtp_n']
#         self.ppi_n = self.id_object_map['ppi_n']
#         self.ppi_c = self.id_object_map['ppi_c']

#         # mrna expression
#         # processing variables------------------------------------------------------------

#         self.pi_n = self.id_object_map['pi_n']
#         self.h_n = self.id_object_map['h_n']
#         self.h2o_n = self.id_object_map['h2o_n']
#         gp = self.gtp_n.elements.copy()
#         gp['O'] -= 6
#         gp['P'] -= 2
#         self.gp = gp
#         self.amet_n = self.id_object_map['amet_n']
#         self.ahcys_n = self.id_object_map['ahcys_n']
#         self.adp_n = self.id_object_map['adp_n']

#         # mrna degradation------------------------------------------------------------
#         self.h_c = self.id_object_map['h_c']
#         self.h2o_c = self.id_object_map['h2o_c']
#         self.pi_c = self.id_object_map['pi_c']
#         self.amet_c = self.id_object_map['amet_c']
#         self.ahcys_c = self.id_object_map['ahcys_c']
#         self.amp_c = self.id_object_map['amp_c']

#         # rrna expression
#         self.atp_c = self.id_object_map['atp_c']

#         # protein expression

#         # nucleus------------------------------------------------------------
#         self.gdp_n = self.id_object_map['gdp_n']

#         # secretory pathway------------------------------------------------------------
#         self.o2_r = self.id_object_map['o2_r']
#         self.h2o2_r = self.id_object_map['h2o2_r']

#         self.hdca_r = self.id_object_map['hdca_r']
#         self.gpi_hs_r = self.id_object_map['gpi_hs_r']
#         self.balanced_gpi = {'C': 71, 'H': 138, 'N': 4, 'O': 41, 'P': 4}  # == dgpi_prot_hs_r.elements

#         self.udpacgal_g = self.id_object_map['udpacgal_g']
#         self.udpgal_g = self.id_object_map['udpgal_g']
#         self.uacgam_g = self.id_object_map['uacgam_g']
#         self.h_g = self.id_object_map['h_g']
#         self.udp_g = self.id_object_map['udp_g']

#         self.udpgal_r = self.id_object_map['udpgal_r']
#         self.uacgam_r = self.id_object_map['uacgam_r']
#         self.udpacgal_r = self.id_object_map['udpacgal_r']
#         self.udp_r = self.id_object_map['udp_r']

#         self.hdca_l = self.id_object_map['hdca_l']
#         self.h_l = self.id_object_map['h_l']
#         self.h2o_l = self.id_object_map['h2o_l']
#         self.o2_l = self.id_object_map['o2_l']
#         self.h2o2_l = self.id_object_map['h2o2_l']

#         self.udpacgal_l = self.id_object_map['udpacgal_l']
#         self.udp_l = self.id_object_map['udp_l']

#         # mrna expression
#         self.seq_metabolite_map = {self.id_object_map['utp_n']: 'U',
#                             self.gtp_n: 'G',
#                             self.id_object_map['ctp_n']: 'C',
#                             self.atp_n: 'A'}

#         # RNA backbone elements
#         self.seq_element_map = dict()
#         for k, v in self.seq_metabolite_map.items():
#             elements = k.elements.copy()
#             elements['O'] = elements['O'] - 7  # lost from incoming ntp
#             elements['P'] = elements['P'] - 2  # lost from incoming ntp
#             elements['H'] = elements['H'] - 1  # lost from 3' end of growing strand
#             self.seq_element_map[v] = elements
#         # lariat degradataion------------------------------------------------------------
#         self.nmp_map_n = {'C': self.id_object_map['cmp_n'],
#                     'U': self.id_object_map['ump_n'],
#                     'G': self.id_object_map['gmp_n'],
#                     'A': self.id_object_map['amp_n']}
#         self.ntp_map_n = {v: k for k, v in self.seq_metabolite_map.items()}

#         # mrna degradation------------------------------------------------------------
#         self.gmp_c = self.id_object_map['gmp_c']
#         self.nmp_map_c = {'C': self.id_object_map['cmp_c'],
#                     'U': self.id_object_map['ump_c'],
#                     'G': self.gmp_c,
#                     'A': self.amp_c}
#         self.gdp_c = self.id_object_map['gdp_c']
#         self.ndp_map_c = {'C': self.id_object_map['cdp_c'],
#                     'U': self.id_object_map['udp_c'],
#                     'G': self.gdp_c,
#                     'A': self.id_object_map['adp_c']}

#         # rrna expression
#         self.ntp_map_c = {'C': self.id_object_map['ctp_c'],
#                     'U': self.id_object_map['utp_c'],
#                     'G': self.id_object_map['gtp_c'],
#                     'A': self.atp_c}

#         # trna expression
#         self.seq_amino_acid_map_c = {
#             'A': self.id_object_map['ala_L_c'],
#             'R': self.id_object_map['arg_L_c'],
#             'N': self.id_object_map['asn_L_c'],
#             'D': self.id_object_map['asp_L_c'],
#             'C': self.id_object_map['cys_L_c'],
#             'E': self.id_object_map['glu_L_c'],
#             'Q': self.id_object_map['gln_L_c'],
#             'G': self.id_object_map['gly_c'],
#             'H': self.id_object_map['his_L_c'],
#             'I': self.id_object_map['ile_L_c'],
#             'L': self.id_object_map['leu_L_c'],
#             'K': self.id_object_map['lys_L_c'],
#             'M': self.id_object_map['met_L_c'],
#             'F': self.id_object_map['phe_L_c'],
#             'P': self.id_object_map['pro_L_c'],
#             'S': self.id_object_map['ser_L_c'],
#             'T': self.id_object_map['thr_L_c'],
#             'W': self.id_object_map['trp_L_c'],
#             'Y': self.id_object_map['tyr_L_c'],
#             'V': self.id_object_map['val_L_c'],
#         }

#         self.seq_amino_acid_map_m = {aa_code: self.id_object_map['_'.join(aa_metabolite.id.split('_')[:-1]) + '_m']
#                                 for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
#         self.seq_amino_acid_map_l = {aa_code: self.id_object_map['_'.join(aa_metabolite.id.split('_')[:-1]) + '_l']
#                                 for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
#         self.seq_amino_acid_map_x = {aa_code: self.id_object_map['_'.join(aa_metabolite.id.split('_')[:-1]) + '_x']
#                                 for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
#         self.seq_amino_acid_map_n = {aa_code: self.id_object_map['_'.join(aa_metabolite.id.split('_')[:-1]) + '_n']
#                                 for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}
#         self.seq_amino_acid_map_r = {aa_code: self.id_object_map['_'.join(aa_metabolite.id.split('_')[:-1]) + '_r']
#                                 for aa_code, aa_metabolite in self.seq_amino_acid_map_c.items()}

#         self.seq_amino_acid_map_compartments = {'c': self.seq_amino_acid_map_c, 'x': self.seq_amino_acid_map_x, 'r': self.seq_amino_acid_map_r,
#                                         'm': self.seq_amino_acid_map_m, 'n': self.seq_amino_acid_map_n, 'l': self.seq_amino_acid_map_l}

#         self.adp_c = self.ndp_map_c['A']

#         self.atp_m = self.id_object_map['atp_m']
#         self.adp_m = self.id_object_map['adp_m']
#         self.h_m = self.id_object_map['h_m']
#         self.pi_m = self.id_object_map['pi_m']
#         self.h2o_m = self.id_object_map['h2o_m']
#         self.h_i = self.id_object_map['h_i']

#         self.h_x = self.id_object_map['h_x']
#         self.h2o_x = self.id_object_map['h2o_x']
#         self.pi_x = self.id_object_map['pi_x']
#         self.atp_x = self.id_object_map['atp_x']
#         self.adp_x = self.id_object_map['adp_x']

#         self.h_r = self.id_object_map['h_r']
#         self.h2o_r = self.id_object_map['h2o_r']
#         self.pi_r = self.id_object_map['pi_r']
#         self.atp_r = self.id_object_map['atp_r']
#         self.adp_r = self.id_object_map['adp_r']

#         self.pi_l = self.id_object_map['pi_l']
#         self.atp_l = self.id_object_map['atp_l']
#         self.adp_l = self.id_object_map['adp_l']

#         self.atp_compartments = {'c': self.atp_c, 'm': self.atp_m, 'i': self.atp_m, 'x': self.atp_x, 'n': self.atp_n, 'r': self.atp_r, 'l': self.atp_l}
#         self.adp_compartments = {'c': self.adp_c, 'm': self.adp_m, 'i': self.adp_m, 'x': self.adp_x, 'n': self.adp_n, 'r': self.adp_r, 'l': self.adp_l}
#         self.h2o_compartments = {'c': self.h2o_c, 'm': self.h2o_m, 'i': self.h2o_m, 'x': self.h2o_x, 'n': self.h2o_n, 'r': self.h2o_r, 'l': self.h2o_l}
#         self.pi_compartments = {'c': self.pi_c, 'm': self.pi_m, 'i': self.pi_m, 'x': self.pi_x, 'n': self.pi_n, 'r': self.pi_r, 'l': self.pi_l}
#         self.h_compartments = {'c': self.h_c, 'm': self.h_m, 'i': self.h_i, 'x': self.h_x, 'n': self.h_n, 'r': self.h_r, 'l': self.h_l}

#         self.datp_n = self.id_object_map['datp_n']
#         self.dctp_n = self.id_object_map['dctp_n']
#         self.dgtp_n = self.id_object_map['dgtp_n']
#         self.dttp_n = self.id_object_map['dttp_n']

#         # carbohydrate
#         self.g6p_c = self.id_object_map['g6p_c']

#         # lipid
#         self.chsterol_c = self.id_object_map['chsterol_c']
#         self.clpn_hs_c = self.id_object_map['clpn_hs_c']
#         self.pail_hs_c = self.id_object_map['pail_hs_c']
#         self.pchol_hs_c = self.id_object_map['pchol_hs_c']
#         self.pe_hs_c = self.id_object_map['pe_hs_c']
#         self.pglyc_hs_c = self.id_object_map['pglyc_hs_c']
#         self.ps_hs_c = self.id_object_map['ps_hs_c']
#         self.sphmyln_hs_c = self.id_object_map['sphmyln_hs_c']

#         # h2o = {'r': h2o_r, 'c': h2o_c, 'l': h2o_l, 'm': h2o_m, 'n': h2o_n, 'x': h2o_x}
#         self.nmp_map = {'n': self.nmp_map_n, 'c': self.nmp_map_c}
#         self.ntp_map = {'n': self.ntp_map_n, 'c': self.ntp_map_c}
