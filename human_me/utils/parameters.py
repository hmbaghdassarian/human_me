import logging
import os
import sys

import cobra
import numpy as np
import pandas as pd
from sympy.parsing.sympy_parser import parse_expr

from human_me.utils.load_environmental_variables import (build_files_path,
                                                         processed_data_path)

logging.basicConfig()
logger = logging.getLogger(cobra.__name__)
logger.setLevel(logging.CRITICAL)

class HiddenPrints:
    '''Supress package print messages.'''
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

with HiddenPrints():
    human_model = cobra.io.read_sbml_model(processed_data_path + 'corrected_model.xml')

psim_me = pd.read_hdf(processed_data_path + 'corrected_psim.h5', key='corrected')

psim_me['SP'] = psim_me['SP'].apply(lambda x: bool(x))
human_model = cobra.io.read_sbml_model(processed_data_path + 'corrected_model.xml')

mu = parse_expr('mu')

compartments = {'c': 'cytoplasm', 'l': 'lysosome', 'r': 'endoplasmic reticulum', 'e': 'extracellular space',
                'm': 'mitochondrion',
                'g': 'Golgi apparatus', 'n': 'nucleus', 'b': 'boundary', 'i': 'mitochondrial intermembrane space',
                'x': 'peroxisome', 'pm': 'plasma membrane'}

allowed_ptms = {'dsb': 'disulfide bond formation', 'gpi': 'GPI Anchor', 'og': 'O-linked glycosylation'}  # ,
# 'ng': 'N-linked glycosylation'}

amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

allowed_trna_modifications = {} 
# number_BiP = len(gene_info.protein_seq)/40

# universal variables and inputs

RATE_INTRON = 10 / 67000  # 10 introns / 67 kbp (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5199132/)
# L_polyA_n = 250 # https://www.nature.com/articles/s41592-019-0503-y
N_UB = 4  # see this - no. of ubiquitins to add to protein

TRANSPORT_TRANSLOCATION_ATP_COST = 0.5  # 1 ATP/2 residues
PROTEOLYSIS_TRANSLOCATION_ATP_COST = 0.5  # 1 ATP/2 residues

PTT_LENGTH = 160  # amino acid length greater than which co-translatioanl translocatoin occurs rather than
# post-translational
NUCLEAR_DIFFUSION_LIMIT = 40  # 40 kDA and less proteins diffuse through nucleus
L_SP = 22  # secretory pathway signal peptide degradation
K_V = 0.7  # secretory pathway vesicle coat coefficients

MEMBRANE_DIFFUSION_LIMIT = 504  # 504 Da includes ATP, uncharged molecules at this diffusion limit are passive, no dummy

# coupling parameters

# enzyme
KEFF_MEDIAN = 3.983 * 3600  # units: hr^-1 (3.983 in s^-1)

# central dogma
ALPHA_M_MEDIAN = 0.06108233261605428  # units: hours (Gregersen et al ) median value
ALPHA_P_MEDIAN = 0.019808138247250934  # units: hours ^-1 (Cambridge et al 2011 + Li et al 2021) median value
PTR_MEDIAN = 65162.83940608428  # (Eraslan et al 2019) median value

ptr = pd.read_csv(build_files_path + 'PTR_Gagneur_processed.tsv', sep='\t', index_col=0)
# don't groupby hgnc ID median, because if tissue option is used, can include unmapped ids in calculation
ptr.drop(columns=['ENSG_ID'], inplace=True)
ptr.columns = pd.Series(ptr.columns).apply(lambda x: x.split('_')[0] if '_PTR' in x else x).tolist()

alpha_p = pd.read_csv(build_files_path + 'protein_turnover.csv', index_col=0)
alpha_p = alpha_p.groupby(alpha_p.HGNC_ID).median().kdeg  # get median across cell lines

alpha_m = pd.read_csv(build_files_path + 'Gregersen_mrna_turnover_processed.tsv', sep='\t', index_col=0)
alpha_m = alpha_m.groupby(alpha_m.HGNC_ID).median().median_turnover  # have true median stored above

turnover = {'alpha_m': alpha_m, 'alpha_p': alpha_p,
            'alpha_m_median': ALPHA_M_MEDIAN, 'alpha_p_median': ALPHA_P_MEDIAN}

# ribosome
RNA_DEGRADATION_CONSTANT = np.log(2) / 72  # bioid 108025
SINGLE_UB_SEQ = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
# ribosomal_degradation_rate = np.log(2)/300 #bioid 110053 # unused

# biomass

# constant fractions
DNA_FRAC = 0.014
CARB_FRAC = 0.071
LIPID_FRAC = 0.097
# other_frac = 0.054

UNMODELED_PROTEIN_FRAC = 1 - 0.12041534186261499
