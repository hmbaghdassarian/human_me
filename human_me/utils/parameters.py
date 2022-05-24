import numpy as np
import pandas as pd
from sympy.parsing.sympy_parser import parse_expr

from human_me.data.file_paths import build_files_url

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
N_UB = 4  # no. of ubiquitins to add to protein

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

ptr = pd.read_csv(build_files_url + 'PTR_Gagneur_processed.tsv', sep='\t', index_col=0)
# don't groupby hgnc ID median, because if tissue option is used, can include unmapped ids in calculation
ptr.drop(columns=['ENSG_ID'], inplace=True)
ptr.columns = pd.Series(ptr.columns).apply(lambda x: x.split('_')[0] if '_PTR' in x else x).tolist()

alpha_p = pd.read_csv(build_files_url + 'protein_turnover.csv', index_col=0)
alpha_p = alpha_p.groupby(alpha_p.HGNC_ID).median().kdeg  # get median across cell lines

alpha_m = pd.read_csv(build_files_url + 'Gregersen_mrna_turnover_processed.tsv', sep='\t', index_col=0)
alpha_m = alpha_m.groupby(alpha_m.HGNC_ID).median().median_turnover  # have true median stored above

turnover = {'alpha_m': alpha_m, 'alpha_p': alpha_p,
            'alpha_m_median': ALPHA_M_MEDIAN, 'alpha_p_median': ALPHA_P_MEDIAN}

# ribosome
RNA_DEGRADATION_CONSTANT = np.log(2) / 72  # bioid 108025
SINGLE_UB_SEQ = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
# ribosomal_degradation_rate = np.log(2)/300 #bioid 110053 # unused
RIBOSOME_TRANSLATION_RATE = 18000 # 5aa/sec/ribosome --> 18000 aa/hr/ribosome (BioID 104598, 107783, https://doi.org/10.1093/nar/gkaa1103)

# biomass
class BiomassParameters:
    """Stores Biomass reactions' asssociated coefficients and mass fractions in one object."""
    def __init__(self):
        """Create all the parameters."""
        self.mass_fraction = {'DNA': 0.014,
                              'carbohydrate': 0.071,
                              'lipid': 0.097,
                              'other': 0.054}
        # metabolite ID : metabolite coefficient (excluding biomass component) from Recon2.2 metabolic model (will be scaled by mass fraction)
        self.coefficients = {'DNA': {'datp_n': -0.941642857142857,
                                     'dctp_n': -0.674428571428572,
                                     'dgtp_n': -0.707,
                                     'dttp_n': -0.935071428571429,
                                     'ppi_n': 3.2581428571428583},
                             'carbohydrate': {'g6p_c': -3.87591549295775}, 
                             'lipid': {'chsterol_c': -0.09580101814936463,
                                        'clpn_hs_c': -0.05474478062767954,
                                        'pail_hs_c': -0.10948486535720975,
                                        'pchol_hs_c': -0.7253284281824849,
                                        'pe_hs_c': -0.26003066413425396,
                                        'pglyc_hs_c': -0.01368384720784515,
                                        'sphmyln_hs_c': -0.08211247504336976,
                                        'ps_hs_c': -0.02737239031383982}
                            }
biomass_parameters = BiomassParameters()
#TODO: implement an unmodeled protein fraction 
UNMODELED_PROTEIN_FRAC = 0 #1 - 0.12041534186261499 
