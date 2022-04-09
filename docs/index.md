# human_me

Python package to generate and analyze human ME Models.

## Protein-Specific Information Matrix (PSIM)

This section provides details on the PSIM. The PSIM is a dataframe which contains all protein-specific features needed to generate the expression module reactions in the ME-Model

### Data Sources
We provide a "gold-standard" PSIM (both in the build/ directory and the inputs/ directory with file name "psim_me.h5"). 
Sources:
* Sequences were generated form MANE Select, RefSeq Select, and APPRIS. 
* The number of exons was found using the Ensembl REST API for the specific transcript isoform. 
* The poly(A)-length was taken from https://doi.org/10.1038/s41592-019-0503-y. 
* Metadata associated to the secretory pathway was taken from human PSIM specified in  https://doi.org/10.1038/s41467-019-13867-y (https://github.com/LewisLabUCSD/MammalianSecretoryRecon/blob/master/JUPYTER_NOTEBOOKS/RECON2s_python/PSIM_HUMAN.tab)
* mrna degradation rates were taken from https://doi.org/10.1016/j.molcel.2014.03.017
* protein degradation rates were taken from https://doi.org/10.1021/pr101183k
* PTRs were taken from https://doi.org/10.15252/msb.20188513

### Formatting
We explain the formatting and values for all columns in the gold-standard PSIM. Any user-provided PSIM can be generated following this formatting (column names must match exactly). 

Legend:
* super required columns<sup>0</sup>: default values unavailable, user must provide
* required columns<sup>1</sup>: default values unavailable, but if not provided or incorrect, will fill in with the gold-standard PSIM values when available (otherwise will error out)
* semi-optional columns<sup>2</sup>: may or may not be required by pipeline
* optional columns<sup>3</sup>: used in pipeline, but default values can be used when not provided 
* optional columns, secretory pathway only<sup>4</sup>: if provided, these are only used for proteins that will be processed via the secretory pathway (ER, Golgi, exctracellular membrane, plasma membrane, lysosome).
* additional information<sup>5</sup>: not used in pipeline, but additional metadata associated with each gene - specific to Recon2.2 and gold-standard PSIM. User provided PSIM does not need these columns. 


The PSIM is read into the pipeline using the preprocess.correct_inputs.correct_psim function, which will check that all values make sense. If the fill_na argument is set to 'select', all NaN values in the user-provided PSIM are filled with values from the gold-standard PSIM when available, otherwise default. If the fill_na argument is set to 'default', all NaN values are replaced with a default value according to the expression.gene_expression.gene_information.GeneInformation class. 

Columns:
* HGNC_ID<sup>0</sup>: The gene ID in HGNC format (HGNC:####). There should be an entry for all genes that are included in the M_Model GPR and in non-machinery. 
    * Datatype: str
* PREMRNA_SEQ<sup>1</sup>: The gene premrna sequence. Requirements include that values can only include 'A', 'C', 'G', 'U', and the sequence length must be >= mrna sequence length.  
    * Datatype: string
    * Default value: Technically none, but preprocess.correct_inputs.correct_psim will fill incorrect values with the gold-standard PSIM values. Requirements include that values can only include 'A', 'C', 'G', 'U', the sequence length must be <= premrna sequence length, and the sequence length must be >= 3*protein sequence length.  
* MRNA_SEQ<sup>1</sup>: The gene mrna sequence (isoform specific). 
    * Datatype: str
    * Default value: Technically none, but preprocess.correct_inputs.correct_psim will fill incorrect values with the gold-standard PSIM values. 
* PROTEIN_SEQ<sup>1</sup>: The gene protein sequence (isoform specific). Requirements include values can only include one-letter amino-acid codes and the sequence length <= (mrna sequence length/3)
    * Datatype: str
    * Default value: Technically none, but preprocess.correct_inputs.correct_psim will fill incorrect values with the gold-standard PSIM values. 
* POLYA_LENGTH<sup>3</sup>: The length of the mature mRNA polyA tail. 
    * Datatype: int 
    * Default value: Randomly draws from a johnsonsu distribution
* N_EXONS<sup>3</sup>: The number of exons in the premrna (isoform specific). Use to estimate the number of introns (as # of exons - 1). 
    * Datatype: int
    * Default value: Estimated as (premrna sequence length)/6700
* TMD<sup>4</sup>: The number of transmembrane domains contained in the sequence.
    * Datatype: int
    * Default value: 0
* SP<sup>4</sup>: Whether the protein contains a secretory pathway signal peptide. This option is not currently implemented as all proteins destined for secretory pathway compartments are assumed to have a signal peptide (SP). In the future, this option will be used for non-canonical secretion
    * Datatype: bool
    * Default value: True for secretory pathway destined proteins, False otherwise
* DSB<sup>4</sup>: The number of disulfide bonds in the protein. 
    * Datatype: int
    * Default value: 0
* GPI<sup>4</sup>: Whether a GPI anchor is present in the protein. 0 if not present, 1 otherwise.  
    * Datatype: int
    * Default value: 0
* OG<sup>4</sup>: The number of utilized O-linked glycosylation sites in the protein. 
    * Datatype: int
    * Default value: 0
* NG<sup>4</sup>: The number of utilized N-linked glycosylation sites in the protein. 
    * Datatype: int
    * Default value: 0
* ALPHA_M<sup>3</sup>: The mrna degradation/turnover rate (hrs^-1). Used in calculating coupling constraints.
    * Datatype: float
    * Default value: 0.06 hrs^-1
* ALPHA_P<sup>3</sup>: The protein degradation/turnover rate (hrs^-1). Used in calculating coupling constraints.
    * Datatype: float
    * Default value: 0.02 hrs^-1
* PTR<sup>3</sup>: The protein to rna ratio, as described in https://doi.org/10.15252/msb.20188513. Used in calculating coupling constraints.
    * Datatype: float
    * Default value: 65163
* LOCATION<sup>2</sup>: The final location of the protein. Required for non-machinery, disregarded for machinery (pipeline infers location from the reaction compartments).  
    * Datatype: str, on of utils.paramaters.compartments.keys()
    
Additional Information in Gold-Standard PSIM:
* Machinery<sup>5</sup>: Whether a protein is considered machinery according to the full Recon2.2 ('Metabolic'), the GPRs for expression reactions ('Expression'), both ('Both'), or neither ('Non-Machinery'). 
* Source<sup>5</sup>: From which database the isoform sequences were attained. 
* Status<sup>5</sup>: 1 for entries that should work with the pipeline, 0 for entries that will cause an error in the pipeline. 
* The remaining columns<sup>5</sup> are various IDs for the gene: 'GENE_SYMBOL', 'GeneID', 'ENSG_ID', 'ENST_ID', 'ENSP_ID', 'REFT_ID', 'REFP_ID','UNIPROT_ID'. 
