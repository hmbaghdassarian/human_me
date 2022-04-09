# human_me
human_me is a python package to generate ME-Models from input context-extracted Recon2.2 M-Models

## Installation
For additional details and alternate install options, see the documentation's [installation instructions](https://hmbaghdassarian.github.io/human_me/install/). 
Requirements: gfortran (>=4.6)


1. Create a [python virtual environment](https://packaging.python.org/guides/) and activate this environment
2. Install human_me:
```console
$pip install human_me
```
3. Set up the qMINOS solver. 

⋅⋅⋅i) download the qminos file into a specified solver directory, which we refer to here as "solver_parent_directory".

⋅⋅⋅ii) the solver can be installed using the human_me Makefile as follows:

```console
$make -C <path/to/human_me/> install-qminos SOLVER_PATH=<path/to/solver/solver_parent_directory>
```

4. Download the files in the "data.zip" and "build_files.zip" folders. If you don't have a PSIM or cobrapy model to input, also download those from input files. If you want to express non_machinery, create that as a text file (list of HGNC IDs separated by \n)

## Downloadable Files
Download directories 1 and 2: 

1. prebuild: not needed to run; these are file inputs/outputs from analyses that helped generate the build files
2. build: all files used in pipeline to building the ME Model; can be downloaded here
3. input: the three files used as inputs to the ME Model building (detailed descriptions below); we have provided our recommended input files here

## Input File Descriptions
1. M_Model: a cobrapy metabolic model in sbml format. We highly recommend Recon2.2 or a context-specific metabolic model generated from Recon2.2, as this is the only model the pipeline has been tested on. Our "inputs" directory provides a version of Recon2.2 with minor modifications to work with the ME-Model building pipeline. Alternatively, you can use the preprocess.correct_inputs.correct_model function on your metabolic model to introduce these modifications. 
2. PSIM: see the [doumentation](https://hmbaghdassarian.github.io/human_me/) for details
3. Non-machinery: an optional list or text file of HGNC IDs for non-machinery proteins to be expressed by the ME-Model. We define machinery as proteins that are utilized in the reaction GPRs. In its current format, it is not possible for a protein to both be machinery and non-machinery (e.g., HGNC:23408 used to catalyze reaction "3HBCDm" in Recon2.2 will be expressed and transported to the mitochondria to catalyze this reaction in the ME-Model, but cannot also be specified as a non-machinery to be secreted).

## About
For instructions on how to use the human ME Model, visit the [documentation](https://hmbaghdassarian.github.io/human_me/) 
