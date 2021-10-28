# human_me
### human_me is a python package to generate ME-Models from input context-extracted Recon2.2 M-Models

Installation
------

1. Create a [python virtual environment](https://packaging.python.org/guides/) using the requirements.txt file:
```console
python3 -m venv <env_name> --python=python3.6.9
source [env_name]/bin/activate
pip install -r path/to/human_me/requirements.txt
```

Note: Conda environments do not work well with solvemepy. 

2. With the environment activated, setup the QMINOS solver as specified in the installation instructions for [solvemepy](https://github.com/SBRG/solvemepy). The QMINOS solver can be obtained for academic use from Prof. Michael Saunders at Stanford University. Make sure gfortran is available to your system, as it is needed for running the solver. 

Note: Disregard the requirements that are delineated in solvmepy's README, with the exception of gfortran (#4). The remainder should have been appropriately installed with the requirements.txt file or are unecessary for human_me. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Step 2a: Untar the QMINOS solver and follow Step 1 of solveme installation guide:

```console
tar -xvf qminos.tar.gz #tar file from Prof. Michael Saunders
cd qminos/
cp Makefile.defs minos56/
cp Makefile.defs qminos56/
cd minos56
make clean
make
cd ../qminos56/
make clean
make
```

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Step 2b: Having exited the qminos directory, clone the solveme github and follow Step 2-3 of the solveme installation guide: 

```console
git clone https://github.com/SBRG/solvemepy.git #@2a2c9c098d5bad957ef41637955fe338a31bac4c
cd solvemepy
cp path/to/qminos/minos56/lib/libminos.a ./
cp path/to/qminos/qminos56/lib/libqminos.a ./
python setup.py develop
```

3. Download the files in the "data.zip" and "build_files.zip" folders. If you don't have a PSIM or cobrapy model to input, also download those from input files. If you want to express non_machinery, create that as a text file (list of HGNC IDs separated by \n)
4. See the Guide for getting started

If jupyter notebook does not load the virtual environment kernel, outside of the environment, try: 
```console
python3 -m ipykernel install --user --name=<env_name>
```

Downloadable Files
------
Download directories 1 and 2: 

0. prebuild: not needed to run; these are file inputs/outputs from analyses that helped generate the build files
1. build: all files used in pipeline to building the ME Model; can be downloaded here
2. input: the three files used as inputs to the ME Model building (detailed descriptions below); we have provided our recommended input files here

Input File Descriptions
------
1. M_Model: a cobrapy metabolic model in sbml format. We highly recommend Recon2.2 or a context-specific metabolic model generated from Recon2.2, as this is the only model the pipeline has been tested on. Our "inputs" directory provides a version of Recon2.2 with minor modifications to work with the ME-Model building pipeline. Alternatively, you can use the preprocess.correct_inputs.correct_model function on your metabolic model to introduce these modifications. 
2. PSIM: see the PSIM README for details
3. Non-machinery: an optional list or text file of HGNC IDs for non-machinery proteins to be expressed by the ME-Model. We define machinery as proteins that are utilized in the reaction GPRs. In its current format, it is not possible for a protein to both be machinery and non-machinery (e.g., HGNC:23408 used to catalyze reaction "3HBCDm" in Recon2.2 will be expressed and transported to the mitochondria to catalyze this reaction in the ME-Model, but cannot also be specified as a non-machinery to be secreted). 






