# Installation Instructions

These provide additional details and alternate installation instructions for setting up the human_me environment and qMINOS solver

## Step 1: Creating the Environment

We highly recommending setting up a [python virtual environment](https://docs.python.org/3/library/venv.html) instead of a conda environment, because it works better with the solver.

This can be done as follows (example with python 3.9):
```console
python3.9 -m venv <env_name>
source <env_name>/bin/activate
```

All remaining steps should be implemented with the virtual environment activated. Make sure this virtual environment has Python 3.8-3.9, setuptools >= 65.4, and pip >= 21. Once the environment is created, if setuptools and pip do not meet these version requirements, you can run the following command:

```console
pip install --upgrade pip setuptools
```

## Step 2: Install the human_me package

human_me can be installed using pip and PyPi:
```console
pip install human_me
```

***Alternatively***, the human_me environment can be setup independent of PyPi using github and setup.py: 
```console
git clone https://github.com/hmbaghdassarian/human_me.git
cd human_me
pip install .

# delete the cloned repo so that the package path is appropriately read:
cd ..
rm -rf human_me
```

If you want to install jupyter notebook: 
```console
pip install human_me[interactive]
python3 -m ipykernel install --user --name=<env_name>
```

For steps 3-4 below, we need to get the path to the package:
```console
PACKAGE_PATH=$(python -c "import human_me; print(human_me.__path__[0])")/
```

## Step 3: Setting up the qMINOS solver
qMINOS is a high precision LP solver necessary for the order-of-magnitude differences in ME Model coefficients.<br>
The QMINOS solver can be obtained for academic use from Prof. Michael Saunders at Stanford University.<br>
gfortran (>=4.6) is required for qMINOS <br>

&emsp;i) download the qminos file into a specified directory, which we refer to here as "solver_parent_directory".

&emsp;ii) the solver can be installed using the human_me Makefile as follows:

```console
make -C $PACKAGE_PATH install-qminos SOLVER_PATH=<path/to/solver_parent_directory>
```

If SOLVER_PATH is not specified, it defaults to human_me/solver.

---

***Alternatively***, instead of using make, you can set up qminos manually as specified in the installation instructions for [solvemepy](https://github.com/SBRG/solvemepy). 

This will require getting both qminos and solvemepy. The solver parent directory specified above stores both qminos and solvemepy, but they do not need to be stored in the same directory.

You can disregard the requirements that are delineated in solvmepy's README, with the exception of gfortran (#4). The remainder should have been appropriately installed with installation of human_me, if necessary.

Ensure that all this is done with the virtual environment activated.

&emsp;i) Untar the QMINOS solver and follow Step 1 of solveme installation guide:

```console
tar -xvf qminos.tar.gz #tar file from Prof. Michael Saunders
cd <path/to/solver/solver_parent_directory/>qminos1114/
cp Makefile.defs minos56/
cp Makefile.defs qminos56/
cd minos56
make clean
make
cd ../qminos56/
make clean
make
```

&emsp;ii) Having exited the qminos directory, clone the solveme github and follow Step 2-3 of the solveme installation guide: 

```console
git clone https://github.com/SBRG/solvemepy.git #@2a2c9c098d5bad957ef41637955fe338a31bac4c
cd <path/to/solvemepy/>
cp <path/to/solver/solver_parent_directory/>qminos1114/minos56/lib/libminos.a ./
cp <path/to/solver/solver_parent_directory/>qminos1114/qminos56/lib/libquadminos.a ./
python setup.py develop
```

## Step 4: Downloading supporting data files

Data for human_me is structured as follows:
```
data
└───build
│      **psim_gold.h5
│      **recon2_2_only_psim.csv
└───inputs
│      **recon2_2.xml
│      psim_user
│      non_machinery
└───*prebuild
```
** = default, downloaded by make build-data

1. inputs: files used as inputs to the ME Model building; the defaul M-Model (Recon2.2) is downloaded. 
   Other main input files are briefly discussed below, and extensively discussed in the tutorials and API documentation. 
   While they do not need to be stored in /data/inputs/, it is a useful organization structure.
2. build: all files used in pipeline to building the ME Model; can be downloaded here
3. prebuild: not needed to run; these are file inputs/outputs from analyses that helped generate the build files


Small files are stored on Github in the [human_me_data repository](https://github.com/hmbaghdassarian/human_me_data) and accessed directly by the package. 
Larger files need to be downloaded directly from public Google Drive files. This is done using the make build-files command as follows:

```console
make -C $PACKAGE_PATH build-data DATA_DIR=</desired/local_data/directory>
```

Downloading the build files will take ~30 min, as the PSIM file is > 10 Gb.

If you also want to download the prebuild files, run the following command:
```console
make -C $PACKAGE_PATH build-data DATA_DIR=</desired/local_data/directory> PREBUILD=1
```
The user-specificed local data directory is stored in a human_me/data/data.ini config file.  

### Input File Descriptions
1. M_Model: a cobrapy metabolic model in sbml format. We highly recommend Recon2.2 or a context-specific metabolic model generated from Recon2.2, as this is the only model the pipeline has been tested on. Our "inputs" directory provides a version of Recon2.2 with minor modifications to work with the ME-Model building pipeline. Alternatively, you can use the preprocess.correct_inputs.correct_model function on your metabolic model to introduce these modifications. 
2. PSIM: see the [doumentation](https://hmbaghdassarian.github.io/human_me/) for details. If a user does not provide a PSIM, either directly into the function as a datatable or via data/inputs/psim_user, the default input PSIM is data/build/psim_gold.h5 (the gold-standard PSIM). 
