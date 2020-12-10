# human_me
### human_me is a python package to generate ME-Models from input context-extracted Recon2.2 M-Models

Installation
------


1v1. DEPRECATED - Create conda environment using environment.yml file

```console
conda create --name <env_name> --file=environment.yml
conda activate <env_name>
```

1v2. Create a [python virtual environment](https://packaging.python.org/guides/) using the requirements.txt file:
```console
python3.6.9 -m venv <env_name>
source [env_name]/bin/activate
pip install -r path/to/human_me/requirements.txt
```

Note: Conda environments do not work well with solvemepy. 

2. With the environment activated, setup the QMINOS solver as specified in the installation instructions for [solvemepy](https://github.com/SBRG/solvemepy). For step 2, you must first clone solvemepy. The QMINOS solver can be obtained for academic use from Prof. Michael Saunders at Stanford University. Make sure gfortran is available to your system, as it is needed for running the solver. 

Note: Disregard the requirements that are delineated in solvmepy's README, with the exception of gfortran (#4). The remainder should have been appropriately installed with the requirements.txt file or are unecessary. 

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

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Step 2b: Having exitted the qminos directory, clone the solveme github and follow Step 2-3 of the solveme installation guide: 

```console
git clone https://github.com/SBRG/solvemepy.git #@0d0ebca585c61ed0f4a559a11ef38706620d7c22
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

