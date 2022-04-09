# Installation Instructions

These provide additional details and alternate installation instructions for setting up the human_me environment and qMINOS solver

## Step 1-2: Creating the Environment

We highly recommending setting up a [python virtual environment](https://packaging.python.org/guides/) instead of a conda environment, because it works better with the solver.

This can be done as follows:
```console
$python3 -m venv <env_name> #--python=python3.6.9
$source [env_name]/bin/activate
```

human_me can be installed using pip and PyPi:
```console
$pip install human_me
```

***Alternatively***, the human_me environment can be setup independent of PyPi using github and requirements.txt: 
```console
$git clone https://github.com/hmbaghdassarian/human_me.git
$pip install -r <path/to/human_me/>requirements.txt
```

or setup.py:
```console
$git clone https://github.com/hmbaghdassarian/human_me.git
$python <path/to/human_me/>setup.py install
```

If jupyter notebook does not load the virtual environment kernel, outside of the environment, try: 
```console
python3 -m ipykernel install --user --name=<env_name>
```

## Step 3: Setting up the qMINOS solver
qMINOS is a high precision LP solver necessary for the order-of-magnitude differences in ME Model coefficients.<br>
The QMINOS solver can be obtained for academic use from Prof. Michael Saunders at Stanford University.<br>
gfortran (>=4.6) is required for qMINOS <br>

&emsp;i) download the qminos file into a specified solver directory, which we refer to here as "solver_parent_directory".

&emsp;ii) the solver can be installed using the human_me Makefile as follows:

```console
$make -C <path/to/human_me/> install-qminos SOLVER_PATH=<path/to/solver/solver_parent_directory>
```

If SOLVER_PATH is not specified, it defaults to path/to/human_me/solver.

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