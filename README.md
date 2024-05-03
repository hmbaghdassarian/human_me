# human_me
human_me is a python package to generate ME-Models from input context-extracted Recon2.2 M-Models

## Installation
For additional details and alternate install options, see the documentation's [installation instructions](https://hmbaghdassarian.github.io/human_me/install/). <br>
Requirements: gfortran (>=4.6)


1. Create a [python virtual environment](https://docs.python.org/3/library/venv.html) and activate this environment. Make sure this virtual environment has Python 3.8-3.9, setuptools >= 65.4, and pip >=21. 
2. Install human_me:
```console
pip install human_me
```
3. Set up the qMINOS solver. You will need the qminos file, which can be obtained for academic use from Prof. Michael Saunders at Stanford University.

&emsp;i) download the qminos file into a specified solver directory, which we refer to here as "solver_parent_directory".

&emsp;ii) the solver can be installed using the human_me Makefile as follows:

```console
PACKAGE_PATH=$(python -c "import human_me; print(human_me.__path__[0])")/
make -C $PACKAGE_PATH install-qminos SOLVER_PATH=<path/to/solver_parent_directory>
```
4. Download the large files for building locally (11 gb of space are needed).

```console
make -C $PACKAGE_PATH build-data DATA_DIR=</desired/local_data/directory>
```

## About
For instructions on how to use the human ME Model, visit the [documentation](https://hmbaghdassarian.github.io/human_me/) 
