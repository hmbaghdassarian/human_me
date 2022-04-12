# human_me
human_me is a python package to generate ME-Models from input context-extracted Recon2.2 M-Models

## Installation
For additional details and alternate install options, see the documentation's [installation instructions](https://hmbaghdassarian.github.io/human_me/install/). <br>
Requirements: gfortran (>=4.6)


1. Create a [python virtual environment](https://packaging.python.org/guides/) and activate this environment
2. Install human_me:
```console
pip install human_me
```
3. Set up the qMINOS solver. You will need the qminos file, which can be obtained for academic use from Prof. Michael Saunders at Stanford University.

&emsp;i) download the qminos file into a specified solver directory, which we refer to here as "solver_parent_directory".

&emsp;ii) the solver can be installed using the human_me Makefile as follows:

```console
make -C <path/to/human_me/> install-qminos SOLVER_PATH=<path/to/solver/solver_parent_directory>
```
4. Download the large files for building locally (11 gb of space are needed).

```console
make -C <path/to/human_me/> build-data DATA_DIR=</desired/local_data/directory>
```

## About
For instructions on how to use the human ME Model, visit the [documentation](https://hmbaghdassarian.github.io/human_me/) 
