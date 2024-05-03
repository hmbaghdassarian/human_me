# SRC_CORE=human_me
#SRC_TEST=tests
#SRC_BENCHMARK=benchmarks

PYTHON=python

SOLVER_PATH=./solver/
SOLVER_TYPE=qminos

DATA_DIR=./data/
PREBUILD=0 # 0 or 1

help:
	@echo "Available Commands:"
	@echo " install-qminos         - Install the qMINOS solver."
	@echo " build-files            - Download the build files and Recon2.2 input."

# make SOLVER_PATH=path/to/install_solver/ install-qminos
# solver path is the parent directory that should contain the qminos.tar.gz file
install-qminos:
	${PYTHON} install_solver.py ${SOLVER_PATH} ${SOLVER_TYPE}
#	${PYTHON} ${SOLVER_PATH}solvemepy/setup.py develop

# make DATA_DIR=/data3/hratch/test_data PREBUILD=0 build-data
build-data:
	${PYTHON} ./data/_download_data.py ${DATA_DIR} ${PREBUILD}

# TODO: the two make commands can be combined (parent dir with solver directory + data directories) if QMINOS can be downloaded