# SRC_CORE=human_me
#SRC_TEST=tests
#SRC_BENCHMARK=benchmarks

PYTHON=python
SOLVER_PATH=./solver/
SOLVER_TYPE=qminos
# VENV=me_env
# PIP=pip

help:
	@echo "Available Commands:"
	@echo " install-qminos         - Install the qMINOS solver."
	@echo " build-files            - Download the build files."
	# @echo " tests-coverage-html    - Run unit tests, code coverage and generate html."

# make SOLVER_PATH=path/to/install_solver/ install-qminos
# solver path is the parent directory that should contain the qminos.tar.gz file
install-qminos:
	${PYTHON} install_solver.py ${SOLVER_PATH} ${SOLVER_TYPE}
	${PYTHON} ${SOLVER_PATH}/solvemepy/setup.py develop

# build-data:
	#YOURCODEHERE
