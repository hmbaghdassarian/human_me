# SRC_CORE=human_me
#SRC_TEST=tests
#SRC_BENCHMARK=benchmarks

# PYTHON=python
# PIP=pip

help:
	@echo "Available Commands:"
	@echo " install-qminos         - Install the qMINOS solver."
	@echo " build-files            - Download the build files."
	# @echo " tests-coverage-html    - Run unit tests, code coverage and generate html."

install-qminos:
# 	$(PYTHON) -m pytest $(SRC_TEST)

build-files:
# 	$(PYTHON) -m pytest --cov=$(SRC_CORE) $(SRC_TEST)
# 	$(PYTHON) -m codecov

# test-coverage-html:
# 	$(PYTHON) -m pytest --cov=$(SRC_CORE) $(SRC_TEST) --cov-report=html
