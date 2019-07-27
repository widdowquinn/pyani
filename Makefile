# Makefile
#
# This file is part of the pyani package distribution
# (https://github.com/widdowquinn/pyani)

# Set up all development dependencies in the current conda environment
setup_env:
	@conda install --file requirements-dev.txt --yes
	@conda install --file requirements.txt --yes
	@conda install --file requirements-thirdparty.txt --yes
	@pip install -r requirements-pip.txt
	@pre-commit install
	@pip install -U -e .

# Run all tests and display coverage report in a browser
test:
	@pytest --cov-report=html --cov=pyani -v tests/ && open htmlcov/index.html

# Build and display documentation
doc:
	@cd docs && make html && open _build/html/index.html

# Clean up test output and coverage
clean:
	@rm -rf htmlcov
	@rm -rf tests/test_output/*