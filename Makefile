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
docs: clean_docs
	@cd docs && make html && open _build/html/index.html

# Clean up test, walkthrough, and coverage output
clean: clean_walkthrough clean_tests clean_coverage

clean_coverage:
	@rm -rf htmlcov

clean_docs:
	@rm -rf docs/_build/html

clean_tests:
	@rm -rf tests/test_output/*

clean_walkthrough:
	@rm -rf C_blochmannia
	@rm -rf C_blochmannia_ANIm
	@rm -rf .pyani/pyanidb

# Run walkthrough
walkthrough: clean_walkthrough
	pyani download --email my.email@my.domain -t 203804 -o C_blochmannia
	pyani createdb -f
	pyani anim -i C_blochmannia -o C_blochmannia_ANIm \
        --name "C. blochmannia run 1" \
        --labels C_blochmannia/labels.txt --classes C_blochmannia/classes.txt
	pyani report --runs -o C_blochmannia_ANIm/ --formats html excel stdout
	pyani report --run_results 1 --formats html excel stdout -o C_blochmannia_ANIm/
	pyani report --run_matrices 1 --formats html excel stdout -o C_blochmannia_ANIm/
	pyani plot --formats png pdf --method seaborn -o C_blochmannia_ANIm --run_ids 1
	# pyani anib C_blochmannia C_blochmannia_ANIb \
    #     --name "C. blochmannia run 2" \
    #     --labels C_blochmannia/labels.txt --classes C_blochmannia/classes.txt
	# pyani report --runs C_blochmannia_ANIb/ --formats html,excel,stdout
	# pyani report --run_results 2 --formats html,excel,stdout C_blochmannia_ANIb/
	# pyani report --run_matrices 2 --formats html,excel,stdout C_blochmannia_ANIb/
	# pyani plot --formats png,pdf --method seaborn C_blochmannia_ANIb 2	

uml:
	pyreverse -o pdf -p pyani pyani
