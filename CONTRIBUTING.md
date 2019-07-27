# Contributing to `pyani`

Contributions, including bugfixes and addition of new features, are welcomed!

The `pyani` package is currently maintained on GitHub under two main branches:

- `master` is the source code underpinning the most recent/current release of pyani. It should always be in sync with the latest release found at [https://github.com/widdowquinn/pyani/releases](https://github.com/widdowquinn/pyani/releases). The only time this code should not be in sync with the release is when there are modifications to documentation, or for a release candidate immediately preceding a release.
- `development` is the current bleeding-edge version of `pyani`. It should (almost) always be in a working and usable condition, but may not be complete and/or some features may be missing or still under development.

Additional branches may also be found on GitHub, with bug fixes and or new features in varying stages of completeness.

## Table of Contents

<!-- TOC -->

- [Table of Contents](#table-of-contents)
- [Licensing](#licensing)
- [How to Contribute](#how-to-contribute)
    - [Getting the source code](#getting-the-source-code)
    - [Setting up a development environment](#setting-up-a-development-environment)
    - [Editing source code](#editing-source-code)
    - [Editing documentation](#editing-documentation)
- [Pull Requests](#pull-requests)
    - [Reviewing changes](#reviewing-changes)
    - [Coding Conventions](#coding-conventions)
    - [Code Style Checks](#code-style-checks)
    - [Pre-commit Checks](#pre-commit-checks)
    - [Testing](#testing)
    - [Local Testing](#local-testing)
    - [Continuous Integration](#continuous-integration)

<!-- /TOC -->

## Licensing

`pyani` is shared under the [MIT License](https://opensource.org/licenses/MIT) (see LICENSE file for details), and it any contributions you make will be licensed under this agreement.

## How to Contribute

### Getting the source code

Please fork the repository to your own GitHub account (you will need to create a GitHub account if you do not already have one), and clone the repository from your fork. Using SSH (recommended):

```bash
git clone git@github.com:<YOUR USERNAME>/pyani.git
cd pyani
```

or HTTPS:

```bash
git clone https://github.com/<YOUR USERNAME>/pyani.git
cd curveball
```

### Setting up a development environment

`pyani` development makes use of the following tools and packages:

- `Anaconda`/`Miniconda` for virtual environment management
- `git` for version control
- `pytest` for testing
- `pre-commit` to manage pre-commit hooks
- `doc8`, `flake8` and `pylint` for code linting
- `black` for code formatting
- `bandit` to identify security issues
- `codecov` to measure code coverage

To set up a local development environment with these tools configured for `pyani` development, first create and activate a new `conda` environment:

```bash
conda create --yes --name pyani_dev python=3.7 && conda activate pyani_dev
```

and then use the command

```bash
make setup_env
```

This will install all dependencies for running and developing `pyani`, as well as pre-commit hooks. Once installation is complete, run the test suite to check installation, availability of dependencies, and code coverage:

```bash
make test
```

#### Cleaning up development environment

You can remove the `conda` development environment with the following commands:

```bash
conda deactivate
conda remove -n pyani_dev --all
```

### Editing source code

The root directory of the repository has subdirectories specific for testing, packaging, and deploying `pyani`:

- `pyani` source code is under the `pyani` subdirectory
    - code for the CLI scripts is in the `pyani/scripts` subdirectory
- tests (written for the `pytest` framework) and test/target data are in the `tests` subdirectory
- documentation is found under the `docs` subdirectory

### Editing documentation

Package documentation for `pyani` is found in the `docs` subdirectory. Both the package documentation and docstrings are written in [`reStructured Text`](http://sphinx-doc.org/rest.html) (i.e. `RST`), and to build the documentation you will need to install [`Sphinx`](https://www.sphinx-doc.org/en/master/).

Although you can build the documentation in multiple formats, the default is the HTML website. Tho build the package docs and view them in your browser you can use the command:

```bash
make docs
```

from the repository root directory

## Pull Requests

If you plan to make a pull request, please begin by forking this repository, and creating a new branch with a short, descriptive name, instead of using the `development` branch.

A workflow might look like this:

1. Identify something you think you can change or add
2. Fork this repository into your own account
3. Obtain source code and set up a development environment as described above
4. Create a new branch with a short, descriptive name (for the thing you're fixing/adding), and work on this branch locally
5. When you're finished, check the changes/additions with `flake8`/`black` and test locally with `pytest -v`.
6. Commit the branch, and submit a pull request
7. Continue the discussion in the [`Pull requests` section](https://github.com/widdowquinn/pyani/pulls) of this repository on GitHub.

**All tests must pass before a pull request will be merged at `GitHub`.**

### Reviewing changes

When you have finished editing source or documentation, please check that everything you modified is committed:

```bash
git status
```

Please also review the differences that have been introduced relative to the `origin/master` or `origin/development` branch (as appropriate):

```bash
git diff origin/master
```

Please check your `git` log and consider rebasing or squashing any commits that have not been pushed:

```bash
git log
```

### Coding Conventions

So far as is possible, we aim to follow the coding conventions as described in [PEP8](http://www.python.org/dev/peps/pep-0008/) and [PEP257](http://www.python.org/dev/peps/pep-0257/), but we have adopted `black` code styling, which does vary from the PEPs in places.

We are moving to writing all docstrings as reStructuredText, for usage with Sphinx.

### Code Style Checks

We use the `flake8`, `doc8`, and `pylint` tools for style checks. These will be installed if you have used `make setup_env` to set up your development environment.

### Pre-commit Checks

The `flake8` and `black` styles are enforced as pre-commit hooks using the `pre-commit` package (included in `requirements-dev.txt` and `requirements-pip.txt`).

The `black` and `flake8` hooks are defined in `.pre-commit-config.yaml`. Custom settings for `flake8` are held in `.flake8`.

To enable pre-commit checks in the codebase on your local machine (once `pre-commit` has been installed), execute the following command in the root directory of this repository:

```bash
pre-commit install
```

This is done automatically if you use `make setup_env` to set up your development environment.

### Testing

New features or functions will not be accepted without tests. Bug fixes should ideally include an additional test to establish that the bug has been squashed. We expect that you will have run tests locally before a pull request is made.

Tests will also be run as part of continuous integration, and changes will not be accepted until continuous integration tests have been passed.

Tests are located in the `tests` subdirectory of this repository.

### Local Testing

We currently use `pytest` for testing, and `codecov` to establish how much of the codebase is covered by tests. These can be installed as follows:

```bash
conda install pytest codecov
```

`pytest` can then be run on the codebase with

```bash
pytest --cov-report=html --cov=pyani -v tests/
```

### Continuous Integration

When you submit a pull request on GitHub, automated tests will be run, and results reported on the pull request. **All tests must pass before a pull request will be merged.**

We currently use [TravisCI](https://travis-ci.org/widdowquinn/pyani) to run tests. The configuration file can be foudn in the repository as `.travis.yml`
