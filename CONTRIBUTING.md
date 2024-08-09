# Contributing to `pyani`

Contributions, including bugfixes and addition of new features, are welcomed!

The `pyani` package is currently maintained on GitHub under two main branches:

- `master` is the source code underpinning the most recent/current release of pyani. It should always be in sync with the latest release found at [https://github.com/widdowquinn/pyani/releases](https://github.com/widdowquinn/pyani/releases). The only time this code should not be in sync with the release is when there are modifications to documentation, or for a release candidate immediately preceding a release.
- `development` is the current bleeding-edge version of `pyani`. It should (almost) always be in a working and usable condition, but may not be complete and/or some features may be missing or still under development.

Additional branches may also be found on GitHub, with bug fixes and or new features in varying stages of completeness.

## Table of Contents

<!-- TOC -->

- [Contributing to `pyani`](#contributing-to-pyani)
  - [Table of Contents](#table-of-contents)
  - [Licensing](#licensing)
  - [How to Contribute](#how-to-contribute)
    - [Getting the source code](#getting-the-source-code)
    - [Setting up a development environment](#setting-up-a-development-environment)
      - [Cleaning up development environment](#cleaning-up-development-environment)
    - [Editing source code](#editing-source-code)
    - [Editing documentation](#editing-documentation)
  - [Pull Requests](#pull-requests)
    - [Naming your branch](#naming-your-branch)
    - [Reviewing changes](#reviewing-changes)
    - [Coding Conventions](#coding-conventions)
    - [Code Style Checks](#code-style-checks)
    - [Commit Message Convention](#commit-message-convention)
    - [Pre-commit Checks](#pre-commit-checks)
    - [Testing](#testing)
    - [Local Testing](#local-testing)
    - [Continuous Integration](#continuous-integration)
  - [Versioning](#versioning)

<!-- /TOC -->

## Licensing

`pyani` is shared under the [MIT License](https://opensource.org/licenses/MIT) (see LICENSE file for details), and any contributions you make will be licensed under this agreement.

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

then set up `conda` channels for the required dependencies:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

and then use the command

```bash
make setup_env
```

This will install all dependencies for running and developing `pyani`, as well as pre-commit hooks. Once installation is complete, run the test suite to check installation, availability of dependencies, and code coverage:

```bash
make test
```

If you want to be able to edit source files and have those changes take immediate effect when calling `pyani` (useful for testing), clone the GitHub repository with:

```bash
git clone https://github.com/widdowquinn/pyani.git
```

then inside the new `pyani` directory run:

```bash
pip install -e .
```

This is the [`pip install --editable`](https://pip.pypa.io/en/stable/cli/pip_install/#install-editable) command, which links the installed package to the specified location (here `.`, i.e. the current directory) rather than the usual package location (`site-packages`). When using this option, edits to the source code are immediately available in the installed package. This allows you to test changes to the source code as you make them, without the need for an additional uninstall/install step.

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

If you plan to make a pull request, please begin by forking this repository, and creating a new branch with a short, descriptive name, instead of using the `development`, `version_0_2`, `version_0_3`, `main`, or `master` branch.

A workflow might look like this:

1. Identify something you think you can change or add
2. Fork this repository into your own account (but [please add write access for repository maintainers](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/working-with-forks/allowing-changes-to-a-pull-request-branch-created-from-a-fork))
3. Obtain source code and set up a development environment as described above
4. Create a new branch with a short, descriptive name (for the thing you're fixing/adding), and work on this branch locally
5. When you're finished, check the changes/additions with `flake8`/`black` and test locally with `pytest -v`.
6. Commit the branch, and submit a pull request
7. Continue the discussion in the [`Pull requests` section](https://github.com/widdowquinn/pyani/pulls) of this repository on GitHub.

**All tests must pass before a pull request will be merged at `GitHub`.**

Please keep each pull request as *atomic* as possible - fixing or adding a single conceptually-complete issue/functionality. If you would like to add several features or fix several issues, please use a separate pull request for each one, where possible.

### Naming your branch

- If you are proposing a fix to an issue, please give your branch a name that reflects this. For example, if you were proposing to fix issue #150, please call your branch `issue_150` (or some similar variation).
- If you are handling a pull request from a fork, with a view to merging, please indicate the PR number somewhere in the local branch name when you pull the forked version into the repository. For example, if you were handling pull request #250 which refers to a fork from `kenny_loggins/pyani`, please call your local branch `pr_250` (or some similar variation).

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

### Coding Conventions

So far as is possible, we aim to follow the coding conventions as described in [PEP8](http://www.python.org/dev/peps/pep-0008/) and [PEP257](http://www.python.org/dev/peps/pep-0257/), but we have adopted `black` code styling, which does vary from the PEPs in places.

We are moving to writing all docstrings as reStructuredText, for usage with Sphinx.

### Code Style Checks

We use the `flake8`, `doc8`, and `pylint` tools for style checks. These will be installed if you have used `make setup_env` to set up your development environment.

### Commit Message Convention

`git` commit messages are an important way to manage a readable revision history. We use the following conventions:

- **If a commit fixes an issue, state this in the commit message**
  - `GitHub` Projects picks up on terms like `fixes #123` and `closes #987`. Using these phrases makes project management much easier.

- **Every commit gets a short description**
  - The short description should be in *imperative form* and *around 50 characters or less*
  - The short description can, but need not, include the name of the file that is modified
  - There should be a short account of what the change does
  - The short discription *should not only contain the name of the file that is modified*

For example, if the commit updates some documentation, the following are good short descriptions:

- `update citations.rst to add new references`
- `update docs to add new references; fixes #567`
- `add new references to citations.rst`

The following are bad short descriptions

- `update citations.rst` (does not say what was done)
- `there were some new references so I added them` (not in imperative form)
- `citations.rst` (does not say what was done)
- `part of some doc updates` (does not say what was done)

- **Some commits get long/extended descriptions**
  - Long descriptions should be included where there is more information to share than can be fit in the short description
  - Long descriptions are free-form
  - Long descriptions should explain the *why* of the change. They may also explain the *what* at a high level, but should not be excessively detailed.
  - They are "long descriptions", but they should be concise, precise and clear
  - Paragraphs are fine
  - Bullet points are fine

A good long description could be

> This fix improves efficiency of the veeblefetzer. The main change is replacing a
> nested loop with asyncio calls to a new function `fetzveebles()`. This commit
> makes affects `veebles.py`, and new tests are added in `test_veeblefetzer.py`.
>
> fixes #246

A bad long description might be

> So, because it was taking too long to fetz all the veebles I looked at the
> structure. It took a while to identify the main problem, which was that there
> were a bunch of nested loops. I timed these and it turned out that one of
> them really hit performance. So I looked into a load of options (like
> multiprocessing, threading and so on) but, inspired by a Stack Overflow
> post (URL HERE:) I decided to try asyncio. After some trial and error I
> managed to get something working with that. I also wrote some tests. I
> also think this might fix one of the issues.

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

We currently use [GitHub Actions](https://github.com/widdowquinn/pyani/actions/) to run tests. The configuration file can be found in the repository under `.github/workflows/`.

Currently, `pyani` is tested under Python 3.6 and 3.7, and coverage is reported at [`CodeCov`](https://codecov.io/gh/widdowquinn/pyani/branch/development).

## Versioning

We aim to comply with [`PEP440`](https://www.python.org/dev/peps/pep-0440/) and [semantic versioning](https://semver.org/).
