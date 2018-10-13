# Contributing to `pyani`

Contributions for bugfixes and addition of new features are welcomed!

The `pyani` package is presented at GitHub under two main branches:

- `master` is the source code underpinning the most recent/current release of pyani. It will (almost) always be in sync with the latest release found at [https://github.com/widdowquinn/pyani/releases](https://github.com/widdowquinn/pyani/releases). The only time this code should not be in sync with the release is when there are modifications to documentation, or immediately preceding a release.
- `development` is the current bleeding-edge version of `pyani`. It should (almost) always be in a working and usable condition, but may not be complete and/or some features may be missing or still under development.

Additional branches may also be available, with bug fixes and or new features in varying stages of completeness.

## Table of Contents

<!-- TOC -->

- [Table of Contents](#table-of-contents)
- [Licensing](#licensing)
- [Git Usage](#git-usage)
- [Coding Conventions](#coding-conventions)
    - [Code Style Checks](#code-style-checks)
    - [Pre-commit Checks](#pre-commit-checks)
- [Testing](#testing)
    - [Local Testing](#local-testing)
    - [Continuous Integration](#continuous-integration)

<!-- /TOC -->

## Licensing

`pyani` is shared under the [MIT License](https://opensource.org/licenses/MIT) (see LICENSE file for details), and it will be presumed that any contributions you make will be licensed under this agreement.

## Git Usage

If you plan to make a pull request, please begin by forking this repository, and creating a new branch with a short, descriptive name, instead of using the `development` branch.

A workflow might look like this:

1. Identify something you think you can change or add
2. Fork this repository into your own account
3. Create a new branch with a short, descriptive name (for the thing you're fixing/adding), and work on this branch locally
4. When you're finished, check the changes/additions with `flake8`/`black` and test locally
5. Commit the branch, and submit a pull request
6. Continue the discussion in the [`Pull requests` section](https://github.com/widdowquinn/pyani/pulls) of this repository on GitHub.

**All tests must pass before a pull request will be merged at `GitHub`.**

## Coding Conventions

So far as is possible, we aim to follow the coding conventions as described in [PEP8](http://www.python.org/dev/peps/pep-0008/) and [PEP257](http://www.python.org/dev/peps/pep-0257/), but we have adopted `black` code styling, which does vary from the PEPs in places.

For now, docstrings are not in any controlled syntax, such as reStructuredText, but this may change.

### Code Style Checks

We use the `flake8` tool for style checks, and this can be installed as follows (with two useful plugins):

```bash
pip install flake8 flake8-docstrings flake8-blind-except
```

`flake8` can then be run directly on the codebase with

```bash
flake8 bin/
flake8 pyani/
```

We use the `black` tool for code style checking, which can be installed with:

```bash
pip install black
```

### Pre-commit Checks

The `flake8` and `black` styles can be enforced as pre-commit hooks using the `pre-commit` package (included in `requirements.txt`).

The `black` and `flake8` hooks are defined in `.pre-commit-config.yaml`. Custom settings for `flake8` are held in `.flake8`.

To enable pre-commit checks in the codebase on your local machine (once `pre-commit` has been installed), execute the following command in the root directory of this repository:

```bash
pre-commit install
```

## Testing

New features or functions will not be accepted without tests. Bug fixes should ideally include an additional test to establish that the bug has been squashed. We expect that you will have run tests locally before any pull request is made. Tests will also be run as part of continuous integration, and changes will not be accepted until continuous integration tests have been passed.

Tests are located in the `tests` subdirectory of this repository.

### Local Testing

We currently use `nose` for testing, and `coverage` to establish how much of the codebase is covered by tests. These can be installed as follows:

```bash
pip install nose coverage
```

`nose` can then be run on the codebase with

```bash
nosetests --cover-erase --cover-package=pyani --cover-html --with-coverage -v
```

### Continuous Integration

When you submit a pull request on GitHub, automated tests will be run, and results reported on the pull request. **All tests must pass before a pull request will be merged.**

We use [TravisCI](https://travis-ci.org/widdowquinn/pyani) to run tests. The configuration file is kept in the repository as `.travis.yml`