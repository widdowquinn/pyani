# Contributing to `pyani`

Contributions for bugfixes and addition of new features are welcomed

## Licensing

`pyani` is shared under the [MIT License](https://opensource.org/licenses/MIT) (see LICENSE file for details), and it will be presumed that any contributions you make will be licensed under this agreement.

## Git Usage

If you plan to make a pull request, please begin by forking this repository, and creating a new branch with a short, descriptive name, instead of using the `master` branch.

A workflow might look like this:

1. Identify something you think you can change or add
2. Fork this repository into your own account
3. Create a new branch with a short, descriptive name (for the thing you're fixing/adding), and work on this branch locally
4. When you're finished, check the changes/additions with `flake8` and test locally
5. Commit the branch, and submit a pull request
6. Continue the discussion in the [`Pull requests` section](https://github.com/widdowquinn/pyani/pulls) of this repository on GitHub.

## Coding Conventions

So far as is possible, we aim to follow the coding conventions as described in [PEP8](http://www.python.org/dev/peps/pep-0008/) and [PEP257](http://www.python.org/dev/peps/pep-0257/)

For now, docstrings are not in any controlled syntax, such as reStructuredText, but this may change.

### Code Style Checks

We use the `flake8` tool for style checks, and this can be installed as follows (with two useful plugins):

```bash
pip install flake8 flake8-docstrings flake8-blind-except
```

`flake8` can then be run directly on the codebase with

```bash
$ flake8 bin/
$ flake8 pyani/
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

We use [TravisCI](https://travis-ci.org/widdowquinn/pyani) to run tests.