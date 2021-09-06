.. _pyani-contributing:

=========================
Contributing to ``pyani``
=========================

-------------------------
Reporting bugs and errors
-------------------------

If you find a bug, or an error in the code or documentation, please report this by raising an issue at the `GitHub issues page`_ for ``pyani``:

- `GitHub issues page`_

----------------------------------
Contributing code or documentation
----------------------------------

We gratefully accept code and other contributions. A list of contributors can be found via the `Github contributors`_ link.

You are welcome to help develop ``pyani``, fix a bug, improve documentation, or contribute in any other way. To make everyone's lives easier in this process, we ask that you please follow the guidelines for developers below:


:::::::::::::::::::::::::::::::::
Pre-commit checks and style guide
:::::::::::::::::::::::::::::::::

So far as is possible, we aim to follow the coding conventions as described in `PEP8`_ and `PEP257`_, but we have adopted ``black`` code styling, which does vary from the PEPs in places.

We use the ``flake8`` tool for style checks, and this can be installed as follows (with two useful plugins):

.. code-block:: bash

    pip install flake8 flake8-docstrings flake8-blind-except

``flake8`` can then be run directly on the codebase with

.. code-block:: bash

    flake8 bin/
    flake8 pyani/

We use the ``black`` tool for code style checking, which can be installed with:

.. code-block:: bash

    pip install black

The ``flake8` and ``black`` styles can be enforced as pre-commit hooks using the `pre-commit`_ package (included in ``requirements.txt``).

The ``black`` and ``flake8`` hooks are defined in ``.pre-commit-config.yaml``; custom settings for ``flake8`` are held in ``.flake8`` (all files are under version control).

To enable pre-commit checks in the codebase on your local machine (once ``pre-commit`` has been installed), execute the following command in the root directory of this repository:

.. code-block:: bash

    pre-commit install


:::::::::::::::::::::::::::::::::::::
Checking changes to the documentation
:::::::::::::::::::::::::::::::::::::

Much of the repository documentation is written in Markdown files, but the main documentation (which you are reading) is prepared for `ReadTheDocs`_, which uses reStructuredText and ``Sphinx``. The ``Sphinx`` configuration is described in ``docs/conf.py`` (under version control).

So long as ``Sphinx`` is installed on your machine, you can check your documentation changes locally, building inplace by changing to the ``docs/`` directory and issuing:

.. code-block:: bash

    make html

This will place a compiled version of the documentation under ``_build/html``, which you can inspect before committing to the repository.

.. TIP::
    To build the documentation in ``ReadTheDocs`` style, you will need to install the corresponding theme with ``pip install sphinx_rtd_theme`` or ``conda install sphinx_rtd_theme``

For now, docstrings in the source code are not required to be in any controlled syntax, such as reStructuredText, but this may change.

::::::::::::::::::::::::::::::::
Making changes and pull requests
::::::::::::::::::::::::::::::::

1. Fork the ``pyani`` `repository`_ under your account at `GitHub`_.
2. Clone your fork to your development machine.
2a. To be able to run `pyani` and have changes you make take effect (useful for testing), you can run `pip install -e .` inside the local cloned repository. 
3. Create a new branch in your forked repository with an informative name like ``fix_issue_107``, using ``git`` (e.g. with the command ``git checkout -b fix_issue_107``).
4. Make the changes you need and commit them to your local branch.
5. Run the repository tests (see the :ref:`pyani-testing` documentation for more details).
6. If the tests all pass, push the changes to your fork, and submit a pull request against the original repository.
7. Indicate one of the ``pyani`` developers as an assignee to review your pull request when you submit your pull request.

The assigned developer will then review your pull request, and merge it or continue the conversation, as appropriate.



---------------------------
Suggestions for improvement
---------------------------

If you would like to make a suggestion for how we could improve ``pyani``, we welcome contributions. If you have a specific problem, or a concrete suggestion, you can submit these at the `GitHub issues page`_. If you would like to discuss an idea with the maintainers and the ``pyani`` community, this can be done at the `Github discussions page`_.

.. _GitHub: https://github.com
.. _Github contributors: https://github.com/widdowquinn/pyani/blob/master/CONTRIBUTORS.md
.. _Github issues page: https://github.com/widdowquinn/pyani/issues
.. _Github discussions page: https//github.com/widdowquinn/pyani/discussions
.. _PEP8: http://www.python.org/dev/peps/pep-0008/
.. _PEP257: http://www.python.org/dev/peps/pep-0257/
.. _pre-commit: https://github.com/pre-commit/pre-commit
.. _ReadTheDocs: https://docs.readthedocs.io/en/latest/#
.. _repository: https://github.com/widdowquinn/pyani
