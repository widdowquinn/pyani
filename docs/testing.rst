.. _pyani-testing:

=======
Testing
=======

We are currently writing tests formatted for the `nosetests`_ package, for testing ``pyani``.

.. WARNING::
    We are aware that ``nosetests`` is in maintenance mode and, while we have no timetable for the move, our plan is to change the test framework to use `pytest`_ at some future date

------------------------
Test directory structure
------------------------

All tests, supporting input and target data, are to be found in the ``tests/`` subdirectory of the ``pyani`` repository.

- Input data for tests is located under the ``tests/test_input`` directory, in subdirectories named for the general operation that is being tested.
- Target data (known correct values) for tests is located under the ``tests/test_targets`` directory, in subdirectories named for the general operation being tested.
- Test output is written to the ``tests/test_output`` directory, in subdirectories named for the general operation being tested.

.. TIP::
    If you wish to write new tests for ``pyani``, we ask that your test data and operations conform to this structure


-------------
Running tests
-------------

To run tests with ``nosetests``, change directory to the root of the ``pyani`` repository, and invoke a ``nosetests`` command.

^^^^^^^^^^^^^^^^^^^^^
Run all tests locally
^^^^^^^^^^^^^^^^^^^^^

To run all tests locally on your machine, issue the following command from the repository root:

.. code-block:: bash

    nosetests -v

This will cause ``nose`` to run all tests under the ``tests/`` subdirectory.


^^^^^^^^^^^^^^^^^^^^
Run individual tests
^^^^^^^^^^^^^^^^^^^^

Tests are grouped in files with filenames that match ``test_*.py``. We aim to write tests as classes that subclass ``unittest.TestCase``, as described in the `nosetests`_ documentation. An example of this style can be found in the ``tests/test_anim.py`` test file.

This style allows us to run tests at several levels of granularity, specifying all tests (see above), all tests within a module (e.g. ``test_anim.py``), or all tests within a class in that test file.

For example, to run all ANIm-related tests, we can issue:

.. code-block:: bash

    nosetests -v tests/test_anim.py


To run all tests of ``nucmer`` command line generation, we can specify a single class within that file using the command:

.. code-block:: bash

    nosetests -v tests/test_anim.py:TestNUCmerCmdline

And to test only "multiple command generation" we can issue the following:

.. code-block:: bash

    nosetests -v tests/test_anim.py:TestNUCmerCmdline.test_multi_cmd_generation



.. _nosetests: https://nose.readthedocs.io/en/latest/
.. _pytest: https://docs.pytest.org/en/latest/