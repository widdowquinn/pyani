.. _pyani-subcmd-createdb:

==================
``pyani createdb``
==================

The ``createdb`` subcommand creates a new, empty database for ``pyani`` to use in subsequent analysis runs.

.. code-block:: text

    usage: pyani createdb [-h] [-l LOGFILE] [-v] [--disable_tqdm]
                         [--dbpath DBPATH] [-f]


-----------------
Flagged arguments
-----------------

``--disable_tqdm``
    Disable the ``tqdm`` progress bar while the download process runs. This is useful when testing to avoid aesthetic problems with test output.

``--dbpath DBPATH``
    Path to the location where the database will be created. Default: ``.pyani/pyanidb``

``-h, --help``
    Display usage information for ``pyani createdb``.

``-l LOGFILE, --logfile LOGFILE``
    Provide the location ``LOGFILE`` to which a logfile of the download process will be written.

``-v, --verbose``
    Provide verbose output to ``STDOUT``
