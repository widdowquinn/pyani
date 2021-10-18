.. _pyani-subcmd-listdeps:

==================
``pyani listdeps``
==================

The ``listdeps`` subcommand writes an account of the local platform, the installed Python version, and installed dependencies (Python packages and third-party software tools) to ``STDOUT``.

.. code-block:: text

    usage: pyani listdeps [-h] [-l LOGFILE] [-v]

-----------------
Flagged arguments
-----------------

``--disable_tqdm``
    Disable the ``tqdm`` progress bar while the download process runs. Does nothing.

``-h, --help``
    Display usage information for ``pyani listdeps``.

``-l LOGFILE, --logfile LOGFILE``
    Provide the location ``LOGFILE`` to which a logfile of the subcommand output will be written.

``-v, --verbose``
    Provide verbose output to ``STDOUT``. This actually provides no additional information.
