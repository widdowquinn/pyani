.. _pyani-subcmd-plot:

==============
``pyani plot``
==============

The ``plot`` subcommand will plot the results of an ANI analysis, specified by ``run_id``, contained within a local `SQLite3`_ database, located at ``dbpath``. Plot files will be written to the ``outdir`` directory.

.. code-block:: text

    usage: pyani.py plot [-h] [-l LOGFILE] [-v] [--debug]
    [--disable_tqdm] [--version] [--citation] [--dbpath DBPATH]
    [--formats FORMATS] [--method {seaborn,mpl,plotly}]
    outdir run_id

.. _SQLite3: https://www.sqlite.org/index.html

--------------------
Positional arguments
--------------------

``outdir``
   Path to a directory where comparison plot files will be written.

``run_id``
  Unique database ID of the run to be plotted.

-----------------
Flagged arguments
-----------------

``--dbpath DBPATH``
    Path to the location of the local ``pyani`` database to be used. Default: ``.pyani/pyanidb``

``--disable_tqdm``
    Disable the ``tqdm`` progress bar while the plotting process runs. This is useful when testing to avoid aesthetic problems with test output.

``--formats FORMATS``
    Graphics output format(s); more than one can be specified. Valid options are: (pdf/png/svg/jpg). (default: png)

``-l LOGFILE, --logfile LOGFILE``
    Provide the location ``LOGFILE`` to which a logfile of the plotting process will be written.

``--method {seaborn,mpl,plotly}``
    Graphics method to use for plotting. (default: seaborn)
