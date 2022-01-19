.. _pyani-subcmd-plot:

==============
``pyani plot``
==============

The ``plot`` subcommand will plot the results of an ANI analysis, specified by ``run_id``, contained within a local `SQLite3`_ database, located at ``dbpath``. Plot files will be written to the ``outdir`` directory.

.. code-block:: text

    usage: pyani plot [-h] [-l LOGFILE] [-v] [--debug] [--disable_tqdm] [--version]
                  [--citation] -o OUTDIR --run_ids RUN_ID [RUN_ID ...]
                  [--dbpath DBPATH] [--formats FORMAT [FORMAT ...]]
                  [--method METHOD] [--workers WORKERS]

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

``--formats FORMAT [FORMAT ...]``
    Graphics output format(s); more than one can be specified. Valid options are: (pdf/png/svg/jpg). (default: png)

``-l LOGFILE, --logfile LOGFILE``
    Provide the location ``LOGFILE`` to which a logfile of the plotting process will be written.

``--method METHOD``
    Graphics method to use for plotting; options (seaborn, mpl, plotly). (default: seaborn)

``-o OUTDIR, --o outdir OUTDIR``
    Path to a directory where comparison plot files will be written.

``--run_ids RUN_ID [RUN_ID ...]``
    Unique database ID of the runs to be plotted.
