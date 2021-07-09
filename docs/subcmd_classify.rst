.. _pyani-subcmd-classify:

==================
``pyani classify``
==================

The ``classify`` subcommand identifies cliques (k-complete) graphs of genomes from a larger set. This is helpful for circumscribing potentially meaningful groups of genomes that can not be described using traditional taxonomy.

.. code-block:: text

    usage: pyani.py classify [-h] [-l LOGFILE] [-v] [--debug]
    [--disable_tqdm] [--version] [--citation] [--dbpath DBPATH]
    [--cov_min COV_MIN] [--id_min ID_MIN] [--min_id MIN_ID]
    [--max_id MAX_ID] [--resolution RESOLUTION] [--show_all]
    outdir run_id

--------------------
Positional arguments
--------------------

``outdir``
   Path to a directory where comparison output files will be written.

``run_id``
  Unique database ID of the run to be plotted.

-----------------
Flagged arguments
-----------------

``--cov_min COV_MIN``
    Minimum percent coverage required to draw an edge between two nodes (genomes). (default: 0.5)

``--dbpath DBPATH``
    Path to the location of the local ``pyani`` database to be used. Default: ``.pyani/pyanidb``

``--disable_tqdm``
    Disable the ``tqdm`` progress bar while the download process runs. This is useful when testing to avoid aesthetic problems with test output. (default: False)

``--formats FORMATS``
    Graphics output format(s); more than one can be specified. Valid options are: (pdf/png/svg/jpg). (default: png)

``--id_min ID_MIN``
    Minimum percent identity required to draw an edge between two nodes (genomes). (default: 0.8)

``-l LOGFILE, --logfile LOGFILE``
    Provide the location ``LOGFILE`` to which a logfile of the download process will be written.

``--max_id MAX_ID``
    Maximum identity threshold to test. (default: None)

``--method {seaborn,mpl,plotly}``
    Graphics method to use for plotting. (default: seaborn)

``--min_id MIN_ID``
    Minimum identity threshold to test. (default: None)

``--resolution RESOLUTION``
    Number of identity thresholds to test. (default: 0.0001)

``--show_all``
    Report all intervals in log. (default: only intervals where all subgraphs are k-complete) (default: False)

``-v, --verbose``
    Provide verbose output to ``STDOUT``.
