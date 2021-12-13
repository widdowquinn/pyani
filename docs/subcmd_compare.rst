.. _pyani-subcmd-anim:

=================
``pyani compare``
=================

The ``compare`` subcommand compares ``pyani`` runs that may involve different methods or parameters, producing plots of the differences and summary files for each pair of runs. Runs for comparison are specified in the form of a set of runs to use as references, and a second set of runs as queries.

.. code-block:: text


    usage: pyani compare [-h] [-l LOGFILE] [-v] [--debug] [--disable_tqdm]
                         [--version] [--citation] -o OUTDIR --ref_ids RUN_ID
                         [RUN_ID ...] --run_ids RUN_ID [RUN_ID ...]
                         [--dbpath DBPATH] [--formats FORMAT [FORMAT ...]]
                         [--method [METHOD]] [--workers WORKERS]





.. _SQLite3: https://www.sqlite.org/index.html

-----------------
Flagged arguments
-----------------

``--dbpath DBPATH``
    Path to the location of the local ``pyani`` database to be used. Default: ``.pyani/pyanidb``

``--debug``
    Turn on debugging output.

``--disable_tqdm``
    Disable the ``tqdm`` progress bar while the download process runs. This is useful when testing to avoid aesthetic problems with test output.

``--formats FORMAT [FORMAT ...]``
    Graphics output format(s); more than one can be specified. Valid options are: (pdf/png/svg/jpg). (default: png)


``-h, --help``
    Display usage information for ``pyani anim``.

``-l LOGFILE, --logfile LOGFILE``
    Provide the location ``LOGFILE`` to which a logfile of the download process will be written.

``--method {seaborn,mpl,plotly}``
    Graphics method to use for plotting. (default: seaborn)

``-o OUTDIR, --outdir OUTDIR``
    Path to a directory where comparison output files will be written.

``--ref_ids RUN_ID [RUN_ID ...]``
    Space-separated list of run_ids to use as reference(s)
    for comparisons (default: None).

``--run_ids RUN_ID [RUN_ID ...]``
    Space-separated list of run_ids to compare to
    reference(s) (default: None).

``-v, --verbose``
    Provide verbose output to ``STDOUT``

``--workers WORKERS``
    Spawn WORKERS worker processes with the ``--scheduler multiprocessing`` option. Default: 0 (use all cores)
