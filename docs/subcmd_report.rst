.. _pyani-subcmd-report:

================
``pyani report``
================

The ``report`` subcommand reports the contents of a local `SQLite3`_ database, located at ``dbpath``. Output files will be written to the ``outdir`` directory.


usage: pyani report [-h] [-l LOGFILE] [-v] [--debug] [--disable_tqdm] [--version]
                    [--citation] -o OUTDIR [--dbpath DBPATH] [--runs] [--genomes]
                    [--runs_genomes] [--genomes_runs]
                    [--run_results RUN_ID [RUN_ID ...]]
                    [--run_matrices RUN_ID [RUN_ID ...]] [--no_matrix_labels]
                    [--formats FORMAT [FORMAT ...]]

.. _SQLite3: https://www.sqlite.org/index.html


-----------------
Flagged arguments
-----------------

``--dbpath DBPATH``
    Path to the location of the local ``pyani`` database to be used. Default: ``.pyani/pyanidb``

``--disable_tqdm``
    Disable the ``tqdm`` progress bar while the report process runs. This is useful when testing to avoid aesthetic problems with test output.

``--formats FORMAT [FORMAT ...]``
    Space-separated list of output formats (in addition to .tab);
    possible values: {html, excel,stdout}. (default: None)

``--genomes``
    Report table of genomes in database. (default: False)

``--genomes_runs``
    Report table of all runs in which each genome in the database participates. (default: False)

``--no_matrix_labels``
    Turn off row/column labels in output matrix files (default: False)

``-o OUTDIR, --outdir OUTDIR``
    Directory where output analysis reports will be written.

``--runs``
    Report table of analysis runs in database. (default: False)

``--runs_genomes``
    Report table of all genomes for each run in database. (default: False)

``--run_matrices RUN_ID [RUN_ID ...]``
    Report matrices of results for space-separated list of runs. (default: None)

``--run_results RUN_ID [RUN_ID ...]``
    Report table of results for space-separated list of runs. (default: None)

``-v, --verbose``
    Provide verbose output to ``STDOUT``.
