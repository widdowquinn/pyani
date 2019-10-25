.. _pyani-subcmd-anim:

==============
``pyani anim``
==============

The ``anim`` subcommand will carry out ANIm analysis using genome files contained in the ``indir`` directory, writing result files to the ``outdir`` directory, and recording data about each comparison and run in a local `SQLite3`_ database.

.. code-block:: text

    usage: pyani.py anim [-h] [-l LOGFILE] [-v] [--disable_tqdm]
                     [--scheduler {multiprocessing,SGE}] [--workers WORKERS]
                     [--SGEgroupsize SGEGROUPSIZE] [--SGEargs SGEARGS]
                     [--jobprefix JOBPREFIX] [--name NAME] [--classes CLASSES]
                     [--labels LABELS] [--recovery] [--dbpath DBPATH]
                     [--nucmer_exe NUCMER_EXE] [--filter_exe FILTER_EXE]
                     [--maxmatch] [--nofilter]
                     indir outdir



.. _SQLite3: https://www.sqlite.org/index.html

--------------------
Positional arguments
--------------------

``indir``
    Path to the directory containing indexed genome files to be used for the analysis.

``outdir``
    Path to a directory where comparison output files will be written.

-----------------
Flagged arguments
-----------------

``--classes CLASSFNAME``
    Use the set of classes (one per genome sequence file) found in the file ``CLASSFNAME`` in ``indir``. Default: ``classes.txt``

``--dbpath DBPATH``
    Path to the location of the local ``pyani`` database to be used. Default: ``.pyani/pyanidb``

``--disable_tqdm``
    Disable the ``tqdm`` progress bar while the download process runs. This is useful when testing to avoid aesthetic problems with test output.

``--filter_exe FILTER_EXE``
    Path to the ``MUMmer`` ``delta-filter`` executable. Default: ``delta-filter``

``-h, --help``
    Display usage information for ``pyani index``.

``--jobprefix JOBPREFIX``
    Use the string ``JOBPREFIX`` as a prefix for SGE job submission names. Default: ``PYANI``

``--labels LABELFNAME``
    Use the set of labels (one per genome sequence file) found in the file ``LABELFNAME`` in ``indir``. Default: ``labels.txt``

``--name NAME``
    Use the string ``NAME`` to identify this ANIm run in the ``pyani`` database.

``--nucmer_exe NUCMER_EXE``
    Path to the ``MUMmer`` ``nucmer`` executable. Default: ``nucmer``

``-l LOGFILE, --logfile LOGFILE``
    Provide the location ``LOGFILE`` to which a logfile of the download process will be written.

``--maxmatch``
    Use the ``MUMmer`` ``--maxmatch`` option to include all ``nucmer`` matches.

``--nofilter``
    Do not use ``delta-filter`` to restrict ``nucmer`` output to 1:1 matches.

``--recovery``
    Use existing ``NUCmer`` comparison output if available, e.g. if recovering from a failed job submission. Using this option will not generate a new comparison if the old output files exist.

``--scheduler {multiprocessing, SGE}``
    Specify the job scheduler to be used when parallelising genome comparisons: one of ``multiprocessing`` (use many cores on the current machine)  or ``SGE`` (use an SGE or OGE job scheduler). Default: ``multiprocessing``.

``--SGEargs SGEARGS``
    Pass additional arguments ``SGEARGS`` to ``qsub`` when running the SGE-distributed jobs.

``--SGEgroupsize SGEGROUPSIZE``
    Create SGE arrays containing SGEGROUPSIZE comparison jobs. Default: 10000

``-v, --verbose``
    Provide verbose output to ``STDOUT``

``--workers WORKERS``
    Spawn WORKERS worker processes with the ``--scheduler multiprocessing`` option. Default: 0 (use all cores)
