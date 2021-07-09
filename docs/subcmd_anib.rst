.. _pyani-subcmd-anib:

==============
``pyani anib``
==============

The ``anib`` subcommand will carry out ANIb analysis using genome files contained in the ``indir`` directory, writing result files to the ``outdir`` directory, and recording data about each comparison and run in a local `SQLite3`_ database.

.. code-block:: text

    usage: pyani.py anib [-h] [-l LOGFILE] [-v] [--disable_tqdm]
                     [--scheduler {multiprocessing,SGE}] [--workers WORKERS]
                     [--SGEgroupsize SGEGROUPSIZE] [--SGEargs SGEARGS]
                     [--jobprefix JOBPREFIX] [--name NAME] [--classes CLASSES]
                     [--labels LABELS] [--recovery] [--dbpath DBPATH]
                     [--blastn_exe BLASTN_EXE] [--format_exe FORMAT_EXE]
                     [--fragsize] -i INDIR -o OUTDIR

.. _SQLite3: https://www.sqlite.org/index.html

-----------------
Flagged arguments
-----------------

``-i, --input INDIR``
    Path to the directory containing indexed genome files to be used for the analysis.

``-o, --outdir OUTDIR``
    Path to a directory where comparison output files will be written.

``--classes CLASSFNAME``
    Use the set of classes (one per genome sequence file) found in the file ``CLASSFNAME`` in ``INDIR``. Default: ``classes.txt``

``--dbpath DBPATH``
    Path to the location of the local ``pyani`` database to be used. Default: ``.pyani/pyanidb``

``--disable_tqdm``
    Disable the ``tqdm`` progress bar while the ANIb process runs. This is useful when testing to avoid aesthetic problems with test output.

``--blastn_exe BLASTN_EXE``
    Path to the ``blastn`` executable. Default: ``blastn``

``--format_exe FORMAT_EXE``
    Path to the ``blastn`` executable. Default: ``makeblastdb``

``--fragsize FRAGSIZE``
    Fragment size to use in analysis. (default: 1020)

``-h, --help``
    Display usage information for ``pyani anib``.

``--jobprefix JOBPREFIX``
    Use the string ``JOBPREFIX`` as a prefix for SGE job submission names. Default: ``PYANI``

``--labels LABELFNAME``
    Use the set of labels (one per genome sequence file) found in the file ``LABELFNAME`` in ``INDIR``. Default: ``labels.txt``

``-l LOGFILE, --logfile LOGFILE``
    Provide the location ``LOGFILE`` to which a logfile of the ANIb process will be written.

``--name NAME``
    Use the string ``NAME`` to identify this ``ANIb`` run in the ``pyani`` database.

``--recovery``
    Use existing ``ANIb`` comparison output if available, e.g. if recovering from a failed job submission. Using this option will not generate a new comparison if the old output files exist.

``--scheduler {multiprocessing, SGE}``
    Specify the job scheduler to be used when parallelising genome comparisons: one of ``multiprocessing`` (use many cores on the current machine)  or ``SGE`` (use an SGE or OGE job scheduler). Default: ``multiprocessing``.

``--SGEargs SGEARGS``
    Pass additional arguments ``SGEARGS`` to ``qsub`` when running the SGE-distributed jobs.

``--SGEgroupsize SGEGROUPSIZE``
    Create SGE arrays containing SGEGROUPSIZE comparison jobs. Default: 10000

``-v, --verbose``
    Provide verbose output to ``STDOUT``.

``--workers WORKERS``
    Spawn WORKERS worker processes with the ``--scheduler multiprocessing`` option. Default: 0 (use all cores)
