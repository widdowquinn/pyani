.. _pyani-subcmd-index:

===============
``pyani index``
===============

The ``index`` subcommand will index the genome files it finds the passed directory ``INDIR``, generating label and class files, and files that contain an MD5 hash of the nucleotide sequence of each genome.

.. code-block:: text

    usage: pyani index [-h] [-l LOGFILE] [-v] [--debug] [--disable_tqdm] [--version]
                      [--citation] -i INDIR [--labels LABELFNAME]
                      [--classes CLASSFNAME]

-----------------
Flagged arguments
-----------------

``--classes CLASSFNAME``
    Write a set of labels (one per genome sequence file) to the file ``CLASSFNAME`` in ``INDIR``. Default: ``classes.txt``

``--disable_tqdm``
    Disable the ``tqdm`` progress bar while the download process runs. This is useful when testing to avoid aesthetic problems with test output.

``-h, --help``
    Display usage information for ``pyani index``.

``-i INDIR, --indir INDIR``
    The ``INDIR`` argument should be the path to a directory containing genome sequence data as FASTA files (one per genome assembly).

``-l LOGFILE, --logfile LOGFILE``
    Provide the location ``LOGFILE`` to which a logfile of the download process will be written.

``--labels LABELFNAME``
    Write a set of labels (one per genome sequence file) to the file ``LABELFNAME`` in ``INDIR``. Default: ``labels.txt``

``-v, --verbose``
    Provide verbose output to ``STDOUT``
