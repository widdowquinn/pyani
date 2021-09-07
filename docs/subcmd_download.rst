.. _pyani-subcmd-download:

==================
``pyani download``
==================

The ``download`` subcommand controls download of genome files from the `NCBI Assembly`_ database for input to ``pyani``.

.. code-block:: text

    usage: pyani.py download [-h] [-l LOGFILE] [-v] [--debug] [--disable_tqdm] [--version]
                         [--citation] [-o OUTDIR] -t TAXON --email EMAIL
                         [--api_key API_KEYPATH] [--retries RETRIES]
                         [--batchsize BATCHSIZE] [--timeout TIMEOUT] [-f]
                         [--noclobber] [--labels LABELFNAME] [--classes CLASSFNAME]
                         [--kraken] [--dry-run]

--------------------
Positional arguments
--------------------

``outdir``
    The ``outdir`` argument should be the path to a directory into which genome files will be downloaded. If the directory exists, a warning will be given and the download will not proceed, to avoid overwriting existing data. To force writing into an existing directory, use the ``-f`` option.

-----------------
Flagged arguments
-----------------

``--api_key PATH_TO_API_KEY``
    The program will attempt to use an NCBI API key (see `here <https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/>`_) located at ``PATH_TO_API_KEY``. Default: ``~/.ncbi/api_key``

``--batchsize BATCHSIZE``
    The download process will attempt to download assemblies in multiples of ``BATCHSIZE``. Default: 10000

``--classes CLASSFNAME``
    Write a set of labels (one per downloaded genome) to the file ``CLASSFNAME`` in ``outdir``. Default: ``classes.txt``

``--disable_tqdm``
    Disable the ``tqdm`` progress bar while the download process runs. This is useful when testing to avoid aesthetic problems with test output.

``--dry-run``
    Perform all actions of the download process except for downloading files.

``--email EMAIL``
    **COMPULSORY**. Provide the email address ``EMAIL`` to NCBI so that they can track problems.

``-f, --force``
    Force use of the ``OUTDIR`` directory when downloaded genomes, even if it already exists.

``-h, --help``
    Display usage information for ``pyani download``.

``--kraken``
    Add taxonomy information to the FASTA file headers of downloaded genomes. This allows the genomes to be readily used to construct databases for the `Kraken`_ software package.

``-l LOGFILE, --logfile LOGFILE``
    Provide the location ``LOGFILE`` to which a logfile of the download process will be written.

``--labels LABELFNAME``
    Write a set of labels (one per downloaded genome) to the file ``LABELFNAME`` in ``outdir``. Default: ``labels.txt``

``--noclobber``
    Do not overwrite individual files in the ``outdir`` directory, when used with ``-f``.

``-o OUTDIR, --outdir OUTDIR``
    The ``OUTDIR`` argument should be the path to a directory into which genome files will be downloaded. If the directory exists, a warning will be given and the download will not proceed, to avoid overwriting existing data. To force writing into an existing directory, use the ``-f`` option.

``--retries RETRIES``
    The download process will attempt to download each batch of assemblies a maximum of ``RETRIES`` times. Default: 20

``-t TAXON, --taxon TAXON``
    **COMPULSORY**. All genomes below taxon ID ``TAXON`` of a node in the `NCBI Taxonomy`_ database will be downloaded to the location specified by ``outdir``.

``--timeout TIMEOUT``
    The download process will wait a amaximum of ``TIMEOUT`` seconds before abandoning a URL connection attempt. Default: 10

``-v, --verbose``
    Provide verbose output to ``STDOUT``


.. _Kraken: https://ccb.jhu.edu/software/kraken/
.. _NCBI Assembly: https://www.ncbi.nlm.nih.gov/assembly
.. _NCBI Taxonomy: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
