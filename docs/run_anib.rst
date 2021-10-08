.. _pyani-run_anib:

=====================
Running ANIb analysis
=====================

``pyani`` implements average nucleotide identity analysis using `NCBI-BLAST+`_ (*ANIb*) as defined in Goris `et al.` (2007) (`doi:10.1099/ijs.0.64483-0`_). To run ANIb on a set of input genomes, use the ``pyani anib`` subcommand.

In brief, the analysis proceeds as follows for a set of input prokaryotic genomes:

1. Each input genome is fragmented into consecutive sequences of a given size (default: 1020bp)
2. A new ``BLAST+`` database is built from each input genome sequence
3. `NCBI-BLAST+`_ is used to perform pairwise comparisons of each input genome fragment set against the databases for each other input genome, to identify homologous (alignable) regions.
4. For each comparison, the alignment output is parsed, and the following values are calculated:

    - total number of aligned bases on each genome
    - fraction of each genome that is aligned (the *coverage*)
    - the proportion of all aligned regions that is identical in each genome (the *ANI*)
    - the number of unaligned or non-identical bases (the *similarity errors*)
    - the product of *coverage* and *ANI*

The output values are recorded in the ``pyani`` database.

.. NOTE::
    The `NCBI-BLAST+`_ comparisons are asymmetric, and performed in both directions for a pair of genomes (i.e. "fragmented A vs complete B" and "fragmented B vs complete A").

.. TIP::
    The `NCBI-BLAST+`_ comparisons are embarrasingly parallel, and can be distributed across cores on an `Open Grid Scheduler`_-compatible cluster, using the ``--scheduler SGE`` option.

.. ATTENTION::
    ``pyani anib`` requires that a working copy of `NCBI-BLAST+`_ is available. Please see :ref:`pyani-installation` for information about installing this package.

For more information about the ``pyani anib`` subcommand, please see the :ref:`pyani-subcmd-anib` page, or issue the command ``pyani anib -h`` to see the inline help.

---------------------
Perform ANIb analysis
---------------------

The basic form of the command is:

.. code-block:: bash

    pyani anib -i <INPUT_DIRECTORY> -o <OUTPUT_DIRECTORY>

This instructs ``pyani`` to perform ANIb on the genome FASTA files in ``<INPUT_DIRECTORY>``, which is passed to the ``-i`` argument, and write any output files to ``<OUTPUT_DIRECTORY>``, which is passed to the ``-o`` argument. For example, the following command performs ANIb on genomes in the directory ``genomes`` and writes output to a new directory ``genomes_ANIb``:

.. code-block:: bash

    pyani anib -i genomes -o genomes_ANIb

.. NOTE::
    While running, ``pyani anib`` will show progress bars unless these are disabled with the option ``--disable_tqdm``

This command will write the intermediate `NCBI-BLAST+`_ output to the directory ``genomes_ANIb``, where the results can be inspected if required.

..
  I am unsure if this is relevant for anib
  .. code-block:: bash

    $ ls genomes_ANIb/
    blastn_output

.. ATTENTION::
    To view the output ANIb results, you will need to use the ``pyani report`` or ``pyani plot`` subcommands. Please see :ref:`pyani-subcmd-report` and :ref:`pyani-subcmd-plot` for more details.

----------------------------------------------
Perform ANIb analysis with Open Grid Scheduler
----------------------------------------------

The `NCBI-BLAST+`_ comparisons are embarrassingly parallel, and `NCBI-BLAST+`_ jobs can be distributed across cores in a cluster using the `Open Grid Scheduler`_. To enable this during the analysis, use the ``--scheduler SGE`` option:

.. code-block:: bash

    pyani anib --scheduler SGE -i genomes -o genomes_ANIb

.. NOTE::
    Jobs are submitted as *array jobs* to keep the scheduler queue short.

.. NOTE::
    If ``--scheduler SGE`` is not specified, all `NCBI-BLAST+`_ jobs are run locally with ``Python``'s ``multiprocessing`` module.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Controlling parameters of Open Grid Scheduler
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to control the following features of `Open Grid Scheduler`_ `via` the ``pyani anib`` subcommand:

- The array job size (by default, comparison jobs are batched in arrays of 10,000)
- The prefix string for the job, as reported in the scheduler queue
- Arguments to the ``qsub`` job submission command

These allow for useful control of job execution. For example, the command:

.. code-block:: bash

    pyani anib --scheduler SGE --SGEgroupsize 5000 -i genomes -o genomes_ANIb

will batch `NCBI-BLAST+`_ jobs in groups of 500 for the scheduler. The command:

.. code-block:: bash

    pyani anib --scheduler SGE --jobprefix My_Ace_Job -i genomes -o genomes_ANIb

will prepend the string ``My_Ace_Job`` to your job in the scheduler queue. And the command:

.. code-block:: bash

    pyani anib --scheduler SGE --SGEargs "-m e -M my.name@my.domain" 5000 -i genomes -o genomes_ANIb

will email ``my.name@my.domain`` when the jobs finish.


----------
References
----------

- Goris`et al.` (2007) `Int J Syst Evol Micr` _57_: 81-91. `doi:10.1099/ijs.0.64483-0`.

.. _doi:10.1099/ijs.0.64483-0: https://dx.doi.org/10.1099/ijs.0.64483-0
.. _NCBI-BLAST+: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
.. _Open Grid Scheduler: http://gridscheduler.sourceforge.net/
