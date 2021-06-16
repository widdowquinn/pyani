.. _pyani-run_fastani:

========================
Running fastANI analysis
========================

``pyani`` implements average nucleotide identity analysis using `fastANI`_(*fastANI*) as defined in __________________. To run fastANI on a set of input genomes, use the ``pyani fastani`` subcommand.

In brief, the analysis proceeds as follows for a set of input prokaryotic genomes:

1. `fastANI`_ is used to perform pairwise comparisons between each possible (directed) pair of input genomes, to identify homologous (alignable) regions.
2. For each comparison, the alignment output is parsed, and the following values are calculated:

    - total number of aligned bases on each genome
    - fraction of each genome that is aligned (the *coverage*)
    - the proportion of all aligned regions that is identical in each genome (the *ANI*)
    - the number of unaligned or non-identical bases (the *similarity errors*)
    - the product of *coverage* and *ANI*

The output values are recorded in the ``pyani`` database.

.. NOTE::
    ``fastANI`` comparisons are not symmetrical, so two comparisons are performed between each pair of genomes; the query and reference genomes from the first comparison are swapped for the second one.

.. ATTENTION::
    ``pyani fastANI`` requires that a working copy of `fastANI` is available. Please see :ref:`pyani-installation` for information about installing this package.

For more information about the ``pyani fastani`` subcommand, please see the :ref:`pyani-subcmd-fastani` page, or issue the command ``pyani fastani -h`` to see the inline help.

------------------------
Perform fastANI analysis
------------------------

The basic form of the command is:

.. code-block:: bash

    pyani fastani -i <INPUT_DIRECTORY> -o <OUTPUT_DIRECTORY>

This instructs ``pyani`` to perform fast ANI on the genome FASTA files in ``<INPUT_DIRECTORY>``, and write any output files to ``<OUTPUT_DIRECTORY>``. For example, the following command performs fastANI on genomes in the directory ``genomes`` and writes output to a new directory ``genoems_fastANI``:

.. code-block:: bash

    pyani fastani -i genomes -o genomes_fastANI

.. NOTE::
    While running, ``pyani fastani`` will show progress bars, unless these are disabled with the option ``--disable_tqdm``.

This command will write the intermediate ``fastANI`` output to the directory ``genomes_fastANI``, where the results can be inspected if required.

.. ATTENTION::
    To view the output fastANI results, you will need to use the ``pyani report`` or ``pyani plot`` subcommands. Please see :ref:`pyani-subcmd-report` and :ref:`pyani-subcd-plot` for more details.

-------------------------------------------------
Perform fastANI analysis with Open Grid Scheduler
-------------------------------------------------

To enable use of a scheduler during the analysis, use the ``--scheduler SGE`` option:

.. code-block:: bash

    pyani fastani --scheduler SGE -i genomes -o genomes_fastANI

.. NOTE::
    Jobs are submitted as *array jobs* to keep the scheduler queue short.

.. NOTE::
    If ``--scheduler SGE`` is not specified, all ``fastANI`` jobs are run locally with ``Python``'s ``multiprocessing`` module.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Controlling parameters of Open Grid Scheduler
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to control the following features of `Open Grid Scheduler`_ `via` the ``pyani fastani`` subcommand:

- The array job size (by default, comparison jobs are batched in arrays of 10,000)
- The prefix string for the job, as reported in the scheduler queue
- Arguments to the ``qsub`` job submission command

These allow for useful control of job execution. For example, the command:

.. code-block:: bash

    pyani fastani --scheduler SGE --SGEgroupsize 5000 -i genomes -o genomes_fastANI

will batch ``fastANI`` jobs in groups of 500 for the scheduler. The command:

.. code-block:: bash

    pyani fastani --scheduler SGE --jobprefix My_Ace_Job -i genomes -o genomes_ANIm

will prepend the string ``My_Ace_Job`` to your job in the scheduler queue. And the command:

.. code-block:: bash

    pyani fastani --scheduler SGE --SGEargs "-m e -M my.name@my.domain" 5000 -i genomes -o genomes_ANIm

will email ``my.name@my.domain`` when the jobs finish.


----------
References
----------

- Jain, C., Rodriguez-R, L.M., Phillippy, A.M. et al. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nat Commun 9, 5114 (2018). `doi.org/10.1038/s41467-018-07641-9`_.

.. _doi.org/10.1038/s41467-018-07641-9: https://doi.org/10.1038/s41467-018-07641-9
.. _Open Grid Scheduler: http://gridscheduler.sourceforge.net/
