.. _pyani-run_anim:

=====================
Running ANIm analysis
=====================

``pyani`` implements average nucleotide identity analysis using `MUMmer3`_ (*ANIm*) as defined in Richter & Rosselló-Móra (2009) (`doi:10.1073/pnas.0906412106`_). To run ANIm on a set of input genomes, use the ``pyani anim`` subcommand.

In brief, the analysis proceeds as follows for a set of input prokaryotic genomes:

1. `MUMmer3`_ is used to perform pairwise comparisons between each possible pair of input genomes, to identify homologous (alignable) regions.
2. For each comparison, the alignment output is parsed, and the following values are calculated:

    - total number of aligned bases on each genome
    - fraction of each genome that is aligned (the *coverage*)
    - the proportion of all aligned regions that is identical in each genome (the *ANI*)
    - the number of unaligned or non-identical bases (the *similarity errors*)
    - the product of *coverage* and *ANI*

The output values are recorded in the ``pyani`` database.

.. NOTE::
    Two ``Mummer`` comparisons are performed between each pair of genomes, so that each genome serves as the reference sequence, and as the subject sequence.
    .. A single ``MUMmer`` comparison is performed between each pair of genomes. Input genomes are sorted into alphabetical order by filename, and the query sequence is the genome that occurs earliest in the list; the subject sequence is the genome that occurs latest in the list.

.. TIP::
    The ``MUMmer`` comparisons are embarrasingly parallel, and can be distributed across cores on an `Open Grid Scheduler`_-compatible cluster, using the ``--scheduler SGE`` option.

.. ATTENTION::
    ``pyani anim`` requires that a working copy of `MUMmer3` is available. Please see :ref:`pyani-installation` for information about installing this package.

For more information about the ``pyani anim`` subcommand, please see the :ref:`pyani-subcmd-anim` page, or issue the command ``pyani anim -h`` to see the inline help.

---------------------
Perform ANIm analysis
---------------------

The basic form of the command is:

.. code-block:: bash

    pyani anim -i <INPUT_DIRECTORY> -o <OUTPUT_DIRECTORY>

This instructs ``pyani`` to perform ANIm on the genome FASTA files in ``<INPUT_DIRECTORY>``, which is passed to the ``-i`` argument, and write any output files to ``<OUTPUT_DIRECTORY>``, which is passed to the ``-o`` argument. For example, the following command performs ANIm on genomes in the directory ``genomes`` and writes output to a new directory ``genomes_ANIm``:

.. code-block:: bash

    pyani anim -i genomes -o genomes_ANIm

.. NOTE::
    While running, ``pyani anim`` will show progress bars unless these are disabled with the option ``--disable_tqdm``

This command will write the intermediate ``nucmer``/``MUMmer`` output to the directory ``genomes_ANIm``, in a subdirectory called ``nucmer_output``, where the results can be inspected if required.

.. code-block:: bash

    $ ls genomes_ANIm/
    nucmer_output

.. ATTENTION::
    To view the output ANIm results, you will need to use the ``pyani report`` or ``pyani plot`` subcommands. Please see :ref:`pyani-subcmd-report` and :ref:`pyani-subcmd-plot` for more details.

----------------------------------------------
Perform ANIm analysis with Open Grid Scheduler
----------------------------------------------

The ``MUMmer`` comparison step of ANIm is embarrassingly parallel, and ``nucmer`` jobs can be distributed across cores in a cluster using the `Open Grid Scheduler`. To enable this during the analysis, use the ``--scheduler SGE`` option:

.. code-block:: bash

    pyani anim --scheduler SGE -i genomes -o genomes_ANIm

.. NOTE::
    Jobs are submitted as *array jobs* to keep the scheduler queue short.

.. NOTE::
    If ``--scheduler SGE`` is not specified, all ``MUMmer`` jobs are run locally with ``Python``'s ``multiprocessing`` module.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Controlling parameters of Open Grid Scheduler
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to control the following features of `Open Grid Scheduler`_ `via` the ``pyani anim`` subcommand:

- The array job size (by default, comparison jobs are batched in arrays of 10,000)
- The prefix string for the job, as reported in the scheduler queue
- Arguments to the ``qsub`` job submission command

These allow for useful control of job execution. For example, the command:

.. code-block:: bash

    pyani anim --scheduler SGE --SGEgroupsize 5000 -i genomes -o genomes_ANIm

will batch ``MUMmer`` jobs in groups of 500 for the scheduler. The command:

.. code-block:: bash

    pyani anim --scheduler SGE --jobprefix My_Ace_Job -i genomes -o genomes_ANIm

will prepend the string ``My_Ace_Job`` to your job in the scheduler queue. And the command:

.. code-block:: bash

    pyani anim --scheduler SGE --SGEargs "-m e -M my.name@my.domain" 5000 -i genomes -o genomes_ANIm

will email ``my.name@my.domain`` when the jobs finish.


----------
References
----------

- Richter & Rosselló-Móra (2009) Proc Natl Acad Sci USA 106: 19126-19131 `doi:10.1073/pnas.0906412106`_.

.. _doi:10.1073/pnas.0906412106: https://dx.doi.org/10.1073/pnas.0906412106
.. _MUMmer3: http://mummer.sourceforge.net/
.. _Open Grid Scheduler: http://gridscheduler.sourceforge.net/
