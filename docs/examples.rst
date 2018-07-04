.. _pyani-examples:

========
Examples
========

----------------------
Using non-NCBI genomes
----------------------

It is usual to want to include or work only with genomes that have been generated locally, or that were not downloaded from NCBI using
`pyani download`. To use these genomes with the `pyani` analysis subcommands, the genomes must be *indexed* [1]_.

To *index* a set of genomes, use the `pyani index` subcommand on the input directory. To index the directory ``mygenomes``, for example:

.. code-block:: bash

    pyani index mygenomes

This will create a ``.md5`` file (containing the *hash*) for each genome, as well as class and label files listing all the input genomes.

All genomes in the ``mygenomes`` directory will now be available for use in `pyani`.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``class`` and ``labels`` files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using the ``pyani download`` command will create two files, by default ``classes.txt`` and ``labels.txt`` containing identifiers for each input genome that are used in later analysis and visualisation. These files are also created when ``pyani index`` is used as above.

The location of the labels and classes files may be changed using the ``--labels`` and ``--classes`` arguments, for example:

.. code-block:: bash

    pyani index mygenomes --classes myclasses.txt --labels mylabels.txt




.. [1] *indexing* here refers to constructing a *hash* of the genome: a short representation of the entire genome's contents that can be used to identify it uniquely