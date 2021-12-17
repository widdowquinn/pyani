.. _pyani-indexing:

================
Indexing Genomes
================

*Indexing* genomes is a necessary step in ``pyani`` to prepare input genomes for analysis.

.. ATTENTION::
    If you use the ``pyani download`` subcommand (see :ref:`pyani-download`) to obtain genomes for analysis, then indexing is carried out automatically. However, if you collect a local set of genomes (e.g. from your own sequencing project), then you will need to index the genomes with the ``pyani index`` subcommand.

For more information about the ``pyani index`` subcommand, please see the :ref:`pyani-subcmd-index` page, or issue the command ``pyani index -h`` to see the inline help.

----------------------
What does indexing do?
----------------------

In the context of ``pyani``, *indexing* refers to generating an index code that is unique to each input genome FASTA file in the input directory. The index code is the MD5 hash for the FASTA file.

This MD5 index code is used to identify each specific input genome sequence (and associated metadata) so that duplicate comparisons can be readily identified, and previous results reused from the ``pyani`` database, if they are available.

Indexing also generates two files (see :ref:`pyani-download`):

- ``classes.txt``: each genome is assigned a *class* which is used to annotate genomes in the graphical output. ``pyani`` attempts to infer genus and species as the default class
- ``labels.txt``: each genome is assigned a text label, which is used to label genomes in the graphical output. ``pyani`` attempts to infer genus, species, and strain ID as the default label

These files are used to associate labels and classes to the genome files in the ``pyani`` database, specific to the analysis run. Both ``classes.txt`` and ``labels.txt`` can be edited to suit the user's classification and labelling scheme.

--------------------------------
Index a directory of FASTA files
--------------------------------

The basic form of the command is:

.. code-block:: bash

    pyani index -i <GENOME_DIRECTORY>

This instructs ``pyani`` to search ``<GENOME_DIRECTORY>`` for files with a standard FASTA suffix (``.fna``, ``.fasta``, ``.fa``, ``.fas``, ``.fsa_nt``). For each file found, it calculates the MD5 hash and writes it to an accompanying file with extension ``.md5``. The hash is then associated with a genome label and a genome class, written to the two files ``labels.txt`` and ``classes.txt`` (see above).

For example, if we have a directory called ``unindexed`` that contains some FASTA format genome sequence files:

.. code-block:: bash

    $ tree unindexed
    unindexed
    ├── GCA_001312105.1_ASM131210v1_genomic.fna
    ├── GCF_000834555.1_ASM83455v1_genomic.fna
    └── GCF_005796105.1_ASM579610v1_genomic.fna

We could run the ``pyani index`` command:

.. code-block:: bash

    $ pyani index -i unindexed/
    $ tree unindexed
    unindexed
    ├── GCA_001312105.1_ASM131210v1_genomic.fna
    ├── GCA_001312105.1_ASM131210v1_genomic.fna.md5
    ├── GCF_000834555.1_ASM83455v1_genomic.fna
    ├── GCF_000834555.1_ASM83455v1_genomic.fna.md5
    ├── GCF_005796105.1_ASM579610v1_genomic.fna
    ├── GCF_005796105.1_ASM579610v1_genomic.fna.md5
    ├── classes.txt
    └── labels.txt

This creates an ``.md5`` file for each genome, and corresponding ``classes.txt`` and ``labels.txt`` files:

.. code-block:: bash

    $ head unindexed/GCA_001312105.1_ASM131210v1_genomic.fna
    >BBCY01000001.1 Pseudomonas tuomuerensis JCM 14085 DNA, contig: JCM14085.contig00001, whole genome shotgun sequence
    ACCAGCATCTGGCGGATCAGGTCGCGGGCCTTCTCGGCCGATTGGCGGATGCGCCCGAGGTAGCGGCCGAGCGGCGCGTC
    GCCGCGCTCGCCCGCCAGCTCCTCGGCCATCTGCGTGTAGCCGAGCATGCTGGTCAGCAGGTTGTTGAAGTCGTGGGCAA
    $ head unindexed/GCA_001312105.1_ASM131210v1_genomic.fna.md5
    e55cd3d913a198ac60afd8d509c02ab4	unindexed/GCA_001312105.1_ASM131210v1_genomic.fna
    $ head unindexed/classes.txt
    527f35b3eb9dd371d8d5309b6043dd9f	GCF_000834555.1_ASM83455v1_genomic	Pseudomonas fulva strain MEJ086 contig_1, whole genome shotgun sequence
    b00c5b1f636b8083b68b128e7ee28a40	GCF_005796105.1_ASM579610v1_genomic	Pseudomonas mosselii strain SC006 Scaffold1, whole genome shotgun sequence
    e55cd3d913a198ac60afd8d509c02ab4	GCA_001312105.1_ASM131210v1_genomic	Pseudomonas tuomuerensis JCM 14085 DNA, contig: JCM14085.contig00001, whole genome shotgun sequence
    $ head unindexed/labels.txt
    527f35b3eb9dd371d8d5309b6043dd9f	GCF_000834555.1_ASM83455v1_genomic	Pseudomonas fulva strain MEJ086 contig_1, whole genome shotgun sequence
    b00c5b1f636b8083b68b128e7ee28a40	GCF_005796105.1_ASM579610v1_genomic	Pseudomonas mosselii strain SC006 Scaffold1, whole genome shotgun sequence
    e55cd3d913a198ac60afd8d509c02ab4	GCA_001312105.1_ASM131210v1_genomic	Pseudomonas tuomuerensis JCM 14085 DNA, contig: JCM14085.contig00001, whole genome shotgun sequence

.. TIP::
    The class and label information produced by ``pyani index`` is different to that generated with ``pyani download``. Genus, species and strain identifiers can reliably be obtained from NCBI metadata when downloading genomes, but with user-provided sequences the information may not be encoded in the sequence description line in a standard manner.

    As a result, when using ``pyani index`` it is often useful to edit the ``classes.txt`` and ``labels.txt`` directly, or generate these files in some other way.
