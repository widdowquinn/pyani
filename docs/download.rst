.. _pyani-download:

=============================
Downloading Genomes from NCBI
=============================

This page describes some typical use cases for downloading genomes from NCBI using the ``pyani download`` subcommand. This command downloads all assembled genomes found beneath a given taxon identifier, at the NCBI GenBank database. You will need to know the NCBI taxon ID for each taxon you wish to download.

.. ATTENTION::
    To use their online resources programmatically, NCBI require that you provide your email address for contact purposes if jobs go wrong, and for their own usage statistics. This should be specified with the ``--email <EMAIL ADDRESS>`` argument of ``pyani download``.

For more information about the ``pyani download`` subcommand, please see the :ref:`pyani-subcmd-download` page, or issue the command ``pyani download -h`` to see the inline help.

--------------------------------------
Download all genomes in a single taxon
--------------------------------------

The basic form of the command is:

.. code-block:: bash

    pyani download --email my.email@my.domain -t <TAXON_ID> <OUTPUT_DIRECTORY>

This instructs ``pyani`` to use the ``download`` subcommand to obtain all available genome assemblies below the taxon ID ``<TAXON_ID>``, passed with the ``-t`` argument, and place the downloaded files - along with label and class information files created by ``pyani`` in the subdirectory ``<OUTPUT_DIRECTORY>``.

For example, if we wished to download all available assemblies for the bacterium *Pseudomonas flexibilis* we would `identify the taxon ID <https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=706570&lvl=3&lin=f&keep=1&srchmode=1&unlock>`_ to be 706570, and use this as the argument to ``-t``, placing the output in a convenient subdirectory (e.g. ``genomes``), with the command:

.. code-block:: bash

    $ pyani download --email my.email@my.domain -t 706570 genomes
    GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_genomic.fna.gz: 2097152it [00:00, 3224293.90it/s]
    GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_hashes.txt: 1048576it [00:00, 110097041.36it/s]
    GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_genomic.fna.gz: 2097152it [00:00, 3125724.89it/s]
    GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_hashes.txt: 1048576it [00:00, 55139621.76it/s]
    GCA_001312105.1_ASM131210v1_genomic.fna.gz: 1048576it [00:00, 1519534.52it/s]
    GCA_001312105.1_ASM131210v1_hashes.txt: 1048576it [00:00, 272426072.29it/s]
    GCF_000806415.1_ASM80641v1_genomic.fna.gz: 2097152it [00:00, 3097667.41it/s]
    GCF_000806415.1_ASM80641v1_hashes.txt: 1048576it [00:00, 129632638.05it/s]
    GCF_000802425.1_ASM80242v1_genomic.fna.gz: 2097152it [00:00, 2973018.02it/s]
    GCF_000802425.1_ASM80242v1_hashes.txt: 1048576it [00:00, 233590743.10it/s]

This displays each assembly as a download is attempted, and places all output in the named subdirectory:

.. code-block:: bash

    $ tree genomes/
    genomes/
    ├── GCA_001312105.1_ASM131210v1_genomic.fna
    ├── GCA_001312105.1_ASM131210v1_genomic.fna.gz
    ├── GCA_001312105.1_ASM131210v1_genomic.md5
    ├── GCA_001312105.1_ASM131210v1_hashes.txt
    ├── GCF_000802425.1_ASM80242v1_genomic.fna
    ├── GCF_000802425.1_ASM80242v1_genomic.fna.gz
    ├── GCF_000802425.1_ASM80242v1_genomic.md5
    ├── GCF_000802425.1_ASM80242v1_hashes.txt
    ├── GCF_000806415.1_ASM80641v1_genomic.fna
    ├── GCF_000806415.1_ASM80641v1_genomic.fna.gz
    ├── GCF_000806415.1_ASM80641v1_genomic.md5
    ├── GCF_000806415.1_ASM80641v1_hashes.txt
    ├── GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_genomic.fna
    ├── GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_genomic.fna.gz
    ├── GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_genomic.md5
    ├── GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_hashes.txt
    ├── GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_genomic.fna
    ├── GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_genomic.fna.gz
    ├── GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_genomic.md5
    ├── GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_hashes.txt
    ├── classes.txt
    └── labels.txt

Each genome is downloaded in compressed format (``.fna.gz`` files) and expanded in-place to give the FASTA file (``.fna``) files. The MD5 hash of each FASTA file is also calculated (``.md5``). This will be used by ``pyani`` to uniquely identify that assembly throughout the analysis process.

.. code-block:: bash

    $ head genomes/GCA_001312105.1_ASM131210v1_genomic.md5
    e55cd3d913a198ac60afd8d509c02ab4	genomes/GCA_001312105.1_ASM131210v1_genomic.fna

``pyani`` also creates two files:

- ``classes.txt``: each genome is assigned a *class* which is used to annotate genomes in the graphical output. ``pyani`` attempts to infer genus and species as the default class
- ``labels.txt``: each genome is assigned a text label, which is used to label genomes in the graphical output. ``pyani`` attempts to infer genus, species, and strain ID as the default label

.. code-block:: bash

    $ head genomes/classes.txt
    1ac4941c569f32f044eba0a8540d4704	GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_genomic	Pseudomonas flexibilis
    8664341798070f1d70b2569a5b3a2320	GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_genomic	Pseudomonas flexibilis
    e55cd3d913a198ac60afd8d509c02ab4	GCA_001312105.1_ASM131210v1_genomic	Pseudomonas flexibilis
    9b9719eb78bf7cf6dd0146a3f9426f60	GCF_000806415.1_ASM80641v1_genomic	Pseudomonas flexibilis
    2bdfffd867d843f970e4dfd388d5332a	GCF_000802425.1_ASM80242v1_genomic	Pseudomonas flexibilis
    $ head genomes/labels.txt
    1ac4941c569f32f044eba0a8540d4704	GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_genomic	P. flexibilis ATCC 29606
    8664341798070f1d70b2569a5b3a2320	GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_genomic	P. flexibilis CGMCC 1.1365
    e55cd3d913a198ac60afd8d509c02ab4	GCA_001312105.1_ASM131210v1_genomic	P. flexibilis JCM 14085
    9b9719eb78bf7cf6dd0146a3f9426f60	GCF_000806415.1_ASM80641v1_genomic	P. flexibilis JCM 14085
    2bdfffd867d843f970e4dfd388d5332a	GCF_000802425.1_ASM80242v1_genomic	P. flexibilis ATCC 29606

These files are used to associate labels and classes to the genome files in the ``pyani`` database, specific to the analysis run.

Both ``classes.txt`` and ``labels.txt`` can be edited to suit the user's classification and labelling scheme.

---------------------------------------
Download all genomes from multiple taxa
---------------------------------------

To download genomes from more than one taxon, you can provide a comma-separated list of taxon IDs to the ``pyani download`` subcommand, e.g.:

.. code-block:: bash

    pyani download --email my.email@my.domain -t <TAXON_ID1>,<TAXON_ID2>,... <OUTPUT_DIRECTORY>

The following command can be used to download assemblies from three different *Pseudomonas* taxa (*P. flexibilis*: 706570, *P. mosselli*: 78327, and *P. fulva*: 47880):

.. code-block:: bash

    $ pyani download --email my.email@my.domain -t 706570,78327,47880 multi_taxa
    GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_genomic.fna.gz: 2097152it [00:00, 3081776.59it/s]
    GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_hashes.txt: 1048576it [00:00, 63489526.95it/s]
    GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_genomic.fna.gz: 2097152it [00:00, 3194885.99it/s]
    GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_hashes.txt: 1048576it [00:00, 77838775.82it/s]

---------------------------------------------------
Dry-run test (identify, but do not download, files)
---------------------------------------------------

If you only want to see which genomes will be downloaded from NCBI with a given ``pyani download`` subcommand, but not download them, then you can use the ``--dry-run`` option. For example:

.. code-block:: bash

    $ pyani download --email my.email@my.domain -t 706570,78327,47880 multi_taxa --dry-run
    WARNING: Dry run only: will not overwrite or download
    WARNING: (dry-run) skipping download of GCF_900155995.1
    WARNING: (dry-run) skipping download of GCF_900101515.1
    WARNING: (dry-run) skipping download of GCA_001312105.1
    WARNING: (dry-run) skipping download of GCF_000806415.1
    [...]

----------------------------------------------------------------
Download genomes for compilation of a custom ``Kraken`` database
----------------------------------------------------------------

`Kraken`_ is a bioinformatics tool that assigns taxonomic identities to short DNA sequences, such as Illumina or Nanopore reads. Several wide-ranging ``Kraken`` databases are available for download, and around the community, but it can sometimes be useful to construct a custom local ``Kraken`` database specific for your organism or taxon of interest (e.g. for filtering out contaminating or suspect reads).

The ``pyani download`` command can prepare downloaded genome files for immediate use with the ``Kraken`` database-building tools, by specifying the ``--kraken`` option:

.. code-block:: bash

    $ pyani download --email my.email@my.domain -t 706570,78327,47880 genomes_kraken --kraken
    GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_genomic.fna.gz: 2097152it [00:00, 3085741.03it/s]
    GCF_900155995.1_IMG-taxon_2681812811_annotated_assembly_hashes.txt: 1048576it [00:00, 140958511.30it/s]
    WARNING: Modifying downloaded sequence for Kraken compatibility
    GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_genomic.fna.gz: 2097152it [00:01, 1902023.68it/s]
    GCF_900101515.1_IMG-taxon_2596583557_annotated_assembly_hashes.txt: 1048576it [00:00, 31428985.47it/s]
    WARNING: Modifying downloaded sequence for Kraken compatibility
    [...]
    $ head multi_taxa/GCA_001312105.1_ASM131210v1_genomic.fna
    >BBCY01000001.1 Pseudomonas tuomuerensis JCM 14085 DNA, contig: JCM14085.contig00001, whole genome shotgun sequence
    ACCAGCATCTGGCGGATCAGGTCGCGGGCCTTCTCGGCCGATTGGCGGATGCGCCCGAGGTAGCGGCCGAGCGGCGCGTC
    GCCGCGCTCGCCCGCCAGCTCCTCGGCCATCTGCGTGTAGCCGAGCATGCTGGTCAGCAGGTTGTTGAAGTCGTGGGCAA
    TGCCGCCGGTCAGGTGGCCGATGGCTTCCATGCGCTGCGCCTGGCGCAGCTGCTGTTCCAGCGCCGCCCGCTCCACCTCG
    $ head genomes_kraken/GCA_001312105.1_ASM131210v1_genomic.fna
    >BBCY01000001.1|kraken:taxid|706570 BBCY01000001.1 Pseudomonas tuomuerensis JCM 14085 DNA, contig: JCM14085.contig00001, whole genome shotgun sequence
    ACCAGCATCTGGCGGATCAGGTCGCGGGCCTTCTCGGCCGATTGGCGGATGCGCCCGAGG
    TAGCGGCCGAGCGGCGCGTCGCCGCGCTCGCCCGCCAGCTCCTCGGCCATCTGCGTGTAG
    CCGAGCATGCTGGTCAGCAGGTTGTTGAAGTCGTGGGCAATGCCGCCGGTCAGGTGGCCG
    ATGGCTTCCATGCGCTGCGCCTGGCGCAGCTGCTGTTCCAGCGCCGCCCGCTCCACCTCG

Using this option does affects downstream performance or use of ``pyani`` only in that the two different input files for the same genome will have distinct hashes:

.. code-block:: bash

    $ head multi_taxa/GCA_001312105.1_ASM131210v1_genomic.md5
    e55cd3d913a198ac60afd8d509c02ab4	multi_taxa/GCA_001312105.1_ASM131210v1_genomic.fna
    $ head genomes_kraken/GCA_001312105.1_ASM131210v1_genomic.md5
    053fd98d8c9ab30de46f56fd601ef529	genomes_kraken/GCA_001312105.1_ASM131210v1_genomic.fna

and so will not be considered to be the "same sequence" when repeating comparisons.

.. _Kraken: https://ccb.jhu.edu/software/kraken/
