.. _pyani-quickstart:

================
QuickStart Guide
================

------------
Installation
------------

To use ``pyani`` you will need to install it. Installation is easiest using one of the two most popular Python package managers:

^^^^^^^^^^^^
1. ``conda``
^^^^^^^^^^^^

``pyani`` is available through the ``bioconda`` channel of Anaconda:

.. code-block:: bash

    conda install -c bioconda pyani 


^^^^^^^^^^^
2. ``PyPI``
^^^^^^^^^^^

``pyani`` is available *via* the `PyPI`_ package manager for Python:

.. code-block:: bash

    pip install pyani 

.. TIP::
    ``pyani`` can also be installed directly from source. Installation instructions can be found on the :ref:`pyani-installation` page.


---------------------
``pyani`` Walkthrough
---------------------

An example ANIm analysis using ``pyani`` is provided as a walkthrough below. The general procedure for any ``pyani`` analysis is:

1. Download genomes for the analysis
2. Create a database to hold genome data and analysis results
3. Perform ANIm analysis
4. Report and visualise analysis results
5. Generate species hypotheses (classify genomes) using the analysis results


To see options available for the ``pyani`` program, use the ``-h``
(help) option:

.. code-block:: bash

    pyani -h

^^^^^^^^^^^^^^^^^^^
1. Download genomes
^^^^^^^^^^^^^^^^^^^

For this walkthrough a new set of genomes will be obtained from GenBank, using the ``pyani download`` command.

.. TIP::
    ``pyani`` requires an *indexed* set of genomes, and the visualisation and classification steps benefit from having ``classes.txt`` and ``labels.txt`` files. These are generated automatically when downloading genomes, but you must create them in other ways when applying ``pyani`` to a set of local files.

.. ATTENTION::
    To use their online resources programmatically, NCBI require that you provide your email address for contact purposes if jobs go wrong, and for their own usage statistics. This is specified with the ``--email <EMAIL ADDRESS>`` argument of ``pyani download``.

.. ATTENTION::
    ``pyani`` can also be used with an existing local set of genomes. This is discussed elsewhere in the documentation.

Using the ``pyani download`` subcommand, we download all available genomes for *Candidatus Blochmannia* from NCBI. The taxon ID for this grouping is ``203804``.

.. code-block:: bash

    pyani download C_blochmannia --email my.email@my.domain -t 203804

This creates a new directory (``C_blochmannia``) with the following contents:

.. code-block:: bash

    $ tree C_blochmannia
    C_blochmannia
    ├── GCF_000011745.1_ASM1174v1_genomic.fna
    ├── GCF_000011745.1_ASM1174v1_genomic.fna.gz
    ├── GCF_000011745.1_ASM1174v1_genomic.md5
    [...]
    ├── GCF_000973545.1_ASM97354v1_hashes.txt
    ├── classes.txt
    └── labels.txt

Each downloaded genome is represented by four files: the genome sequence (expanded: ``.fna``, compressed: ``.fna.gz``), an NCBI hashes file (``_hashes.txt``) and an MD5 hash of the genome sequence file (``.md5``).

Two additional files are created: 

- ``classes.txt``: defines a *class* to which each input genome belongs. This is used for determining membership of groups and annotating graphical output.
- ``labels.txt``: provides text which will be used to label each input genome in the graphical output from ``pyani``

^^^^^^^^^^^^^^^^^^
2. Create database
^^^^^^^^^^^^^^^^^^

``pyani`` uses a database to store genome data and analysis results. For this walkthrough, we create a new, empty database by executing the command:

.. code-block:: bash

    pyani createdb

.. TIP::
    This creates a new database in the default location (``.pyani/pyanidb``), but the name and location of this database can be controlled with the ``pyani createdb`` command, and specified in each of the subsequent commands.

^^^^^^^^^^^^^^^^^^^^^^^^
3. Conduct ANIm analysis
^^^^^^^^^^^^^^^^^^^^^^^^

We run ANIm on the downloaded genomes, by issuing the command:

.. code-block:: bash

    pyani anim C_blochmannia C_blochmannia_ANIm \
        --name "C. blochmannia run 1" \
        --labels C_blochmannia/labels.txt --classes C_blochmannia/classes.txt

This will run an ANIm analysis on the genomes in the ``C_blochmannia`` directory. The analysis results will be stored in the database we created earlier, identified by the name ``C. blochmannia run 1``. The comparison files will be written to the ``C_blochmannia_ANIm`` directory.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
4. Reporting Analyses and Analysis Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can list all the runs contained in the (default) database by using the command:

.. code-block:: bash

    pyani report --runs C_blochmannia_ANIm/ --formats html,excel,stdout

This will report the relevant information to new files in the ``C_blochmannia_ANIm`` directory.

.. code-block:: bash

    $ tree -L 1 C_blochmannia_ANIm/
    C_blochmannia_ANIm/
    ├── nucmer_output
    ├── runs.html
    ├── runs.tab
    └── runs.xlsx

.. TIP::
    By default the ``pyani report`` command will create a tab-separated text file with the ``.tab`` suffix, but by using the ``--formats`` option, we have also created an HTML file, and an Excel file with the same data. The ``stdout`` option prints the output table to the terminal window.

By inspecting the ``runs.tab` file (or any of the other ``runs.*`` files) we see that our walkthrough analysis has run ID ``1``. So we can use this ID to get tables of specific information for that run, such as:

**the genomes that were analysed in the run**

.. code-block:: bash

    pyani report --runs_genomes --formats html,excel,stdout C_blochmannia_ANIm/

**the complete set of pairwise comparison results** (listed by comparison)

.. code-block:: bash

    pyani report --run_results 1 --formats html,excel,stdout C_blochmannia_ANIm/

**comparison results as matrices** (percentage identity and coverage, number of aligned bases and "similarity errors", and a Hadamard matrix of identity multiplied by coverage).

.. code-block:: bash

    pyani report --run_matrices 1 --formats html,excel,stdout C_blochmannia_ANIm/

.. ATTENTION::
    The ``--run_results`` and ``--run_matrices`` options take the run ID (or a comma-separated list of IDs, such as ``1,3,4,5,9``) as an argument, and produce output for each run.

Graphical output is obtained by executing the ``pyani plot`` subcommand and specifying the output directory and run ID.

.. code-block:: bash

    pyani plot --formats png,pdf --method seaborn C_blochmannia_ANIm 1

Supported output methods are:

- ``seaborn``
- ``mpl`` (``matplotlib``)
- ``plotly``

and each generates five plots corresponding to the matrices that ``pyani report`` produces: percentage identity and coverage, number of aligned bases and "similarity errors", and a Hadamard matrix of identity multiplied by coverage.

.. figure:: images/matrix_identity_1.png
    :alt: percentage identity matrix for *Candidatus Blochmannia* ANIm analysis
    
    Percentage identity matrix for *Candidatus Blochmannia* ANIm analysis

    Each cell represents a pairwise comparison between the named genomes on rows and columns, and the number in the cell is the pairwise identity *of aligned regions*. The dendrograms are single-linkage clustering trees generated from the matrix of pairwise identity results. The default colour scheme colours cells with identity > 0.95 as red, and those with < 0.95 as blue. This division corresponds to the convention for bacterial species boundaries.

.. figure:: images/matrix_coverage_1.png
    :alt: percentage coverage matrix for *Candidatus Blochmannia* ANIm analysis
    
    Percentage coverage matrix for *Candidatus Blochmannia* ANIm analysis

    Each cell represents a pairwise comparison between the named genomes on rows and columns, and the number in the cell is pairwise coverage of the alignment. The dendrograms are single-linkage clustering trees generated from the matrix of pairwise coverage results. The default colour scheme colours cells with identity > 0.50 as red, and those with < 0.50 as blue. This division corresponds to a strict majority of each genome in the comparison being alignable (a plausible minimum requirement for being considered "the same thing").

Several graphics output formats are available, including ``.png``, ``.pdf`` and ``.svg``.



.. _NCBI Taxonomy database: https://www.ncbi.nlm.nih.gov/taxonomy
.. _PyPI: https://pypi.python.org/pypi