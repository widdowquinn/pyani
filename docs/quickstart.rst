.. _pyani-quickstart:

================
QuickStart Guide
================

------------
Installation
------------

Installing `pyani` is easiest using one of the two most popular Python package managers:

^^^^^^^^^^^^
1. ``conda``
^^^^^^^^^^^^

``pyani`` is available through the ``bioconda`` channel of `Anaconda`_:

.. code-block:: bash

    conda install -c bioconda pyani

Installing pyani this way will also install all the required python packages and third-party tools, such as ``blast`` and ``mummer``.


^^^^^^^^^^^
2. ``PyPI``
^^^^^^^^^^^

``pyani`` is available *via* the `PyPI`_ package manager for Python:

.. code-block:: bash

    pip install pyani

.. ATTENTION::
    If you install ``pyani`` with ``pip``, then you will also need to install third-party tools, such as ``blast`` and ``mummer``, yourself. More detailed installation instructions can be found on the :ref:`pyani-installation` page.

.. TIP::
    ``pyani`` can also be installed directly from source. More detailed installation instructions can be found on the :ref:`pyani-installation` page.


---------------------
``pyani`` Walkthrough
---------------------

An example ANIm analysis using ``pyani`` is provided as a walkthrough to help you get started, below. The general procedure for any ``pyani`` analysis is:

1. Collect genomes for analysis (and index them)43
2. Create a database to hold genome data and analysis results
3. Perform ANIm analysis
4. Report and visualise analysis results
5. Generate species hypotheses (classify the input genomes) using the analysis results


To see options available for the ``pyani`` program, use the ``-h``
(help) option:

.. code-block:: bash

    pyani -h

^^^^^^^^^^^^^^^^^^
1. Collect genomes
^^^^^^^^^^^^^^^^^^

It is possible to use genomes you already have locally on your computer with ``pyani`` but, to illustrate downloading genomes for this walkthrough, a new set of files will be obtained from GenBank using the ``pyani download`` command.

.. TIP::
    To read more about using local files with ``pyani``, please see the :ref:`pyani-indexing` documentation. To read more about downloading genomes from NCBI, please see the :ref:`pyani-download` documentation.

.. ATTENTION::
    To use their online resources programmatically, NCBI require that you provide your email address for contact purposes if jobs go wrong, and for their own usage statistics. This can be specified with the ``--email <EMAIL ADDRESS>`` argument of ``pyani download``.

Using the ``pyani download`` subcommand, we download all available genomes for *Candidatus Blochmannia* from `NCBI <https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=203804&lvl=3&lin=f&keep=1&srchmode=1&unlock>`_. The taxon ID for this grouping is ``203804``, and this ID is passed as the ``-t`` argument. The final (compulsory) argument is the path to the directory into which the genome data will be downloaded.

.. code-block:: bash

    pyani download --email my.email@my.domain -t 203804 C_blochmannia

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

Each downloaded genome is represented by four files: the genome sequence (expanded: ``*.fna``, compressed: ``*.fna.gz``), an NCBI hashes file (``*_hashes.txt``) and an MD5 hash of the genome sequence file (``*.md5``). Two additional files are created:

- ``classes.txt``: defines a *class* to which each input genome belongs. This is used for determining membership of groups for each genome, and annotating graphical output.
- ``labels.txt``: provides text which will be used to label each input genome in the graphical output from ``pyani``

^^^^^^^^^^^^^^^^^^
2. Create database
^^^^^^^^^^^^^^^^^^

``pyani`` uses a local `SQLite3`_ database to store genome data and analysis results. For this walkthrough, we create a new, empty database by executing the command:

.. code-block:: bash

    pyani createdb

.. TIP::
    This creates the new database in a default location (``.pyani/pyanidb``), but the name and location of this database can be controlled with the ``pyani createdb`` command (see the :ref:`pyani-createdb` documentation). The path to the database can be specified in each of the subsequent commands, to enable maintenance and sharing of multiple analysis runs.

^^^^^^^^^^^^^^^^^^^^^^^^
3. Conduct ANIm analysis
^^^^^^^^^^^^^^^^^^^^^^^^

We can now run ANIm on the downloaded genomes and place the output in the database, by specifying paths first to the directory containing the genome data (here, ``C_blochmannia``), then to a directory which will contain the analysis results (``C_blochmannia_ANIm`` for this walkthrough.

We also provide a name for the analysis (``--name``) for later human-readable reference, and optional files defining labels for each genome to be used when plotting output (``--labels``) and a set of classes to which each genome belongs (``--classes``) for downstream analysis.:

.. code-block:: bash

    pyani anim C_blochmannia C_blochmannia_ANIm \
        --name "C. blochmannia run 1" \
        --labels C_blochmannia/labels.txt --classes C_blochmannia/classes.txt

This will run ANIm analysis on all genomes in the ``C_blochmannia`` directory. Analysis results will be stored in the database we created earlier (``.pyani/pyanidb``), where they will be identified by the name ``C. blochmannia run 1``. The comparison result files will be written to the ``C_blochmannia_ANIm`` directory.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
4. Reporting Analyses and Analysis Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can list all the runs contained in the (default) database, using the command:

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
    By default the ``pyani report`` command will create a tab-separated text file with the ``.tab`` suffix but, by using the ``--formats`` option, we have also created an HTML file, and an Excel file with the same data. The ``stdout`` option prints the output table to the terminal window.

By inspecting the ``runs.tab` file (or any of the other ``runs.*`` files) we see that our walkthrough analysis has run ID ``1``. So we can use this ID to get tables of specific information for that run, such as:

**the genomes that were analysed in all runs**

.. code-block:: bash

    pyani report --runs_genomes --formats html,excel,stdout C_blochmannia_ANIm/

**the complete set of pairwise comparison results for a single run** (listed by comparison)

.. code-block:: bash

    pyani report --run_results 1 --formats html,excel,stdout C_blochmannia_ANIm/

**comparison results as matrices** (percentage identity and coverage, number of aligned bases and "similarity errors", and a Hadamard matrix of identity multiplied by coverage).

.. code-block:: bash

    pyani report --run_matrices 1 --formats html,excel,stdout C_blochmannia_ANIm/

.. ATTENTION::
    The ``--run_results`` and ``--run_matrices`` options take the run ID (or a comma-separated list of IDs, such as ``1,3,4,5,9``) as an argument, and produce output for each run.

""""""""""""""""
Graphical output
""""""""""""""""

Graphical output is obtained by executing the ``pyani plot`` subcommand, specifying the output directory and run ID. Optionally, output file formats and the graphics drawing method can be specified.

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

    Each cell represents a pairwise comparison between the named genomes on rows and columns, and the number in the cell is the pairwise identity *of aligned regions* (i.e. regions with detectable homology). The dendrograms are single-linkage clustering trees generated from the matrix of pairwise identity results. The default colour scheme colours cells with identity > 0.95 as red, and those with < 0.95 as blue. This division corresponds to a widely-used convention for bacterial species boundaries.

.. figure:: images/matrix_coverage_1.png
    :alt: percentage coverage matrix for *Candidatus Blochmannia* ANIm analysis

    Percentage coverage matrix for *Candidatus Blochmannia* ANIm analysis

    Each cell represents a pairwise comparison between the named genomes on rows and columns, and the number in the cell is pairwise coverage of the alignment. The dendrograms are single-linkage clustering trees generated from the matrix of pairwise coverage results. The default colour scheme colours cells with identity > 0.50 as red, and those with < 0.50 as blue. This division corresponds to a strict majority of each genome in the comparison being alignable (i.e. homologous: a plausible minimum requirement for being considered "the same thing").

Several graphics output formats are available, including ``.png``, ``.pdf`` and ``.svg``.


.. _Anaconda: https://www.anaconda.com/what-is-anaconda/
.. _NCBI Taxonomy database: https://www.ncbi.nlm.nih.gov/taxonomy
.. _PyPI: https://pypi.python.org/pypi
.. _SQLite3: https://www.sqlite.org/index.html
