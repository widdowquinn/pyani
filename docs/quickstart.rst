.. _pyani-quickstart:

================
QuickStart Guide
================

------------
Installation
------------

To use ``pyani`` you will need to install it. Installation is possible *via* two popular Python package managers:

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
    ``pyani`` can be installed from source. Installation instructions can be found on the :ref:`pyani-installation` page.


-----------------
``pyani`` Example
-----------------

An example ANIm analysis using ``pyani`` is provided as a walkthrough in the stages below:

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

``pyani`` can be used with an existing local set of genomes. For this walkthrough a new set of genomes will be obtained from GenBank, using the ``pyani download`` command.

.. TIP::
    ``pyani`` requires an *indexed* set of genomes, and the visualisation and classification steps benefit from having ``classes.txt`` and ``labels.txt`` files. These are generated automatically when downloading genomes, but you must create them in other ways when applying ``pyani`` to a set of local files.

.. ATTENTION::
    To use their online resources programmatically, NCBI require that you provide your email address for contact purposes if jobs go wrong, and for their own usage statistics. This is specified with the ``--email <EMAIL ADDRESS>`` argument of ``pyani download``.

Use the pyani.py download subcommand to download all available genomes for *Candidatus Blochmannia* from NCBI. The taxon ID for this grouping is ``203804``.

* `NCBI Taxonomy database`_

.. code-block:: bash

    pyani download C_blochmannia --email my.email@my.domain -t 203804

This produces a new ``directory`` (``C_blochmannia``) with the following contents:

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


^^^^^^^^^^^^^^^^^^
2. Create database
^^^^^^^^^^^^^^^^^^

``pyani`` uses a database to store genome data and analysis results. Create a new, clean, database with the command:

.. code-block:: bash

    pyani createdb

.. TIP::
    This creates a new database in the default location (``.pyani/pyanidb``), but the name and location of this database can be controlled with the ``pyani createdb`` command, and specified in each of the subsequent commands.

^^^^^^^^^^^^^^^^^^^^^^^^
3. Conduct ANIm analysis
^^^^^^^^^^^^^^^^^^^^^^^^

Run ANIm on the downloaded genomes, using the command:

.. code-block:: bash

    pyani anim C_blochmannia C_blochmannia_ANIm \
        --name "C. blochmannia run 1" \
        --labels C_blochmannia/labels.txt --classes C_blochmannia/classes.txt

This will run an ANIm analysis on the genomes in the ``C_blochmannia`` directory. The analysis results will be stored in the database you created, identified with the name ``C. blochmannia run 1``, but the comparison files will be written to the ``C_blochmannia_ANIm`` directory.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
4. Reporting Analyses and Analysis Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Show all the runs contained in the (default) database with the command:

.. code-block:: bash

    pyani report --runs C_blochmannia_ANIm/ --formats html,excel,stdout

This will report the relevant information to new files the ``pyani report`` command creates in the ``C_blocahmannia_ANIm`` directory.

.. code-block:: bash

    $ tree -L 1 C_blochmannia_ANIm/
    C_blochmannia_ANIm/
    ├── nucmer_output
    ├── runs.html
    ├── runs.tab
    └── runs.xlsx


.. _NCBI Taxonomy database: https://www.ncbi.nlm.nih.gov/taxonomy
.. _PyPI: https://pypi.python.org/pypi