.. _pyani-quickstart:

================
QuickStart Guide
================

------------
Installation
------------

^^^^^^^^^^^^^^^
Using ``conda``
^^^^^^^^^^^^^^^

``pyani`` is available through the ``bioconda`` channel of Anaconda:

.. code-block:: bash

    conda install -c bioconda pyani 


^^^^^^^^^^^^^^
Using ``PyPI``
^^^^^^^^^^^^^^

``pyani`` is available *via* the `PyPI`_ package manager for Python:

.. code-block:: bash

    pip install pyani 


^^^^^^^^^^^
From source
^^^^^^^^^^^

At the command-line, use ``git`` to clone the current version of the ``pyani`` repository:

.. code-block:: bash

    git clone git@github.com:widdowquinn/pyani.git

Change to the newly-created ``pyani`` subdirectory:

.. code-block:: bash

    cd pyani

Install the package and program, using the ``setup.py`` file:

.. code-block:: bash

    python setup.py install

(other installation methods can be found on the :ref:`pyani-installation` page)

-----------------
``pyani`` Example
-----------------

To see options available for the ``pyani`` program, use the ``-h``
(help) option:

.. code-block:: bash

    pyani -h

An example ``pyani`` analysis is provided as a walkthrough in the stages below:

1. Download genomes
2. Create a database to hold genome data and analysis results
3. Perform ANI analysis
4. Report and visualise analysis results
5. Generate species hypotheses (classify genomes) using the analysis results

^^^^^^^^^^^^^^^^^^^
1. Download genomes
^^^^^^^^^^^^^^^^^^^

``pyani`` can be used with an existing local set of genomes, but for this walkthrough we will download a new set of genomes with the ``pyani download`` command.

.. TIP::
    ``pyani`` requires an *indexed* set of genomes, and the visualisation and classification steps benefit from having ``classes.txt`` and ``labels.txt`` files. These are generated automatically when downloading genomes, but you must create these in other ways when applying ``pyani`` to a set of local files.

.. ATTENTION::
    To use their online resources programmatically, NCBI require that you provide your email address for contact purposes if jobs go wrong, and for their own usage statistics. This is specified with the ``--email <EMAIL ADDRESS>`` argument of ``pyani download``.

Use the pyani.py download subcommand to download all available genomes for Candidatus Blochmannia from NCBI. The taxon ID for this grouping is ``203804``.

* `NCBI Taxonomy database`_

.. code-block:: bash

    pyani.py download C_blochmannia --email my.email@my.domain -t 203804

This produces a new subdirectory (C_blochmannia) with the following contents:

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


.. _NCBI Taxonomy database: https://www.ncbi.nlm.nih.gov/taxonomy
.. _PyPI: https://pypi.python.org/pypi