.. _pyani-installation:

==================
Installation Guide
==================

We support four ways to install and run ``pyani`` on your system:

1. Installation with `Anaconda`_ (i.e. the ``conda`` package manager) **[Recommended]**
2. Installation *via* ``pip`` (i.e. from `PyPI`_)
3. Installation from source (i.e. download from `GitHub`_)
4. Installation of a `Docker`_ image

.. NOTE::
    If you wish to contribute to development of ``pyani``, you will require developer tools and packages that are not included in these installation instructions. To set up a local environment suitable for developing ``pyani``, please refer to the :ref:`pyani-contributing` page.

---------------------------------
1. Installation with ``Anaconda``
---------------------------------

`Anaconda <https://www.anaconda.com/>` is a ``Python`` language distribution that includes the ``conda`` package manager. Several *channels*, themed collections of packages, are available throught he ``conda`` package manager. The latest release of ``pyani`` should always be available from the `bioconda <https://anaconda.org/bioconda>`_ channel.

We recommend installation of ``pyani`` using the `Anaconda3 <https://www.anaconda.com/>`_ or `Miniconda3 <https://docs.conda.io/en/latest/miniconda.html>`_ distributions, and the ``conda`` package manager. This provides a straightforward way of managing third-party tool (e.g. ``MUMmer`` and ``NCBI-BLAST``) installation, and allows for installation of ``pyani`` in a virtual environment, protected from other installations on the system.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1.1 Install ``Anaconda`` or ``Miniconda``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If not already available on your system, install ``Anaconda3`` or ``Miniconda3`` on your system, following the instructions for your system on their respective websites:

- `Anaconda3 download <https://www.anaconda.com/distribution/>`_
- `Miniconda3 download <https://docs.conda.io/en/latest/miniconda.html>`_

After installation, you can check whether ``conda`` is available by issuing the following command in your terminal:

.. code-block:: bash

    $ conda -V
    conda 4.7.10

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1.2 (optional) Create a new ``conda`` environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Creation of a new ``conda`` environment is not necessary for installation of ``pyani``, but it can be sometimes useful to specify the available version of ``Python``, or to keep the versions of third party tools and Python packages separate from the system-level install. This may be particularly useful if, for example, the default system-level ``Python`` is ``Python2``, as ``Python3`` is required for ``pyani``.

The following commands will create and activate a new ``conda`` environment called ``pyani``, using ``Python 3.7``:

.. code-block:: bash

    conda create --name pyani --yes python=3.7  # --yes accepts all suggestions
    conda activate pyani

If successful, you should see the prefix ``(pyani)`` before your terminal prompt. To exit the current ``conda`` environment, issue the following command:

.. code-block:: bash

    conda deactivate

- `Managing conda environments <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1.3 Add required ``conda`` channels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``pyani`` and many of the packages it requires are provided in the ``conda`` *channel* called `bioconda <https://anaconda.org/bioconda>`_. *Channels* are locations where ``conda`` will look for packages, and they are typically grouped by topic area.

To install ``pyani``, you will need to make available the ``defaults``, ``conda-forge`` and ``bioconda`` channels, with the following commands:

.. code-block:: bash

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

- `Managing conda channels <https://docs.anaconda.com/anaconda/navigator/tutorials/manage-channels/>`_

^^^^^^^^^^^^^^^^^^^^^
1.4 Install ``pyani``
^^^^^^^^^^^^^^^^^^^^^

The ``bioconda`` distribution of ``pyani`` will install all necessary packages and software required to run the software, including ``NCBI-BLAST+`` and ``MUMmer``, with the following command:

.. code-block:: bash

    conda install --yes pyani

When installation is complete, you can check for the availability of the ``pyani`` programs with the following commands:

.. code-block::  bash

    $ pyani --version
    pyani 0.2.9
    $ average_nucleotide_identity.py --version
    average_nucleotide_identity.py: pyani 0.2.9
    $ genbank_get_genomes_by_taxon.py --version
    genbank_get_genomes_by_taxon.py: pyani 0.2.9

.. ATTENTION::
    If you wish to use the ``ANIblastall`` legacy ANIb method, then the legacy ``NCBI-BLAST`` tools need to be installed. These are not available through ``conda`` and must be installed manually for your system, as described below.


----------------------------
2. Installation with ``pip``
----------------------------

`PyPI is the Python Packaging Index <https://pypi.org/>`_, a repository of software for the ``Python`` language. Packages from ``PyPI`` can be installed using the ``pip`` package installer, which should come preinstalled with your system's ``Python3``. The latest release of ``pyani`` should always be available from ``PyPI``.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
2.1 (optional) Create and activate a new virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As with the ``conda`` installation route above, it can be useful to separate the version of ``Python`` you use for ``pyani``, and any installed packages and tools, from the system-level ``Python`` installation. There are multiple tools available to do this, and for convenience we list some below:

- `Anaconda distribution <https://anaconda.org/>`_
- `pipenv <https://docs.pipenv.org/en/latest/>`_
- `pyenv <https://github.com/pyenv/pyenv>`_
- `virtualenv <https://virtualenv.pypa.io/en/latest/>`_

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
2.2 Install third-party tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two third-party software tools are needed to perform ANIm and ANIb analysis:

- `MUMmer3 <http://mummer.sourceforge.net/>`_ for ANIm
- `NCBI-BLAST+ <https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_ for ANIb

These tools are not part of the ``PyPI`` distribution of ``pyani``, and should be installed according to the instructions on their respective websites.

.. TIP::
    If you are using a ``conda`` environment for ``pyani``, you can install both tools with a single command: ``conda install --yes blast mummer``.

.. ATTENTION::
    If you wish to use the ``ANIblastall`` legacy ANIb method, then the legacy ``NCBI-BLAST`` tools need to be installed. These are not available through ``conda`` and must be installed manually for your system, as described below.

^^^^^^^^^^^^^^^^^^^^^
2.3 Install ``pyani``
^^^^^^^^^^^^^^^^^^^^^

The ``pyani`` programs, and their ``Python`` dependencies, can be installed with the command:

.. code-block:: bash

    pip install pyani

When installation is complete, you can check for the availability of the ``pyani`` programs with the following commands:

.. code-block::  bash

    $ pyani --version
    pyani 0.2.9
    $ average_nucleotide_identity.py --version
    average_nucleotide_identity.py: pyani 0.2.9
    $ genbank_get_genomes_by_taxon.py --version
    genbank_get_genomes_by_taxon.py: pyani 0.2.9

---------------------------
3. Installation from source
---------------------------

The source code for the current ``pyani`` release can always be found on `GitHub`_:

- `Current pyani release source <https://github.com/widdowquinn/pyani/releases>`_

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
3.1 (optional) Create and activate a new virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As with the ``conda`` installation route above, it can be useful to separate the version of ``Python`` you use for ``pyani``, and any installed packages and tools, from the system-level ``Python`` installation. There are multiple tools available to do this, and for convenience we list some below:

- `Anaconda distribution <https://anaconda.org/>`_
- `pipenv <https://docs.pipenv.org/en/latest/>`_
- `pyenv <https://github.com/pyenv/pyenv>`_
- `virtualenv <https://virtualenv.pypa.io/en/latest/>`_

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
3.2 Install third-party tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two third-party software tools are needed to perform ANIm and ANIb analysis:

- `MUMmer3 <http://mummer.sourceforge.net/>`_ for ANIm
- `NCBI-BLAST+ <https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_ for ANIb

These tools are not part of the ``PyPI`` distribution of ``pyani``, and should be installed according to the instructions on their respective websites.

.. TIP::
    If you are using a ``conda`` environment for ``pyani``, you can install both tools with a single command: ``conda install --yes blast mummer``.

.. ATTENTION::
    If you wish to use the ``ANIblastall`` legacy ANIb method, then the legacy ``NCBI-BLAST`` tools need to be installed. These are not available through ``conda`` and must be installed manually for your system, as described below.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
3.3 Download and extract the ``pyani`` source code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Click on one of the two links on the ``pyani`` `releases page <https://github.com/widdowquinn/pyani/releases>`_, or use a download tool such as ``curl`` or ``wget`` to download the link to a convenient location. For example:

.. code-block:: bash

    wget https://github.com/widdowquinn/pyani/archive/v0.2.9.tar.gz

Then extract the source code archive. This will create a new directory called ``pyani-<CURRENT_VERSION>`` in your current location:

.. code-block:: bash

    $ tar -zxf v0.2.9.tar.gz
    $ ls
    [...]
    pyani-0.2.9/
    [...]

^^^^^^^^^^^^^^^^^^^^^
3.4 Install ``pyani``
^^^^^^^^^^^^^^^^^^^^^

Change directory in the terminal to the source code directory, e.g.:

.. code-block:: bash

    $ cd pyani-0.2.9
    $ ls
    CHANGES.md                               LICENSE                                  _config.yml                              requirements.txt
    CITATIONS                                MANIFEST.in                              bin/                                     setup.cfg
    CONTRIBUTORS.md                          Makefile                                 pyani/                                   setup.py
    Dockerfile-average_nucleotide_identity   README-docker.md                         requirements-dev.txt                     test-requirements.txt
    Dockerfile-genbank_get_genomes_by_taxon  README.md                                requirements-pip.txt                     tests/

``pyani`` can be installed using ``Python3``'s setup tools, using the command:

.. code-block:: bash

    python setup.py install

This will download and install the ``Python`` packages that ``pyani`` needs to run. When installation is complete, you can check for the availability of the ``pyani`` programs with the following commands:

.. code-block::  bash

    $ pyani --version
    pyani 0.2.9
    $ average_nucleotide_identity.py --version
    average_nucleotide_identity.py: pyani 0.2.9
    $ genbank_get_genomes_by_taxon.py --version
    genbank_get_genomes_by_taxon.py: pyani 0.2.9


-------------------------------
4. Installation with ``Docker``
-------------------------------

``Docker`` is a platform for running applications that uses *containerisation* to share virtual machines that come pre-installed with all required tools and dependencies. The latest ``pyani`` programs should always be available as separate ``Docker`` containers from ``DockerHub``:

- `average_nucleotide_identity.py <https://cloud.docker.com/repository/docker/leightonpritchard/average_nucleotide_identity>`_
- `genbank_get_genomes_by_taxon.py <https://cloud.docker.com/repository/docker/leightonpritchard/genbank_get_genomes_by_taxon>`_

In order to use the containerised versions of ``pyani``, you must have ``Docker`` installed and working on your system. To do so, please follow the instructions at the ``Docker`` website:

- `Docker website <https://www.docker.com>`_

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
4.1 Running ``pyani`` with ``Docker``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To pull (if necessary) and run the ``pyani`` programs in a ``Docker`` image on your local system, use the following commands (with the ``Docker`` daemon running):

.. code-block:: bash

    docker run -v ${PWD}:/host_dir leightonpritchard/average_nucleotide_identity:v<REQUIRED_VERSION>
    docker run -v ${PWD}:/host_dir leightonpritchard/genbank_get_genomes_by_taxon:v<REQUIRED_VERSION>
    docker run -v ${PWD}:/host_dir leightonpritchard/pyani:v<REQUIRED_VERSION>

.. TIP::
    If no tag is specified, then ``Docker`` will attempt to use the ``:latest`` tag, which may not exist.

.. NOTE::
    The ``-v ${PWD}:/host_dir`` links the ``Docker`` image to the current working directory, and enables the ``pyani`` programs to see files below the current execution point in your filesystem. This is necessary for the analysis to proceed.


.. _Anaconda: https://anaconda.org/bioconda/pyani
.. _Docker: https://cloud.docker.com/repository/docker/leightonpritchard/average_nucleotide_identity
.. _GitHub: https://github.com/widdowquinn/pyani
.. _PyPI: https://pypi.org/project/pyani/

