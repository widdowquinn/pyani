.. _pyani-requirements:

============
Requirements
============

""""""
Python
""""""

``pyani`` is written in, and runs on, `Python 3`_. This is an essential requirement to run ``pyani``. We use several features specific to Python 3, and ``pyani`` will not work with Python 2.

We recommend the `Anaconda`_ or `Minicodan`_ distribution for Python, as this provides the ``conda`` package manager and also allows for installation of ``pyani`` and its dependencies using the `Bioconda`_ channel.

""""""
BLAST+
""""""

``BLAST+`` is the NCBI biological sequence search program, and is required for the ``pyani anib`` command (average nucleotide identity with ``BLAST+``).

""""""""""""
Legacy BLAST
""""""""""""

Legacy ``BLAST`` is the deprecated  version of NCBI biological sequence search program. It is no longer maintained, but is required for the ``pyani aniblastall`` command (average nucleotide identity with ``BLAST``).


""""""
MUMmer
""""""

``MUMmer`` is a biological sequence alignment program, and is required for the ``pyani anim`` subcommand (average nucleotide identity with ``MUMmer``)



.. _Anaconda: https://www.anaconda.com/distribution/
.. _Bioconda: https://anaconda.org/bioconda
.. _Python 3: https://www.python.org/