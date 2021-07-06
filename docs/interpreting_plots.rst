=================================
Interpreting the Graphical Output
=================================
..
    Graphical output is obtained by executing the ``pyani plot`` subcommand, specifying the output directory and run ID. Optionally, output file formats and the graphics drawing method can be specified.

    .. code-block:: bash

    pyani plot --formats png,pdf --method seaborn C_blochmannia_ANIm 1

    Supported output methods are:

    - ``seaborn``
    - ``mpl`` (``matplotlib``)
    - ``plotly``

----------
The Output
----------

``pyani plot`` generates five plots corresponding to the matrices that ``pyani report`` produces:

  - percentage identity of aligned regions
  - percentage coverage of each genome by aligned regions
  - number of aligned bases on each genome
  - number of "similarity errors" on each genome
  - a Hadamard matrix of percentage identity multiplied by percentage coverage for each comparison

.. figure:: images/matrix_identity_1.png
    :alt: percentage identity matrix for *Candidatus Blochmannia* ANIm analysis

    Percentage identity matrix for *Candidatus Blochmannia* ANIm analysis

    Each cell represents a pairwise comparison between the named genomes on rows and columns, and the number in the cell is the pairwise identity *of aligned regions*. The dendrograms are single-linkage clustering trees generated from the matrix of pairwise identity results. The default colour scheme colours cells with identity > 0.95 as red, and those with < 0.95 as blue. This division corresponds to a widely-used convention for bacterial species boundaries.

.. figure:: images/matrix_coverage_1.png
    :alt: percentage coverage matrix for *Candidatus Blochmannia* ANIm analysis

    Percentage coverage matrix for *Candidatus Blochmannia* ANIm analysis

    Each cell represents a pairwise comparison between the named genomes on rows and columns, and the number in the cell is pairwise coverage of each genome by aligned regions in the comparison. The dendrograms are single-linkage clustering trees generated from the matrix of pairwise coverage results. The default colour scheme colours cells with identity > 0.50 as red, and those with < 0.50 as blue. This division corresponds to a strict majority of each genome in the comparison being alignable (a plausible minimum requirement for two sequences being considered "the same thing").

-----------------------
Understanding the plots
-----------------------

For each of the ANI methods provided by `pyani` the values shown in the plots are calculated differently — as a result of the differences in method/output of the third-party tool being called.

**Average nucleotide identity** is a measure of similarity between two genomes. Depending on the ANI method used, this may be symmetrical: comparing genome A to genome B is the same as comparing genome B to genome A; or asymmetrical: A --> B is different from B --> A. `pyani` provides both symmetrical and assymmetrical ANI methods:

  - ANIm — symmetrical
  - FastANI — asymmetrical (only available in version 0.3.0-alpha)
  - ANIb — asymmetrical (only available in version 0.2)
  - ANIblastall — asymmetrical (only available in version 0.2)
  - TETRA — (only available in version 0.2)

**Alignment coverage** is the proportion of the query genome that aligns against the reference genome (this can be asymmetrical: comparing genome A to genome B may give different coverage values for A and B).

  - in ANIm this is ``alignment_length / genome_length``
  - in fastANI this is ``matched_fragments / all_fragments``


**Alignment length** is the count of bases in the pairwise alignment of two genomes.

  - in ANIm this is calculated as ``reference_length + insertions - deletions``
  - in fastANI this is ``matched_fragments * fragment_length``

The **similarity errors** graph shows a measure of the number of bases/positions that do not match.

  - in ANIm this is ``non-identities + insertions + deletions``
  - in fastANI this is ``all_fragments - matched_fragents``

The **Hadamard** ouptut is the product (identity x coverage), as described at `Hadamard product`_ of identity and coverage. It's meant to provide a measure that allows you to interpret identity and coverage simultaneously.

  - this is always ``ANI * coverage``, but the plot is not symmetric coverage may differ for query and reference genomes

.. _Hadamard product: https://en.wikipedia.org/wiki/Hadamard_product_(matrices)
