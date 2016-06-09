# CHANGES.md

## v0.2.0.dev

* `pyani` now requires `rpy2` v2.8.0 in order to satisfy running under Anaconda (see issue #26)
* `pyani` now checks for presence of `rpy2` and - when run from source - if `rpy2` is not available, `pyani` doesn't throw an error until R graphical output is requested. If installed *via* `pip`, then `pyani` still raises `pkg_resources.DistributionNotFound` if `rpy2` is missing.
* Updated `genbank_get_genomes_by_taxon.py` script to use the new FTP locations at NCBI for each assembly.
* Fixed bug where `ANIb` would not go to completion if empty BLASTN files were generated
* Fixed bug where `ANIm` would not finish under `multiprocessing` if input sequences were highly divergent.


## v0.2.0

* Merged [pull request](https://github.com/widdowquinn/pyani/pull/17) from peterjc to make printing from tests Python3-friendly.
* Merged [pull request](https://github.com/widdowquinn/pyani/pull/21) from peterjc to use `open()` for opening files.
* Merged [pull request](https://github.com/widdowquinn/pyani/pull/24) from peterjc to cope with missing labels/classes more gracefully
* Fixed `-s`/`--fragsize` option in `average_nucleotide_identity.py` (thanks to Joseph Adelskov for hte report).
* BLAST and `nucmer` results are now written to a subdirectory of the output folder. By default, these sequence search output files are compressed, but this behaviour can be suppressed using the `--nocompress` option.
* Added `genbank_get_genomes_by_taxon.py` as an aid to downloading publicly-available genome files from GenBank, for analysis.
