# CHANGES.md

## v0.2.7
* Fix issue #97 where valid input arguments were not recognised in the download script

## v0.2.6
* Add Dockerfiles for making Docker images

## v0.2.5
* Version bump to fix `PyPI` wheel/script shebangs

## v0.2.4

* `ANIm` now uses `delta-filter` to remove alignments of repeat regions (issue #91)
* added `--filter_exe` option to specify location of `delta-filter` utility (issue #91)
* fixed `--format` option so that GenBank downloads work again (issue #89)
* add `--SGEargs` option to `average_nucleotide_identity.py` for custom qsub settings
* `README.md` badges now clickable
* `--version` switch added to `average_nucleotide_identity.py`
* FTP timeouts are now caught differently in `genbank_get_genomes_by_taxon.py`
* Additional characters in NCBI FTP URIs now escaped in `genbank_get_genomes_by_taxon.py` - should be fewer failed downloads
* Modified error messaging when `NUCmer` alignment fails
* `average_nucleotide_identity.py` argument documentation improvements
* Script now fails immediately if label or class files missing (issue #78)
* Changes to `--noclobber` log behaviour (issue #79)
* fixed `--rerender` code (issue #85)


## v0.2.3
* fixes a bug in the installed scripts where the shebang (`#!`) in wheel and egg packages pointed to a development Python

## v0.2.2

* fix for issue #53 (--maxmatch has no effect)
* fix to `genbank_get_genomes_by_taxon.py` to account for NCBI FTP location changes
* fixed issue #52 (local variable bug)
* fixed issued #49 (TETRA failure) and #51 (matplotlib bug)
* add several tests and support for `codecov.io`, `landscape.io` and `Travis-CI`
* removed requirement for `rpy2`
* moved scripts to `bin/` subdirectory


## v0.2.1

* `pyani` now requires `rpy2` v2.8.0 in order to satisfy running under Anaconda (see issue #26)
* `pyani` now checks for presence of `rpy2` and - when run from source - if `rpy2` is not available, `pyani` doesn't throw an error until R graphical output is requested. If installed *via* `pip`, then `pyani` still raises `pkg_resources.DistributionNotFound` if `rpy2` is missing.
* Updated `genbank_get_genomes_by_taxon.py` script to use the new FTP locations at NCBI for each assembly.
* Fixed bug where `ANIb` would not go to completion if empty BLASTN files were generated (see issue #27)
* Fixed bug where `ANIm` would not finish under `multiprocessing` if input sequences were highly divergent.
* Added Hadamard product of percentage identity and alignment coverage as output.
* Fixed bug where label/classes are out of sync with new NCBI downloaded filenames
* Added --rerender option to draw (new) graphics from old output, without recalculation
* Corrected matplotlib row dendrogram orientation
* Seaborn output no longer dumps core on large (ca. 500 genome) datasets
* `genbank_get_genomes_by_taxon.py` attempts to identify cause for failed downloads and correct, where nomenclature/versions are at fault
* graceful replacement of classes that are not present in `classes.txt`
* add `pyani` version to log file

## v0.2.0

* Merged [pull request](https://github.com/widdowquinn/pyani/pull/17) from peterjc to make printing from tests Python3-friendly.
* Merged [pull request](https://github.com/widdowquinn/pyani/pull/21) from peterjc to use `open()` for opening files.
* Merged [pull request](https://github.com/widdowquinn/pyani/pull/24) from peterjc to cope with missing labels/classes more gracefully
* Fixed `-s`/`--fragsize` option in `average_nucleotide_identity.py` (thanks to Joseph Adelskov for hte report).
* BLAST and `nucmer` results are now written to a subdirectory of the output folder. By default, these sequence search output files are compressed, but this behaviour can be suppressed using the `--nocompress` option.
* Added `genbank_get_genomes_by_taxon.py` as an aid to downloading publicly-available genome files from GenBank, for analysis.
