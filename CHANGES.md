# CHANGES.md

## v0.2

* Merged [pull request](https://github.com/widdowquinn/pyani/pull/17) from peterjc to make printing from tests Python3-friendly.
* Merged [pull request](https://github.com/widdowquinn/pyani/pull/21) from peterjc to use `open()` for opening files.
* Merged [pull request](https://github.com/widdowquinn/pyani/pull/24) from peterjc to cope with missing labels/classes more gracefully
* Fixed `-s`/`--fragsize` option in `average_nucleotide_identity.py` (thanks to Joseph Adelskov for hte report).
* BLAST and `nucmer` results are now written to a subdirectory of the output folder. By default, these sequence search output files are compressed, but this behaviour can be suppressed using the `--nocompress` option.
* Added `genbank_get_genomes_by_taxon.py` as an aid to downloading publicly-available genome files from GenBank, for analysis.
