# README.md (pyani)

[![pyani PyPi version](https://img.shields.io/pypi/v/pyani.svg "PyPi version")](https://pypi.python.org/pypi/pyani)
[![pyani licence](https://img.shields.io/pypi/l/pyani.svg "PyPi licence")](https://github.com/widdowquinn/pyani/blob/master/LICENSE)
[![pyani TravisCI build status](https://api.travis-ci.org/widdowquinn/pyani.svg?branch=master)](https://travis-ci.org/widdowquinn/pyani/branches)
[![pyani codecov.io coverage](https://img.shields.io/codecov/c/github/widdowquinn/pyani/master.svg)](https://codecov.io/github/widdowquinn/pyani)
[![Code Health](https://landscape.io/github/widdowquinn/pyani/master/landscape.svg?style=flat)](https://landscape.io/github/widdowquinn/pyani/master) 

## Overview
`pyani` is a Python3 module and script that provides support for calculating average nucleotide identity (ANI) and related measures for whole genome comparisons, and rendering relevant graphical summary output. Where available, it takes advantage of multicore systems, and can integrate with [SGE/OGE](http://gridscheduler.sourceforge.net/)-type job schedulers for the sequence comparisons.

`pyani` also installs the script `pyani.py`, which enables command-line based analysis of genomes.

## Installation

The easiest way to install `pyani` is to use `pip3`:

```
pip3 install pyani
```

From version 0.1.3.2 onwards, this should also install all the required Python package dependencies. Prior to this version (i.e. 0.1.3.1 and earlier), you can acquire these dependencies with `pip -r`, and pointing at `requirements.txt` from this repository:

```
pip3 install -r requirements.txt
```

## Walkthrough: A First Analysis

The command-line interface to `pyani` uses subcommands. These separate individual steps of an analysis into separate actions.

1. Download genomes
2. Create a database to hold genome data and analysis results
3. Perform ANI analysis
4. Report and visualise analysis results
5. Generate species hypotheses (classify genomes) using the analysis results

The steps are described in detail with examples, below.

### 1. Download genome data

The first step is to obtain genome data for analysis. `pyani` expects to find each individual genome in its own FASTA file (that file can contain multiple sequences - chromosomes and plasmids; sequenced scaffolds, etc). All the FASTA files for an analysis are expected to be located in a single subdirectory (with optional `labels` and `classes` files). You can arrange your data manually, but `pyani` provides a subcommand that downloads all genomes in a taxon subtree from NCBI, and organises them ready for use with `pyani`.

We'll use the `pyani.py download` subcommand to download all available genomes for *Candidatus Blochmannia* from NCBI. The taxon ID for this grouping is [203804](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=203804&lvl=3&lin=f&keep=1&srchmode=1&unlock).

```bash
pyani.py download C_blochmannia --email my.email@my.domain -t 203804 -v -l C_blochmannia_dl.log 
```

The first argument is the output directory into which the downloaded genomes will be written (`C_blochmannia`). To download anything from NCBI we must provide an email address (`--email my.email@my.domain`), and to specify which taxon subtree we want to download we provide the taxon ID (`-t 203804`).

Here we have also requested verbose output (`-v`), and write a log file for reproducible research/diagnosing bugs and errors (`-l C_blochmannia_dl.log`). 

This produces a new subdirectory (`C_blochmannia`) with the following contents:

```bash
$ tree C_blochmannia
C_blochmannia
├── GCF_000011745.1_ASM1174v1_genomic.fna
├── GCF_000011745.1_ASM1174v1_genomic.fna.gz
├── GCF_000011745.1_ASM1174v1_genomic.md5
├── GCF_000011745.1_ASM1174v1_hashes.txt
├── GCF_000043285.1_ASM4328v1_genomic.fna
├── GCF_000043285.1_ASM4328v1_genomic.fna.gz
├── GCF_000043285.1_ASM4328v1_genomic.md5
├── GCF_000043285.1_ASM4328v1_hashes.txt
├── GCF_000185985.2_ASM18598v2_genomic.fna
├── GCF_000185985.2_ASM18598v2_genomic.fna.gz
├── GCF_000185985.2_ASM18598v2_genomic.md5
├── GCF_000185985.2_ASM18598v2_hashes.txt
├── GCF_000331065.1_ASM33106v1_genomic.fna
├── GCF_000331065.1_ASM33106v1_genomic.fna.gz
├── GCF_000331065.1_ASM33106v1_genomic.md5
├── GCF_000331065.1_ASM33106v1_hashes.txt
├── GCF_000973505.1_ASM97350v1_genomic.fna
├── GCF_000973505.1_ASM97350v1_genomic.fna.gz
├── GCF_000973505.1_ASM97350v1_genomic.md5
├── GCF_000973505.1_ASM97350v1_hashes.txt
├── GCF_000973545.1_ASM97354v1_genomic.fna
├── GCF_000973545.1_ASM97354v1_genomic.fna.gz
├── GCF_000973545.1_ASM97354v1_genomic.md5
├── GCF_000973545.1_ASM97354v1_hashes.txt
├── classes.txt
└── labels.txt
```

Seven genomes have been downloaded, and each is represented by four files:

- `_genomic.fna.gz`: the compressed genome sequence
- `_genomic.fna`: the uncompressed genome sequence
- `_genomic.md5`: an MD5 hash/checksum of the (uncompressed) genome sequence; this was generated during the download
- `_hashes.txt`: a list of MD5 hashes; this is provided by NCBI and is a reference to be sure that the download did not corrupt the genome sequence

There are two additional plain text files: `classes.txt` and `labels.txt`, which provide alternative labels for use in the analysis. These files are generated during the download.

### 2. Create an analysis database

`pyani` uses a database to store genome data and analysis results. This is convenient for data sharing and developing custom analyses, but also makes it easier to extend an existing ANI analysis with new genomes, without having to repeat calculations that were already performed. 

To create a new, clean, database in the default location (`.pyani/pyanidb`) issue the command:

```bash
pyani.py createdb -v -l C_blochmannia_createdb.log
```

As above, the verbose (`-v`) and log file (`-l C_blochmannia_createdb.log`) options allow for reproducible work. The default database location is in the hidden directory (`.pyani`):

```bash
$ tree .pyani
.pyani
└── pyanidb
```

Subsequent `pyani` commands will assume this location for the database, but you can specify the location when creating a database, or using an existing database.

### 3. Conduct ANI analysis

`pyani` provides four subcommands to run ANI analyses:

- `anim`: ANIm
- `anib`: ANIb, using BLAST+
- `aniblastall`: ANIb, using legacy BLAST
- `tetra`: TETRA

In this walkthrough, we'll run ANIm on the downloaded genomes, using the command:

```bash
pyani.py anim C_blochmannia C_blochmannia_ANIm -v -l C_blochmannia_ANIm.log
```

All four analysis commands operate in a similar way. The first two arguments are paths to directories: the first path is to a directory containing input genomes, and the second is the path to an output directory for storing intermediate results. The `-v` and `-l` arguments work as above, specifying verbose output and logging output to a file. 

You will probably notice that the verbose output is very verbose, to enable informative identification of any problems. In particular, the verbose output (which is also written to the log file) writes out the command-lines used for the pairwise comparisons so, if something goes wrong, you can test whether a specific comparison can be run at the command-line *at all*, to aid diagnosis of any problems.

#### Rerunning the same analysis

One reason for using a database backend for analysis results is so that, for very large analyses, we do not ever need to recalculate a pairwise genome comparison. All the analysis subcommands check whether input genomes have been used before (using the unique MD5 hash for each genome to identify whether it's been used previously), and whether the comparison of two genomes has been run, with the particular analysis settings that were used. If either genome was not seen before, or if the analysis settings are different, the comparison is performed.

You can test this for yourself by running the analysis command again, as below. You will see a number of messages indicating that genomes have been seen before, and that analyses performed before were skipped:

```bash
$ pyani.py anim C_blochmannia C_blochmannia_ANIm -v -l C_blochmannia_ANIm.log
INFO: command-line: pyani.py anim C_blochmannia C_blochmannia_ANIm -v -l C_blochmannia_ANIm.log
INFO: Running ANIm analysis
INFO: Adding analysis information to database .pyani/pyanidb
INFO: Current analysis has ID 2 in this database
INFO: Identifying input genome/hash files:
[…]
INFO: Adding genome data to database...
WARNING: Genome already in database with this hash and path!
WARNING: Using existing genome from database, row 1
[…]
INFO: Complete pairwise comparison list:
	[(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 3), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)]
INFO: Excluding pre-calculated comparisons
INFO: Comparisons still to be performed:
	[]
INFO: All comparison results already present in database (skipping comparisons)
INFO: Completed. Time taken: 0.211
```

### 4. Reporting Analyses and Analysis Results


### 5. Classifying Genomes from Analysis Results



## Testing `pyani`

`pyani` includes tests that can be run with `nosetest` (including coverage measurement using `coverage.py`) with the following command, executed from the repository root directory:

```
nosetests --cover-erase --cover-package=pyani --cover-html --with-coverage
```

Coverage output will be placed (by default) in the `cover` subdirectory, and can be loaded into the web browser.

## Running `pyani`

### Script: <a name="average_nucleotide_identity.py">`average_nucleotide_identity.py`</a>

The `average_nucleotide_identity.py` script - installed as part of this package - enables straightforward ANI analysis at the command-line, and uses the `pyani` module behind the scenes.

You can get a summary of available command-line options with `average_nucleotide_identity.py -h`

```
$ ./average_nucleotide_identity.py -h
usage: average_nucleotide_identity.py [-h] [-o OUTDIRNAME] [-i INDIRNAME] [-v]
                                      [-f] [-s FRAGSIZE] [-l LOGFILE]
                                      [--skip_nucmer] [--skip_blastn]
                                      [--noclobber] [--nocompress] [-g]
                                      [--gformat GFORMAT] [--gmethod GMETHOD]
                                      [--labels LABELS] [--classes CLASSES]
                                      [-m METHOD] [--scheduler SCHEDULER]
                                      [--workers WORKERS]
                                      [--SGEgroupsize SGEGROUPSIZE]
                                      [--maxmatch] [--nucmer_exe NUCMER_EXE]
                                      [--blastn_exe BLASTN_EXE]
                                      [--makeblastdb_exe MAKEBLASTDB_EXE]
                                      [--blastall_exe BLASTALL_EXE]
                                      [--formatdb_exe FORMATDB_EXE]
                                      [--write_excel] [--subsample SUBSAMPLE]
                                      [--seed SEED] [--jobprefix JOBPREFIX]


[…]
```

Example data and output can be found in the directory `test_ani_data`. The data are chromosomes of four isolates of *Caulobacter*. Basic analyses can be performed with the command lines:

```
$ ./average_nucleotide_identity.py -i tests/test_ani_data/ -o tests/test_ANIm_output -m ANIm -g
$ ./average_nucleotide_identity.py -i tests/test_ani_data/ -o tests/test_ANIb_output -m ANIb -g
$ ./average_nucleotide_identity.py -i tests/test_ani_data/ -o tests/test_ANIblastall_output -m ANIblastall -g
$ ./average_nucleotide_identity.py -i tests/test_ani_data/ -o tests/test_TETRA_output -m TETRA -g
```

The graphical output below, supporting assignment of `NC_002696` and `NC_011916` to the same species (*C.crescentus*), and the other two isolates to distinct species (`NC_014100`:*C.segnis*; `NC_010338`:*C.* sp K31), was generated with the command-line:

```
./average_nucleotide_identity.py -v -i tests/test_ani_data/ \
    -o tests/test_ANIm_output/ -g --gformat png,pdf,eps \
    --classes tests/test_ani_data/classes.tab \
    --labels tests/test_ani_data/labels.tab
```


![ANIm percentage identity for *Caulobacter* test data](tests/test_ani_data/ANIm_percentage_identity.png "ANIm percentage identity")
![ANIm alignment coverage for *Caulobacter* test data](tests/test_ani_data/ANIm_alignment_coverage.png "ANIm alignment coverage")
![ANIm alignment length for *Caulobacter* test data](tests/test_ani_data/ANIm_alignment_lengths.png "ANIm alignment length")
![ANIm alignment similarity errors for *Caulobacter* test data](tests/test_ani_data/ANIm_similarity_errors.png "ANIm alignment similarity")

### Script: <a name="genbank_get_genomes_by_taxon.py">`genbank_get_genomes_by_taxon.py`</a>

The script `genbank_get_genomes_by_taxon.py`, installed by this package, enables download of genomes from NCBI, specified by taxon ID. The script will download all available assemblies for taxa at or below the specified node in the NCBI taxonomy tree.

Command-line options can be viewed using:

```
$ genbank_get_genomes_by_taxon.py -h
usage: genbacnk_get_genomes_by_taxon.py [-h] [-o OUTDIRNAME] [-t TAXON] [-v]
                                        [-f] [--noclobber] [-l LOGFILE]
                                        [--format FORMAT] [--email EMAIL]
                                        [--retries RETRIES]
                                        [--batchsize BATCHSIZE]
[…]
```

For example, the NCBI taxonomy ID for *Caulobacter* is 75, so all publicly-available *Caulobacter* sequences can be obtained using the command-line:

```
$ genbank_get_genomes_by_taxon.py -o Caulobacter_downloads -v -t 75 -l Caulobacter_downloads.log --email me@my.email.domain
INFO: genbank_get_genomes_by_taxon.py: Mon Apr 18 17:22:54 2016
INFO: command-line: /Users/lpritc/Virtualenvs/pyani3/bin/genbank_get_genomes_by_taxon.py -o Caulobacter_downloads -v -t 75 -l Caulobacter_downloads.log --email me@my.email.domain
INFO: Namespace(batchsize=10000, email='me@my.email.domain', force=False, format='gbk,fasta', logfile='Caulobacter_downloads.log', noclobber=False, outdirname='Caulobacter_downloads', retries=20, taxon='75', verbose=True)
INFO: Set NCBI contact email to me@my.email.domain
INFO: Creating directory Caulobacter_downloads
INFO: Output directory: Caulobacter_downloads
INFO: Passed taxon IDs: 75
INFO: Entrez ESearch with query: txid75[Organism:exp]
INFO: Entrez ESearch returns 29 assembly IDs
INFO: Identified 29 unique assemblies
INFO: Taxon 75: 29 assemblies
[…]
INFO: Assembly 639581: 271 contigs
INFO: Assembly 233261: 17 contigs
INFO: Assembly 575291: 48 contigs
INFO: Mon Apr 18 17:25:46 2016
INFO: Done.
```

**NOTE:** You must provide a valid email to identify yourself to NCBI for troubleshooting.

The number of attempted retries for each download, and the size of a batch download can be modified. By default, the script will attempt 20 download retries, and obtain sequences in batches of 10000.


## DEPENDENCIES

Note that Python package dependencies should automatically be installed if you are using version 0.1.3.2 or greater, and installing with `pip install pyani`.

For earlier versions, you can satisfy dependencies by using `pip install -r requirements.txt` (using the `requirements.txt` file in this repository).

### For ANI analysis

* **Biopython** <http://www.biopython.org>
* **NumPy** <http://www.numpy.org/>
* **pandas** <http://pandas.pydata.org/>
* **SciPy** <http://www.scipy.org/>

#### Alignment tools

* **BLAST+** executable in the `$PATH`, or available on the command line (required for **ANIb** analysis) <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>
* **legacy BLAST** executable in the `$PATH` or available on the command line (required for **ANIblastall** analysis) <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/>
* **MUMmer** executables in the $PATH, or available on the command line (required for **ANIm** analysis) <http://mummer.sourceforge.net/>

### For graphical output

* **matplotlib** <http://matplotlib.org/>
* **seaborn** <https://github.com/mwaskom/seaborn>

## Method and Output Description

### Average Nucleotide Identity (ANI)

This module calculates Average Nucleotide Identity (ANI) according to one of a number of alternative methods described in, e.g.

* Richter M, Rossello-Mora R (2009) Shifting the genomic gold standard for the prokaryotic species definition. Proc Natl Acad Sci USA 106: 19126-19131. doi:10.1073/pnas.0906412106. (ANI1020, ANIm, ANIb)
* Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, et al. (2007) DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. Int J Syst Evol Micr 57: 81-91. doi:10.1099/ijs.0.64483-0.

ANI is proposed to be the appropriate *in silico* substitute for DNA-DNA 
hybridisation (DDH), and so useful for delineating species boundaries. A 
typical percentage threshold for species boundary in the literature is 95% 
ANI (e.g. Richter et al. 2009).

All ANI methods follow the basic algorithm:

- Align the genome of organism 1 against that of organism 2, and identify the matching regions
- Calculate the percentage nucleotide identity of the matching regions, as an average for all matching regions

Methods differ on: (1) what alignment algorithm is used, and the choice of parameters (this affects the aligned region boundaries); (2) what the input is for alignment (typically either fragments of fixed size, or the most complete assembly available).

* **ANIm**: uses MUMmer (NUCmer) to align the input sequences.
* **ANIb**: uses BLASTN+ to align 1020nt fragments of the input sequences
* **ANIblastall**: uses legacy BLASTN to align 1020nt fragments of the input sequences
* **TETRA**: calculates tetranucleotide frequencies of each input sequence

The algorithms takes as input correctly-formatted FASTA multiple sequence files. All sequences for a single organism should be contained in only one sequence file. Although it is possible to provide new labels for each input genome, for rendering graphical output, the names of these files are used for identification so it is best to name 
them sensibly.

Output is written to a named directory. The output files differ depending on the chosen ANI method.

* **ANIm**: MUMmer/NUCmer .delta files, describing each pairwise sequence alignment. Output as tab-separated plain text format tables describing: alignment coverage; total alignment lengths; similarity errors; and percentage identity (ANIm).
* **ANIb** and **ANIblastall**: FASTA sequences describing 1020nt fragments of each input sequence; BLAST nucleotide databases - one for each set of fragments; and BLASTN output files (tab-separated tabular format plain text) - one for each pairwise comparison of input sequences. Output as tab-separated plain text tables describing: alignment coverage; total alignment lengths; similarity errors; and percentage identity (ANIb or ANIblastall).
* **TETRA**: Tab-separated plain text files describing the Pearson correlations between Z-score distributions for each tetranucleotide in each input sequence (TETRA).

If graphical output is chosen, the output directory will also contain PDF, PNG and EPS files representing the various output measures as a heatmap with row and column dendrograms. Other output formats (e.g. SVG) can be specified with the `--gformat` argument.



## Licensing

Unless otherwise indicated, all code is subject to the following agreement:

    (c) The James Hutton Institute 2014, 2015
    Author: Leighton Pritchard

    Contact: leighton.pritchard@hutton.ac.uk

    Address: 
    Leighton Pritchard,
    Information and Computational Sciences,
    James Hutton Institute,
    Errol Road,
    Invergowrie,
    Dundee,
    DD6 9LH,
    Scotland,
    UK

The MIT License

Copyright (c) 2014-2015 The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
