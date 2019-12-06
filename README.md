# README.md (pyani)

[![pyani PyPi version](https://img.shields.io/pypi/v/pyani.svg "PyPi version")](https://pypi.python.org/pypi/pyani)
[![pyani licence](https://img.shields.io/pypi/l/pyani.svg "PyPi licence")](https://github.com/widdowquinn/pyani/blob/master/LICENSE)
[![CircleCI](https://circleci.com/gh/widdowquinn/pyani/tree/master.svg?style=shield)](https://circleci.com/gh/widdowquinn/pyani/tree/master)
[![pyani codecov.io coverage](https://img.shields.io/codecov/c/github/widdowquinn/pyani/development.svg)](https://codecov.io/github/widdowquinn/pyani)

[![GitHub Issues](https://img.shields.io/github/issues-closed/widdowquinn/pyani.svg)](https://github.com/widdowquinn/pyani/issues)
[![GitHub Stars](https://img.shields.io/github/stars/widdowquinn/pyani.svg)](https://github.com/widdowquinn/pyani/stargazers)

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/f3e56b2bf118471aabf09514a3e6af75)](https://www.codacy.com/app/widdowquinn/pyani?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=widdowquinn/pyani&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/widdowquinn/pyani/badge)](https://www.codefactor.io/repository/github/widdowquinn/pyani)
[![pyani sourcerank](https://img.shields.io/librariesio/sourcerank/pypi/pyani.svg?logo=koding&logoColor=white)](https://libraries.io/pypi/pyani)

<!-- TOC -->

- [Overview](#overview)
- [Installation](#installation)
    - [`pip3`](#pip3)
    - [`bioconda`](#bioconda)
    - [Installing legacy BLAST](#installing-legacy-blast)
- [Citing `pyani`](#citing-pyani)
- [Documentation](#documentation)
    - [Older documentation](#older-documentation)
- [Bugs, issues, problems and questions](#bugs-issues-problems-and-questions)
- [Contributing](#contributing)
- [Walkthrough: A First Analysis](#walkthrough-a-first-analysis)
    - [1. Download genome data](#1-download-genome-data)
    - [2. Create an analysis database](#2-create-an-analysis-database)
    - [3. Conduct ANI analysis](#3-conduct-ani-analysis)
    - [4. Reporting Analyses and Analysis Results](#4-reporting-analyses-and-analysis-results)
    - [5. Generating graphical output for ANI](#5-generating-graphical-output-for-ani)
    - [6. Classifying Genomes from Analysis Results](#6-classifying-genomes-from-analysis-results)
- [Testing `pyani`](#testing-pyani)
- [Running `pyani`](#running-pyani)
    - [Script: `average_nucleotide_identity.py`](#script-average_nucleotide_identitypy)
    - [Script: `genbank_get_genomes_by_taxon.py`](#script-genbank_get_genomes_by_taxonpy)
- [DEPENDENCIES](#dependencies)
    - [For ANI analysis](#for-ani-analysis)
    - [For graphical output](#for-graphical-output)
- [Method and Output Description](#method-and-output-description)
    - [Average Nucleotide Identity (ANI)](#average-nucleotide-identity-ani)
- [Licensing](#licensing)

<!-- /TOC -->

## Overview

`pyani` is a Python3 module and script that provides support for calculating average nucleotide identity (ANI) and related measures for whole genome comparisons, and rendering relevant graphical summary output. Where available, it takes advantage of multicore systems, and can integrate with [SGE/OGE](http://gridscheduler.sourceforge.net/)-type job schedulers for the sequence comparisons.

`pyani` also installs the prgram `pyani`, which enables command-line based analysis of genomes.

## Installation

The easiest way to install `pyani` is to use `pip3` or `bioconda`:

### `pip3`

```bash
pip3 install pyani
```

From version 0.1.3.2 onwards, this should also install all the required Python package dependencies. Prior to this version (i.e. 0.1.3.1 and earlier), you can acquire these dependencies with `pip -r`, and pointing at `requirements.txt` from this repository:

```bash
pip3 install -r requirements.txt
pip3 install -r requirements-pip.txt
```

### `bioconda`

With a working `anaconda` installation, install the `bioconda` and `conda-forge` channels:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Then install `pyani`:

```bash
conda install pyani
```

**NOTE:** The current `conda` version is *behind* the development version and has a different command-line interface. Documentation for the `conda` version (0.2.x) can be found at the link below:

- [https://github.com/widdowquinn/pyani/blob/master/README_v_0_2_x.md](https://github.com/widdowquinn/pyani/blob/master/README_v_0_2_x.md)


### Installing legacy BLAST

The NCBI legacy `BLAST` package is deprecated and not recommended. However, if you wish to use `pyani blastall` or the `ANIblastall` method with the legacy `pyani` interface, you will require a locally-installed copy of this package. This is one of the packages linked in the `requirements-thirdparty.txt` file.

## Citing `pyani`

A guide to citing `pyani` is included in the file [`CITATIONS`](CITATIONS)

## Documentation

This `README.md` file provides a quick overview and walkthrough for **THE DEVELOPMENT VERSION (v0.3+) OF `pyani`**, and full documentation can be found at the link below:

- [`pyani` v0.3+ Documentation](https://readthedocs.org/projects/pyani/)

### Older documentation

If you are using an older version of `pyani` (v0.2.x), then please note that the command-line API has changed, and documentation for this version can be found at the following page:

- [`pyani` v0.2.x Documentation](./README_v_0_2_x.md)

## Bugs, issues, problems and questions

If wou would like to report a bug or problem with `pyani`, or ask a question of the developer(s), please raise an issue at the link below:

- [`pyani` Issues page](https://github.com/widdowquinn/pyani/issues)

## Contributing

If you would like to contribute code to the `pyani` project (as a bug fix or new feature), please refer to the [`CONTRIBUTING.md`](CONTRIBUTING.md) guide for more details.

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

We'll use the `pyani download` subcommand to download all available genomes for *Candidatus Blochmannia* from NCBI. The taxon ID for this grouping is [203804](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=203804&lvl=3&lin=f&keep=1&srchmode=1&unlock).

```bash
pyani download C_blochmannia --email my.email@my.domain -t 203804 -v -l C_blochmannia_dl.log
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
pyani createdb -v -l C_blochmannia_createdb.log
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
pyani anim C_blochmannia C_blochmannia_ANIm -v -l C_blochmannia_ANIm.log \
    --name "C. blochmannia run 1" \
    --labels C_blochmannia/labels.txt --classes C_blochmannia/classes.txt
```

All four analysis commands operate in a similar way. The first two arguments are paths to directories: the first path is to a directory containing input genomes, and the second is the path to an output directory for storing intermediate results. The `-v` and `-l` arguments work as above, specifying verbose output and logging output to a file.

You will probably notice that the verbose output is very verbose, to enable informative identification of any problems. In particular, the verbose output (which is also written to the log file) writes out the command-lines used for the pairwise comparisons so, if something goes wrong, you can test whether a specific comparison can be run at the command-line *at all*, to aid diagnosis of any problems.

#### Rerunning the same analysis

One reason for using a database backend for analysis results is so that, for very large analyses, we do not ever need to recalculate a pairwise genome comparison. All the analysis subcommands check whether input genomes have been used before (using the unique MD5 hash for each genome to identify whether it's been used previously), and whether the comparison of two genomes has been run, with the particular analysis settings that were used. If either genome was not seen before, or if the analysis settings are different, the comparison is performed.

You can test this for yourself by running the analysis command again, as below. You will see a number of messages indicating that genomes have been seen before, and that analyses performed before were skipped:

```bash
$ pyani anim C_blochmannia C_blochmannia_ANIm -v -l C_blochmannia_ANIm.log \
    --name "C. blochmannia run 2" \
    --labels C_blochmannia/labels.txt --classes C_blochmannia/classes.txt
INFO: command-line: pyani anim C_blochmannia C_blochmannia_ANIm -v -l C_blochmannia_ANIm.log
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

Once an analysis is run, the results are placed in a local `SQLite` database, which can be queried for information about the analyses that have been run. You can request information about:

- `--runs`: show all analysis runs with results stored in the database
- `--runs_genomes`: show all the analysis runs with results in the database, and all the genomes analysed in each run
- `--genomes`: show all the genomes used for any analysis in the database
- `--genomes_runs`: for each genome in the database, also list the analysis results it participates in
- `--run_results`: show all the pairwise comparison results for a named run (run IDs can be obtained with the `--runs` argument

The report tables are written to a named directory (compulsory argument), and are written by default `.tab` plain-text format, but HTML and Excel format can alos be requested with the `--formats` argument:

```bash
$ pyani report -v --runs C_blochmannia_ANIm/ --formats html,excel,stdout
INFO: Processed arguments: Namespace(cmdline='./pyani report -v --runs C_blochmannia_ANIm/ --formats html,excel', dbpath='.pyani/pyanidb', formats='html,excel', func=<function subcmd_report at 0x10c674a60>, logfile=None, outdir='C_blochmannia_ANIm/', run_results=False, show_genomes=False, show_genomes_runs=False, show_runs=True, show_runs_genomes=False, verbose=True)
INFO: command-line: ./pyani report -v --runs C_blochmannia_ANIm/ --formats html,excel
INFO: Creating output in formats: ['excel', 'tab', 'html']
INFO: Using database: .pyani/pyanidb
INFO: Writing table of pyani runs from the database to C_blochmannia_ANIm/runs.*
INFO: Completed. Time taken: 0.937
$ tree -L 1 C_blochmannia_ANIm/
C_blochmannia_ANIm/
├── nucmer_output
├── runs.html
├── runs.tab
└── runs.xlsx
```

To see all the pairwise results for an individual run, the run ID must be provided. It is possible to get results for more than one run ID by providing a comma-separated list of run IDs (though each run's results will be provided in a separate file):

```bash
$ pyani report -v --runs C_blochmannia_ANIm/ --formats html,excel --run_results 1,2,3,4
INFO: Processed arguments: Namespace(cmdline='./pyani report -v --runs C_blochmannia_ANIm/ --formats html,excel --run_results 1,2,3,4', dbpath='.pyani/pyanidb', formats='html,excel', func=<function subcmd_report at 0x108616a60>, logfile=None, outdir='C_blochmannia_ANIm/', run_results='1,2,3,4', show_genomes=False, show_genomes_runs=False, show_runs=True, show_runs_genomes=False, verbose=True)
INFO: command-line: ./pyani report -v --runs C_blochmannia_ANIm/ --formats html,excel --run_results 1,2,3,4
INFO: Creating output in formats: ['tab', 'excel', 'html']
INFO: Using database: .pyani/pyanidb
INFO: Writing table of pyani runs from the database to C_blochmannia_ANIm/runs.*
INFO: Attempting to write results tables for runs: ['1', '2', '3', '4']
INFO: Collecting data for run with ID: 1
INFO: Collecting data for run with ID: 2
INFO: Collecting data for run with ID: 3
INFO: Collecting data for run with ID: 4
INFO: Completed. Time taken: 1.285
```

You can see a run's results in the terminal by specifying the `stdout` format. For example, to see the identity, coverage, and other output matrices, you would specify `--run_matrices <RUN>` and `--formats=stdout` as below:

```bash
$ pyani report C_blochmannia_ANIm --formats=stdout --run_matrices 1
TABLE: C_blochmannia_ANIm/matrix_identity_1
                                                    C. Blochmannia pennsylvanicus BPEN  C. Blochmannia floridanus  C. Blochmannia vafer BVAF  C. Blochmannia chromaiodes 640  B. endosymbiont of Polyrhachis (Hedomyrma) turneri 675  B. endosymbiont of Camponotus (Colobopsis) obliquus 757
C. Blochmannia pennsylvanicus BPEN                                            1.000000                   0.834866                   0.836903                        0.980244                                           0.843700                                                0.829509
C. Blochmannia floridanus                                                     0.834866                   1.000000                   0.828733                        0.834916                                           0.847060                                                0.857859
C. Blochmannia vafer BVAF                                                     0.836903                   0.828733                   1.000000                        0.837811                                           0.866015                                                0.844438
C. Blochmannia chromaiodes 640                                                0.980244                   0.834916                   0.837811                        1.000000                                           0.849834                                                0.834769
B. endosymbiont of Polyrhachis (Hedomyrma) turn...                            0.843700                   0.847060                   0.866015                        0.849834                                           1.000000                                                0.844228
B. endosymbiont of Camponotus (Colobopsis) obli...                            0.829509                   0.857859                   0.844438                        0.834769                                           0.844228                                                1.000000

TABLE: C_blochmannia_ANIm/matrix_coverage_1
                                                    C. Blochmannia pennsylvanicus BPEN  C. Blochmannia floridanus  C. Blochmannia vafer BVAF  C. Blochmannia chromaiodes 640  B. endosymbiont of Polyrhachis (Hedomyrma) turneri 675  B. endosymbiont of Camponotus (Colobopsis) obliquus 757
C. Blochmannia pennsylvanicus BPEN                                            1.000000                   0.045736                   0.041404                        1.000306                                           0.017263                                                0.021027
C. Blochmannia floridanus                                                     0.051317                   1.000000                   0.152609                        0.054930                                           0.016366                                                0.010749
C. Blochmannia vafer BVAF                                                     0.045362                   0.149012                   1.000000                        0.046520                                           0.008356                                                0.014706
C. Blochmannia chromaiodes 640                                                1.000856                   0.048983                   0.042485                        1.000000                                           0.014056                                                0.016140
B. endosymbiont of Polyrhachis (Hedomyrma) turn...                            0.018238                   0.015410                   0.008058                        0.014841                                           1.000000                                                0.020416
B. endosymbiont of Camponotus (Colobopsis) obli...                            0.021508                   0.009799                   0.013730                        0.016500                                           0.019766                                                1.000000

TABLE: C_blochmannia_ANIm/matrix_aln_lengths_1
                                                    C. Blochmannia pennsylvanicus BPEN  C. Blochmannia floridanus  C. Blochmannia vafer BVAF  C. Blochmannia chromaiodes 640  B. endosymbiont of Polyrhachis (Hedomyrma) turneri 675  B. endosymbiont of Camponotus (Colobopsis) obliquus 757
C. Blochmannia pennsylvanicus BPEN                                            791654.0                    36207.0                    32778.0                        791896.0                                            13666.0                                                 16646.0
C. Blochmannia floridanus                                                      36207.0                   705557.0                   107674.0                         38756.0                                            11547.0                                                  7584.0
C. Blochmannia vafer BVAF                                                      32778.0                   107674.0                   722585.0                         33615.0                                             6038.0                                                 10626.0
C. Blochmannia chromaiodes 640                                                791896.0                    38756.0                    33615.0                        791219.0                                            11121.0                                                 12770.0
B. endosymbiont of Polyrhachis (Hedomyrma) turn...                             13666.0                    11547.0                     6038.0                         11121.0                                           749321.0                                                 15298.0
B. endosymbiont of Camponotus (Colobopsis) obli...                             16646.0                     7584.0                    10626.0                         12770.0                                            15298.0                                                     NaN

TABLE: C_blochmannia_ANIm/matrix_sim_errors_1
                                                    C. Blochmannia pennsylvanicus BPEN  C. Blochmannia floridanus  C. Blochmannia vafer BVAF  C. Blochmannia chromaiodes 640  B. endosymbiont of Polyrhachis (Hedomyrma) turneri 675  B. endosymbiont of Camponotus (Colobopsis) obliquus 757
C. Blochmannia pennsylvanicus BPEN                                                 0.0                     5979.0                     5346.0                         15645.0                                             2136.0                                                  2838.0
C. Blochmannia floridanus                                                       5979.0                        0.0                    18441.0                          6398.0                                             1766.0                                                  1078.0
C. Blochmannia vafer BVAF                                                       5346.0                    18441.0                        0.0                          5452.0                                              809.0                                                  1653.0
C. Blochmannia chromaiodes 640                                                 15645.0                     6398.0                     5452.0                             0.0                                             1670.0                                                  2110.0
B. endosymbiont of Polyrhachis (Hedomyrma) turn...                              2136.0                     1766.0                      809.0                          1670.0                                                0.0                                                  2383.0
B. endosymbiont of Camponotus (Colobopsis) obli...                              2838.0                     1078.0                     1653.0                          2110.0                                             2383.0                                                     0.0

TABLE: C_blochmannia_ANIm/matrix_hadamard_1
                                                    C. Blochmannia pennsylvanicus BPEN  C. Blochmannia floridanus  C. Blochmannia vafer BVAF  C. Blochmannia chromaiodes 640  B. endosymbiont of Polyrhachis (Hedomyrma) turneri 675  B. endosymbiont of Camponotus (Colobopsis) obliquus 757
C. Blochmannia pennsylvanicus BPEN                                            1.000000                   0.038183                   0.034652                        0.980543                                           0.014564                                                0.017442
C. Blochmannia floridanus                                                     0.042843                   1.000000                   0.126472                        0.045862                                           0.013863                                                0.009221
C. Blochmannia vafer BVAF                                                     0.037964                   0.123491                   1.000000                        0.038975                                           0.007237                                                0.012418
C. Blochmannia chromaiodes 640                                                0.981082                   0.040896                   0.035594                        1.000000                                           0.011945                                                0.013473
B. endosymbiont of Polyrhachis (Hedomyrma) turn...                            0.015387                   0.013053                   0.006978                        0.012613                                           1.000000                                                0.017236
B. endosymbiont of Camponotus (Colobopsis) obli...                            0.017841                   0.008406                   0.011594                        0.013774                                           0.016687                                                1.000000
```

### 5. Generating graphical output for ANI

The output of a `pyani` run can also be represented graphically, using the `plot` subcommand. For example, the command:

```bash
pyani plot C_blochmannia_ANIm 1 -v --formats png,pdf
```

will place `.pdf` and `.png` format output in the `C_blochmannia_ANIm` output directory for the run wuth ID 1, generated above. Five heatmaps are generated:

- percentage identity
- percentage coverage (for both query and subject)
- alignment length (total aligned length)
- similarity errors (total number of mismatches, not including indels)
- hadamard matrix (dot product of identity and coverage matrices)

The heatmaps also include dendrograms, clustering the rows and columns by overall similarity.

### 6. Classifying Genomes from Analysis Results

## Testing `pyani`

`pyani` includes tests that can be run with `nosetest` (including coverage measurement using `coverage.py`) with the following command, executed from the repository root directory:

```bash
nosetests --cover-erase --cover-package=pyani --cover-html --with-coverage
```

Coverage output will be placed (by default) in the `cover` subdirectory, and can be loaded into the web browser.

## Running `pyani`

### Script: `average_nucleotide_identity.py`

The `average_nucleotide_identity.py` script - installed as part of this package - enables straightforward ANI analysis at the command-line, and uses the `pyani` module behind the scenes.

You can get a summary of available command-line options with `average_nucleotide_identity.py -h`

```bash
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

```bash
./average_nucleotide_identity.py -i tests/test_ani_data/ -o tests/test_ANIm_output -m ANIm -g
./average_nucleotide_identity.py -i tests/test_ani_data/ -o tests/test_ANIb_output -m ANIb -g
./average_nucleotide_identity.py -i tests/test_ani_data/ -o tests/test_ANIblastall_output -m ANIblastall -g
./average_nucleotide_identity.py -i tests/test_ani_data/ -o tests/test_TETRA_output -m TETRA -g
```

The graphical output below, supporting assignment of `NC_002696` and `NC_011916` to the same species (*C.crescentus*), and the other two isolates to distinct species (`NC_014100`:*C.segnis*; `NC_010338`:*C.* sp K31), was generated with the command-line:

```bash
./average_nucleotide_identity.py -v -i tests/test_ani_data/ \
    -o tests/test_ANIm_output/ -g --gformat png,pdf,eps \
    --classes tests/test_ani_data/classes.tab \
    --labels tests/test_ani_data/labels.tab
```

![ANIm percentage identity for *Caulobacter* test data](tests/test_ani_data/ANIm_percentage_identity.png "ANIm percentage identity")
![ANIm alignment coverage for *Caulobacter* test data](tests/test_ani_data/ANIm_alignment_coverage.png "ANIm alignment coverage")
![ANIm alignment length for *Caulobacter* test data](tests/test_ani_data/ANIm_alignment_lengths.png "ANIm alignment length")
![ANIm alignment similarity errors for *Caulobacter* test data](tests/test_ani_data/ANIm_similarity_errors.png "ANIm alignment similarity")

### Script: `genbank_get_genomes_by_taxon.py`

The script `genbank_get_genomes_by_taxon.py`, installed by this package, enables download of genomes from NCBI, specified by taxon ID. The script will download all available assemblies for taxa at or below the specified node in the NCBI taxonomy tree.

Command-line options can be viewed using:

```bash
$ genbank_get_genomes_by_taxon.py -h
usage: genbacnk_get_genomes_by_taxon.py [-h] [-o OUTDIRNAME] [-t TAXON] [-v]
                                        [-f] [--noclobber] [-l LOGFILE]
                                        [--format FORMAT] [--email EMAIL]
                                        [--retries RETRIES]
                                        [--batchsize BATCHSIZE]
[…]
```

For example, the NCBI taxonomy ID for *Caulobacter* is 75, so all publicly-available *Caulobacter* sequences can be obtained using the command-line:

```bash
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

- **Biopython** <http://www.biopython.org>
- **NumPy** <http://www.numpy.org/>
- **pandas** <http://pandas.pydata.org/>
- **SciPy** <http://www.scipy.org/>

#### Alignment tools

- **BLAST+** executable in the `$PATH`, or available on the command line (required for **ANIb** analysis) <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>
- **legacy BLAST** executable in the `$PATH` or available on the command line (required for **ANIblastall** analysis) <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/>
- **MUMmer** executables in the $PATH, or available on the command line (required for **ANIm** analysis) <http://mummer.sourceforge.net/>

### For graphical output

- **matplotlib** <http://matplotlib.org/>
- **seaborn** <https://github.com/mwaskom/seaborn>

## Method and Output Description

### Average Nucleotide Identity (ANI)

This module calculates Average Nucleotide Identity (ANI) according to one of a number of alternative methods described in, e.g.

- Richter M, Rossello-Mora R (2009) Shifting the genomic gold standard for the prokaryotic species definition. Proc Natl Acad Sci USA 106: 19126-19131. doi:10.1073/pnas.0906412106. (ANI1020, ANIm, ANIb)
- Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, et al. (2007) DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. Int J Syst Evol Micr 57: 81-91. doi:10.1099/ijs.0.64483-0.

ANI is proposed to be the appropriate *in silico* substitute for DNA-DNA hybridisation (DDH), and so useful for delineating species boundaries. A typical percentage threshold for species boundary in the literature is 95% ANI (e.g. Richter et al. 2009).

All ANI methods follow the basic algorithm:

- Align the genome of organism 1 against that of organism 2, and identify the matching regions
- Calculate the percentage nucleotide identity of the matching regions, as an average for all matching regions

Methods differ on: (1) what alignment algorithm is used, and the choice of parameters (this affects the aligned region boundaries); (2) what the input is for alignment (typically either fragments of fixed size, or the most complete assembly available).

- **ANIm**: uses MUMmer (NUCmer) to align the input sequences.
- **ANIb**: uses BLASTN+ to align 1020nt fragments of the input sequences
- **ANIblastall**: uses legacy BLASTN to align 1020nt fragments of the input sequences
- **TETRA**: calculates tetranucleotide frequencies of each input sequence

The algorithms takes as input correctly-formatted FASTA multiple sequence files. All sequences for a single organism should be contained in only one sequence file. Although it is possible to provide new labels for each input genome, for rendering graphical output, the names of these files are used for identification so it is best to name them sensibly.

Output is written to a named directory. The output files differ depending on the chosen ANI method.

- **ANIm**: MUMmer/NUCmer .delta files, describing each pairwise sequence alignment. Output as tab-separated plain text format tables describing: alignment coverage; total alignment lengths; similarity errors; and percentage identity (ANIm).
- **ANIb** and **ANIblastall**: FASTA sequences describing 1020nt fragments of each input sequence; BLAST nucleotide databases - one for each set of fragments; and BLASTN output files (tab-separated tabular format plain text) - one for each pairwise comparison of input sequences. Output as tab-separated plain text tables describing: alignment coverage; total alignment lengths; similarity errors; and percentage identity (ANIb or ANIblastall).
- **TETRA**: Tab-separated plain text files describing the Pearson correlations between Z-score distributions for each tetranucleotide in each input sequence (TETRA).

If graphical output is chosen, the output directory will also contain PDF, PNG and EPS files representing the various output measures as a heatmap with row and column dendrograms. Other output formats (e.g. SVG) can be specified with the `--gformat` argument.

## Licensing

Unless otherwise indicated, all code is subject to the following agreement:

```text
    (c) The James Hutton Institute 2014-2019
    (c) The University of Strathclyde 2019
    Author: Leighton Pritchard

    Contact: leighton.pritchard@strath.ac.uk

    Address:
    Leighton Pritchard,
    Strathclyde Institute of Pharmacy and Biomedical Sciences
    161 Cathedral Street
    Glasgow
    G4 0RE,
    Scotland,
    UK

The MIT License

Copyright (c) 2014-2019 The James Hutton Institute

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
```
