# README.md (pyani)

## Overview
`pyani` is a Python module that provides support for calculating average nucleotide identity (ANI) and related measures for whole genome comparisons, and rendering relevant graphical summary output. Where available, it takes advantage of multicore systems, and can integrate with [SGE/OGE](http://gridscheduler.sourceforge.net/)-type job schedulers for the sequence comparisons.

`pyani` also installs a script: `average_nucleotide_identity.py` that enables command-line ANI analysis.

## Script: <a name="average_nucleotide_identity.py">`average_nucleotide_identity.py`</a>

The `average_nucleotide_identity.py` script - part of this module - enables ANI analysis at the command-line, and uses the `pyani` module behind the scenes.

### Script Usage

```
average_nucleotide_identity.py [-h] [-o OUTDIRNAME] [-i INDIRNAME] [-v]
                                      [-f] [-s] [-l LOGFILE] [--skip_nucmer]
                                      [--skip_blastn] [--noclobber] [-g]
                                      [--gformat GFORMAT] [--gmethod GMETHOD]
                                      [--labels LABELS] [--classes CLASSES]
                                      [-m METHOD] [--scheduler SCHEDULER]
                                      [--maxmatch] [--nucmer_exe NUCMER_EXE]
                                      [--blastn_exe BLASTN_EXE]
                                      [--makeblastdb_exe MAKEBLASTDB_EXE]
                                      [--blastall_exe BLASTALL_EXE]
                                      [--formatdb_exe FORMATDB_EXE]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIRNAME, --outdir OUTDIRNAME
                        Output directory
  -i INDIRNAME, --indir INDIRNAME
                        Input directory name
  -v, --verbose         Give verbose output
  -f, --force           Force file overwriting
  -s, --fragsize        Sequence fragment size for ANIb
  -l LOGFILE, --logfile LOGFILE
                        Logfile location
  --skip_nucmer         Skip NUCmer runs, for testing (e.g. if output already
                        present)
  --skip_blastn         Skip BLASTN runs, for testing (e.g. if output already
                        present)
  --noclobber           Don't nuke existing files
  -g, --graphics        Generate heatmap of ANI
  --gformat GFORMAT     Graphics output format [pdf|png|jpg|svg]
  --gmethod GMETHOD     Graphics output method [mpl|R]
  --labels LABELS       Path to file containing sequence labels
  --classes CLASSES     Path to file containing sequence classes
  -m METHOD, --method METHOD
                        ANI method [ANIm|ANIb|ANIblastall|TETRA]
  --scheduler SCHEDULER
                        Job scheduler [multiprocessing|SGE]
  --maxmatch            Override MUMmer to allow all NUCmer matches
  --nucmer_exe NUCMER_EXE
                        Path to NUCmer executable
  --blastn_exe BLASTN_EXE
                        Path to BLASTN+ executable
  --makeblastdb_exe MAKEBLASTDB_EXE
                        Path to BLAST+ makeblastdb executable
  --blastall_exe BLASTALL_EXE
                        Path to BLASTALL executable
  --formatdb_exe FORMATDB_EXE
                        Path to BLAST formatdb executable
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
    -o tests/test_ANIm_output/ -g --gformat png \
    --classes tests/test_ani_data/classes.tab \
    --labels tests/test_ani_data/labels.tab
```


![ANIm percentage identity for *Caulobacter* test data](tests/test_ani_data/ANIm_percentage_identity.png "ANIm percentage identity")
![ANIm alignment coverage for *Caulobacter* test data](tests/test_ani_data/ANIm_alignment_coverage.png "ANIm alignment coverage")
![ANIm alignment length for *Caulobacter* test data](tests/test_ani_data/ANIm_alignment_lengths.png "ANIm alignment length")
![ANIm alignment similarity errors for *Caulobacter* test data](tests/test_ani_data/ANIm_similarity_errors.png "ANIm alignment similarity")

## DEPENDENCIES

### For ANI analysis

* **Biopython** <http://www.biopython.org>
* **NumPy** <http://www.numpy.org/>
* **pandas** <http://pandas.pydata.org/>
* **SciPy** <http://www.scipy.org/>

* **BLAST+** executable in the `$PATH`, or available on the command line (required for **ANIb** analysis) <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>
* **legacy BLAST** executable in the `$PATH` or available on the command line (required for **ANIblastall** analysis) <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/>
* **MUMmer** executables in the $PATH, or available on the command line (required for **ANIm** analysis) <http://mummer.sourceforge.net/>

### For graphical output

* **matplotlib** <http://matplotlib.org/>

and/or

* **R** with shared libraries installed on the system <http://cran.r-project.org/>
* **Rpy2** <http://rpy.sourceforge.net/rpy2.html>


# INSTALLATION

If you have downloaded v0.1.0 or greater, and the dependencies above are satisfied, then installation should be as simple as downloading the latest release and uncompressing it, or cloning the repository as below:

```
$ git clone https://github.com/widdowquinn/pyani
```

then changing to the appropriate directory:

```
$ cd pyani
```

then issuing:

```
$ python setup.py install
```

(or whatever variant you wish, e.g. for a home directory-local installation) from the top directory in the repository (with root permissions, if necessary).

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

* **ANIm**: MUMmer/NUCmer .delta files, describing each pairwise sequence alignment. Output as tab-separated plain text and Excel format tables describing: alignment coverage; total alignment lengths; similarity errors; and percentage identity (ANIm).
* **ANIb** and **ANIblastall**: FASTA sequences describing 1020nt fragments of each input sequence; BLAST nucleotide databases - one for each set of fragments; and BLASTN output files (tab-separated tabular format plain text) - one for each pairwise comparison of input sequences. Output as tab-separated plain text and Excel format tables describing: alignment coverage; total alignment lengths; similarity errors; and percentage identity (ANIb or ANIblastall).
* **TETRA**: Tab-separated plain text and Excel format files describing the Pearson correlations between Z-score distributions for each tetranucleotide in each input sequence (TETRA).

If graphical output is chosen, the output directory will also contain PDF files representing the similarity between sequences as a heatmap with row and column dendrograms.



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
