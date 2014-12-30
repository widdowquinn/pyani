#!/usr/bin/env python
#
# calculate_ani.py
#
# This script calculates Average Nucleotide Identity (ANI) according to one of
# a number of alternative methods described in, e.g.
#
# Richter M, Rossello-Mora R (2009) Shifting the genomic gold standard for the
# prokaryotic species definition. Proc Natl Acad Sci USA 106: 19126-19131.
# doi:10.1073/pnas.0906412106. (ANI1020, ANIm, ANIb)
#
# Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, et al.
# (2007) DNA-DNA hybridization values and their relationship to whole-genome
# sequence similarities. Int J Syst Evol Micr 57: 81-91.
# doi:10.1099/ijs.0.64483-0.
#
# ANI is proposed to be the appropriate in silico substitute for DNA-DNA
# hybridisation (DDH), and so useful for delineating species boundaries. A
# typical percentage threshold for species boundary in the literature is 95%
# ANI (e.g. Richter et al. 2009).
#
# All ANI methods follow the basic algorithm:
# - Align the genome of organism 1 against that of organism 2, and identify
#   the matching regions
# - Calculate the percentage nucleotide identity of the matching regions, as
#   an average for all matching regions
# Methods differ on: (1) what alignment algorithm is used, and the choice of
# parameters (this affects the aligned region boundaries); (2) what the input
# is for alignment (typically either fragments of fixed size, or the most
# complete assembly available).
#
# ANIm: uses MUMmer (NUCmer) to align the input sequences.
# ANIb: uses BLASTN to align 1000nt fragments of the input sequences
# TETRA: calculates tetranucleotide frequencies of each input sequence
#
# This script takes as input a directory containing a set of
# correctly-formatted FASTA multiple sequence files. All sequences for a
# single organism should be contained in only one sequence file. The names of
# these files are used for identification, so it would be advisable to name
# them sensibly.
#
# Output is written to a named directory. The output files differ depending on
# the chosen ANI method.
#
# ANIm: MUMmer/NUCmer .delta files, describing the sequence
#       alignment; tab-separated format plain text tables describing total
#       alignment lengths, and total alignment percentage identity
#
# ANIb: FASTA sequences describing 1000nt fragments of each input sequence;
#       BLAST nucleotide databases - one for each set of fragments; and BLASTN
#       output files (tab-separated tabular format plain text) - one for each
#       pairwise comparison of input sequences. There are potentially a lot of
#       intermediate files.
#
# TETRA: Tab-separated text file describing the Z-scores for each
#        tetranucleotide in each input sequence.
#
# In addition, all methods produce a table of output percentage identity (ANIm
# and ANIb) or correlation (TETRA), between each sequence.
#
# If graphical output is chosen, the output directory will also contain PDF
# files representing the similarity between sequences as a heatmap with
# row and column dendrograms.
#
# DEPENDENCIES
# ============
#
# o Biopython (http://www.biopython.org)
#
# o BLAST+ executable in the $PATH, or available on the command line (ANIb)
#       (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
#
# o MUMmer executables in the $PATH, or available on the command line (ANIm)
#       (http://mummer.sourceforge.net/)
#
# For graphical output
# --------------------
#
# o R with shared libraries installed on the system, for graphical output
#       (http://cran.r-project.org/)
#
# o Rpy2 (http://rpy.sourceforge.net/rpy2.html)
#
#
# USAGE
# =====
#
# calculate_ani.py [options]
#
# Options:
#   -h, --help            show this help message and exit
#   -o OUTDIRNAME, --outdir=OUTDIRNAME
#                         Output directory
#   -i INDIRNAME, --indir=INDIRNAME
#                         Input directory name
#   -v, --verbose         Give verbose output
#   -f, --force           Force file overwriting
#   -s, --fragsize        Sequence fragment size for ANIb
#   --skip_nucmer         Skip NUCmer runs, for testing (e.g. if output already
#                         present)
#   --skip_blast          Skip BLAST runs, for testing (e.g. if output already
#                         present)
#   --noclobber           Don't nuke existing files
#   -g, --graphics        Generate heatmap of ANI
#   -m METHOD, --method=METHOD
#                         ANI method
#   --maxmatch            Override MUMmer settings and allow all matches in
#                         NUCmer
#   --nucmer_exe=NUCMER_EXE
#                         Path to NUCmer executable
#   --blast_exe=BLAST_EXE
#                         Path to BLASTN+ executable
#   --makeblastdb_exe=MAKEBLASTDB_EXE
#                         Path to BLAST+ makeblastdb executable
#
# (c) The James Hutton Institute 2013
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@hutton.ac.uk
#
# Leighton Pritchard,
# Information and Computing Sciences,
# James Hutton Institute,
# Errol Road,
# Invergowrie,
# Dundee,
# DD6 9LH,
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2010-2014 The James Hutton Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

#=============
# GLOBALS

# Store the cumulative sum of returned values from multiprocessing runs.
# If external apps are well-behaved, they should return 0 on successful
# execution, and some other integer otherwise. So, by keeping a running total
# of returned values from the multiprocessing callback, we can establish
# whether one or more of the external applications had an issue.
# We use this in conjunction with thrown errors to suggest possible causes
# and solutions for the error.
CUM_RETVALS = 0

#=============
# IMPORTS

import collections
import logging
import logging.handlers
import math
import multiprocessing
import os
import shutil
import subprocess
import sys
import time
import traceback

try:
    from Bio import SeqIO
except ImportError:
    print "Biopython required for script, but not found (exiting)"
    sys.exit(1)
from optparse import OptionParser

try:
    import rpy2.robjects as robjects
    rpy2_import = True
except ImportError:
    print "Could not import rpy2: graphical output " +\
        "and TETRA method not available"
    rpy2_import = False

#=============
# FUNCTIONS


# Parse command-line
def parse_cmdline(args):
    """ Parse command-line arguments for the script
    """
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outdir", dest="outdirname",
                      action="store", default=None,
                      help="Output directory")
    parser.add_option("-i", "--indir", dest="indirname",
                      action="store", default=None,
                      help="Input directory name")
    parser.add_option("-v", "--verbose", dest="verbose",
                      action="store_true", default=False,
                      help="Give verbose output")
    parser.add_option("-f", "--force", dest="force",
                      action="store_true", default=False,
                      help="Force file overwriting")
    parser.add_option("-s", "--fragsize", dest="fragsize",
                      action="store_const", default=1020,
                      help="Sequence fragment size for ANIb")
    parser.add_option("-l", "--logfile", dest="logfile",
                      action="store", default=None,
                      help="Logfile location")
    parser.add_option("--skip_nucmer", dest="skip_nucmer",
                      action="store_true", default=False,
                      help="Skip NUCmer runs, for testing " +
                      "(e.g. if output already present)")
    parser.add_option("--skip_blast", dest="skip_blast",
                      action="store_true", default=False,
                      help="Skip BLAST runs, for testing " +
                      "(e.g. if output already present)")
    parser.add_option("--noclobber", dest="noclobber",
                      action="store_true", default=False,
                      help="Don't nuke existing files")
    parser.add_option("-g", "--graphics", dest="graphics",
                      action="store_true", default=False,
                      help="Generate heatmap of ANI")
    parser.add_option("--format", dest="gformat",
                      action="store", default="pdf",
                      help="Graphics output format [pdf|png|jpg]")
    parser.add_option("-m", "--method", dest="method",
                      action="store", default="ANIm",
                      help="ANI method [ANIm|ANIb|TETRA]")
    parser.add_option("--maxmatch", dest="maxmatch",
                      action="store_true", default=False,
                      help="Override MUMmer to allow all NUCmer matches")
    parser.add_option("--nucmer_exe", dest="nucmer_exe",
                      action="store", default="nucmer",
                      help="Path to NUCmer executable")
    parser.add_option("--blast_exe", dest="blast_exe",
                      action="store", default="blastn",
                      help="Path to BLASTN+ executable")
    parser.add_option("--makeblastdb_exe", dest="makeblastdb_exe",
                      action="store", default="makeblastdb",
                      help="Path to BLAST+ makeblastdb executable")
    return parser.parse_args()


# Report last exception as string
def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# METHOD: ANIm
# This method uses NUCmer to calculate pairwise alignments for the input
# organisms, without chopping sequences into fragments. We follow the method
# of Richter et al. (2009)
def calculate_anim():
    """ Calculate ANI by the ANIm method, as described in Richter et al (2009)
        Proc Natl Acad Sci USA 106: 19126-19131 doi:10.1073/pnas.0906412106.

        All FASTA format files (selected by suffix) in the input directory
        are compared against each other, pairwise, using NUCmer (which must
        be in the path). NUCmer output is stored in the output directory.

        The NUCmer .delta file output is parsed to obtain an alignment length
        and similarity error count for every unique region alignment between
        the two organisms, as represented by the sequences in the FASTA files.

        These are processed to give matrices of aligned sequence lengths,
        similarity error counts, average nucleotide identity (ANI) percentages,
        and minimum aligned percentage (of whole genome) for each pairwise
        comparison.

        The matrices are written to file in a plain text tab-separated format.
    """
    logger.info("Running ANIm method")
    infiles = get_fasta_files(options.indirname)
    org_lengths = get_org_lengths()
    pairwise_nucmer(infiles)
    lengths, sim_errors, perc_ids, perc_aln = process_delta(org_lengths)
    # Sanity check print for organisms of same species
    #for k, v in sorted(perc_ids.items()):
    #    if v > 0.95:
    #        print k, v
    # Write output to file
    write_table('aln_lengths.tab', org_lengths.keys(), lengths,
                "Aligment lengths")
    write_table('sim_errors.tab', org_lengths.keys(), sim_errors,
                "Similarity errors")
    write_table('perc_ids.tab', org_lengths.keys(), perc_ids,
                "ANIm", "float")
    write_table('perc_aln.tab', org_lengths.keys(), perc_aln,
                "Minimum % aligned nt", "float")


# METHOD: ANIb
# This method uses BLAST to calculate pairwise alignments for input organisms,
# fragmented into consecutive 1020bp fragments, as described in Goris et al.
# (2007).
def calculate_anib():
    """ Calculate ANI by the ANIb method, as described in Goris et al. (2007)
        Int J Syst Evol Micr 57: 81-91. doi:10.1099/ijs.0.64483-0.

        All FASTA format files (selected by suffix) in the input directory are
        used to construct BLAST databases, placed in the output directory.
        Each file's contents are also split into sequence fragments of length
        options.fragsize, and the multiple FASTA file that results written to
        the output directory. These are BLASTNed, pairwise, against the
        databases.

        The BLAST output is interrogated for all fragment matches that cover
        at least 70% of the query sequence, with at least 30% nucleotide
        identity over the full length of the query sequence. This is an odd
        choice and doesn't correspond to the twilight zone limit as implied by
        Goris et al. We persist with their definition, however.  Only these
        qualifying matches contribute to the total aligned length, and total
        aligned sequence identity used to calculate ANI.

        The results are processed to give matrices of aligned sequence length
        (aln_lengths.tab), similarity error counts (sim_errors.tab), ANIs
        (perc_ids.tab), and minimum aligned percentage (perc_aln.tab) of
        each genome, for each pairwise comparison. These are written to the
        output directory in plain text tab-separated format.
    """
    logger.info("Running ANIb method")
    fragment_input_files()
    org_lengths = get_org_lengths()
    make_blast_dbs()
    pairwise_blast(get_input_files(options.outdirname, '.fasta'))
    lengths, sim_errors, perc_ids, perc_aln = process_blast(org_lengths)
    # Sanity check print for organisms of same species
    for k, v in sorted(perc_ids.items()):
        if v > 0.95:
            print k, v
    # Write output to file
    write_table('aln_lengths.tab', org_lengths.keys(), lengths,
                "Aligment lengths")
    write_table('sim_errors.tab', org_lengths.keys(), sim_errors,
                "Similarity errors")
    write_table('perc_ids.tab', org_lengths.keys(), perc_ids,
                "ANIb", "float")
    write_table('perc_aln.tab', org_lengths.keys(), perc_aln,
                "Minimum % aligned nt", "float")


# METHOD: TETRA
# This method calculates tetranucleotide frequencies for the input organisms,
# as used by JSpecies, and described in Richter et al (2009) and Teeling et
# al. (2004) and Teeling et al. (2004)
def calculate_tetra():
    """ Calculate tetranucleotide frequencies for each input sequence, and
        their Pearson correlation, as described in Teeling et al. (2004a)
        Env. Micro. 6 938-947 doi:10.1111/j.1462-2920.2004.00624.x;
        Teeling et al. (2004b) BMC Bioinf. 5 163 doi:10.1186/1471-2105-5-163;
        and Richter & Rossello-Mora (2009) PNAS 106 19126-19131
        doi:10.1073/pnas.0906412106.

        FASTA format files (selected by suffix) in the input directory are used
        to construct frequency tables of tetranucleotide occurrence. These
        frequencies are converted into Z-scores, and tabulated
        (rows=tetranucleotides, columns=input sequences, written as
        tetra_z_scores.tab). The Pearson correlation between Z-scores for
        input sequences is then used as a measure of sequence similarity and
        written to tetra_corr.tab.
    """
    logger.info("Running TETRA method")
    tetra_z = calc_org_tetra()   # Calculate Z-scores for tetranucleotides
    write_tetraz('tetra_z_scores.tab', tetra_z)
    corr_z = calc_tetra_corr(tetra_z)  # Calculate Pearson correlation
    write_table('tetra_corr.tab', tetra_z.keys(), corr_z, "TETRA", "float")


# SUPPORT FUNCTIONS
# Write the set of tetranucleotide frequency Z scores to a plain text
# tab-separated table, in the output directory
def write_tetraz(filename, tetra_z):
    """ Writes the Z score for each tetranucleotide to a plain text
        tab-separated table with one row for each tetranucleotide, and one
        column for each input sequence.

        - filename is the location of the file to which the tetranucleotide
              Z-scores should be written

        - tetra_z is a dictionary containing the tetranucleotide Z-scores
              for each input sequence
    """
    try:
        fh = open(os.path.join(options.outdirname, filename), 'w')
        logger.info("Writing Z-scores to %s" % fh.name)
    except:
        logger.error("Could not open %s for writing (exiting)")
        sys.exit(1)
    orgs = sorted(tetra_z.keys())
    tets = sorted(tetra_z.values()[0].keys())
    # Write headers
    print >> fh, "# calculate_ani.py %s" % time.asctime()
    print >> fh, "# tetranucleotide frequency Z-scores"
    print >> fh, '\t' + '\t'.join(orgs)
    for tet in tets:
        outstr = [tet] + ["%.2f" % tetra_z[org][tet] for org in orgs]
        print >> fh, '\t'.join(outstr)
    fh.close()


# Calculate Pearson's correlation coefficient from the Z-scores for each
# tetranucleotide. If we're forcing rpy2, might as well use that, though...
def calc_tetra_corr(tetra_z):
    """ Calculate Pearson correlation coefficient from Z scores for each
        tetranucleotide. This is done longhand here, which is fast enough,
        but for robustness we might want to just hand this over to R,
        as it is a dependency for this method anyway (TODO).

        Note that we report a correlation by this method, rather than a
        percentage identity.

        - tetra_z is a dictionary of tetranucleotide Z-scores, for each
              input sequence
    """
    corrs = {}
    orgs = sorted(tetra_z.keys())
    for idx, o1 in enumerate(orgs[:-1]):
        for o2 in orgs[idx+1:]:
            assert sorted(tetra_z[o1].keys()) == sorted(tetra_z[o2].keys())
            tets = sorted(tetra_z[o1].keys())
            z1 = [tetra_z[o1][t] for t in tets]
            z2 = [tetra_z[o2][t] for t in tets]
            z1_mean = sum(z1) / len(z1)
            z2_mean = sum(z2) / len(z2)
            z1diffs = [z - z1_mean for z in z1]
            z2diffs = [z - z2_mean for z in z2]
            diffprods = sum([z1diffs[i] * z2diffs[i] for i in
                             range(len(z1diffs))])
            z1diffs2 = sum([z * z for z in z1diffs])
            z2diffs2 = sum([z * z for z in z2diffs])
            corrs[(o1, o2)] = diffprods/math.sqrt(z1diffs2 * z2diffs2)
    return corrs


# Calculate tetranucleotide values for each input sequence
def calc_org_tetra():
    """ We calculate the mono-, di-, tri- and tetranucleotide frequencies
        for each sequence, on each strand, and follow Teeling et al. (2004)
        in calculating a corresponding Z-score for each observed
        tetranucleotide frequency, dependent on the mono-, di- and tri-
        nucleotide frequencies for that input sequence.
    """
    org_tetraz = {}
    for fn in get_fasta_files(options.indirname):
        org = os.path.splitext(os.path.split(fn)[-1])[0]
        logger.info("Calculating tetranucleotide frequencies for %s" % fn)
        # For the Teeling et al. method, the Z-scores require us to count
        # mono, di, tri and tetranucleotide sequences
        monocnt, dicnt, tricnt, tetracnt = (collections.defaultdict(int),
                                            collections.defaultdict(int),
                                            collections.defaultdict(int),
                                            collections.defaultdict(int))
        for rec in SeqIO.parse(fn, 'fasta'):
            for s in [str(rec.seq).upper(),
                      str(rec.seq.reverse_complement()).upper()]:
                # Since the Teeling et al. algorithm requires us to consider
                # both strand orientations, monocounts are easy
                monocnt['G'] += s.count('G')
                monocnt['C'] += s.count('C')
                monocnt['T'] += s.count('T')
                monocnt['A'] += s.count('A')
                # For di, tri and tetranucleotide counts, we loop over the
                # sequence and its reverse complement, until near the end:
                for i in range(len(s[:-4])):
                    di, tri, tetra = s[i:i+2], s[i:i+3], s[i:i+4]
                    dicnt[str(di)] += 1
                    tricnt[str(tri)] += 1
                    tetracnt[str(tetra)] += 1
                # We clean up the straggling bit at the end:
                tricnt[str(s[-4:-1])] += 1
                tricnt[str(s[-3:])] += 1
                dicnt[str(s[-4:-2])] += 1
                dicnt[str(s[-3:-1])] += 1
                dicnt[str(s[-2:])] += 1
        logger.info("%d mono, %d di, %d tri, %d tetranucleotides found" %
                    (len(monocnt), len(dicnt), len(tricnt), len(tetracnt)))
        # Following Teeling (2004), we calculate expected frequencies for each
        # tetranucleotide; we ignore ambiguity symbols
        tetra_exp = {}
        for t in [tet for tet in tetracnt if tet_clean(tet)]:
            tetra_exp[t] = 1.*tricnt[t[:3]]*tricnt[t[1:]]/dicnt[t[1:3]]
        logger.info("%d non-ambiguous tetranucleotides" % len(tetra_exp))
        # Following Teeling (2004) we approximate the std dev of each
        # tetranucleotide
        tetra_sd = {}
        for t, exp in tetra_exp.items():
            den = dicnt[t[1:3]]
            tetra_sd[t] = math.sqrt(exp * (den - tricnt[t[:3]]) *
                                    (den - tricnt[t[1:]]) / (den * den))
        # Following Teeling (2004) we calculate the Z-score for each
        # tetranucleotide
        tetra_z = {}
        for t, exp in tetra_exp.items():
            try:
                tetra_z[t] = (tetracnt[t] - exp)/tetra_sd[t]
            except ZeroDivisionError:
                # We hit a zero in the estimation of variance
                zeroes = [k for k, v in tetra_sd.items() if v == 0]
                logger.warning("Zero variance for tetranucleotides %s" %
                               zeroes)
                tetra_z[t] = 1 / (dicnt[t[1:3]] * dicnt[t[1:3]])
        #print len(tetra_z), sorted(tetra_z.items())
        #for tet in sorted(tetra_exp.keys()):
        #    print tet, tetracnt[tet], tetra_exp[tet], tetra_sd[tet], \
        #        tetra_z[tet]
        org_tetraz[org] = tetra_z
    return org_tetraz


# Returns true if the passed string contains only A, C, G or T
def tet_clean(s):
    """ Checks that a passed string contains only unambiguous IUPAC nucleotide
        symbols. We are assuming that a low frequency of IUPAC ambiguity
        symbols doesn't affect our calculation.
    """
    if not len(set(s) - set('ACGT')):
        return True
    return False


# Divide the input FASTA sequences into fragments, and place multiple sequence
# FASTA files into the output directory
def fragment_input_files():
    """ Takes every sequence from every FASTA file in the input directory,
        splits them into consecutive fragments of length options.fragsize,
        (with any trailing sequences being included, even if shorter), and
        writes the resulting set of sequences to a file with the same name
        in the output directory. All fragments are named consecutively and
        uniquely as fragNNNNN.
    """
    logger.info("Fragmenting input FASTA files")
    for fn in get_fasta_files(options.indirname):
        logger.info("Processing %s" % fn)
        ostem = os.path.splitext(os.path.split(fn)[-1])[0]
        ofn = os.path.join(options.outdirname, ostem) + '.fasta'
        logger.info("Writing fragments to %s" % ofn)
        outseqs = []
        i, count = 0, 0
        for s in SeqIO.parse(fn, 'fasta'):
            while i < len(s):
                count += 1
                newseq = s[i:i+options.fragsize]
                newseq.id = "frag%05d" % count
                outseqs.append(newseq)
                i += options.fragsize
        SeqIO.write(outseqs, ofn, 'fasta')


# Make BLAST databases for each of the fragmented input files
def make_blast_dbs():
    """ Use local makeblastdb to build BLAST a nucleotide database for each
        input sequence.

        For ANIb, the input sequence has been split into consecutive fragments.
    """
    global CUM_RETVALS
    logger.info("Making BLAST databases for fragment files")
    cmdlines = []
    for fn in get_fasta_files(options.indirname):
        cmdlines.append(make_makeblastdb_cmd(fn))
    logger.info("BLAST makeblastdb command lines:\n\t%s" %
                '\n\t'.join(cmdlines))
    multiprocessing_run(cmdlines)
    if 0 < CUM_RETVALS:           # At least one error
        logger.error("At least one BLAST database creation step had a " +
                     "nonzero return value. ANIb may fail.")
        CUM_RETVALS = 0           # Reset error


# Get the list of FASTA files from the input directory
def get_fasta_files(dirname=None):
    """ Return a list of FASTA files in the input directory
    """
    infiles = get_input_files(dirname,  # '.fna')
                              '.fasta', '.fas', '.fa', '.fna')
    logger.info("Input files:\n\t%s" % '\n\t'.join(infiles))
    return infiles


# Get lengths of sequence for each organism
def get_org_lengths():
    """ Returns a dictionary of total input sequence lengths, keyed by
        organism.

        Biopython's SeqIO module is used to parse all sequences in the FASTA
        file corresponding to each organism, and the total base count in each
        is obtained.

        NOTE: ambiguity symbols are not discounted.
    """
    logger.info("Processing input organism sequence lengths")
    tot_lengths = {}
    for fn in get_fasta_files(options.indirname):
        tot_lengths[os.path.splitext(os.path.split(fn)[-1])[0]] = \
            sum([len(s) for s in SeqIO.parse(fn, 'fasta')])
    return tot_lengths


# Write a table of values to file, organised as a square matrix with row/col
# headers, in tab-separated format
def write_table(filename, org_names, values, comment='', fmt=None):
    """ Writes a tab-separated plain text square matrix file, with row and
        column headers, describing the passed data.

        - filename is the full path for the output file

        - org_names describes the row and column headers, as a list of strings

        - values describes the values in each cell, as a dictionary keyed by
             (row, col) tuple, with identifiers from org_names

        - comment is an optional comment string for the output file
    """
    # Open file and write file header
    fname = os.path.join(options.outdirname, filename)
    try:
        logger.info("Opening %s for writing" % fname)
        fh = open(fname, 'w')
        print >> fh, "# calculate_ani.py %s" % time.asctime()
        if len(comment):
            print >> fh, "# %s" % comment
    except:
        logger.error("Could not open file %s for output (exiting)" % fname)
        logger.error(last_exception())
        sys.exit(1)
    names = sorted(list(set(org_names)))  # Set name order
    print >> fh, '\t'.join([''] + names)
    for n1 in names:
        outrow = [n1]
        for n2 in names:
            if n1 == n2:
                outrow.append('NA')
                continue
            try:
                val = values[(n1, n2)]
            except KeyError:  # This error is not thrown for a square matrix
                val = values[(n2, n1)]
            if fmt == "float":
                outrow.append("%.4f" % val)
            else:
                outrow.append(str(val))
        print >> fh, '\t'.join(outrow)
    logger.info("Wrote data to %s" % fname)


# Parse NUCmer delta output to store alignment total length, sim_error,
# and percentage identity, for each pairwise comparison
def process_delta(org_lengths):
    """ Returns a tuple containing a list and four dictionaries. The list
        describes the names of all organisms (derived from their filenames).
        The dictionaries describe results for pairwise comparisons: total
        aligned lengths; similarity errors in those alignments; the percentage
        of aligned length that matches (ANIm); and the percentage of the
        pairwise comparison that is aligned.

        For the total aligned length, similarity error, and ANIm dictionaries,
        as these are triangular/symmetrical matrices we only key them by
        (query, subject), but as the percentage aligned measure depends on the
        sequence we calculate it against, we report (query, subject) and
        (subject, query) values.

        - org_lengths is a dictionary of total sequence lengths for each
              input sequence
    """
    global CUM_RETVALS
    infiles = get_input_files(options.outdirname, '.delta')
    logger.info("Delta files:\n\t%s" % '\n\t'.join(infiles))
    logger.info("Processing .delta files")
    # We store pairwise comparison lengths in dictionaries, keyed by organism
    # ID pair tuples:
    # perc_aln is useful, as it is a matrix of the minimum percentage of an
    # organism's genome involved in a pairwise alignment
    lengths, sim_errors, perc_ids, perc_aln = {}, {}, {}, {}
    for dname in infiles:
        logger.info("Processing %s" % dname)
        qname, sname = \
            os.path.splitext(os.path.split(dname)[-1])[0].split('_vs_')
        logger.info("Query organism: %s; Subject organism: %s" %
                    (qname, sname))
        tot_length, tot_sim_error = parse_delta(dname)
        try:
            perc_id = 1 - 1. * tot_sim_error/tot_length
        except ZeroDivisionError:
            # If this is thrown, the proximate cause is an empty NUCmer output
            # The root cause may be a failed NUCmer run (when CUM_RETVALS>0)
            # or that one or more of the sequences is too distant for NUCmer
            # to identify a similarity.
            logger.error("One or more of the NUCmer output files contains " +
                         "no useable output.")
            if 0 < CUM_RETVALS:
                logger.error("One or more NUCmer runs failed. " +
                             "Please investigate.")
                logger.error("Please retry the NUCmer comparison for " +
                             "%s vs %s manually" % (qname, sname))
            else:
                logger.error("The NUCmer comparison between %s and %s " %
                             (qname, sname) +
                             "has no usable output. The comparison may be " +
                             "too distant for use. Consider using --maxmatch.")
            logger.error(last_exception())
        lengths[(qname, sname)] = tot_length
        sim_errors[(qname, sname)] = tot_sim_error
        perc_ids[(qname, sname)] = perc_id
        perc_aln[(qname, sname)] = 1.*tot_length/org_lengths[qname]
        perc_aln[(sname, qname)] = 1.*tot_length/org_lengths[sname]
    return lengths, sim_errors, perc_ids, perc_aln


# Parse BLAST tabular output and store total alignment length, similarity
# counts and percentage identity for each pairwise comparison
def process_blast(org_lengths):
    """ Read in the BLASTN comparison output files, and calculate alignment
        lengths, similarity errors, and percentage identity and alignment
        coverage for each input sequence comparison.

        - org_lengths is a dictionary of total sequence lengths for each
              input sequence.
    """
    infiles = get_input_files(options.outdirname, '.blast_tab')
    logger.info("BLAST files:\n\t%s" % '\n\t'.join(infiles))
    logger.info("Processing .blast_tab files")
    lengths, sim_errors, perc_ids, perc_aln = {}, {}, {}, {}
    org_names = set()
    for tname in infiles:
        logger.info("Processing %s" % tname)
        qname, sname = \
            os.path.splitext(os.path.split(tname)[-1])[0].split('_vs_')
        org_names.add(qname)
        org_names.add(sname)
        logger.info("Query organism: %s; Subject organism: %s" %
                    (qname, sname))
        tot_length, tot_sim_error = parse_blast(tname)
        if tot_length:
            perc_id = 1 - 1. * tot_sim_error/tot_length
        else:
            perc_id = 0
        lengths[(qname, sname)] = tot_length
        sim_errors[(qname, sname)] = tot_sim_error
        perc_ids[(qname, sname)] = perc_id
        perc_aln[(qname, sname)] = 1.*tot_length/org_lengths[qname]
        perc_aln[(sname, qname)] = 1.*tot_length/org_lengths[sname]
    return lengths, sim_errors, perc_ids, perc_aln


# Parse NUCmer delta file to get total alignment length and total sim_errors
def parse_delta(filename):
    """ Reads a NUCmer output .delta file, extracting the aligned length and
        number of similarity errors for each aligned uniquely-matched region,
        and returns the cumulative total for each as a tuple.

        - filename is the path to the input .delta file
    """
    aln_length, sim_errors = 0, 0
    for line in [l.strip().split() for l in open(filename, 'rU').readlines()]:
        if line[0] == 'NUCMER' or line[0].startswith('>'):  # Skip headers
            continue
        # We only want lines with seven columns:
        if len(line) == 7:
            aln_length += abs(int(line[1]) - int(line[0]))
            sim_errors += int(line[4])
    return aln_length, sim_errors


# Parse custom BLASTN output to get total alignment length and mismatches
def parse_blast(filename):
    """ Calculate the alignment length and total number of similarity errors
        for the passed BLASTN alignment file generated by comparing fragmented
        input sequences.

        - filename is the location of the BLASTN output for a pairwise
              comparison between input sequences
    """
    aln_length, sim_errors = 0, 0,
    qname, sname = \
        os.path.splitext(os.path.split(filename)[-1])[0].split('_vs_')
    qalnlen, qnumid, qlen, qerr = (collections.defaultdict(float),
                                   collections.defaultdict(float),
                                   collections.defaultdict(float),
                                   collections.defaultdict(float))
    seen = set()  # IDs of queries that have been processed
    for line in [l.strip().split() for l in open(filename, 'rU').readlines()
                 if len(l) and not l.startswith('#')]:
        # We need to collate matches by query ID, to determine whether the
        # match has > 30% identity and > 70% coverage.
        # Following Goris et al (2007) we only use matches that contribute to
        # a total match identity of at least 30% and a total match coverage
        # of at least 70% of either query or reference length
        # As of BLASTN 2.2.29+, the max_target_seqs/num_alignments still
        # will not restrict to the best hit. We rely on BLAST output giving
        # the best hit for any query first in the table, and check whether we
        # have seen it before. If not, we carry on processing the data.
        qid = line[0]
        if qid in seen:
            continue
        seen.add(qid)
        qalnlen[qid] += int(line[2])
        qnumid[qid] += int(line[5])
        qlen[qid] = int(line[6])
        qerr[qid] += int(line[3])
    for qid, ql in qlen.items():
        if 1.*qalnlen[qid]/ql > 0.7 and 1.*qnumid[qid]/ql > 0.3:
            aln_length += int(qalnlen[qid])
            sim_errors += int(qerr[qid])
    return aln_length, sim_errors


# Run BLASTN pairwise on input files, using multiprocessing
def pairwise_blast(filenames):
    """ Run BLASTN for each pairwise comparison of fragmented input sequences,
        using multiprocessing to take advantage of multiple cores where
        possible, and writing results to the nominated output directory.

        - filenames is an iterable of locations of input FASTA files, from
              which BLASTN command lines are constructed.

        We loop over all FASTA files in the input directory, generating
        BLASTN command line for each pairwise comparison, and then pass those
        command lines to be run using multiprocessing.
    """
    global CUM_RETVALS
    logger.info("Running pairwise BLASTN to generate *.blast_tab")
    cmdlines = []
    for idx, f1 in enumerate(filenames[:-1]):
        cmdlines.extend([make_blast_cmd(f1, f2)
                         for f2 in filenames[idx+1:]])
    logger.info("BLASTN command lines:\n\t%s" % '\n\t'.join(cmdlines))
    if not options.skip_blast:
        multiprocessing_run(cmdlines)
        if 0 < CUM_RETVALS:
            logger.error("At least one BLAST comparison returned a nonzero " +
                         "value. ANIb may fail.")
    else:
        logger.warning("BLASTN run skipped!")


# Run NUCmer pairwise on the input files, using multiprocessing
def pairwise_nucmer(filenames):
    """ Run NUCmer to generate pairwise alignment data for each of the
        input FASTA files.

        - filenames is a list of input FASTA filenames, from which NUCmer
              command lines are constructed

        We loop over all FASTA files in the input directory, generating NUCmer
        command lines for each pairwise comparison, and then pass those
        command lines to be run using multiprocessing.
    """
    global CUM_RETVALS
    logger.info("Running pairwise NUCmer comparison to generate *.delta")
    cmdlines = []
    for idx, f1 in enumerate(filenames[:-1]):
        cmdlines.extend([make_nucmer_cmd(f1, f2)
                         for f2 in filenames[idx+1:]])
    logger.info("NUCmer command lines:\n\t%s" % '\n\t'.join(cmdlines))
    if not options.skip_nucmer:
        multiprocessing_run(cmdlines)
        if 0 < CUM_RETVALS:
            logger.error("At least one NUCmer comparison failed. ANIm " +
                         "may fail.")
    else:
        logger.warning("NUCmer run skipped!")


# Run a set of command lines using multiprocessing
def multiprocessing_run(cmdlines):
    """ Distributes the passed command-line jobs using multiprocessing.

        - cmdlines is an iterable of command line strings
    """
    logger.info("Running %d jobs with multiprocessing" % len(cmdlines))
    pool = multiprocessing.Pool()
    completed = []
    if options.verbose:
        callback_fn = logger_callback
    else:
        callback_fn = completed.append
    pool_outputs = [pool.apply_async(subprocess.call,
                                     (str(cline), ),
                                     {'stderr': subprocess.PIPE,
                                      'shell': sys.platform != "win32"},
                                     callback=callback_fn)
                    for cline in cmdlines]
    pool.close()        # Run jobs
    pool.join()         # Collect output
    logger.info("Multiprocessing jobs completed:\n%s" % completed)


# Multiprocessing callback to logger
def logger_callback(val):
    """ Basic callback for multiprocessing just to log status of each job

        - val is an integer returned by multiprocessing, describing the run
             status
    """
    global CUM_RETVALS
    logger.info("Multiprocessing run completed with status: %s" % val)
    # Keep track of returned values, as these help diagnose problems
    # for ANIm analyses
    CUM_RETVALS += val


# Construct a command-line for NUCmer
def make_nucmer_cmd(f1, f2):
    """ Construct a command-line for NUCmer pairwise comparison, and return as
        a string

        - f1, f2 are the locations of two input FASTA files for analysis

        We use the -mum option so that we consider matches that are unique in
        both the reference and the query. -mumreference gives us matches
        unique only in the reference and -maxmatch gives us matches to all
        regions, regardless of uniqueness. We may want to make this an option.
    """
    prefix = os.path.join(options.outdirname, "%s_vs_%s" %
                          (os.path.splitext(os.path.split(f1)[-1])[0],
                           os.path.splitext(os.path.split(f2)[-1])[0]))
    # Do we use the --maxmatch option?
    if options.maxmatch:
        mode = "-maxmatch"
    else:
        mode = "-mum"
    cmd = "%s %s -p %s %s %s" % (options.nucmer_exe, mode, prefix, f1, f2)
    return cmd


# Construct a command-line for BLASTN
def make_blast_cmd(f1, f2):
    """ Construct a BLASTN command line to conduct sequence comparison between
        two fragmented input sequences, for ANIb.

        - f1, f2 are the locations of two input FASTA files for analysis
    """
    prefix = os.path.join(options.outdirname, "%s_vs_%s" %
                          (os.path.splitext(os.path.split(f1)[-1])[0],
                           os.path.splitext(os.path.split(f2)[-1])[0]))
    blastdb = os.path.splitext(f2)[0]
    cmd = "%s -out %s.blast_tab -query %s -db %s " % \
        (options.blast_exe, prefix, f1, blastdb) + \
        "-xdrop_gap_final 150 -dust no -evalue 1e-15 " +\
        "-max_target_seqs 1 -outfmt '6 qseqid sseqid length mismatch " +\
        "pident nident qlen slen qstart qend sstart send positive " +\
        "ppos gaps' -task blastn"
    return cmd


# Construct a command line for BLAST makeblastdb
def make_makeblastdb_cmd(filename):
    """ Construct a makeblastdb command line to make a BLAST nucleotide
        database from the passed fragmented input sequence FASTA file.

        - filename is the location of the fragmented input FASTA sequence
              file, for constructing the database
    """
    db_prefix = os.path.join(options.outdirname,
                             os.path.split(os.path.splitext(filename)[0])[-1])
    cmd = "%s -out %s -dbtype nucl -in %s" % (options.makeblastdb_exe,
                                              db_prefix, filename)
    return cmd


# Get list of FASTA files in a directory
def get_input_files(dir, *ext):
    """ Returns a list of files in the input directory with the passed
        extension

        - dir is the location of the directory containing the input files

        - *ext is a list of arguments describing permissible file extensions
    """
    filelist = [f for f in os.listdir(dir)
                if os.path.splitext(f)[-1] in ext]
    return [os.path.join(dir, f) for f in filelist]


# Create output directory if it doesn't exist
def make_outdir():
    """ Make the output directory, if required.

        This is a little involved.  If the output directory already exists,
        we take the safe option by default, and stop with an error.  We can,
        however, choose to force the program to go on, in which case we can
        either clobber or not the existing directory.  The options turn out
        as the following, if the directory exists:

        DEFAULT: stop
        FORCE: continue, and remove the existing output directory
        NOCLOBBER+FORCE: continue, but do not remove the existing output
    """
    if os.path.exists(options.outdirname):
        if not options.force:
            logger.error("Output directory %s would " % options.outdirname +
                         "overwrite existing files (exiting)")
            sys.exit(1)
        else:
            logger.info("Removing directory %s and everything below it" %
                        options.outdirname)
            if options.noclobber:
                logger.warning("NOCLOBBER: not actually deleting directory")
            else:
                shutil.rmtree(options.outdirname)
    logger.info("Creating directory %s" % options.outdirname)
    try:
        os.makedirs(options.outdirname)   # We make the directory recursively
    except OSError:
        # This gets thrown if the directory exists. If we've forced overwrite/
        # delete and we're not clobbering, we let things slide
        if options.noclobber and options.force:
            logger.info("NOCLOBBER+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            sys.exit(1)


#=============
# SCRIPT

if __name__ == '__main__':

    # Parse command-line
    # options are all options - no arguments
    options, args = parse_cmdline(sys.argv)

    # We set up logging, and modify loglevel according to whether we need
    # verbosity or not
    # err_handler points to sys.stderr
    # err_handler_file points to a logfile, if named
    logger = logging.getLogger('calculate_ani.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if options.logfile is not None:
        try:
            logstream = open(options.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging" %
                         options.logfile)
            sys.exit(1)
    if options.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)
    logger.info('# calculate_ani.py logfile')
    logger.info('# Run: %s' % time.asctime())

    # Report arguments, if verbose
    logger.info(options)
    logger.info(args)

    # Have we got an input and output directory? If not, exit.
    if options.indirname is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s" % options.indirname)
    if options.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    make_outdir()
    logger.info("Output directory: %s" % options.outdirname)

    # Have we got a valid method choice?
    methods = {"ANIm": calculate_anim,
               "ANIb": calculate_anib,
               "TETRA": calculate_tetra}
    if options.method not in methods:
        logger.error("ANI method %s not recognised (exiting)" % options.method)
        logger.error("Valid methods are: %s" % methods.keys())
        sys.exit(1)
    logger.info("Using ANI method: %s" % options.method)

    # Run method on the contents of the input directory, writing out
    # to the named output directory
    methods[options.method]()

    # If graphics have been selected, use R to generate a heatmap of the ANI
    # scores from the perc_id.tab output
    if options.graphics:
        logger.info("Rendering heatmap with R")
        if options.method == "TETRA":
            filename = "tetra_corr.tab"
        else:
            filename = "perc_ids.tab"
        if options.gformat.lower() not in ('pdf', 'png', 'jpg'):
            options.gformat = 'pdf'
        if not rpy2_import:
            logger.error("No rpy2 module: graphics are unavailable")
        else:
            rstr = ["library(gplots)",
                    "ani <- read.table('%s', " %
                    os.path.join(options.outdirname, filename) +
                    "header=T, sep='\\t', " +
                    "row.names=1)",
                    "ani[] <- lapply(ani, " +
                    "function(x){replace(x, x == 0, NA)})",
                    "%s('%s')" %
                    (options.gformat.lower(),
                     os.path.join(options.outdirname, '%s.%s' %
                                  (options.method,
                                   options.gformat.lower()))),
                    "heatmap.2(as.matrix(ani), col=bluered, " +
                    "breaks=seq(min(0.9, " +
                    "max(ani[!is.na(ani)]))-0.01,1,0.001), " +
                    "margins=c(15,12), " +
                    "cexCol=1/log10(ncol(ani)), " +
                    "cexRow=1/log10(nrow(ani)), main='%s')" % options.method,
                    "dev.off()"]
            rstr = '\n'.join(rstr)
            logger.info("R command:\n%s" % rstr)
            robjects.r(rstr)
