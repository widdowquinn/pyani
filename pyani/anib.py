# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence. 
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to implement the ANIb average nucleotide identity method.

Calculates ANI by the ANIb method, as described in Goris et al. (2007)
Int J Syst Evol Micr 57: 81-91. doi:10.1099/ijs.0.64483-0.

All input FASTA format files are used to construct BLAST databases.
Each file's contents are also split into sequence fragments of length
options.fragsize, and the multiple FASTA file that results written to
the output directory. These are BLASTNed, pairwise, against the
databases.

BLAST output is interrogated for all fragment matches that cover
at least 70% of the query sequence, with at least 30% nucleotide
identity over the full length of the query sequence. This is an odd
choice and doesn't correspond to the twilight zone limit as implied by
Goris et al. We persist with their definition, however.  Only these
qualifying matches contribute to the total aligned length, and total
aligned sequence identity used to calculate ANI.
"""

import pandas as pd

import collections
import os

import pyani_config, pyani_files

from Bio import SeqIO


# Divide input FASTA sequences into fragments
def fragment_FASTA_files(infiles, outdirname, fragsize):
    """Chops sequences of the passed files into fragments, returns filenames.

    - infiles - paths to each input sequence file
    - outdirname - path to output directory
    - fragsize - the size of sequence fragments

    Takes every sequence from every file in infiles, and splits them into
    consecutive fragments of length fragsize, (with any trailing sequences
    being included, even if shorter than fragsize), and writes the resulting
    set of sequences to a file with the same name in the output directory.
    All fragments are named consecutively and uniquely (within a file) as
    fragNNNNN. Sequence description fields are retained.
    """
    outfnames = []
    for fname in infiles:
        outstem = os.path.splitext(os.path.split(fname)[-1])[0]
        outfname = os.path.join(outdirname, outstem) + '.fasta'
        outseqs = []
        idx, count = 0, 0
        for seq in SeqIO.parse(fname, 'fasta'):
            while idx < len(seq):
                count += 1
                newseq = seq[idx:idx+fragsize]
                newseq.id = "frag%05d" % count
                outseqs.append(newseq)
                idx += fragsize
        outfnames.append(outfname)
        SeqIO.write(outseqs, outfname, 'fasta')
    return outfnames


# Generate list of makeblastdb command lines from passed filenames
def generate_blastdb_commands(filenames, outdir,
                              blastdb_exe=pyani_config.MAKEBLASTDB_DEFAULT):
    """Return a list of makeblastdb command-lines for ANIm

    - filenames - a list of paths to input FASTA files
    - outdir - path to output directory
    - blastdb_exe - path to the makeblastdb executable
    """
    cmdlines = [construct_makeblastdb_cmdline(fname, blastdb_exe) for 
                fname in filenames]
    return cmdlines


# Generate single makeblastdb command line
def construct_makeblastdb_cmdline(filename,
                                  blastdb_exe=pyani_config.MAKEBLASTDB_DEFAULT):
    """Returns a single makeblastdb command.

    - filename - input filename
    - blastdb_exe - path to the makeblastdb executable
    """
    title = os.path.splitext(os.path.split(filename)[-1])[0] + '-fragments'
    return "{0} -dbtype nucl -in {1} -title {2}".format(blastdb_exe,
                                                        filename, title)

# Generate list of BLASTN command lines from passed filenames
def generate_blastn_commands(filenames, outdir,
                             blastn_exe=pyani_config.BLASTN_DEFAULT):
    """Return a list of makeblastdb command-lines for ANIm

    - filenames - a list of paths to input FASTA files
    - outdir - path to output directory
    - blastn_exe - path to BLASTN executable
    """
    cmdlines = []
    for idx, fname1 in enumerate(filenames[:-1]):
        cmdlines.extend([construct_makeblastn_cmdline(fname1, fname2, outdir,
                                                      blastn_exe) for 
                         fname2 in filenames[idx+1:]])
    return cmdlines


# Generate single BLASTN command line
def construct_makeblastn_cmdline(fname1, fname2, outdir,
                                 blastn_exe=pyani_config.BLASTN_DEFAULT):
    """Returns a single makeblastdb command.

    - filename - input filename
    - blastn_exe - path to BLASTN executable
    """
    fstem1 = os.path.splitext(os.path.split(fname1)[-1])[0]
    fstem2 = os.path.splitext(os.path.split(fname2)[-1])[0]
    prefix = os.path.join(outdir, "%s_vs_%s" % (fstem1, fstem2))
    cmd = "{0} -out {1}.blast_tab -query {2} -db {3} " +\
        "-xdrop_gap_final 150 -dust no -evalue 1e-15 " +\
        "-max_target_seqs 1 -outfmt '6 qseqid sseqid length mismatch " +\
        "pident nident qlen slen qstart qend sstart send positive " +\
        "ppos gaps' -task blastn"
    return cmd.format(blastn_exe, prefix, fname1, fname2) 


# Process pairwise BLASTN output
def process_blastn(blast_dir, org_lengths):
    """Returns a tuple of ANIb results for .blast_tab files in the output dir.

    - blast_dir - path to the directory containing .blast_tab files
    - org_lengths - the base count for each input sequence
    
    Returns the following pandas dataframes in a tuple:

    - alignment_lengths - symmetrical: total length of alignment
    - percentage_identity - symmetrical: percentage identity of alignment
    - alignment_coverage - non-symmetrical: coverage of query and subject
    - similarity_errors - symmetrical: count of similarity errors

    May throw a ZeroDivisionError if one or more BLAST runs failed, or a
    very distant sequence was included in the analysis.
    """
    # Process directory to identify input files
    blastfiles = pyani_files.get_input_files(blast_dir, '.blast_tab')
    labels = org_lengths.keys()
    # Hold data in pandas dataframe
    alignment_lengths = pd.DataFrame(index=labels, columns=labels,
                                     dtype=float)
    similarity_errors = pd.DataFrame(index=labels, columns=labels,
                                     dtype=float).fillna(0)
    percentage_identity = pd.DataFrame(index=labels, columns=labels,
                                       dtype=float).fillna(1.0)
    alignment_coverage = pd.DataFrame(index=labels, columns=labels,
                                      dtype=float).fillna(1.0)
    # Fill diagonal NA values for alignment_length with org_lengths
    for org, length in org_lengths.items():
        alignment_lengths[org][org] = length
    # Process .blast_tab files assuming that the filename format holds:
    # org1_vs_org2.blast_tab:
    for blastfile in blastfiles:
        qname, sname = \
            os.path.splitext(os.path.split(blastfile)[-1])[0].split('_vs_')
        tot_length, tot_sim_error = parse_blast(blastfile)
        query_cover = float(tot_length) / org_lengths[qname]
        sbjct_cover = float(tot_length) / org_lengths[sname]
        # Calculate percentage ID of aligned length. This may fail if 
        # total length is zero.
        # The ZeroDivisionError that would arise should be handled 
        # Common causes are that a BLASTN run failed, or that a very 
        # distant sequence was included in the analysis.
        perc_id = 1 - float(tot_sim_error) / tot_length
        # Populate dataframes: when assigning data, pandas dataframes 
        # take column, index order, i.e. df['column']['row'] - this only
        # matters for asymmetrical data
        alignment_lengths[qname][sname] = tot_length
        alignment_lengths[sname][qname] = tot_length
        similarity_errors[qname][sname] = tot_sim_error
        similarity_errors[sname][qname] = tot_sim_error
        percentage_identity[qname][sname] = perc_id
        percentage_identity[sname][qname] = perc_id
        alignment_coverage[sname][qname] = query_cover
        alignment_coverage[qname][sname] = sbjct_cover
    return(alignment_lengths, percentage_identity, alignment_coverage,
           similarity_errors)
        
# Parse BLASTN output to get total alignment length and mismatches
def parse_blast(filename):
    """Returns (alignment length, similarity errors) tuple from .blast_tab

    - filename - path to .blast_tab file

    Calculate the alignment length and total number of similarity errors
    for the passed BLASTN alignment .blast_tab file.
    """
    aln_length, sim_errors = 0, 0
    # Assuming that the filename format holds org1_vs_org2.blast_tab:
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
