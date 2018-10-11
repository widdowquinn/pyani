# -*- coding: utf-8 -*-
"""Code to implement the ANIb average nucleotide identity method.

Calculates ANI by the ANIb method, as described in Goris et al. (2007)
Int J Syst Evol Micr 57: 81-91. doi:10.1099/ijs.0.64483-0.

From Goris et al.

'''The genomic sequence from one of the genomes in a pair (the query)
was cut into consecutive 1020 nt fragments. The 1020 nt cut-off was used
to correspond with the fragmentation of the genomic DNA to approximately
1 kb fragments during the DDH experiments. [...] The 1020 nt fragments
were then used to search against the whole genomic sequence of the other
genome in the pair (the reference) by using the BLASTN algorithm;
the best BLASTN match was saved for further analysis. The BLAST
algorithm was run using the following settings: X=150 (where X is the
drop-off value for gapped alignment), q=-1 (where q is the penalty
for nucleotide mismatch) and F=F (where F is the filter for repeated
sequences); the rest of the parameters were used at the default settings.
These settings give better sensitivity than the default settings when
more distantly related genomes are being compared, as the latter
target sequences that are more similar to each other.
[...]
The ANI between the query genome and the reference genome was
calculated as the mean identity of all BLASTN matches that showed more
than 30% overall sequence identity (recalculated to an identity along
the entire sequence) over an alignable region of at least 70% of their
length. This cut-off is above the 'twilight zone' of similarity searches in
which an inference of homology is error prone because of low levels of
Reverse searching, i.e. in which the reference genome is used as the
query, was also performed to provide reciprocal values.'''

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

(c) The James Hutton Institute 2013-2018
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2013-2018 The James Hutton Institute

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
"""

import os
import shutil

import pandas as pd

from Bio import SeqIO

from . import pyani_config
from . import pyani_files
from . import pyani_jobs
from .pyani_tools import ANIResults, BLASTcmds, BLASTexes, BLASTfunctions


# Divide input FASTA sequences into fragments
def fragment_fasta_files(infiles, outdirname, fragsize):
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
        outstem, outext = os.path.splitext(os.path.split(fname)[-1])
        outfname = os.path.join(outdirname, outstem) + "-fragments" + outext
        outseqs = []
        count = 0
        for seq in SeqIO.parse(fname, "fasta"):
            idx = 0
            while idx < len(seq):
                count += 1
                newseq = seq[idx : idx + fragsize]
                newseq.id = "frag%05d" % count
                outseqs.append(newseq)
                idx += fragsize
        outfnames.append(outfname)
        SeqIO.write(outseqs, outfname, "fasta")
    return outfnames, get_fraglength_dict(outfnames)


# Get lengths of all sequences in all files
def get_fraglength_dict(fastafiles):
    """Returns dictionary of sequence fragment lengths, keyed by query name.

    - fastafiles - list of FASTA input whole sequence files

    Loops over input files and, for each, produces a dictionary with fragment
    lengths, keyed by sequence ID. These are returned as a dictionary with
    the keys being query IDs derived from filenames.
    """
    fraglength_dict = {}
    for filename in fastafiles:
        qname = os.path.split(filename)[-1].split("-fragments")[0]
        fraglength_dict[qname] = get_fragment_lengths(filename)
    return fraglength_dict


# Get lengths of all sequences in a file
def get_fragment_lengths(fastafile):
    """Returns dictionary of sequence fragment lengths, keyed by fragment ID.

    Biopython's SeqIO module is used to parse all sequences in the FASTA
    file.

    NOTE: ambiguity symbols are not discounted.
    """
    fraglengths = {}
    for seq in SeqIO.parse(fastafile, "fasta"):
        fraglengths[seq.id] = len(seq)
    return fraglengths


# Create dictionary of database building commands, keyed by dbname
def build_db_jobs(infiles, blastcmds):
    """Returns dictionary of db-building commands, keyed by dbname."""
    dbjobdict = {}  # Dict of database construction jobs, keyed by filename
    # Create dictionary of database building jobs, keyed by db name
    # defining jobnum for later use as last job index used
    for idx, fname in enumerate(infiles):
        dbjobdict[blastcmds.get_db_name(fname)] = pyani_jobs.Job(
            "%s_db_%06d" % (blastcmds.prefix, idx), blastcmds.build_db_cmd(fname)
        )
    return dbjobdict


def make_blastcmd_builder(
    mode, outdir, format_exe=None, blast_exe=None, prefix="ANIBLAST"
):
    """Returns BLASTcmds object for construction of BLAST commands."""
    if mode == "ANIb":  # BLAST/formatting executable depends on mode
        blastcmds = BLASTcmds(
            BLASTfunctions(construct_makeblastdb_cmd, construct_blastn_cmdline),
            BLASTexes(
                format_exe or pyani_config.MAKEBLASTDB_DEFAULT,
                blast_exe or pyani_config.BLASTN_DEFAULT,
            ),
            prefix,
            outdir,
        )
    else:
        blastcmds = BLASTcmds(
            BLASTfunctions(construct_formatdb_cmd, construct_blastall_cmdline),
            BLASTexes(
                format_exe or pyani_config.FORMATDB_DEFAULT,
                blast_exe or pyani_config.BLASTALL_DEFAULT,
            ),
            prefix,
            outdir,
        )
    return blastcmds


# Make a dependency graph of BLAST commands
def make_job_graph(infiles, fragfiles, blastcmds):
    """Return a job dependency graph, based on the passed input sequence files.

    - infiles - a list of paths to input FASTA files
    - fragfiles - a list of paths to fragmented input FASTA files

    By default, will run ANIb - it *is* possible to make a mess of passing the
    wrong executable for the mode you're using.

    All items in the returned graph list are BLAST executable jobs that must
    be run *after* the corresponding database creation. The Job objects
    corresponding to the database creation are contained as dependencies.
    How those jobs are scheduled depends on the scheduler (see
    run_multiprocessing.py, run_sge.py)
    """
    joblist = []  # Holds list of job dependency graphs

    # Get dictionary of database-building jobs
    dbjobdict = build_db_jobs(infiles, blastcmds)

    # Create list of BLAST executable jobs, with dependencies
    jobnum = len(dbjobdict)
    for idx, fname1 in enumerate(fragfiles[:-1]):
        for fname2 in fragfiles[idx + 1 :]:
            jobnum += 1
            jobs = [
                pyani_jobs.Job(
                    "%s_exe_%06d_a" % (blastcmds.prefix, jobnum),
                    blastcmds.build_blast_cmd(fname1, fname2.replace("-fragments", "")),
                ),
                pyani_jobs.Job(
                    "%s_exe_%06d_b" % (blastcmds.prefix, jobnum),
                    blastcmds.build_blast_cmd(fname2, fname1.replace("-fragments", "")),
                ),
            ]
            jobs[0].add_dependency(dbjobdict[fname1.replace("-fragments", "")])
            jobs[1].add_dependency(dbjobdict[fname2.replace("-fragments", "")])
            joblist.extend(jobs)

    # Return the dependency graph
    return joblist


# Generate list of makeblastdb command lines from passed filenames
def generate_blastdb_commands(filenames, outdir, blastdb_exe=None, mode="ANIb"):
    """Return a list of makeblastdb command-lines for ANIb/ANIblastall

    - filenames - a list of paths to input FASTA files
    - outdir - path to output directory
    - blastdb_exe - path to the makeblastdb executable
    """
    if mode == "ANIb":
        construct_db_cmdline = construct_makeblastdb_cmd
    else:
        construct_db_cmdline = construct_formatdb_cmd
    if blastdb_exe is None:
        cmdlines = [construct_db_cmdline(fname, outdir) for fname in filenames]
    else:
        cmdlines = [
            construct_db_cmdline(fname, outdir, blastdb_exe) for fname in filenames
        ]
    return cmdlines


# Generate single makeblastdb command line
def construct_makeblastdb_cmd(
    filename, outdir, blastdb_exe=pyani_config.MAKEBLASTDB_DEFAULT
):
    """Returns a single makeblastdb command.

    - filename - input filename
    - blastdb_exe - path to the makeblastdb executable
    """
    title = os.path.splitext(os.path.split(filename)[-1])[0]
    outfilename = os.path.join(outdir, os.path.split(filename)[-1])
    return (
        "{0} -dbtype nucl -in {1} -title {2} -out {3}".format(
            blastdb_exe, filename, title, outfilename
        ),
        outfilename,
    )


# Generate single makeblastdb command line
def construct_formatdb_cmd(filename, outdir, blastdb_exe=pyani_config.FORMATDB_DEFAULT):
    """Returns a single formatdb command.

    - filename - input filename
    - blastdb_exe - path to the formatdb executable
    """
    title = os.path.splitext(os.path.split(filename)[-1])[0]
    newfilename = os.path.join(outdir, os.path.split(filename)[-1])
    shutil.copy(filename, newfilename)
    return (
        "{0} -p F -i {1} -t {2}".format(blastdb_exe, newfilename, title),
        newfilename,
    )


# Generate list of BLASTN command lines from passed filenames
def generate_blastn_commands(filenames, outdir, blast_exe=None, mode="ANIb"):
    """Return a list of blastn command-lines for ANIm

    - filenames - a list of paths to fragmented input FASTA files
    - outdir - path to output directory
    - blastn_exe - path to BLASTN executable

    Assumes that the fragment sequence input filenames have the form
    ACCESSION-fragments.ext, where the corresponding BLAST database filenames
    have the form ACCESSION.ext. This is the convention followed by the
    fragment_FASTA_files() function above.
    """
    if mode == "ANIb":
        construct_blast_cmdline = construct_blastn_cmdline
    else:
        construct_blast_cmdline = construct_blastall_cmdline
    cmdlines = []
    for idx, fname1 in enumerate(filenames[:-1]):
        dbname1 = fname1.replace("-fragments", "")
        for fname2 in filenames[idx + 1 :]:
            dbname2 = fname2.replace("-fragments", "")
            if blast_exe is None:
                cmdlines.append(construct_blast_cmdline(fname1, dbname2, outdir))
                cmdlines.append(construct_blast_cmdline(fname2, dbname1, outdir))
            else:
                cmdlines.append(
                    construct_blast_cmdline(fname1, dbname2, outdir, blast_exe)
                )
                cmdlines.append(
                    construct_blast_cmdline(fname2, dbname1, outdir, blast_exe)
                )
    return cmdlines


# Generate single BLASTN command line
def construct_blastn_cmdline(
    fname1, fname2, outdir, blastn_exe=pyani_config.BLASTN_DEFAULT
):
    """Returns a single blastn command.

    - filename - input filename
    - blastn_exe - path to BLASTN executable
    """
    fstem1 = os.path.splitext(os.path.split(fname1)[-1])[0]
    fstem2 = os.path.splitext(os.path.split(fname2)[-1])[0]
    fstem1 = fstem1.replace("-fragments", "")
    prefix = os.path.join(outdir, "%s_vs_%s" % (fstem1, fstem2))
    cmd = (
        "{0} -out {1}.blast_tab -query {2} -db {3} "
        + "-xdrop_gap_final 150 -dust no -evalue 1e-15 "
        + "-max_target_seqs 1 -outfmt '6 qseqid sseqid length mismatch "
        + "pident nident qlen slen qstart qend sstart send positive "
        + "ppos gaps' -task blastn"
    )
    return cmd.format(blastn_exe, prefix, fname1, fname2)


# Generate single BLASTALL command line
def construct_blastall_cmdline(
    fname1, fname2, outdir, blastall_exe=pyani_config.BLASTALL_DEFAULT
):
    """Returns a single blastall command.

    - blastall_exe - path to BLASTALL executable
    """
    fstem1 = os.path.splitext(os.path.split(fname1)[-1])[0]
    fstem2 = os.path.splitext(os.path.split(fname2)[-1])[0]
    fstem1 = fstem1.replace("-fragments", "")
    prefix = os.path.join(outdir, "%s_vs_%s" % (fstem1, fstem2))
    cmd = (
        "{0} -p blastn -o {1}.blast_tab -i {2} -d {3} "
        + "-X 150 -q -1 -F F -e 1e-15 "
        + "-b 1 -v 1 -m 8"
    )
    return cmd.format(blastall_exe, prefix, fname1, fname2)


# Process pairwise BLASTN output
def process_blast(
    blast_dir,
    org_lengths,
    fraglengths=None,
    mode="ANIb",
    identity=0.3,
    coverage=0.7,
    logger=None,
):
    """Returns a tuple of ANIb results for .blast_tab files in the output dir.

    - blast_dir - path to the directory containing .blast_tab files
    - org_lengths - the base count for each input sequence
    - fraglengths - dictionary of query sequence fragment lengths, only
    needed for BLASTALL output
    - mode - parsing BLASTN+ or BLASTALL output?
    - logger - a logger for messages

    Returns the following pandas dataframes in an ANIResults object;
    query sequences are rows, subject sequences are columns:

    - alignment_lengths - non-symmetrical: total length of alignment
    - percentage_identity - non-symmetrical: ANIb (Goris) percentage identity
    - alignment_coverage - non-symmetrical: coverage of query
    - similarity_errors - non-symmetrical: count of similarity errors

    May throw a ZeroDivisionError if one or more BLAST runs failed, or a
    very distant sequence was included in the analysis.
    """
    # Process directory to identify input files
    blastfiles = pyani_files.get_input_files(blast_dir, ".blast_tab")
    # Hold data in ANIResults object
    results = ANIResults(list(org_lengths.keys()), mode)

    # Fill diagonal NA values for alignment_length with org_lengths
    for org, length in list(org_lengths.items()):
        results.alignment_lengths[org][org] = length

    # Process .blast_tab files assuming that the filename format holds:
    # org1_vs_org2.blast_tab:
    for blastfile in blastfiles:
        qname, sname = os.path.splitext(os.path.split(blastfile)[-1])[0].split("_vs_")

        # We may have BLAST files from other analyses in the same directory
        # If this occurs, we raise a warning, and skip the file
        if qname not in list(org_lengths.keys()):
            if logger:
                logger.warning(
                    "Query name %s not in input " % qname
                    + "sequence list, skipping %s" % blastfile
                )
            continue
        if sname not in list(org_lengths.keys()):
            if logger:
                logger.warning(
                    "Subject name %s not in input " % sname
                    + "sequence list, skipping %s" % blastfile
                )
            continue
        resultvals = parse_blast_tab(blastfile, fraglengths, identity, coverage, mode)
        query_cover = float(resultvals[0]) / org_lengths[qname]

        # Populate dataframes: when assigning data, we need to note that
        # we have asymmetrical data from BLAST output, so only the
        # upper triangle is populated
        results.add_tot_length(qname, sname, resultvals[0], sym=False)
        results.add_sim_errors(qname, sname, resultvals[1], sym=False)
        results.add_pid(qname, sname, 0.01 * resultvals[2], sym=False)
        results.add_coverage(qname, sname, query_cover)
    return results


# Parse BLASTALL output to get total alignment length and mismatches
def parse_blast_tab(filename, fraglengths, identity, coverage, mode="ANIb"):
    """Returns (alignment length, similarity errors, mean_pid) tuple
    from .blast_tab

    - filename - path to .blast_tab file

    Calculate the alignment length and total number of similarity errors (as
    we would with ANIm), as well as the Goris et al.-defined mean identity
    of all valid BLAST matches for the passed BLASTALL alignment .blast_tab
    file.

    '''ANI between the query genome and the reference genome was calculated as
    the mean identity of all BLASTN matches that showed more than 30% overall
    sequence identity (recalculated to an identity along the entire sequence)
    over an alignable region of at least 70% of their length.
    '''
    """
    # Assuming that the filename format holds org1_vs_org2.blast_tab:
    qname = os.path.splitext(os.path.split(filename)[-1])[0].split("_vs_")[0]
    # Load output as dataframe
    if mode == "ANIblastall":
        qfraglengths = fraglengths[qname]
        columns = [
            "sid",
            "blast_pid",
            "blast_alnlen",
            "blast_mismatch",
            "blast_gaps",
            "q_start",
            "q_end",
            "s_start",
            "s_end",
            "e_Value",
            "bit_score",
        ]
    else:
        columns = [
            "sbjct_id",
            "blast_alnlen",
            "blast_mismatch",
            "blast_pid",
            "blast_identities",
            "qlen",
            "slen",
            "q_start",
            "q_end",
            "s_start",
            "s_end",
            "blast_pos",
            "ppos",
            "blast_gaps",
        ]
    # We may receive an empty BLASTN output file, if there are no significant
    # regions of homology. This causes pandas to throw an error on CSV import.
    # To get past this, we create an empty dataframe with the appropriate
    # columns.
    try:
        data = pd.read_csv(filename, header=None, sep="\t", index_col=0)
        data.columns = columns
    except pd.io.common.EmptyDataError:
        data = pd.DataFrame(columns=columns)
    # Add new column for fragment length, only for BLASTALL
    if mode == "ANIblastall":
        data["qlen"] = pd.Series(
            [qfraglengths[idx] for idx in data.index], index=data.index
        )
    # Add new columns for recalculated alignment length, proportion, and
    # percentage identity
    data["ani_alnlen"] = data["blast_alnlen"] - data["blast_gaps"]
    data["ani_alnids"] = data["ani_alnlen"] - data["blast_mismatch"]
    data["ani_coverage"] = data["ani_alnlen"] / data["qlen"]
    data["ani_pid"] = data["ani_alnids"] / data["qlen"]
    # Filter rows on 'ani_coverage' > 0.7, 'ani_pid' > 0.3
    filtered = data[(data["ani_coverage"] > coverage) & (data["ani_pid"] > identity)]
    # Dedupe query hits, so we only take the best hit
    filtered = filtered.groupby(filtered.index).first()
    # Replace NaNs with zero
    filtered = filtered.fillna(value=0)  # Needed if no matches
    # The ANI value is then the mean percentage identity.
    # We report total alignment length and the number of similarity errors
    # (mismatches and gaps), as for ANIm
    # NOTE: We report the mean of 'blast_pid' for concordance with JSpecies
    # Despite this, the concordance is not exact. Manual inspection during
    # development indicated that a handful of fragments are differentially
    # filtered out in JSpecies and this script. This is often on the basis
    # of rounding differences (e.g. coverage being close to 70%).
    # NOTE: If there are no hits, then ani_pid will be nan - we replace this
    # with zero if that happens
    ani_pid = filtered["blast_pid"].mean()
    if pd.isnull(ani_pid):  # Happens if there are no matches in ANIb
        ani_pid = 0
    aln_length = filtered["ani_alnlen"].sum()
    sim_errors = filtered["blast_mismatch"].sum() + filtered["blast_gaps"].sum()
    filtered.to_csv(filename + ".dataframe", sep="\t")
    return aln_length, sim_errors, ani_pid
