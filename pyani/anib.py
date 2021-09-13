# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019-2020
# Author: Leighton Pritchard
#
# Contact: leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016-2019 The James Hutton Institute
# Copyright (c) 2019-2020 University of Strathclyde
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
"""

import platform
import re
import shutil
import subprocess

from logging import Logger
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import pandas as pd  # type: ignore

from Bio import SeqIO  # type: ignore

from . import pyani_config
from . import pyani_files
from . import pyani_jobs
from .pyani_tools import ANIResults, BLASTcmds, BLASTexes, BLASTfunctions


def get_version(blast_exe: Path = pyani_config.BLASTN_DEFAULT) -> str:
    """Return BLAST+ blastn version as a string.

    :param blast_exe:  path to blastn executable

    We expect blastn to return a string as, for example

    .. code-block:: bash

        $ blastn -version
        blastn: 2.9.0+
        Package: blast 2.9.0, build Jun 10 2019 09:40:53

    This is concatenated with the OS name.
    """
    cmdline = [blast_exe, "-version"]
    result = subprocess.run(
        cmdline,  # type: ignore
        shell=False,
        stdout=subprocess.PIPE,  # type: ignore
        stderr=subprocess.PIPE,
        check=True,
    )
    version = re.search(  # type: ignore
        r"(?<=blastn:\s)[0-9\.]*\+", str(result.stdout, "utf-8")
    ).group()
    return f"{platform.system()}_{version}"


# Divide input FASTA sequences into fragments
def fragment_fasta_files(
    infiles: List[Path], outdirname: Path, fragsize: int
) -> Tuple[List, Dict]:
    """Chop sequences of the passed files into fragments, return filenames.

    :param infiles:  collection of paths to each input sequence file
    :param outdirname:  Path, path to output directory
    :param fragsize:  Int, the size of sequence fragments

    Takes every sequence from every file in infiles, and splits them into
    consecutive fragments of length fragsize, (with any trailing sequences
    being included, even if shorter than fragsize), writing the resulting
    set of sequences to a file with the same name in the specified
    output directory.

    All fragments are named consecutively and uniquely (within a file) as
    fragNNNNN. Sequence description fields are retained.

    Returns a tuple ``(filenames, fragment_lengths)`` where ``filenames`` is a
    list of paths to the fragment sequence files, and ``fragment_lengths`` is
    a dictionary of sequence fragment lengths, keyed by the sequence files,
    with values being a dictionary of fragment lengths, keyed by fragment
    IDs.
    """
    outfnames = []
    for fname in infiles:
        outfname = outdirname / f"{fname.stem}-fragments{fname.suffix}"
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
def get_fraglength_dict(fastafiles: List[Path]) -> Dict:
    """Return dictionary of sequence fragment lengths, keyed by query name.

    :param fastafiles:  list of paths to FASTA input whole sequence files

    Loops over input files and, for each, produces a dictionary with fragment
    lengths, keyed by sequence ID. These are returned as a dictionary with
    the keys being query IDs derived from filenames.
    """
    fraglength_dict = {}
    for filename in fastafiles:
        qname = filename.stem.split("-fragments")[0]
        fraglength_dict[qname] = get_fragment_lengths(filename)
    return fraglength_dict


# Get lengths of all sequences in a file
def get_fragment_lengths(fastafile: Path) -> Dict:
    """Return dictionary of sequence fragment lengths, keyed by fragment ID.

    :param fastafile:

    Biopython's SeqIO module is used to parse all sequences in the FASTA
    file.

    NOTE: ambiguity symbols are not discounted.
    """
    fraglengths = {}
    for seq in SeqIO.parse(fastafile, "fasta"):
        fraglengths[seq.id] = len(seq)
    return fraglengths


# Create dictionary of database building commands, keyed by dbname
def build_db_jobs(infiles: List[Path], blastcmds: BLASTcmds) -> Dict:
    """Return dictionary of db-building commands, keyed by dbname.

    :param infiles:
    :param blastcmds:
    """
    dbjobdict = {}  # Dict of database construction jobs, keyed by filename
    # Create dictionary of database building jobs, keyed by db name
    # defining jobnum for later use as last job index used
    for idx, fname in enumerate(infiles):
        dbjobdict[blastcmds.get_db_name(fname)] = pyani_jobs.Job(
            f"{blastcmds.prefix}_db_{idx:06}", blastcmds.build_db_cmd(fname)
        )
    return dbjobdict


def make_blastcmd_builder(
    method: str,
    outdir: Path,
    format_exe: Optional[Path] = None,
    blast_exe: Optional[Path] = None,
    prefix: str = "ANIBLAST",
) -> BLASTcmds:
    """Return BLASTcmds object for construction of BLAST commands.

    :param method:  str, the kind of ANI analysis (ANIb)
    :param outdir:
    :param format_exe:
    :param blast_exe:
    :param prefix:
    """
    blastcmds = BLASTcmds(
        BLASTfunctions(construct_makeblastdb_cmd, construct_blastn_cmdline),
        BLASTexes(
            format_exe or pyani_config.MAKEBLASTDB_DEFAULT,
            blast_exe or pyani_config.BLASTN_DEFAULT,
        ),
        prefix,
        outdir,
    )
    return blastcmds


# Make a dependency graph of BLAST commands
def make_job_graph(
    infiles: List[Path], fragfiles: List[Path], blastcmds: BLASTcmds
) -> List[pyani_jobs.Job]:
    """Return job dependency graph, based on the passed input sequence files.

    :param infiles:  list of paths to input FASTA files
    :param fragfiles:  list of paths to fragmented input FASTA files
    :param blastcmds:

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
                    f"{blastcmds.prefix}_exe_{jobnum:06d}_a",
                    blastcmds.build_blast_cmd(
                        fname1, fname2.parent / fname2.name.replace("-fragments", "")
                    ),
                ),
                pyani_jobs.Job(
                    f"{blastcmds.prefix}_exe_{jobnum:06d}_b",
                    blastcmds.build_blast_cmd(
                        fname2, fname1.parent / fname1.name.replace("-fragments", "")
                    ),
                ),
            ]
            jobs[0].add_dependency(
                dbjobdict[fname1.parent / fname1.name.replace("-fragments", "")]
            )
            jobs[1].add_dependency(
                dbjobdict[fname2.parent / fname2.name.replace("-fragments", "")]
            )
            joblist.extend(jobs)

    # Return the dependency graph
    return joblist


# Generate list of makeblastdb command lines from passed filenames
def generate_blastdb_commands(
    filenames: List[Path],
    outdir: Path,
    blastdb_exe: Optional[Path] = None,
) -> List[Tuple[str, Path]]:
    """Return list of makeblastdb command-lines for ANIb.

    :param filenames:  a list of paths to input FASTA files
    :param outdir:  path to output directory
    :param blastdb_exe:  path to the makeblastdb executable
    :param method:  str, ANIb analysis type (ANIb)
    """
    if blastdb_exe is None:
        cmdlines = [construct_makeblastdb_cmd(fname, outdir) for fname in filenames]
    else:
        cmdlines = [
            construct_makeblastdb_cmd(fname, outdir, blastdb_exe) for fname in filenames
        ]
    return cmdlines


# Generate single makeblastdb command line
def construct_makeblastdb_cmd(
    filename: Path, outdir: Path, blastdb_exe: Path = pyani_config.MAKEBLASTDB_DEFAULT
) -> Tuple[str, Path]:
    """Return makeblastdb command and path to output file.

    :param filename:  Path, input filename
    :param outdir:  Path, directory for output
    :param blastdb_exe:  Path, path to the makeblastdb executable
    """
    outfilename = outdir / filename.name.replace("-fragments", "")
    return (
        f"{blastdb_exe} -dbtype nucl -in {filename} -title {filename.stem} -out {outfilename}",
        outfilename,
    )


# Generate list of BLASTN command lines from passed filenames
def generate_blastn_commands(
    query: Path,
    subject: Path,
    outdir: Path,
    blast_exe: Optional[Path] = None,
) -> List[str]:
    """Return a list of blastn command-lines for ANIb.

    :param query:  a paths to the query's fragmented input FASTA file
    :param subject:  a paths to the subject's fragmented input FASTA file
    :param outdir:  path to output directory
    :param blastn_exe:  path to BLASTN executable
    :param method:  str, analysis type (ANIb)

    Assumes that the fragment sequence input filenames have the form
    ACCESSION-fragments.ext, where the corresponding BLAST database filenames
    have the form ACCESSION.ext. This is the convention followed by the
    fragment_FASTA_files() function above.
    """
    subj_db = outdir / "blastdbs" / Path(str(subject.name))
    if blast_exe is None:
        cmdline = construct_blastn_cmdline(query, subj_db, outdir)
    else:
        cmdline = construct_blastn_cmdline(query, subj_db, outdir, blast_exe)
    return cmdline


# Generate single BLASTN command line
def construct_blastn_cmdline(
    query: Path,
    subj_db: Path,
    outdir: Path,
    blastn_exe: Path = pyani_config.BLASTN_DEFAULT,
) -> str:
    """Return a single blastn command.

    :param fname1:  Path, FASTA file for query genome
    :param fname2:  Path, database of fragments for subject genome
    :param outdir:  Path, to the output directory
    :param blastn_exe:  str, path to blastn executable
    """
    prefix = outdir / f"{query.stem.replace('-fragments', '')}_vs_{subj_db.stem}"
    return (
        f"{blastn_exe} -out {prefix}.blast_tab -query {query} -db {subj_db} "
        "-xdrop_gap_final 150 -dust no -evalue 1e-15 -max_target_seqs 1 -outfmt "
        "'6 qseqid sseqid length mismatch pident nident qlen slen "
        "qstart qend sstart send positive ppos gaps' "
        "-task blastn"
    )


# Process pairwise BLASTN output
def process_blast(
    blast_dir: Path,
    org_lengths: Dict,
    fraglengths: Dict,
    method: str = "ANIb",
    logger: Optional[Logger] = None,
) -> ANIResults:
    """Return tuple of ANIb results for .blast_tab files in the output dir.

    :param blast_dir:  Path, path to the directory containing .blast_tab files
    :param org_lengths:  Dict, the base count for each input sequence
    :param fraglengths:  dictionary of query sequence fragment lengths, only
        needed for BLASTALL output
    :param method:  str, analysis type (ANIb)
    :param logger:  a logger for messages

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
    results = ANIResults(list(org_lengths.keys()), method)

    # Fill diagonal NA values for alignment_length with org_lengths
    for org, length in list(org_lengths.items()):
        results.alignment_lengths[org][org] = length

    # Process .blast_tab files assuming that the filename format holds:
    # org1_vs_org2.blast_tab:
    for blastfile in blastfiles:
        qname, sname = blastfile.stem.split("_vs_")

        # We may have BLAST files from other analyses in the same directory
        # If this occurs, we raise a warning, and skip the file
        if qname not in list(org_lengths.keys()):
            if logger:
                logger.warning(
                    "Query name %s not in input sequence list, skipping %s",
                    qname,
                    blastfile,
                )
            continue
        if sname not in list(org_lengths.keys()):
            if logger:
                logger.warning(
                    "Subject name %s not in input sequence list, skipping %s",
                    sname,
                    blastfile,
                )
            continue
        resultvals = parse_blast_tab(blastfile, fraglengths)
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
def parse_blast_tab(filename: Path, fraglengths: Dict) -> Tuple[int, int, int]:
    """Return (alignment length, similarity errors, mean_pid) tuple.

    :param filename:  Path, path to .blast_tab file
    :param fraglengths:  Optional[Dict], dictionary of fragment lengths for each
        genome.
    :param method:  str, analysis type (ANIb)

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
    # Load output as dataframe
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
    # Add new columns for recalculated alignment length, proportion, and
    # percentage identity
    data["ani_alnlen"] = data["blast_alnlen"] - data["blast_gaps"]
    data["ani_alnids"] = data["ani_alnlen"] - data["blast_mismatch"]
    data["ani_coverage"] = data["ani_alnlen"] / data["qlen"]
    data["ani_pid"] = data["ani_alnids"] / data["qlen"]
    # Filter rows on 'ani_coverage' > 0.7, 'ani_pid' > 0.3
    filtered = data[(data["ani_coverage"] > 0.7) & (data["ani_pid"] > 0.3)]
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
    filtered.to_csv(Path(filename).with_suffix(".blast_tab.dataframe"), sep="\t")
    return aln_length, sim_errors, ani_pid
