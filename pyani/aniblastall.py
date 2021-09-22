# -*- coding: utf-8 -*-
# (c) University of Strathclyde 2021
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
# Copyright (c) 2021 University of Strathclyde
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
"""Code to implement the ANIblastall average nucleotide identity method."""

import logging
import os
import platform
import re
import shutil
import subprocess
from logging import Logger

import pandas as pd

from pathlib import Path
from typing import List, Tuple, Dict, Optional

from Bio import SeqIO

from . import pyani_config
from . import pyani_jobs
from .pyani_tools import BLASTcmds, BLASTfunctions, BLASTexes, ANIResults
from pyani import pyani_files


def get_version(blastall_exe: Path = pyani_config.BLASTALL_DEFAULT) -> str:
    r"""Return BLAST blastall version as a string.

    :param blast_exe:  path to blastall executable

    We expect blastall to return a string as, for example

    .. code-block:: bash

        $ blastall -version
        [blastall 2.2.26] ERROR: Number of database sequences to show \
        one-line descriptions for (V) [ersion] is bad or out of range [? to ?]

    This is concatenated with the OS name.

    The following circumstances are explicitly reported as strings

    - no executable at passed path
    - non-executable file at passed path
    - no version info returned
    - executable cannot be run on this OS
    """
    logger = logging.getLogger(__name__)

    try:
        blastall_path = Path(shutil.which(blastall_exe))  # type:ignore
    except TypeError:
        return f"{blastall_exe} is not found in $PATH"

    if not blastall_path.is_file():  # no executable
        return f"No blastall at {blastall_path}"

    # This should catch cases when the file can't be executed by the user
    if not os.access(blastall_path, os.X_OK):  # file exists but not executable
        return f"blastall exists at {blastall_path} but not executable"

    if platform.system() == "Darwin":
        cmdline = [blastall_exe, "-version"]
    else:
        cmdline = [blastall_exe]

    try:
        result = subprocess.run(
            cmdline,  # type: ignore
            shell=False,
            stdout=subprocess.PIPE,  # type: ignore
            stderr=subprocess.PIPE,
            check=False,  # blastall doesn't return 0
        )

    except OSError:
        logger.warning("blastall executable will not run", exc_info=True)
        return f"blastall exists at {blastall_path} but could not be executed"

    version = re.search(  # type: ignore
        r"(?<=blastall\s)[0-9\.]*", str(result.stderr, "utf-8")
    ).group()

    if 0 == len(version.strip()):
        return f"blastall exists at {blastall_path} but could not retrieve version"

    return f"{platform.system()}_{version} ({blastall_path})"


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
    # method: str,
    outdir: Path,
    format_exe: Optional[Path] = None,
    blast_exe: Optional[Path] = None,
    prefix: str = "ANIBLAST",
) -> BLASTcmds:
    """Return BLASTcmds object for construction of BLAST commands.

    :param outdir:
    :param format_exe:
    :param blast_exe:
    :param prefix:
    """
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


## Generate list of makeblastdb command lines from passed filenames
# def generate_blastdb_commands(
#     filenames: List[Path],
#     outdir: Path,
#     blastdb_exe: Optional[Path] = None,
# ) -> List[Tuple[str, Path]]:
#     """Return list of makeblastdb command-lines for ANIblastall.
#
#     :param filenames:  a list of paths to input FASTA files
#     :param outdir:  path to output directory
#     :param blastdb_exe:  path to the makeblastdb executable
#     :param method:  str, ANI analysis type (ANIblastall)
#     """
#     if blastdb_exe is None:
#         cmdlines = [construct_formatdb_cmd(fname, outdir) for fname in filenames]
#     else:
#         cmdlines = [
#             construct_formatdb_cmd(fname, outdir, blastdb_exe) for fname in filenames
#         ]
#     return cmdlines


# Generate single makeblastdb command line
def construct_formatdb_cmd(
    filename: Path, outdir: Path, blastdb_exe: Path = pyani_config.FORMATDB_DEFAULT
) -> Tuple[str, Path]:
    """Return formatdb command and path to output file.

    :param filename:  Path, input filename
    :param outdir:  Path, path to output directory
    :param blastdb_exe:  Path, path to the formatdb executable
    """
    newfilename = (
        Path(outdir) / Path(filename).name
    )  # Path(filename.name.replace("-fragments", ""))
    shutil.copy(filename, newfilename)
    return (f"{blastdb_exe} -p F -i {newfilename} -t {filename.stem}", newfilename)


# Generate list of BLASTN command lines from passed filenames
def generate_blastall_commands(
    # filenames: List[Path],
    query: Path,
    subject: Path,
    outdir: Path,
    blastall_exe: Optional[Path] = None,
) -> List[str]:
    """Return a list of blastn command-lines for ANIblastall.

    :param query:  a paths to the query's fragmented input FASTA file
    :param subject:  a paths to the subject's fragmented input FASTA file
    :param outdir:  path to output directory
    :param blastall_exe:  path to BLASTALL executable
    :param method:  str, analysis type (ANIblastall)

    Assumes that the fragment sequence input filenames have the form
    ACCESSION-fragments.ext, where the corresponding BLAST database filenames
    have the form ACCESSION.ext. This is the convention followed by the
    fragment_FASTA_files() function above.
    """
    subj_db = outdir / "blastalldbs" / Path(str(subject.name))
    if blastall_exe is None:
        cmdline = construct_blastall_cmdline(query, subj_db, outdir)
    else:
        cmdline = construct_blastall_cmdline(query, subj_db, outdir, blastall_exe)
        return cmdline


# Generate single BLASTALL command line
def construct_blastall_cmdline(
    query: Path,
    subj_db: Path,
    outdir: Path,
    blastall_exe: Path = pyani_config.BLASTALL_DEFAULT,
) -> str:
    """Return single blastall command.

    :param query:  Path, FASTA file for query genome
    :param subj_db:  Path, database of fragments for subject genome
    :param outdir:  Path, to the output directory
    :param blastall_exe:  str, path to blastall executable
    """
    prefix = Path(outdir) / f"{query.stem.replace('-fragments', '')}_vs_{subj_db.stem}"
    return (
        f"{blastall_exe} -p blastn -o {prefix}.blast_tab -i {query} -d {subj_db} "
        "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
    )


# Parse BLASTALL output to get total alignment length and mismatches
def parse_blast_tab(filename: Path, fraglengths: Dict) -> Tuple[int, int, int]:
    """Return (alignment length, similarity errors, mean_pid) tuple.

    :param filename:  Path, path to .blast_tab file
    :param fraglengths:  Optional[Dict], dictionary of fragment lengths for each
        genome.
    :param method:  str, analysis type (ANIblastall)

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
    qname = filename.stem.split("_vs_")[0]
    # Load output as dataframe
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


def process_blast(
    blast_dir: Path,
    org_lengths: Dict,
    fraglengths: Dict,
    mode: str = "ANIblastall",
    logger: Optional[Logger] = None,
) -> ANIResults:
    """Return tuple of ANIb results for .blast_tab files in the output dir.
    :param blast_dir:  Path, path to the directory containing .blast_tab files
    :param org_lengths:  Dict, the base count for each input sequence
    :param fraglengths:  dictionary of query sequence fragment lengths, only
        needed for BLASTALL output
    :param mode:  str, analysis type (ANIb or ANIblastall)
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
    results = ANIResults(list(org_lengths.keys()), mode)

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
