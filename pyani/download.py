# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G1 1XQ
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
"""Module providing functions useful for downloading genomes from NCBI."""

import hashlib
import re
import shlex
import subprocess
import sys
import traceback
import urllib.request

from pathlib import Path
from subprocess import CompletedProcess
from typing import Any, List, NamedTuple, Tuple
from urllib.error import HTTPError, URLError

from Bio import Entrez  # type: ignore
from tqdm import tqdm  # type: ignore
from namedlist import namedlist  # type: ignore


# Regular expression for NCBI taxon numbers
TAXONREGEX = re.compile(r"([0-9]\,?){1,}")


# Custom exceptions
class NCBIDownloadException(Exception):

    """General exception for failed NCBI download."""

    def __init__(self, msg: str = "Error downloading file from NCBI"):
        """Instantiate class."""
        Exception.__init__(self, msg)


class FileExistsException(Exception):

    """A specified file exists."""

    def __init__(self, msg: str = "Specified file exists"):
        """Instantiate class."""
        Exception.__init__(self, msg)


class ASMIDs(NamedTuple):

    """Matching Assembly ID information for a query taxID."""

    query: str
    result_count: int
    asm_ids: List[str]


class Classification(NamedTuple):

    """Taxonomic classification for an isolate."""

    organism: str
    genus: str
    species: str
    strain: str


class Hashstatus(NamedTuple):

    """Status report on file hash comparison."""

    passed: bool
    localhash: str
    filehash: str


def last_exception() -> str:
    """Return last exception as a string."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))


def set_ncbi_email(email: str) -> None:
    """Set contact email for NCBI.

    :param email:  str, email address to give to Entrez at NCBI
    """
    Entrez.email = email
    Entrez.tool = "pyani.py"


def parse_api_key(api_path: Path) -> str:
    """Parse the contents of the file for NCBI API key.

    :param api_path:  Path, location of NCBI API key

    This function expects that the indicated file contains a single line
    (perhaps with a newline) containing no data other than an NCBI API key.
    """
    with api_path.open() as ifh:
        api_key = ifh.readline().strip()
    return api_key


# Get results from NCBI web history, in batches
def entrez_batch_webhistory(
    record, expected, batchsize, retries, *fnargs, **fnkwargs
) -> List[str]:
    """Recover the Entrez data from a prior NCBI webhistory search.

    :param record:  Entrez webhistory record
    :param expected:  number of expected search returns
    :param batchsize:  how many search returns to retrieve in a batch
    :param retries:  int
    :param *fnargs:  arguments to Efetch
    :param **fnkwargs:  keyword arguments to Efetch

    Recovers results in batches of defined size, using Efetch.
    Returns all results as a list.
    """
    results = []  # type: List[Any]
    for start in range(0, expected, batchsize):
        batch_handle = entrez_retry(
            Entrez.efetch,
            retries,
            retstart=start,
            retmax=batchsize,
            webenv=record["WebEnv"],
            query_key=record["QueryKey"],
            *fnargs,
            **fnkwargs,
        )
        batch_record = Entrez.read(batch_handle, validate=False)
        results.extend(batch_record)
    return results


# Retry an Entrez query a specified number of times
def entrez_retry(func, retries, *fnargs, **fnkwargs):
    """Retry the passed function up to the number of times specified.

    :param func:  function to be executed
    :param retries:  int, number of times to retry function execution
    :param *fnargs:  optional arguments to passed function
    :param **fnkwargs:  optional keyword arguments to passed function
    """
    tries, success = 0, False
    while not success and tries < retries:
        try:
            output = func(*fnargs, **fnkwargs)
            success = True
        except (HTTPError, URLError):
            tries += 1
    if not success:
        raise NCBIDownloadException("Too many Entrez failures")
    return output


# Split a list of taxon ids into components, checking for correct formatting
def split_taxa(taxa: str) -> List[str]:
    """Return list of taxon ids from the passed comma-separated list.

    :param taxa:  str, comma-separated list of valid NCBI taxonomy IDs

    The function checks the passed taxon argument against a regular expression
    that permits comma-separated numerical symbols only.
    """
    # Check format of passed taxa
    match = TAXONREGEX.match(taxa)
    if match is None or len(match.group()) != len(taxa):
        raise ValueError("invalid taxon string: {0}".format(taxa))
    return [taxon for taxon in taxa.split(",") if len(taxon)]


# Get assembly UIDs for the subtree rooted at the passed taxon
def get_asm_uids(taxon_uid: str, retries: int) -> ASMIDs:
    """Return set of NCBI UIDs associated with the passed taxon UID.

    :param taxon_uid:  str, NCBI taxID for taxon to download
    :param retries:  int, number of download retry attempts

    This query at NCBI returns all assemblies for the taxon subtree
    rooted at the passed taxon_uid.
    """
    query = "txid%s[Organism:exp]" % taxon_uid

    # Perform initial search for assembly UIDs with taxon ID as query.
    # Use NCBI history for the search.
    handle = entrez_retry(
        Entrez.esearch, retries, db="assembly", term=query, format="xml", usehistory="y"
    )
    record = Entrez.read(handle, validate=False)
    result_count = int(record["Count"])

    # Recover assembly UIDs from the web history
    asm_ids = entrez_batch_webhistory(
        record, result_count, 250, retries, db="assembly", retmode="xml"
    )

    return ASMIDs(query, result_count, asm_ids)


# Get a filestem from Entrez eSummary data
def extract_filestem(esummary) -> str:
    """Extract filestem from Entrez eSummary data.

    :param esummary:

    Function expects esummary['DocumentSummarySet']['DocumentSummary'][0]

    Some illegal characters may occur in AssemblyName - for these, a more
    robust regex replace/escape may be required. Sadly, NCBI don't just
    use standard percent escapes, but instead replace certain
    characters with underscores: white space, slash, comma, hash, brackets.
    """
    escapes = re.compile(r"[\s/,#\(\)]")
    escname = re.sub(escapes, "_", esummary["AssemblyName"])
    return "_".join([esummary["AssemblyAccession"], escname])


# Get eSummary data for a single assembly UID
def get_ncbi_esummary(asm_uid, retries, api_key=None) -> Tuple:
    """Obtain full eSummary info for the passed assembly UID.

    :param asm_uid:
    :param retries:
    :param api_key:
    """
    # Obtain full eSummary data for the assembly
    summary = Entrez.read(
        entrez_retry(
            Entrez.esummary,
            retries,
            db="assembly",
            id=asm_uid,
            report="full",
            api_key=api_key,
        ),
        validate=False,
    )

    # Extract filestem from assembly data
    data = summary["DocumentSummarySet"]["DocumentSummary"][0]
    filestem = extract_filestem(data)

    return (data, filestem)


# Get the taxonomic classification strings for eSummary data
def get_ncbi_classification(esummary) -> Classification:
    """Return organism, genus, species, strain info from eSummary data.

    :param esummary:
    """
    # Extract species/strain info
    organism = esummary["SpeciesName"]
    try:
        strain = esummary["Biosource"]["InfraspeciesList"][0]["Sub_value"]
    except (KeyError, IndexError):
        # we consider this an error/incompleteness in the NCBI metadata
        strain = ""
    genus, species = organism.split(" ", 1)

    return Classification(organism, genus, species, strain)


# Given a remote filestem, generate URIs for download
def compile_url(filestem: str, suffix: str, ftpstem: str) -> Tuple[str, str]:
    """Compile download URLs given a passed filestem.

    :param filestem:
    :param suffix:
    :param ftpstem:

    The filestem corresponds to <AA>_<AN>, where <AA> and <AN> are
    AssemblyAccession and AssemblyName: data fields in the eSummary record.
    These correspond to downloadable files for each assembly at
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GC[AF]/nnn/nnn/nnn/<AA>_<AN>/
    where <AA> is AssemblyAccession, and <AN> is AssemblyName. The choice
    of GCA vs GCF, and the values of nnn, are derived from <AA>

    The files in this directory all have the stem <AA>_<AN>_<suffix>, where
    suffixes are:
    assembly_report.txt
    assembly_stats.txt
    feature_table.txt.gz
    genomic.fna.gz
    genomic.gbff.gz
    genomic.gff.gz
    protein.faa.gz
    protein.gpff.gz
    rm_out.gz
    rm.run
    wgsmaster.gbff.gz
    """
    gcstem, acc, _ = tuple(filestem.split("_", 2))
    aaval = acc.split(".")[0]
    subdirs = "/".join([acc[i : i + 3] for i in range(0, len(aaval), 3)])

    url = f"{ftpstem}/{gcstem}/{subdirs}/{filestem}/{filestem}_{suffix}"
    hashurl = f"{ftpstem}/{gcstem}/{subdirs}/{filestem}/md5checksums.txt"
    return (url, hashurl)


# Download a remote file to the specified directory
def download_url(
    url: str, outfname: Path, disable_tqdm: bool = False
) -> None:
    """Download remote URL to a local directory.

    :param url:  URL of remote file for download
    :param outfname: Path, path to write output
    :param disable_tqdm:  Boolean, show tqdm progress bar?

    This function downloads the contents of the passed URL to the passed
    filename, in buffered chunks
    """
    # Open connection, and get expected filesize
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as response:
        fsize = int(response.info().get("Content-length"))

        # Define buffer sizes
        bsize = 1048576  # buffer size
        fsize_dl = 0  # bytes downloaded

        # Download file
        with open(outfname, "wb") as ofh:
            with tqdm(total=fsize, disable=disable_tqdm, desc=outfname.name) as pbar:
                while True:
                    buffer = response.read(bsize)
                    if not buffer:
                        break
                    fsize_dl += len(buffer)
                    ofh.write(buffer)
                    pbar.update(bsize)


# Construct filepaths for downloaded files and their hashes
def construct_output_paths(
    filestem: str, suffix: str, outdir: Path
) -> Tuple[Path, Path]:
    """Construct paths to output files for genome and hash.

    :param filestem:  str, output filename stem
    :param suffix:  str, output filename suffix
    :param outdir:  Path, path to output directory
    """
    outfname = outdir / "_".join([filestem, suffix])
    outfhash = outdir / "_".join([filestem, "hashes.txt"])
    return (outfname, outfhash)


# Download a remote genome from NCBI and its MD5 hash
def retrieve_genome_and_hash(
    filestem: str,
    suffix: str,
    ftpstem: str,
    outdir: Path,
    timeout: int,
    disable_tqdm: bool = False,
) -> namedlist:
    """Download genome contigs and MD5 hash data from NCBI.

    :param filestem:
    :param suffix:
    :param ftpstem:
    :param outdir:
    :param timeout:
    :param disable_tqdm:  Boolean, show tqdm progress bar?
    """
    DLStatus = namedlist("DLStatus", "url hashurl outfname outfhash skipped error")
    skipped = False  # Flag - set True if we skip download for existing file
    error = None  # Text of last-raised error

    # Construct remote URLs and output filenames
    url, hashurl = compile_url(filestem, suffix, ftpstem)
    outfname, outfhash = construct_output_paths(filestem, suffix, outdir)

    # Download the genome sequence and corresponding hash file
    try:
        download_url(url, outfname, timeout, disable_tqdm)
        download_url(hashurl, outfhash, timeout, disable_tqdm)
    except IOError:
        error = last_exception()

    return DLStatus(url, hashurl, outfname, outfhash, skipped, error)


# Check the file hash against the downloaded hash
def check_hash(fname: Path, hashfile: Path) -> Hashstatus:
    """Check MD5 of passed file against downloaded NCBI hash file.

    :param fname:  Path, path to local hash file
    :param hashfile:  Path, path to NCBI hash file
    """
    filehash = ""
    passed = False  # Flag - set to True if the hash matches

    # Generate MD5 hash
    localhash = create_hash(fname)

    # Get hash from file
    filehash = extract_hash(hashfile, fname.name)

    # Check for match
    if filehash == localhash:
        passed = True

    return Hashstatus(passed, localhash, filehash)


# Extract contigs from a compressed file, using gunzip
def extract_contigs(fname: Path, ename: Path) -> CompletedProcess:
    """Extract contents of fname to ename using gunzip.

    :param fname:  str, path to input compressed file
    :param ename:  str, path to output uncompressed file

    Returns status of subprocess call
    """
    cmd = ["gunzip", "-c", shlex.quote(str(fname))]
    with open(ename, "w") as efh:
        return subprocess.run(cmd, stdout=efh, check=True, shell=False)


# Using a genomes UID, create class and label text files
def create_labels(
    classification: Classification, filestem: str, genomehash: str
) -> Tuple[str, str]:
    r"""Return class and label text from UID classification.

    :param classification:  Classification named tuple (org, genus, species, strain)
    :param filestem:  str, filestem of input genome file
    :param genomehash:  str, MD5 hash of genome data

    The 'class' data is the organism as provided in the passed Classification
    named tuple; the 'label' data is genus, species and strain information
    from the same tuple. The label is intended to be human-readable, the class
    data to be a genuine class identifier.

    Returns a tuple of two strings: (label, class).

    The two strings are tab-separated strings: <HASH>\t<FILE>\t<CLASS/LABEL>.
    The hash is used to help uniquely identify the genome in the database
    (label/class is unique by a combination of hash and run ID).
    """
    return (
        (
            f"{genomehash}\t{filestem}_genomic\t{classification.genus[0] + '.'} "
            f"{classification.species} {classification.strain}"
        ),
        f"{genomehash}\t{filestem}_genomic\t{classification.organism}",
    )


# Create an MD5 hash for the passed genome
def create_hash(fname: Path) -> str:
    """Return MD5 hash of the passed file contents.

    :param fname:  Path, path to file for hashing

    We can ignore the Bandit B303 error as we're not using the hash for
    cryptographic purposes.
    """
    hash_md5 = hashlib.md5()  # nosec
    with open(fname, "rb") as fhandle:
        for chunk in iter(lambda: fhandle.read(65536), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


# Create an MD5 hash for the passed genome
def extract_hash(hashfile: Path, name: str) -> str:
    """Return MD5 hash from file of name:MD5 hashes.

    :param hashfile:  Path, path to file containing name:MD5 pairs
    :param name:  str, name associated with hash
    """
    filehash = None
    with open(hashfile, "r") as hhandle:
        for line in [_.strip().split() for _ in hhandle if len(_.strip())]:
            if Path(line[1]).name == name:  # hash filename
                filehash = line[0]
    return str(filehash)
