#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# 161 Cathedral Street,
# Glasgow,
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2019 The James Hutton Institute
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
"""Script to download from NCBI all genomes in a specified taxon subtree.

This script takes an NCBI taxonomy identifier (or string, though this is
not always reliable for taxonomy tree subgraphs...) and downloads all genomes
it can find from NCBI in the corresponding taxon subgraph that has
the passed argument as root.
"""

import logging
import re
import shutil
import subprocess  # nosec
import sys
import time
import traceback

from argparse import ArgumentParser, Namespace
from collections import defaultdict
from pathlib import Path
from socket import timeout
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from Bio import Entrez, SeqIO

from pyani import __version__
from pyani.download import create_hash
from pyani.scripts.logger import config_logger


class NCBIDownloadException(Exception):

    """General exception for failed NCBI download."""

    def __init__(self):
        """Instantiate exception."""
        Exception.__init__(self, "Error downloading file from NCBI")


# Parse command-line
def parse_cmdline(argv=None):
    """Parse command-line arguments.

    :param argv:  list of command-line arguments
    """
    parser = ArgumentParser(prog="genbank_get_genomes_by_taxon.py")
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdirname",
        required=True,
        action="store",
        default=None,
        type=Path,
        help="Output directory (required)",
    )
    parser.add_argument(
        "-t",
        "--taxon",
        dest="taxon",
        action="store",
        default=None,
        help="NCBI taxonomy ID",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Give verbose output",
    )
    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        default=False,
        help="Report debugging output",
    )
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force file overwriting",
    )
    parser.add_argument(
        "--noclobber",
        dest="noclobber",
        action="store_true",
        default=False,
        help="Don't nuke existing files",
    )
    parser.add_argument(
        "-l",
        "--logfile",
        dest="logfile",
        action="store",
        default=None,
        type=Path,
        help="Logfile location",
    )
    parser.add_argument(
        "--format",
        dest="format",
        action="store",
        default="fasta",
        help="Output file format [gbk|fasta]",
    )
    parser.add_argument(
        "--email",
        dest="email",
        required=True,
        action="store",
        default=None,
        help="Email associated with NCBI queries (required)",
    )
    parser.add_argument(
        "--retries",
        dest="retries",
        action="store",
        default=20,
        help="Number of Entrez retry attempts per request.",
    )
    parser.add_argument(
        "--batchsize",
        dest="batchsize",
        action="store",
        default=10000,
        help="Entrez record return batch size",
    )
    parser.add_argument(
        "--timeout",
        dest="timeout",
        action="store",
        default=10,
        help="Timeout for URL connection (s)",
    )
    # Parse arguments
    if argv is None:
        argv = sys.argv[1:]
    return parser.parse_args(argv)


# Report last exception as string
def last_exception():
    """Return last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))


# Set contact email for NCBI
def set_ncbi_email(args: Namespace) -> None:
    """Set contact email for NCBI.

    :param args:  Namespace, command-line arguments
    """
    logger = logging.getLogger(__name__)

    Entrez.email = args.email
    logger.info("Set NCBI contact email to %s", args.email)
    Entrez.tool = "genbank_get_genomes_by_taxon.py"


# Create output directory if it doesn't exist
def make_outdir(args: Namespace) -> None:
    """Make the output directory, if required.

    :param args:  Namespace, command-line arguments

    This is a little involved.  If the output directory already exists,
    we take the safe option by default, and stop with an error.  We can,
    however, choose to force the program to go on, in which case we can
    either clobber the existing directory, or not.  The options turn out
    as the following, if the directory exists:
    DEFAULT: stop and report the collision
    FORCE: continue, and remove the existing output directory
    NOCLOBBER+FORCE: continue, but do not remove the existing output
    """
    logger = logging.getLogger(__name__)

    if args.outdirname.exists():
        if not args.force:
            logger.error(
                "Output directory %s would overwrite existing files (exiting)",
                args.outdirname,
            )
            raise SystemExit(1)
        logger.info("--force output directory use")
        if args.noclobber:
            logger.warning("--noclobber: existing output directory kept")
        else:
            logger.info(
                "Removing directory %s and everything below it", args.outdirname
            )
            shutil.rmtree(args.outdirname)
    logger.info("Creating directory %s", args.outdirname)
    try:
        args.outdirname.mkdir(exist_ok=True)  # We make the directory recursively
    except OSError:
        # This gets thrown if the directory exists. If we've forced overwrite/
        # delete and we're not clobbering, we let things slide
        if args.noclobber and args.force:
            logger.info("NOCLOBBER+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            raise SystemExit(1)


# Retry Entrez requests (or any other function)
def entrez_retry(args, func, *fnargs, **fnkwargs):
    """Retry the passed function a defined number of times.

    :param args:  Namespace, command-line arguments
    :param func:  func, Entrez function to attempt
    :param *fnargs:  tuple, arguments to the Entrez function
    :param **fnkwargs:  dict, keyword arguments to the Entrez function
    """
    logger = logging.getLogger(__name__)

    tries, success = 0, False
    while not success and tries < args.retries:
        try:
            output = func(*fnargs, **fnkwargs)
            success = True
        except (HTTPError, URLError):
            tries += 1
            logger.warning(
                "Entrez query %s(%s, %s) failed (%d/%d)",
                func,
                fnargs,
                fnkwargs,
                tries + 1,
                args.retries,
            )
            logger.warning(last_exception())
    if not success:
        logger.error("Too many Entrez failures (exiting)")
        raise SystemExit(1)
    return output


# Get results from NCBI web history, in batches
def entrez_batch_webhistory(args, record, expected, batchsize, *fnargs, **fnkwargs):
    """Recover Entrez data from a prior NCBI webhistory search.

    :param args:  Namespace, command-line arguments
    :param record:  Entrez webhistory record
    :param expected:  int, number of expected search returns
    :param batchsize:  int, number of search returns to retrieve in each batch
    :param *fnargs:  tuple, arguments to Efetch
    :param **fnkwargs:  dict, keyword arguments to Efetch

    Recovery is performed in in batches of defined size, using Efetch.
    Returns all results as a list.
    """
    results = []
    for start in range(0, expected, batchsize):
        batch_handle = entrez_retry(
            args,
            Entrez.efetch,
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


# Get assembly UIDs for the root taxon
def get_asm_uids(args, taxon_uid):
    """Return a set of NCBI UIDs associated with the passed taxon.

    :param args:  Namespace, command-line arguments
    :param taxon_uid:  str, NCBI taxon ID

    This query at NCBI returns all assemblies for the taxon subtree
    rooted at the passed taxon_uid.
    """
    logger = logging.getLogger(__name__)

    query = f"txid{taxon_uid}[Organism:exp]"
    logger.info("Entrez ESearch with query: %s", query)

    # Perform initial search for assembly UIDs with taxon ID as query.
    # Use NCBI history for the search.
    handle = entrez_retry(
        args, Entrez.esearch, db="assembly", term=query, format="xml", usehistory="y",
    )
    record = Entrez.read(handle, validate=False)
    result_count = int(record["Count"])
    logger.info("Entrez ESearch returns %d assembly IDs", result_count)

    # Recover assembly UIDs from the web history
    asm_ids = entrez_batch_webhistory(
        args, record, result_count, 250, db="assembly", retmode="xml"
    )
    logger.info("Identified %d unique assemblies", len(asm_ids))
    return asm_ids


# Extract filestem from Entrez eSummary
def extract_filestem(data):
    """Extract filestem from Entrez eSummary data.

    :param data:  Entrez eSummary

    Function expects esummary['DocumentSummarySet']['DocumentSummary'][0]

    Some illegal characters may occur in AssemblyName - for these, a more
    robust regex replace/escape may be required. Sadly, NCBI don't just
    use standard percent escapes, but instead replace certain
    characters with underscores: white space, slash, comma, hash, brackets.
    """
    escapes = re.compile(r"[\s/,#\(\)]")
    escname = re.sub(escapes, "_", data["AssemblyName"])
    return "_".join([data["AssemblyAccession"], escname])


# Download NCBI assembly file for a passed Assembly UID
def get_ncbi_asm(args, asm_uid, fmt="fasta"):
    """Return the NCBI AssemblyAccession and AssemblyName for an assembly.

    :param args:  Namespace, command-line arguments
    :param asm_uid:  NCBI assembly UID
    :param fmt:  str, format to retrieve assembly information

    Returns organism data for class/label files also, as well
    as accession, so we can track whether downloads fail because only the
    most recent version is available..

    AssemblyAccession and AssemblyName are data fields in the eSummary record,
    and correspond to downloadable files for each assembly at
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GC[AF]/nnn/nnn/nnn/<AA>_<AN>
    where <AA> is AssemblyAccession, and <AN> is AssemblyName, and the choice
    of GCA vs GCF, and the three values of nnn are taken from <AA>
    """
    logger = logging.getLogger(__name__)

    logger.info("Identifying assembly information from NCBI for %s", asm_uid)

    # Obtain full eSummary data for the assembly
    summary = Entrez.read(
        entrez_retry(args, Entrez.esummary, db="assembly", id=asm_uid, report="full"),
        validate=False,
    )

    # Extract filestem from assembly data
    data = summary["DocumentSummarySet"]["DocumentSummary"][0]
    filestem = extract_filestem(data)

    # Report interesting things from the summary for those interested
    logger.info("\tOrganism: %s", data["Organism"])
    logger.info("\tTaxid: %s", data["SpeciesTaxid"])
    logger.info("\tAccession: %s", data["AssemblyAccession"])
    logger.info("\tName: %s", data["AssemblyName"])
    # NOTE: Maybe parse out the assembly stats here, in future?

    # Get class and label text
    organism = data["SpeciesName"]
    try:
        strain = data["Biosource"]["InfraspeciesList"][0]["Sub_value"]
    except (KeyError, IndexError):
        # we consider this an error/incompleteness in the NCBI metadata
        strain = ""

    # Download and extract genome assembly
    hash_md5 = None
    try:
        fastafname = retrieve_asm_contigs(args, filestem, fmt=fmt)
        hash_md5 = create_hash(fastafname)
    except NCBIDownloadException:
        # This is a little hacky. Sometimes, RefSeq assemblies are
        # suppressed (presumably because they are non-redundant),
        # but the GenBank assembly persists. In those cases, we
        # *assume* (because it may not be true) that the corresponding
        # genbank sequence shares the same accession number, except
        # that GCF is replaced by GCA
        gbfilestem = re.sub("^GCF_", "GCA_", filestem)
        logger.warning("Could not download %s, trying %s", filestem, gbfilestem)
        try:
            fastafname = retrieve_asm_contigs(args, gbfilestem, fmt=fmt)
            hash_md5 = create_hash(fastafname)
        except NCBIDownloadException:
            fastafname = None

    # Create label and class strings
    genus, species = organism.split(" ", 1)
    lbltxt = "%s\t%s_genomic\t%s %s %s" % (
        hash_md5,
        filestem,
        genus[0] + ".",
        species,
        strain,
    )
    clstxt = "%s\t%s_genomic\t%s" % (hash_md5, filestem, organism)
    logger.info("\tLabel: %s", lbltxt)
    logger.info("\tClass: %s", clstxt)

    return (fastafname, clstxt, lbltxt, data["AssemblyAccession"])


# Download and extract an NCBI assembly file, given a filestem
def retrieve_asm_contigs(
    args, filestem, ftpstem="ftp://ftp.ncbi.nlm.nih.gov/genomes/all", fmt="fasta",
):
    """Download assembly sequence to a local directory.

    :param args:  Namespace, command-line arguments
    :param filestem:  str, filestem for output file
    :param ftpstem:  str, URI stem for NCBI FTP site
    :param fmt:  str, format for output file

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

    This function downloads the genomic_fna.gz file, and extracts it in the
    output directory name specified when the script is called.
    """
    logger = logging.getLogger(__name__)

    logger.info("Retrieving assembly sequence for %s", filestem)

    # Define format suffix
    logger.info("%s format requested", fmt)
    if fmt == "fasta":
        suffix = "genomic.fna.gz"
    elif fmt == "gbk":
        suffix = "genomic.gbff.gz"

    # Compile URL
    fnameparts = tuple(filestem.split("_", 2))  # three elements: GC*, AA, discard
    subdirs = "/".join(
        [
            fnameparts[1][i : i + 3]
            for i in range(0, len(fnameparts[1].split(".")[0]), 3)
        ]
    )

    asmurl = "{0}/{1}/{2}/{3}/{3}_{4}".format(
        ftpstem, fnameparts[0], subdirs, filestem, suffix
    )
    logger.info("Using URL: %s", asmurl)

    # Get data info
    try:
        response = urlopen(asmurl, timeout=args.timeout)
    except HTTPError:
        logger.error("Download failed for URL: %s", asmurl, exc_info=True)
        raise NCBIDownloadException()
    except URLError as err:
        if isinstance(err.reason, timeout):
            logger.error("Download timed out for URL: %s", asmurl, exc_info=True)
        else:
            logger.error("Download failed for URL: %s", asmurl, exc_info=True)
        raise NCBIDownloadException()
    except timeout:
        logger.error("Download timed out for URL: %s", asmurl, exc_info=True)
        raise NCBIDownloadException()
    fsize = int(response.info().get("Content-length"))
    logger.info("Opened URL and parsed metadata.")

    # Download data
    outfname = args.outdirname / f"{filestem}_{suffix}"
    if outfname.exists():
        logger.warning("Output file %s exists, not downloading", outfname)
    else:
        logger.info("Downloading %s (%s bytes)", asmurl, fsize)
        bsize = 1_048_576  # buffer size
        fsize_dl = 0  # bytes downloaded
        try:
            with open(outfname, "wb") as outfh:
                while True:
                    buffer = response.read(bsize)
                    if not buffer:
                        break
                    fsize_dl += len(buffer)
                    outfh.write(buffer)
                    status = r"%10d  [%3.2f%%]" % (fsize_dl, fsize_dl * 100.0 / fsize)
                    logger.info(status)
        except IOError:
            logger.error("Download failed for %s", asmurl)
            logger.error(last_exception())
            raise NCBIDownloadException()

        # Extract gzip archive and return path to extracted file
        return extract_archive(outfname)


def extract_archive(archivepath):
    """Return path to extracted gzip file.

    :param archivepath:  Path, path to gzipped file with ".tar.gz" suffix
    """
    logger = logging.getLogger(__name__)

    # Extract data from targzed file.
    if archivepath.suffix == ".gz":
        ename = archivepath.with_suffix("")  # Strips only .gz from filename
    else:
        logger.info("Expected .gz file, got %s - not extracting", ename)
    if ename.exists():
        logger.warning("Output file %s exists - not extracting", ename)
    else:
        logger.info("Extracting archive %s to %s", archivepath, ename)
        try:
            with open(ename, "w") as efh:
                subprocess.call(
                    ["gunzip", "-c", archivepath], stdout=efh
                )  # can be subprocess.run in Py3.5
                logger.info("Archive extracted to %s", ename)
        except IOError:
            logger.error("Extracting archive %s failed", archivepath, exc_info=True)
            raise NCBIDownloadException()

    return ename


# Write contigs for a single assembly out to file
def write_contigs(args, asm_uid, contig_uids, batchsize=10000):
    """Write assembly contigs to a single FASTA file.

    :param args:  Namespace, command-line arguments
    :param asm_uid:  str, NCBI assembly UID
    :param contig_uids:
    :param batchsize:  int

    FASTA records are returned, as GenBank and even GenBankWithParts format
    records don't reliably give correct sequence in all cases.

    The script returns two strings for each assembly, a 'class' and a 'label'
    string - this is for use with, e.g. pyani.
    """
    logger = logging.getLogger(__name__)

    # Has duplicate code with get_class_label_info() - needs refactoring
    logger.info("Collecting contig data for %s", asm_uid)
    # Assembly record - get binomial and strain names
    asm_record = Entrez.read(
        entrez_retry(args, Entrez.esummary, db="assembly", id=asm_uid, rettype="text"),
        validate=False,
    )
    asm_smry = asm_record["DocumentSummarySet"]["DocumentSummary"][0]
    asm_organism = asm_smry["SpeciesName"]
    try:
        asm_strain = asm_smry["Biosource"]["InfraspeciesList"][0]["Sub_value"]
    except KeyError:
        asm_strain = ""
    # Assembly UID (long form) for the output filename
    outfilename = f"args.outdirname, {asm_smry['AssemblyAccession']}.fasta"

    # Create label and class strings
    genus, species = asm_organism.split(" ", 1)

    # Get FASTA records for contigs
    logger.info(
        "Downloading FASTA records for assembly %s (%s)",
        asm_uid,
        " ".join([genus[0] + ".", species, asm_strain]),
    )
    # We're doing an explicit outer retry loop here because we want to confirm
    # we have the correct data, as well as test for Entrez connection errors,
    # which is all the entrez_retry function does.
    tries, success = 0, False
    while not success and tries < args.retries:
        records = []  # Holds all return records
        # We may need to batch contigs
        query_uids = ",".join(contig_uids)
        try:
            for start in range(0, len(contig_uids), batchsize):
                logger.info("Batch: %d-%d", start, start + batchsize)
                records.extend(
                    list(
                        SeqIO.parse(
                            entrez_retry(
                                args,
                                Entrez.efetch,
                                db="nucleotide",
                                id=query_uids,
                                rettype="fasta",
                                retmode="text",
                                retstart=start,
                                retmax=batchsize,
                            ),
                            "fasta",
                        )
                    )
                )
            tries += 1
            # Check only that correct number of records returned.
            if len(records) == len(contig_uids):
                success = True
            else:
                logger.warning(
                    "%d contigs expected, %d contigs returned",
                    len(contig_uids),
                    len(records),
                )
                logger.warning("FASTA download for assembly %s failed", asm_uid)
                logger.warning("try %d/20", tries)
            # Could also check expected assembly sequence length?
            logger.info("Downloaded genome size: %d", sum([len(r) for r in records]))
        except HTTPError:
            logger.warning("FASTA download for assembly %s failed", asm_uid)
            logger.warning(last_exception())
            logger.warning("try %d/20", tries)
    if not success:
        # Could place option on command-line to stop or continue here.
        logger.error("Failed to download records for %s (continuing)", asm_uid)

    # Write contigs to file
    retval = SeqIO.write(records, outfilename, "fasta")
    logger.info("Wrote %d contigs to %s", retval, outfilename)


# Function to report whether an accession has been downloaded
def logreport_downloaded(accn, skiplist, accndict, uidaccndict):
    """Report to logger if alternative assemblies were downloaded.

    :param accn:
    :param skiplist:
    :param accndict:
    :param uidaccndict:
    """
    logger = logging.getLogger(__name__)

    for vid in accndict[accn.split(".")[0]]:
        if vid in skiplist:
            status = "NOT DOWNLOADED"
        else:
            status = "DOWNLOADED"
        logger.warning("\t\t%s: %s - %s", vid, uidaccndict[vid], status)


# Run as script
def run_main(args=None):
    """Run main process for average_nucleotide_identity.py script.

    :param args:  Namespace, command-line arguments
    """
    logger = logging.getLogger(__name__)

    # If we need to (i.e. a namespace isn't passed), parse the command-line
    if args is None:
        args = parse_cmdline()
        config_logger(args)

    # Catch execution with no arguments
    if len(sys.argv) == 1:
        sys.stderr.write("pyani version: {0}\n".format(__version__))
        return 0

    # Have we got an email address? If not, exit.
    if args.email is None:
        logger.error("No email contact address provided (exiting)")
        raise SystemExit(1)
    set_ncbi_email(args)

    # Have we got an output directory? If not, exit.
    if args.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    make_outdir(args)
    logger.info("Output directory: %s", args.outdirname)

    # We might have more than one taxon in a comma-separated list
    taxon_ids = args.taxon.split(",")
    logger.info("Passed taxon IDs: %s", ", ".join(taxon_ids))

    # Get all NCBI assemblies for each taxon UID
    asm_dict = defaultdict(set)
    for tid in taxon_ids:
        asm_dict[tid] = get_asm_uids(args, tid)
    for tid, asm_uids in list(asm_dict.items()):
        logger.info("Taxon %s: %d assemblies", tid, len(asm_uids))

    # Download contigs for each assembly UID
    classes, labels = [], []
    contig_dict = defaultdict(set)
    accessiondict = defaultdict(list)  # UIDs, keyed by accession
    uidaccdict = {}  # accessions, keyed by UID
    skippedlist = []
    for tid, asm_uids in list(asm_dict.items()):
        for uid in asm_uids:
            (fastafilename, classtxt, labeltxt, accession) = get_ncbi_asm(
                args, uid, args.format
            )
            # fastafilename is None if there was an error thrown
            if fastafilename is not None:
                contig_dict[uid] = fastafilename
            else:
                logger.error("Skipping download for %s", uid)
                skippedlist.append(uid)
            # Populate dictionaries for all attempted downloads
            classes.append(classtxt)
            labels.append(labeltxt)
            accessiondict[accession.split(".")[0]].append(uid)
            uidaccdict[uid] = accession

    # Write class and label files
    classfilename = args.outdirname / "classes.txt"
    labelfilename = args.outdirname / "labels.txt"
    logger.info("Writing classes file to %s", classfilename)
    with open(classfilename, "w") as ofh:
        ofh.write("\n".join(classes) + "\n")
    logger.info("Writing labels file to %s", labelfilename)
    with open(labelfilename, "w") as ofh:
        ofh.write("\n".join(labels) + "\n")

    # How many downloads did we do/have to skip?
    logger.info("Obtained %d assemblies", len(contig_dict))
    if skippedlist:
        logger.warning("Skipped %d downloads through error", len(skippedlist))
        for uid in sorted(skippedlist):
            logger.warning("Assembly UID %s skipped", uid)
            acc = uidaccdict[uid]
            logger.warning("\tUID: %s - accession: %s", uid, acc)
            # Has another version of this genome been successfully dl'ed
            logger.warning("\tAccession %s has versions:", acc.split(".")[0])
            logreport_downloaded(acc, skippedlist, accessiondict, uidaccdict)
            url = "http://www.ncbi.nlm.nih.gov/assembly/%s" % uid
            # Is this a GenBank sequence with no RefSeq counterpart?
            # e.g. http://www.ncbi.nlm.nih.gov/assembly/196191/
            if acc.startswith("GCA"):
                logger.warning("\tAccession is GenBank: does RefSeq exist?")
                logger.warning("\tCheck under 'history' at %s", url)
                # Check for RefSeq counterparts
                rsacc = re.sub("^GCA_", "GCF_", uidaccdict[uid])
                logger.warning(
                    "\tAlternative RefSeq candidate accession: %s", rsacc.split(".")[0]
                )
                logger.warning("\tWere alternative assemblies downloaded?")
                logreport_downloaded(rsacc, skippedlist, accessiondict, uidaccdict)
            # Is this a suppressed RefSeq sequence?
            if acc.startswith("GCF"):
                logger.warning("\tAccession is RefSeq: is it suppressed?")
                logger.warning("\tCheck under 'history' at %s", url)
                # Check for GenBank counterparts
                gbacc = re.sub("^GCF_", "GCA_", uidaccdict[uid])
                logger.warning(
                    "\tAlternative GenBank candidate accession: %s", gbacc.split(".")[0]
                )
                logger.warning("\tWere alternative assemblies downloaded?")
                logreport_downloaded(gbacc, skippedlist, accessiondict, uidaccdict)
    logger.info("Skipped assembly UIDs: %s", skippedlist)

    # Let the user know we're done
    logger.info(time.asctime())
    logger.info("Done.")

    # Exit
    return 0
