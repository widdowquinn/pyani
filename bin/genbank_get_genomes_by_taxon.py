#!/usr/bin/env python3
#
# genbank_get_genomes_by_taxon.py
#
# Copyright 2015-2016, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""A script to download sequence/class/label data from NCBI

This script takes an NCBI taxonomy identifier (or string, though this is
not always reliable for taxonomy tree subgraphs...) and downloads all genomes
it can find from NCBI in the corresponding taxon subgraph that has
the passed argument as root.
"""

import logging
import logging.handlers
import os
import re
import shutil
import subprocess
import sys
import time
import traceback

from argparse import ArgumentParser
from collections import defaultdict
from socket import timeout
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from Bio import Entrez, SeqIO


class NCBIDownloadException(Exception):
    """General exception for failed NCBI download."""
    def __init__(self):
        Exception.__init__(self, "Error downloading file from NCBI")


# Parse command-line
def parse_cmdline():
    """Parse command-line arguments"""
    parser = ArgumentParser(prog="genbank_get_genomes_by_taxon.py")
    parser.add_argument("-o", "--outdir", dest="outdirname",
                        action="store", default=None,
                        help="Output directory")
    parser.add_argument("-t", "--taxon", dest="taxon",
                        action="store", default=None,
                        help="NCBI taxonomy ID")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Give verbose output")
    parser.add_argument("-f", "--force", dest="force",
                        action="store_true", default=False,
                        help="Force file overwriting")
    parser.add_argument("--noclobber", dest="noclobber",
                        action="store_true", default=False,
                        help="Don't nuke existing files")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None,
                        help="Logfile location")
    parser.add_argument("--format", dest="format",
                        action="store", default="gbk,fasta",
                        help="Output file format [gbk|fasta]")
    parser.add_argument("--email", dest="email",
                        action="store", default=None,
                        help="Email associated with NCBI queries")
    parser.add_argument("--retries", dest="retries",
                        action="store", default=20,
                        help="Number of Entrez retry attempts per request.")
    parser.add_argument("--batchsize", dest="batchsize",
                        action="store", default=10000,
                        help="Entrez record return batch size")
    parser.add_argument("--timeout", dest="timeout",
                        action="store", default=10,
                        help="Timeout for URL connection (s)")
    return parser.parse_args()


# Report last exception as string
def last_exception():
    """ Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))

# Set contact email for NCBI
def set_ncbi_email():
    """Set contact email for NCBI."""
    Entrez.email = args.email
    logger.info("Set NCBI contact email to %s", args.email)
    Entrez.tool = "genbank_get_genomes_by_taxon.py"


# Create output directory if it doesn't exist
def make_outdir():
    """Make the output directory, if required.
    This is a little involved.  If the output directory already exists,
    we take the safe option by default, and stop with an error.  We can,
    however, choose to force the program to go on, in which case we can
    either clobber the existing directory, or not.  The options turn out
    as the following, if the directory exists:
    DEFAULT: stop and report the collision
    FORCE: continue, and remove the existing output directory
    NOCLOBBER+FORCE: continue, but do not remove the existing output
    """
    if os.path.exists(args.outdirname):
        if not args.force:
            logger.error("Output directory %s would overwrite existing " +
                         "files (exiting)", args.outdirname)
            sys.exit(1)
        else:
            logger.info("--force output directory use")
            if args.noclobber:
                logger.warning("--noclobber: existing output directory kept")
            else:
                logger.info("Removing directory %s and everything below it",
                            args.outdirname)
                shutil.rmtree(args.outdirname)
    logger.info("Creating directory %s", args.outdirname)
    try:
        os.makedirs(args.outdirname)   # We make the directory recursively
    except OSError:
        # This gets thrown if the directory exists. If we've forced overwrite/
        # delete and we're not clobbering, we let things slide
        if args.noclobber and args.force:
            logger.info("NOCLOBBER+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            sys.exit(1)

# Retry Entrez requests (or any other function)
def entrez_retry(func, *fnargs, **fnkwargs):
    """Retries the passed function up to the number of times specified
    by args.retries
    """
    tries, success = 0, False
    while not success and tries < args.retries:
        try:
            output = func(*fnargs, **fnkwargs)
            success = True
        except (HTTPError, URLError):
            tries += 1
            logger.warning("Entrez query %s(%s, %s) failed (%d/%d)",
                           func, fnargs, fnkwargs, tries+1, args.retries)
            logger.warning(last_exception())
    if not success:
        logger.error("Too many Entrez failures (exiting)")
        sys.exit(1)
    return output


# Get results from NCBI web history, in batches
def entrez_batch_webhistory(record, expected, batchsize, *fnargs, **fnkwargs):
    """Recovers the Entrez data from a prior NCBI webhistory search, in
    batches of defined size, using Efetch. Returns all results as a list.

    - record: Entrez webhistory record
    - expected: number of expected search returns
    - batchsize: how many search returns to retrieve in a batch
    - *fnargs: arguments to Efetch
    - **fnkwargs: keyword arguments to Efetch
    """
    results = []
    for start in range(0, expected, batchsize):
        batch_handle = entrez_retry(Entrez.efetch,
                                    retstart=start, retmax=batchsize,
                                    webenv=record["WebEnv"],
                                    query_key=record["QueryKey"],
                                    *fnargs, **fnkwargs)
        batch_record = Entrez.read(batch_handle, validate=False)
        results.extend(batch_record)
    return results


# Get assembly UIDs for the root taxon
def get_asm_uids(taxon_uid):
    """Returns a set of NCBI UIDs associated with the passed taxon.

    This query at NCBI returns all assemblies for the taxon subtree
    rooted at the passed taxon_uid.
    """
    query = "txid%s[Organism:exp]" % taxon_uid
    logger.info("Entrez ESearch with query: %s", query)

    # Perform initial search for assembly UIDs with taxon ID as query.
    # Use NCBI history for the search.
    handle = entrez_retry(Entrez.esearch, db="assembly", term=query,
                          format="xml", usehistory="y")
    record = Entrez.read(handle, validate=False)
    result_count = int(record['Count'])
    logger.info("Entrez ESearch returns %d assembly IDs", result_count)

    # Recover assembly UIDs from the web history
    asm_ids = entrez_batch_webhistory(record, result_count, 250,
                                      db="assembly", retmode="xml")
    logger.info("Identified %d unique assemblies", len(asm_ids))
    return asm_ids


# Extract filestem from Entrez eSummary
def extract_filestem(data):
    """Extract filestem from Entrez eSummary data.

    Function expects esummary['DocumentSummarySet']['DocumentSummary'][0]

    Some illegal characters may occur in AssemblyName - for these, a more
    robust regex replace/escape may be required. Sadly, NCBI don't just
    use standard percent escapes, but instead replace certain
    characters with underscores.
    """
    escapes = re.compile(r"[\s/,]")
    escname = re.sub(escapes, '_', data['AssemblyName'])
    return '_'.join([data['AssemblyAccession'], escname])


# Download NCBI assembly file for a passed Assembly UID
def get_ncbi_asm(asm_uid):
    """Returns the NCBI AssemblyAccession and AssemblyName for the assembly
    with passed UID, and organism data for class/label files also, as well
    as accession, so we can track whether downloads fail because only the
    most recent version is available..

    AssemblyAccession and AssemblyName are data fields in the eSummary record,
    and correspond to downloadable files for each assembly at
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GC[AF]/nnn/nnn/nnn/<AA>_<AN>
    where <AA> is AssemblyAccession, and <AN> is AssemblyName, and the choice
    of GCA vs GCF, and the three values of nnn are taken from <AA>
    """
    logger.info("Identifying assembly information from NCBI for %s",
                asm_uid)

    # Obtain full eSummary data for the assembly
    summary = Entrez.read(entrez_retry(Entrez.esummary, db="assembly",
                                       id=asm_uid, report="full"),
                          validate=False)

    # Extract filestem from assembly data
    data = summary['DocumentSummarySet']['DocumentSummary'][0]
    filestem = extract_filestem(data)

    # Report interesting things from the summary for those interested
    logger.info("\tOrganism: %s", data['Organism'])
    logger.info("\tTaxid: %s", data['SpeciesTaxid'])
    logger.info("\tAccession: %s", data['AssemblyAccession'])
    logger.info("\tName: %s", data['AssemblyName'])
    # NOTE: Maybe parse out the assembly stats here, in future?

    # Get class and label text
    organism = data['SpeciesName']
    try:
        strain = data['Biosource']['InfraspeciesList'][0]['Sub_value']
    except (KeyError, IndexError):
        # we consider this an error/incompleteness in the NCBI metadata
        strain = ""

    # Create label and class strings
    genus, species = organism.split(' ', 1)
    labeltxt = "%s_genomic\t%s %s %s" % (filestem, genus[0] + '.',
                                         species, strain)
    classtxt = "%s_genomic\t%s" % (filestem, organism)
    logger.info("\tLabel: %s", labeltxt)
    logger.info("\tClass: %s", classtxt)

    # Download and extract genome assembly
    try:
        fastafilename = retrieve_asm_contigs(filestem)
    except NCBIDownloadException:
        # This is a little hacky. Sometimes, RefSeq assemblies are
        # suppressed (presumably because they are non-redundant),
        # but the GenBank assembly persists. In those cases, we
        # *assume* (because it may not be true) that the corresponding
        # genbank sequence shares the same accession number, except
        # that GCF is replaced by GCA
        gbfilestem = re.sub('^GCF_', 'GCA_', filestem)
        logger.warning("Could not download %s, trying %s",
                       filestem, gbfilestem)
        try:
            fastafilename = retrieve_asm_contigs(gbfilestem)
        except NCBIDownloadException:
            fastafilename = None

    return (fastafilename, classtxt, labeltxt, data['AssemblyAccession'])


# Download and extract an NCBI assembly file, given a filestem
def retrieve_asm_contigs(filestem,
                         ftpstem="ftp://ftp.ncbi.nlm.nih.gov/genomes/all",
                         suffix="genomic.fna.gz"):
    """Downloads an assembly sequence to a local directory.

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
    logger.info("Retrieving assembly sequence for %s", filestem)

    # Compile URL
    gc, aa, an = tuple(filestem.split('_', 2))
    aaval = aa.split('.')[0]
    subdirs = '/'.join([aa[i:i+3] for i in range(0, len(aaval), 3)])
               
    url = "{0}/{1}/{2}/{3}/{3}_{4}".format(ftpstem, gc, subdirs,
                                           filestem, suffix)
    logger.info("Using URL: %s", url)

    # Get data info
    try:
        response = urlopen(url, timeout=args.timeout)
    except (HTTPError, URLError):
        logger.error("Download failed for URL: %s\n%s",
                     url, last_exception())
        raise NCBIDownloadException()
    except timeout:
        logger.error("Download timed out for URL: %s\n%s",
                     url, last_exception())
        raise NCBIDownloadException()
    else:
        fsize = int(response.info().get("Content-length"))
        logger.info("Opened URL and parsed metadata.")

    # Download data
    outfname = os.path.join(args.outdirname, '_'.join([filestem, suffix]))
    if os.path.exists(outfname):
        logger.warning("Output file %s exists, not downloading", outfname)
    else:
        logger.info("Downloading %s (%d bytes)", url, fsize)
        bsize = 1048576  # buffer size
        fsize_dl = 0     # bytes downloaded
        try:
            with open(outfname, "wb") as ofh:
                while True:
                    buffer = response.read(bsize)
                    if not buffer:
                        break
                    fsize_dl += len(buffer)
                    ofh.write(buffer)
                    status = r"%10d  [%3.2f%%]" % (fsize_dl,
                                                   fsize_dl * 100. / fsize)
                    logger.info(status)
        except:
            logger.error("Download failed for %s", url)
            logger.error(last_exception())
            raise NCBIDownloadException()

    # Extract data
    ename = os.path.splitext(outfname)[0]  # Strips only .gz from filename
    # The code below would munge the extracted filename to suit the expected
    # class/label from the old version of this script.
    # The .gz file downloaded from NCBI has format
    # <assembly UID>_<string>_genomic.fna.gz - which we would extract to
    # <assembly UID>.fna
    #regex = ".{3}_[0-9]{9}.[0-9]"
    #outparts = os.path.split(outfname)
    #print(outparts[0])
    #print(re.match(regex, outparts[-1]).group())
    #ename = os.path.join(outparts[0],
    #                     re.match(regex, outparts[-1]).group() + '.fna')
    if os.path.exists(ename):
        logger.warning("Output file %s exists, not extracting", ename)
    else:
        try:
            logger.info("Extracting archive %s to %s",
                        outfname, ename)
            with open(ename, 'w') as efh:
                subprocess.call(['gunzip', '-c', outfname],
                                stdout=efh)  # can be subprocess.run in Py3.5
                logger.info("Archive extracted to %s", ename)
        except:
            logger.error("Extracting archive %s failed", outfname)
            logger.error(last_exception())
            raise NCBIDownloadException()

    return ename


# Write contigs for a single assembly out to file
def write_contigs(asm_uid, contig_uids, batchsize=10000):
    """Writes assembly contigs out to a single FASTA file in the script's
    designated output directory.

    FASTA records are returned, as GenBank and even GenBankWithParts format
    records don't reliably give correct sequence in all cases.

    The script returns two strings for each assembly, a 'class' and a 'label'
    string - this is for use with, e.g. pyani.
    """
    # Has duplicate code with get_class_label_info() - needs refactoring
    logger.info("Collecting contig data for %s", asm_uid)
    # Assembly record - get binomial and strain names
    asm_record = Entrez.read(entrez_retry(Entrez.esummary, db='assembly',
                                          id=asm_uid, rettype='text'),
                             validate=False)
    asm_organism = asm_record['DocumentSummarySet']['DocumentSummary']\
                   [0]['SpeciesName']
    try:
        asm_strain = asm_record['DocumentSummarySet']['DocumentSummary']\
                     [0]['Biosource']['InfraspeciesList'][0]['Sub_value']
    except KeyError:
        asm_strain = ""
    # Assembly UID (long form) for the output filename
    outfilename = "%s.fasta" % os.path.join(args.outdirname,
                                            asm_record['DocumentSummarySet']\
                                            ['DocumentSummary']\
                                            [0]['AssemblyAccession'])

    # Create label and class strings
    genus, species = asm_organism.split(' ', 1)

    # Get FASTA records for contigs
    logger.info("Downloading FASTA records for assembly %s (%s)",
                asm_uid, ' '.join([genus[0] + '.', species, asm_strain]))
    # We're doing an explicit outer retry loop here because we want to confirm
    # we have the correct data, as well as test for Entrez connection errors,
    # which is all the entrez_retry function does.
    tries, success = 0, False
    while not success and tries < args.retries:
        records = []  # Holds all return records
        # We may need to batch contigs
        query_uids = ','.join(contig_uids)
        try:
            for start in range(0, len(contig_uids), batchsize):
                logger.info("Batch: %d-%d", start, start+batchsize)
                records.extend(list(SeqIO.parse(entrez_retry(Entrez.efetch,
                                                             db='nucleotide',
                                                             id=query_uids,
                                                             rettype='fasta',
                                                             retmode='text',
                                                             retstart=start,
                                                             retmax=batchsize),
                                                'fasta')))
            tries += 1
            # Check only that correct number of records returned.
            if len(records) == len(contig_uids):
                success = True
            else:
                logger.warning("%d contigs expected, %d contigs returned",
                               len(contig_uids), len(records))
                logger.warning("FASTA download for assembly %s failed",
                               asm_uid)
                logger.warning("try %d/20", tries)
            # Could also check expected assembly sequence length?
            logger.info("Downloaded genome size: %d",
                        sum([len(r) for r in records]))
        except:
            logger.warning("FASTA download for assembly %s failed", asm_uid)
            logger.warning(last_exception())
            logger.warning("try %d/20", tries)
    if not success:
        # Could place option on command-line to stop or continue here.
        logger.error("Failed to download records for %s (continuing)",
                     asm_uid)

    # Write contigs to file
    retval = SeqIO.write(records, outfilename, 'fasta')
    logger.info("Wrote %d contigs to %s", retval, outfilename)


# Function to report whether an accession has been downloaded
def logreport_downloaded(accession, skippedlist, accessiondict, uidaccdict):
    """Reports to logger whether alternative assemblies for an accession that
    was missing have been downloaded
    """
    for vid in accessiondict[accession.split('.')[0]]:
        if vid in skippedlist:
            status = "NOT DOWNLOADED"
        else:
            status = "DOWNLOADED"
        logger.warning("\t\t%s: %s - %s",
                       vid, uidaccdict[vid], status)

# Run as script
if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline()

    # Set up logging
    logger = logging.getLogger('genbank_get_genomes_by_taxon.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)

    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging",
                         args.logfile)
            sys.exit(1)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info("genbank_get_genomes_by_taxon.py: %s", time.asctime())
    logger.info("command-line: %s", ' '.join(sys.argv))
    logger.info(args)

    # Have we got an output directory? If not, exit.
    if args.email is None:
        logger.error("No email contact address provided (exiting)")
        sys.exit(1)
    set_ncbi_email()

    # Have we got an output directory? If not, exit.
    if args.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    make_outdir()
    logger.info("Output directory: %s", args.outdirname)

    # We might have more than one taxon in a comma-separated list
    taxon_ids = args.taxon.split(',')
    logger.info("Passed taxon IDs: %s", ', '.join(taxon_ids))

    # Get all NCBI assemblies for each taxon UID
    asm_dict = defaultdict(set)
    for tid in taxon_ids:
        asm_dict[tid] = get_asm_uids(tid)
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
            fastafilename, classtxt, labeltxt, accession = get_ncbi_asm(uid)
            # fastafilename is None if there was an error thrown
            if fastafilename is not None:
                contig_dict[uid] = fastafilename
            else:
                logger.error("Skipping download for %s", uid)
                skippedlist.append(uid)
            # Populate dictionaries for all attempted downloads
            classes.append(classtxt)
            labels.append(labeltxt)
            accessiondict[accession.split('.')[0]].append(uid)
            uidaccdict[uid] = accession

    # Write class and label files
    classfilename = os.path.join(args.outdirname, 'classes.txt')
    labelfilename = os.path.join(args.outdirname, 'labels.txt')
    logger.info("Writing classes file to %s", classfilename)
    with open(classfilename, 'w') as ofh:
        ofh.write('\n'.join(classes) + '\n')
    logger.info("Writing labels file to %s", labelfilename)
    with open(labelfilename, 'w') as ofh:
        ofh.write('\n'.join(labels) + '\n')

    # How many downloads did we do/have to skip?
    logger.info("Obtained %d assemblies", len(contig_dict))
    if len(skippedlist):
        logger.warning("Skipped %d downloads through error", len(skippedlist))
        for uid in sorted(skippedlist):
            logger.warning("Assembly UID %s skipped", uid)
            acc = uidaccdict[uid]
            logger.warning("\tUID: %s - accession: %s", uid, acc)
            # Has another version of this genome been successfully dl'ed
            logger.warning("\tAccession %s has versions:", acc.split('.')[0])
            logreport_downloaded(acc, skippedlist, accessiondict, uidaccdict)
            url = "http://www.ncbi.nlm.nih.gov/assembly/%s" % uid
            # Is this a GenBank sequence with no RefSeq counterpart?
            # e.g. http://www.ncbi.nlm.nih.gov/assembly/196191/
            if acc.startswith('GCA'):
                logger.warning("\tAccession is GenBank: does RefSeq exist?")
                logger.warning("\tCheck under 'history' at %s", url)
                # Check for RefSeq counterparts
                rsacc = re.sub('^GCA_', 'GCF_', uidaccdict[uid])
                logger.warning("\tAlternative RefSeq candidate accession: %s",
                               rsacc.split('.')[0])
                logger.warning("\tWere alternative assemblies downloaded?")
                logreport_downloaded(rsacc, skippedlist,
                                     accessiondict, uidaccdict)
            # Is this a suppressed RefSeq sequence?
            if acc.startswith('GCF'):
                logger.warning("\tAccession is RefSeq: is it suppressed?")
                logger.warning("\tCheck under 'history' at %s", url)
                # Check for GenBank counterparts
                gbacc = re.sub('^GCF_', 'GCA_', uidaccdict[uid])
                logger.warning("\tAlternative GenBank candidate accession: %s",
                               gbacc.split('.')[0])
                logger.warning("\tWere alternative assemblies downloaded?")
                logreport_downloaded(gbacc, skippedlist,
                                     accessiondict, uidaccdict)
    logger.info("Skipped assembly UIDs: %s", skippedlist)

    # Let the user know we're done
    logger.info(time.asctime())
    logger.info("Done.")
