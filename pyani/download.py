# Copyright 2017, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

import hashlib
import os
import re
import sys
import subprocess
import traceback

from collections import namedtuple
from socket import timeout
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from Bio import Entrez
from tqdm import tqdm
from namedlist import namedlist


taxonregex = re.compile('([0-9]\,?){1,}')


# Custom exceptions
class NCBIDownloadException(Exception):
    """General exception for failed NCBI download."""
    def __init__(self, msg="Error downloading file from NCBI"):
        Exception.__init__(self, msg)

class FileExistsException(Exception):
    """A specified file exists."""
    def __init__(self, msg="Specified file exists"):
        Exception.__init__(self, msg)


def last_exception():
    """ Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


def set_ncbi_email(email):
    """Set contact email for NCBI."""
    Entrez.email = email
    Entrez.tool = "pyani.py"


# Get results from NCBI web history, in batches
def entrez_batch_webhistory(record, expected, batchsize, retries,
                            *fnargs, **fnkwargs):
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
        batch_handle = entrez_retry(Entrez.efetch, retries,
                                    retstart=start, retmax=batchsize,
                                    webenv=record["WebEnv"],
                                    query_key=record["QueryKey"],
                                    *fnargs, **fnkwargs)
        batch_record = Entrez.read(batch_handle, validate=False)
        results.extend(batch_record)
    return results


# Retry an Entrez query a specified number of times
def entrez_retry(func, retries, *fnargs, **fnkwargs):
    """Retries the passed function up to the number of times specified
    by retries
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
def split_taxa(taxa):
    """Returns a list of taxon ids from the passed comma-separated list.

    The function checks the passed taxon argument against a regular expression
    that permits comma-separated numerical symbols only.
    """
    # Check format of passed taxa
    match = taxonregex.match(taxa)
    if match is None or len(match.group()) != len(taxa):
        raise ValueError("invalid taxon string: {0}".format(taxa))
    return [taxon for taxon in taxa.split(',') if len(taxon)]


# Get assembly UIDs for the subtree rooted at the passed taxon
def get_asm_uids(taxon_uid, retries):
    """Returns a set of NCBI UIDs associated with the passed taxon UID.

    This query at NCBI returns all assemblies for the taxon subtree
    rooted at the passed taxon_uid.
    """
    Results = namedtuple("ASM_UIDs", "query count asm_ids")
    query = "txid%s[Organism:exp]" % taxon_uid

    # Perform initial search for assembly UIDs with taxon ID as query.
    # Use NCBI history for the search.
    handle = entrez_retry(Entrez.esearch, retries, db="assembly",
                          term=query, format="xml", usehistory="y")
    record = Entrez.read(handle, validate=False)
    result_count = int(record['Count'])

    # Recover assembly UIDs from the web history
    asm_ids = entrez_batch_webhistory(record, result_count, 250, retries,
                                      db="assembly", retmode="xml")

    return Results(query, result_count, asm_ids)
        

# Get a filestem from Entrez eSummary data
def extract_filestem(esummary):
    """Extract filestem from Entrez eSummary data.

    Function expects esummary['DocumentSummarySet']['DocumentSummary'][0]

    Some illegal characters may occur in AssemblyName - for these, a more
    robust regex replace/escape may be required. Sadly, NCBI don't just
    use standard percent escapes, but instead replace certain
    characters with underscores: white space, slash, comma, hash, brackets.
    """
    escapes = re.compile(r"[\s/,#\(\)]")
    escname = re.sub(escapes, '_', esummary['AssemblyName'])
    return '_'.join([esummary['AssemblyAccession'], escname])


# Get eSummary data for a single assembly UID
def get_ncbi_esummary(asm_uid, retries):
    """Obtain full eSummary info for the passed assembly UID."""
    # Obtain full eSummary data for the assembly
    summary = Entrez.read(entrez_retry(Entrez.esummary, retries,
                                       db="assembly",
                                       id=asm_uid, report="full"),
                          validate=False)

    # Extract filestem from assembly data
    data = summary['DocumentSummarySet']['DocumentSummary'][0]
    filestem = extract_filestem(data)

    return(data, filestem)


# Get the taxonomic classification strings for eSummary data
def get_ncbi_classification(esummary):
    """Return organism, genus, species, strain info from eSummary data."""
    Classification = namedtuple("Classsification",
                                "organism genus species strain")

    # Extract species/strain info
    organism = esummary['SpeciesName']
    try:
        strain = esummary['Biosource']['InfraspeciesList'][0]\
                 ['Sub_value']
    except (KeyError, IndexError):
        # we consider this an error/incompleteness in the NCBI metadata
        strain = ""
    genus, species = organism.split(' ', 1)

    return Classification(organism, genus, species, strain)


# Given a remote filestem, generate URIs for download
def compile_url(filestem, suffix, ftpstem):
    """Compile download URLs given a passed filestem.

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
    gc, aa, an = tuple(filestem.split('_', 2))
    aaval = aa.split('.')[0]
    subdirs = '/'.join([aa[i:i+3] for i in range(0, len(aaval), 3)])
               
    url = "{0}/{1}/{2}/{3}/{3}_{4}".format(ftpstem, gc, subdirs,
                                           filestem, suffix)
    hashurl = "{0}/{1}/{2}/{3}/{4}".format(ftpstem, gc, subdirs,
                                           filestem, "md5checksums.txt")
    return (url, hashurl)


# Download a remote file to the specified directory
def download_url(url, outfname, timeout):
    """Downloads a remote URL to a local directory.

    This function downloads the contents of the passed URL to the passed
    filename, in buffered chunks
    """
    # Open connection, and get expected filesize
    response = urlopen(url, timeout=timeout)
    fsize = int(response.info().get("Content-length"))

    # Define buffer sizes
    bsize = 1048576  # buffer size
    fsize_dl = 0     # bytes downloaded
 
    # Download file
    with open(outfname, "wb") as ofh:
        with tqdm(total=fsize) as pbar:
            while True:
                buffer = response.read(bsize)
                if not buffer:
                    break
                fsize_dl += len(buffer)
                ofh.write(buffer)
                pbar.update(bsize)


# Construct filepaths for downloaded files and their hashes
def construct_output_paths(filestem, suffix, outdir):
    """Construct paths to output files for genome and hash."""
    outfname = os.path.join(outdir, '_'.join([filestem, suffix]))
    outfhash = os.path.join(outdir, '_'.join([filestem, "hashes.txt"]))
    return (outfname, outfhash)




# Download a remote genome from NCBI and its MD5 hash
def retrieve_genome_and_hash(filestem, suffix, ftpstem, outdir, timeout):
    """Download genome contigs and MD5 hash data from NCBI."""
    DLStatus = namedlist("DLStatus",
                         "url hashurl outfname outfhash skipped error")
    skipped = False  # Flag - set True if we skip download for existing file
    error = None     # Text of last-raised error

    # Construct remote URLs and output filenames
    url, hashurl = compile_url(filestem, suffix, ftpstem)
    outfname, outfhash = construct_output_paths(filestem, suffix, outdir)

    # Download the genome sequence and corresponding hash file
    try:
        download_url(url, outfname, timeout)
        download_url(hashurl, outfhash, timeout)
    except:
        error = last_exception()

    return DLStatus(url, hashurl, outfname, outfhash, skipped, error)


# Check the file hash against the downloaded hash
def check_hash(fname, hashfile):
    """Checks the MD5 hash of the passed file against an entry in the
    downloaded NCBI hash file."""
    Hashstatus = namedtuple("Hashstatus", "passed localhash filehash")
    filehash = ""
    passed = False   # Flag - set to True if the hash matches

    # Generate MD5 hash
    localhash = create_hash(fname)

    # Get hash from file
    localfname = os.path.split(fname)[-1]
    filehash = extract_hash(hashfile, localfname)
    
    # Check for match
    if filehash == localhash:
        passed = True

    return Hashstatus(passed, localhash, filehash)


# Extract contigs from a compressed file, using gunzip
def extract_contigs(fname, ename):
    """Extract contents of fname to ename using gunzip"""
    with open(ename, 'w') as efh:
        subprocess.run(['gunzip', '-c', fname],
                        stdout=efh)  # can be subprocess.run


# Using a genomes UID, create class and label text files
def create_labels(classification, filestem):
    """Constructs class and label text from UID classification."""
    class_data = (filestem, classification.genus[0] + '.',
                  classification.species, classification.strain)
    labeltxt = "{0}_genomic\t{1} {2} {3}".format(*class_data)
    classtxt = "{0}_genomic\t{1}".format(filestem, classification.organism)
    
    return (labeltxt, classtxt)


# Create an MD5 hash for the passed genome
def create_hash(fname):
    """Return MD5 hash of the passed file contents."""
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as fhandle:
        for chunk in iter(lambda: fhandle.read(65536), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
    

# Create an MD5 hash for the passed genome
def extract_hash(hashfile, name):
    """Return MD5 hash from file of name:MD5 hashes

    - hashfile             path to file containing name:MD5 pairs
    - name                 name asspcoated with hash
    """
    filehash = None
    with open(hashfile, "r") as hhandle:
        for l in [l.strip().split() for l in hhandle if len(l.strip())]:
            hashfname = os.path.split(l[1])[-1]
            if os.path.split(l[1])[-1] == name:
                filehash = l[0]
    return filehash
    


