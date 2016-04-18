#!/usr/bin/env python3
#
# genbank_get_genomes_by_taxon.py
#
# A script that takes an NCBI taxonomy identifier (or string, though this is
# not reliable for taxonomy tree subgraphs...) and downloads all genomes it 
# can find from NCBI in the corresponding taxon subgraph with the passed
# argument as root.
#
# (c) The James Hutton Institute 2015
# Author: Leighton Pritchard

import logging
import logging.handlers
import os
import shutil
import subprocess
import sys
import time
import traceback
from urllib.request import urlopen

from argparse import ArgumentParser
from collections import defaultdict

from Bio import Entrez, SeqIO

# Parse command-line
def parse_cmdline(args):
    """Parse command-line arguments"""
    parser = ArgumentParser(prog="genbacnk_get_genomes_by_taxon.py")
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
    logger.info("Set NCBI contact email to %s" % args.email)
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
            logger.error("Output directory %s would " % args.outdirname +
                         "overwrite existing files (exiting)")
            sys.exit(1)
        else:
            logger.info("--force output directory use")
            if args.noclobber:
                logger.warning("--noclobber: existing output directory kept")
            else:
                logger.info("Removing directory %s and everything below it" %
                            args.outdirname)
                shutil.rmtree(args.outdirname)
    logger.info("Creating directory %s" % args.outdirname)
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
def entrez_retry(fn, *fnargs, **fnkwargs):
    """Retries the passed function up to the number of times specified
    by args.retries
    """
    tries, success = 0, False
    while not success and tries < args.retries:
        try:
            output = fn(*fnargs, **fnkwargs)
            success = True
        except:
            tries += 1
            logger.warning("Entrez query %s(%s, %s) failed (%d/%d)" %
                           (fn, fnargs, fnkwargs, tries+1, args.retries))
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
    logger.info("Entrez ESearch with query: %s" % query)
    
    # Perform initial search for assembly UIDs with taxon ID as query.
    # Use NCBI history for the search.
    handle = entrez_retry(Entrez.esearch, db="assembly", term=query,
                          format="xml", usehistory="y")
    record = Entrez.read(handle, validate=False)
    result_count = int(record['Count'])
    logger.info("Entrez ESearch returns %d assembly IDs" % result_count)
    
    # Recover assembly UIDs from the web history
    asm_ids = entrez_batch_webhistory(record, result_count, 250,
                                      db="assembly", retmode="xml")
    logger.info("Identified %d unique assemblies" % len(asm_ids))
    return asm_ids


# Get contig UIDs for a specified assembly UID
def get_contig_uids(asm_uid):
    """Returns a set of NCBI UIDs for each contig belonging to the assembly
    with passsed UID.

    Assemblies are linked to contigs through Elink in three ways:
    assembly_nuccore_insdc; assembly_nuccore_refseq; and
    assembly_nuccore_wgsmaster. For most insdc and refseq links, we can
    use the results of an Elink query directly, recovering contigs with
    a (history/batch) Entrez ESearch. For wgsmaster-only links, we need
    to download from a URL that is not given in the Esummary, and which we
    have to munge together from information in the summary.

    For insdc/refseq Elink queries that report 100000 contig links, this is
    a capped return value result, and forces us to download contigs as though
    we use the wgsmaster links.
    """
    logger.info("Finding contig UID links for assembly %s" % asm_uid)
    # The Elink search returns a single result, so we don't need to batch
    linklist = entrez_retry(Entrez.elink, dbfrom="assembly", db="nucleotide",
                            retmode="gb", from_uid=asm_uid)
    links = Entrez.read(linklist, validate=False) 
    # Assemblies may be in 'assembly_nuccore_insdc', 'nuccore', or
    # 'assembly_nuccore_wgsmaster' databases. For a list of link names:
    # http://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html
    # We prefer insdc > refseq > wgsmaster, but if there are 100000 contig 
    # UIDs returned, this indicates that not all contigs will be downloaded,
    # so we revert to wgsmaster
    try:
        wgsmaster, fastafilename = False, None
        linktype = guess_linktype(links)
        logger.info("Using %s links" % linktype)
        if linktype in ['assembly_nuccore_insdc', 'assembly_nuccore_inrefseq']:
            contig_uids = retrieve_contig_uids(links, linktype)
            logger.info("Recovered %d contig UIDs" % len(contig_uids))
            if len(contig_uids) == 100000:
                logger.info("insdc/refseq result cap: reverting to wgsmaster")
                wgsmaster = True
                contig_uids, fastafilename = download_wgsmaster_contigs(links)
            else:
                write_contigs(asm_uid, contig_uids)
        if linktype == 'assembly_nuccore_wgsmaster':
            wgsmaster = True
            contig_uids, fastafilename = download_wgsmaster_contigs(links)
    except:
        logger.error("Could not recover contigs (exiting)")
        logger.error(last_exception())
        sys.exit(1)
    logger.info("Identified %d contig UIDs" % len(contig_uids))
    return {'wgsmaster': wgsmaster,
            'filename': fastafilename,
            'contig_uids': contig_uids}

# Download contigs from wgsmaster record indicated in Elink result
def download_wgsmaster_contigs(links):
    """Download assembly contigs from the wgsmaster record in the passed
    Elink results, and extract the archive. Returns a list of IDs for the
    contigs, and path to the extracted FASTA file.
    """
    wgsmaster_uid = [l['Link'][0]['Id'] for l in links[0]['LinkSetDb']
                     if l['LinkName']=='assembly_nuccore_wgsmaster'][0]
    return retrieve_wgsmaster_contigs(wgsmaster_uid)

# Retrieves contigs from Entrez by links from assembly entries
def retrieve_contig_uids(links, linktype):
    """Passed the returned 'assembly_nuccore_insdc' or
    'assembly_nuccore_refseq' links for an assembly, retrieves the contig
    UIDs.

    - links: returned links from an Elink query
    - linktype: the kind of links we want to use for contig UIDs
    """
    contigs = [l for l in links[0]['LinkSetDb'] 
               if l['LinkName'] == linktype][0]
    contig_uids = set([e['Id'] for e in contigs['Link']])
    return contig_uids

# Returns the linktype to be used for downloading contigs
def guess_linktype(links):
    """Returns the type of link to be used for downloading contigs, from an
    Entrez Elink list. We prefer insdc > refseq > wgsmaster
    """
    for linktype in ['assembly_nuccore_insdc', 'assembly_nuccore_inrefseq',
                     'assembly_nuccore_wgsmaster']:
        if linktype in [l['LinkName'] for l in links[0]['LinkSetDb']]:
            return linktype
    # If we're here, then there's no valid linktype
    logger.error("No valid linktype found (exiting)")
    sys.exit(1)

# Download and uncompress the wgsmaster contig archive
def retrieve_wgsmaster_contigs(uid):
    """Munges a download URL from the passed UID and downloads the
    corresponding archive from NCBI, extracting it to the output
    directory.
    """
    logger.info("Processing wgsmaster UID: %s" % uid)
    summary = Entrez.read(Entrez.esummary(db='nuccore', id=uid,
                                          rettype='text', validate=False))
    # Assume that the 'Extra' field is present and is well-formatted.
    # Which means that the first six characters of the last part of
    # the 'Extra' string correspond to the download archive filestem.
    dlstem = summary[0]['Extra'].split('|')[-1][:6]
    dlver = summary[0]['Extra'].split('|')[3].split('.')[-1]
    # Download archive to output directory
    # Establish download size; if version number not in sync with
    # download, try again with version number decremented by 1. This 
    # may be necessary because genome sequence version and genome/
    # assembly version numbers are not synchronised. 
    fsize = None
    while str(dlver) != '0' and not fsize:
        try:
            fname = "%s.%s.fsa_nt.gz" % (dlstem, dlver)
            outfname = os.path.join(args.outdirname, fname)
            url = "http://www.ncbi.nlm.nih.gov/Traces/wgs/?download=%s" % \
                fname
            logger.info("Trying URL: %s" % url)
            response = urlopen(url)
            meta = response.info()
            fsize = int(meta.getheaders("Content-length")[0])
            logger.info("Downloading: %s Bytes: %s" % (fname, fsize))
            fsize_dl = 0
            bsize = 1048576
        except:  # Download didn't work. Assuming it's because of version
            fsize = None
            logger.error("Download failed for (%s)" % url)
            if str(dlver) != '0':
                dlver = int(dlver) - 1
                logger.info("Retrying download with version = %s" % dlver)
            else:
                logger.error("No more versions to try (exiting)")
                sys.exit(1)
    # Download data
    try:
        with open(outfname, 'wb') as fh:
            while True:
                buffer = response.read(bsize)
                if not buffer:
                    break
                fsize_dl += len(buffer)
                fh.write(buffer)
                status = r"%10d  [%3.2f%%]" % (fsize_dl,
                                               fsize_dl * 100. / fsize)
                logger.info(status)
    except:
        logger.error("Download failed for %s (exiting)" % fname)
        logger.error(last_exception())
        sys.exit(1)
    # Extract archive
    asm_summary = entrez_retry(Entrez.esummary, db='assembly', id=asm_uid,
                               rettype='text')
    asm_record = Entrez.read(asm_summary, validate=False)
    gname = asm_record['DocumentSummarySet']['DocumentSummary']\
            [0]['AssemblyAccession']
    extractfname = os.path.join(args.outdirname,
                                '.'.join([gname, 'fasta']))
    try:
        logger.info("Extracting archive %s to %s" %
                    (outfname, extractfname))
        with open(extractfname, 'w') as efh:
            subprocess.call(['gunzip', '-c', outfname],
                            stdout=efh)  # can be subprocess.run in Py3.5
        logger.info("Archive extracted to %s" % extractfname)
    except:
        logger.error("Extracting archive %s failed (exiting)" % outfname)
        logger.error(last_exception())
        sys.exit(1)
    # Get contig_uids
    contig_uids = [s.description for s in SeqIO.parse(extractfname, 'fasta')]
    return contig_uids, extractfname


def get_class_label_info(asm_uid):
    """Returns class and label strings for an assembly indicated by asm_uid."""
    # Currently duplicated in write_contigs() - needs refactoring
    logger.info("Recovering class and label info for assembly %s" % asm_uid)
    # Assembly record - get binomial and strain names
    asm_summary = entrez_retry(Entrez.esummary, db='assembly', id=asm_uid,
                               rettype='text')
    asm_record = Entrez.read(asm_summary, validate=False)
    asm_organism = asm_record['DocumentSummarySet']['DocumentSummary']\
                   [0]['SpeciesName']
    try:
        asm_strain = asm_record['DocumentSummarySet']['DocumentSummary']\
                     [0]['Biosource']['InfraspeciesList'][0]['Sub_value']
    except:
        asm_strain = ""
    # Assembly UID (long form) for the output filename 
    gname = asm_record['DocumentSummarySet']['DocumentSummary']\
            [0]['AssemblyAccession']
    outfilename = "%s.fasta" % os.path.join(args.outdirname, gname)

    # Create label and class strings
    genus, species = asm_organism.split(' ', 1)
    ginit = genus[0] + '.'
    labeltxt = "%s\t%s %s %s" % (gname, ginit, species, asm_strain)
    classtxt = "%s\t%s" % (gname, asm_organism)
    logger.info("UID: %s Label: %s" % (asm_uid, labeltxt))
    logger.info("UID: %s Class: %s" % (asm_uid, classtxt))

    # Return labels
    return classtxt, labeltxt


# Write contigs for a single assembly out to file
def write_contigs(asm_uid, contig_uids):
    """Writes assembly contigs out to a single FASTA file in the script's
    designated output directory.

    FASTA records are returned, as GenBank and even GenBankWithParts format
    records don't reliably give correct sequence in all cases.

    The script returns two strings for each assembly, a 'class' and a 'label'
    string - this is for use with, e.g. pyani.
    """
    # Has duplicate code with get_class_label_info() - needs refactoring
    logger.info("Collecting contig data for %s" % asm_uid)
    # Assembly record - get binomial and strain names
    asm_summary = entrez_retry(Entrez.esummary, db='assembly', id=asm_uid,
                               rettype='text')
    asm_record = Entrez.read(asm_summary, validate=False)
    asm_organism = asm_record['DocumentSummarySet']['DocumentSummary']\
                   [0]['SpeciesName']
    try:
        asm_strain = asm_record['DocumentSummarySet']['DocumentSummary']\
                     [0]['Biosource']['InfraspeciesList'][0]['Sub_value']
    except:
        asm_strain = ""
    # Assembly UID (long form) for the output filename
    gname = asm_record['DocumentSummarySet']['DocumentSummary']\
            [0]['AssemblyAccession']
    outfilename = "%s.fasta" % os.path.join(args.outdirname, gname)

    # Create label and class strings
    genus, species = asm_organism.split(' ', 1)
    ginit = genus[0] + '.'
    labeltxt = "%s\t%s %s %s" % (gname, ginit, species, asm_strain)
    classtxt = "%s\t%s" % (gname, asm_organism)

    # Get FASTA records for contigs
    logger.info("Downloading FASTA records for assembly %s (%s)" %
                (asm_uid, ' '.join([ginit, species, asm_strain])))
    # We're doing an explicit retry loop here because we want to confirm we
    # have the correct data, as well as test for Entrez connection errors,
    # which is all the entrez_retry function does.
    tries, success = 0, False
    while not success and tries < args.retries:
        try:
            records = []  # Holds all return records
            # We may need to batch contigs
            query_uids = ','.join(contig_uids)
            batch_size = 10000
            for start in range(0, len(contig_uids), batch_size):
                logger.info("Batch: %d-%d" % (start, start+batch_size))
                seqdata = entrez_retry(Entrez.efetch, db='nucleotide',
                                       id=query_uids,
                                       rettype='fasta', retmode='text',
                                       retstart=start, retmax=batch_size)
                records.extend(list(SeqIO.parse(seqdata, 'fasta')))
            tries += 1
            # Check only that correct number of records returned.
            if len(records) == len(contig_uids):  
                success = True
            else:  
                logger.warning("%d contigs expected, %d contigs returned" %
                               (len(contig_uids), len(records)))
                logger.warning("FASTA download for assembly %s failed" %
                               asm_uid)
                logger.warning("try %d/20" % tries)
            # Could also check expected assembly sequence length?
            totlen = sum([len(r) for r in records])
            logger.info("Downloaded genome size: %d" % totlen)
        except:
            logger.warning("FASTA download for assembly %s failed" % asm_uid)
            logger.warning(last_exception())
            logger.warning("try %d/20" % tries)            
    if not success:
        # Could place option on command-line to stop or continue here.
        #logger.error("Failed to download records for %s (exiting)" % asm_uid)
        #sys.exit(1)
        logger.error("Failed to download records for %s (continuing)" % asm_uid)

    # Write contigs to file
    retval = SeqIO.write(records, outfilename, 'fasta')
    logger.info("Wrote %d contigs to %s" % (retval, outfilename))


# Run as script
if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline(sys.argv)

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
            logger.error("Could not open %s for logging" %
                         args.logfile)
            sys.exit(1)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info("genbank_get_genomes_by_taxon.py: %s" % time.asctime())
    logger.info("command-line: %s" % ' '.join(sys.argv))
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
    logger.info("Output directory: %s" % args.outdirname)

    # We might have more than one taxon in a comma-separated list
    taxon_ids = args.taxon.split(',')
    logger.info("Passed taxon IDs: %s" % ', '.join(taxon_ids))    

    # Get all NCBI assemblies for each taxon UID
    asm_dict = defaultdict(set)
    for tid in taxon_ids:
        asm_dict[tid] = get_asm_uids(tid)
    for tid, asm_uids in list(asm_dict.items()):
        logger.info("Taxon %s: %d assemblies" % (tid, len(asm_uids)))

    # Download contigs for each assembly UID
    classes, labels = [], []
    contig_dict = defaultdict(set)
    for tid, asm_uids in list(asm_dict.items()):
        for asm_uid in asm_uids:
            contig_dict[asm_uid] = get_contig_uids(asm_uid)
            classtxt, labeltxt = get_class_label_info(asm_uid)
            classes.append(classtxt)
            labels.append(labeltxt)

    # Write class and label files
    classfilename = os.path.join(args.outdirname, 'classes.txt')
    labelfilename = os.path.join(args.outdirname, 'labels.txt')
    logger.info("Writing classes file to %s" % classfilename)
    with open(classfilename, 'w') as fh:
        fh.write('\n'.join(classes))
    logger.info("Writing labels file to %s" % labelfilename)
    with open(labelfilename, 'w') as fh:
        fh.write('\n'.join(labels))

    # Report the number of contigs we got for each assembly
    for asm_uid, contig_uids in list(contig_dict.items()):
        logger.info("Assembly %s: %d contigs" %
                    (asm_uid, len(contig_uids['contig_uids'])))

    # Let the user know we're done
    logger.info(time.asctime())
    logger.info("Done.")
