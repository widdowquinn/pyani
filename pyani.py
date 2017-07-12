#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pyani.py

This script uses the pyani module to produce ANI analyses and classifications
of prokaryotic genome sequences (draft or complete).

(c) The James Hutton Institute 2016-2017
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

Copyright (c) 2016-2017 The James Hutton Institute

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

import logging
import os
import re
import shutil
import sys
import time
import traceback

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict, namedtuple

from pyani import __version__, download

class PyaniDownloadException(Exception):
    """General exception for downloading"""
    def __init__(self, msg="Error in download subcommand"):
        Exception.__init__(self, msg)


# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# Process command-line
def parse_cmdline(args):
    """Parse command-line arguments for script.
    The script offers a single main parser, with subcommands for the actions:
    
    classify - produce a graph-based classification of each input genome
    """
    # Main parent parser
    parser_main = ArgumentParser(prog="pyani.py",
                                 formatter_class=ArgumentDefaultsHelpFormatter)
    subparsers = parser_main.add_subparsers(title="subcommands",
                                            description="valid subcommands",
                                            help="additional help")

    # A 'common' parser, with shared options for all subcommands
    # Not all commands require input or output directories, but all commands
    # support verbose output, and logging.
    parser_common = ArgumentParser(add_help=False)
    parser_common.add_argument('-l', '--logfile', dest='logfile',
                               action='store', default=None,
                               help='logfile location')
    parser_common.add_argument('-v', '--verbose', action='store_true',
                               dest='verbose', default=False,
                               help='report verbose progress to log')

    # Subcommand parsers
    # Download genomes from NCBI
    parser_download = subparsers.add_parser('download',
                                            parents=[parser_common],
                                            formatter_class=\
                                            ArgumentDefaultsHelpFormatter)
    # Classify pyani output into genomotypes
    parser_classify = subparsers.add_parser('classify',
                                            parents=[parser_common],
                                            formatter_class=\
                                            ArgumentDefaultsHelpFormatter)

    # DOWNLOAD: Genome download options
    # Output directory, positional and required
    parser_download.add_argument(action='store',
                                 dest='outdir', default=None,
                                 help='output directory')
    # Required arguments for NCBI query
    parser_download.add_argument("-t", "--taxon", dest="taxon",
                                 action="store", default=None,
                                 help="NCBI taxonomy IDsr (required)",
                                 required=True)
    parser_download.add_argument("--email", dest="email", required=True,
                                 action="store", default=None,
                                 help="Email associated with NCBI queries " +\
                                 "(required)")
    # Arguments controlling connection to NCBI for download
    parser_download.add_argument("--retries", dest="retries",
                                 action="store", default=20,
                                 help="Number of Entrez retry attempts per " +\
                                 "request")
    parser_download.add_argument("--batchsize", dest="batchsize",
                                 action="store", default=10000,
                                 help="Entrez record return batch size")
    parser_download.add_argument("--timeout", dest="timeout",
                                 action="store", default=10,
                                 help="Timeout for URL connection (s)")
    # Arguments controlling local filehandling
    parser_download.add_argument("-f", "--force", dest="force",
                                 action="store_true", default=False,
                                 help="Allow download to existing directory")
    parser_download.add_argument("--noclobber", dest="noclobber",
                                 action="store_true", default=False,
                                 help="Don't replace existing files")
    # Names for output files
    parser_download.add_argument("--labels", dest="labelfname",
                                 action="store", default="labels.txt",
                                 help="Filename for labels file")
    parser_download.add_argument("--classes", dest="classfname",
                                 action="store", default="classes.txt",
                                 help="Filename for classes file")

    # CLASSIFY: Genome classification options
    # Input directory, positional and required
    parser_classify.add_argument(action='store',
                                 dest='indir', default=None,
                                 help='input directory')
    # Output directory, defaults to input directory indir
    parser_classify.add_argument('-o', '--outdir', action='store',
                                 dest='outdir', default=None,
                                 help='output directory')
    # Label file, defaults to indir/labels.txt
    parser_classify.add_argument('--labels', dest='labelfile',
                                 action='store', default=None,
                                 help='file with labels for input genomes')
    # Parameters for classification: minimum %coverage, %identity,
    # and the resolution of thresholds to test
    parser_classify.add_argument('--cov_min', dest='cov_min',
                                 action='store', type=float, default=0.5,
                                 help='minimum %%coverage for an edge')
    parser_classify.add_argument('--id_min', dest='id_min',
                                 action='store', type=float, default=0.8,
                                 help='minimum %%identity for an edge')
    parser_classify.add_argument('--resolution', dest='resolution',
                                 action='store', type=int, default=1500,
                                 help='number of identity thresholds to test')
    
    # Parse arguments
    return parser_main


def download_genome_and_hash(filestem, suffix, ftpstem, outdir, timeout,
                             logger):
    """Download genome and accompanying MD5 hash from NCBI.

    This function tries the (assumed to be passed) RefSeq FTP URL first and,
    if that fails, then attempts to download the corresponding GenBank data.

    We attempt to gracefully skip genomes with download errors.
    """
    # First attempt: RefSeq download
    dlstatus = download.retrieve_genome_and_hash(filestem, suffix,
                                                 ftpstem, outdir, timeout)
    if dlstatus.error is not None:  # Something went awry
        logger.warning("RefSeq download failed: skipping!\n%s", dlstatus.error)
        # Second attempt: GenBank download
        logger.warning("Trying GenBank alternative assembly")
        gbfilestem = re.sub('^GCF_', 'GCA_', filestem)
        logger.info("Retrieving URLs for %s", gbfilestem)
        gbdlstatus = download.retrieve_genome_and_hash(gbfilestem, suffix,
                                                       ftpstem, outdir,
                                                       timeout)
        if gbdlstatus.error:  # Something went awry again
            logger.error("GenBank download failed: skipping!\n%s",
                         gbdlstatus.error)
            dlstatus = gbdlstatus
            dlstatus.skipped = True

    return dlstatus


def write_classes_labels(classes, labels, outdir, classfname, labelfname,
                         noclobber, logger):
    """Write classes and labels files for the downloads."""
    # Write classes
    classfname = os.path.join(outdir, classfname)
    logger.info("Writing classes file to %s", classfname)
    if os.path.exists(classfname) and noclobber:
        logger.warning("Class file %s exists, not overwriting", classfname)
    else:
        with open(classfname, "w") as cfh:
            cfh.write('\n'.join(classes) + '\n')

    # Write labels
    labelfname = os.path.join(outdir, labelfname)
    logger.info("Writing labels file to %s", labelfname) 
    if os.path.exists(labelfname) and noclobber:
        logger.warning("Label file %s exists, not overwriting", labelfname)
    else:
        with open(labelfname, "w") as lfh:
            lfh.write('\n'.join(labels) + '\n')     


def make_outdir(outdir, force, noclobber, logger):
    """Create output directory (allows for force and noclobber).

    The intended outcomes are:
    outdir doesn't exist: create outdir
    outdir exists: raise exception
    outdir exists, --force only: remove the directory tree
    outdir exists, --force --noclobber: continue with existing directory tree
                                        but do not overwrite files

    So long as the outdir is created with this function, we need only check
    for args.noclobber elsewhere to see how to proceed when a file exists.
    """
    if os.path.isdir(outdir):
        logger.warning("Output directory %s exists", outdir)
        if not force:
            raise PyaniDownloadException("Will not modify existing " +
                                         "directory %s", outdir)
        elif force and not noclobber:
            # Delete old directory and start again
            logger.warning("Overwrite forced. Removing %s and everything " +
                           "below it (--force)", outdir)
            shutil.rmtree(outdir)
        else:
            logger.warning("Keeping existing directory, skipping existing " +
                           "files (--force --noclobber).")
    os.makedirs(outdir, exist_ok=True)
    

def make_asm_dict(taxon_ids, retries):
    """Return a dict of assembly UIDs, keyed by each passed taxon ID."""
    asm_dict = dict()

    for tid in taxon_ids:
        asm_uids = download.get_asm_uids(tid, retries)
        asm_dict[tid] = asm_uids.asm_ids

    return asm_dict



# DOWNLOAD
# Download sequence/class/label data from NCBI
def subcmd_download(args, logger):
    """Download all assembled genomes beneath a passed taxon ID from NCBI."""
    # Create output directory
    make_outdir(args.outdir, args.force, args.noclobber, logger)
    
    # Set Entrez email
    download.set_ncbi_email(args.email)
    logger.info("Setting Entrez email address: %s", args.email)
    
    # Get list of taxon IDs to download
    taxon_ids = download.split_taxa(args.taxon)
    logger.info("Taxon IDs received: %s", taxon_ids)

    # Get assembly UIDs for each taxon
    asm_dict = make_asm_dict(taxon_ids, args.retries)
    for tid, uids in asm_dict.items():
        logger.info("Taxon ID summary\n\tQuery: " +\
                    "%s\n\tasm count: %s\n\tUIDs: %s", tid, len(uids), uids)


    # Compile list of outputs for class and label files, and a list of
    # skipped downloads (and a helper tuple for collating skipped genome
    # information)
    classes = []
    labels = []
    skippedlist = []
    Skipped = namedtuple("Skipped",
                         "taxon_id accession organism strain " +
                         "refseq_url genbank_url")

    # Download contigs and hashes for each assembly UID
    for tid, uids in asm_dict.items():
        logger.info("Downloading contigs for Taxon ID %s", tid)
        for uid in uids:
            # Obtain eSummary            
            logger.info("Get eSummary information for UID %s", uid)
            esummary, filestem = download.get_ncbi_esummary(uid, args.retries)
            uid_class = download.get_ncbi_classification(esummary)

            # Report summary
            outstr = '\n\t'.join(["Taxid: %s" % esummary['SpeciesTaxid'],
                                  "Accession: %s" %
                                  esummary['AssemblyAccession'],
                                  "Name: %s" % esummary['AssemblyName'],
                                  "Organism: %s" % uid_class.organism,
                                  "Genus: %s" % uid_class.genus,
                                  "Species: %s" % uid_class.species,
                                  "Strain: %s" % uid_class.strain])
            logger.info("eSummary information:\n\t%s", outstr)

            # Make label/class text
            labeltxt, classtxt = download.create_labels(uid_class, filestem)
            classes.append(classtxt)
            labels.append(labeltxt)
            logger.info("Label and class file entries\n" + 
                        "\tLabel: %s\n\tClass: %s", labeltxt, classtxt)
    
            # Obtain URLs
            ftpstem="ftp://ftp.ncbi.nlm.nih.gov/genomes/all"
            suffix="genomic.fna.gz"
            logger.info("Retrieving URLs for %s", filestem)
            dlstatus = download_genome_and_hash(filestem, suffix, ftpstem,
                                                args.outdir, args.timeout,
                                                logger)
            if dlstatus.skipped:
                skippedlist.append(Skipped(tid, uid,
                                           uid_class.organism,
                                           uid_class.strain,
                                           dlstatus.url, gbdlstatus.url))
                continue  # Move straight on to the next download

            # Report the working download
            logger.info("Downloaded from URL: %s", dlstatus.url)
            logger.info("Wrote assembly to: %s", dlstatus.outfname)
            logger.info("Wrote MD5 hashes to: %s", dlstatus.outfhash)

            # Check hash for the download
            hashstatus = download.check_hash(dlstatus.outfname,
                                             dlstatus.outfhash)
            logger.info("Local MD5 hash: %s", hashstatus.localhash)
            logger.info("NCBI MD5 hash: %s", hashstatus.filehash)
            if hashstatus.passed:
                logger.info("MD5 hash check passed")
            else:
                logger.warning("MD5 hash check failed.")

            # Extract downloaded files
            ename = os.path.splitext(dlstatus.outfname)[0]
            if os.path.exists(ename) and args.noclobber:
                logger.warning("Output file %s exists, not extracting", ename)
            else:
                logger.info("Extracting archive %s to %s",
                            dlstatus.outfname, ename)
                download.extract_contigs(dlstatus.outfname, ename)
        
    # Write class and label files
    write_classes_labels(classes, labels, args.outdir, args.classfname,
                         args.labelfname, args.noclobber, logger)

    # Show skipped genome list
    if len(skippedlist):
        logger.warning("%d genome downloads were skipped", len(skippedlist))
        for skipped in skippedlist:
            outstr = '\n\t'.join(["taxon id: %s" % skipped.taxon_id,
                                   "accession: %s" % skipped.accession,
                                   "RefSeq URL: %s" % skipped.refseq_url,
                                   "GenBank URL: %s" % skipped.genbank_url])
            logger.warning("%s %s:\n\t%s", skipped.organism, skipped.strain,
                           outstr)


# CLASSIFY
# Classify input genomes on basis of ANI coverage and identity output
def subcmd_classify(args, logger):
    """Take pyani output, and generate a series of classifications of the
    input data.
    """
    raise NotImplementedError


# Set up logger
def build_logger(args):
    """Return a logger for the script."""
    # Set up logger
    logger = logging.getLogger('pyani.py: %s' % time.asctime())
    t0 = time.time()

    # Default logging level is DEBUG
    logger.setLevel(logging.DEBUG)

    # Create an error handler in STDERR
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)    

    # STDERR gets INFO and above when verbose, WARNING and above otherwise
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)

    # Was a logfile specified in the passed args? If so, open a stream
    # to this file for logging, and always record at INFO and above
    if args.logfile is not None:
        logstream = open(args.logfile, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)

    return logger


# MAIN
# The pyani main script function
def pyani_main(scriptargs=None):
    """The main routine."""
    # Get arguments from command-line if none provided directly
    if scriptargs is None:
        scriptargs = sys.argv

    # Parse command-line
    parser = parse_cmdline(scriptargs)
    args = parser.parse_args()

    # Get a logger - if we can't, drop out.
    try:
        logger = build_logger(args)
    except:
        sys.stderr.write("Could not create logger (exiting)")
        sys.stderr.write(last_exception())
        sys.exit(1)

    # If no arguments provided, show usage and drop out
    if len(sys.argv) == 1:
        print("pyani version: {0}".format(__version__))
        parser.print_help()
        sys.exit(1)

    # Report arguments, if verbose
    logger.info('Processed arguments: %s' % args)
    logger.info('command-line: %s' % ' '.join(sys.argv))

    # Define subcommand main functions, and distribute on basis of subcommand
    # NOTE: Halting errors/exceptions raised during running subcommands
    #       bubble up to be caught and logged here - they do not generally
    #       need to be caught before this point.
    subcmd = sys.argv[1]
    subcmds = {'classify': subcmd_classify,
               'download': subcmd_download}
    try:
        subcmds[subcmd](args, logger)
    except KeyError:
        logger.error("Subcommand %s not recognised (exiting)", subcmd)
        sys.exit(1)
    except NotImplementedError:
        logger.error("Subcommand %s not yet implemented (exiting)", subcmd)
        sys.exit(1)
    except:
        logger.error("Could not execute subcommand %s", subcmd)
        logger.error(last_exception())
        sys.exit(1)

    # Let the user know we're done
    logger.info(time.asctime())
    logger.info("Done.")    

    # Exit cleanly (POSIX)
    sys.exit(0)

    
###
# Run as script
if __name__ == '__main__':
    # Call main function
    pyani_main()

