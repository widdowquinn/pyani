#!/usr/bin/env python3
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
# 161 Cathedral Street,
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
"""Provides the download subcommand for pyani."""

import logging

from argparse import Namespace
from collections import namedtuple

from Bio import SeqIO

from pyani import download
from pyani.pyani_tools import termcolor
from pyani.scripts import tools


def subcmd_download(args: Namespace) -> int:
    """Download assembled genomes in subtree of passed NCBI taxon ID.

    :param args:  Namespace, command-line arguments
    """
    # Create logger
    logger = logging.getLogger(__name__)
    logger.info(termcolor("Downloading genomes from NCBI", "red"))

    # Create output directory, respecting force/noclobber
    if not args.dryrun:
        tools.make_outdir(args.outdir, args.force, args.noclobber)
    else:
        logger.warning(
            termcolor("Dry run only: will not overwrite or download", "cyan")
        )

    # Set Entrez email
    download.set_ncbi_email(args.email)
    logger.info("Setting Entrez email address: %s", args.email)

    # Parse Entrez API key, if provided
    api_path = args.api_keypath.expanduser()
    if not api_path.is_file():
        logger.warning("API path %s not a valid file. Not using API key.", api_path)
        api_key = None
    else:
        api_key = download.parse_api_key(api_path)
        logger.info("API key recovered from %s", api_path)

    # Get list of taxon IDs to download
    taxon_ids = download.split_taxa(args.taxon)
    logger.info(termcolor("Taxon IDs received: %s", "blue"), taxon_ids)

    # Get assembly UIDs for each taxon
    asm_dict = tools.make_asm_dict(taxon_ids, args.retries)
    for tid, uids in asm_dict.items():
        logger.debug(
            "Taxon ID summary\n\tQuery: %s\n\tasm count: %s\n\tUIDs: %s",
            tid,
            len(uids),
            uids,
        )

    # Compile outputs to write class and label files, and a list of
    # skipped downloads (and define a helper tuple for collating skipped
    # genome information)
    classes = []
    labels = []
    skippedlist = []
    Skipped = namedtuple("Skipped", "taxon_id accession organism strain url dltype")

    # Download contigs and hashes for each assembly UID in the list
    # On completion of this loop, each assembly in the list will either be
    # downloaded or skipped (with skipped genome information preserved in
    # skippedlist), and class/label info will be collated, ready for writing
    # to file.
    # Summary information is reported to the logger for each eSummary that
    # can be recovered
    for tid, uids in asm_dict.items():
        logger.info(termcolor("Downloading contigs for Taxon ID %s", "blue"), uids)
        for uid in uids:
            # Obtain eSummary
            logger.info(
                termcolor("Retrieving eSummary information for UID %s", "cyan"), uid
            )
            logger.debug(
                "NCBI eSummary:\n%s",
                download.get_ncbi_esummary(uid, args.retries, api_key),
            )
            esummary, filestem = download.get_ncbi_esummary(uid, args.retries, api_key)
            uid_class = download.get_ncbi_classification(esummary)

            # Report summary
            outstr = "\n\t".join(
                [
                    f"Species Taxid: {esummary['SpeciesTaxid']}",
                    f"TaxID: {esummary['Taxid']}",
                    f"Accession: {esummary['AssemblyAccession']}",
                    f"Name: {esummary['AssemblyName']}",
                    f"Organism: {uid_class.organism}",
                    f"Genus: {uid_class.genus}",
                    f"Species: {uid_class.species}",
                    f"Strain: {uid_class.strain}",
                ]
            )
            logger.debug("eSummary information:\n\t%s", outstr)
            if args.dryrun:
                logger.warning(
                    "(dry-run) skipping download of %s", esummary["AssemblyAccession"]
                )
                continue

            # Obtain URLs, trying the RefSeq filestem first, then GenBank if
            # there's a failure
            dlfiledata = tools.DLFileData(
                filestem, "ftp://ftp.ncbi.nlm.nih.gov/genomes/all", "genomic.fna.gz"
            )
            logger.info("Retrieving URLs for %s", filestem)
            # Try RefSeq first
            dlstatus = tools.download_genome_and_hash(
                args, dlfiledata, dltype="RefSeq", disable_tqdm=args.disable_tqdm,
            )
            # RefSeq failed, try GenBank
            # Pylint is confused by the content of dlstatus (a namedlist)
            if dlstatus.skipped:  # pylint: disable=no-member
                skippedlist.append(
                    Skipped(
                        tid,
                        uid,
                        uid_class.organism,
                        uid_class.strain,
                        dlstatus.url,  # pylint: disable=no-member
                        "RefSeq",
                    )
                )
                logger.warning("RefSeq failed. Trying GenBank alternative assembly")
                # Try GenBank assembly
                dlstatus = tools.download_genome_and_hash(
                    args, dlfiledata, dltype="GenBank", disable_tqdm=args.disable_tqdm,
                )
                # Pylint is confused by the content of dlstatus (a namedlist)
                if dlstatus.skipped:  # pylint: disable=no-member
                    skippedlist.append(
                        Skipped(
                            tid,
                            uid,
                            uid_class.organism,
                            uid_class.strain,
                            dlstatus.url,
                            "GenBank",
                        )
                    )
                    logger.warning("GenBank failed.")
                    continue  # Move straight on to the next download

            # One of the downloads worked: report information
            logger.debug("Downloaded from URL: %s", dlstatus.url)
            logger.debug("Wrote assembly to: %s", dlstatus.outfname)
            logger.debug("Wrote MD5 hashes to: %s", dlstatus.outfhash)

            # Check hash for the download
            hashstatus = download.check_hash(dlstatus.outfname, dlstatus.outfhash)
            logger.debug("Local MD5 hash: %s", hashstatus.localhash)
            logger.debug("NCBI MD5 hash: %s", hashstatus.localhash)
            if hashstatus.passed:
                logger.info(termcolor("MD5 hash check passed", "green"))
            else:
                logger.warning("MD5 hash check failed. Please check and retry.")

            # Extract downloaded files
            ename = dlstatus.outfname.with_suffix("")  # should strip only last suffix
            if ename.exists() and args.noclobber:
                logger.warning("Output file %s exists, not extracting", ename)
            else:
                logger.debug("Extracting archive %s to %s", dlstatus.outfname, ename)
                download.extract_contigs(dlstatus.outfname, ename)

            # Modify sequence ID header if Kraken option active
            if args.kraken:
                logger.warning("Modifying downloaded sequence for Kraken compatibility")
                seqdata = list(SeqIO.parse(ename, "fasta"))
                logger.debug("Modifying %s", ename)
                for seq in seqdata:
                    seq.id = "|".join(
                        [seq.id, "kraken:taxid", esummary["SpeciesTaxid"]]
                    )
                SeqIO.write(seqdata, ename, "fasta")

            # Create MD5 hash for the downloaded contigs
            logger.debug("Creating local MD5 hash for %s", ename)
            hashfname = ename.with_suffix(".md5")
            datahash = download.create_hash(ename)
            logger.debug("Writing hash to %s", hashfname)
            with open(hashfname, "w") as hfh:
                hfh.write("\t".join([datahash, str(ename)]) + "\n")
            # Make label/class text
            labeltxt, classtxt = download.create_labels(uid_class, filestem, datahash)
            classes.append(classtxt)
            labels.append(labeltxt)
            logger.info(
                "Label and class file entries\n\tLabel: %s\n\tClass: %s",
                labeltxt,
                classtxt,
            )

    # Write class and label files
    classfname = args.outdir / args.classfname
    logger.info("Writing classes file to %s", classfname)
    if classfname.exists() and args.noclobber:
        logger.warning("Class file %s exists, not overwriting", classfname)
    else:
        with open(classfname, "w") as ofh:
            ofh.write("\n".join(classes) + "\n")

    labelfname = args.outdir / args.labelfname
    logger.info("Writing labels file to %s", labelfname)
    if labelfname.exists() and args.noclobber:
        logger.warning("Labels file %s exists, not overwriting", labelfname)
    else:
        with open(labelfname, "w") as ofh:
            ofh.write("\n".join(labels) + "\n")

    # Report skipped genome list
    if skippedlist:
        logger.warning(
            termcolor("%s genome downloads were skipped", "red"), len(skippedlist)
        )
        for skipped in skippedlist:
            outstr = "\n\t".join(
                [
                    f"taxon id: {skipped.taxon_id}",
                    f"accession: {skipped.accession}",
                    f"URL: {skipped.url}",
                    f"source: {skipped.dltype}",
                ]
            )
            logger.warning("%s %s:\n\t%s", skipped.organism, skipped.strain, outstr)

    return 0
