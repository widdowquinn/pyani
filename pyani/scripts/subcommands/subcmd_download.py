#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""subcmd_download.py

Provides the download subcommand for pyani

(c) The James Hutton Institute 2017-18

Author: Leighton Pritchard
Contact: leighton.pritchard@hutton.ac.uk

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

Copyright (c) 2017-18 The James Hutton Institute

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

from collections import namedtuple

from Bio import SeqIO

from pyani import download
from pyani.scripts import tools

# Convenience struct for file download data
DLFileData = namedtuple("DLFileData", "filestem ftpstem suffix")


def subcmd_download(args, logger):
    """Download assembled genomes in subtree of passed NCBI taxon ID."""
    # Create output directory, respecting force/noclobber
    if not args.dryrun:
        tools.make_outdir(args.outdir, args.force, args.noclobber, logger)
    else:
        logger.warning("Dry run only: will not overwrite or download")

    # Set Entrez email
    download.set_ncbi_email(args.email)
    logger.info("Setting Entrez email address: %s", args.email)

    # Get list of taxon IDs to download
    taxon_ids = download.split_taxa(args.taxon)
    logger.info("Taxon IDs received: %s", taxon_ids)

    # Get assembly UIDs for each taxon
    asm_dict = tools.make_asm_dict(taxon_ids, args.retries)
    for tid, uids in asm_dict.items():
        logger.info(
            "Taxon ID summary\n\tQuery: " + "%s\n\tasm count: %s\n\tUIDs: %s",
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
    Skipped = namedtuple(
        "Skipped", "taxon_id accession organism strain " + "url dltype"
    )

    # Download contigs and hashes for each assembly UID in the list
    # On completion of this loop, each assembly in the list will either be
    # downloaded or skipped (with skipped genome information preserved in
    # skippedlist), and class/label info will be collated, ready for writing
    # to file.
    for tid, uids in asm_dict.items():
        logger.info("Downloading contigs for Taxon ID %s", tid)
        for uid in uids:
            # Obtain eSummary
            logger.info("Get eSummary information for UID %s", uid)
            esummary, filestem = download.get_ncbi_esummary(uid, args.retries)
            uid_class = download.get_ncbi_classification(esummary)

            # Report summary
            outstr = "\n\t".join(
                [
                    "Species Taxid: %s" % esummary["SpeciesTaxid"],
                    "TaxID: %s" % esummary["Taxid"],
                    "Accession: %s" % esummary["AssemblyAccession"],
                    "Name: %s" % esummary["AssemblyName"],
                    "Organism: %s" % uid_class.organism,
                    "Genus: %s" % uid_class.genus,
                    "Species: %s" % uid_class.species,
                    "Strain: %s" % uid_class.strain,
                ]
            )
            logger.info("eSummary information:\n\t%s", outstr)
            if args.dryrun:
                logger.warning(
                    "(dry-run) skipping download of %s", esummary["AssemblyAccession"]
                )
                continue

            # Obtain URLs, trying the RefSeq filestem first, then GenBank if
            # there's a failure
            dlfiledata = DLFileData(
                filestem, "ftp://ftp.ncbi.nlm.nih.gov/genomes/all", "genomic.fna.gz"
            )
            logger.info("Retrieving URLs for %s", filestem)
            # Try RefSeq first
            dlstatus = tools.download_genome_and_hash(
                args,
                logger,
                dlfiledata,
                dltype="RefSeq",
                disable_tqdm=args.disable_tqdm,
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
                logger.warning(
                    "RefSeq failed. Trying GenBank alternative " + "assembly"
                )
                # Try GenBank assembly
                dlstatus = tools.download_genome_and_hash(
                    args,
                    logger,
                    dlfiledata,
                    dltype="GenBank",
                    disable_tqdm=args.disable_tqdm,
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
            logger.info("Downloaded from URL: %s", dlstatus.url)
            logger.info("Wrote assembly to: %s", dlstatus.outfname)
            logger.info("Wrote MD5 hashes to: %s", dlstatus.outfhash)

            # Check hash for the download
            hashstatus = download.check_hash(dlstatus.outfname, dlstatus.outfhash)
            logger.info("Local MD5 hash: %s", hashstatus.localhash)
            logger.info("NCBI MD5 hash: %s", hashstatus.filehash)
            if hashstatus.passed:
                logger.info("MD5 hash check passed")
            else:
                logger.warning("MD5 hash check failed. " + "Please check and retry.")

            # Extract downloaded files
            ename = os.path.splitext(dlstatus.outfname)[0]
            if os.path.exists(ename) and args.noclobber:
                logger.warning("Output file %s exists, not extracting", ename)
            else:
                logger.info("Extracting archive %s to %s", dlstatus.outfname, ename)
                download.extract_contigs(dlstatus.outfname, ename)

            # Modify sequence ID header if Kraken option active
            if args.kraken:
                logger.warning("Modifying downloaded sequence for Kraken compatibility")
                seqdata = list(SeqIO.parse(ename, "fasta"))
                logger.info("Modifying %s", ename)
                for seq in seqdata:
                    seq.id = "|".join(
                        [seq.id, "kraken:taxid", esummary["SpeciesTaxid"]]
                    )
                SeqIO.write(seqdata, ename, "fasta")

            # Create MD5 hash for the downloaded contigs
            logger.info("Creating local MD5 hash for %s" % ename)
            hashfname = os.path.splitext(ename)[0] + ".md5"
            datahash = download.create_hash(ename)
            logger.info("Writing hash to %s" % hashfname)
            with open(hashfname, "w") as hfh:
                hfh.write("\t".join([datahash, ename]) + "\n")
            # Make label/class text
            labeltxt, classtxt = download.create_labels(uid_class, filestem, datahash)
            classes.append(classtxt)
            labels.append(labeltxt)
            logger.info(
                "Label and class file entries\n" + "\tLabel: %s\n\tClass: %s",
                labeltxt,
                classtxt,
            )

    # Write class and label files
    classfname = os.path.join(args.outdir, args.classfname)
    logger.info("Writing classes file to %s", classfname)
    if os.path.exists(classfname) and args.noclobber:
        logger.warning("Class file %s exists, not overwriting", classfname)
    else:
        with open(classfname, "w") as ofh:
            ofh.write("\n".join(classes) + "\n")

    labelfname = os.path.join(args.outdir, args.labelfname)
    logger.info("Writing labels file to %s", labelfname)
    if os.path.exists(labelfname) and args.noclobber:
        logger.warning("Labels file %s exists, not overwriting", labelfname)
    else:
        with open(labelfname, "w") as ofh:
            ofh.write("\n".join(labels) + "\n")

    # Report skipped genome list
    if skippedlist:
        logger.warning("%d genome downloads were skipped", len(skippedlist))
        for skipped in skippedlist:
            outstr = "\n\t".join(
                [
                    "taxon id: %s" % skipped.taxon_id,
                    "accession: %s" % skipped.accession,
                    "URL: %s" % skipped.url,
                    "source: %s" % skipped.dltype,
                ]
            )
            logger.warning("%s %s:\n\t%s", skipped.organism, skipped.strain, outstr)
