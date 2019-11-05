#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019
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
# Copyright (c) 2019 University of Strathclyde
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

from collections import namedtuple

from Bio import SeqIO

from pyani import download
from pyani.scripts import tools

# Convenience struct for file download data
DLFileData = namedtuple("DLFileData", "filestem ftpstem suffix")


def subcmd_download(args, logger):
    """Download assembled genomes in subtree of passed NCBI taxon ID.

    :param args:  Namespace, command-line arguments
    :param logger:  logging object
    """
    # Create output directory, respecting force/noclobber
    if not args.dryrun:
        tools.make_outdir(args.outdir, args.force, args.noclobber, logger)
    else:
        logger.warning("Dry run only: will not overwrite or download")

    # Set Entrez email
    download.set_ncbi_email(args.email)
    logger.info(f"Setting Entrez email address: {args.email}")

    # Parse Entrez API key, if provided
    api_path = args.api_keypath.expanduser()
    if not api_path.is_file():
        logger.warning(f"API path {api_path} not a valid file. Not using API key.")
        api_key = None
    else:
        api_key = download.parse_api_key(api_path)
        logger.info(f"API key recovered from {api_path}")

    # Get list of taxon IDs to download
    taxon_ids = download.split_taxa(args.taxon)
    logger.info(f"Taxon IDs received: {taxon_ids}")

    # Get assembly UIDs for each taxon
    asm_dict = tools.make_asm_dict(taxon_ids, args.retries)
    for tid, uids in asm_dict.items():
        logger.info(
            f"Taxon ID summary\n\tQuery: {tid}\n\tasm count: {len(uids)}\n\tUIDs: {uids}"
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
        logger.info(f"Downloading contigs for Taxon ID {tid}")
        for uid in uids:
            # Obtain eSummary
            logger.info(f"Get eSummary information for UID {uid}")
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
            logger.info(f"eSummary information:\n\t{outstr}")
            if args.dryrun:
                logger.warning(
                    f"(dry-run) skipping download of {esummary['AssemblyAccession']}"
                )
                continue

            # Obtain URLs, trying the RefSeq filestem first, then GenBank if
            # there's a failure
            dlfiledata = DLFileData(
                filestem, "ftp://ftp.ncbi.nlm.nih.gov/genomes/all", "genomic.fna.gz"
            )
            logger.info(f"Retrieving URLs for {filestem}")
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
                logger.warning("RefSeq failed. Trying GenBank alternative assembly")
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
            logger.info(f"Downloaded from URL: {dlstatus.url}")
            logger.info(f"Wrote assembly to: {dlstatus.outfname}")
            logger.info(f"Wrote MD5 hashes to: {dlstatus.outfhash}")

            # Check hash for the download
            hashstatus = download.check_hash(dlstatus.outfname, dlstatus.outfhash)
            logger.info(f"Local MD5 hash: {hashstatus.localhash}")
            logger.info(f"NCBI MD5 hash: {hashstatus.filehash}")
            if hashstatus.passed:
                logger.info("MD5 hash check passed")
            else:
                logger.warning("MD5 hash check failed. Please check and retry.")

            # Extract downloaded files
            ename = dlstatus.outfname.with_suffix("")  # should strip only last suffix
            if ename.exists() and args.noclobber:
                logger.warning(f"Output file {ename} exists, not extracting")
            else:
                logger.info(f"Extracting archive {dlstatus.outfname} to {ename}")
                download.extract_contigs(dlstatus.outfname, ename)

            # Modify sequence ID header if Kraken option active
            if args.kraken:
                logger.warning("Modifying downloaded sequence for Kraken compatibility")
                seqdata = list(SeqIO.parse(ename, "fasta"))
                logger.info(f"Modifying {ename}")
                for seq in seqdata:
                    seq.id = "|".join(
                        [seq.id, "kraken:taxid", esummary["SpeciesTaxid"]]
                    )
                SeqIO.write(seqdata, ename, "fasta")

            # Create MD5 hash for the downloaded contigs
            logger.info(f"Creating local MD5 hash for {ename}")
            hashfname = ename.stem + ".md5"
            datahash = download.create_hash(ename)
            logger.info("Writing hash to %s" % hashfname)
            with open(hashfname, "w") as hfh:
                hfh.write("\t".join([datahash, ename]) + "\n")
            # Make label/class text
            labeltxt, classtxt = download.create_labels(uid_class, filestem, datahash)
            classes.append(classtxt)
            labels.append(labeltxt)
            logger.info(
                f"Label and class file entries\n\tLabel: {labeltxt}\n\tClass: {classtxt}"
            )

    # Write class and label files
    classfname = args.outdir / args.classfname
    logger.info(f"Writing classes file to {classfname}")
    if classfname.exists() and args.noclobber:
        logger.warning(f"Class file {classfname} exists, not overwriting")
    else:
        with open(classfname, "w") as ofh:
            ofh.write("\n".join(classes) + "\n")

    labelfname = args.outdir / args.labelfname
    logger.info(f"Writing labels file to {labelfname}")
    if labelfname.exists() and args.noclobber:
        logger.warning(f"Labels file {labelfname} exists, not overwriting")
    else:
        with open(labelfname, "w") as ofh:
            ofh.write("\n".join(labels) + "\n")

    # Report skipped genome list
    if skippedlist:
        logger.warning(f"{len(skippedlist)} genome downloads were skipped")
        for skipped in skippedlist:
            outstr = "\n\t".join(
                [
                    f"taxon id: {skipped.taxon_id}",
                    f"accession: {skipped.accession}",
                    f"URL: {skipped.url}",
                    f"source: {skipped.dltype}",
                ]
            )
            logger.warning(f"{skipped.organism} {skipped.strain}:\n\t{outstr}")
