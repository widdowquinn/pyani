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
import subprocess

from argparse import Namespace
from typing import Dict, List, NamedTuple, Optional, Tuple

from Bio import SeqIO

from pyani import download
from pyani.pyani_tools import termcolor
from pyani.scripts import make_outdir


class Skipped(NamedTuple):

    """Convenience struct for holding information about skipped genomes."""

    taxon_id: str
    accession: str
    organism: str
    strain: str
    url: str
    dltype: str


def configure_entrez(args: Namespace) -> Optional[str]:
    """Configure Entrez email, return API key.

    :param args:  Namespace, command-line arguments

    Returns None if no API key found
    """
    logger = logging.getLogger(__name__)

    download.set_ncbi_email(args.email)
    logger.info("Setting Entrez email address: %s", args.email)
    return parse_api_key(args)


def dl_info_to_str(esummary, uid_class) -> str:
    """Return descriptive string for passed download data.

    :param esummary:
    :param uid_class:
    """
    return "\n\t".join(
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


def download_data(
    args: Namespace,
    api_key: Optional[str],
    asm_dict: Dict[str, List],
) -> Tuple[List, List, List]:
    """Download the accessions indicated in the passed dictionary.

    :param args:  Namespace of command-line arguments
    :param api_key:  str, API key for NCBI downloads
    :param asm_dict:  dictionary of assembly UIDs to download, keyed by taxID

    Returns lists of information about downloaded genome classes and labels, and a
    list of skipped downloads (as Skipped objects).
    """
    logger = logging.getLogger(__name__)

    classes, labels, skippedlist = [], [], []

    for tid, uids in asm_dict.items():
        logger.info(termcolor("Downloading contigs for Taxon ID %s", "blue"), uids)
        for uid in uids:
            # Obtain eSummary for each assembly UID
            logger.info(
                termcolor("Retrieving eSummary information for UID %s", "cyan"), uid
            )
            try:
                esummary, filestem = download.get_ncbi_esummary(
                    uid, args.retries, api_key
                )
            except download.NCBIDownloadException:
                logger.warning(
                    termcolor("Skipping download of record for UID %s", "red"),
                    uid,
                    exc_info=True,
                )
                skippedlist.append(
                    Skipped(
                        tid,
                        uid,
                        "",
                        "",
                        None,
                        "RefSeq",
                    )
                )  # pylint: disable=no-member
                continue

            uid_class = download.get_ncbi_classification(esummary)
            logger.debug(
                "eSummary information (%s):\n\t%s",
                filestem,
                dl_info_to_str(esummary, uid_class),
            )

            if args.dryrun:
                logger.warning(
                    "(dry-run) skipping download of %s", esummary["AssemblyAccession"]
                )
                continue

            # Download genome for UID, and extract compressed files
            dlstatus, skipped_genomes = download_genome(
                args, filestem, tid, uid, uid_class
            )
            skippedlist.extend(skipped_genomes)
            if not dlstatus.skipped:
                extract_genomes(args, dlstatus, esummary)
                labeltxt, classtxt = hash_genomes(args, dlstatus, filestem, uid_class)
                classes.append(classtxt)
                labels.append(labeltxt)
                logger.info(
                    "Label and class file entries\n\tLabel: %s\n\tClass: %s",
                    labeltxt,
                    classtxt,
                )

    return classes, labels, skippedlist


def extract_genomes(args: Namespace, dlstatus: download.DLStatus, esummary) -> None:
    """Extract genome files in passed dlstatus.

    :param args:  Namespace of command-line arguments
    :param dlstatus:
    :param esummary:
    """
    logger = logging.getLogger(__name__)

    # Extract downloaded files
    ename = dlstatus.outfname.with_suffix("")  # should strip only last suffix
    if ename.exists() and args.noclobber:
        logger.warning("Output file %s exists, not extracting", ename)
    else:
        logger.debug("Extracting archive %s to %s", dlstatus.outfname, ename)
        try:
            download.extract_contigs(dlstatus.outfname, ename)
        except subprocess.CalledProcessError:
            logger.warning("Could not extract %s; continuing", dlstatus.outfname)
            pass

    # Modify sequence ID header if Kraken option active
    if args.kraken:
        logger.warning("Modifying downloaded sequence for Kraken compatibility")
        seqdata = list(SeqIO.parse(ename, "fasta"))
        logger.debug("Modifying %s", ename)
        for seq in seqdata:
            seq.id = "|".join([seq.id, "kraken:taxid", esummary["SpeciesTaxid"]])
        SeqIO.write(seqdata, ename, "fasta")


def hash_genomes(
    args: Namespace, dlstatus: download.DLStatus, filestem: str, uid_class
) -> Tuple[str, str]:
    """Hash genome files in passed dlstatus.

    :param args:  Namespace of command-line arguments
    :param dlstatus:
    :param filestem:  str, filestem for output
    :param uid_class:
    """
    logger = logging.getLogger(__name__)

    # Create MD5 hash for the downloaded contigs
    ename = dlstatus.outfname.with_suffix("")  # should strip only last suffix
    logger.debug("Creating local MD5 hash for %s", ename)
    hashfname = ename.with_suffix(".md5")
    datahash = download.create_hash(ename)
    logger.debug("Writing hash to %s", hashfname)
    with open(hashfname, "w") as hfh:
        hfh.write("\t".join([datahash, str(ename)]) + "\n")
    # Make label/class text
    labeltxt, classtxt = download.create_labels(uid_class, filestem, datahash)
    return labeltxt, classtxt


def download_genome(args: Namespace, filestem: str, tid: str, uid: str, uid_class):
    """Download single genome data to output directory.

    :param args:  Namespace, command-line arguments
    :param filestem:  str, output filestem
    :param tid:  str, taxonID
    :param uid:  str, assembly UID
    :param uid_class:
    """
    logger = logging.getLogger(__name__)

    skippedlist = []
    refseq_status, genbank_status = True, True  # set False if skipped

    dlfiledata = download.DLFileData(
        filestem, "ftp://ftp.ncbi.nlm.nih.gov/genomes/all", "genomic.fna.gz"
    )
    logger.info("Retrieving URLs for %s", filestem)
    # Try RefSeq first
    dlstatus = download.download_genome_and_hash(
        args.outdir,
        args.timeout,
        dlfiledata,
        dltype="RefSeq",
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
                "RefSeq",  # pylint: disable=no-member
            )
        )
        refseq_status = False

    # RefSeq fails, so try GenBank
    if refseq_status is False:
        logger.warning(
            termcolor("RefSeq failed. Trying GenBank alternative assembly", "magenta")
        )
        # Try GenBank assembly
        dlstatus = download.download_genome_and_hash(
            args.outdir,
            args.timeout,
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
            genbank_status = False
            logger.warning(termcolor("GenBank failed.", "magenta"))

    if genbank_status or refseq_status:
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

    return dlstatus, skippedlist


def get_tax_asm_dict(args: Namespace) -> Dict[str, List]:
    """Return dictionary of assembly UIDs to download, keyed by taxID.

    :param args:  Namespace of command-line arguments
    """
    logger = logging.getLogger(__name__)

    taxon_ids = download.split_taxa(args.taxon)
    logger.info(termcolor("Taxon IDs received: %s", "blue"), taxon_ids)
    asm_dict = download.make_asm_dict(taxon_ids, args.retries)
    for tid, uids in asm_dict.items():
        logger.debug(
            "Taxon ID summary\n\tQuery: %s\n\tasm count: %s\n\tUIDs: %s",
            tid,
            len(uids),
            uids,
        )
    return asm_dict


def parse_api_key(args: Namespace) -> Optional[str]:
    """Returns NCBI API key if present, None otherwise.

    :param args:  Namespace of command-line arguments

    Checks for key in args.api_keypath.
    """
    logger = logging.getLogger(__name__)

    api_path = args.api_keypath.expanduser()
    if not api_path.is_file():
        logger.warning("API path %s not a valid file. Not using API key.", api_path)
        api_key = None
    else:
        with api_path.open() as ifh:
            api_key = ifh.readline().strip()
        logger.info("API key recovered from %s", api_path)

    return api_key


def subcmd_download(args: Namespace) -> int:
    """Download assembled genomes in subtree of passed NCBI taxon ID.

    :param args:  Namespace, command-line arguments
    """
    # Create logger
    logger = logging.getLogger(__name__)
    logger.info(termcolor("Downloading genomes from NCBI", "red"))

    # Create output directory, respecting force/noclobber
    if args.dryrun:
        logger.warning(
            termcolor("Dry run only: will not overwrite or download", "cyan")
        )
    else:
        make_outdir(args.outdir, args.force, args.noclobber)

    api_key = configure_entrez(args)  # set up email/get API key
    asm_dict = get_tax_asm_dict(args)  # dictionary of assembly UIDs for download

    # Download contigs and hashes for each assembly UID in the dictionary
    # Collect class and label information for each downloaded genome, plus a list
    # of skipped genome data
    classes, labels, skippedlist = download_data(args, api_key, asm_dict)

    # Write class and label files
    if not args.dryrun:
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
