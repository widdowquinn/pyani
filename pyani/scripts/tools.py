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
"""Module providing functions in support of the pyani command-line scripts."""

import re
import shutil

from pyani import download


# EXCEPTIONS
# ==========
# General errors that occur in a script
class PyaniScriptException(Exception):

    """General exception for pyani.py script."""

    def __init__(self, msg="Error in pyani.py script"):
        """Insantiate exception."""
        Exception.__init__(self, msg)


# UTILITY FUNCTIONS
# =================

# Create a directory (handling force/noclobber options)
def make_outdir(outdir, force, noclobber, logger):
    """Create output directory (allows for force and noclobber).

    :param outdir:  Path, path to output directory
    :param force:
    :param noclobber:
    :param logger:

    The intended outcomes are:
    outdir doesn't exist: create outdir
    outdir exists: raise exception
    outdir exists, --force only: remove the directory tree
    outdir exists, --force --noclobber:
        continue with existing directory tree but do not overwrite files

    So long as the outdir is created with this function, we need only check
    for args.noclobber elsewhere to see how to proceed when a file exists.
    """
    if outdir.is_dir():
        logger.warning("Output directory %s exists", outdir)
        if force and not noclobber:  # user forces directory reuse, and overwrite
            # Delete old directory and start again
            logger.warning(
                "Overwrite forced. Removing %s and everything " + "below it (--force)",
                outdir,
            )
            shutil.rmtree(outdir)
        elif force and noclobber:  # user forces directory reuse, not file overwrite
            logger.warning(
                "Keeping existing directory, skipping existing "
                + "files (--force --noclobber)."
            )
        else:  # user is not forcing directory reuse
            raise PyaniScriptException(
                "Will not modify existing directory " + "%s" % outdir
            )
    outdir.mkdir(exist_ok=True)


# Make a dictionary of assembly download info
def make_asm_dict(taxon_ids, retries):
    """Return a dict of assembly UIDs, keyed by each passed taxon ID.

    :param taxon_ids:
    :param retries:
    """
    asm_dict = dict()

    for tid in taxon_ids:
        asm_uids = download.get_asm_uids(tid, retries)
        asm_dict[tid] = asm_uids.asm_ids

    return asm_dict


# Download the RefSeq genome and MD5 hash from NCBI
def download_genome_and_hash(
    args, logger, dlfiledata, dltype="RefSeq", disable_tqdm=False
):
    """Download genome and accompanying MD5 hash from NCBI.

    :param args:  Namespace for command-line arguments
    :param logger:  logging object
    :param dlfiledata:  namedtuple of info for file to download
    :param dltype:  reference database to use: RefSeq or GenBank
    :param disable_tqdm:  disable progress bar

    This function tries the (assumed to be passed) RefSeq FTP URL first and,
    if that fails, then attempts to download the corresponding GenBank data.

    We attempt to gracefully skip genomes with download errors.
    """
    if dltype == "GenBank":
        filestem = re.sub("^GCF_", "GCA_", dlfiledata.filestem)
    else:
        filestem = dlfiledata.filestem
    dlstatus = download.retrieve_genome_and_hash(
        filestem,
        dlfiledata.suffix,
        dlfiledata.ftpstem,
        args.outdir,
        args.timeout,
        disable_tqdm,
    )
    # Pylint is confused by the content of dlstatus (a namedlist)
    if dlstatus.error is not None:  # pylint: disable=no-member
        logger.warning(
            "%s download failed: skipping!\n%s",
            dltype,
            dlstatus.error,  # pylint: disable=no-member
        )
        dlstatus.skipped = True

    return dlstatus
