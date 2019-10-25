#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
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
# Copyright (c) 2017-2019 The James Hutton Institute
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
"""Provides the createdb subcommand for pyani."""

import os

# from pyani import pyani_db
from pyani import pyani_orm


def subcmd_createdb(args, logger):
    """Create an empty pyani database."""
    # If the database exists, raise an error rather than overwrite
    if os.path.isfile(args.dbpath):
        if not args.force:
            logger.error(f"Database {args.dbpath} already exists (exiting)")
            raise SystemError(1)
        logger.warning(f"Database {args.dbpath} already exists - overwriting")
        os.remove(args.dbpath)

    # If the path to the database doesn't exist, create it
    dbdir = os.path.split(args.dbpath)[0]
    if not os.path.isdir(dbdir):
        logger.info("Creating database directory %s", dbdir)
        os.makedirs(dbdir, exist_ok=True)

    # Create the empty database
    logger.info("Creating pyani database at %s", args.dbpath)
    pyani_orm.create_db(args.dbpath)
