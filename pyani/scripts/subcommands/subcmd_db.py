#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""subcmd_db.py

Provides the db subcommand for pyani

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

from pyani import pyani_db


def subcmd_db(args, logger):
    """Perform operations on an existing database."""
    # If the database exists, raise an error rather than overwrite
    if not os.path.isfile(args.dbpath):
        logger.error("Database %s does not exist (exiting)", args.dbpath)
        raise SystemError(1)

    logger.info("Working with pyani database at %s", args.dbpath)
    logger.info("Working with run ID %d", args.run_id)

    # If requested, carry out relabelling of genomes
    if args.relabelfname is not None:
        logger.info("Relabelling genomes from file %s", args.relabelfname)
        pyani_db.relabel_genomes_from_file(
            args.dbpath, args.relabelfname, args.run_id, args.force
        )

    # If requested, change classes of genomes
    if args.reclassfname is not None:
        logger.info("Changing classes of genomes from file %s", args.reclassfname)
        pyani_db.reclass_genomes_from_file(
            args.dbpath, args.reclassfname, args.run_id, args.force
        )
