#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""subcmd_plot.py

Provides the plot subcommand for pyani

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

from pyani import pyani_config, pyani_db, pyani_graphics


def subcmd_plot(args, logger):
    """Produce graphical output for an analysis.

    This is graphical output for representing the ANI analysis results, and
    takes the form of a heatmap, or heatmap with dendrogram.
    """
    # Announce what's going on to the user
    logger.info("Generating graphical output for analyses")
    logger.info("Writing output to: %s", args.outdir)
    logger.info("Rendering method: %s", args.method)

    # Distribution dictionary of graphics methods
    gmethod = {
        "mpl": pyani_graphics.heatmap_mpl,
        "seaborn": pyani_graphics.heatmap_seaborn,
    }

    # Work on each run:
    run_ids = [int(run) for run in args.run_id.split(",")]
    logger.info("Generating graphics for runs: %s", run_ids)
    for run_id in run_ids:

        # Get results matrices for the run
        logger.info("Acquiring results for run %d", run_id)
        results = pyani_db.ANIResults(args.dbpath, run_id)

        # Parse output formats
        outfmts = args.formats.split(",")
        logger.info("Requested output formats: %s", outfmts)

        # Generate filestems
        for matname in [
            "identity",
            "coverage",
            "aln_lengths",
            "sim_errors",
            "hadamard",
        ]:
            df = getattr(results, matname)  # results matrix
            cmap = pyani_config.get_colormap(df, matname)
            for fmt in outfmts:
                outfname = os.path.join(
                    args.outdir, "matrix_{0}_{1}.{2}".format(matname, run_id, fmt)
                )
                logger.info("Writing graphics to %s", outfname)
                params = pyani_graphics.Params(cmap, results.labels, results.classes)
                # Draw figure
                gmethod[args.method](
                    df,
                    outfname,
                    title="matrix_{0}_{1}".format(matname, run_id),
                    params=params,
                )
