#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""subcmd_classify.py

Provides the classify subcommand for pyani

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

from pyani import pyani_classify, pyani_db


def subcmd_classify(args, logger):
    """Generate classifications for an analysis."""
    # Tell the user what's going on
    logger.info("Generating classification for ANI run: %s", args.run_id)
    logger.info("Writing output to: %s", args.outdir)
    logger.info("Coverage threshold: %s", args.cov_min)
    logger.info("Initial minimum identity threshold: %s", args.id_min)

    # Get results data for the specified run
    logger.info("Acquiring results for run: %s", args.run_id)
    results = pyani_db.ANIResults(args.dbpath, args.run_id)

    # Generate initial graph on basis of results
    logger.info("Constructing graph from results.")
    initgraph = pyani_classify.build_graph_from_results(
        results, args.cov_min, args.id_min
    )
    logger.info(
        "Returned graph has %d nodes:\n\t%s",
        len(initgraph),
        "\n\t".join([n for n in initgraph]),
    )
    logger.info(
        "Initial graph clique information:\n\t%s",
        pyani_classify.analyse_cliques(initgraph),
    )

    # Report all subgraphs, thresholding by identity
    logger.info(
        "Summarising cliques at all identity thresholds:\n\t%s",
        "\n\t".join(
            [
                "\t".join([str(gdata[0]), str(gdata[2])])
                for gdata in pyani_classify.trimmed_graph_sequence(initgraph)
            ]
        ),
    )

    # Report 'natural breaks' in the identity subgraphs
    logger.info(
        "Identifying 'natural breaks' with no clique-confusion:\n\t%s",
        "\n\t".join(
            [
                "\t".join([str(gdata[0]), str(gdata[2])])
                for gdata in pyani_classify.unconfused_graphs(initgraph)
            ]
        ),
    )
