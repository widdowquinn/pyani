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

from pyani import pyani_classify, pyani_orm


def subcmd_classify(args, logger):
    """Generate classifications for an analysis."""
    # Tell the user what's going on
    logger.info(f"Generating classification for ANI run: {args.run_id}")
    logger.info(f"\tWriting output to: {args.outdir}")
    logger.info(f"\tCoverage threshold: {args.cov_min}")
    logger.info(f"\tInitial minimum identity threshold: {args.id_min}")

    # Get results data for the specified run
    logger.info(f"Acquiring results for run: {args.run_id}")
    logger.info(f"Connecting to database: {args.dbpath}")
    session = pyani_orm.get_session(args.dbpath)
    logger.info("Retrieving results matrices")
    results = (
        session.query(pyani_orm.Run).filter(pyani_orm.Run.run_id == args.run_id).first()
    )
    result_label_dict = pyani_orm.get_matrix_labels_for_run(session, args.run_id)

    # Generate initial graph on basis of results
    logger.info("Constructing graph from results.")
    initgraph = pyani_classify.build_graph_from_results(
        results, result_label_dict, args.cov_min, args.id_min
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
