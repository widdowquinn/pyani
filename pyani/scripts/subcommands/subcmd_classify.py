#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
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
# Copyright (c) 2017-2019 The James Hutton Institute
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
"""Provides the classify subcommand for pyani."""

import logging

from argparse import Namespace
from typing import Generator, NamedTuple

import networkx as nx  # pylint: disable=E0401
import numpy as np

from tqdm import tqdm

from pyani import pyani_classify, pyani_orm
from pyani.pyani_tools import termcolor


class SubgraphData(NamedTuple):

    """Subgraph clique/classification output."""

    interval: float  # the trimming threshold for this subgraph
    graph: nx.Graph  # the trimmed subgraph
    cliqueinfo: pyani_classify.Cliquesinfo


def subcmd_classify(args: Namespace) -> int:
    """Generate classifications for an analysis.

    :param args:  Namespace, command-line arguments
    """
    logger = logging.getLogger(__name__)

    # Tell the user what's going on
    logger.info(
        termcolor("Generating classification for ANI run: %s", "red"), args.run_id
    )
    logger.info("\tWriting output to: %s", args.outdir)
    logger.info(termcolor("\tCoverage threshold: %s", "cyan"), args.cov_min)
    logger.info(
        termcolor("\tInitial minimum identity threshold: %s", "cyan"), args.id_min
    )

    # Get results data for the specified run
    logger.info("Acquiring results for run: %s", args.run_id)
    logger.debug("Connecting to database: %s", args.dbpath)
    session = pyani_orm.get_session(args.dbpath)
    logger.debug("Retrieving results matrices")
    results = (
        session.query(pyani_orm.Run).filter(pyani_orm.Run.run_id == args.run_id).first()
    )
    result_label_dict = pyani_orm.get_matrix_labels_for_run(session, args.run_id)

    # Generate initial graph on basis of results
    logger.info("Constructing graph from results.")
    initgraph = pyani_classify.build_graph_from_results(
        results, result_label_dict, args.cov_min, args.id_min
    )
    logger.debug(
        "Returned graph has %d nodes:\n\t%s",
        len(initgraph),
        "\n\t".join(n for n in initgraph),
    )
    logger.debug(
        "Initial graph clique information:\n\t%s",
        pyani_classify.analyse_cliques(initgraph),
    )

    # Obtain all subgraph splits, thresholding by identity
    subgraphs = trimmed_graph_sequence(initgraph, args)
    special_intervals = [_ for _ in subgraphs if _.cliqueinfo.all_k_complete]
    outstr = "\n\t".join([f"{_.interval}\t{_.cliqueinfo}" for _ in special_intervals])
    logger.info(
        "%s intervals with special property:\n\t%s", len(special_intervals), outstr
    )
    if args.show_all:
        outstr = "\n\t".join([f"{_.interval}\t{_.cliqueinfo}" for _ in subgraphs])
        logger.debug("Subgraphs at all identity thresholds:\n\t%s", outstr)

    return 0


# Generate a list of graphs from lowest to highest pairwise identity threshold
def trimmed_graph_sequence(
    ingraph: nx.Graph, args: Namespace, attribute: str = "identity"
) -> Generator:
    """Return graphs trimmed from lowest to highest attribute value.

    :param ingraph:  nx.Graph of genomes as nodes, having edges weighted by
                     the property named in attribute
    :param args:  Namespace, parsed command-line arguments
    :param attribute:  str, name of the property by which the graph edges
                       should be trimmed

    A generator which, starting from the initial graph, yields in sequence a
    series of graphs from which the edge(s) with the lowest threshold value
    attribute were removed. The generator returns a tuple of:

    (threshold, graph, analyse_cliques(graph))

    This will be slow with moderate-large graphs
    """
    graph = ingraph.copy()

    # Trim the graph now, removing edges at the minimum allowed identity and below
    graph, edgelist = pyani_classify.remove_low_weight_edges(
        graph, args.min_id, attribute
    )

    #  There's no sense resolving to small intervals if there would be more
    # of these than there are edges to remove from the graph. If the
    # resolution leads to more breakpoints than there are edges, just trim
    # edge by edge. Otherwise work over a range of property value breaks.
    if len(edgelist) < 1 / args.resolution:
        breaks = [_[-1] for _ in edgelist]
    else:
        breaks = np.arange(
            args.min_id or edgelist[0][-1], args.max_id or 1, args.resolution
        )

    # Remove edges at each breakpoint and yield the resulting subgraph
    for threshold in tqdm(breaks, disable=args.disable_tqdm):
        yield SubgraphData(threshold, graph, pyani_classify.analyse_cliques(graph))
        while edgelist and edgelist[0][-1] <= threshold:
            edge = edgelist.pop(0)
            graph.remove_edge(edge[0], edge[1])

    # Yield final edge/graph (identity)
    yield SubgraphData(
        1,
        pyani_classify.remove_low_weight_edges(graph, 1, attribute)[0],
        pyani_classify.analyse_cliques(graph),
    )
