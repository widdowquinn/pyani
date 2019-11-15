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
"""Module providing functions to generate clusters/species hypotheses."""

from typing import List, NamedTuple, Tuple

import networkx as nx
import pandas as pd

from pyani.pyani_tools import label_results_matrix


# Holds summary information about a graph's cliques
class Cliquesinfo(NamedTuple):

    """Summary of clique structure."""

    n_nodes: int
    n_subgraphs: int
    all_k_complete: bool


# Build an undirected graph from an ANIResults object
def build_graph_from_results(
    results, label_dict, cov_min: float = 0, id_min: float = 0
) -> nx.Graph:
    """Return undirected graph representing the passed ANIResults object.

    The passed ANIResults object is converted to an undirected graph where
    nodes on the graph represent genomes, and edges represent pairwise
    comparisons having the minimum coverage and identity indicated.

    :param results:     - Run object from pyani_orm
    :param label_dict:  dictionary of genome labels for result matrices
        the dict is keyed by the index/column values for the results
        matrices
    :param cov_min:     - minimum coverage for an edge
    :param id_min:      - minimum identity for an edge
    """
    # Parse identity and coverage matrices
    mat_identity = label_results_matrix(pd.read_json(results.df_identity), label_dict)
    mat_coverage = label_results_matrix(pd.read_json(results.df_coverage), label_dict)

    node_names = mat_coverage.columns
    rows = []
    # Loop over dataframe entries and add edges to the base graph for classification
    # We loop over only the lower-triangle indices of each matrix, but check the
    # identity and coverage in both directions. Coverage is not expected to be
    # symmetrical, and while identity should be symmetrical for ANIm, it will not
    # likely be so for ANIb and related methods.
    # Where results are non-symmetrical, we take the minimum value when creating the
    # base graph.
    for idx, node_from in enumerate(node_names[:-1]):
        for node_to in node_names[idx + 1 :]:
            datadict = {
                "from": node_from,
                "to": node_to,
                "coverage": min(
                    mat_coverage[node_from][node_to], mat_coverage[node_to][node_from]
                ),
                "identity": min(
                    mat_identity[node_from][node_to], mat_identity[node_to][node_from]
                ),
            }
            if datadict["identity"] > id_min and datadict["coverage"] > cov_min:
                rows.append(datadict)
    node_data = pd.DataFrame(rows, columns=["from", "to", "coverage", "identity"])

    # Convert reordered data to undirected graph and return
    graph = nx.from_pandas_edgelist(node_data, "from", "to", ["coverage", "identity"])
    for node in node_names[:-1]:
        if node not in graph.nodes:
            graph.add_node(node)
    return graph


# Report clique info for a graph
def analyse_cliques(graph: nx.Graph) -> Cliquesinfo:
    """Return Cliquesinfo NamedTuple describing clique data for a graph.

    :param graph:  NetworkX Graph object
    """
    subgraphs = [graph.subgraph(_).copy() for _ in nx.connected_components(graph)]
    return Cliquesinfo(len(graph), len(subgraphs), all_components_k_complete(graph))


def all_components_k_complete(graph: nx.Graph) -> bool:
    """Return True if all components in passed graph are k-complete.

    :param graph:  NetworkX Graph object
    """
    return all(k_complete_component_status(graph))


def k_complete_component_status(graph: nx.Graph) -> List[bool]:
    """Return list of Booleans of whether connected components of the graph are k-complete.

    :param graph:  NetworkX Graph object

    For each component in the passed graph, a list of Booleans is calculated,
    representing whether each node has property P: the degree of the node
    is equal to the number of nodes in that component, minus 1.

    The all() gives a Boolean indicating whether all nodes in that
    component have property P.
    """
    node_degrees = {key: val for key, val in graph.degree}
    return [
        all([(node_degrees[n] == len(c) - 1) for n in c])
        for c in nx.connected_components(graph)
    ]


def remove_low_weight_edges(
    graph: nx.Graph, threshold: float, attribute: str = "identity"
) -> Tuple[nx.Graph, List]:
    """Return graph and edgelist where edges having weight < threshold are removed.

    :param graph:  NetworkX Graph
    :param threshold:  float, minimum edge weight
    :param attribute:  String, attribute to use as weight
    """
    edgelist = sorted(graph.edges(data=attribute), key=lambda x: x[-1])
    threshold = threshold or edgelist[0][-1]
    while edgelist and edgelist[0][-1] <= threshold:
        edge = edgelist.pop(0)
        graph.remove_edge(edge[0], edge[1])
    return graph, edgelist
