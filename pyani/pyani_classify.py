#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module providing functions to generate clusters/species hypotheses.

(c) The James Hutton Institute 2016-2018
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD2 5DA,
Scotland,
UK

The MIT License

Copyright (c) 2016-2018 The James Hutton Institute

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

from collections import namedtuple

import networkx as nx
import pandas as pd

from pyani.pyani_tools import label_results_matrix

# Holds summary information about a graph's cliques
Cliquesinfo = namedtuple(
    "Cliquesinfo", "n_nodes n_subgraphs n_cliques n_cliquenodes confused"
)


# Build an undirected graph from an ANIResults object
def build_graph_from_results(results, label_dict, cov_min=0, id_min=0):
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
def analyse_cliques(graph):
    """Return Cliquesinfo namedtuple describing clique data for a graph."""
    cliques = list(nx.find_cliques(graph))
    tot_clique_members = sum([len(c) for c in cliques])
    subgraphs = [graph.subgraph(_).copy() for _ in nx.connected_components(graph)]
    return Cliquesinfo(
        len(graph),
        len(subgraphs),
        len(cliques),
        tot_clique_members,
        tot_clique_members - len(graph),
    )


# Generate a list of graphs from lowest to highest pairwise identity threshold
def trimmed_graph_sequence(ingraph, attribute="identity"):
    """Return graphs trimmed from lowest to highest attribute value

    A generator which, starting from the initial graph, yields in sequence a
    series of graphs from which the edge(s) with the lowest threshold value
    attribute were removed. The generator returns a tuple of:

    (threshold, graph, analyse_cliques(graph))

    ingraph          - the initial graph to work from
    attribute  - string describing the attribute to work on

    This will be slow with moderate-large graphs
    """
    graph = ingraph.copy()
    edgelist = sorted(graph.edges(data=attribute), key=lambda x: x[-1])
    while len(edgelist) > 1:
        threshold = edgelist[0][-1]
        yield (threshold, graph, analyse_cliques(graph))
        while edgelist[0][-1] <= threshold:
            edge = edgelist.pop(0)
            graph.remove_edge(edge[0], edge[1])
    # For last edge/graph
    threshold = edgelist[0][-1]
    yield (threshold, graph, analyse_cliques(graph))


# Generate a list of graphs with no clique confusion
def unconfused_graphs(ingraph, attribute="identity"):
    """Return graphs having no clique-confused nodes

    A generator which, starting from the initial graph, yields in sequence
    a series of graphs from lowest to highest threshold edge where no node
    in the graph participates in more than one clique.

    ingraph         - the initial graph to start from
    attribute  - the attribute to use for thresholds

    This will be slow with moderate-large graphs
    """
    graph = ingraph.copy()
    for subgraph in trimmed_graph_sequence(graph, attribute):
        if subgraph[-1].confused:
            continue
        yield subgraph
