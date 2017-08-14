# -*- coding: utf-8 -*-
"""Module providing functions to generate clusters/species hypotheses.

(c) The James Hutton Institute 2016-2017
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

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

Copyright (c) 2016-2017 The James Hutton Institute

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

import networkx as nx
import pandas as pd


# Build an undirected graph from an ANIResults object
def build_graph_from_results(results, cov_min, id_min):
    """Return undirected graph representing the passed ANIResults object.

    The passed ANIResults object is converted to an undirected graph where
    nodes on the graph represent genomes, and edges represent pairwise
    comparisons having the minimum coverage and identity indicated.

    results     - ANIResults object with result of ANI analysis
    cov_min     - minimum coverage for an edge
    id_min      - minimum identity for an edge
    """
    # Generate reordered data as dataframe
    node_names = results.coverage.columns
    rows = []
    for idx, node_from in enumerate(node_names[:-1]):
        for node_to in node_names[idx+1:]:
            datadict = {'from': node_from,
                        'to': node_to,
                        'coverage': min(results.coverage[node_from][node_to],
                                        results.coverage[node_to][node_from]),
                        'identity': results.identity[node_from][node_to]}
            rows.append(datadict)
    node_data = pd.DataFrame(rows,
                             columns=['from', 'to', 'coverage', 'identity'])

    # Convert reordered data to undirected graph and return
    G = nx.from_pandas_dataframe(node_data, 'from', 'to',
                                 ['coverage', 'identity'])
    return G
