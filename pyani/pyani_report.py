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
"""Module providing functions for presenting analysis/db output."""

import sys

import seaborn as sns

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import pandas as pd  # type: ignore


def cluster_data(dfr):
    """Reorder contents of matrices to match the clustering performed
    by pyani plot.

    :param dfr:  the dataframe to be reordered

    """
    fig = sns.clustermap(dfr)
    cols = [f"Genome_id:{_+1}" for _ in fig.dendrogram_col.reordered_ind]
    rows = [f"Genome_id:{_+1}" for _ in fig.dendrogram_row.reordered_ind]
    return dfr.loc[rows, cols]


def colour_rows(
    series: pd.Series, even_colour: str = "#DDECF5", odd_colour: str = "#6CB6E4"
) -> List[str]:
    """Return alternating colours for rows in a dataframe.

    :param series:  pd.Series
    :param even_colour:  str, hex colour for even rows
    :param odd_colour:  str, hex colour for odd rows
    """
    is_odd = [idx % 2 for idx, row in enumerate(series.index)]
    return [
        "background-color: %s" % odd_colour
        if v
        else "background-color: %s" % even_colour
        for v in is_odd
    ]


def table_padding() -> Dict[str, Any]:
    """Return HTML for table cell padding."""
    return dict(selector="td", props=[("padding", "15px")])


def hover_highlight(hover_colour: str = "#FFFF99") -> Dict[str, Any]:
    """Return HTML style to colour dataframe row when hovering.

    :param hover_colour:  str, hex colour for hover highlight
    """
    return dict(selector="tr:hover", props=[("background-color", "%s" % hover_colour)])


def header_font() -> Dict[str, Any]:
    """Return header HTML font style."""
    return dict(
        selector="th",
        props=[
            ("text-align", "center"),
            ("font-family", "Helvetica"),
            ("font-size", "small"),
        ],
    )


def colour_identity(
    series: pd.Series, threshold: float = 0.95, colour: str = "#FF2222"
) -> List[str]:
    """Highlight percentage identities over a threshold.

    :param series:
    :param threshold:  float, threshold for cell highlighting
    :param colour:  str, hex colour for highlighted cells
    """
    if series.name == "percentage identity":
        mask = series >= threshold
        return ["color: %s" % colour if v else "" for v in mask]
    return ["" for v in series]


def colour_coverage(
    series: pd.Series, threshold: float = 0.95, colour: str = "#FF2222"
) -> List[str]:
    """Highlight percent coverage over a threshold.

    :param series:
    :param threshold:  float, threshold for cell highlighting
    :param colour:  str, hex colour for highlighted cells
    """
    if "coverage" in str(series.name):
        mask = series >= threshold
        return ["color: %s" % colour if v else "" for v in mask]
    return ["" for v in series]


def colour_numeric(val: float, threshold: float = 0.95, colour: str = "#FF2222") -> str:
    """Highlight numeric values over a threshold.

    :param val:
    :param threshold:  float, threshold for cell highlighting
    :param colour:  str, hex colour for highlighted cell
    """
    if val < threshold:
        colour = "black"
    return "color: %s" % colour


# Write a dataframe in pyani-styled HTML
def write_styled_html(
    path: Path, dfm: pd.DataFrame, index: Optional[str] = None, colour_num: bool = False
) -> None:
    """Add CSS styling to a dataframe and write as HTML.

    :param path:       path to write output file
    :param dfm:         dataframe to be written out
    :param index:      column to be set as index (if necessary)
    """
    # Reset the index to a specified column
    if index is not None and index in dfm.columns:
        dfm.set_index(index, inplace=True)

    # Colour rows in alternating shades of blue
    styled = dfm.style.apply(colour_rows)

    # Colour percentage identity threshold/coverage values > 95% in red
    styled = styled.apply(colour_identity).apply(colour_coverage)

    # Colour numbers over a given threshold
    if colour_num:
        styled = styled.applymap(colour_numeric)

    # Apply styles
    styled = styled.set_table_styles(
        [hover_highlight(), header_font(), table_padding()]
    )

    # Set font to Helvetica
    styled = styled.set_properties(**{"font-family": "Helvetica", "font-size": "small"})

    # Write styled HTML to path
    html = styled.render()
    with open(path, "w") as ofh:
        ofh.write(html)


# Write a dataframe to STDOUT
def write_to_stdout(
    stem: str, dfm: pd.DataFrame, show_index: bool = False, line_width: float = None
) -> None:
    """Write dataframe in tab-separated form to STDOUT.

    :param stem:  str
    :param dfm:  pd.Dataframe
    :param show_index:  Boolean, include index in output table
    :param line_width:
    """
    sys.stdout.write(f"TABLE: {stem}\n")
    sys.stdout.write(dfm.to_string(index=show_index, line_width=line_width) + "\n\n")


# Write a table returned from the pyani database in the requested format
def write_dbtable(
    dfm: pd.DataFrame,
    path: Path,
    formats: Sequence[str] = ("tab",),
    show_index: bool = True,
    colour_num: bool = False,
) -> None:
    """Write database result table to output file in named format.

    :param dfm:  pd.Dataframe
    :param path:  Path to output file
    :param formats:  tuple of str, output file formats
    :param show_index:  output row and column labels
    :param colour_num:  use colours for values in HTML output

    colours are used for identity/coverage tables
    """
    formatdict = {
        "tab": (dfm.to_csv, {"sep": "\t", "index": show_index}, ".tab"),
        "excel": (dfm.to_excel, {"index": show_index}, ".xlsx"),
        "html": (
            write_styled_html,
            {"dfm": dfm, "index": show_index, "colour_num": colour_num},
            ".html",
        ),
        "stdout": (write_to_stdout, {"dfm": dfm, "show_index": show_index}, ""),
    }
    for fmt in formats:
        func, args, ext = formatdict[fmt]
        ofname = path.with_suffix(ext)
        func(ofname, **args)
