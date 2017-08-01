# -*- coding: utf-8 -*-
"""Module providing useful functions for presenting analysis/db output

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

import os

import numpy as np
import pandas as pd

def colour_rows(series, even="#DDECF5", odd="#6CB6E4"):
    """Colour rows in a dataframe."""
    is_odd = [row % 2 for row in list(series.index)]
    return ['background-color: %s' % odd if v
            else 'background-color: %s' % even for v in is_odd]


def table_padding():
    """Return HTML for table cell padding."""
    return dict(selector="td",
                props=[('padding', '15px')])


def hover_highlight(hover_colour="#FFFF99"):
    """Return HTML style to colour dataframe row when hovering."""
    return dict(selector="tr:hover",
                props=[('background-color', '%s' % hover_colour)])


def header_font():
    """Return header HTML font style."""
    return dict(selector="th",
                props=[('text-align', 'center'),
                       ('font-family', 'Helvetica'),
                       ('font-size', 'small')])


def colour_identity(series, threshold=0.95, colour="#FF2222"):
    """Highlight percentage identities over a threshold."""
    if series.name == 'percentage identity':
        mask = series >= threshold
        return ['color: %s' % colour if v else '' for v in mask]
    else:
        return ['' for v in series]


def colour_coverage(series, threshold=0.95, colour="#FF2222"):
    """Highlight percent coverage over a threshold."""
    if 'coverage' in series.name:
        mask = series >= threshold
        return ['color: %s' % colour if v else '' for v in mask]
    else:
        return ['' for v in series]


# Write a dataframe in pyani-styled HTML
def write_styled_html(path, df, index=None):
    """Add CSS styling to a dataframe and write as HTML.

    path       path to write output file
    df         dataframe to be written out
    index      column to be set as index (if necessary)
    """
    # Reset the index to a specified column
    if index is not None and index in df.columns:
        df.set_index(index, inplace=True)

    # Colour rows in alternating shades of blue
    styled = df.style.apply(colour_rows)

    # Colour percentage identity threshold/coverage values > 95% in red
    styled = styled.apply(colour_identity).apply(colour_coverage)

    # Apply styles
    styled = styled.set_table_styles([hover_highlight(),
                                      header_font(),
                                      table_padding()])

    # Set font to Helvetica
    styled = styled.set_properties(**{'font-family': 'Helvetica',
                                      'font-size': 'small'})

    # Write styled HTML to path
    html = styled.render()
    with open(path, 'w') as ofh:
        ofh.write(html)


# Write a table returned from the pyani database in the requested format
def write_dbtable(data, headers, path=None, formats=('tab',), index=False):
    """Write database result table to output file in named format."""
    df = pd.DataFrame(data)
    df.columns = headers
    formatdict = {'tab': (df.to_csv, {'sep': '\t', 'index': False}, '.tab'),
                  'excel': (df.to_excel, {'index': False}, '.xlsx'),
                  'html': (write_styled_html, {'df': df, 'index': index},
                           '.html')}
    for format in formats:
        func, args, ext = formatdict[format]
        ofname = path + ext
        func(ofname, **args)
