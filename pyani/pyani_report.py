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


def highlight_odd_rows(s):
    """Highlight odd rows in a dataframe."""
    print(s)


# Write a dataframe in pyani-styled HTML
def write_styled_html(path, index, df):
    """Add CSS styling to a dataframe."""
    numeric_col_mask = df.dtypes.apply(lambda d:
                                       issubclass(np.dtype(d).type,
                                                  np.number))
    # Header styles
    header_styles = dict(selector="th",
                         props=[('text-align', 'center'),
                                ('font-family', 'Helvetica')])
    # Pre-render HTML with styling for dataframe
    # Set global background and font colour for cells, also cell font
    df = df.style.set_properties(**{'background-color': 'black',
                                    'color': 'lawngreen',
                                    'border-color': 'white',
                                    'font-family': 'Helvetica',
                                    'font-size': 'small'})
    df = df.set_properties(subset=df.columns[numeric_col_mask], 
                                 **{'width':'10em', 'text-align':'right'})
    df = df.set_properties(subset=df.columns[~numeric_col_mask],
                                 **{'width':'10em', 'text-align':'left'})
    df = df.format(lambda x: '{:,.0f}'.format(x) if
                   x > 1e3 else '{:,.2f}'.format(x),
                   subset=pd.IndexSlice[:,
                                        df.columns[numeric_col_mask]])
    df = df.set_table_styles([header_styles])
    html = df.render()
    # Write styled HTML to path
    with open(path, 'w') as ofh:
        ofh.write(html)


# Write a table returned from the pyani database in the requested format
def write_dbtable(data, headers, path=None, formats=('tab',)):
    """Write database result table to output file in named format."""
    df = pd.DataFrame(data)
    df.columns = headers
    formatdict = {'tab': (df.to_csv, {'sep': '\t'}, '.tab'),
                  'excel': (df.to_excel, {}, '.xlsx'),
                  'html': (write_styled_html, {'df': df}, '.html')}
    for format in formats:
        func, args, ext = formatdict[format]
        ofname = path + ext
        func(ofname, index=False, **args)
