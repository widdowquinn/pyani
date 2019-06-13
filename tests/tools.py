#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""tools.py

Provides tools to support tests in the pyani package

(c) The James Hutton Institute 2019
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

Copyright (c) 2019 The James Hutton Institute

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

import copy


def modify_namespace(namespace, args):
    """Update arguments in a passed Namespace.

    :param namespace:       argparse.Namespace object
    :param args:            dict of argument: value pairs

    The expected usage pattern is, for a command-line application with many
    or complex arguments, to define a base argparse.Namespace object, then
    change only a few arguments, specific to a test. This function takes
    a base namespace and a dictionary of argument: value pairs, and
    returns the modified namespace.
    """
    new_namespace = copy.deepcopy(namespace)
    for argname, argval in args.items():
        setattr(new_namespace, argname, argval)
    return new_namespace
