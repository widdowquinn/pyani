#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""tools.py

Provides tools to support tests in the pyani package

The pylint error C0103 is disabled as we follow the camelCase example
of unittest for our tests, rather than maintaining snakecase.

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
import json
import os
import unittest

import pandas as pd

from pyani import blast, nucmer


class PyaniFileEqualityTests(unittest.TestCase):
    """Tests for equality of filetypes used in pyani.

    Each test defines a comparison for a specific filetype, with contents that are
    expected to be equal in some way.
    """

    def assertJsonEqual(self, json1, json2):  # pylint: disable=C0103
        """Assert that two passed JSON files are equal.
        :param json1:  path to reference JSON file
        :param json2:  path to comparator JSON file
        As we can't always be sure that JSON elements are in the same order in otherwise
        equal files, we compare ordered components.
        """
        with open(json1, "r") as fh1:
            with open(json2, "r") as fh2:
                self.assertEqual(
                    ordered(json.load(fh1)),
                    ordered(json.load(fh2)),
                    msg="JSON files {} and {} do not contain equal contents".format(
                        json1, json2
                    ),
                )

    def assertNucmerEqual(self, fname1, fname2):  # pylint: disable=C0103
        """Assert that two passed nucmer output files are equal.
        :param fname1:  path to reference .delta/.filter file
        :param fname2:  path to comparator .delta/.filter file
        This is a standard file comparison, skipping the first line.
        """
        with open(fname1, "r") as fh1:
            with open(fname2, "r") as fh2:
                fdata1 = nucmer.DeltaData("fh1", fh1)
                fdata2 = nucmer.DeltaData("fh2", fh2)
                self.assertEqual(
                    fdata1,
                    fdata2,
                    msg="Nucmer files {} and {} are not equivalent".format(
                        fname1, fname2
                    ),
                )

    def assertBlasttabEqual(self, fname1, fname2):  # pylint: disable=C0103
        """Assert that two passed BLAST+ .tab output files contain the same data.
        This is not a simple comparison, as we can't rely on the same ordering,
        so we parse the files and compare objects.
        """
        with open(fname1, "r") as fh1:
            with open(fname2, "r") as fh2:
                data1 = blast.parse_blasttab(fh1)
                data2 = blast.parse_blasttab(fh2)
                for line1, line2 in zip(data1, data2):
                    self.assertEqual(line1, line2)

    def assertBinaryEqual(self, fname1, fname2):  # pylint: disable=C0103
        """Assert that two passed binary files contain the same data"""
        with open(fname1, "rb") as fh1:
            with open(fname2, "rb") as fh2:
                data1 = bytearray(fh1.read())
                data2 = bytearray(fh2.read())
        self.assertEqual(data1, data2)

    def assertFilesEqual(self, fname1, fname2):  # pylint: disable=C0103
        """Assert that the two passed files have the same contents."""
        with open(fname1, "r") as fh1:
            with open(fname2, "r") as fh2:
                self.assertEqual(
                    fh1.read(),
                    fh2.read(),
                    msg="Files {} and {} do not have the same contents".format(
                        fname1, fname2
                    ),
                )

    def assertTabEqual(self, fname1, fname2):  # pylint: disable=C0103
        """Assert that two tabular files are essentially equal

        To do this, we need to load each .tab file as a dataframe, order rows
        and columns, then compare.
        """
        df1 = pd.read_csv(fname1, sep="\t", index_col=0)
        df1 = df1.reindex(columns=sorted(df1.columns)).reindex(index=sorted(df1.index))
        df2 = pd.read_csv(fname2, sep="\t", index_col=0)
        df2 = df2.reindex(columns=sorted(df2.columns)).reindex(index=sorted(df2.index))
        if not df1.equals(df2):
            raise AssertionError(
                f"Dataframes {fname1} and {fname2} are not equal.\n{df1}\n{df2}"
            )


class PyaniTestCase(PyaniFileEqualityTests, unittest.TestCase):
    """Specific pyani unit tests."""

    def assertDirsEqual(self, dir1, dir2, filt=None):  # pylint: disable=C0103
        """Assert that two passed directories have the same contents.

        :param dir1:  Reference directory for comparison
        :param dir2:  Comparator direcotry for comparison

        Files to be compared can be restricted using the filter argument. For
        instance:

        assertDirsEqual(d1, d2, filter=".tab") will only compare files with
        the tab extension.

        Directories are compared recursively.
        """
        # List directories and skip hidden files
        dir1files = ordered([_ for _ in os.listdir(dir1) if not _.startswith(".")])
        dir2files = ordered([_ for _ in os.listdir(dir2) if not _.startswith(".")])
        self.assertEqual(
            dir1files,
            dir2files,
            msg="{} and {} do not have same file listings".format(dir1, dir2),
        )
        #  Compare contents of directories; descend through directories, but
        # filter file extensions if needed
        if filt is not None:
            dir1files = [
                _
                for _ in dir1files
                if (os.path.isdir(_) is False) and (os.path.splitext(_)[-1] == filt)
            ]
        for fpath in dir1files:
            if os.path.isdir(os.path.join(dir1, fpath)):  # Compare dictionaries
                self.assertDirsEqual(
                    os.path.join(dir1, fpath), os.path.join(dir2, fpath)
                )
            else:  # Compare files
                ext = os.path.splitext(fpath)[-1]
                fname1 = os.path.join(dir1, fpath)
                fname2 = os.path.join(dir2, fpath)
                if ext.lower() in (".gz", ".pdf", ".png"):  # skip these files
                    continue
                elif ext.lower() == ".tab":  # Compare tab-separated tabular files
                    self.assertTabEqual(fname1, fname2)
                elif ext.lower() == ".json":  # Compare JSON files
                    self.assertJsonEqual(fname1, fname2)
                elif ext.lower() == ".blasttab":  # Compare BLAST+ .tab output
                    self.assertBlasttabEqual(fname1, fname2)
                elif ext.lower() in (
                    ".delta",
                    ".filter",
                ):  # Compare nucmer/delta-filter output
                    self.assertNucmerEqual(fname1, fname2)
                # elif ext.lower() in (".png",):  # Compare binary files
                #     self.assertBinaryEqual(fname1, fname2)
                else:  # Compare standard files
                    self.assertFilesEqual(fname1, fname2)


def ordered(obj):
    """Return ordered version of the passed object
    Dictionaries are not ordered in all Python versions, and the
    implementation of sort_keys() in the the JSON library seems
    erratic in terms of effect
    """
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        try:
            return sorted(ordered(x) for x in obj)
        except TypeError:  # list contains non-comparable types
            return obj
    else:
        return obj


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
