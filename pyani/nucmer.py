# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2019
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
"""Code for handling NUCmer output files."""

import os


class DeltaData:

    """Class to hold MUMmer/nucmer output "delta" data.

    This is required because the ordering within files differs depending on MUMmer
    build, for the same version (as evidenced by differences between OSX and Linux
    builds), and a means of testing for equality of outputs is necessary.

    The output file structure and format is described at
    http://mummer.sourceforge.net/manual/#nucmeroutput

    Each file is represented as:

    - header: first line of the .delta file, naming the two input comparison files; stored
              as a tuple (path1, path2), returned as the combined string; the individual
              files are stored as self._query and self._subject
    - program: name of the MUMmer program that produced the output
    - query: path to the query sequence file
    - subject: path to the subject sequence file
    """

    def __init__(self, name, handle=None):
        """Initialise DeltaData object.

        :param name:
        :param handle:
        """
        self.name = name
        self._metadata = None
        self._comparisons = []
        if handle is not None:
            self.from_delta(handle)

    def from_delta(self, handle):
        """Populate the object from the passed .delta or .filter filehandle."""
        parser = DeltaIterator(handle)
        for element in parser:
            if isinstance(element, DeltaMetadata):
                self._metadata = element
            if isinstance(element, DeltaComparison):
                self._comparisons.append(element)

    @property
    def comparisons(self):
        """Comparisons in the .delta file."""
        return self._comparisons

    @property
    def metadata(self):
        """Metadata from the .delta file."""
        return self._metadata

    @property
    def reference(self):
        """Reference file for the MUMmer comparison."""
        return self._metadata.reference

    @property
    def program(self):
        """The MUMmer program used for the comparison."""
        return self._metadata.program

    @property
    def query(self):
        """Query file for the MUMmer comparison."""
        return self._metadata.query

    def __eq__(self, other):
        # We do not enforce equality of metadata, as the full path to both query and reference is
        # written in the .delta file, and we care only about the alignment data, and the program
        # that was used.
        if not isinstance(other, DeltaData):
            return False
        return (self.program == other.program) and (
            self._comparisons == other._comparisons
        )

    def __len__(self):
        return len(self._comparisons)

    def __str__(self):
        """Return the object in .delta format output."""
        outstr = os.linesep.join(
            [str(self._metadata)] + [str(_) for _ in self._comparisons]
        )
        return outstr


class DeltaHeader:

    """Represents a single sequence comparison header from a MUMmer .delta file."""

    def __init__(self, reference, query, reflen, querylen):
        """Initialise DeltaHeader object.

        :param reference:
        :param query:
        :param reflen:
        :param querylen:
        """
        self.reference = reference
        self.query = query
        self.referencelen = int(reflen)
        self.querylen = int(querylen)

    def __eq__(self, other):
        if not isinstance(other, DeltaHeader):
            return False
        return (self.reference, self.query, self.referencelen, self.querylen) == (
            other.reference,
            other.query,
            other.referencelen,
            other.querylen,
        )

    def __str__(self):
        return ">{} {} {} {}".format(
            self.reference, self.query, self.referencelen, self.querylen
        )


class DeltaAlignment:

    """Represents a single alignment region and scores for a pairwise comparison."""

    def __init__(self, refstart, refend, qrystart, qryend, errs, simerrs, stops):
        """Initialise DeltaAlignment object.

        :param refstart:
        :param refend:
        :param qrystart:
        :param qryend:
        :param errs:
        :param simerrs:
        :param stops:
        """
        self.refstart = int(refstart)
        self.refend = int(refend)
        self.querystart = int(qrystart)
        self.queryend = int(qryend)
        self.errs = int(errs)
        self.simerrs = int(simerrs)
        self.stops = int(stops)
        self.indels = []

    def __lt__(self, other):
        return (self.refstart, self.refend, self.querystart, self.queryend) < (
            other.refstart,
            other.refend,
            other.querystart,
            other.queryend,
        )

    def __eq__(self, other):
        return (self.refstart, self.refend, self.querystart, self.queryend) == (
            other.refstart,
            other.refend,
            other.querystart,
            other.queryend,
        )

    def __str__(self):
        outstr = [
            "{} {} {} {} {} {} {}".format(
                self.refstart,
                self.refend,
                self.querystart,
                self.queryend,
                self.errs,
                self.simerrs,
                self.stops,
            )
        ] + [str(_) for _ in self.indels]
        return os.linesep.join(outstr)


class DeltaMetadata:

    """Represents the metadata header for a MUMmer .delta file."""

    def __init__(self):
        """Initialise DeltaMetadata object."""
        self.reference = None
        self.query = None
        self.program = None

    def __eq__(self, other):
        if not isinstance(other, DeltaMetadata):
            return False
        return (self.reference, self.query, self.program) == (
            other.reference,
            other.query,
            other.program,
        )

    def __str__(self):
        return "{} {}{}{}".format(self.reference, self.query, os.linesep, self.program)


class DeltaComparison:

    """Represents a comparison between two sequences in a .delta file."""

    def __init__(self, header, alignments):
        """Initialise DeltaComparison object.

        :param header:
        :param alignments:
        """
        self.header = header
        self.alignments = alignments

    def add_alignment(self, aln):
        """Add passed alignment to this object.

        :param aln:  DeltaAlignment object
        """
        self.alignments.append(aln)

    def __eq__(self, other):
        if not isinstance(other, DeltaComparison):
            return False
        return (self.header == other.header) and (
            sorted(self.alignments) == sorted(other.alignments)
        )

    def __len__(self):
        return len(self.alignments)

    def __str__(self):
        outstr = os.linesep.join([str(self.header)] + [str(_) for _ in self.alignments])
        return outstr


class DeltaIterator:

    """Iterator for MUMmer .delta files.

    Returns a stream of DeltaMetadata, DeltaComparison and DeltaAlignment
    objects when iterated over a filehandle

    The .delta file structure and format is described at
    http://mummer.sourceforge.net/manual/#nucmeroutput
    """

    def __init__(self, handle):
        """Instantiate DeltaIterator object with the passed filehandle.

        :param handle:
        """
        self._handle = handle
        self._metadata = None  # metadata for a .delta file
        self._header = None  # header information for a pairwise comparison
        self._comparison = None  # current comparison region

    def __iter__(self):
        """Iterate over elements of the .delta file as DeltaHeader and DeltaAlignment objects."""
        return iter(self.__next__, None)

    def __next__(self):
        """Parse the next element from the .delta file."""
        # Parse .delta file metadata
        if self._metadata is None:
            self._metadata = DeltaMetadata()
            self._metadata.reference, self._metadata.query = (
                self._handle.readline().strip().split()
            )
            self._metadata.program = self._handle.readline().strip()
            return self._metadata

        # Parse remaining lines into a DeltaHeader for each comparison, and corresponding
        # DeltaAlignments
        line = self._handle.readline()
        while line:
            # If we're at the start of an alignment, create a new DeltaAlignment
            if line.startswith(">"):
                if self._comparison is not None:
                    return self._comparison
                self._header = DeltaHeader(*(line[1:].split()))
                self._comparison = DeltaComparison(self._header, [])
            # Populate the current pairwise alignment with each individual alignment
            else:
                alndata = line.rstrip().split()
                if len(alndata) > 1:  # alignment header
                    alignment = DeltaAlignment(*alndata)
                elif alndata[0] == "0":
                    alignment.indels.append(alndata[0])
                    self._comparison.add_alignment(alignment)
                else:
                    alignment.indels.append(alndata[0])
            # Get the next line and return the final comparison if we're at the end of file
            line = self._handle.readline()
            if not line:
                return self._comparison
