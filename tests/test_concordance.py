#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""test_concordance.py

Test for concordance of pyani package output with JSpecies

These tests are intended to be run from the repository root using:

nosetests -v

print() statements will be caught by nosetests unless there is an
error. They can also be recovered with the -s option.

(c) The James Hutton Institute 2017
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

Copyright (c) 2017 The James Hutton Institute

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
import shutil
import subprocess
import sys
import unittest

import pandas as pd

from nose.tools import assert_equal, assert_less, nottest

from pyani import run_multiprocessing as run_mp
from pyani import (anib, anim, tetra, pyani_files, pyani_config)


def parse_jspecies(infile):
    """Parse JSpecies output into Pandas dataframes.

    The function expects a single file containing (legacy) ANIb,
    ANIm, and TETRA output.

    - infile        path to JSpecies output file

    This is an ugly function!
    """
    dfs = dict()
    methods = ('ANIm', 'ANIb', 'Tetra')
    with open(infile, 'r') as ifh:
        header, in_table = False, False
        for line in [l.strip() for l in ifh.readlines() + ['\n']]:
            if line in methods and not in_table:
                method, header = line, True
            elif header:
                columns = line.split('\t')
                data = pd.DataFrame(index=columns, columns=columns)
                in_table, header = True, False
            elif in_table:
                if not len(line):
                    dfs[method] = data.sort_index(axis=0).sort_index(axis=1)
                    in_table = False
                else:
                    ldata = line.split('\t')
                    row = ldata[0]
                    for idx, val in enumerate(ldata[1:]):
                        if val != '---':
                            data[columns[idx]][row] = float(val)
                        elif method.startswith('ANI'):
                            data[columns[idx]][row] = 100.
                        else:
                            data[columns[idx]][row] = 1.
            else:
                pass
    return dfs


class TestConcordance(unittest.TestCase):

    """Class defining tests of pyani concordance with JSpecies."""

    def setUp(self):
        """Set values and parameters for tests."""
        self.indir = os.path.join('tests', 'test_input', 'concordance')
        self.outdir = os.path.join('tests', 'test_output', 'concordance')
        self.tgtdir = os.path.join('tests', 'test_targets', 'concordance')
        self.deltadir = os.path.join(self.outdir, 'nucmer_output')
        self.infiles = pyani_files.get_fasta_files(self.indir)
        self.orglengths = pyani_files.get_sequence_lengths(self.infiles)
        self.target = parse_jspecies(os.path.join(self.tgtdir, 'jspecies_output.tab'))
        self.tolerance = {'ANIm': 0.1, 'ANIb_lo': 5, 'ANIb_hi': 0.1,
                          'ANIblastall': 0.1, 'TETRA': 0.1}
        self.fragsize = 1020
        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(self.deltadir, exist_ok=True)

    def test_anim_concordance(self):
        """ANIm results concordant with JSpecies."""
        # Perform ANIm on the input directory contents
        # We have to separate nucmer/delta-filter command generation
        # because Travis-CI doesn't play nicely with changes we made
        # for local SGE/OGE integration.
        # This might be avoidable with a scheduler flag passed to
        # jobgroup generation in the anim.py module. That's a TODO.
        ncmds, fcmds = anim.generate_nucmer_commands(self.infiles,
                                                     self.outdir)
        run_mp.multiprocessing_run(ncmds)

        # delta-filter commands need to be treated with care for
        # Travis-CI. Our cluster won't take redirection or semicolon
        # separation in individual commands, but the wrapper we wrote
        # for this (delta_filter_wrapper.py) can't be called under
        # Travis-CI. So we must deconstruct the commands below
        dfcmds = [' > '.join([' '.join(fcmd.split()[1:-1]),
                              fcmd.split()[-1]]) for fcmd in fcmds]
        run_mp.multiprocessing_run(dfcmds)

        results = anim.process_deltadir(self.deltadir, self.orglengths)
        result_pid = results.percentage_identity
        result_pid.to_csv(os.path.join(self.outdir, 'pyani_anim.tab'),
                          sep="\t")

        # Compare JSpecies output to results
        result_pid = result_pid.sort_index(axis=0).sort_index(axis=1) * 100.
        diffmat = result_pid.as_matrix() - self.target['ANIm'].as_matrix()
        anim_diff = pd.DataFrame(diffmat, index=result_pid.index,
                                 columns=result_pid.columns)
        anim_diff.to_csv(os.path.join(self.outdir, 'pyani_anim_diff.tab'),
                         sep="\t")
        assert_less(anim_diff.abs().values.max(), self.tolerance['ANIm'])

    def test_anib_concordance(self):
        """ANIb results concordant with JSpecies.

        We expect ANIb results to be quite different, as the BLASTN
        algorithm changed substantially between BLAST and BLAST+
        """
        # Perform ANIb on the input directory contents
        outdir = os.path.join(self.outdir, 'blastn')
        os.makedirs(outdir, exist_ok=True)
        fragfiles, fraglengths = anib.fragment_fasta_files(self.infiles,
                                                           outdir,
                                                           self.fragsize)
        jobgraph = anib.make_job_graph(self.infiles, fragfiles,
                                       anib.make_blastcmd_builder('ANIb',
                                                                  outdir))
        assert_equal(0, run_mp.run_dependency_graph(jobgraph))
        results = anib.process_blast(outdir, self.orglengths, fraglengths,
                                     mode="ANIb")
        result_pid = results.percentage_identity
        result_pid.to_csv(os.path.join(self.outdir, 'pyani_anib.tab'),
                          sep="\t")

        # Compare JSpecies output to results. We do this in two blocks,
        # masked according to whether the expected result is greater than
        # 90% identity, or less than that threshold.
        # The complete difference matrix is written to output, though
        result_pid = result_pid.sort_index(axis=0).sort_index(axis=1) * 100.
        lo_result = result_pid.mask(result_pid >= 90).fillna(0)
        hi_result = result_pid.mask(result_pid < 90).fillna(0)
        lo_target = self.target['ANIb'].mask(self.target['ANIb'] >= 90).fillna(0)
        hi_target = self.target['ANIb'].mask(self.target['ANIb'] < 90).fillna(0)
        lo_diffmat = lo_result.as_matrix() - lo_target.as_matrix()
        hi_diffmat = hi_result.as_matrix() - hi_target.as_matrix()
        diffmat = result_pid.as_matrix() - self.target['ANIb'].as_matrix()
        lo_diff = pd.DataFrame(lo_diffmat, index=result_pid.index,
                               columns=result_pid.columns)
        hi_diff = pd.DataFrame(hi_diffmat, index=result_pid.index,
                               columns=result_pid.columns)
        anib_diff = pd.DataFrame(diffmat, index=result_pid.index,
                                 columns=result_pid.columns)
        anib_diff.to_csv(os.path.join(self.outdir, 'pyani_anib_diff.tab'),
                         sep="\t")
        assert_less(lo_diff.abs().values.max(), self.tolerance['ANIb_lo'])
        assert_less(hi_diff.abs().values.max(), self.tolerance['ANIb_hi'])

    def test_aniblastall_concordance(self):
        """ANIblastall results concordant with JSpecies."""
        # Perform ANIblastall on the input directory contents
        outdir = os.path.join(self.outdir, 'blastall')
        os.makedirs(outdir, exist_ok=True)
        fragfiles, fraglengths = anib.fragment_fasta_files(self.infiles,
                                                           outdir,
                                                           self.fragsize)
        jobgraph = anib.make_job_graph(self.infiles, fragfiles,
                                       anib.make_blastcmd_builder('ANIblastall',
                                                                  outdir))
        assert_equal(0, run_mp.run_dependency_graph(jobgraph))
        results = anib.process_blast(outdir, self.orglengths, fraglengths,
                                     mode="ANIblastall")
        result_pid = results.percentage_identity
        result_pid.to_csv(os.path.join(self.outdir, 'pyani_aniblastall.tab'),
                          sep="\t")

        # Compare JSpecies output to results
        result_pid = result_pid.sort_index(axis=0).sort_index(axis=1) * 100.
        diffmat = result_pid.as_matrix() - self.target['ANIb'].as_matrix()
        aniblastall_diff = pd.DataFrame(diffmat, index=result_pid.index,
                                        columns=result_pid.columns)
        aniblastall_diff.to_csv(os.path.join(self.outdir, 'pyani_aniblastall_diff.tab'),
                                sep="\t")
        assert_less(aniblastall_diff.abs().values.max(),
                    self.tolerance['ANIblastall'])

    def test_tetra_concordance(self):
        """TETRA results concordant with JSpecies."""
        # Perform TETRA analysis
        zscores = dict()
        for filename in self.infiles:
            org = os.path.splitext(os.path.split(filename)[-1])[0]
            zscores[org] = tetra.calculate_tetra_zscore(filename)
        results = tetra.calculate_correlations(zscores)
        results.to_csv(os.path.join(self.outdir, 'pyani_tetra.tab'), sep="\t")

        # Compare JSpecies output
        diffmat = results.as_matrix() - self.target['Tetra'].as_matrix()
        tetra_diff = pd.DataFrame(diffmat, index=results.index,
                                  columns=results.columns)
        tetra_diff.to_csv(os.path.join(self.outdir, 'pyani_tetra_diff.tab'),
                          sep="\t")
        assert_less(tetra_diff.abs().values.max(),
                    self.tolerance['TETRA'])
