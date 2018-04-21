#!/usr/bin/env python

"""test_anib.py

Test anib.py module.

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
DD2 5DA,
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
import unittest

import pandas as pd

from nose.tools import (assert_equal,)
from pandas.util.testing import (assert_frame_equal,)

from pyani import (anib, pyani_files)


class TestBLASTCmdline(unittest.TestCase):

    """Class defining tests of BLAST command-line generation."""

    def setUp(self):
        """Set parameters for tests."""
        self.indir = os.path.join('tests', 'test_input', 'anib')
        self.outdir = os.path.join('tests', 'test_output', 'anib')
        self.seqdir = os.path.join('tests', 'test_input', 'sequences')
        self.infiles = [os.path.join(self.seqdir, fname) for fname in
                        os.listdir(self.seqdir)]
        self.fraglen = 1000
        self.fmtdboutdir = os.path.join(self.outdir, 'formatdb')
        self.fmtdbcmd = ' '.join(["formatdb -p F -i",
                                  "tests/test_output/anib/formatdb/NC_002696.fna",
                                  "-t NC_002696"])
        self.makeblastdbdir = os.path.join(self.outdir, 'makeblastdb')
        self.makeblastdbcmd = ' '.join(["makeblastdb -dbtype nucl -in",
                                        "tests/test_input/sequences/NC_002696.fna",
                                        "-title NC_002696 -out",
                                        os.path.join('tests', 'test_output',
                                                     'anib', 'makeblastdb',
                                                     'NC_002696.fna')])
        self.blastdbfnames = [os.path.join(self.seqdir, fname) for fname in
                              ('NC_002696.fna', 'NC_010338.fna')]
        self.blastdbtgt = [(' '.join(['makeblastdb -dbtype nucl -in',
                                      'tests/test_input/sequences/NC_002696.fna',
                                      '-title NC_002696 -out',
                                      'tests/test_output/anib/NC_002696.fna']),
                            'tests/test_output/anib/NC_002696.fna'),
                           (' '.join(['makeblastdb -dbtype nucl -in',
                                      'tests/test_input/sequences/NC_010338.fna',
                                      '-title NC_010338 -out',
                                      'tests/test_output/anib/NC_010338.fna']),
                            'tests/test_output/anib/NC_010338.fna')]
        self.blastdbtgtlegacy = [(' '.join(['formatdb -p F -i',
                                            'tests/test_output/anib/NC_002696.fna',
                                            '-t NC_002696']),
                                  'tests/test_output/anib/NC_002696.fna'),
                                 (' '.join(['formatdb -p F -i',
                                            'tests/test_output/anib/NC_010338.fna',
                                            '-t NC_010338']),
                                  'tests/test_output/anib/NC_010338.fna')]
        self.blastncmd = ' '.join(["blastn -out",
                                   os.path.join("tests", "test_output", "anib",
                                                "NC_002696_vs_NC_010338.blast_tab"),
                                   "-query tests/test_input/sequences/NC_002696.fna",
                                   "-db tests/test_input/sequences/NC_010338.fna",
                                   "-xdrop_gap_final 150 -dust no -evalue 1e-15",
                                   "-max_target_seqs 1 -outfmt '6 qseqid sseqid",
                                   "length mismatch pident nident qlen slen qstart",
                                   "qend sstart send positive ppos gaps' -task blastn"])
        self.blastallcmd = ' '.join(['blastall -p blastn -o',
                                     os.path.join('tests', 'test_output', 'anib',
                                                  'NC_002696_vs_NC_010338.blast_tab'),
                                     '-i tests/test_input/sequences/NC_002696.fna',
                                     '-d tests/test_input/sequences/NC_010338.fna',
                                     '-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8'])
        self.blastntgt = [' '.join(["blastn -out",
                                    os.path.join("tests", "test_output", "anib",
                                                 "NC_002696_vs_NC_010338.blast_tab"),
                                    "-query tests/test_input/sequences/NC_002696.fna",
                                    "-db tests/test_input/sequences/NC_010338.fna",
                                    "-xdrop_gap_final 150 -dust no -evalue 1e-15",
                                    "-max_target_seqs 1 -outfmt '6 qseqid sseqid",
                                    "length mismatch pident nident qlen slen qstart",
                                    "qend sstart send positive ppos gaps' -task blastn"]),
                          ' '.join(["blastn -out",
                                    os.path.join("tests", "test_output", "anib",
                                                 "NC_010338_vs_NC_002696.blast_tab"),
                                    "-query tests/test_input/sequences/NC_010338.fna",
                                    "-db tests/test_input/sequences/NC_002696.fna",
                                    "-xdrop_gap_final 150 -dust no -evalue 1e-15",
                                    "-max_target_seqs 1 -outfmt '6 qseqid sseqid length",
                                    "mismatch pident nident qlen slen qstart qend",
                                    "sstart send positive ppos gaps' -task blastn"])]
        self.blastalltgt = [' '.join(['blastall -p blastn -o',
                                      os.path.join('tests', 'test_output', 'anib',
                                                   'NC_002696_vs_NC_010338.blast_tab'),
                                      '-i tests/test_input/sequences/NC_002696.fna',
                                      '-d tests/test_input/sequences/NC_010338.fna',
                                      '-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8']),
                            ' '.join(['blastall -p blastn -o',
                                      os.path.join('tests', 'test_output', 'anib',
                                                   'NC_010338_vs_NC_002696.blast_tab'),
                                      '-i tests/test_input/sequences/NC_010338.fna',
                                      '-d tests/test_input/sequences/NC_002696.fna',
                                      '-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8'])]
        self.blastnjobdict = sorted(
            [('tests/test_output/anib/NC_002696.fna',
              'makeblastdb -dbtype nucl ' +
              '-in tests/test_input/sequences/NC_002696.fna ' +
              '-title NC_002696 -out tests/test_output/anib/NC_002696.fna'),
             ('tests/test_output/anib/NC_010338.fna',
              'makeblastdb -dbtype nucl ' +
              '-in tests/test_input/sequences/NC_010338.fna ' +
              '-title NC_010338 -out tests/test_output/anib/NC_010338.fna'),
             ('tests/test_output/anib/NC_011916.fna',
              'makeblastdb -dbtype nucl ' +
              '-in tests/test_input/sequences/NC_011916.fna ' +
              '-title NC_011916 -out tests/test_output/anib/NC_011916.fna'),
             ('tests/test_output/anib/NC_014100.fna',
              'makeblastdb -dbtype nucl ' +
              '-in tests/test_input/sequences/NC_014100.fna ' +
              '-title NC_014100 -out tests/test_output/anib/NC_014100.fna')])
        self.blastalljobdict = sorted(
            [('tests/test_output/anib/NC_002696.fna',
              'formatdb -p F -i tests/test_output/anib/NC_002696.fna ' +
              '-t NC_002696'),
             ('tests/test_output/anib/NC_010338.fna',
              'formatdb -p F -i tests/test_output/anib/NC_010338.fna ' +
              '-t NC_010338'),
             ('tests/test_output/anib/NC_011916.fna',
              'formatdb -p F -i tests/test_output/anib/NC_011916.fna ' +
              '-t NC_011916'),
             ('tests/test_output/anib/NC_014100.fna',
              'formatdb -p F -i tests/test_output/anib/NC_014100.fna ' +
              '-t NC_014100')])
        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(self.fmtdboutdir, exist_ok=True)
        os.makedirs(self.makeblastdbdir, exist_ok=True)

    def test_formatdb_generation(self):
        """generate formatdb command-line."""
        cmd = anib.construct_formatdb_cmd(os.path.join(self.seqdir,
                                                       "NC_002696.fna"),
                                          self.fmtdboutdir)
        assert_equal(cmd[0], self.fmtdbcmd)  # correct command
        assert(os.path.isfile(cmd[1]))       # creates new file

    def test_makeblastdb_generation(self):
        """generate makeblastdb command-line."""
        cmd = anib.construct_makeblastdb_cmd(os.path.join(self.seqdir,
                                                          "NC_002696.fna"),
                                             self.makeblastdbdir)
        assert_equal(cmd[0], self.makeblastdbcmd)  # correct command

    def test_blastdb_commands(self):
        """generate both BLAST+ and legacy BLAST db commands."""
        # BLAST+
        cmds = anib.generate_blastdb_commands(self.blastdbfnames, self.outdir,
                                              mode="ANIb")
        assert_equal(cmds, self.blastdbtgt)
        # legacy
        cmds = anib.generate_blastdb_commands(self.blastdbfnames, self.outdir,
                                              mode="ANIblastall")
        assert_equal(cmds, self.blastdbtgtlegacy)

    def test_blastn_generation(self):
        """generate BLASTN+ command-line."""
        cmd = anib.construct_blastn_cmdline(self.blastdbfnames[0],
                                            self.blastdbfnames[1],
                                            self.outdir)
        assert_equal(cmd, self.blastncmd)

    def test_blastall_generation(self):
        """generate legacy BLASTN command-line."""
        cmd = anib.construct_blastall_cmdline(self.blastdbfnames[0],
                                              self.blastdbfnames[1],
                                              self.outdir)
        assert_equal(cmd, self.blastallcmd)

    def test_blastn_commands(self):
        """generate both BLASTN+ and legacy BLASTN commands."""
        # BLAST+
        cmds = anib.generate_blastn_commands(self.blastdbfnames, self.outdir,
                                             mode="ANIb")
        assert_equal(cmds, self.blastntgt)
        cmds = anib.generate_blastn_commands(self.blastdbfnames, self.outdir,
                                             mode="ANIblastall")
        assert_equal(cmds, self.blastalltgt)

    def test_blastall_dbjobdict(self):
        """generate dictionary of legacy BLASTN database jobs."""
        blastcmds = anib.make_blastcmd_builder("ANIblastall", self.outdir)
        jobdict = anib.build_db_jobs(self.infiles, blastcmds)
        assert_equal(sorted([(k, v.script) for (k, v) in jobdict.items()]),
                     self.blastalljobdict)

    def test_blastn_dbjobdict(self):
        """generate dictionary of BLASTN+ database jobs."""
        blastcmds = anib.make_blastcmd_builder("ANIb", self.outdir)
        jobdict = anib.build_db_jobs(self.infiles, blastcmds)
        assert_equal(sorted([(k, v.script) for (k, v) in jobdict.items()]),
                     self.blastnjobdict)

    def test_blastn_graph(self):
        """create jobgraph for BLASTN jobs."""
        fragresult = anib.fragment_fasta_files(self.infiles, self.outdir,
                                               self.fraglen)
        blastcmds = anib.make_blastcmd_builder("ANIb", self.outdir)
        jobgraph = anib.make_job_graph(self.infiles, fragresult[0],
                                       blastcmds)
        # We check that the main script job is a blastn job, and that there
        # is a single dependency, which is a makeblastdb job
        for job in jobgraph:
            assert(job.script.startswith('blastn'))
            assert_equal(1, len(job.dependencies))
            dep = job.dependencies[0]
            assert(dep.script.startswith('makeblastdb'))

    def test_blastall_graph(self):
        """create jobgraph for legacy BLASTN jobs."""
        fragresult = anib.fragment_fasta_files(self.infiles, self.outdir,
                                               self.fraglen)
        blastcmds = anib.make_blastcmd_builder("ANIblastall", self.outdir)
        jobgraph = anib.make_job_graph(self.infiles, fragresult[0],
                                       blastcmds)
        # We check that the main script job is a blastn job, and that there
        # is a single dependency, which is a makeblastdb job
        for job in jobgraph:
            assert(job.script.startswith('blastall -p blastn'))
            assert_equal(1, len(job.dependencies))
            dep = job.dependencies[0]
            assert(dep.script.startswith('formatdb'))


class TestFragments(unittest.TestCase):

    """Class defining tests of ANIb FASTA fragmentation"""

    def setUp(self):
        """Initialise parameters for tests."""
        self.outdir = os.path.join('tests', 'test_output', 'anib')
        self.seqdir = os.path.join('tests', 'test_input', 'sequences')
        self.infnames = [os.path.join(self.seqdir, fname) for fname in
                         ('NC_002696.fna', 'NC_010338.fna',
                          'NC_011916.fna', 'NC_014100.fna')]
        self.outfnames = [os.path.join(self.outdir, fname) for fname in
                          ('NC_002696-fragments.fna', 'NC_010338-fragments.fna',
                           'NC_011916-fragments.fna', 'NC_014100-fragments.fna')]
        self.fraglen = 1000
        os.makedirs(self.outdir, exist_ok=True)

    def test_fragment_files(self):
        """fragment files for ANIb/ANIblastall."""
        result = anib.fragment_fasta_files(self.infnames, self.outdir,
                                           self.fraglen)
        # Are files created?
        for outfname in self.outfnames:
            assert os.path.isfile(outfname)

        # Test fragment lengths
        for accession, fragdict in result[-1].items():
            for fragname, fraglen in fragdict.items():
                assert fraglen <= self.fraglen


class TestParsing(unittest.TestCase):

    """Class defining tests of BLAST output parsing."""

    def setUp(self):
        self.indir = os.path.join('tests', 'test_input', 'anib')
        self.seqdir = os.path.join('tests', 'test_input', 'sequences')
        self.fragdir = os.path.join('tests', 'test_input', 'anib', 'fragfiles')
        self.anibdir = os.path.join('tests', 'test_input', 'anib', 'blastn')
        self.aniblastalldir = os.path.join(
            'tests', 'test_input', 'anib', 'blastall')
        self.fname_legacy = os.path.join(self.indir,
                                         "NC_002696_vs_NC_010338.blast_tab")
        self.fname = os.path.join(self.indir,
                                  "NC_002696_vs_NC_011916.blast_tab")
        self.fragfname = os.path.join(self.indir,
                                      "NC_002696-fragments.fna")
        self.fraglens = 1000
        self.infnames = [os.path.join(self.seqdir, fname) for fname in
                         ('NC_002696.fna', 'NC_010338.fna',
                          'NC_011916.fna', 'NC_014100.fna')]
        self.fragfiles = [os.path.join(self.fragdir, fname) for fname in
                          ('NC_002696-fragments.fna', 'NC_010338-fragments.fna',
                           'NC_011916-fragments.fna', 'NC_014100-fragments.fna')]
        self.anibtgt = pd.DataFrame([[1.000000, 0.796974, 0.999977, 0.837285],
                                     [0.795958, 1.000000, 0.795917, 0.798250],
                                     [0.999922, 0.795392, 1.000000, 0.837633],
                                     [0.836780, 0.798704, 0.836823, 1.000000]],
                                    columns=['NC_002696', 'NC_010338',
                                             'NC_011916', 'NC_014100'],
                                    index=['NC_002696', 'NC_010338', 'NC_011916',
                                           'NC_014100'])
        self.aniblastalltgt = pd.DataFrame([[1.000000, 0.785790, 0.999977, 0.830641],
                                            [0.781319, 1.000000,
                                                0.781281, 0.782723],
                                            [0.999937, 0.782968,
                                                1.000000, 0.830431],
                                            [0.828919, 0.784533, 0.828853, 1.000000]],
                                           columns=['NC_002696', 'NC_010338',
                                                    'NC_011916', 'NC_014100'],
                                           index=['NC_002696', 'NC_010338', 'NC_011916',
                                                  'NC_014100'])

    def test_parse_blasttab(self):
        """parses ANIblastall .blast_tab output."""
        fragdata = anib.get_fraglength_dict([self.fragfname])
        # ANIb output
        result = anib.parse_blast_tab(self.fname, fragdata,
                                      mode="ANIb")
        assert_equal(result, (4016551, 93, 99.997693577050029))
        # ANIblastall output
        result = anib.parse_blast_tab(self.fname_legacy, fragdata,
                                      mode="ANIblastall")
        assert_equal(result, (1966922, 406104, 78.578978313253018))

    def test_blastdir_processing(self):
        """parses directory of .blast_tab output."""
        orglengths = pyani_files.get_sequence_lengths(self.infnames)
        fraglengths = anib.get_fraglength_dict(self.fragfiles)
        # ANIb
        result = anib.process_blast(self.anibdir, orglengths,
                                    fraglengths, mode="ANIb")
        assert_frame_equal(result.percentage_identity.sort_index(1).sort_index(),
                           self.anibtgt.sort_index(1).sort_index())
        # ANIblastall
        result = anib.process_blast(self.aniblastalldir, orglengths,
                                    fraglengths, mode="ANIblastall")
        assert_frame_equal(result.percentage_identity.sort_index(1).sort_index(),
                           self.aniblastalltgt.sort_index(1).sort_index())
