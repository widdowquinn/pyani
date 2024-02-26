#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019-2024
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# 161 Cathedral Street,
# Glasgow,
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016-2019 The James Hutton Institute
# Copyright (c) 2019-2024 University of Strathclyde
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
"""Test for concordance of pyani package output with JSpecies.

These tests are intended to be run from the repository root using:

pytest -v

print() statements will be caught by nosetests unless there is an
error. They can also be recovered with the -s option.
"""

import unittest

from pathlib import Path

import pandas as pd
import pytest

from pyani import run_multiprocessing as run_mp
from pyani import anib, aniblastall, anim, tetra, pyani_files


def parse_jspecies(infile):
    """Parse JSpecies output into Pandas dataframes.

    The function expects a single file containing (legacy) ANIb,
    ANIm, and TETRA output.

    :param infile:  path to JSpecies output file

    This is an ugly function!
    """
    dfs = dict()
    methods = ("ANIm", "ANIb", "Tetra")
    with open(infile, "r") as ifh:
        header, in_table = False, False
        for line in [_.strip() for _ in ifh.readlines() + ["\n"]]:
            if line in methods and not in_table:
                mth, header = line, True
            elif header:
                columns = line.split("\t")
                data = pd.DataFrame(index=columns, columns=columns)
                in_table, header = True, False
            elif in_table:
                if not line:
                    dfs[mth] = data.sort_index(axis=0).sort_index(axis=1)
                    in_table = False
                else:
                    ldata = line.split("\t")
                    row = ldata[0]
                    for idx, val in enumerate(ldata[1:]):
                        if val != "---":
                            data.loc[columns[idx], row] = float(val)
                        elif mth.startswith("ANI"):
                            data.loc[columns[idx], row] = 100.0
                        else:
                            data.loc[columns[idx], row] = 1.0
            else:
                pass
    return dfs


@pytest.fixture
def paths_concordance_fna(path_fixtures_base):
    """Paths to FASTA inputs for concordance analysis."""
    return [
        _
        for _ in (path_fixtures_base / "concordance").iterdir()
        if _.is_file() and _.suffix == ".fna"
    ]


@pytest.fixture
def path_concordance_jspecies(path_fixtures_base):
    """Path to JSpecies analysis output."""
    return path_fixtures_base / "concordance/jspecies_output.tab"


@pytest.fixture
def threshold_anib_lo_hi():
    """Threshold for concordance comparison split between high and low identity.

    When comparing ANIb results with ANIblastall results, we need to account for
    the differing performances of BLASTN and BLASTN+ on more distantly-related
    sequences. On closely-related sequences both methods give similar results;
    for more distantly-related sequences, the results can be quite different. This
    threshold is the percentage identity we consider to separate "close" from
    "distant" related sequences.
    """
    return 90


@pytest.fixture
def tolerance_anib_hi():
    """Tolerance for ANIb concordance comparisons.

    This tolerance is for comparisons between "high identity" comparisons, i.e.
    genomes having identity greater than threshold_anib_lo_hi in homologous regions.

    These "within-species" level comparisons need to be more accurate
    """
    return 0.2


@pytest.fixture
def tolerance_anib_lo():
    """Tolerance for ANIb concordance comparisons.

    This tolerance is for comparisons between "low identity" comparisons, i.e.
    genomes having identity less than threshold_anib_lo_hi in homologous regions.

    These "intra-species" level comparisons vary more as a result of the change of
    algorithm from BLASTN to BLASTN+ (megablast).
    """
    return 5


@pytest.fixture
def tolerance_anim():
    """Tolerance for ANIm concordance comparisons."""
    return 0.1


@pytest.fixture
def tolerance_tetra():
    """Tolerance for TETRA concordance comparisons."""
    return 0.1


@pytest.mark.skip_if_exe_missing("nucmer")
def test_anim_concordance(
    paths_concordance_fna, path_concordance_jspecies, tolerance_anim, tmp_path
):
    """Check ANIm results are concordant with JSpecies."""
    # Perform ANIm on the input directory contents
    # We have to separate nucmer/delta-filter command generation
    # because Travis-CI doesn't play nicely with changes we made
    # for local SGE/OGE integration.
    # This might be avoidable with a scheduler flag passed to
    # jobgroup generation in the anim.py module. That's a TODO.
    ncmds, fcmds = anim.generate_nucmer_commands(paths_concordance_fna, tmp_path)
    (tmp_path / "nucmer_output").mkdir(exist_ok=True, parents=True)
    run_mp.multiprocessing_run(ncmds)

    # delta-filter commands need to be treated with care for
    # Travis-CI. Our cluster won't take redirection or semicolon
    # separation in individual commands, but the wrapper we wrote
    # for this (delta_filter_wrapper.py) can't be called under
    # Travis-CI. So we must deconstruct the commands below
    dfcmds = [
        " > ".join([" ".join(fcmd.split()[1:-1]), fcmd.split()[-1]]) for fcmd in fcmds
    ]
    run_mp.multiprocessing_run(dfcmds)

    orglengths = pyani_files.get_sequence_lengths(paths_concordance_fna)

    results = anim.process_deltadir(tmp_path / "nucmer_output", orglengths)
    result_pid = results.percentage_identity
    result_pid.to_csv(tmp_path / "pyani_anim.tab", sep="\t")

    # Compare JSpecies output to results
    result_pid = (result_pid.sort_index(axis=0).sort_index(axis=1) * 100.0).values
    tgt_pid = parse_jspecies(path_concordance_jspecies)["ANIm"].values

    assert result_pid - tgt_pid == pytest.approx(0, abs=tolerance_anim)


@pytest.mark.skip_if_exe_missing("blastn")
def test_anib_concordance(
    paths_concordance_fna,
    path_concordance_jspecies,
    tolerance_anib_hi,
    tolerance_anib_lo,
    threshold_anib_lo_hi,
    fragment_length,
    tmp_path,
):
    """Check ANIb results are concordant with JSpecies.

    We expect ANIb results to be quite different, as the BLASTN
    algorithm changed substantially between BLAST and BLAST+ (the
    megaBLAST algorithm is now the default for BLASTN)
    """
    # Get lengths of input genomes
    orglengths = pyani_files.get_sequence_lengths(paths_concordance_fna)

    # Build and run BLAST jobs
    fragfiles, fraglengths = anib.fragment_fasta_files(
        paths_concordance_fna, tmp_path, fragment_length
    )
    jobgraph = anib.make_job_graph(
        paths_concordance_fna, fragfiles, anib.make_blastcmd_builder(tmp_path)
    )
    assert 0 == run_mp.run_dependency_graph(jobgraph)  # Jobs must run correctly

    # Process BLAST output
    result_pid = anib.process_blast(
        tmp_path, orglengths, fraglengths
    ).percentage_identity

    # Compare JSpecies output to results. We do this in two blocks,
    # masked according to whether the expected result is greater than
    # a threshold separating "low" from "high" identity comparisons.
    result_pid = result_pid.sort_index(axis=0).sort_index(axis=1) * 100.0
    lo_result = result_pid.mask(result_pid >= threshold_anib_lo_hi).fillna(0).values
    hi_result = result_pid.mask(result_pid < threshold_anib_lo_hi).fillna(0).values

    tgt_pid = parse_jspecies(path_concordance_jspecies)["ANIb"]
    lo_target = tgt_pid.mask(tgt_pid >= threshold_anib_lo_hi).fillna(0).values
    hi_target = tgt_pid.mask(tgt_pid < threshold_anib_lo_hi).fillna(0).values

    assert (lo_result - lo_target, hi_result - hi_target) == (
        pytest.approx(0, abs=tolerance_anib_lo),
        pytest.approx(0, abs=tolerance_anib_hi),
    )


@pytest.mark.skip_if_exe_missing("blastall")
def test_aniblastall_concordance(
    paths_concordance_fna,
    path_concordance_jspecies,
    tolerance_anib_hi,
    fragment_length,
    tmp_path,
):
    """Check ANIblastall results are concordant with JSpecies."""
    # Get lengths of input genomes
    orglengths = pyani_files.get_sequence_lengths(paths_concordance_fna)

    # Perform ANIblastall on the input directory contents
    fragfiles, fraglengths = aniblastall.fragment_fasta_files(
        paths_concordance_fna, tmp_path, fragment_length
    )
    jobgraph = aniblastall.make_job_graph(
        paths_concordance_fna,
        fragfiles,
        aniblastall.make_blastcmd_builder(tmp_path),
    )
    assert 0 == run_mp.run_dependency_graph(jobgraph)  # Jobs must run correctly

    # Process BLAST output
    result_pid = aniblastall.process_blast(
        tmp_path, orglengths, fraglengths
    ).percentage_identity

    # Compare JSpecies output to results
    result_pid = (result_pid.sort_index(axis=0).sort_index(axis=1) * 100.0).values
    tgt_pid = parse_jspecies(path_concordance_jspecies)["ANIb"].values
    assert result_pid - tgt_pid == pytest.approx(0, abs=tolerance_anib_hi)


def test_tetra_concordance(
    paths_concordance_fna, path_concordance_jspecies, tolerance_tetra, tmp_path
):
    """Check TETRA results are concordant with JSpecies."""
    # Perform TETRA analysis
    zscores = dict()
    for filename in paths_concordance_fna:
        zscores[filename.stem] = tetra.calculate_tetra_zscore(filename)
    results = tetra.calculate_correlations(zscores).values

    # Compare JSpecies output
    tgt_mat = parse_jspecies(path_concordance_jspecies)["Tetra"].values
    assert results - tgt_mat == pytest.approx(0, abs=tolerance_tetra)
