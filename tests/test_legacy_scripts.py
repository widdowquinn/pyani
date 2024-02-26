#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019-2020
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
# Copyright (c) 2019-2020 University of Strathclyde
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
"""Test legacy scripts.

The test suite is intended to be run from the repository root using:

pytest -v

The two legacy scripts download genomes, then carry out ANI analysis. The
pytest ordering plug-in is used to guarantee that the download script
tests are conducted first.
"""

from argparse import Namespace
from pathlib import Path

import pytest

from pyani.pyani_config import (
    BLASTALL_DEFAULT,
    BLASTN_DEFAULT,
    FILTER_DEFAULT,
    FORMATDB_DEFAULT,
    MAKEBLASTDB_DEFAULT,
    NUCMER_DEFAULT,
)
from pyani.scripts import average_nucleotide_identity, genbank_get_genomes_by_taxon
from tools import modify_namespace


@pytest.fixture
def legacy_ani_namespace(path_fixtures_base, tmp_path):
    """Base namespace for legacy average_nucleotide_identity.py tests."""
    return Namespace(
        outdirname=tmp_path,
        indirname=path_fixtures_base / "legacy" / "ANI_input",
        verbose=False,
        debug=False,
        force=True,
        fragsize=1020,
        logfile=Path("test_ANIm.log"),
        skip_nucmer=False,
        skip_blastn=False,
        noclobber=False,
        nocompress=False,
        graphics=True,
        gformat="pdf,png",
        gmethod="seaborn",
        labels=path_fixtures_base / "legacy" / "ANI_input" / "labels.txt",
        classes=path_fixtures_base / "legacy" / "ANI_input" / "classes.txt",
        method="ANIm",
        scheduler="multiprocessing",
        workers=None,
        sgeargs=None,
        sgegroupsize=10000,
        maxmatch=False,
        nucmer_exe=NUCMER_DEFAULT,
        filter_exe=FILTER_DEFAULT,
        blastn_exe=BLASTN_DEFAULT,
        blastall_exe=BLASTALL_DEFAULT,
        makeblastdb_exe=MAKEBLASTDB_DEFAULT,
        formatdb_exe=FORMATDB_DEFAULT,
        write_excel=False,
        rerender=False,
        subsample=None,
        seed=None,
        jobprefix="ANI",
        tree=False,
    )


@pytest.fixture
def legacy_anib_sns_namespace(tmp_path, legacy_ani_namespace):
    """Namespace for legacy ANIm script tests.

    Uses the base namespace to run ANIm with seaborn output
    """
    return modify_namespace(legacy_ani_namespace, method="ANIb")


@pytest.fixture
def legacy_anib_mpl_namespace(tmp_path, legacy_ani_namespace):
    """Namespace for legacy ANIm script tests.

    Runs ANIm with matplotlib output
    """
    return modify_namespace(legacy_ani_namespace, gmethod="mpl", method="ANIb")


@pytest.fixture
def legacy_anim_sns_namespace(tmp_path, legacy_ani_namespace):
    """Namespace for legacy ANIm script tests.

    Uses the base namespace to run ANIm with seaborn output
    """
    return legacy_ani_namespace


@pytest.fixture
def legacy_anim_mpl_namespace(tmp_path, legacy_ani_namespace):
    """Namespace for legacy ANIm script tests.

    Runs ANIm with matplotlib output
    """
    return modify_namespace(legacy_ani_namespace, gmethod="mpl")


@pytest.fixture
def legacy_download_namespace(tmp_path):
    """Namespace for legacy download script tests."""
    return Namespace(
        outdirname=tmp_path,
        taxon="203804",
        verbose=False,
        force=True,
        noclobber=False,
        logfile=None,
        format="fasta",
        email="pyani@pyani.tests",
        retries=20,
        batchsize=10000,
        timeout=10,
        debug=False,
    )


@pytest.fixture
def legacy_tetra_sns_namespace(tmp_path, legacy_ani_namespace):
    """Namespace for legacy ANIm script tests.

    Uses the base namespace to run ANIm with seaborn output
    """
    return modify_namespace(legacy_ani_namespace, method="TETRA")


@pytest.fixture
def legacy_tetra_mpl_namespace(tmp_path, legacy_ani_namespace):
    """Namespace for legacy ANIm script tests.

    Uses the base namespace to run ANIm with mpl output
    """
    return modify_namespace(legacy_ani_namespace, method="TETRA", gmethod="mpl")


@pytest.mark.skip_if_exe_missing("nucmer")
@pytest.mark.slow
def test_legacy_anim_sns(legacy_anim_sns_namespace):
    r"""Use legacy script to run ANIm (seaborn output).

    average_nucleotide_identity.py \
        -l test_ANIm.log \
        -i tests/fixtures/legacy/ANI_input \
        -o tests/test_output/legacy_scripts/ANIm_seaborn \
        -g --gmethod seaborn --gformat pdf,png \
        -f --jobprefix ANI
    """
    average_nucleotide_identity.run_main(legacy_anim_sns_namespace)


@pytest.mark.skip_if_exe_missing("nucmer")
@pytest.mark.slow
def test_legacy_anim_mpl(legacy_anim_mpl_namespace):
    r"""Use legacy script to run ANIm (mpl output).

    average_nucleotide_identity.py \
        -l test_ANIm.log \
        -i tests/fixtures/legacy/ANI_input \
        -o tests/test_output/legacy_scripts/ANIm_mpl \
        -g --gmethod mpl --gformat pdf,png \
        -f --jobprefix ANI
    """
    average_nucleotide_identity.run_main(legacy_anim_mpl_namespace)


@pytest.mark.skip_if_exe_missing("blastn")
@pytest.mark.slow
def test_legacy_anib_sns(legacy_anib_sns_namespace):
    r"""Use legacy script to run ANIb (seaborn output).

    average_nucleotide_identity.py \
        -l test_ANIb.log \
        -i tests/test_output/legacy_scripts/C_blochmannia \
        -o tests/test_output/legacy_scripts/ANIb_seaborn \
        -g --gmethod seaborn --gformat pdf,png \
        -f --jobprefix ANI
    """
    average_nucleotide_identity.run_main(legacy_anib_sns_namespace)


@pytest.mark.skip_if_exe_missing("blastn")
@pytest.mark.slow
def test_legacy_anib_mpl(legacy_anib_mpl_namespace):
    r"""Use legacy script to run ANIb (mpl output).

    average_nucleotide_identity.py \
        -l test_ANIb.log \
        -i tests/test_output/legacy_scripts/C_blochmannia \
        -o tests/test_output/legacy_scripts/ANIb_mpl \
        -g --gmethod mpl --gformat pdf,png \
        -f --jobprefix ANI
    """
    average_nucleotide_identity.run_main(legacy_anib_mpl_namespace)


@pytest.mark.slow
def test_legacy_tetra_sns(legacy_tetra_sns_namespace):
    r"""Use legacy script to run TETRA (seaborn output)."""
    average_nucleotide_identity.run_main(legacy_tetra_sns_namespace)


@pytest.mark.slow
def test_legacy_tetra_mpl(legacy_tetra_mpl_namespace):
    r"""Use legacy script to run TETRA (mpl output)."""
    average_nucleotide_identity.run_main(legacy_tetra_mpl_namespace)


def test_legacy_genome_downloads(
    legacy_download_namespace, mock_legacy_single_genome_dl
):
    r"""Use legacy script to download genomes; mocks file downloading.

    Otherwise emulates a command such as:

    genbank_get_genomes_by_taxon.py \
        -o tests/fixtures/legacy/ANI_input \
        --email emailme@my.email.domain \
        -t 203804 -f
    """
    genbank_get_genomes_by_taxon.run_main(legacy_download_namespace)
