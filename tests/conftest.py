#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The University of Strathclude 2019-2021
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
# 161 Cathedral Street
# Glasgow
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# (c) The University of Strathclude 2019-2021
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
"""Pytest configuration file."""

import subprocess
import shutil
import os
import re

from pathlib import Path
from typing import NamedTuple

import pytest

from pyani import download
from pyani.download import ASMIDs, DLStatus
from pyani.pyani_config import (
    BLASTALL_DEFAULT,
    BLASTN_DEFAULT,
    NUCMER_DEFAULT,
    FRAGSIZE,
)
from pyani.scripts import genbank_get_genomes_by_taxon

# Path to tests, contains tests and data subdirectories
# This conftest.py file should be found in the top directory of the tests
# module. The fixture data should be in a subdirectory named fixtures
TESTSPATH = Path(__file__).parents[0]
FIXTUREPATH = TESTSPATH / "fixtures"


# Convenience structs to emulate returned objects
class MockProcess(NamedTuple):
    """Mock process object."""

    stdout: str
    stderr: str


class MockMatch(NamedTuple):
    """Mock match object."""

    def group(self):
        return ""


@pytest.fixture
def blastall_available():
    """Returns True if blastall can be run, False otherwise."""
    cmd = str(BLASTALL_DEFAULT)
    # Can't use check=True, as blastall without arguments returns 1!
    try:
        result = subprocess.run(
            cmd,
            shell=False,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError:
        return False
    return result.stdout[1:9] == b"blastall"


@pytest.fixture
def blastn_available():
    """Returns True if blastn can be run, False otherwise."""
    cmd = [str(BLASTN_DEFAULT), "-version"]
    try:
        result = subprocess.run(
            cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
        )
    except OSError:
        return False
    return result.stdout[:6] == b"blastn"


@pytest.fixture
def dir_anib_in():
    """Input files for ANIb tests."""
    return FIXTUREPATH / "anib"


@pytest.fixture
def dir_aniblastall_in():
    """Input files for ANIblastall tests."""
    return FIXTUREPATH / "anib"


@pytest.fixture
def dir_anim_in():
    """Input files for ANIm tests."""
    return FIXTUREPATH / "anim"


@pytest.fixture
def dir_graphics_in():
    """Input files for graphics tests."""
    return FIXTUREPATH / "graphics"


@pytest.fixture
def dir_seq():
    """Sequence files for tests."""
    return FIXTUREPATH / "sequences"


@pytest.fixture
def dir_targets():
    """Target files for output comparisons."""
    return FIXTUREPATH / "targets"


@pytest.fixture
def dir_tgt_fragments(dir_targets):
    """Target files for FASTA file fragmentation."""
    return dir_targets / "fragments"


@pytest.fixture
def email_address():
    """Dummy email address."""
    return "pyani.tests@pyani.org"


@pytest.fixture
def executable_incompatible_with_os(monkeypatch):
    """
    Mocks an executable file that is incompatible with the OS.

    (This situation likely only applies to blastall.)
    """

    def mock_which(*args, **kwargs):
        """Mock an absolute file path."""
        return args[0]

    def mock_isfile(*args, **kwargs):
        """Mock a call to `os.path.isfile()`."""
        return True

    def mock_access(*args, **kwargs):
        """Mock a call to `os.access()`."""
        return True

    def mock_subprocess(*args, **kwargs):
        """Mock a call to `subprocess.run()` with an incompatible program."""
        raise OSError

    monkeypatch.setattr(shutil, "which", mock_which)
    monkeypatch.setattr(Path, "is_file", mock_isfile)
    monkeypatch.setattr(os.path, "isfile", mock_isfile)
    monkeypatch.setattr(os, "access", mock_access)
    monkeypatch.setattr(subprocess, "run", mock_subprocess)


@pytest.fixture
def executable_missing(monkeypatch):
    """Mocks an executable path that does not point to a file."""

    def mock_which(*args, **kwargs):
        """Mock a call to `shutil.which()`, which produces an absolute file path."""
        return args[0]

    def mock_isfile(*args, **kwargs):
        """Mock a call to `os.path.isfile()`."""
        return False

    monkeypatch.setattr(shutil, "which", mock_which)  # Path(test_file_1))
    monkeypatch.setattr(Path, "is_file", mock_isfile)
    monkeypatch.setattr(os.path, "isfile", mock_isfile)


@pytest.fixture
def executable_not_executable(monkeypatch):
    """
    Mocks an executable path that does not point to an executable file,
    but does point to a file.
    """

    def mock_which(*args, **kwargs):
        """Mock an absolute file path."""
        return args[0]

    def mock_isfile(*args, **kwargs):
        """Mock a call to `os.path.isfile()`."""
        return True

    def mock_access(*args, **kwargs):
        """Mock a call to `os.access()`."""
        return False

    monkeypatch.setattr(shutil, "which", mock_which)
    monkeypatch.setattr(Path, "is_file", mock_isfile)
    monkeypatch.setattr(os.path, "isfile", mock_isfile)
    monkeypatch.setattr(os, "access", mock_access)


@pytest.fixture
def executable_without_version(monkeypatch):
    """
    Mocks an executable file for which the version can't be obtained, but
    which runs without incident.
    """

    def mock_which(*args, **kwargs):
        """Mock an absolute file path."""
        return args[0]

    def mock_isfile(*args, **kwargs):
        """Mock a call to `os.path.isfile()`."""
        return True

    def mock_access(*args, **kwargs):
        """Mock a call to `os.access()`."""
        return True

    def mock_subprocess(*args, **kwargs):
        """Mock a call to `subprocess.run()`."""
        return MockProcess(b"mock bytes", b"mock bytes")

    def mock_search(*args, **kwargs):
        """Mock a call to `re.search()`."""
        return MockMatch()

    monkeypatch.setattr(shutil, "which", mock_which)
    monkeypatch.setattr(Path, "is_file", mock_isfile)
    monkeypatch.setattr(os.path, "isfile", mock_isfile)
    monkeypatch.setattr(os, "access", mock_access)
    monkeypatch.setattr(subprocess, "run", mock_subprocess)
    monkeypatch.setattr(re, "search", mock_search)


@pytest.fixture
def fragment_length():
    """Fragment size for ANIb-related analyses."""
    return FRAGSIZE


@pytest.fixture
def mock_legacy_single_genome_dl(monkeypatch):
    """Mocks remote database calls for single-genome downloads.

    This masks calls to functions in genbank_get_genomes_by_taxon, for safe testing.

    This will be deprecated once the genbank_get_genomes_by_taxon.py script is
    converted to use the pyani.download module.
    """

    def mock_asmuids(*args, **kwargs):
        """Mock genbank_get_genomes_by_taxon.get_asm_uids()."""
        return ["32728"]

    def mock_ncbi_asm(*args, **kwargs):
        """Mock genbank_get_genomes_by_taxon.get_ncbi_asm()."""
        return (
            Path(
                "tests/test_output/legacy_scripts/C_blochmannia_legacy/GCF_000011605.1_ASM1160v1_genomic.fna"
            ),
            "8b0cab310cb638c977d453ff06eceb64\tGCF_000011605.1_ASM1160v1_genomic\tPectobacterium atrosepticum",
            "8b0cab310cb638c977d453ff06eceb64\tGCF_000011605.1_ASM1160v1_genomic\tP. atrosepticum SCRI1043",
            "GCF_000011605.1",
        )

    monkeypatch.setattr(genbank_get_genomes_by_taxon, "get_asm_uids", mock_asmuids)
    monkeypatch.setattr(genbank_get_genomes_by_taxon, "get_ncbi_asm", mock_ncbi_asm)


@pytest.fixture
def mock_single_genome_dl(monkeypatch):
    """Mocks remote database calls for single-genome downloads.

    This masks calls to the download module, for safe testing.
    """

    def mock_asmuids(*args, **kwargs):
        """Mock download.get_asm_uids()."""
        return ASMIDs("txid218491[Organism:exp]", 1, ["32728"])

    def mock_ncbi_esummary(*args, **kwargs):
        """Mock download.get_ncbi_esummary()."""
        return (
            {
                "Taxid": "218491",
                "SpeciesTaxid": "29471",
                "AssemblyAccession": "GCF_000011605.1",
                "AssemblyName": "ASM1160v1",
                "SpeciesName": "Pectobacterium atrosepticum",
            },
            "GCF_000011605.1_ASM1160v1",
        )

    def mock_genome_hash(*args, **kwargs):
        """Mock download.retrieve_genome_and_hash()."""
        return DLStatus(
            "ftp://ftp.ncbi.nlm.nih.gov/dummy_genomic.fna.gz",
            "ftp://ftp.ncbi.nlm.nih.gov/dummy/md5checksums.txt",
            FIXTUREPATH
            / "single_genome_download"
            / "GCF_000011605.1_ASM1160v1_genomic.fna.gz",
            FIXTUREPATH / "single_genome_download/GCF_000011605.1_ASM1160v1_hashes.txt",
            False,
            None,
        )

    monkeypatch.setattr(download, "get_asm_uids", mock_asmuids)
    monkeypatch.setattr(download, "get_ncbi_esummary", mock_ncbi_esummary)
    monkeypatch.setattr(download, "retrieve_genome_and_hash", mock_genome_hash)


@pytest.fixture
def nucmer_available():
    """Test that nucmer is available."""
    cmd = [str(NUCMER_DEFAULT), "--version"]
    try:
        result = subprocess.run(
            cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
        )
    except OSError:
        return False
    return result.stderr[:6] == b"nucmer"


@pytest.fixture
def path_file_two():
    """Path to two arbitrary filenames."""
    return [Path(f"file{_:d}.fna") for _ in range(1, 3)]


@pytest.fixture
def path_file_four():
    """Path to four arbitrary filenames."""
    return [Path(f"file{_:d}.fna") for _ in range(1, 5)]


@pytest.fixture
def path_fixtures_base():
    """Base path to fixture data folders."""
    return FIXTUREPATH


@pytest.fixture
def path_fna(dir_seq):
    """Path to one .fna sequence file from dir_seq."""
    fnapaths = [_ for _ in dir_seq.iterdir() if _.is_file() and _.suffix == ".fna"]
    return fnapaths[0]


@pytest.fixture
def path_fna_two(dir_seq):
    """Paths to two .fna sequence file in dir_seq."""
    fnapaths = [_ for _ in dir_seq.iterdir() if _.is_file() and _.suffix == ".fna"]
    return fnapaths[:2]


@pytest.fixture
def path_fna_all(dir_seq):
    """Paths to all .fna sequence file in dir_seq."""
    return [_ for _ in dir_seq.iterdir() if _.is_file() and _.suffix == ".fna"]


@pytest.fixture(autouse=True)
def skip_by_unavailable_executable(
    request, blastall_available, blastn_available, nucmer_available
):
    """Skip test if executable is unavailable.

    Use with @pytest.mark.skip_if_exe_missing("executable") decorator.
    """
    if request.node.get_closest_marker("skip_if_exe_missing"):
        exe_name = request.node.get_closest_marker("skip_if_exe_missing").args[0]
        tests = {
            "blastall": blastall_available,
            "blastn": blastn_available,
            "nucmer": nucmer_available,
        }
        try:
            if not tests[exe_name]:
                pytest.skip(f"Skipped as {exe_name} not available")
        except KeyError:  # Unknown executables are ignored
            pytest.skip(f"Executable {exe_name} not recognised")
