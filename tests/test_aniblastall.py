"""Test anib.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path
from typing import List, NamedTuple

import pandas as pd
import pytest  # noqa: F401  # pylint: disable=unused-import

from pandas.util.testing import assert_frame_equal

from pyani import anib, pyani_files


class ANIbOutput(NamedTuple):

    """Convenience struct for ANIb output."""

    fragfile: Path
    tabfile: Path
    legacytabfile: Path


class ANIbOutputDir(NamedTuple):

    """Convenience struct for ANIb output."""

    infiles: List[Path]
    fragfiles: List[Path]
    blastdir: Path
    legacyblastdir: Path
    blastresult: pd.DataFrame
    legacyblastresult: pd.DataFrame


@pytest.fixture
def anib_output(dir_anib_in):
    """Namedtuple of example ANIb output.

    fragfile - fragmented FASTA query file
    tabfile  - BLAST+ tabular output
    legacytabfile - blastall tabular output
    """
    return ANIbOutput(
        dir_anib_in / "NC_002696-fragments.fna",
        dir_anib_in / "NC_002696_vs_NC_011916.blast_tab",
        dir_anib_in / "NC_002696_vs_NC_010338.blast_tab",
    )


@pytest.fixture
def anib_output_dir(dir_anib_in):
    """Namedtuple of example ANIb output - full directory.

    infiles - list of FASTA query files
    fragfiles - list of fragmented FASTA query files
    blastdir - path to BLAST+ output data
    legacyblastdir - path to blastall output data
    blastresult - pd.DataFrame result for BLAST+
    legacyblastresult - pd.DataFrame result for blastall
    """
    return ANIbOutputDir(
        [
            _
            for _ in (dir_anib_in / "sequences").iterdir()
            if _.is_file() and _.suffix == ".fna"
        ],
        [
            _
            for _ in (dir_anib_in / "fragfiles").iterdir()
            if _.is_file() and _.suffix == ".fna"
        ],
        dir_anib_in / "blastn",
        dir_anib_in / "blastall",
        pd.read_csv(dir_anib_in / "dataframes" / "blastn_result.csv", index_col=0),
        pd.read_csv(dir_anib_in / "dataframes" / "blastall_result.csv", index_col=0),
    )


# Test legacy BLAST (blastall) command generation
def test_blastall_dbjobdict(path_fna_all, tmp_path):
    """Generate dictionary of legacy BLASTN database jobs."""
    blastcmds = anib.make_blastcmd_builder("ANIblastall", tmp_path)
    jobdict = anib.build_db_jobs(path_fna_all, blastcmds)
    expected = [
        (tmp_path / _.name, f"formatdb -p F -i {tmp_path / _.name} -t {_.stem}")
        for _ in path_fna_all
    ]
    assert sorted([(k, v.script) for (k, v) in jobdict.items()]) == sorted(expected)


def test_blastall_graph(path_fna_all, tmp_path, fragment_length):
    """Create jobgraph for legacy BLASTN jobs."""
    fragresult = anib.fragment_fasta_files(path_fna_all, tmp_path, fragment_length)
    blastcmds = anib.make_blastcmd_builder("ANIblastall", tmp_path)
    jobgraph = anib.make_job_graph(path_fna_all, fragresult[0], blastcmds)
    # We check that the main script job is a blastn job, and that there
    # is a single dependency, which is a makeblastdb job
    for job in jobgraph:
        assert job.script.startswith("blastall -p blastn")
        assert len(job.dependencies) == 1
        assert job.dependencies[0].script.startswith("formatdb")


def test_blastall_multiple(path_fna_two, tmp_path):
    """Generate legacy BLASTN commands."""
    cmds = anib.generate_blastn_commands(path_fna_two, tmp_path, mode="ANIblastall")
    expected = [
        (
            "blastall -p blastn -o "
            f"{tmp_path / str(path_fna_two[0].stem + '_vs_' + path_fna_two[1].stem + '.blast_tab')} "
            f"-i {path_fna_two[0]} "
            f"-d {path_fna_two[1]} "
            "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
        ),
        (
            "blastall -p blastn -o "
            f"{tmp_path / str(path_fna_two[1].stem + '_vs_' + path_fna_two[0].stem + '.blast_tab')} "
            f"-i {path_fna_two[1]} "
            f"-d {path_fna_two[0]} "
            "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
        ),
    ]
    assert cmds == expected


def test_blastall_single(path_fna_two, tmp_path):
    """Generate legacy BLASTN command-line."""
    cmd = anib.construct_blastall_cmdline(path_fna_two[0], path_fna_two[1], tmp_path)
    expected = (
        f"blastall -p blastn -o {tmp_path / str(path_fna_two[0].stem + '_vs_' + path_fna_two[1].stem + '.blast_tab')} "
        f"-i {path_fna_two[0]} "
        f"-d {path_fna_two[1]} "
        "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
    )
    assert cmd == expected


# Test legacy BLAST database formatting (formatdb) command generation
def test_formatdb_multiple(path_fna_two, tmp_path):
    """Generate legacy BLAST db creation commands."""
    cmds = anib.generate_blastdb_commands(path_fna_two, tmp_path, mode="ANIblastall")
    expected = [
        (
            f"formatdb -p F -i {tmp_path / path_fna_two[0].name} -t {path_fna_two[0].stem}",
            tmp_path / path_fna_two[0].name,
        ),
        (
            f"formatdb -p F -i {tmp_path / path_fna_two[1].name} -t {path_fna_two[1].stem}",
            tmp_path / path_fna_two[1].name,
        ),
    ]
    assert cmds == expected


def test_formatdb_single(path_fna, tmp_path):
    """Generate legacy BLAST formatdb command-line."""
    cmd = anib.construct_formatdb_cmd(path_fna, tmp_path)
    expected = f"formatdb -p F -i {tmp_path / path_fna.name} -t {path_fna.stem}"
    assert cmd[0] == expected


# Test output file parsing for ANIb methods
def test_parse_legacy_blastdir(anib_output_dir):
    """Parses directory of legacy BLAST output."""
    orglengths = pyani_files.get_sequence_lengths(anib_output_dir.infiles)
    fraglengths = anib.get_fraglength_dict(anib_output_dir.fragfiles)
    result = anib.process_blast(
        anib_output_dir.legacyblastdir, orglengths, fraglengths, mode="ANIblastall"
    )
    assert_frame_equal(
        result.percentage_identity.sort_index(1).sort_index(),
        anib_output_dir.legacyblastresult.sort_index(1).sort_index(),
    )


def test_parse_legacy_blasttab(anib_output):
    """Parses ANIB legacy .blast_tab output."""
    fragdata = anib.get_fraglength_dict([anib_output.fragfile])
    result = anib.parse_blast_tab(
        anib_output.legacytabfile, fragdata, mode="ANIblastall"
    )
    assert (
        a == b for a, b in zip(result, [1_966_922, 406_104, 78.578_978_313_253_018])
    )
