#!/usr/bin/env python

"""Tests for pyani package scripts

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/), from the repository root
directory.

If the test is run directly at the command-line, the output obtained by each
test is returned to STDOUT.
"""


from nose.tools import assert_equal, assert_less

import os
import pandas as pd
import shutil
import subprocess
import sys

# Thresholds for testing
ANIM_THRESHOLD = 1e-3
ANIB_THRESHOLD = 1e-1
ANIBLASTALL_THRESHOLD = 1e-1
TETRA_THRESHOLD = 1e-3


# Work out where we are. We need to do this to find related data files
# for testing
curdir = os.path.dirname(os.path.abspath(__file__))


def test_anim_script():
    """Test average_nucleotide_identity.py script with ANIm"""
    # Remove output directory
    outdir = os.path.join("tests", "test_script_ANIm")
    try:
        shutil.rmtree(outdir, ignore_errors=True)
    except FileNotFoundError:
        print("Did not find %s to delete it (not an error, continuing)")
        pass    

    # Run ANIm script
    indir = os.path.join("tests", "test_ani_data")
    logfile = os.path.join("tests", "test_ANIm_script.log")
    cmd = ' '.join(["python",
                    "average_nucleotide_identity.py",
                    "-v",
                    "-i %s" % indir,
                    "-o %s" % outdir,
                    "-l %s" % logfile,
                    "--classes %s" % os.path.join(indir, "classes.tab"),
                    "--labels %s" % os.path.join(indir, "labels.tab"),
                    "-g",
                    "--gformat png,pdf,svg",
                    "-m ANIm"])
    print(cmd)
    result = subprocess.run(cmd, shell=sys.platform != "win32",
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    
    # Test result
    expected_file = os.path.join("tests", "target_ANIm_output",
                                 "ANIm_percentage_identity.tab")
    returned_file = os.path.join(outdir, "ANIm_percentage_identity.tab")
    expected = pd.read_csv(expected_file, index_col=0,
                           sep="\t").sort_index(axis=0).sort_index(axis=1)
    index, columns = expected.index, expected.columns
    expected = expected.as_matrix() * 100
    print(expected)
    returned = pd.read_csv(returned_file, index_col=0,
                           sep="\t").sort_index(axis=0).sort_index(axis=1)
    returned = returned.as_matrix() * 100
    print(returned)
    diffmat = pd.DataFrame(expected - returned,
                           index=index, columns=columns)
    max_diff = diffmat.abs().values.max()
    print("Maximum observed difference: %e" % max_diff)
    assert_less(max_diff, ANIM_THRESHOLD)


def test_anib_script():
    """Test average_nucleotide_identity.py script with ANIb"""
    # Remove output directory
    outdir = os.path.join("tests", "test_script_ANIb")
    try:
        shutil.rmtree(outdir, ignore_errors=True)
    except FileNotFoundError:
        print("Did not find %s to delete it (not an error, continuing)")
        pass    

    # Run ANIb script
    indir = os.path.join("tests", "test_ani_data")
    logfile = os.path.join("tests", "test_ANIb_script.log")
    cmd = ' '.join(["python",
                    "average_nucleotide_identity.py",
                    "-v",
                    "-i %s" % indir,
                    "-o %s" % outdir,
                    "-l %s" % logfile,
                    "--classes %s" % os.path.join(indir, "classes.tab"),
                    "--labels %s" % os.path.join(indir, "labels.tab"),
                    "-g",
                    "--gmethod seaborn",
                    "-m ANIb"])
    print(cmd)
    result = subprocess.run(cmd, shell=sys.platform != "win32",
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    
    # Test result
    expected_file = os.path.join("tests", "target_ANIb_output",
                                 "ANIb_percentage_identity.tab")
    returned_file = os.path.join(outdir, "ANIb_percentage_identity.tab")
    expected = pd.read_csv(expected_file, index_col=0,
                           sep="\t").sort_index(axis=0).sort_index(axis=1)
    index, columns = expected.index, expected.columns
    expected = expected.as_matrix() * 100
    print(expected)
    returned = pd.read_csv(returned_file, index_col=0,
                           sep="\t").sort_index(axis=0).sort_index(axis=1)
    returned = returned.as_matrix() * 100
    print(returned)
    diffmat = pd.DataFrame(expected - returned,
                           index=index, columns=columns)
    max_diff = diffmat.abs().values.max()
    print("Maximum observed difference: %e" % max_diff)
    assert_less(max_diff, ANIB_THRESHOLD)


def test_aniblastall_script():
    """Test average_nucleotide_identity.py script with ANIblastall"""
    # Remove output directory
    outdir = os.path.join("tests", "test_script_ANIblastall")
    try:
        shutil.rmtree(outdir, ignore_errors=True)
    except FileNotFoundError:
        print("Did not find %s to delete it (not an error, continuing)")
        pass    

    # Run ANIblastall script
    indir = os.path.join("tests", "test_ani_data")
    logfile = os.path.join("tests", "test_ANIblastall_script.log")
    cmd = ' '.join(["python",
                    "average_nucleotide_identity.py",
                    "-v",
                    "-i %s" % indir,
                    "-o %s" % outdir,
                    "-l %s" % logfile,
                    "--classes %s" % os.path.join(indir, "classes.tab"),
                    "--labels %s" % os.path.join(indir, "labels.tab"),
                    "-m ANIblastall"])
    print(cmd)
    result = subprocess.run(cmd, shell=sys.platform != "win32",
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    
    # Test result
    expected_file = os.path.join("tests", "target_ANIblastall_output",
                                 "ANIblastall_percentage_identity.tab")
    returned_file = os.path.join(outdir, "ANIblastall_percentage_identity.tab")
    expected = pd.read_csv(expected_file, index_col=0,
                           sep="\t").sort_index(axis=0).sort_index(axis=1)
    index, columns = expected.index, expected.columns
    expected = expected.as_matrix() * 100
    print(expected)
    returned = pd.read_csv(returned_file, index_col=0,
                           sep="\t").sort_index(axis=0).sort_index(axis=1)
    returned = returned.as_matrix() * 100
    print(returned)
    diffmat = pd.DataFrame(expected - returned,
                           index=index, columns=columns)
    max_diff = diffmat.abs().values.max()
    print("Maximum observed difference: %e" % max_diff)
    assert_less(max_diff, ANIBLASTALL_THRESHOLD)
    

def test_tetra_script():
    """Test average_nucleotide_identity.py script with TETRA"""
    # Remove output directory
    outdir = os.path.join("tests", "test_script_TETRA")
    try:
        shutil.rmtree(outdir, ignore_errors=True)
    except FileNotFoundError:
        print("Did not find %s to delete it (not an error, continuing)")
        pass    

    # Run TETRA script
    indir = os.path.join("tests", "test_ani_data")
    logfile = os.path.join("tests", "test_tetra_script.log")
    cmd = ' '.join(["python",
                    "average_nucleotide_identity.py",
                    "-v",
                    "-i %s" % indir,
                    "-o %s" % outdir,
                    "-l %s" % logfile,
                    "--classes %s" % os.path.join(indir, "classes.tab"),
                    "--labels %s" % os.path.join(indir, "labels.tab"),
                    "-m TETRA"])
    print(cmd)
    result = subprocess.run(cmd, shell=sys.platform != "win32",
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    
    # Test result
    expected_file = os.path.join("tests", "target_TETRA_output",
                                 "TETRA_correlations.tab")
    returned_file = os.path.join(outdir, "TETRA_correlations.tab")
    expected = pd.read_csv(expected_file, index_col=0,
                           sep="\t").sort_index(axis=0).sort_index(axis=1)
    index, columns = expected.index, expected.columns
    expected = expected.as_matrix() * 100
    print(expected)
    returned = pd.read_csv(returned_file, index_col=0,
                           sep="\t").sort_index(axis=0).sort_index(axis=1)
    returned = returned.as_matrix() * 100
    print(returned)
    diffmat = pd.DataFrame(expected - returned,
                           index=index, columns=columns)
    max_diff = diffmat.abs().values.max()
    print("Maximum observed difference: %e" % max_diff)
    assert_less(max_diff, TETRA_THRESHOLD)


def test_genbank_download():
    """Test GenBank download script."""
    # Remove output directory
    outdir = os.path.join("tests", "test_genbank_download")
    try:
        shutil.rmtree(outdir, ignore_errors=True)
    except FileNotFoundError:
        print("Did not find %s to delete it (not an error, continuing)")
        pass

    # Run GenBank download script
    cmd = ' '.join(["python",
                    "genbank_get_genomes_by_taxon.py",
                    "-o %s" % outdir,
                    "--email pyani@github.com.dev.null",
                    "-l tests/test_genbank_download.log", 
                    "-t 1224150",
                    "-v"])
    print(cmd)
    result = subprocess.run(cmd, shell=sys.platform != "win32",
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # Test result
    expected = sorted(["GCF_000400505.1_DPA2511_1.0_genomic.fna",
                       "GCF_000400505.1_DPA2511_1.0_genomic.fna.gz",
                       "classes.txt",
                       "labels.txt"])
    returned = sorted(os.listdir(outdir))
    assert_equal(expected, returned)
