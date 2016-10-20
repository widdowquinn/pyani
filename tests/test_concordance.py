#!/usr/bin/env python

"""Tests for pyani package concordance with JSpecies output

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/).

If the test is run directly at the command-line, the output obtained by each
test is returned to STDOUT.
"""
from nose.tools import assert_equal, assert_less, nottest
from pyani.run_multiprocessing import run_dependency_graph, \
    multiprocessing_run

import os
import pandas as pd
import shutil
import subprocess
import sys

from pyani import anib, anim, tetra, pyani_files, pyani_config

# Work out where we are. We need to do this to find related data files
# for testing
curdir = os.path.dirname(os.path.abspath(__file__))

# Path to JSpecies output data. This data is pre-prepared. If you replace
# the test data with your own data, you will need to replace this file,
# or change the file path.
JSPECIES_OUTFILE = os.path.relpath(os.path.join(curdir, 'test_JSpecies',
                                               'jspecies_results.tab'))

# Path to test input data
INDIRNAME = os.path.relpath(os.path.join(curdir, 'test_ani_data'))

# Path to directory for concordance test output
OUTDIRNAME = os.path.relpath(os.path.join(curdir, 'test_concordance'))

# Thresholds for allowable difference
TETRA_THRESHOLD = 0.1
ANIB_THRESHOLD = 2  # This threshold higher because BLASTN+ != legacy BLASTN
ANIBLASTALL_THRESHOLD = 0.1
ANIM_THRESHOLD = 0.1

# Helper function to parse tables from JSpecies output
def parse_table(filename, title):
    """Parse method output from JSpecies tabular output, returns dataframe.

    - filename - path to JSpecies output file
    - title - header of ANI method output

    The JSpecies output may have one or more tables; which table is required
    is indicated by the passed title.
    """
    assert title in ('Tetra', 'ANIb', 'ANIm'), "Invalid header: %s" % title
    with open(filename, 'rU') as fh:
        header, in_table = False, False
        for line in iter(fh):
            if line.strip() == title:  # We have the table we want
                header = True
            elif header:
                columns = line.strip().split('\t')
                data = pd.DataFrame(index=columns, columns=columns)
                in_table, header = True, False
            elif in_table:
                if not len(line.strip()):
                    return data.sort_index(axis=0).sort_index(axis=1)
                else:
                    ldata = line.strip().split('\t')
                    row = ldata[0]
                    for idx, val in enumerate(ldata[1:]):
                        if val != '---':
                            data[columns[idx]][row] = float(val)
                        elif title in ("ANIb", "ANIm"):
                            data[columns[idx]][row] = 100.0
                        else:
                            data[columns[idx]][row] = 1.0
            else:
                pass
    return data.sort_index(axis=0).sort_index(axis=1)


# Make output directory if necessary
def delete_and_remake_outdir(mode):
    """Make concordance test output directories."""
    outdirname = '_'.join([OUTDIRNAME, mode])
    try:
        shutil.rmtree(outdirname, ignore_errors=True)
    except FileNotFoundError:
        print("Did not find %s to delete it (not an error, continuing)")
        pass
    os.makedirs(outdirname, exist_ok=True)
    return outdirname


# Test concordance of ANIb with JSpecies output
@nottest
def test_anib_concordance():
    """Test concordance of ANIb method with JSpecies output.

    This may take some time. Please be patient.
    """
    # Make/check output directories
    mode = "ANIb"
    outdirname = delete_and_remake_outdir(mode)

    # Get dataframes of JSpecies output
    anib_jspecies = parse_table(JSPECIES_OUTFILE, 'ANIb')

    # Identify our input files, and the total lengths of each organism seq
    infiles = pyani_files.get_fasta_files(INDIRNAME)
    org_lengths = pyani_files.get_sequence_lengths(infiles)

    # Test ANIb concordance:
    # Make fragments
    fragfiles, fraglengths = anib.fragment_FASTA_files(infiles, outdirname,
                                                       pyani_config.FRAGSIZE)

    
    # Build jobgraph
    jobgraph = anib.make_job_graph(infiles, fragfiles, outdirname,
                                   mode="ANIb")
    print("\nJobgraph:\n", jobgraph)

    # Run jobgraph with multiprocessing
    run_dependency_graph(jobgraph)
    print("Ran multiprocessing jobs")

    # Process BLAST; the pid data is in anib_data[1]
    anib_data = anib.process_blast(outdirname, org_lengths, fraglengths,
                                   mode="ANIb")
    anib_pid = anib_data[1].sort_index(axis=0).sort_index(axis=1) * 100.

    index, columns = anib_pid.index, anib_pid.columns
    diffmat = anib_pid.as_matrix() - anib_jspecies.as_matrix()
    anib_diff = pd.DataFrame(diffmat, index=index, columns=columns)

    # Write dataframes to file, for reference
    anib_pid.to_csv(os.path.join(outdirname,
                                'ANIb_pid.tab'),
                   sep='\t')
    anib_jspecies.to_csv(os.path.join(outdirname,
                                      'ANIb_jspecies.tab'),
                         sep='\t')
    anib_diff.to_csv(os.path.join(outdirname,
                                  'ANIb_diff.tab'),
                     sep='\t')
    print("ANIb concordance test output placed in %s" % outdirname)
    print("ANIb PID: \n", anib_pid)
    print("ANIb JSpecies: \n", anib_jspecies)
    print("ANIb diff: \n", anib_diff)

    # We'd like the absolute difference reported to be < ANIB_THRESHOLD
    max_diff = anib_diff.abs().values.max()
    print("Maximum difference for ANIb: %e" % max_diff)
    assert_less(max_diff, ANIB_THRESHOLD)


# Test concordance of ANIblastall with JSpecies output
@nottest
def test_aniblastall_concordance():
    """Test concordance of ANIblastall method with JSpecies output."""
    # Make/check output directory
    mode = "ANIblastall"
    outdirname = delete_and_remake_outdir(mode)

    # Get dataframes of JSpecies output
    aniblastall_jspecies = parse_table(JSPECIES_OUTFILE, 'ANIb')

    # Identify our input files, and the total lengths of each organism seq
    infiles = pyani_files.get_fasta_files(INDIRNAME)
    org_lengths = pyani_files.get_sequence_lengths(infiles)

    # Test ANIblastall concordance:
    # Make fragments
    fragfiles, fraglengths = anib.fragment_FASTA_files(infiles, outdirname,
                                                       pyani_config.FRAGSIZE)

    # Build jobgraph
    jobgraph = anib.make_job_graph(infiles, fragfiles, outdirname,
                                   mode="ANIblastall")
    print("\nJobgraph:\n", jobgraph)
    print("\nJob 0:\n", jobgraph[0].script)

    # Run jobgraph with multiprocessing
    run_dependency_graph(jobgraph)
    print("Ran multiprocessing jobs")

    # Process BLAST; the pid data is in anib_data[1]
    aniblastall_data = anib.process_blast(outdirname, org_lengths,
                                          fraglengths,
                                          mode="ANIblastall")
    aniblastall_pid = \
        aniblastall_data[1].sort_index(axis=0).sort_index(axis=1) * 100.

    index, columns = aniblastall_pid.index, aniblastall_pid.columns
    diffmat = aniblastall_pid.as_matrix() - aniblastall_jspecies.as_matrix()
    aniblastall_diff = pd.DataFrame(diffmat, index=index, columns=columns)

    # Write dataframes to file, for reference
    aniblastall_pid.to_csv(os.path.join(outdirname,
                                        'ANIblastall_pid.tab'),
                           sep='\t')
    aniblastall_jspecies.to_csv(os.path.join(outdirname,
                                             'ANIblastall_jspecies.tab'),
                                sep='\t')
    aniblastall_diff.to_csv(os.path.join(outdirname,
                                  'ANIblastall_diff.tab'),
                     sep='\t')
    print("ANIblastall concordance test output placed in %s" % outdirname)
    print("ANIblastall PID:\n", aniblastall_pid)
    print("ANIblastall JSpecies:\n", aniblastall_jspecies)
    print("ANIblastall diff:\n", aniblastall_diff)

    # We'd like the absolute difference reported to be < ANIBLASTALL_THRESHOLD
    max_diff = aniblastall_diff.abs().values.max()
    print("Maximum difference for ANIblastall: %e" % max_diff)
    assert_less(max_diff, ANIB_THRESHOLD)


# Test concordance of ANIm with JSpecies output
def test_anim_concordance():
    """Test concordance of ANIm method with JSpecies output."""
    # Make/check output directory
    mode = "ANIm"
    outdirname = delete_and_remake_outdir(mode)
    nucmername = os.path.join(outdirname, 'nucmer_output') 
    os.makedirs(nucmername, exist_ok=True)

    # Get dataframes of JSpecies output
    anim_jspecies = parse_table(JSPECIES_OUTFILE, 'ANIm')

    # Identify our input files, and the total lengths of each organism seq
    infiles = pyani_files.get_fasta_files(INDIRNAME)
    org_lengths = pyani_files.get_sequence_lengths(infiles)

    # Test ANIm concordance:
    # Run pairwise NUCmer
    cmdlist = anim.generate_nucmer_commands(infiles, outdirname,
                                            pyani_config.NUCMER_DEFAULT)
    print('\n'.join(cmdlist))
    multiprocessing_run(cmdlist, verbose=False)
    # Process .delta files
    anim_data = anim.process_deltadir(nucmername, org_lengths)
    anim_pid = anim_data[1].sort_index(axis=0).sort_index(axis=1) * 100.

    print("ANIm data\n", anim_data)

    index, columns = anim_pid.index, anim_pid.columns
    diffmat = anim_pid.as_matrix() - anim_jspecies.as_matrix()
    anim_diff = pd.DataFrame(diffmat, index=index, columns=columns)

    # Write dataframes to file, for reference
    anim_pid.to_csv(os.path.join(outdirname,
                                'ANIm_pid.tab'),
                   sep='\t')
    anim_jspecies.to_csv(os.path.join(outdirname,
                                      'ANIm_jspecies.tab'),
                         sep='\t')
    anim_diff.to_csv(os.path.join(outdirname,
                                  'ANIm_diff.tab'),
                     sep='\t')
    print("ANIm concordance test output placed in %s" % outdirname)
    print("ANIm PID\n", anim_pid)
    print("ANIm JSpecies\n", anim_jspecies)
    print("ANIm diff\n", anim_diff)

    # We'd like the absolute difference reported to be < ANIB_THRESHOLD
    max_diff = anim_diff.abs().values.max()
    print("Maximum difference for ANIm: %e" % max_diff)
    assert_less(max_diff, ANIM_THRESHOLD)


# Test concordance of TETRA code with JSpecies output
@nottest
def test_tetra_concordance():
    """Test concordance of TETRA method with JSpecies output."""
    # Make/check output directory
    mode = "TETRA"
    outdirname = delete_and_remake_outdir(mode)

    # Get dataframes of JSpecies output
    tetra_jspecies = parse_table(JSPECIES_OUTFILE, 'Tetra')

    # Identify our input files, and the total lengths of each organism seq
    infiles = pyani_files.get_fasta_files(INDIRNAME)
    org_lengths = pyani_files.get_sequence_lengths(infiles)

    # Test TETRA concordance
    tetra_zscores = {}
    for filename in infiles:
        org = os.path.splitext(os.path.split(filename)[-1])[0]
        tetra_zscores[org] = tetra.calculate_tetra_zscore(filename)
    tetra_correlations = tetra.calculate_correlations(tetra_zscores)
    index, columns = tetra_correlations.index, tetra_correlations.columns
    tetra_diff = pd.DataFrame(tetra_correlations.as_matrix() -\
                              tetra_jspecies.as_matrix(),
                              index=index, columns=columns)

    # Write dataframes to file, for reference
    tetra_correlations.to_csv(os.path.join(outdirname,
                                           'tetra_correlations.tab'),
                              sep='\t')
    tetra_jspecies.to_csv(os.path.join(outdirname,
                                       'tetra_jspecies.tab'),
                          sep='\t')
    tetra_diff.to_csv(os.path.join(outdirname,
                                   'tetra_diff.tab'),
                      sep='\t')
    print("TETRA concordance test output placed in %s" % outdirname)
    print("TETRA correlations:\n", tetra_correlations)
    print("TETRA JSpecies:\n", tetra_jspecies)
    print("TETRA diff:\n", tetra_diff)

    # We'd like the absolute difference reported to be < TETRA_THRESHOLD
    max_diff = tetra_diff.abs().values.max()
    print("Maximum difference for TETRA: %e" % max_diff)
    assert_less(max_diff, TETRA_THRESHOLD)
    
