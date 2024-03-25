from scripts.parsing_delta import parse_delta
from scripts.simple_parsing_LP import parse_coords
import pytest
from pathlib import Path




"""The below tests were written to check if the returned value for preasumed total number of aligned bases
returned by parse_coords function (mcoords) are the same as the ones returned by parse_coords.
"""

#Test 1
def test_1_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 1 - Reference)
    The test should PASS."""

    parse_coords_test_1_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_1/test_1.mcoords"))[0]
    parse_delta_test_1_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_1/test_1.mdelta"))[0]

    assert parse_coords_test_1_REF == parse_delta_test_1_REF

def test_1_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 1 - Query)
    The test should PASS."""

    parse_coords_test_1_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_1/test_1.mcoords"))[1]
    parse_delta_test_1_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_1/test_1.mdelta"))[1]

    assert parse_coords_test_1_QRY == parse_delta_test_1_QRY

#Test 2
def test_2_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 2 - Reference)
    The test should PASS."""

    parse_coords_test_2_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_2/test_2.mcoords"))[0]
    parse_delta_test_2_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_2/test_2.mdelta"))[0]

    assert parse_coords_test_2_REF == parse_delta_test_2_REF

def test_2_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 2 - Query)
    The test should PASS."""

    parse_coords_test_2_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_2/test_2.mcoords"))[1]
    parse_delta_test_2_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_2/test_2.mdelta"))[1]

    assert parse_coords_test_2_QRY == parse_delta_test_2_QRY


#Test 3
def test_3_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 3 - Reference)
    The test should PASS."""

    parse_coords_test_3_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_3/test_3.mcoords"))[0]
    parse_delta_test_3_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_3/test_3.mdelta"))[0]

    assert parse_coords_test_3_REF == parse_delta_test_3_REF

def test_3_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 3 - Query)
    The test should PASS."""

    parse_coords_test_3_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_3/test_3.mcoords"))[1]
    parse_delta_test_3_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_3/test_3.mdelta"))[1]

    assert parse_coords_test_3_QRY == parse_delta_test_3_QRY


#Test 4
def test_4_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 4 - Reference)
    The test should PASS."""

    parse_coords_test_4_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_4/test_4.mcoords"))[0]
    parse_delta_test_4_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_4/test_4.mdelta"))[0]

    assert parse_coords_test_4_REF == parse_delta_test_4_REF

def test_4_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 4 - Query)
    The test should PASS."""

    parse_coords_test_4_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_4/test_4.mcoords"))[1]
    parse_delta_test_4_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_4/test_4.mdelta"))[1]

    assert parse_coords_test_4_QRY == parse_delta_test_4_QRY


#Test 5
def test_5_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 5 - Reference)
    The test should PASS."""

    parse_coords_test_5_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_5/test_5.mcoords"))[0]
    parse_delta_test_5_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_5/test_5.mdelta"))[0]

    assert parse_coords_test_5_REF == parse_delta_test_5_REF

def test_5_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 5 - Query)
    The test should PASS."""

    parse_coords_test_5_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_5/test_5.mcoords"))[1]
    parse_delta_test_5_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_5/test_5.mdelta"))[1]

    assert parse_coords_test_5_QRY == parse_delta_test_5_QRY


    #Test 6
def test_6_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 6 - Reference)
    The test should PASS."""

    parse_coords_test_6_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_6/donovan_dnadiff.mcoords"))[0]
    parse_delta_test_6_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_6/donovan_dnadiff.mdelta"))[0]

    assert parse_coords_test_6_REF == parse_delta_test_6_REF

def test_6_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 6 - Query)
    The test should PASS."""

    parse_coords_test_6_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_6/donovan_dnadiff.mcoords"))[1]
    parse_delta_test_6_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_6/donovan_dnadiff.mdelta"))[1]

    assert parse_coords_test_6_QRY == parse_delta_test_6_QRY

    #Test 7
def test_7_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 7 - Reference)
    The test should PASS."""

    parse_coords_test_7_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_7/test_7.mcoords"))[0]
    parse_delta_test_7_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_7/test_7.mdelta"))[0]

    assert parse_coords_test_7_REF == parse_delta_test_7_REF

def test_7_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 7 - Query)
    The test should PASS."""

    parse_coords_test_7_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_7/test_7.mcoords"))[1]
    parse_delta_test_7_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_7/test_7.mdelta"))[1]

    assert parse_coords_test_7_QRY == parse_delta_test_7_QRY


    #Test 8
def test_8_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 8 - Reference)
    The test should PASS."""

    parse_coords_test_8_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_8/test_8.mcoords"))[0]
    parse_delta_test_8_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_8/test_8.mdelta"))[0]

    assert parse_coords_test_8_REF == parse_delta_test_8_REF

def test_8_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 8 - Query)
    The test should PASS."""

    parse_coords_test_8_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_8/test_8.mcoords"))[1]
    parse_delta_test_8_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_8/test_8.mdelta"))[1]

    assert parse_coords_test_8_QRY == parse_delta_test_8_QRY

    #Test 9
def test_9_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 9 - Reference)
    The test should PASS."""

    parse_coords_test_9_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_9/test_9.mcoords"))[0]
    parse_delta_test_9_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_9/test_9.mdelta"))[0]

    assert parse_coords_test_9_REF == parse_delta_test_9_REF

def test_9_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 9 - Query)
    The test should PASS."""

    parse_coords_test_9_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_9/test_9.mcoords"))[1]
    parse_delta_test_9_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_9/test_9.mdelta"))[1]

    assert parse_coords_test_9_QRY == parse_delta_test_9_QRY


    #Test 10
def test_10_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 10 - Reference)
    The test should PASS."""

    parse_coords_test_10_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_10/test_10.mcoords"))[0]
    parse_delta_test_10_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_10/test_10.mdelta"))[0]

    assert parse_coords_test_10_REF == parse_delta_test_10_REF

def test_10_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 10 - Query)
    The test should PASS."""

    parse_coords_test_10_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_10/test_10.mcoords"))[1]
    parse_delta_test_10_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_10/test_10.mdelta"))[1]

    assert parse_coords_test_10_QRY == parse_delta_test_10_QRY


    #Test 11
def test_11_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 11 - Reference)
    The test should PASS."""

    parse_coords_test_11_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_11/test_11.mcoords"))[0]
    parse_delta_test_11_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_11/test_11.mdelta"))[0]

    assert parse_coords_test_11_REF == parse_delta_test_11_REF

def test_11_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 11 - Query)
    The test should PASS."""

    parse_coords_test_11_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_11/test_11.mcoords"))[1]
    parse_delta_test_11_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_11/test_11.mdelta"))[1]

    assert parse_coords_test_11_QRY == parse_delta_test_11_QRY


    #Test 12
def test_12_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 12 - Reference)
    The test should PASS."""

    parse_coords_test_12_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_12/test_12.mcoords"))[0]
    parse_delta_test_12_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_12/test_12.mdelta"))[0]

    assert parse_coords_test_12_REF == parse_delta_test_12_REF

def test_12_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 12 - Query)
    The test should PASS."""

    parse_coords_test_12_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_12/test_12.mcoords"))[1]
    parse_delta_test_12_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_12/test_12.mdelta"))[1]

    assert parse_coords_test_12_QRY == parse_delta_test_12_QRY

    #Test 13
def test_13_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 13 - Reference)
    The test should PASS."""

    parse_coords_test_13_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_13/test_13.mcoords"))[0]
    parse_delta_test_13_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_13/test_13.mdelta"))[0]

    assert parse_coords_test_13_REF == parse_delta_test_13_REF

def test_13_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 13 - Query)
    The test should PASS."""

    parse_coords_test_13_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_13/test_13.mcoords"))[1]
    parse_delta_test_13_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_13/test_13.mdelta"))[1]

    assert parse_coords_test_13_QRY == parse_delta_test_13_QRY


    #Test 14
def test_14_parse_REF():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 14 - Reference)
    The test should PASS."""

    parse_coords_test_14_REF = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_14/test_14.mcoords"))[0]
    parse_delta_test_14_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_14/test_14.mdelta"))[0]

    assert parse_coords_test_14_REF == parse_delta_test_14_REF

def test_14_parse_QRY():
    """Comparisions of AlignedBases from mdelta and mcoords. (Test 14 - Query)
    The test should PASS."""

    parse_coords_test_14_QRY = parse_coords(Path("issue_340_tests_AK/outputs_dnadiff/test_14/test_14.mcoords"))[1]
    parse_delta_test_14_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_14/test_14.mdelta"))[1]

    assert parse_coords_test_14_QRY == parse_delta_test_14_QRY