from scripts.mummer_alignedbases import AlignedBases
import pytest
from pathlib import Path


#Getting fixtures. Here, we pass the known AlignedBases values provided by dnadiff for both reference and query for all 14 tests
@pytest.fixture
def Test_1_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_1 sample (reference) """
    known = 4480707
    return known


@pytest.fixture
def Test_2_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_2 sample (reference) """
    known = 5264558
    return known


@pytest.fixture
def Test_3_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_3 sample (reference) """
    known = 444097
    return known


@pytest.fixture
def Test_4_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_4 sample (reference) """
    known = 6287002
    return known


@pytest.fixture
def Test_5_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_5 sample (reference) """
    known = 4832959
    return known


@pytest.fixture
def Test_6_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_6 sample (reference) """
    known = 39169
    return known


@pytest.fixture
def Test_7_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_7 sample (reference) """
    known = 29049
    return known


@pytest.fixture
def Test_8_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_8 sample (reference) """
    known = 11687
    return known


@pytest.fixture
def Test_9_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_9 sample (reference) """
    known = 24955
    return known


@pytest.fixture
def Test_10_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_10 sample (reference) """
    known = 18958
    return known


@pytest.fixture
def Test_11_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_11 sample (reference) """
    known = 7844783
    return known


@pytest.fixture
def Test_12_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_12 sample (reference) """
    known = 1574169
    return known


@pytest.fixture
def Test_13_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_13 sample (reference) """
    known = 1811679
    return known


@pytest.fixture
def Test_14_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_14 sample (reference) """
    known = 5612279
    return known






def test_AlignedBases_REF_test_1(Test_1_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_1 reference)
    The test should PASS."""

    our_AlignedBases_Test_1_REF = AlignedBases("issue_340_tests_AK/inputs/test_1/GCF_003369795.1_ASM336979v1_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_1/test_1.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_1/test_1.rdiff")

    assert our_AlignedBases_Test_1_REF == Test_1_AlignedBases_REF



def test_AlignedBases_REF_test_2(Test_2_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_2 reference)
    The test should PASS."""

    our_AlignedBases_Test_2_REF = AlignedBases("issue_340_tests_AK/inputs/test_2/GCF_002802945.1_ASM280294v1_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_2/test_2.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_2/test_2.rdiff")

    assert our_AlignedBases_Test_2_REF == Test_2_AlignedBases_REF


def test_AlignedBases_REF_test_3(Test_3_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_3 reference)
    The test should PASS."""

    our_AlignedBases_Test_3_REF = AlignedBases("issue_340_tests_AK/inputs/test_3/GCF_900095885.1_Daq1742_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_3/test_3.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_3/test_3.rdiff")

    assert our_AlignedBases_Test_3_REF == Test_3_AlignedBases_REF


def test_AlignedBases_REF_test_4(Test_4_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_4 reference)
    The test should PASS."""

    our_AlignedBases_Test_4_REF = AlignedBases("issue_340_tests_AK/inputs/test_4/genome1.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_4/test_4.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_4/test_4.rdiff")

    assert our_AlignedBases_Test_4_REF == Test_4_AlignedBases_REF


def test_AlignedBases_REF_test_5(Test_5_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_5 reference)
    The test should PASS."""

    our_AlignedBases_Test_5_REF = AlignedBases("issue_340_tests_AK/inputs/test_5/GCA_000011605.1_ASM1160v1_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_5/test_5.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_5/test_5.rdiff")

    assert our_AlignedBases_Test_5_REF == Test_5_AlignedBases_REF


def test_AlignedBases_REF_test_6(Test_6_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_6 reference)
    The test should PASS."""

    our_AlignedBases_Test_6_REF = AlignedBases("issue_340_tests_AK/inputs/test_6/MGV-GENOME-0264574.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_6/donovan_dnadiff.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_6/donovan_dnadiff.rdiff")

    assert our_AlignedBases_Test_6_REF == Test_6_AlignedBases_REF


def test_AlignedBases_REF_test_7(Test_7_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_7 reference)
    The test should PASS."""

    our_AlignedBases_Test_7_REF = AlignedBases("issue_340_tests_AK/inputs/test_7/GCA_031122125.1_ASM3112212v1_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_7/test_7.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_7/test_7.rdiff")

    assert our_AlignedBases_Test_7_REF == Test_7_AlignedBases_REF


def test_AlignedBases_REF_test_8(Test_8_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_8 reference)
    The test should PASS."""

    our_AlignedBases_Test_8_REF = AlignedBases("issue_340_tests_AK/inputs/test_8/GCA_000866645.1_ViralMultiSegProj15620_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_8/test_8.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_8/test_8.rdiff")

    assert our_AlignedBases_Test_8_REF == Test_8_AlignedBases_REF



def test_AlignedBases_REF_test_9(Test_9_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_9 reference)
    The test should PASS."""

    our_AlignedBases_Test_9_REF = AlignedBases("issue_340_tests_AK/inputs/test_9/GCA_032918235.1_ASM3291823v1_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_9/test_9.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_9/test_9.rdiff")

    assert our_AlignedBases_Test_9_REF == Test_9_AlignedBases_REF


def test_AlignedBases_REF_test_10(Test_10_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_10 reference)
    The test should PASS."""

    our_AlignedBases_Test_10_REF = AlignedBases("issue_340_tests_AK/inputs/test_10/GCA_034098425.1_ASM3409842v1_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_10/test_10.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_10/test_10.rdiff")

    assert our_AlignedBases_Test_10_REF == Test_10_AlignedBases_REF



def test_AlignedBases_REF_test_11(Test_11_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_11 reference)
    The test should PASS."""

    our_AlignedBases_Test_11_REF = AlignedBases("issue_340_tests_AK/inputs/test_11/GCF_016906225.1_ASM1690622v1_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_11/test_11.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_11/test_11.rdiff")

    assert our_AlignedBases_Test_11_REF == Test_11_AlignedBases_REF


def test_AlignedBases_REF_test_12(Test_12_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_12 reference)
    The test should PASS."""

    our_AlignedBases_Test_12_REF = AlignedBases("issue_340_tests_AK/inputs/test_12/GCF_000253235.1_ASM25323v1_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_12/test_12.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_12/test_12.rdiff")

    assert our_AlignedBases_Test_12_REF == Test_12_AlignedBases_REF


def test_AlignedBases_REF_test_13(Test_13_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_13 reference)
    The test should PASS."""

    our_AlignedBases_Test_13_REF = AlignedBases("issue_340_tests_AK/inputs/test_13/GCF_000719035.1_ASM71903v1_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_13/test_13.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_13/test_13.rdiff")

    assert our_AlignedBases_Test_13_REF == Test_13_AlignedBases_REF


def test_AlignedBases_REF_test_14(Test_14_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_14 reference)
    The test should PASS."""

    our_AlignedBases_Test_14_REF = AlignedBases("issue_340_tests_AK/inputs/test_14/GCF_002803175.1_ASM280317v1_genomic.fna", 
                                               "issue_340_tests_AK/outputs_dnadiff/test_14/test_14.mcoords",
                                               "issue_340_tests_AK/outputs_dnadiff/test_14/test_14.rdiff")

    assert our_AlignedBases_Test_14_REF == Test_14_AlignedBases_REF

