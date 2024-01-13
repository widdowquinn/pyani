from scripts.parsing_delta import parse_delta
import pytest
from pathlib import Path


#Getting fixtures. Here, we pass the known AlignedBases values provided by dnadiff for both reference and query for all 14 tests
@pytest.fixture
def Test_1_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_1 sample (reference) """
    known = 4480707
    return known

@pytest.fixture
def Test_1_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_1 sample (query) """
    known = 4483285
    return known

@pytest.fixture
def Test_2_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_2 sample (reference) """
    known = 5264558
    return known

@pytest.fixture
def Test_2_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_2 sample (query) """
    known = 5251795
    return known

@pytest.fixture
def Test_3_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_3 sample (reference) """
    known = 444097
    return known

@pytest.fixture
def Test_3_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_3 sample (query) """
    known = 415667
    return known

@pytest.fixture
def Test_4_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_4 sample (reference) """
    known = 6287002
    return known

@pytest.fixture
def Test_4_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_4 sample (query) """
    known = 6274559
    return known

@pytest.fixture
def Test_5_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_5 sample (reference) """
    known = 4832959
    return known

@pytest.fixture
def Test_5_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_5 sample (query) """
    known = 4829346
    return known

@pytest.fixture
def Test_6_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_6 sample (reference) """
    known = 39169
    return known

@pytest.fixture
def Test_6_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_6 sample (query) """
    known = 39176
    return known

@pytest.fixture
def Test_7_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_7 sample (reference) """
    known = 29049
    return known

@pytest.fixture
def Test_7_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_7 sample (query) """
    known = 29071
    return known

@pytest.fixture
def Test_8_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_8 sample (reference) """
    known = 11687
    return known

@pytest.fixture
def Test_8_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_8 sample (query) """
    known = 11687
    return known

@pytest.fixture
def Test_9_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_9 sample (reference) """
    known = 24955
    return known

@pytest.fixture
def Test_9_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_9 sample (query) """
    known = 24955
    return known

@pytest.fixture
def Test_10_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_10 sample (reference) """
    known = 18958
    return known

@pytest.fixture
def Test_10_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_10 sample (query) """
    known = 18958
    return known

@pytest.fixture
def Test_11_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_11 sample (reference) """
    known = 7844783
    return known

@pytest.fixture
def Test_11_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_11 sample (query) """
    known = 7838472
    return known

@pytest.fixture
def Test_12_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_12 sample (reference) """
    known = 1574169
    return known

@pytest.fixture
def Test_12_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_12 sample (query) """
    known = 1583270
    return known

@pytest.fixture
def Test_13_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_13 sample (reference) """
    known = 1811679
    return known

@pytest.fixture
def Test_13_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_13 sample (query) """
    known = 1811473
    return known

@pytest.fixture
def Test_14_AlignedBases_REF():
    """Known TotalLength provided by dnadiff for test_14 sample (reference) """
    known = 5612279
    return known

@pytest.fixture
def Test_14_AlignedBases_QRY():
    """Known TotalLength provided by dnadiff for test_14 sample (query) """
    known = 5669103
    return known




def test_AlignedBases_REF_test_1(Test_1_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_1 reference)
    The test should PASS."""

    our_AlignedBases_Test_1_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_1/test_1.mdelta"))[0]

    assert our_AlignedBases_Test_1_REF == Test_1_AlignedBases_REF

def test_AlignedBases_QRY_test_1(Test_1_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_1 query)
    The test should PASS."""

    our_AlignedBases_Test_1_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_1/test_1.mdelta"))[1]

    assert our_AlignedBases_Test_1_QRY == Test_1_AlignedBases_QRY


def test_AlignedBases_REF_test_2(Test_2_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_2 reference)
    The test should PASS."""

    our_AlignedBases_Test_2_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_2/test_2.mdelta"))[0]

    assert our_AlignedBases_Test_2_REF == Test_2_AlignedBases_REF

def test_AlignedBases_QRY_test_2(Test_2_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_2 query)
    The test should PASS."""

    our_AlignedBases_Test_2_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_2/test_2.mdelta"))[1]

    assert our_AlignedBases_Test_2_QRY == Test_2_AlignedBases_QRY

def test_AlignedBases_REF_test_3(Test_3_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_3 reference)
    The test should PASS."""

    our_AlignedBases_Test_3_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_3/test_3.mdelta"))[0]

    assert our_AlignedBases_Test_3_REF == Test_3_AlignedBases_REF

def test_AlignedBases_QRY_test_3(Test_3_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_3 query)
    The test should PASS."""

    our_AlignedBases_Test_3_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_3/test_3.mdelta"))[1]

    assert our_AlignedBases_Test_3_QRY == Test_3_AlignedBases_QRY


def test_AlignedBases_REF_test_4(Test_4_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_4 reference)
    The test should PASS."""

    our_AlignedBases_Test_4_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_4/test_4.mdelta"))[0]

    assert our_AlignedBases_Test_4_REF == Test_4_AlignedBases_REF

def test_AlignedBases_QRY_test_4(Test_4_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_4 query)
    The test should PASS."""

    our_AlignedBases_Test_4_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_4/test_4.mdelta"))[1]

    assert our_AlignedBases_Test_4_QRY == Test_4_AlignedBases_QRY

def test_AlignedBases_REF_test_5(Test_5_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_5 reference)
    The test should PASS."""

    our_AlignedBases_Test_5_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_5/test_5.mdelta"))[0]

    assert our_AlignedBases_Test_5_REF == Test_5_AlignedBases_REF

def test_AlignedBases_QRY_test_5(Test_5_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_5 query)
    The test should PASS."""

    our_AlignedBases_Test_5_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_5/test_5.mdelta"))[1]

    assert our_AlignedBases_Test_5_QRY == Test_5_AlignedBases_QRY

def test_AlignedBases_REF_test_6(Test_6_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_6 reference)
    The test should PASS."""

    our_AlignedBases_Test_6_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_6/donovan_dnadiff.mdelta"))[0]

    assert our_AlignedBases_Test_6_REF == Test_6_AlignedBases_REF

def test_AlignedBases_QRY_test_6(Test_6_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_6 query)
    The test should PASS."""

    our_AlignedBases_Test_6_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_6/donovan_dnadiff.mdelta"))[1]

    assert our_AlignedBases_Test_6_QRY == Test_6_AlignedBases_QRY

def test_AlignedBases_REF_test_7(Test_7_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_7 reference)
    The test should PASS."""

    our_AlignedBases_Test_7_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_7/test_7.mdelta"))[0]

    assert our_AlignedBases_Test_7_REF == Test_7_AlignedBases_REF

def test_AlignedBases_QRY_test_7(Test_7_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_7 query)
    The test should PASS."""

    our_AlignedBases_Test_7_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_7/test_7.mdelta"))[1]

    assert our_AlignedBases_Test_7_QRY == Test_7_AlignedBases_QRY

def test_AlignedBases_REF_test_8(Test_8_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_8 reference)
    The test should PASS."""

    our_AlignedBases_Test_8_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_8/test_8.mdelta"))[0]

    assert our_AlignedBases_Test_8_REF == Test_8_AlignedBases_REF

def test_AlignedBases_QRY_test_8(Test_8_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_8 query)
    The test should PASS."""

    our_AlignedBases_Test_8_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_8/test_8.mdelta"))[1]

    assert our_AlignedBases_Test_8_QRY == Test_8_AlignedBases_QRY


def test_AlignedBases_REF_test_9(Test_9_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_9 reference)
    The test should PASS."""

    our_AlignedBases_Test_9_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_9/test_9.mdelta"))[0]

    assert our_AlignedBases_Test_9_REF == Test_9_AlignedBases_REF

def test_AlignedBases_QRY_test_9(Test_9_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_9 query)
    The test should PASS."""

    our_AlignedBases_Test_9_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_9/test_9.mdelta"))[1]

    assert our_AlignedBases_Test_9_QRY == Test_9_AlignedBases_QRY

def test_AlignedBases_REF_test_10(Test_10_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_10 reference)
    The test should PASS."""

    our_AlignedBases_Test_10_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_10/test_10.mdelta"))[0]

    assert our_AlignedBases_Test_10_REF == Test_10_AlignedBases_REF

def test_AlignedBases_QRY_test_10(Test_10_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_10 query)
    The test should PASS."""

    our_AlignedBases_Test_10_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_10/test_10.mdelta"))[1]

    assert our_AlignedBases_Test_10_QRY == Test_10_AlignedBases_QRY

def test_AlignedBases_REF_test_11(Test_11_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_11 reference)
    The test should PASS."""

    our_AlignedBases_Test_11_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_11/test_11.mdelta"))[0]

    assert our_AlignedBases_Test_11_REF == Test_11_AlignedBases_REF

def test_AlignedBases_QRY_test_11(Test_11_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_11 query)
    The test should PASS."""

    our_AlignedBases_Test_11_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_11/test_11.mdelta"))[1]

    assert our_AlignedBases_Test_11_QRY == Test_11_AlignedBases_QRY

def test_AlignedBases_REF_test_12(Test_12_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_12 reference)
    The test should PASS."""

    our_AlignedBases_Test_12_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_12/test_12.mdelta"))[0]

    assert our_AlignedBases_Test_12_REF == Test_12_AlignedBases_REF

def test_AlignedBases_QRY_test_12(Test_12_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_12 query)
    The test should PASS."""

    our_AlignedBases_Test_12_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_12/test_12.mdelta"))[1]

    assert our_AlignedBases_Test_12_QRY == Test_12_AlignedBases_QRY

def test_AlignedBases_REF_test_13(Test_13_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_13 reference)
    The test should PASS."""

    our_AlignedBases_Test_13_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_13/test_13.mdelta"))[0]

    assert our_AlignedBases_Test_13_REF == Test_13_AlignedBases_REF

def test_AlignedBases_QRY_test_13(Test_13_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_13 query)
    The test should PASS."""

    our_AlignedBases_Test_13_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_13/test_13.mdelta"))[1]

    assert our_AlignedBases_Test_13_QRY == Test_13_AlignedBases_QRY

def test_AlignedBases_REF_test_14(Test_14_AlignedBases_REF):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_14 reference)
    The test should PASS."""

    our_AlignedBases_Test_14_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_14/test_14.mdelta"))[0]

    assert our_AlignedBases_Test_14_REF == Test_14_AlignedBases_REF

def test_AlignedBases_QRY_test_14(Test_14_AlignedBases_QRY):
    """AlignedBases values from mdelta and dnadiff values comparision. (test_14 query)
    The test should PASS."""

    our_AlignedBases_Test_14_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_14/test_14.mdelta"))[1]

    assert our_AlignedBases_Test_14_QRY == Test_14_AlignedBases_QRY
