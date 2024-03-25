from scripts.parsing_delta import parse_delta
import pytest
from pathlib import Path



#Getting fixtures. Here, we pass the known TotalLength values provided by dnadiff for both reference and query for all 14 tests
@pytest.fixture
def Test_1_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_1 sample (reference) """
    known = 4508050
    return known

@pytest.fixture
def Test_1_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_1 sample (query) """
    known = 4509104
    return known
@pytest.fixture
def Test_2_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_2 sample (reference) """
    known = 5284424
    return known

@pytest.fixture
def Test_2_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_2 sample (query) """
    known = 5264616
    return known

@pytest.fixture
def Test_3_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_3 sample (reference) """
    known = 460558
    return known

@pytest.fixture
def Test_3_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_3 sample (query) """
    known = 460939
    return known

@pytest.fixture
def Test_4_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_4 sample (reference) """
    known = 6375058
    return known

@pytest.fixture
def Test_4_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_4 sample (query) """
    known = 6374916
    return known

@pytest.fixture
def Test_5_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_5 sample (reference) """
    known = 4908215
    return known

@pytest.fixture
def Test_5_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_5 sample (query) """
    known = 4908459
    return known

@pytest.fixture
def Test_6_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_6 sample (reference) """
    known = 59174
    return known

@pytest.fixture
def Test_6_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_6 sample (query) """
    known = 59187
    return known

@pytest.fixture
def Test_7_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_7 sample (reference) """
    known = 29049
    return known

@pytest.fixture
def Test_7_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_7 sample (query) """
    known = 29071
    return known

@pytest.fixture
def Test_8_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_8 sample (reference) """
    known = 11687
    return known

@pytest.fixture
def Test_8_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_8 sample (query) """
    known = 11687
    return known

@pytest.fixture
def Test_9_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_9 sample (reference) """
    known = 24955
    return known

@pytest.fixture
def Test_9_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_9 sample (query) """
    known = 24955
    return known

@pytest.fixture
def Test_10_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_10 sample (reference) """
    known = 18958
    return known

@pytest.fixture
def Test_10_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_10 sample (query) """
    known = 18958
    return known

@pytest.fixture
def Test_11_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_11 sample (reference) """
    known = 8483896
    return known

@pytest.fixture
def Test_11_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_11 sample (query) """
    known = 8484036
    return known

@pytest.fixture
def Test_12_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_12 sample (reference) """
    known = 1626091
    return known

@pytest.fixture
def Test_12_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_12 sample (query) """
    known = 1628081
    return known

@pytest.fixture
def Test_13_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_13 sample (reference) """
    known = 1827407
    return known

@pytest.fixture
def Test_13_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_13 sample (query) """
    known = 1825865
    return known

@pytest.fixture
def Test_14_TotalLength_REF():
    """Known TotalLength provided by dnadiff for test_14 sample (reference) """
    known = 5688413
    return known

@pytest.fixture
def Test_14_TotalLength_QRY():
    """Known TotalLength provided by dnadiff for test_14 sample (query) """
    known = 5681712
    return known

#Comparing TotalLengths values from dnadiff and the values obtained from our code

def test_TotalLength_REF_test_1(Test_1_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_1 reference)
    The test should PASS."""

    TotalLength_Test_1_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_1/test_1.mdelta"))[2]

    assert TotalLength_Test_1_REF == Test_1_TotalLength_REF


def test_TotalLength_QRY_test_1(Test_1_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_2 query)
    The test should PASS."""

    TotalLength_Test_1_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_1/test_1.mdelta"))[3]

    assert TotalLength_Test_1_QRY == Test_1_TotalLength_QRY

def test_TotalLength_REF_test_2(Test_2_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_2 reference)
    The test should PASS."""

    TotalLength_Test_2_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_2/test_2.mdelta"))[2]

    assert TotalLength_Test_2_REF == Test_2_TotalLength_REF


def test_TotalLength_QRY_test_2(Test_2_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_2 query)
    The test should PASS."""

    TotalLength_Test_2_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_2/test_2.mdelta"))[3]

    assert TotalLength_Test_2_QRY == Test_2_TotalLength_QRY


def test_TotalLength_REF_test_3(Test_3_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_3 reference)
    The test should PASS."""

    TotalLength_Test_3_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_3/test_3.mdelta"))[2]

    assert TotalLength_Test_3_REF == Test_3_TotalLength_REF


def test_TotalLength_QRY_test_3(Test_3_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_3 query)
    The test should PASS."""

    TotalLength_Test_3_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_3/test_3.mdelta"))[3]

    assert TotalLength_Test_3_QRY == Test_3_TotalLength_QRY


def test_TotalLength_REF_test_4(Test_4_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_4 reference)
    The test should PASS."""

    TotalLength_Test_4_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_4/test_4.mdelta"))[2]

    assert TotalLength_Test_4_REF == Test_4_TotalLength_REF


def test_TotalLength_QRY_test_4(Test_4_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_4 query)
    The test should PASS."""

    TotalLength_Test_4_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_4/test_4.mdelta"))[3]

    assert TotalLength_Test_4_QRY == Test_4_TotalLength_QRY


def test_TotalLength_REF_test_5(Test_5_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_5 reference)
    The test should PASS."""

    TotalLength_Test_5_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_5/test_5.mdelta"))[2]

    assert TotalLength_Test_5_REF == Test_5_TotalLength_REF


def test_TotalLength_QRY_test_5(Test_5_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_5 query)
    The test should PASS."""

    TotalLength_Test_5_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_5/test_5.mdelta"))[3]

    assert TotalLength_Test_5_QRY == Test_5_TotalLength_QRY


def test_TotalLength_REF_test_6(Test_6_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_6 reference)
    The test should PASS."""

    TotalLength_Test_6_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_6/donovan_dnadiff.mdelta"))[2]

    assert TotalLength_Test_6_REF == Test_6_TotalLength_REF


def test_TotalLength_QRY_test_6(Test_6_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_6 query)
    The test should PASS."""

    TotalLength_Test_6_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_6/donovan_dnadiff.mdelta"))[3]

    assert TotalLength_Test_6_QRY == Test_6_TotalLength_QRY


def test_TotalLength_REF_test_7(Test_7_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_7 reference)
    The test should PASS."""

    TotalLength_Test_7_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_7/test_7.mdelta"))[2]

    assert TotalLength_Test_7_REF == Test_7_TotalLength_REF


def test_TotalLength_QRY_test_7(Test_7_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_7 query)
    The test should PASS."""

    TotalLength_Test_7_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_7/test_7.mdelta"))[3]

    assert TotalLength_Test_7_QRY == Test_7_TotalLength_QRY


def test_TotalLength_REF_test_8(Test_8_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_8 reference)
    The test should PASS."""

    TotalLength_Test_8_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_8/test_8.mdelta"))[2]

    assert TotalLength_Test_8_REF == Test_8_TotalLength_REF


def test_TotalLength_QRY_test_8(Test_8_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_8 query)
    The test should PASS."""

    TotalLength_Test_8_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_8/test_8.mdelta"))[3]

    assert TotalLength_Test_8_QRY == Test_8_TotalLength_QRY


def test_TotalLength_REF_test_9(Test_9_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_9 reference)
    The test should PASS."""

    TotalLength_Test_9_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_9/test_9.mdelta"))[2]

    assert TotalLength_Test_9_REF == Test_9_TotalLength_REF


def test_TotalLength_QRY_test_9(Test_9_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_9 query)
    The test should PASS."""

    TotalLength_Test_9_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_9/test_9.mdelta"))[3]

    assert TotalLength_Test_9_QRY == Test_9_TotalLength_QRY


def test_TotalLength_REF_test_10(Test_10_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_10 reference)
    The test should PASS."""

    TotalLength_Test_10_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_10/test_10.mdelta"))[2]

    assert TotalLength_Test_10_REF == Test_10_TotalLength_REF


def test_TotalLength_QRY_test_10(Test_10_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_10 query)
    The test should PASS."""

    TotalLength_Test_10_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_10/test_10.mdelta"))[3]

    assert TotalLength_Test_10_QRY == Test_10_TotalLength_QRY


def test_TotalLength_REF_test_11(Test_11_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_11 reference)
    The test should PASS."""

    TotalLength_Test_11_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_11/test_11.mdelta"))[2]

    assert TotalLength_Test_11_REF == Test_11_TotalLength_REF


def test_TotalLength_QRY_test_11(Test_11_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_11 query)
    The test should PASS."""

    TotalLength_Test_11_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_11/test_11.mdelta"))[3]

    assert TotalLength_Test_11_QRY == Test_11_TotalLength_QRY


def test_TotalLength_REF_test_12(Test_12_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_12 reference)
    The test should PASS."""

    TotalLength_Test_12_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_12/test_12.mdelta"))[2]

    assert TotalLength_Test_12_REF == Test_12_TotalLength_REF


def test_TotalLength_QRY_test_12(Test_12_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_12 query)
    The test should PASS."""

    TotalLength_Test_12_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_12/test_12.mdelta"))[3]

    assert TotalLength_Test_12_QRY == Test_12_TotalLength_QRY


def test_TotalLength_REF_test_13(Test_13_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_13 reference)
    The test should PASS."""

    TotalLength_Test_13_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_13/test_13.mdelta"))[2]

    assert TotalLength_Test_13_REF == Test_13_TotalLength_REF


def test_TotalLength_QRY_test_13(Test_13_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_13 query)
    The test should PASS."""

    TotalLength_Test_13_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_13/test_13.mdelta"))[3]

    assert TotalLength_Test_13_QRY == Test_13_TotalLength_QRY


def test_TotalLength_REF_test_14(Test_14_TotalLength_REF):
    """TotalLength values from mdelta and dnadiff values comparision. (test_14 reference)
    The test should PASS."""

    TotalLength_Test_14_REF = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_14/test_14.mdelta"))[2]

    assert TotalLength_Test_14_REF == Test_14_TotalLength_REF


def test_TotalLength_QRY_test_14(Test_14_TotalLength_QRY):
    """TotalLength values from mdelta and dnadiff values comparision. (test_14 query)
    The test should PASS."""

    TotalLength_Test_14_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_14/test_14.mdelta"))[3]

    assert TotalLength_Test_14_QRY == Test_14_TotalLength_QRY

#Running additional tests on false data. All should fail.
def test_TotalLength_test_15(Test_1_TotalLength_REF):
    """Passing parse_delta function on a file that does not exists.
    The test should FAIL."""

    TotalLength_Test_15_QRY = parse_delta(Path("this_path_does_not_exists"))[3]

    assert TotalLength_Test_15_QRY == Test_1_TotalLength_REF

def test_TotalLength_QRY_test_16(Test_1_TotalLength_QRY):
    """Comparing the values, but comparing two diffrent types (eg. string to int)
    The test should FAIL."""

    TotalLength_Test_1_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_1/test_1.mdelta"))[3]

    assert TotalLength_Test_1_QRY == str(Test_1_TotalLength_QRY)

def test_TotalLength_QRY_test_17(Test_2_TotalLength_QRY):
    """Comparing wrong values. Here, TotalLength from test 2 (query) from dnadiff to our TotalLength obtained from dataset test 1 (query). 
    The test should FAIL."""

    TotalLength_Test_1_QRY = parse_delta(Path("issue_340_tests_AK/outputs_dnadiff/test_1/test_1.mdelta"))[3]

    assert TotalLength_Test_1_QRY == Test_2_TotalLength_QRY