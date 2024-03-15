"""This script was used to investigate discrepancies between
average %ID reported by dnadiff and our current approach.
"""

# Set UP
from pathlib import Path


def parse_delta(filename):

    """Return weighedt %IDs calculated by
    diffrent approaches (approach 1, approach 2, approach 3).

    :param filename: Path to the deltafile

    Approach 1:
    Follows approach implemented by dnadiff.

    These calculations are repeated for each aligment:
    sequence ID =
        ((reference aln length + query aln length) -
        (2 * similarity error))
        / (reference aln length + query aln length) * 100


    weighted seq ID =
        (sequence ID rounded to two decimal places / 100) *
        ((refernce aln length + query aln length))

    Final calculations of %ID for all alignments
    %ID = sum of all weighted seq IDs / sum of all alignment lengths


    Approach 2:
    This is calulated as in approach 1, but no sequence ID values
    are rounded to two decimal places.

    Approach 3:
    Here, we skip the intermediate calculation of identity, as follow

    Calculations for each individual aligment
    weighted seq ID = (reference aln length + query aln length) -
                    (2 * similarity error)

    Final calculations of %ID for all alignments
    %ID = sum of all weighted seq IDs / sum of all alignment lengths
    """
    aligned_bases = []  # Hold a list for aligned bases for each sequence

    weighted_identical_AP_1 = []
    weighted_identical_AP_2 = []
    weighted_identical_AP_3 = []  # Hold a list for weighted identical bases

    for line in [_.strip().split() for _ in filename.open("r").readlines()]:
        if line[0] == "NUCMER":  # Skip headers
            continue
        # Lines with seven columns are alignment region headers:
        if len(line) == 7:

            # Calculate aligned bases for each sequence
            ref_aln_lengths = abs(int(line[1]) - int(line[0])) + 1
            qry_aln_lengths = abs(int(line[3]) - int(line[2])) + 1
            aligned_bases.append(ref_aln_lengths)
            aligned_bases.append(qry_aln_lengths)

            # Calculating Sequence IDs
            sequence_perecentage_id = (
                ((ref_aln_lengths + qry_aln_lengths) - (2 * int(line[4])))
                / (ref_aln_lengths + qry_aln_lengths)
                * 100
            )

            # Calculating Weighted IDs
            weighted_identical_AP_1.append(
                (round(sequence_perecentage_id, 2) / 100)
                * ((ref_aln_lengths + qry_aln_lengths))
            )
            weighted_identical_AP_2.append(
                (sequence_perecentage_id / 100) * ((ref_aln_lengths + qry_aln_lengths))
            )
            weighted_identical_AP_3.append(
                (ref_aln_lengths + qry_aln_lengths) - (2 * int(line[4]))
            )

    # Calculating overall %ID
    ID_AP_1 = round(sum(weighted_identical_AP_1) / sum(aligned_bases) * 100, 2)

    ID_AP_2 = sum(weighted_identical_AP_2) / sum(aligned_bases) * 100

    ID_AP_3 = sum(weighted_identical_AP_3) / sum(aligned_bases) * 100

    return ID_AP_1, ID_AP_2, ID_AP_3


print(parse_delta(Path("../data/donovan_test/output/dnadiff/output.1delta")))
print(parse_delta(Path("../data/donovan_AF_bug/output/dnadiff/output.1delta")))
print(parse_delta(Path("../data/streptomyces_genomes/output/dnadiff/output.1delta")))
