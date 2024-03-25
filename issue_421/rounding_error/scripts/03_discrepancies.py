"""This script was used to investigate, discrepancies
discrepancies in our and dnadiff value for %ID.

Comparing of two Streptomyces genomes with dnadiff
reports %ID of 85.35%, where we report 85.37%,
regardless of approach used to calulate these values.


HOW DOES DNADIFF CALCULATES %ID?
We know that dnadiff to calulate the average %ID,
uses .1coords file. The calulations are as follow:

Assigns Combined alignment lengths and weighted
alignment length sum to 0 (for both 1-to-1
and M-to-M aligments).
    my ($rqSumLen1, $rqSumLenM) = (0,0);
    my ($rqSumIdy1, $rqSumIdyM) = (0,0);

The weighted aligment lengths are calulated, by
deviding the average %ID of each sequence by a 100,
and multiplying it by the sum of the length of both
the reference and query.
    $rqSumIdy1 += ($A[6] / 100.0) * ($A[4] + $A[5]);

The length of aligments are calulated as the sum of
the length of both the reference and query:
    $rqSumLen1 += ($A[4] + $A[5]);

The Averge Identity is calulates by deviding the sum
of all weighted alignment lengths, by the sum of all
combined alignment lengths, and multiplying that value
by a 100.
        ($rqSumLen1 ? $rqSumIdy1 / $rqSumLen1 * 100.0 : 0)
"""
# Set Up
import pandas as pd
from pathlib import Path


# Calulating values from coord files
coords1 = pd.read_csv(
    "../data/streptomyces_genomes/output/dnadiff/output.1coords", sep="\t", header=None
)
column_names = [f"col_{i}" for i in range(0, len(coords1.columns))]
coords1.columns = column_names

reqSumIdy1 = 0
rqSumLen1 = 0


for index, row in coords1.iterrows():
    reqSumIdy1 += (row["col_6"] / 100) * (row["col_4"] + row["col_5"])

    rqSumLen1 += row["col_4"] + row["col_5"]

AverageID = reqSumIdy1 / rqSumLen1 * 100

print(f"Average %ID from 1coords file is: {AverageID}")


# Knowing, that we can replicate the same Average %ID with
# coord files and wich values are exaclty used to calculate
# these values, we can identify which values do not match.

# Firslty, let's check if the individual aligment %IDs
# match the ones we have calculated.


check = []

for line in [
    _.strip().split()
    for _ in Path("../data/streptomyces_genomes/output/dnadiff/output.1delta")
    .open("r")
    .readlines()
]:
    if line[0] == "NUCMER":  # Skip headers
        continue
        # Lines with seven columns are alignment region headers:
    if len(line) == 7:

        # Calculate aligned bases for each sequence
        ref_aln_lengths = abs(int(line[1]) - int(line[0])) + 1
        qry_aln_lengths = abs(int(line[3]) - int(line[2])) + 1

        # Calculating Sequence IDs
        sequence_perecentage_id = (
            ((ref_aln_lengths + qry_aln_lengths) - (2 * int(line[4])))
            / (ref_aln_lengths + qry_aln_lengths)
            * 100
        )

        print(sequence_perecentage_id)

        check.append(
            [ref_aln_lengths, qry_aln_lengths, round(sequence_perecentage_id, 2)]
        )


for index, row in coords1.iterrows():
    if [row["col_4"], row["col_5"], row["col_6"]] not in check:
        print([row["col_4"], row["col_5"], row["col_6"]])
