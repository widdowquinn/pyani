"""This script was used to replicate the value of AlignedBases provided by dnadiff using the same
approach as implemented by dnadiff.

Full code explanation is avaliable in issue_340_info.md.
"""

#Set Up
from Bio import SeqIO
import pandas as pd
from pathlib import Path


def AlignedBases(genome_seq, mcoords_file, rdiff_file):
    """Calculate and return AlignedBases value
    using the same approach as implemented by dnadiff.

    :param genome_seq: Path to a genome's sequence used in the analysis
    :param mcoords_file: Path to the mcoords file outputed by dnadiff
    :param rdiff_file: Path to the rdiff file outputed by dnadiff
    """

    #Step 1. Getting basic sequence information

    rnSeqs = 0
    rnBases = 0
    rnABases = 0

    records = list(SeqIO.parse(Path(genome_seq), "fasta"))
    refs = {record.id:len(record.seq) for record in records}

    for sequence_id, sequence_length in refs.items():
        rnSeqs += 1
        rnBases += sequence_length

    #Step 2. Retriving information from M-to-M alignments (mcoords)
    rqnAlignsM = 0
    rSumLenM = 0
    rqSumIdyM = 0
    rqSumLenM = 0
    rnASeqs = 0

    mcoords = pd.read_csv(Path(mcoords_file), sep='\t', names=[f"col_{_}" for _ in range(0,13)])

    for index, row in mcoords.iterrows():
        rqnAlignsM += 1
        rSumLenM += row['col_4']
        rqSumIdyM += (row['col_6']/100) * (row['col_4'] + row['col_5'])
        rqSumLenM += (row['col_4'] + row['col_5'])

    seen_ref_seq = []
    for index, row in mcoords.iterrows():
        if row['col_11'] not in seen_ref_seq:
            rnASeqs += 1
            rnABases += refs[row['col_11']]
            seen_ref_seq.append(row['col_11'])

    #Step 3: Retrive information from `.rdiff` files, and updating the AlignedBases values
    rdiff = pd.read_csv(Path(rdiff_file), sep='\t', names=[f"col_{_}" for _ in range(0,7)])

    rnTIns = 0
    rSumTIns = 0
    rnIns = 0
    rSumIns = 0

    for index, row in rdiff.iterrows():
        gap = row['col_4']
        ins = gap
        if row['col_1'] == 'GAP':
            if int(row['col_6']) > gap:
                ins = row['col_6']
            if int(row['col_4']) <=0 and int(row['col_5']) <=0 and int(row['col_6']) >0:
                rnTIns +=1
                rSumTIns += int(row['col_6'])
        if row['col_1'] != 'DUP' and gap >0:
            rnABases -= gap
        if int(ins) >0:
            rnIns +=1
            rSumIns += int(ins)

    return rnABases

