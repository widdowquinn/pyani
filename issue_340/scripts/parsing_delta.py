"""This script was written in regards to pyANI issue #340.

The main aim of this script is to replicate the value of AlignedBases given by dnadiff by parsing .mdelta file. 
"""

#Set Up
from pathlib import Path
from collections import defaultdict
import intervaltree


def parse_delta(infname):
    """Parse delta files.

    :param infname: Path to delta file
    """

    current_ref, current_qry = None, None

    regions_ref = defaultdict(list) #Hold a dictionary for refence regions
    regions_qry = defaultdict(list) #Hold a dictionary for query regions

    for line in [_.strip().split() for _ in Path(infname).open("r").readlines()]:

            if line[0] == "NUCMER":  # Skip headers
                    continue
            if line[0].startswith(">"):
                current_ref = line[0].strip('>')
                current_qry = line[1]
            if len(line) == 7:
                regions_ref[current_ref].append(tuple(sorted(list([int(line[0]), int(line[1])])))) #aligned regions reference
                regions_qry[current_qry].append(tuple(sorted(list([int(line[2]), int(line[3])])))) #aligned regions qry

    #Getting aligned based for reference sequence
    ref_total_aligned_size = 0
    for key in regions_ref:
        ref_tree = intervaltree.IntervalTree.from_tuples(regions_ref[key])
        ref_tree.merge_overlaps(strict=False)
        print(ref_tree)
        ref_aligned_size = 0
        for interval in ref_tree:
            ref_aligned_size += interval.end - interval.begin + 1
        ref_total_aligned_size += ref_aligned_size

    #Getting aligned bases for quey sequence
    query_total_aligned_size = 0
    for key in regions_qry:
        qry_tree = intervaltree.IntervalTree.from_tuples(regions_qry[key])
        qry_tree.merge_overlaps(strict=False)
        qry_aligned_size = 0
        for interval in qry_tree:
            qry_aligned_size += interval.end - interval.begin + 1
        query_total_aligned_size += qry_aligned_size
    
    
    return ref_total_aligned_size, query_total_aligned_size


print(parse_delta(Path("../test_1_incongruencies/out.mdelta")))
