#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv

from collections import defaultdict
from pathlib import Path

import intervaltree

# Parse .coords files and try to replicate AlignedBases
def parse_coords(infname):
    aln_ref = 0
    aln_query = 0
    cov_ref = 0
    cov_query = 0
    num_rows = 0
    last_row = None
    overlaps = 0
    overlap_len = 0

    ref_intervals = defaultdict(list)
    query_intervals = defaultdict(list)

    with infname.open() as ifh:
        fieldnames = [
            "start_ref",
            "end_ref",
            "start_query",
            "end_query",
            "aln_len_ref",
            "aln_len_query",
            "aln_id",
            "len_ref",
            "len_query",
            "cov_ref",
            "cov_query",
            "id_ref",
            "id_query",
        ]
        reader = csv.DictReader(ifh, delimiter="\t", fieldnames=fieldnames)
        for row in reader:
            num_rows += 1
            aln_ref += int(row["aln_len_ref"])
            aln_query += int(row["aln_len_query"])
            cov_ref += float(row["cov_ref"])
            cov_query += float(row["cov_query"])
            key_query = row["id_query"]
            key_ref = row["id_ref"]

            ref_intervals[key_ref].append(sorted((int(row["start_ref"]), int(row["end_ref"]))))
            query_intervals[key_query].append(
                sorted((int(row["start_query"]), int(row["end_query"])))
            )

    ref_total_aligned_size = 0
    for key in ref_intervals:
        ref_tree = intervaltree.IntervalTree.from_tuples(ref_intervals[key])
    #     # print(f"{len(ref_tree)=}")
        ref_tree.merge_overlaps()
        ref_aligned_size = 0
        for interval in ref_tree:
            ref_aligned_size += interval.end - interval.begin + 1
        ref_total_aligned_size += ref_aligned_size




    query_total_aligned_size = 0
    for key in query_intervals:
        query_tree = intervaltree.IntervalTree.from_tuples(query_intervals[key])
        # print(f"{len(query_tree)=}")
        query_tree.merge_overlaps()
        query_aligned_size = 0
        for interval in query_tree:
            query_aligned_size += interval.end - interval.begin + 1
        query_total_aligned_size += query_aligned_size



    return ref_total_aligned_size, query_total_aligned_size


