from pyani.evolve import MutatableRecord, MutationEvent
import random
from pathlib import Path
from typing import List, Tuple, NamedTuple
import argparse

outdir = "../sequences"

ntides = ["A", "C", "G", "T"]


def parse_args():
    description = "Creates a set of fake genome sequences."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--version", action="version", version="%(prog)s dev")
    parser.add_argument(
        "-o", "--outdir", dest="outdir", help=("output directory location")
    )
    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    sequences = {}

    sequence1 = MutatableRecord(
        id="sequence1", description=": base sequence with base sequence"
    )
    base_sequence = sequence1.seq
    sequences["sequence1"] = sequence1

    sequence2 = MutatableRecord(
        base_sequence,
        id="sequence2",
        description=": base sequence with one 60bp mutation",
    )
    sequence2.mutate([MutationEvent(60, 1)])
    sequences["sequence2"] = sequence2

    sequence2a = MutatableRecord(
        base_sequence,
        id="sequence2a",
        description=": base sequence with sixty 1bp mutations",
    )
    sequence2a.mutate([MutationEvent(1, 60)])
    sequences["sequence2a"] = sequence2a

    # sequence2b = MutatableRecord(base_sequence, id="sequence2b", description=": base sequence with twenty 3bp mutations")
    # sequence2b.mutate([MutationEvent(3, 20)])
    # sequences["sequence2b"] = sequence2b

    sequence3 = MutatableRecord(
        base_sequence,
        id="sequence3",
        description=": base sequence with one 60bp insertion",
    )
    sequence3.insert([MutationEvent(60, 1)])
    sequences["sequence3"] = sequence3

    # sequence3a = MutatableRecord(base_sequence, id="sequence3a", description=": base sequence with sixty 1bp insertions")
    # sequence3a.insert([MutationEvent(1, 60)])
    # sequences["sequence3a"] = sequence3a

    sequence3b = MutatableRecord(
        base_sequence,
        id="sequence3b",
        description=": base sequence with twenty 3bp insertions",
    )
    sequence3b.insert([MutationEvent(3, 20)])
    sequences["sequence3b"] = sequence3b

    sequence4 = MutatableRecord(
        base_sequence,
        id="sequence4",
        description=": base sequence with one 60bp deletion",
    )
    sequence4.delete([MutationEvent(60, 1)])
    sequences["sequence4"] = sequence4

    # sequence4a = MutatableRecord(base_sequence, id="sequence4a", description=": base sequence with sixty 1bp deletions")
    # sequence4a.delete([MutationEvent(1, 60)])
    # sequences["sequence4a"] = sequence4a

    sequence4b = MutatableRecord(
        base_sequence,
        id="sequence4b",
        description=": base sequence with twenty 3bp deletions",
    )
    sequence4b.delete([MutationEvent(3, 20)])
    sequences["sequence4b"] = sequence4b

    sequence5 = MutatableRecord(
        base_sequence,
        id="sequence5",
        description=": base sequence with one 60bp repeat",
    )
    sequence5.tandem_repeat([MutationEvent(60, 1)])
    sequences["sequence5"] = sequence5

    sequence5a = MutatableRecord(
        base_sequence,
        id="sequence5a",
        description=": base sequence with two 30bp repeats",
    )
    sequence5a.tandem_repeat([MutationEvent(30, 2)])
    sequences["sequence5a"] = sequence5a

    sequence6 = MutatableRecord(
        base_sequence,
        id="sequence6",
        description=": base sequence with one 60bp tandem repeat",
    )
    sequence6.tandem_repeat([MutationEvent(60, 1)])
    sequences["sequence6"] = sequence6

    sequence6a = MutatableRecord(
        base_sequence,
        id="sequence6a",
        description=": base sequence with two 30bp tandem repeats",
    )
    sequence6a.tandem_repeat([MutationEvent(30, 2)])
    sequences["sequence6a"] = sequence6a

    for id, record in sequences.items():
        record.write_to_file(f"{args.outdir}/{id}.fasta")


if __name__ == "__main__":
    main()


# A version as functions that don't use object-oriented programming.


def generate_kmer(length: int):
    """Generate a kmer of a given length using randomly-selected 4-mers.

    :param length:   the length of the desired kmer

    """
    # Builds the genome
    kmer = "".join(random.choices(ntides, k=length))

    return kmer


def mutate(sequence, mutations: List[NamedTuple]):

    for pair in mutations:
        for event in range(pair.number):
            start = random.choice(range(len(sequence) - pair.length))
            mut = generate_kmer(pair.length)
            sequence = sequence[:start] + mut + sequence[start + pair.length :]
    return sequence


def insert(sequence, insertions):
    for pair in insertions:
        for event in range(pair.number):
            start = random.choice(range(len(sequence)))
            mut = generate_kmer(pair.length)
            sequence = sequence[:start] + mut + sequence[start:]
    return sequence


def delete(sequence, deletions):
    for pair in deletions:
        for event in range(pair.number):
            start = random.choice(range(len(sequence)))
            sequence = sequence[:start] + sequence[start + pair.length :]
    return sequence


def repeat(sequence, repetitions):
    for pair in repetitions:
        for event in range(pair.number):
            start = random.choice(range(len(sequence)))
            sequence = (
                sequence[:start]
                + sequence[start : start + pair.length]
                + sequence[start:]
            )
    return sequence
