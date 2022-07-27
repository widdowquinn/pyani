import itertools as it
import random
import sys
from typing import List, Tuple, NamedTuple
from pathlib import Path
from Bio import SeqIO, Seq
from . import PyaniException

# Create a list of the nucleotides
ntides = ["A", "C", "G", "T"]


class PyaniEvolveException(PyaniException):
    """Evolve-specific exception for pyani."""


class MutationEvent(NamedTuple):
    length: int
    number: int


class MutatableRecord(SeqIO.SeqRecord):
    def __init__(self, sequence=None, length=1200, **kwargs):
        """Create a custom type of SeqIO.SeqRecord object. If a sequence is
        passed, it is used to create the record; if not, one is randomly
        generated. The length of this can be set, or will default to 1200.

        :param sequence:   a string to be used as the sequence
        :param length:     an int that determines the length of the record's sequence, if one is generated

        Additional keyword arguments derive from SeqIO.SeqRecord:

        Main attributes:
         - id          - Identifier such as a locus tag (string)
         - seq         - The sequence itself (Seq object or similar)

        Additional attributes:
         - name        - Sequence name, e.g. gene name (string)
         - description - Additional text (string)
         - dbxrefs     - List of database cross references (list of strings)
         - features    - Any (sub)features defined (list of SeqFeature objects)
         - annotations - Further information about the whole sequence (dictionary).
        Most entries are strings, or lists of strings.
         - letter_annotations - Per letter/symbol annotation (restricted
        dictionary). This holds Python sequences (lists, strings
        or tuples) whose length matches that of the sequence.
        A typical use would be to hold a list of integers
        representing sequencing quality scores, or a string
        representing the secondary structure.

        """
        # Ensure sequence length is evenly divisible by 12, and
        # in the range x <= length < x + 12, where x is the
        # original value that is passed.
        # We use 12 because itis 3 * 4; 3 for the size of a codon;
        # 4 for the size of the kmers used to generate the sequence
        self.length = (length + 11) // 12 * 12

        # Set sequence (generating, if necessary)
        if not sequence:
            sequence = self.generate_kmer(self.length)
        else:
            sequence = sequence

        # Ensure any other kwargs are set
        SeqIO.SeqRecord.__init__(self, Seq.Seq(sequence), **kwargs)

    def generate_kmer(self, length):
        """Generate a kmer of a given length using randomly-selected 4-mers.

        :param length:   the length of the desired kmer

        """
        # Builds the genome
        kmer = "".join(random.choices(ntides, k=length))

        return kmer

    def delete(self, deletions):
        """Deletes one or more segments from the sequence; the deletions are given as MutationEvent objects, specifying the length of the deletion and number of occurrences.

        :param deletions:  a list of MutationEvent tuples

        """
        sequence = self.seq
        for pair in deletions:
            for event in range(pair.number):
                start = random.choice(range(len(sequence) - pair.length))
                sequence = sequence[:start] + sequence[start + pair.length :]
        self.seq = Seq.Seq(sequence)
        return self

    def mutate(self, mutations: List[NamedTuple]):
        """Mutates one or more segments in the sequence; the mutations are given as MutationEvent objects, specifying the length of the mutation and number of occurrences.

        :param mutations:  a list of MutationEvent tuples

        """
        sequence = self.seq
        for pair in mutations:
            for event in range(pair.number):
                start = random.choice(range(len(sequence) - pair.length))
                mut = self.generate_kmer(pair.length)
                sequence = sequence[:start] + mut + sequence[start + pair.length :]
        self.seq = Seq.Seq(sequence)
        return self

    def insert(self, insertions):
        """Inserts one or more segments in the sequence; the insertions are given as MutationEvent objects, specifying the length of the insertion and number of occurrences.

        :param insertions:  a list of MutationEvent tuples

        """
        sequence = self.seq
        for pair in insertions:
            for event in range(pair.number):
                start = random.choice(range(len(sequence)))
                mut = self.generate_kmer(pair.length)
                sequence = sequence[:start] + mut + sequence[start:]
        self.seq = Seq.Seq(sequence)
        return self

    def invert(self, inversions):
        """Inverts one or more segments in the sequence; the inversions are given as MutationEvent objects, specifying the length of the inversion and number of occurrences.

        :param inversions:  a list of MutationEvent tuples

        """
        sequence = self.seq
        for pair in inversions:
            for event in range(pair.number):
                start = random.choice(range(len(sequence) - pair.length))
                inv = sequence[start : start + pair.length]
                sequence = (
                    sequence[:start] + inv[::-1] + sequence[start + pair.length :]
                )
        self.seq = Seq.Seq(sequence)
        return self

    def repeat(self, repetitions):
        """Repeats one or more segments in the sequence; the repetitions are given as MutationEvent objects, specifying the length of the repetition and number of occurrences.

        :param repetitions:  a list of MutationEvent tuples

        Extra steps are taken to ensure the repeated portions of code do not overlap. This is accoplished by dividing the sequence into partitions and not operating within the same partition more than once.
        """
        sequence = str(self.seq)
        for pair in repetitions:
            sys.stdout.write(f"{sequence}")
            if not len(sequence) > pair.length * pair.number * pair.number:
                return "Can not (currently) guarantee non-overlapping repeats of this length. Using a longer sequence, or requesting shorter repeats will solve this."
            for event in range(pair.number):
                section_length = len(sequence) // (pair.number * pair.number)
                section_start = section_length * event
                section_end = section_length * (event + 1) - pair.length
                start = random.choice(range(section_start, section_end))
                insert = random.choice(
                    range(
                        section_start + (2 * section_length),
                        section_end + (2 * section_length),
                    )
                )
                sequence = (
                    sequence[:insert]
                    + sequence[start : start + pair.length]
                    + sequence[insert:]
                )
        self.seq = Seq.Seq(sequence)
        return self

    def tandem_repeat(self, repetitions):
        """Repeats one or more segments in the sequence, adjacent to their original occurence; the repetitions are given as MutationEvent objects, specifying the length of the repetition and number of occurrences.

        :param repetitions:  a list of MutationEvent tuples

        """
        sequence = self.seq
        for pair in repetitions:
            for event in range(pair.number):
                start = random.choice(range(len(sequence) - pair.length))
                sequence = (
                    sequence[:start]
                    + sequence[start : start + pair.length]
                    + sequence[start:]
                )
        self.seq = Seq.Seq(sequence)
        return self

    def write_to_file(self, file, format="fasta"):
        """Writes the sequence to a file; defaults to `fasta`.

        :param file:  the output file name
        :param format:  the desired output format

        """
        with open(file, "w") as output_handle:
            SeqIO.write(self, output_handle, format)


# def read_from_file(file, format="fasta"):

#    # Check if file exists
#    if not Path(file).resolve().is_file():
#        return f"File not found at {file}"

#    # Try reading from file
#    try:
#        # assume only one sequence in file to start with
#        try:
#            sequence = SeqIO.read(file, format)
#        # try a different option that accommodates multiple sequences
#        except ValueError:
#            sequences = list(SeqIO.parse(file, format))
#    except ValueError:
#        raise PyaniEvolveException(f"{file} is not a valid {format} format")
#    # What should this return? There are conflicting options.


def make_mutatable_record(sequence: str):
    record = MutatableRecord(sequence)
    return record
