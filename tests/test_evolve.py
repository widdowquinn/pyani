from typing import NamedTuple
import sys
import pytest

from pyani.evolve import MutatableRecord, MutationEvent
from Bio import SeqIO, Seq


@pytest.fixture
# fixture
def record():
    return MutatableRecord(
        "CAGAATCACGTTGAGACGAGTACGCCGAGGGTCACAGAACCGGTCCATGTTCCGATTAGG", id="60bp"
    )


@pytest.fixture
def short_record():
    return MutatableRecord("TGACTGTAGGGAAAACATAAAGTG", id="24bp")


@pytest.fixture
def large_mutation_event():
    return [MutationEvent(6, 1)]


@pytest.fixture
def small_mutation_events():
    return [MutationEvent(3, 2)]


@pytest.fixture
def repeats_too_large():
    return [MutationEvent(6, 2)]


# Test mutation (length doesn't change)
def test_mutate_large(record, large_mutation_event):
    record.mutate(large_mutation_event)
    assert len(record.seq) == 60


def test_mutate_small(record, small_mutation_events):
    record.mutate(small_mutation_events)
    assert len(record.seq) == 60


# Test insertion (length increases)
def test_insert_large(record, large_mutation_event):
    record.insert(large_mutation_event)
    assert len(record.seq) == 66


def test_insert_small(record, small_mutation_events):
    record.insert(small_mutation_events)
    assert len(record.seq) == 66


# Test deletion (length decreases)
def test_delete_large(record, large_mutation_event):
    record.delete(large_mutation_event)
    assert len(record.seq) == 54


def test_delete_small(record, small_mutation_events):
    record.delete(small_mutation_events)
    assert len(record.seq) == 54


# Test tandem repeat (length increases)
def test_repeat_large(record, large_mutation_event):
    record.repeat(large_mutation_event)
    assert len(record.seq) == 66


def test_repeat_small(record, small_mutation_events):
    record.repeat(small_mutation_events)
    assert len(record.seq) == 66


# Test repeat (length increases)
def test_tandem_large(record, large_mutation_event):
    record.tandem_repeat(large_mutation_event)
    assert len(record.seq) == 66


def test_tandem_small(record, small_mutation_events):
    record.tandem_repeat(small_mutation_events)
    assert len(record.seq) == 66


# Test sequence too short for repeat size
def test_repeat_failure(short_record, repeats_too_large):
    result = short_record.repeat(repeats_too_large)
    assert (
        result
        == "Can not (currently) guarantee non-overlapping repeats of this length. Using a longer sequence, or requesting shorter repeats will solve this."
    )
