# -*- coding: utf-8 -*-
"""Module providing useful functions for manipulating pyani's SQLite3 db.

This SQLAlchemy-based ORM replaces the previous SQL-based module

(c) The James Hutton Institute 2018
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD2 5DA,
Scotland,
UK

The MIT License

Copyright (c) 2018 The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os

from collections import namedtuple

import numpy as np
import pandas as pd

from sqlalchemy import UniqueConstraint, create_engine, Table
from sqlalchemy import Column, DateTime, ForeignKey, Integer, String, Float, Boolean
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker

from pyani import PyaniException
from pyani.pyani_files import (
    get_fasta_and_hash_paths,
    load_classes_labels,
    read_fasta_description,
    read_hash_string,
)
from pyani.pyani_tools import get_genome_length


class PyaniORMException(PyaniException):
    """Exception raised when ORM or database interaction fails."""

    def __init__(self, msg="Error in pyani ORM/database interface"):
        PyaniException.__init__(self, msg)


# Using the declarative system
Base = declarative_base()
Session = sessionmaker()

rungenome = Table(
    "runs_genomes",
    Base.metadata,
    Column("genome_id", Integer, ForeignKey("genomes.genome_id")),
    Column("run_id", Integer, ForeignKey("runs.run_id")),
)

runcomparison = Table(
    "runs_comparisons",
    Base.metadata,
    Column("comparison_id", Integer, ForeignKey("comparisons.comparison_id")),
    Column("run_id", Integer, ForeignKey("runs.run_id")),
)


# convenience namedtuples
label_tuple = namedtuple("ClassData", "label class_label")


class Label(Base):
    """Describes relationship between genome, run and genome label

    Each genome and run combination can be assigned a single label
    """

    __tablename__ = "labels"

    label_id = Column(Integer, primary_key=True)
    genome_id = Column(Integer, ForeignKey("genomes.genome_id"))
    run_id = Column(Integer, ForeignKey("runs.run_id"))
    label = Column(String)
    class_label = Column(String)

    genome = relationship("Genome", back_populates="labels")
    run = relationship("Run", back_populates="labels")

    def __str__(self):
        return str(
            "Genome ID: {}, Run ID: {}, Label ID: {}, Label: {}, Class: {}".format(
                self.genome_id, self.run_id, self.label_id, self.label, self.class_label
            )
        )

    def __repr__(self):
        return "<Label(key=({}, {}, {}))>".format(
            self.label_id, self.run_id, self.genome_id
        )


class BlastDB(Base):
    """Describes relationship between genome, run and BLAST database

    Each genome and run combination can be assigned a single BLAST database
    for the comparisons

    path      path to database files
    size      number of fragment sequences in database
    dbcmd     command used to generate database
    """

    __tablename__ = "blastdb"

    blastdb_id = Column(Integer, primary_key=True)
    genome_id = Column(Integer, ForeignKey("genomes.genome_id"))
    run_id = Column(Integer, ForeignKey("runs.run_id"))
    path = Column(String)
    size = Column(Integer)
    dbcmd = Column(String)

    genome = relationship("Genome", back_populates="blastdbs")
    run = relationship("Run", back_populates="blastdbs")

    def __str__(self):
        return str(
            "Genome ID: {}, Run ID: {}, Label ID: {}, Label: {}, Class: {}".format(
                self.genome_id, self.run_id, self.label_id, self.label, self.class_label
            )
        )

    def __repr__(self):
        return "<Label(key=({}, {}, {}))>".format(
            self.label_id, self.run_id, self.genome_id
        )


class Genome(Base):
    """Describes an input genome for a pyani run"""

    __tablename__ = "genomes"
    __table_args__ = (UniqueConstraint("genome_hash"),)

    genome_id = Column(Integer, primary_key=True)
    genome_hash = Column(String)
    path = Column(String)
    length = Column(Integer)
    description = Column(String)

    labels = relationship("Label", back_populates="genome", lazy="dynamic")
    blastdbs = relationship("BlastDB", back_populates="genome", lazy="dynamic")
    runs = relationship(
        "Run", secondary=rungenome, back_populates="genomes", lazy="dynamic"
    )
    query_comparisons = relationship(
        "Comparison",
        back_populates="query",
        primaryjoin="Genome.genome_id == Comparison.query_id",
    )
    subject_comparisons = relationship(
        "Comparison",
        back_populates="subject",
        primaryjoin="Genome.genome_id == Comparison.subject_id",
    )

    def __str__(self):
        return str("Genome {}: {}".format(self.genome_id, self.description))

    def __repr__(self):
        return "<Genome(id='{}',desc='{}')>".format(self.genome_id, self.description)


class Run(Base):
    """Describes a single pyani run"""

    __tablename__ = "runs"

    run_id = Column(Integer, primary_key=True)
    method = Column(String)
    cmdline = Column(String)
    date = Column(DateTime)
    status = Column(String)
    name = Column(String)
    df_identity = Column(String)  # JSON-encoded Pandas dataframe
    df_coverage = Column(String)  # JSON-encoded Pandas dataframe
    df_alnlength = Column(String)  # JSON-encoded Pandas dataframe
    df_simerrors = Column(String)  # JSON-encoded Pandas dataframe
    df_hadamard = Column(String)  # JSON-encoded Pandas dataframe

    genomes = relationship(
        "Genome", secondary=rungenome, back_populates="runs", lazy="dynamic"
    )
    comparisons = relationship(
        "Comparison", secondary=runcomparison, back_populates="runs", lazy="dynamic"
    )
    labels = relationship("Label", back_populates="run", lazy="dynamic")
    blastdbs = relationship("BlastDB", back_populates="run", lazy="dynamic")

    def __str__(self):
        return str("Run {}: {} ({})".format(self.run_id, self.name, self.date))

    def __repr__(self):
        return "<Run(run_id={})>".format(self.run_id)


class Comparison(Base):
    """Describes a single pairwise comparison between two genomes"""

    __tablename__ = "comparisons"
    __table_args__ = (
        UniqueConstraint(
            "query_id", "subject_id", "program", "version", "fragsize", "maxmatch"
        ),
    )

    comparison_id = Column(Integer, primary_key=True)
    query_id = Column(Integer, ForeignKey("genomes.genome_id"), nullable=False)
    subject_id = Column(Integer, ForeignKey("genomes.genome_id"), nullable=False)
    aln_length = Column(Integer)
    sim_errs = Column(Integer)
    identity = Column(Float)
    cov_query = Column(Float)
    cov_subject = Column(Float)
    program = Column(String)
    version = Column(String)
    fragsize = Column(Integer)
    maxmatch = Column(Boolean)

    query = relationship(
        "Genome", foreign_keys=[query_id], back_populates="query_comparisons"
    )
    subject = relationship(
        "Genome", foreign_keys=[subject_id], back_populates="subject_comparisons"
    )
    runs = relationship(
        "Run", secondary=runcomparison, back_populates="comparisons", lazy="dynamic"
    )

    def __str__(self):
        return str(
            "Query: {}, Subject: {}, %%ID={}, ({} {})".format(
                self.query_id,
                self.subject_id,
                self.identity,
                self.program,
                self.version,
            )
        )

    def __repr__(self):
        return "<Comparison(comparison_id={})>".format(self.comparison_id)


def create_db(dbpath):
    """Create an empty pyani SQLite3 database at the passed path"""
    engine = create_engine("sqlite:///{}".format(dbpath), echo=False)
    Base.metadata.create_all(engine)


def get_session(dbpath):
    """Connect to an existing pyani SQLite3 database and return a session"""
    engine = create_engine("sqlite:///{}".format(dbpath), echo=False)
    Session.configure(bind=engine)
    return Session()


def get_comparison_dict(session):
    """Return a dictionary of comparisons in the session database

    session      live SQLAlchemy session of pyani database

    Returns Comparison objects, keyed by (_.query_id, _.subject_id,
    _.program, _.version, _.fragsize, _.maxmatch) tuple
    """
    return {
        (_.query_id, _.subject_id, _.program, _.version, _.fragsize, _.maxmatch): _
        for _ in session.query(Comparison).all()
    }


def filter_existing_comparisons(
    session, run, comparisons, program, version, fragsize=None, maxmatch=None
):
    """Filter list of (Genome, Genome) comparisons for those not in the session db.

    session       live SQLAlchemy session of pyani database
    run           Run object describing parent pyani run
    comparisons   list of (Genome, Genome) query vs subject comparisons
    program       program used for comparison
    version       version of program for comparison
    fragsize      fragment size for BLAST databases
    maxmatch      maxmatch used with nucmer comparison

    When passed a list of (Genome, Genome) comparisons as comparisons, check whether
    the comparison exists in the database and, if so, associate it with the passed run.
    If not, then add the (Genome, Genome) pair to a list for returning as the
    comparisons that still need to be run.
    """
    existing_comparisons = get_comparison_dict(session)
    comparisons_to_run = []
    for (qgenome, sgenome) in comparisons:
        try:
            # Associate run with existing comparisons
            run.comparisons.append(
                existing_comparisons[
                    (
                        qgenome.genome_id,
                        sgenome.genome_id,
                        program,
                        version,
                        fragsize,
                        maxmatch,
                    )
                ]
            )
            session.commit()
        except KeyError:
            comparisons_to_run.append((qgenome, sgenome))
    return comparisons_to_run


def add_run_genomes(session, run, indir, classpath=None, labelpath=None):
    """Add genomes for a run to the database

    session       live SQLAlchemy session of pyani database
    run           Run object describing the parent pyani run
    indir         path to the directory containing genomes
    classpath     path to the file containing class information for each genome
    labelpath     path to the file containing class information for each genome

    This function expects a single directory (indir) containing all FASTA files
    for a run, and optional paths to plain text files that contain information
    on class and label strings for each genome.

    The function will attempt to associate new Genome objects with the passed
    Run object.

    The session changes are committed once genomes and labels are added to the
    database without error, as a single transaction.
    """
    # Get list of genome files and paths to class and labels files
    infiles = get_fasta_and_hash_paths(indir)  # paired FASTA/hash files
    class_data, label_data, all_keys = None, None, []
    if classpath:
        class_data = load_classes_labels(classpath)
        all_keys += list(class_data.keys())
    if labelpath:
        label_data = load_classes_labels(labelpath)
        all_keys += list(label_data.keys())

    # Make dictionary of labels and/or classes
    new_keys = set(all_keys)
    label_dict = {}
    for key in new_keys:
        label_dict[key] = label_tuple(label_data[key] or None, class_data[key] or None)

    # Get hash and sequence description for each FASTA/hash pair, and add
    # to current session database
    for fastafile, hashfile in infiles:
        try:
            inhash, _ = read_hash_string(hashfile)
            indesc = read_fasta_description(fastafile)
        except Exception:
            raise PyaniORMException("Could not read genome files for database import")
        abspath = os.path.abspath(fastafile)
        genome_len = get_genome_length(abspath)

        # If the genome is not already in the database, add it as a Genome object
        genome = session.query(Genome).filter(Genome.genome_hash == inhash).first()
        if not isinstance(genome, Genome):
            try:
                genome = Genome(
                    genome_hash=inhash,
                    path=abspath,
                    length=genome_len,
                    description=indesc,
                )
                session.add(genome)
            except Exception:
                raise PyaniORMException(
                    "Could not add genome {} to database".format(genome)
                )

        # Associate this genome with the current run
        try:
            genome.runs.append(run)
        except Exception:
            raise PyaniORMException(
                "Could not associate genome {} with run {}".format(genome, run)
            )

        # If there's an associated class or label for the genome, add it
        if inhash in label_dict:
            try:
                session.add(
                    Label(
                        genome=genome,
                        run=run,
                        label=label_dict[inhash].label,
                        class_label=label_dict[inhash].class_label,
                    )
                )
            except Exception:
                raise PyaniORMException("Could not add new genome labels to database.")
        try:
            session.commit()
        except Exception:
            raise PyaniORMException("Could not commit new genomes in database.")


def update_comparison_matrices(session, run):
    """Update the Run table with summary matrices for the analysis.

    session       active pyanidb session via ORM
    run           Run ORM object for the current ANIm run
    """
    # Create dataframes for storing in the Run table
    # Rows and columns are the (ordered) list of genome IDs
    genome_ids = sorted([_.genome_id for _ in run.genomes.all()])
    df_identity = pd.DataFrame(index=genome_ids, columns=genome_ids, dtype=float)
    df_coverage = pd.DataFrame(index=genome_ids, columns=genome_ids, dtype=float)
    df_alnlength = pd.DataFrame(index=genome_ids, columns=genome_ids, dtype=float)
    df_simerrors = pd.DataFrame(index=genome_ids, columns=genome_ids, dtype=float)
    df_hadamard = pd.DataFrame(index=genome_ids, columns=genome_ids, dtype=float)

    # Set appropriate diagonals for each matrix
    np.fill_diagonal(df_identity.values, 1.0)
    np.fill_diagonal(df_coverage.values, 1.0)
    np.fill_diagonal(df_simerrors.values, 1.0)
    np.fill_diagonal(df_hadamard.values, 1.0)
    for genome in run.genomes.all():
        df_alnlength.loc[genome.genome_id, genome.genome_id] = genome.length

    # Loop over all comparisons for the run and fill in result matrices
    for cmp in run.comparisons.all():
        qid, sid = cmp.query_id, cmp.subject_id
        df_identity.loc[qid, sid] = cmp.identity
        df_identity.loc[sid, qid] = cmp.identity
        df_coverage.loc[qid, sid] = cmp.cov_query
        df_coverage.loc[sid, qid] = cmp.cov_subject
        df_alnlength.loc[qid, sid] = cmp.aln_length
        df_alnlength.loc[sid, qid] = cmp.aln_length
        df_simerrors.loc[qid, sid] = cmp.sim_errs
        df_simerrors.loc[sid, qid] = cmp.sim_errs
        df_hadamard.loc[qid, sid] = cmp.identity * cmp.cov_query
        df_hadamard.loc[sid, qid] = cmp.identity * cmp.cov_subject

    # Add matrices to the database
    run.df_identity = df_identity.to_json()
    run.df_coverage = df_coverage.to_json()
    run.df_alnlength = df_alnlength.to_json()
    run.df_simerrors = df_simerrors.to_json()
    run.df_hadamard = df_hadamard.to_json()
    session.commit()


if __name__ == "__main__":
    # Create test database if run as script

    import datetime
    import time

    dbpath = os.path.join(".", "test.sqlite")
    if os.path.isfile(dbpath):
        os.remove(dbpath)
    engine = create_engine("sqlite:///{}".format(dbpath), echo=False)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)

    session = Session()

    # Create test run
    run = Run(
        method="ANIm",
        cmdline="pyani anim C_blochmannia C_blochmannia_ANIm --name C. blochmannia run 1 --labels C_blochmannia/labels.txt --classes C_blochmannia/classes.txt",
        date=datetime.datetime.fromtimestamp(time.time()),
        status="started",
        name="C. blochmannia run 1",
    )
    session.add(run)
    session.commit()

    # Create test genomes
    genome1 = Genome(
        genome_hash="c8c15b0f79742e14afeee07cd63cead1",
        path="C_blochmannia/GCF_000973545.1_ASM97354v1_genomic.fna",
        length=773940,
        description="NZ_CP010049.1 Blochmannia endosymbiont of Camponotus (Colobopsis) obliquus strain 757, complete genome",
    )
    genome1.runs.append(run)
    session.add(genome1)
    genome2 = Genome(
        genome_hash="562fb00b9325563bf352512530897161",
        path="C_blochmannia/GCF_000011745.1_ASM1174v1_genomic.fna",
        length=791654,
        description="NC_007292.1 Candidatus Blochmannia pennsylvanicus str. BPEN, complete genome",
    )
    genome2.runs.append(run)
    session.add(genome2)
    session.commit()

    # Add label and class for each genome
    for (genome, run, glabel, gclass) in [
        (genome1, run, "757", "C. blochmannia"),
        (genome2, run, "BPEN", "C. blochmannia"),
    ]:
        session.add(Label(genome=genome, run=run, label=glabel, class_label=gclass))
    session.commit()

    # Add a comparison for the two genomes
    comparison = Comparison(
        query=genome1,
        subject=genome2,
        aln_length=16627,
        sim_errs=2838,
        identity=0.8293137667649,
        cov_query=0.0214835775383105,
        cov_subject=0.0210028623615873,
        program="nucmer",
        version="3.1",
        fragsize=0,
        maxmatch=False,
    )
    session.add(comparison)
    session.commit()

    # Test that we can use the objects/tables
    our_genomes = session.query(Genome).all()
    print("Genomes:")
    print("\n".join([str(_) for _ in our_genomes]))

    print("\nComparisons by genome")
    for genome in our_genomes:
        print("{}:".format(genome.description))
        print("Query: {}".format(genome.query_comparisons))
        print("Subject: {}".format(genome.subject_comparisons))
