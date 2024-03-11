# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2018-2019
# (c) The University of Strathclyde 2019-2024
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute of Pharmacy and Biomedical Sciences
# University of Strathclyde
# 161 Cathedral Street
# Glasgow
# Scotland,
# G4 0RE
# UK
#
# The MIT License
#
# Copyright (c) 2018-2019 The James Hutton Institute
# Copyright (c) 2019-2024 The University of Strathclyde
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""Module providing useful functions for manipulating pyani's SQLite3 db.

This SQLAlchemy-based ORM replaces the previous SQL-based module
"""

import logging

from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Optional, Tuple

import numpy as np  # type: ignore
import pandas as pd  # type: ignore

from sqlalchemy import and_  # type: ignore
import sqlalchemy
from sqlalchemy import UniqueConstraint, create_engine, Table
from sqlalchemy import Column, DateTime, ForeignKey, Integer, String, Float, Boolean

# from sqlalchemy.ext.declarative import declarative_base  # type: ignore
from sqlalchemy.orm import declarative_base, relationship, sessionmaker  # type: ignore

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


# Using the declarative system
# We follow Flask-like naming conventions, so override some of the pyline errors
# Mypy doesn't like dynamic base classes, see https://github.com/python/mypy/issues/2477
Base = declarative_base()  # type: Any
Session = sessionmaker()  # pylint: disable=C0103

# Linker table between genomes and runs tables
rungenome = Table(  # pylint: disable=C0103
    "runs_genomes",
    Base.metadata,
    Column("genome_id", Integer, ForeignKey("genomes.genome_id")),
    Column("run_id", Integer, ForeignKey("runs.run_id")),
)

# Linker table between comparisons and runs tables
runcomparison = Table(  # pylint: disable=C0103
    "runs_comparisons",
    Base.metadata,
    Column("comparison_id", Integer, ForeignKey("comparisons.comparison_id")),
    Column("run_id", Integer, ForeignKey("runs.run_id")),
)


# Convenience struct for labels and classes
class LabelTuple(NamedTuple):
    """Label and Class for each file."""

    label: str
    class_label: str


class Label(Base):
    """Describes relationship between genome, run and genome label.

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

    def __str__(self) -> str:
        """Return string representation of Label table row."""
        return str(
            "Genome ID: {}, Run ID: {}, Label ID: {}, Label: {}, Class: {}".format(
                self.genome_id, self.run_id, self.label_id, self.label, self.class_label
            )
        )

    def __repr__(self) -> str:
        """Return string representation of Label table object."""
        return "<Label(key=({}, {}, {}))>".format(
            self.label_id, self.run_id, self.genome_id
        )


class BlastDB(Base):
    """Describes relationship between genome, run, source BLAST database and query fragments.

    Each genome and run combination can be assigned a single BLAST database
    for the comparisons

    - fragpath      path to fragmented genome (query in ANIb)
    - dbpath        path to source genome database (subject in ANIb)
    - fragsizes     JSONified dict of fragment sizes
    - dbcmd         command used to generate database
    """

    __tablename__ = "blastdbs"

    blastdb_id = Column(Integer, primary_key=True)
    genome_id = Column(Integer, ForeignKey("genomes.genome_id"))
    run_id = Column(Integer, ForeignKey("runs.run_id"))
    fragpath = Column(String)
    dbpath = Column(String)
    fragsizes = Column(String)
    dbcmd = Column(String)

    genome = relationship("Genome", back_populates="blastdbs")
    run = relationship("Run", back_populates="blastdbs")

    def __str__(self) -> str:
        """Return string representation of BlastDB table row."""
        return str(
            "BlastDB: {}, Run ID: {}, Label ID: {}, Label: {}, Class: {}".format(
                self.genome_id, self.run_id, self.label_id, self.label, self.class_label
            )
        )

    def __repr__(self) -> str:
        """Return string representation of BlastDB table object."""
        return "<BlastDB(key=({}, {}, {}))>".format(
            self.label_id, self.run_id, self.genome_id
        )


class Genome(Base):
    """Describes an input genome for a pyani run.

    - genome_id
        primary key
    - genome_hash
        MD5 hash of input genome file (in ``path``)
    - path
        path to FASTA genome file
    - length
        length of genome (total bases)
    - description
        genome description
    """

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

    def __str__(self) -> str:
        """Return string representation of Genome table row."""
        return str("Genome {}: {}".format(self.genome_id, self.description))

    def __repr__(self) -> str:
        """Return string representation of Genome table object."""
        return "<Genome(id='{}',desc='{}')>".format(self.genome_id, self.description)


class Run(Base):
    """Describes a single pyani run."""

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

    def __str__(self) -> str:
        """Return string representation of Run table row."""
        return str("Run {}: {} ({})".format(self.run_id, self.name, self.date))

    def __repr__(self) -> str:
        """Return string representation of Run table object."""
        return "<Run(run_id={})>".format(self.run_id)


class Comparison(Base):
    """Describes a single pairwise comparison between two genomes."""

    __tablename__ = "comparisons"
    __table_args__ = (
        UniqueConstraint(
            "query_id",
            "subject_id",
            "program",
            "version",
            "fragsize",
            "maxmatch",
            "kmersize",
            "minmatch",
        ),
    )

    comparison_id = Column(Integer, primary_key=True)
    query_id = Column(Integer, ForeignKey("genomes.genome_id"), nullable=False)
    subject_id = Column(Integer, ForeignKey("genomes.genome_id"), nullable=False)
    query_aln_length = Column(Integer)  # in fastANI this is matchedfrags * fragLength
    subject_aln_length = Column(Integer)
    sim_errs = Column(Integer)  # in fastANI this is allfrags - matchedfrags
    perc_id = Column(Float)
    cov_query = Column(Float)  # in fastANI this is matchedfrags/allfrags
    cov_subject = Column(Float)  # in fastANI this is Null
    program = Column(String)
    version = Column(String)
    fragsize = Column(Integer)  # in fastANI this is fragLength
    maxmatch = Column(Boolean)  # in fastANi this is Null
    kmersize = Column(Integer)
    minmatch = Column(Float)

    query = relationship(
        "Genome", foreign_keys=[query_id], back_populates="query_comparisons"
    )
    subject = relationship(
        "Genome", foreign_keys=[subject_id], back_populates="subject_comparisons"
    )
    runs = relationship(
        "Run", secondary=runcomparison, back_populates="comparisons", lazy="dynamic"
    )

    def __str__(self) -> str:
        """Return string representation of Comparison table row."""
        return str(
            "Query: {}, Subject: {}, %%ID={}, ({} {}), FragSize: {}, MaxMatch: {}, KmerSize: {}, MinMatch: {}".format(
                self.query_id,
                self.subject_id,
                self.perc_id,
                self.program,
                self.version,
                self.fragsize,
                self.maxmatch,
                self.kmersize,
                self.minmatch,
            )
        )

    def __repr__(self) -> str:
        """Return string representation of Comparison table object."""
        return "<Comparison(comparison_id={})>".format(self.comparison_id)


def create_db(dbpath: Path) -> None:
    """Create an empty pyani SQLite3 database at the passed path.

    :param dbpath:  path to pyani database
    """
    engine = create_engine("sqlite:///{}".format(dbpath), echo=False)
    Base.metadata.create_all(engine)


def get_session(dbpath: Path) -> Any:
    """Connect to an existing pyani SQLite3 database and return a session.

    :param dbpath: path to pyani database
    """
    engine = create_engine("sqlite:///{}".format(dbpath), echo=False)
    Session.configure(bind=engine)
    return Session()


def get_comparison_dict(session: Any) -> Dict[Tuple, Any]:
    """Return a dictionary of comparisons in the session database.

    :param session:      live SQLAlchemy session of pyani database

    Returns Comparison objects, keyed by (_.query_id, _.subject_id,
    _.program, _.version, _.fragsize, _.maxmatch) tuple
    """

    return {
        (
            _.query_id,
            _.subject_id,
            _.program,
            _.version,
            _.fragsize,
            _.maxmatch,
            _.kmersize,
            _.minmatch,
        ): _
        for _ in session.query(Comparison).all()
    }


def get_matrix_labels_for_run(session: Any, run_id: int) -> Dict:
    """Return dictionary of genome labels, keyed by row/column ID.

    :param session:  live SQLAlchemy session
    :param run_id:   the Run.run_id value for matrices

    The labels should be valid for identity, coverage and other complete
    matrix results accessed via the .df_* attributes of a run.

    Labels are returned keyed by the string of the genome ID, for compatibility with
    matplotlib.
    """
    results = (
        session.query(Genome.genome_id, Label.label)
        .join(rungenome, Genome.genome_id == rungenome.c.genome_id)
        .join(
            Label, and_(Genome.genome_id == Label.genome_id, Run.run_id == Label.run_id)
        )
        .filter(Run.run_id == run_id)
        .all()
    )
    return {str(_.genome_id): _.label for _ in results}


def get_matrix_classes_for_run(session: Any, run_id: int) -> Dict[str, List]:
    """Return dictionary of genome classes, keyed by row/column ID.

    :param session:  live SQLAlchemy session
    :param run_id:   the Run.run_id value for matrices

    The class labels should be valid for identity, coverage and other complete
    matrix results accessed via the .df_* attributes of a run

    Labels are returned keyed by the string of the genome ID, for compatibility with
    matplotlib.
    """
    results = (
        session.query(Genome.genome_id, Label.class_label)
        .join(rungenome, Genome.genome_id == rungenome.c.genome_id)
        .join(
            Label, and_(Genome.genome_id == Label.genome_id, Run.run_id == Label.run_id)
        )
        .filter(Run.run_id == run_id)
        .all()
    )
    return {str(_.genome_id): _.class_label for _ in results}


def filter_existing_comparisons(
    session,
    run,
    comparisons,
    program,
    version,
    fragsize: Optional[int] = None,
    maxmatch: Optional[bool] = False,
    kmersize: Optional[int] = None,
    minmatch: Optional[float] = None,
) -> List:
    """Filter list of (Genome, Genome) comparisons for those not in the session db.

    :param session:       live SQLAlchemy session of pyani database
    :param run:           Run object describing parent pyani run
    :param comparisons:   list of (Genome, Genome) query vs subject comparisons
    :param program:       program used for comparison
    :param version:       version of program for comparison
    :param fragsize:      fragment size for BLAST databases
    :param maxmatch:      maxmatch used with nucmer comparison

    When passed a list of (Genome, Genome) comparisons as comparisons, check whether
    the comparison exists in the database and, if so, associate it with the passed run.
    If not, then add the (Genome, Genome) pair to a list for returning as the
    comparisons that still need to be run.
    """
    logger = logging.getLogger(__name__)

    existing_comparisons = get_comparison_dict(session)
    logger.debug("Existing comparisons\n\t%s", existing_comparisons)
    comparisons_to_run = []
    logger.debug(
        (
            "Checking for existing comparisons, with unique constraints \n"
            "\tprogram: %s\n"
            "\tversion: %s\n"
            "\tfragsize: %s\n"
            "\tmaxmatch: %s\n"
            "\tkmersize: %s\n"
            "\tminmatch: %s\n"
        ),
        program,
        version,
        fragsize,
        maxmatch,
        kmersize,
        minmatch,
    )
    for qgenome, sgenome in comparisons:
        logger.debug(
            "Checking for existing comparison: %s (%s) vs %s (%s)",
            qgenome,
            qgenome.genome_id,
            sgenome,
            sgenome.genome_id,
        )
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
                        kmersize,
                        minmatch,
                    )
                ]
            )
            session.commit()
        except KeyError:
            comparisons_to_run.append((qgenome, sgenome))
    return comparisons_to_run


def add_run(session, method, cmdline, date, status, name):
    """Create a new Run and add it to the session.

    :param session:      live SQLAlchemy session of pyani database
    :param method:       string describing analysis run type
    :param cmdline:      string describing pyani command-line for run
    :param date:         datetime object describing analysis start time
    :param status:       string describing status of analysis
    :param name:         string - name given to the analysis run

    Creates a new Run object with the passed parameters, and returns it.
    """
    try:
        run = Run(method=method, cmdline=cmdline, date=date, status=status, name=name)
    except Exception:
        raise PyaniORMException(
            f"Could not create {method} run with command line: {cmdline}"
        )
    try:
        session.add(run)
        session.commit()
    except Exception:
        raise PyaniORMException(f"Could not add run {run} to the database")
    return run, run.run_id


def add_run_genomes(
    session, run, indir: Path, classpath: Path, labelpath: Path, **kwargs
) -> List:
    """Add genomes for a run to the database.

    :param session:       live SQLAlchemy session of pyani database
    :param run:           Run object describing the parent pyani run
    :param indir:         path to the directory containing genomes
    :param classpath:     path to the file containing class information for each genome
    :param labelpath:     path to the file containing class information for each genome

    This function expects a single directory (indir) containing all FASTA files
    for a run, and optional paths to plain text files that contain information
    on class and label strings for each genome.

    If the genome already exists in the database, then a Genome object is recovered
    from the database. Otherwise, a new Genome object is created. All Genome objects
    will be associated with the passed Run object.

    The session changes are committed once all genomes and labels are added to the
    database without error, as a single transaction.
    """
    # Get list of genome files and paths to class and labels files
    infiles = get_fasta_and_hash_paths(indir)  # paired FASTA/hash files
    class_data = {}  # type: Dict[str,str]
    label_data = {}  # type: Dict[str,str]
    all_keys = []  # type: List[str]
    if classpath:
        class_data = load_classes_labels(classpath)
        all_keys += list(class_data.keys())
    if labelpath:
        label_data = load_classes_labels(labelpath)
        all_keys += list(label_data.keys())

    # Make dictionary of labels and/or classes
    new_keys = set(all_keys)
    label_dict = {}  # type: Dict
    for key in new_keys:
        label_dict[key] = LabelTuple(label_data[key] or "", class_data[key] or "")

    # Get hash and sequence description for each FASTA/hash pair, and add
    # to current session database
    genome_ids = []
    for fastafile, hashfile in infiles:
        try:
            inhash, _ = read_hash_string(hashfile)
            indesc = read_fasta_description(fastafile)
        except Exception:
            raise PyaniORMException("Could not read genome files for database import")
        abspath = fastafile.absolute()
        genome_len = get_genome_length(abspath)
        # If the genome is not already in the database, add it as a Genome object
        genome = session.query(Genome).filter(Genome.genome_hash == inhash).first()
        if not isinstance(genome, Genome):
            try:
                genome = Genome(
                    genome_hash=inhash,
                    path=str(abspath),
                    length=genome_len,
                    description=indesc,
                )
                session.add(genome)
            except Exception:
                raise PyaniORMException(f"Could not add genome {genome} to database")

        # Associate this genome with the current run
        try:
            genome.runs.append(run)
        except Exception:
            raise PyaniORMException(
                f"Could not associate genome {genome} with run {run}"
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
                raise PyaniORMException(
                    f"Could not add labels for {genome} to database."
                )
        genome_ids.append(genome.genome_id)

    try:
        session.commit()
    except Exception:
        raise PyaniORMException("Could not commit new genomes in database.")

    return genome_ids


def update_comparison_matrices(session, run) -> None:
    logger = logging.getLogger(__name__)

    """Update the Run table with summary matrices for the analysis.

    :param session:       active pyanidb session via ORM
    :param run:           Run ORM object for the current ANIm run
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

    logger.debug("Existing comparisons\n%s", run.comparisons.all())
    for cmp in run.comparisons.all():
        qid, sid = cmp.query_id, cmp.subject_id
        df_identity.loc[qid, sid] = cmp.perc_id
        df_coverage.loc[qid, sid] = cmp.cov_query
        df_alnlength.loc[qid, sid] = cmp.query_aln_length
        df_simerrors.loc[qid, sid] = cmp.sim_errs
        df_hadamard.loc[qid, sid] = cmp.perc_id * cmp.cov_query
        if (qid, sid) == (sid, qid):
            df_hadamard.loc[sid, qid] = cmp.perc_id * cmp.cov_subject
            df_simerrors.loc[sid, qid] = cmp.sim_errs
            df_alnlength.loc[sid, qid] = cmp.query_aln_length
            df_coverage.loc[sid, qid] = cmp.cov_query
            df_identity.loc[sid, qid] = cmp.perc_id

    # Add matrices to the database
    run.df_identity = df_identity.to_json()
    run.df_coverage = df_coverage.to_json()
    run.df_alnlength = df_alnlength.to_json()
    run.df_simerrors = df_simerrors.to_json()
    run.df_hadamard = df_hadamard.to_json()
    session.commit()
