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

from sqlalchemy import UniqueConstraint, create_engine, Table
from sqlalchemy import Column, DateTime, ForeignKey, Integer, String, Float, Boolean
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker

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


class Class(Base):
    """Describes the classification of an input genome"""

    __tablename__ = "classes"
    __table_args__ = (UniqueConstraint("genome_id", "run_id"),)

    class_id = Column(Integer, primary_key=True)
    genome_id = Column(Integer, ForeignKey("genomes.genome_id"), nullable=False)
    run_id = Column(Integer, ForeignKey("runs.run_id"), nullable=False)
    genome_class = Column(String)

    genome = relationship("Genome", back_populates="genome_classes")
    run = relationship("Run", back_populates="run_classes")

    def __str__(self):
        return str(
            "Class, genome {} - run {}: {}".format(
                self.genome_id, self.run_id, self.genome_class
            )
        )

    def __repr__(self):
        return "<Class(class_id={})>".format(self.class_id)


class Label(Base):
    """Describes the label for an input genome to be used in visualisation"""

    __tablename__ = "labels"
    __table_args__ = (UniqueConstraint("genome_id", "run_id"),)

    label_id = Column(Integer, primary_key=True)
    genome_id = Column(Integer, ForeignKey("genomes.genome_id"), nullable=False)
    run_id = Column(Integer, ForeignKey("runs.run_id"), nullable=False)
    genome_label = Column(String)

    genome = relationship("Genome", back_populates="genome_labels")
    run = relationship("Run", back_populates="run_labels")

    def __str__(self):
        return str(
            "Label, genome {} - run {}: {}".format(
                self.genome_id, self.run_id, self.genome_label
            )
        )

    def __repr__(self):
        return "<Label(label_id=[}])>".format(self.label_id)


class Genome(Base):
    """Describes an input genome for a pyani run"""

    __tablename__ = "genomes"
    __table_args__ = (UniqueConstraint("genome_hash"),)

    genome_id = Column(Integer, primary_key=True)
    genome_hash = Column(String)
    path = Column(String)
    length = Column(Integer)
    description = Column(String)

    genome_classes = relationship("Class", back_populates="genome")
    genome_labels = relationship("Label", back_populates="genome")
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
    run_labels = relationship("Label", back_populates="run", lazy="dynamic")
    run_classes = relationship("Class", back_populates="run", lazy="dynamic")

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


if __name__ == "__main__":
    # Create test database if run as script

    import datetime
    import time
    import os

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
        glabel = Label(genome=genome, run=run, genome_label=glabel)
        session.add(glabel)
        gclass = Class(genome=genome, run=run, genome_class=gclass)
        session.add(gclass)
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
