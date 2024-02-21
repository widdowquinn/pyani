#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019-2024
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# 161 Cathedral Street,
# Glasgow,
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2019 The James Hutton Institute
# Copyright (c) 2019-2024 University of Strathclyde
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
"""Provides the report subcommand for pyani."""

import logging
import warnings

from argparse import Namespace
from typing import List, NamedTuple
from pathlib import Path

import pandas as pd

from sqlalchemy import and_
from sqlalchemy.orm import aliased

from pyani import pyani_orm, pyani_report
from pyani.pyani_orm import (
    Run,
    Genome,
    Comparison,
    Label,
    get_matrix_labels_for_run,
    rungenome,
)
from pyani.pyani_tools import label_results_matrix, termcolor, MatrixData


class ReportParams(NamedTuple):
    """Report query/header data."""

    name: str  # name of table being reported
    statement: str  # SQL statement of query
    headers: List[str]  # Column headers for table


def subcmd_report(args: Namespace) -> int:
    """Present report on ANI results and/or database contents.

    :param args:  Namespace, command-line arguments

    The report subcommand takes any of several long options that do one of two
    things:

    1. perform a single action.
    2. set a parameter/format

    These will typically take an output path to a file or directory into which
    the report will be written (whatever form it takes). By default, text
    output is written in plain text format, but for some outputs this can
    be modified by an 'excel' or 'html' format specifier, which writes outputs
    in that format, where possible.
    """
    logger = logging.getLogger(__name__)

    # Output formats will apply across all tabular data requested
    # Expect comma-separated format arguments, and turn them into an iterable
    formats = process_formats(args)
    logger.info(termcolor("Creating report output in formats: %s", "red"), formats)

    # Declare which database is being used, and connect to session
    logger.debug("Using database: %s", args.dbpath)
    session = pyani_orm.get_session(args.dbpath)

    # Report runs in the database
    if args.show_runs:
        statement = session.query(
            Run.run_id, Run.name, Run.method, Run.date, Run.cmdline
        ).statement
        headers = ["run ID", "name", "method", "date run", "command-line"]
        report(args, session, formats, ReportParams("runs", statement, headers))

    # Report genomes in the database
    if args.show_genomes:
        statement = session.query(
            Genome.genome_id,
            Genome.description,
            Genome.path,
            Genome.genome_hash,
            Genome.length,
        ).statement
        headers = ["genome ID", "description", "path", "MD5 hash", "genome length"]
        report(args, session, formats, ReportParams("genomes", statement, headers))

    # Report table of all genomes used for each run
    # This is essentially an expanded equivalent of:
    # SELECT genomes.genome_id, runs_genomes.run_id, labels.label_id
    #  FROM genomes
    #    JOIN runs_genomes ON genomes.genome_id == runs_genomes.genome_id
    #    JOIN labels ON (genomes.genome_id = labels.genome_id
    #                    AND runs_genomes.run_id = labels.run_id);
    if args.show_runs_genomes:
        statement = (
            session.query(
                Run.run_id,
                Run.name,
                Run.method,
                Run.date,
                Genome.genome_id,
                Genome.description,
                Genome.path,
                Genome.genome_hash,
                Label.label,
                Label.class_label,
            )
            .join(rungenome, Genome.genome_id == rungenome.c.genome_id)
            .join(
                Label,
                and_(
                    Genome.genome_id == Label.genome_id,
                    rungenome.c.run_id == Label.run_id,
                ),
            )
            .join(Run, Run.run_id == rungenome.c.run_id)
            .order_by(Run.run_id, Genome.genome_id)
            .statement
        )
        headers = [
            "run ID",
            "run name",
            "method",
            "date run",
            "genome ID",
            "genome description",
            "genome path",
            "genome hash",
            "genome label",
            "genome class",
        ]
        report(args, session, formats, ReportParams("runs_genomes", statement, headers))

    # Report table of all runs in which a genome is involved
    # Compiles the equivalent of:
    # SELECT genomes.genome_id, runs.run_id, genomes.description, genomes.path,
    #        genomes.genome_hash, labels.label, labels.class_label,
    #        runs.name, runs.method, runs.date
    #   FROM genomes
    #     JOIN runs_genomes ON genomes.genome_id = runs_genomes.genome_id
    #     JOIN labels ON genomes.genome_id = labels.genome_id
    #                    AND runs_genomes.run_id = labels.run_id
    #     JOIN runs ON runs.run_id = runs_genomes.run_id
    #  ORDER BY genomes.genome_id, runs.run_id
    if args.show_genomes_runs:
        statement = (
            session.query(
                Genome.genome_id,
                Run.run_id,
                Genome.description,
                Genome.path,
                Genome.genome_hash,
                Label.label,
                Label.class_label,
                Run.name,
                Run.method,
                Run.date,
            )
            .join(rungenome, Genome.genome_id == rungenome.c.genome_id)
            .join(
                Label,
                and_(
                    Genome.genome_id == Label.genome_id,
                    rungenome.c.run_id == Label.run_id,
                ),
            )
            .join(Run, Run.run_id == rungenome.c.run_id)
            .order_by(Genome.genome_id, Run.run_id)
            .statement
        )
        headers = [
            "genome ID",
            "run ID",
            "genome description",
            "genome path",
            "genome hash",
            "genome label",
            "genome class",
            "run name",
            "method",
            "date run",
        ]
        report(args, session, formats, ReportParams("genomes_runs", statement, headers))

    # Report table of comparison results for the indicated runs
    if args.run_results:
        run_ids = [run_id for run_id in args.run_results]
        logger.debug("Attempting to write results tables for runs: %s", run_ids)
        for run_id in run_ids:
            logger.debug("Processing run ID %s", run_id)
            genome_query = aliased(Genome, name="genome_query")
            genome_subject = aliased(Genome, name="genome_subject")
            statement = (
                session.query(
                    Comparison.comparison_id,
                    Comparison.query_id,
                    genome_query.description,
                    Comparison.subject_id,
                    genome_subject.description,
                    Comparison.identity,
                    Comparison.cov_query,
                    Comparison.cov_subject,
                    Comparison.aln_length,
                    Comparison.sim_errs,
                    Comparison.program,
                    Comparison.version,
                    Comparison.fragsize,
                    Comparison.maxmatch,
                    Run.run_id,
                )
                .join(genome_query, Comparison.query_id == genome_query.genome_id)
                .join(genome_subject, Comparison.subject_id == genome_subject.genome_id)
                .filter(Run.run_id == run_id)
                .statement
            )
            logger.debug("Results query: %s", statement)
            headers = [
                "Comparison ID",
                "Query ID",
                "Query description",
                "Subject ID",
                "Subject description",
                "% identity",
                "% query coverage",
                "% subject coverage",
                "alignment length",
                "similarity errors",
                "program",
                "version",
                "fragment size",
                "maxmatch",
                "Run ID",
            ]
            report(
                args,
                session,
                formats,
                ReportParams(f"results_{run_id}", statement, headers),
            )

    # Report matrices of comparison results for the indicated runs
    # For ANIm, all results other than coverage are symmetric matrices,
    # so we only get results in the forward direction.
    # As we need to pull down the matrices as Pandas dataframes by reading from
    # JSON, we don't bother with a helper function like report(), and write out
    # our matrices directly, here
    if args.run_matrices:
        show_index = not args.no_matrix_labels
        for run_id in [run_id for run_id in args.run_matrices]:
            logger.debug("Extracting matrices for run %s", run_id)
            run = session.query(Run).filter(Run.run_id == run_id).first()
            matlabel_dict = get_matrix_labels_for_run(session, run_id)
            for matdata in [
                MatrixData(*_)
                for _ in [
                    ("identity", run.df_identity, {"colour_num": 0.95}),
                    ("coverage", run.df_coverage, {"colour_num": 0.95}),
                    ("aln_lengths", run.df_alnlength, {}),
                    ("sim_errors", run.df_simerrors, {}),
                    ("hadamard", run.df_hadamard, {}),
                ]
            ]:
                logger.debug("Writing %s results", matdata.name)
                matrix = pd.read_json(matdata.data)
                # Matrix rows and columns are labelled if there's a label dictionary,
                # and take the dataframe index otherwise
                matrix = label_results_matrix(matrix, matlabel_dict)
                pyani_report.write_dbtable(
                    matrix,
                    Path(
                        "_".join(
                            [str(args.outdir / "matrix"), matdata.name, str(run_id)]
                        )
                    ),
                    formats,
                    show_index=show_index,
                    **matdata.graphic_args,
                )

    return 0


def report(args: Namespace, session, formats: List[str], params: ReportParams) -> None:
    """Write tabular report of pyani runs from database.

    :param args:  Namespace of command-line arguments
    :param session:  SQLAlchemy database session
    :param formats:  list of output formats
    :param params:  ReportParams namedtuple
    """
    logger = logging.getLogger(__name__)

    outfname = args.outdir / params.name
    logger.debug(
        "Writing table of pyani %s from the database to %s.*", params.name, outfname
    )

    # With newer versions of SQLAlchemy, Pandas may throw a warning due to the composition
    # of our statement including a Cartesian product, even though this is what we want:
    #   SAWarning: SELECT statement has a cartesian product between FROM element(s)
    #   "runs" and FROM element "genome_query".  Apply join condition(s) between each
    #   element to resolve.
    # We could use SQLAlchemy's true() function to force the join condition, but this has
    # to be done from within Pandas, and is an issue for them to fix.
    # We suppress the warning, instead.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=(
                "SELECT statement has a cartesian product between FROM "
                'element\\(s\\) "runs" and FROM element "genome_query"'
            ),
        )
        warnings.filterwarnings(
            "ignore",
            message=(
                "SELECT statement has a cartesian product between FROM "
                'element\\(s\\) "runs" and FROM element "comparisons"'
            ),
        )
        warnings.filterwarnings(
            "ignore",
            message=(
                "SELECT statement has a cartesian product between FROM "
                'element\\(s\\) "runs" and FROM element "genome_subject"'
            ),
        )
        warnings.filterwarnings(
            "ignore",
            message=(
                "SELECT statement has a cartesian product between FROM "
                'element\\(s\\) "genome_query", "genome_subject", "comparisons" '
                'and FROM element "runs"'
            ),
        )
        data = pd.read_sql(params.statement, session.bind)

    data.columns = params.headers
    pyani_report.write_dbtable(data, outfname, formats)


def process_formats(args: Namespace) -> List[str]:
    """Return processed list of output formats for writing reports.

    :param args:  Namespace of command-line arguments
    """
    formats = ["tab"]
    if args.formats:
        formats += [fmt for fmt in args.formats]
    return list(set(formats))  # remove duplicates
