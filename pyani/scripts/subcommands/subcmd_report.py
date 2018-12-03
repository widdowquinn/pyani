#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""subcmd_report.py

Provides the report subcommand for pyani

(c) The James Hutton Institute 2017-18

Author: Leighton Pritchard
Contact: leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2017-18 The James Hutton Institute

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

import pandas as pd

from sqlalchemy import and_
from sqlalchemy.orm import aliased

from pyani import pyani_orm, pyani_report, pyani_db
from pyani.pyani_orm import Run, Genome, Comparison, Label, rungenome, runcomparison


def subcmd_report(args, logger):
    """Present report on ANI results and/or database contents.

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
    # Output formats will apply across all tabular data requested
    # Expect comma-separated format arguments, and turn them into an iterable
    formats = ["tab"]
    if args.formats:
        formats += [fmt.strip() for fmt in args.formats.split(",")]
    formats = list(set(formats))  # remove duplicates
    logger.info("Creating output in formats: %s", formats)

    # Declare which database is being used, and connect to the session
    logger.info("Using database: %s", args.dbpath)
    session = pyani_orm.get_session(args.dbpath)

    # Report runs in the database
    if args.show_runs:
        outfname = os.path.join(args.outdir, "runs")
        logger.info("Writing table of pyani runs from the database to %s.*", outfname)
        statement = session.query(
            Run.run_id, Run.name, Run.method, Run.date, Run.cmdline
        ).statement
        data = pd.read_sql(statement, session.bind)
        headers = ["run ID", "name", "method", "date run", "command-line"]
        data.columns = headers
        pyani_report.write_dbtable(data, outfname, formats)

    # Report genomes in the database
    if args.show_genomes:
        outfname = os.path.join(args.outdir, "genomes")
        logger.info("Writing table of genomes from the database to %s.*", outfname)
        statement = session.query(
            Genome.genome_id,
            Genome.description,
            Genome.path,
            Genome.genome_hash,
            Genome.length,
        ).statement
        data = pd.read_sql(statement, session.bind)
        headers = ["genome ID", "description", "path", "MD5 hash", "genome length"]
        data.columns = headers
        pyani_report.write_dbtable(data, outfname, formats)

    # Report table of all genomes used for each run
    if args.show_runs_genomes:
        outfname = os.path.join(args.outdir, "runs_genomes")
        logger.info(
            "Writing table of pyani runs, with associated genomes to %s.*", outfname
        )
        # We query on the LabelMembership table, then join out to Run, Genome
        # and Label, so that we avoid multiple left joins
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
            .join(rungenome, Genome)
            .join(
                Label,
                and_(Genome.genome_id == Label.genome_id, Run.run_id == Label.run_id),
            )
            .order_by(Run.run_id, Genome.genome_id)
            .statement
        )
        data = pd.read_sql(statement, session.bind)
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
        data.columns = headers
        pyani_report.write_dbtable(data, outfname, formats)

    # Report table of all runs in which a genome is involved
    if args.show_genomes_runs:
        outfname = os.path.join(args.outdir, "genomes_runs")
        logger.info(
            "Writing table of genomes, with associated pyani runs to %s.*", outfname
        )
        # We query on the LabelMembership table, then join out to Run, Genome
        # and Label, so that we avoid multiple left joins
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
            .join(rungenome, Run)
            .join(
                Label,
                and_(Genome.genome_id == Label.genome_id, Run.run_id == Label.run_id),
            )
            .order_by(Genome.genome_id, Run.run_id)
            .statement
        )
        data = pd.read_sql(statement, session.bind)
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
        data.columns = headers
        pyani_report.write_dbtable(data, outfname, formats)

    # Report table of comparison results for the indicated runs
    if args.run_results:
        outfstem = os.path.join(args.outdir, "results")
        run_ids = [run_id.strip() for run_id in args.run_results.split(",")]
        logger.info("Attempting to write results tables for runs: %s", run_ids)
        for run_id in run_ids:
            logger.info("Processing run %d", run_id)
            outfname = "_".join([outfstem, str(run_id)])
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
            data = pd.read_sql(statement, session.bind)
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
            data.columns = headers
            pyani_report.write_dbtable(data, outfname, formats)

    # Report matrices of comparison results for the indicated runs
    # For ANIm, all results other than coverage are symmetric matrices,
    # so we only get results in the forward direction.
    if args.run_matrices:
        outfstem = os.path.join(args.outdir, "matrix")
        run_ids = [run_id.strip() for run_id in args.run_matrices.split(",")]
        logger.info("Writing result matrices for runs: %s", run_ids)
        for run_id in run_ids:
            logger.info("Extracting matrices for run %s", run_id)
            run = session.query(Run).filter(Run.run_id == run_id).first()
            #  Make dictionary of labels for each genome, keyed by genome ID
            results = (
                session.query(Genome.genome_id, Label.label)
                .join(rungenome, Run)
                .join(
                    Label,
                    and_(
                        Genome.genome_id == Label.genome_id, Run.run_id == Label.run_id
                    ),
                )
                .filter(Run.run_id == run_id)
                .all()
            )
            label_dict = {_.genome_id: _.label for _ in results}
            for matname, data, args in [
                ("identity", run.df_identity, {"colour_num": 0.95}),
                ("coverage", run.df_coverage, {"colour_num": 0.95}),
                ("aln_lengths", run.df_alnlength, {}),
                ("sim_errors", run.df_simerrors, {}),
                ("hadamard", run.df_hadamard, {}),
            ]:
                logger.info("Writing %s results", matname)
                outfname = "_".join([outfstem, matname, str(run_id)])
                matrix = pd.read_json(data)
                # Matrix rows and columns are labelled if there is a label
                # dictionary, and take the dataframe index otherwise
                matrix.columns = [
                    "{}:{}".format(label_dict.get(col, col), col)
                    for col in matrix.columns
                ]
                matrix.index = [
                    "{}:{}".format(label_dict.get(idx, idx), idx)
                    for idx in matrix.index
                ]
                pyani_report.write_dbtable(
                    matrix, outfname, formats, show_index=True, **args
                )
