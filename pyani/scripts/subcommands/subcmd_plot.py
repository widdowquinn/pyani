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
"""Provides the plot subcommand for pyani."""

import logging
import os
import multiprocessing

from argparse import Namespace
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import pandas as pd

from pyani import pyani_config, pyani_orm, pyani_graphics
from pyani.pyani_tools import termcolor, MatrixData

# Distribution dictionary of matrix graphics methods
GMETHODS = {"mpl": pyani_graphics.mpl.heatmap, "seaborn": pyani_graphics.sns.heatmap}
SMETHODS = {"mpl": pyani_graphics.mpl.scatter, "seaborn": pyani_graphics.sns.scatter}
# Distribution dictionary of distribution graphics methods
DISTMETHODS = {
    "mpl": pyani_graphics.mpl.distribution,
    "seaborn": pyani_graphics.sns.distribution,
}


def subcmd_plot(args: Namespace) -> int:
    """Produce graphical output for an analysis.

    :param args:  Namespace of command-line arguments

    This is graphical output for representing the ANI analysis results, and
    takes the form of a heatmap, or heatmap with dendrogram.
    """
    logger = logging.getLogger(__name__)

    # Announce what's going on to the user
    logger.info(termcolor("Generating graphical output for analyses", "red"))
    logger.info("Writing output to: %s", args.outdir)
    os.makedirs(args.outdir, exist_ok=True)
    logger.info("Rendering method: %s", args.method)
    logger.info("Rendering method: %s", args)

    # Connect to database session
    logger.debug("Activating session for database: %s", args.dbpath)
    session = pyani_orm.get_session(args.dbpath)

    # Parse output formats
    outfmts = args.formats  # .formats.split(",")
    logger.debug("Requested output formats: %s", outfmts)
    logger.debug("Type of formats variable: %s", type(outfmts))

    # Work on each run:
    run_ids = [int(run) for run in args.run_ids]
    logger.debug("Generating graphics for runs: %s", run_ids)
    for run_id in run_ids:
        write_run_plots(run_id, session, outfmts, args)

    return 0


def write_run_plots(run_id: int, session, outfmts: List[str], args: Namespace) -> None:
    """Write all heatmaps for a specified run to file.

    :param run_id:  int, run identifier in database session
    :param session:  Session, active SQLite session
    :param outfmts:  list of output format types
    :param args:  Namespace, command line arguments
    """
    logger = logging.getLogger(__name__)

    # Get results matrices for the run
    logger.debug("Retrieving results matrices for run %s", run_id)
    results = (
        session.query(pyani_orm.Run).filter(pyani_orm.Run.run_id == run_id).first()
    )
    result_label_dict = pyani_orm.get_matrix_labels_for_run(session, run_id)
    result_class_dict = pyani_orm.get_matrix_classes_for_run(session, run_id)
    logger.debug(
        f"Have {len(result_label_dict)} labels and {len(result_class_dict)} classes"
    )

    # Write heatmap and distribution plot for each results matrix

    # Create worker pool and empty command list
    pool = multiprocessing.Pool(processes=args.workers)
    plotting_commands = []

    # Build and collect the plotting commands
    for matdata in [
        MatrixData(*_)
        for _ in [
            ("identity", pd.read_json(results.df_identity), {}),
            ("coverage", pd.read_json(results.df_coverage), {}),
            ("aln_lengths", pd.read_json(results.df_alnlength), {}),
            ("sim_errors", pd.read_json(results.df_simerrors), {}),
            ("hadamard", pd.read_json(results.df_hadamard), {}),
        ]
    ]:
        plotting_commands.append(
            (
                write_heatmap,
                [run_id, matdata, result_label_dict, result_class_dict, outfmts, args],
            )
        )

        plotting_commands.append((write_distribution, [run_id, matdata, outfmts, args]))

    id_matrix = MatrixData("identity", pd.read_json(results.df_identity), {})
    cov_matrix = MatrixData("coverage", pd.read_json(results.df_coverage), {})
    plotting_commands.append(
        (
            write_scatter,
            [
                run_id,
                id_matrix,
                cov_matrix,
                result_label_dict,
                result_class_dict,
                outfmts,
                args,
            ],
        )
    )

    # Run the plotting commands
    logger.debug("Running plotting commands")
    for func, options in plotting_commands:
        logger.debug("Running %s with options %s", func, options)
        pool.apply_async(func, args=options)

    # Close worker pool
    pool.close()
    pool.join()


def write_distribution(
    run_id: int, matdata: MatrixData, outfmts: List[str], args: Namespace
) -> None:
    """Write distribution plots for each matrix type.

    :param run_id:  int, run_id for this run
    :param matdata:  MatrixData object for this distribution plot
    :param args:  Namespace for command-line arguments
    :param outfmts:  list of output formats for files
    """
    logger = logging.getLogger(__name__)

    logger.info("Writing distribution plot for %s matrix", matdata.name)
    for fmt in outfmts:
        outfname = Path(args.outdir) / f"distribution_{matdata.name}_run{run_id}.{fmt}"
        logger.debug("\tWriting graphics to %s", outfname)
        DISTMETHODS[args.method[0]](
            matdata.data,
            outfname,
            matdata.name,
            title=f"matrix_{matdata.name}_run{run_id}",
        )
    # Be tidy with matplotlib caches
    plt.close("all")


def write_heatmap(
    run_id: int,
    matdata: MatrixData,
    result_labels: Dict,
    result_classes: Dict,
    outfmts: List[str],
    args: Namespace,
) -> None:
    """Write a single heatmap for a pyani run.

    :param run_id:  int, run_id for this run
    :param matdata:  MatrixData object for this heatmap
    :param result_labels:  dict of result labels
    :param result_classes: dict of result classes
    :param args:  Namespace for command-line arguments
    :param outfmts:  list of output formats for files
    """
    logger = logging.getLogger(__name__)

    logger.info("Writing %s matrix heatmaps", matdata.name)
    cmap = pyani_config.get_colormap(matdata.data, matdata.name)
    for fmt in outfmts:
        outfname = Path(args.outdir) / f"matrix_{matdata.name}_run{run_id}.{fmt}"
        logger.debug("\tWriting graphics to %s", outfname)
        params = pyani_graphics.Params(cmap, result_labels, result_classes)
        # Draw heatmap
        GMETHODS[args.method[0]](
            matdata.data,
            outfname,
            title=f"matrix_{matdata.name}_run{run_id}",
            params=params,
        )

    # Be tidy with matplotlib caches
    plt.close("all")


def write_scatter(
    run_id: int,
    matdata1: MatrixData,
    matdata2: MatrixData,
    result_labels: Dict,
    result_classes: Dict,
    outfmts: List[str],
    args: Namespace,
) -> None:
    """Write a single scatterplot for a pyani run.

    :param run_id:  int, run_id for this run
    :param matdata1:  MatrixData object for this scatterplot
    :param matdata2:  MatrixData object for this scatterplot
    :param result_labels:  dict of result labels
    :param result_classes: dict of result classes
    :param args:  Namespace for command-line arguments
    :param outfmts:  list of output formats for files
    """
    logger = logging.getLogger(__name__)

    logger.info("Writing %s vs %s scatterplot", matdata1.name, matdata2.name)
    cmap = pyani_config.get_colormap(matdata1.data, matdata1.name)
    for fmt in outfmts:
        outfname = (
            Path(args.outdir)
            / f"scatter_{matdata1.name}_vs_{matdata2.name}_run{run_id}.{fmt}"
        )
        logger.debug("\tWriting graphics to %s", outfname)

        params = pyani_graphics.Params(cmap, result_labels, result_classes)
        # Draw scatterplot
        SMETHODS[args.method[0]](
            matdata1.data,
            matdata2.data,
            outfname,
            matdata1.name,
            matdata2.name,
            title=f"{matdata1.name.title()} vs {matdata2.name.title()}",
            params=params,
        )

        # Be tidy with matplotlib caches
        plt.close("all")
