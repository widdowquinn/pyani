import logging
import os
import sys
import multiprocessing

from argparse import Namespace
from pathlib import Path
from typing import Dict, List
import pandas as pd

from pyani import pyani_config, pyani_orm, pyani_graphics
from pyani.pyani_tools import termcolor, MatrixData

# TREEMETHODS = {}
TREEMETHODS = {"ete3": pyani_graphics.tree.tree}

NEWICKS = {}


def subcmd_tree(args: Namespace) -> int:
    """Produce tree output for an analysis.

    :param args:  Namespace of command-line arguments

    This is graphical output for representing the ANI analysis results, and
    takes the form of a tree, or dendrogram.
    """
    logger = logging.getLogger(__name__)

    # Announce what's going on to the user
    logger.info(termcolor("Generating tree output for analyses", "red"))
    logger.info("Writing output to: %s", args.outdir)
    os.makedirs(args.outdir, exist_ok=True)
    logger.info("Rendering method: %s", args.method)

    # Connect to database session
    logger.debug("Activating session for database: %s", args.dbpath)
    session = pyani_orm.get_session(args.dbpath)

    # Parse output formats
    outfmts = args.formats
    logger.debug("Requested output formats: %s", outfmts)
    logger.debug("Type of formats variable: %s", type(outfmts))

    # Work on each run:
    run_ids = args.run_ids
    logger.debug("Generating trees for runs: %s", run_ids)
    for run_id in run_ids:
        write_run_trees(run_id, session, outfmts, args)

        if NEWICKS:
            write_newicks(args, run_id)
        NEWICKS.clear()

    return 0


def write_run_trees(
    run_id: int,
    session,
    outfmts: List[str],
    args: Namespace,
) -> None:
    """Write tree plots for each matrix type.

    :param run_id:  int, run_id for this run
    :param matdata:  MatrixData object for this distribution plot
    :param args:  Namespace for command-line arguments
    :param outfmts:  list of output formats for files
    """
    logger = logging.getLogger(__name__)
    logger.debug("Retrieving results matrices for run %s", run_id)

    results = (
        session.query(pyani_orm.Run).filter(pyani_orm.Run.run_id == run_id).first()
    )
    result_label_dict = pyani_orm.get_matrix_labels_for_run(session, run_id)
    result_class_dict = pyani_orm.get_matrix_classes_for_run(session, run_id)
    logger.debug(
        f"Have {len(result_label_dict)} labels and {len(result_class_dict)} classes"
    )

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
        if _[0] in args.trees
    ]:
        logger.info("Writing tree plot for %s matrix", matdata.name)
        plotting_commands.append(
            (
                write_tree,
                [run_id, matdata, result_label_dict, result_class_dict, outfmts, args],
            )
        )

    sys.stdout.write(str(plotting_commands))

    # Run the plotting commands
    for func, options in plotting_commands:
        result = pool.apply_async(func, options, {}, callback=logger.debug)
        result.get()

    # Close worker pool
    pool.close()
    pool.join()


def write_tree(
    run_id: int,
    matdata: MatrixData,
    result_labels: Dict,
    result_classes: Dict,
    outfmts: List[str],
    args: Namespace,
) -> None:
    """Write a single tree for a pyani run.

    :param run_id:  int, run_id for this run
    :param matdata:  MatrixData object for this heatmap
    :param result_labels:  dict of result labels
    :param result_classes: dict of result classes
    :param args:  Namespace for command-line arguments
    :param outfmts:  list of output formats for files
    """
    # logger = logging.getLogger(__name__)
    cmap = pyani_config.get_colormap(matdata.data, matdata.name)

    for fmt in outfmts:
        outfname = Path(args.outdir) / f"distribution_{matdata.name}_run{run_id}.{fmt}"

        params = pyani_graphics.Params(cmap, result_labels, result_classes)

        TREEMETHODS[args.method](
            matdata.data,
            outfname,
            title=f"matrix_{matdata.name}_run{run_id}",
            params=params,
            format=fmt,
            args=args,
        )


def write_newicks(args: Namespace, run_id):
    # If Newick strings were generated, write them out.
    newick_file = Path(args.outdir) / f"newicks_run{run_id}.nw"
    with open(newick_file, "w") as nfh:
        for name, nw in NEWICKS.items():
            nfh.write(f"{name}\t{nw}\n")
