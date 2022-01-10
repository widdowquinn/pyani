import logging
import os
import multiprocessing

from argparse import Namespace
from pathlib import Path
from typing import Dict, List


from pyani import pyani_config, pyani_orm, pyani_graphics
from pyani.pyani_tools import termcolor, MatrixData

TREEMETHODS = {}
# TREEMETHODS = {"ete3": pyani_graphics.tree.tree}

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
        write_trees(run_id, session, outfmts, args)

        if NEWICKS:
            write_newicks(args, run_id)
        NEWICKS.clear()

    return 0


def write_trees(
    run_id: int,
    matdata: MatrixData,
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

    logger.info("Writing tree plot for %s matrix", matdata.name)

    for fmt in outfmts:
        outfname = Path(args.outdir) / f"distribution_{matdata.name}_run{run_id}.{fmt}"

        TREEMETHODS[args.method](
            matdata.data,
            outfname,
            matdata.name,
            title=f"matrix_{matdata.name}_run{run_id}",
        )


def write_newicks(args: Namespace, run_id):
    # If Newick strings were generated, write them out.
    newick_file = Path(args.outdir) / f"newicks_run{run_id}.nw"
    with open(newick_file, "w") as nfh:
        for name, nw in NEWICKS.items():
            nfh.write(f"{name}\t{nw}\n")
