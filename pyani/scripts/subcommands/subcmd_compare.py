import pandas as pd
from argparse import Namespace
import logging
import matplotlib.pyplot as plt
from pyani import pyani_config, pyani_graphics
from pyani.pyani_tools import termcolor, MatrixData
from pyani.pyani_orm import (
    PyaniORMException,
    get_session,
    get_genome_pair_dict,
    get_matrix_classes_for_run,
    get_matrix_labels_for_run,
    filter_uncommon_comparisons,
    get_df_of_scores,
    Comparison,
    Run,
    rungenome,
)
from typing import Any, NamedTuple, Dict, List, Set
from itertools import combinations, permutations
from pyani.pyani_graphics.sns import get_clustermap, get_colorbar

# Distribution dictionary of matrix graphics methods
GMETHODS = {"mpl": pyani_graphics.mpl.heatmap, "seaborn": pyani_graphics.sns.heatmap}
# Distribution dictionary of distribution graphics methods
DISTMETHODS = {
    "mpl": pyani_graphics.mpl.distribution,
    "seaborn": pyani_graphics.sns.distribution,
}


class RunMetadata(NamedTuple):
    run_id: int
    method: str
    cmdline: str
    genomes: Set
    identity: MatrixData
    coverage: MatrixData
    aln_length: MatrixData
    sim_errors: MatrixData
    hadamard: MatrixData


class SubData(NamedTuple):
    identity: MatrixData
    coverage: MatrixData
    aln_length: MatrixData
    sim_errors: MatrixData
    hadamard: MatrixData


def subcmd_compare(args: Namespace):
    # Setup
    # Create logger
    logger = logging.getLogger(__name__)

    # Get run ids
    run_a, run_b = int(args.run_a), int(args.run_b)

    # Announce the analysis
    logger.info(termcolor(f"Comparing runs {run_a} and {run_b}", bold=True))

    # ¶ No need to get version info here

    # Get connection to existing database. This may or may not have data
    logger.debug(f"Connecting to database {args.dbpath}")
    try:
        session = get_session(args.dbpath)
    except PyaniORMException:
        logger.error(
            f"Could not connect to database {args.dbpath} (exiting)", exc_info=True
        )
        raise SystemExit(1)

    # Get information on the runs
    run_data = {}

    logger.debug(f"Getting run metadata from database {args.dbpath}")
    for run in run_a, run_b:
        try:
            run_data.update(
                {f"run_{run}": RunMetadata(*x) for x in get_metadata(session, run)}
            )
        except PyaniORMException:
            logger.error(f"Run {run} not found in the database {args.dbpath} (exiting)")
            raise SystemExit(1)
    for value in run_data.values():
        logger.info(f"Run_{value.run_id} command line: {value.cmdline}")

    # Get comparisons in run A
    logger.debug(f"Getting comparisons for run {run_a} from database {args.dbpath}")
    try:
        run_A_comps = get_genome_pair_dict(session, run_a)
    except PyaniORMException:
        logger.error(
            f"Unable to retrieve comparison information for run {run_a} (exiting)"
        )
        raise SystemExit(1)
    logger.debug(
        f"\t...run {run_a} is present in the database with {len(run_A_comps)} comparisons"
    )

    # Get comparisons in run B
    logger.debug(f"Getting comparisons for run {run_b} from database {args.dbpath}")
    try:
        run_B_comps = get_genome_pair_dict(session, run_b)
    except PyaniORMException:
        logger.error(
            f"Unable to retrieve comparison information for run {run_b} (exiting)"
        )
        raise SystemExit(1)
    logger.debug(
        f"\t...run {run_b} is present in the database with {len(run_B_comps)} comparisons"
    )

    # Get overlapping pairs of compared genomes
    # these are the only ones we can compare
    run_A_common, run_B_common = filter_uncommon_comparisons(
        session, run_A_comps, run_B_comps
    )
    if not run_A_common:
        logger.error(f"No comparisons in common between {run_a} and {run_b}")
        raise SystemExit(1)
    logger.debug(
        f"\t...{len(run_A_common)}, {len(run_B_common)} comparisons in common between {run_a} and {run_b}."
    )

    # Retrieve comparison
    # overlap_comps = session.query(pyani_orm.Comparison).filter(pyani_orm.Comparison.comparison_id.in_(overlap))

    # Get dataframes of scores for each run
    # List of columns to extract
    scores = [
        Comparison.query_id,
        Comparison.subject_id,
        Comparison.aln_length,
        Comparison.sim_errs,
        Comparison.identity,
        Comparison.cov_query,
        Comparison.cov_subject,
    ]

    run_A_scores = get_df_of_scores(session, run_A_common, scores)
    run_B_scores = get_df_of_scores(session, run_B_common, scores)

    # List of columns; needs to vary if Tetra is used
    columns = [
        "aln_length",
        "sim_errs",
        "identity",
        "cov_query",
        "cov_subject",
        "hadamard",
    ]

    # Add hadamard columns
    run_A_scores["hadamard"] = run_A_scores["identity"] * run_A_scores["cov_query"]
    run_B_scores["hadamard"] = run_B_scores["identity"] * run_B_scores["cov_query"]
    logger.debug(f"\t...run_A_scores has dimensions: {run_A_scores.shape}.")
    logger.debug(f"\t...run_A_scores.head(): \n{run_A_scores.head()}")
    logger.debug(f"\t...run_B_scores has dimensions: {run_B_scores.shape}.")
    logger.debug(f"\t...run_B_scores.head(): \n{run_B_scores.head()}")

    # Merge dataframes; this ensures all columns are exactly the same length
    merged_scores = run_A_scores.merge(
        run_B_scores, how="inner", on=["query_id", "subject_id"], suffixes=["_A", "_B"]
    )
    logger.debug(f"\t...merged_scores has dimensions: {merged_scores.shape}.")
    logger.debug(f"\t...merged_scores.head(): \n{merged_scores.head()}")

    # Generate dataframes of differences for each measure
    diffs = pd.DataFrame(
        index=merged_scores.index, columns=[col + "_diff" for col in columns]
    )

    abs_diffs = pd.DataFrame(
        index=merged_scores.index, columns=[col + "_abs" for col in columns]
    )

    # Tetra doesn't report all of the same things
    for col in run_A_scores.columns:
        diffs[f"{col}_diff"] = merged_scores[f"{col}_A"] - merged_scores[f"{col}_B"]
        # abs_diffs[f"{col}_abs"] = [abs(_) for _ in diffs[f"{col}_diff"]]
        abs_diffs[f"{col}_abs"] = abs(
            merged_scores[f"{col}_A"] - merged_scores[f"{col}_B"]
        )

    # Send dataframes for heatmaps, scatterplots
    # Heatmaps will use... something

    # Plot distributions of differences to look at normality

    # Generate summary report

    # raise NotImplementedError


def get_metadata(session: Any, run_id: int) -> RunMetadata:
    """Get metadata for a run in the database.

    :param session:  live SQLAlchemy session of pyani database
    :param run_id:  unique identifier for the run in question

    """
    return session.query(Run.run_id, Run.method, Run.cmdline).filter_by(run_id=run_id)
