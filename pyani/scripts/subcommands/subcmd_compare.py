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
    runs = [run_a, run_b]
    run_dict = {}
    logger.debug(f"Getting run data from database {args.dbpath}")
    try:
        run_data = session.query(
            Run.run_id,
            Run.method,
            Run.cmdline,
            Run.df_identity,
            Run.df_coverage,
            Run.df_alnlength,
            Run.df_simerrors,
            Run.df_hadamard,
        ).filter(Run.run_id.in_(runs))
    except PyaniORMException:
        logger.error(
            f"At least one specified run not found in the database {args.dbpath} (exiting)"
        )
        raise SystemExit(1)

    # Get sets of genomes for each run and parse run data
    for run in run_data:
        genome_query = session.query(rungenome).filter_by(run_id=run.run_id)
        genome_set = set(gen for (gen, run) in genome_query)
        run_dict.update({f"{run.run_id}": parse_data(run, genome_set)})

    # Loop over pairs of runs
    for ref, query in permutations(runs, 2):
        ref, query = run_dict[str(ref)], run_dict[str(query)]
        # Find common genomes
        common = ref.genomes & query.genomes

        if not common:
            logger.error(f"No genomes in common between {ref} and {query}")
            raise SystemExit(1)
        logger.debug(f"\t...{len(common)} genomes in common between {ref} and {query}.")

        # Subset matrices based on common genomes
        sub_ref = subset_matrix(common, ref)
        sub_query = subset_matrix(common, query)

        # Generate dataframes of differences for each measure
        diffs = {
            a.name: MatrixData(a.name, a.data - b.data, {})
            for a, b in zip(sub_ref, sub_query)
        }
        abs_diffs = {
            a.name: MatrixData(a.name, abs(a.data - b.data), {})
            for a, b in zip(sub_ref, sub_query)
        }

        logger.debug(f"Difference matrix:\n {diffs}")
        logger.debug(f"Absolute difference matrix:\n {abs_diffs}")

        # Tetra doesn't report all of the same things

        # Get matrix labels, classes
        labels, classes = {}, {}
        logger.info(f"{ref.run_id}")
        labels.update(get_matrix_labels_for_run(session, ref.run_id))
        labels.update(get_matrix_labels_for_run(session, query.run_id))
        classes.update(get_matrix_classes_for_run(session, ref.run_id))
        classes.update(get_matrix_classes_for_run(session, query.run_id))

        logger.info(f"Labels: {labels}")

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


def parse_data(run, genome_set) -> RunMetadata:
    """Get metadata for a run in the database.

    :param run:  data pertaining to one run in the database
    :param genome_set:  a set of genome_ids included in the run

    """
    return RunMetadata(
        run.run_id,
        run.method,
        run.cmdline,
        genome_set,
        MatrixData("identity", pd.read_json(run.df_identity), {}),
        MatrixData("coverage", pd.read_json(run.df_coverage), {}),
        MatrixData("aln_length", pd.read_json(run.df_alnlength), {}),
        MatrixData("sim_errors", pd.read_json(run.df_simerrors), {}),
        MatrixData("hadamard", pd.read_json(run.df_hadamard), {}),
    )


def subset_matrix(common, run):
    """Subsets a score matrix based on a set of indices."""
    return SubData(
        MatrixData("identity", run.identity.data.loc[common, common], {}),
        MatrixData("coverage", run.coverage.data.loc[common, common], {}),
        MatrixData("aln_length", run.aln_length.data.loc[common, common], {}),
        MatrixData("sim_errors", run.sim_errors.data.loc[common, common], {}),
        MatrixData("hadamard", run.hadamard.data.loc[common, common], {}),
    )
