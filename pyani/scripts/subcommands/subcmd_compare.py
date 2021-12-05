from argparse import Namespace
from pathlib import Path
from typing import Any, NamedTuple, Dict, List, Set
from itertools import combinations, permutations

import logging
import multiprocessing

import pandas as pd
import matplotlib.pyplot as plt
from itertools import product

from pyani import pyani_config, pyani_graphics
from pyani.pyani_tools import termcolor, MatrixData
from pyani.pyani_orm import (
    PyaniORMException,
    get_session,
    get_matrix_classes_for_run,
    get_matrix_labels_for_run,
    Comparison,
    Run,
    rungenome,
)
from pyani.pyani_graphics.sns import get_clustermap, get_colorbar
import sys

# Distribution dictionary of matrix graphics methods
GMETHODS = {"mpl": pyani_graphics.mpl.heatmap, "seaborn": pyani_graphics.sns.heatmap}
# Distribution dictionary of distribution graphics methods
DISTMETHODS = {
    "mpl": pyani_graphics.mpl.distribution,
    "seaborn": pyani_graphics.sns.distribution,
}
# Dictionary of scatter graphics methods
SMETHODS = {"seaborn": pyani_graphics.sns.scatter, "mpl": pyani_graphics.mpl.scatter}


# Convenience struct for run data
class RunData(NamedTuple):
    run_id: int
    method: str
    cmdline: str
    genomes: Set
    identity: MatrixData
    coverage: MatrixData
    aln_length: MatrixData
    sim_errors: MatrixData
    hadamard: MatrixData


# Convenience struct for a set of matrices
class RunMatrices(NamedTuple):
    identity: MatrixData
    coverage: MatrixData
    aln_length: MatrixData
    sim_errors: MatrixData
    hadamard: MatrixData


def subcmd_compare(args: Namespace):
    """Performs a comparison between two completed pyani runs. Runs
    may differ in method used, parameter values, or both.

    :param args:  Namespace, command-line arguments

    Produces a series of scatter, heatpmap, and distribution plots, and a
    summary report.

    """
    # Setup
    # Create logger
    logger = logging.getLogger(__name__)

    # Get connection to existing database. This may or may not have data
    logger.debug("Connecting to database %s", args.dbpath)
    try:
        session = get_session(args.dbpath)
    except PyaniORMException:
        logger.error(
            "Could not connect to database %s (exiting)", args.dbpath, exc_info=True
        )
        raise SystemExit(1)

    # Get run ids
    # run_a, run_b = int(args.run_a), int(args.run_b)
    comparisons = [_ for _ in product(args.ref_ids, args.run_ids)]

    # Announce the analysis
    logger.info(
        "Running %d comparisons between refs: %s and runs: %s",
        len(comparisons),
        args.ref_ids,
        args.run_ids,
    )
    logger.debug("Comparisons to run: %s", comparisons)
    # Parse output formats
    outfmts = args.formats
    logger.debug("Requested output formats: %s", outfmts)

    # Get information on the runs
    # runs = [run_a, run_b]
    runs = args.ref_ids + args.run_ids
    run_dict = {}
    logger.debug("Getting run data from database %s", args.dbpath)
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
            "At least one specified run not found in the database %s (exiting)",
            args.dbpath,
        )
        raise SystemExit(1)

    # Get sets of genomes for each run and parse run data
    for run in run_data:
        genome_query = session.query(rungenome).filter_by(run_id=run.run_id)
        genome_set = set(gen for (gen, run) in genome_query)
        run_dict.update({f"{run.run_id}": parse_data(run, genome_set)})

    # Loop over pairs of runs
    for ref, query in comparisons:
        ref, query = run_dict[str(ref)], run_dict[str(query)]
        # Find common genomes
        common = ref.genomes & query.genomes

        if not common:
            logger.error("No genomes in common between %s and %s", ref, query)
            raise SystemExit(1)
        logger.debug(
            "\t...%d genomes in common between %s and %s.", len(common), ref, query
        )

        # Get matrix labels, classes
        labels = get_labels(session, ref, query)
        classes = get_classes(session, ref, query)

        # Subset matrices based on common genomes
        sub_ref = subset_matrix(common, ref)
        sub_query = subset_matrix(common, query)

        # Generate dataframes of differences for each measure
        difference_matrices = get_difference_matrices(sub_ref, sub_query)

        # Tetra doesn't report all of the same things

        # Create worker pool and empty command list
        pool = multiprocessing.Pool(processes=args.workers)
        plotting_commands = []

        for A, B in zip(sub_ref, sub_query):
            # Plot scatter plots for each score
            logger.info(f"{A.name}, {B.name}")
            logger.info(f"{ref.run_id}, {query.run_id}")
            plotting_commands.append(
                (get_scatter, [ref.run_id, query.run_id, A, B, outfmts, args])
            )
            # Create Bland-Altman plots

        # Send dataframes for heatmaps, scatterplots
        for matdata in difference_matrices.values():
            # Write heatmap for each results matrix
            plotting_commands.append(
                (
                    get_heatmap,
                    [ref.run_id, query.run_id, matdata, labels, classes, outfmts, args],
                )
            )
            # Plot distributions of differences to look at normality
            plotting_commands.append(
                (get_distribution, [ref.run_id, query.run_id, matdata, outfmts, args])
            )

        # Run the plotting commands
        for func, args in plotting_commands:
            # pool.apply_async(func, args, {})
            logger.info(f"{func}")

        # Close worker pool
        pool.close()
        pool.join()

        # Generate summary report


def get_metadata(session: Any, run_id: int) -> RunData:
    """Get metadata for a run in the database.

    :param session:  live SQLAlchemy session of pyani database
    :param run_id:  unique identifier for the run in question

    """
    return session.query(Run.run_id, Run.method, Run.cmdline).filter_by(run_id=run_id)


def parse_data(run, genome_set) -> RunData:
    """Get metadata for a run in the database.

    :param run:  data pertaining to one run in the database
    :param genome_set:  a set of genome_ids included in the run

    """
    return RunData(
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


def get_labels(session, reference, query):
    """Grab labels information for the runs.

    :param session:  live SQLAlchemy session of pyani database
    :param reference:  reference run
    :param query:  query run

    """
    labels = {}
    for run in (reference, query):
        labels.update(get_matrix_labels_for_run(session, run.run_id))
    return labels


def get_classes(session, reference, query):
    """Grab class information for the runs.

    :param session:  live SQLAlchemy session of pyani database
    :param reference:  reference run
    :param query:  query run

    """
    classes = {}
    for run in (reference, query):
        classes.update(get_matrix_classes_for_run(session, run.run_id))
    return classes


def subset_matrix(common, run):
    """Subsets a score matrix based on a set of indices."""
    return RunMatrices(
        MatrixData("identity", run.identity.data.loc[common, common], {}),
        MatrixData("coverage", run.coverage.data.loc[common, common], {}),
        MatrixData("aln_length", run.aln_length.data.loc[common, common], {}),
        MatrixData("sim_errors", run.sim_errors.data.loc[common, common], {}),
        MatrixData("hadamard", run.hadamard.data.loc[common, common], {}),
    )


def get_difference_matrices(reference, query) -> Dict:
    """Calculate difference matrices between two sets of runs.

    :param reference:  a set of score matrices to be used as the references
    :param query:  a set of score matrices to subtract from the references

    Returns a dictionary object with absolute and relative difference matrices
    as values, keyed by e.g., 'identity_diffs', 'identity_absdiffs'.
    """
    difference_matrices = {}
    difference_matrices.update(
        {
            f"{a.name}_diffs": MatrixData(f"{a.name}_diffs", a.data - b.data, {})
            for a, b in zip(reference, query)
        }
    )
    difference_matrices.update(
        {
            f"{a.name}_absdiffs": MatrixData(
                f"{a.name}_absdiffs", abs(a.data - b.data), {}
            )
            for a, b in zip(reference, query)
        }
    )
    return difference_matrices


def get_heatmap(
    run_a: int,
    run_b: int,
    matdata: MatrixData,
    result_labels: Dict,
    result_classes: Dict,
    outfmts: List[str],
    args: Namespace,
):
    """Write a single heatmap for a pyani run.

    :param run_id:  int, run_id for this run
    :param matdata:  MatrixData object for this heatmap
    :param result_labels:  dict of result labels
    :param result_classes: dict of result classes
    :param args:  Namespace for command-line arguments
    :param outfmts:  list of output formats for files
    """
    cmap = ("BuRd", matdata.data.values.min(), matdata.data.values.max())

    logger = logging.getLogger(__name__)

    logger.info("Creating %s matrix heatmaps", matdata.name)
    logger.info(
        "Cmap min: %s, cmap max: %s",
        matdata.data.values.min(),
        matdata.data.values.max(),
    )

    for fmt in outfmts:
        outfname = (
            Path(args.outdir) / f"compare_{matdata.name}_run{run_a}_run{run_b}.{fmt}"
        )
        logger.debug("\tWriting graphics to %s", outfname)
        params = pyani_graphics.Params(cmap, result_labels, result_classes)

        # Draw heatmap
        # get_clustermap(matdata.data, params, "")
        GMETHODS["seaborn"](
            matdata.data,
            outfname,
            title=f"compare_{matdata.name}_run{run_a}_run{run_b}",
            params=params,
        )

    # Be tidy with matplotlib caches
    plt.close("all")


def get_distribution(
    run_a: int,
    run_b: int,
    matdata: MatrixData,
    outfmts: List[str],
    args: Namespace,
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
        outfname = (
            Path(args.outdir)
            / f"distribution_{matdata.name}_run{run_a}_run{run_b}.{fmt}"
        )
        logger.debug("\tWriting graphics to %s", outfname)
        DISTMETHODS[args.method](
            matdata.data,
            outfname,
            matdata.name,
            title=f"compare_{matdata.name}_run{run_a}_run{run_b}",
        )


def get_scatter(
    run_a: int,
    run_b: int,
    matdata1: MatrixData,
    matdata2: MatrixData,
    outfmts: List[str],
    args: Namespace,
) -> None:
    """Write a single scatterplot for a pyani run.

    :param run_a:  int, run_id for the reference
    :param run_b:  int, run_id for the query
    :param matdata1:  MatrixData object for this scatterplot
    :param matdata2:  MatrixData object for this scatterplot
    :param result_labels:  dict of result labels
    :param result_classes: dict of result classes
    :param args:  Namespace for command-line arguments
    :param outfmts:  list of output formats for files
    """
    logger = logging.getLogger(__name__)
    sys.stderr.write("Grapevine Fires\n")
    sys.stderr.write(f"Writing {matdata1.name} vs {matdata2.name} scatterplot\n")
    extreme = max(abs(matdata1.data.values.min()), abs(matdata1.data.values.max()))
    cmap = ("BuRd", extreme * -1, extreme)
    for fmt in outfmts:
        outfname = (
            Path(args.outdir)
            / f"scatter_{matdata1.name}_run{run_a}_vs_{matdata2.name}_run{run_b}.{fmt}"
        )
        sys.stderr.write(f"{outfname}\n")
        logger.debug("\tWriting graphics to %s", outfname)
        params = pyani_graphics.Params(cmap, {}, {})
        # Draw scatterplot
        sys.stderr.write(f"{SMETHODS[args.method]}")
        SMETHODS[args.method](
            matdata1.data,
            matdata2.data,
            outfname,
            f"{matdata1.name}_{run_a}",
            f"{matdata2.name}_{run_b}",
            title=f"{matdata1.name.title()} run {run_a} vs {matdata2.name.title()} run {run_b}",
            params=params,
        )

        # Be tidy with matplotlib caches
        # plt.close("all")
