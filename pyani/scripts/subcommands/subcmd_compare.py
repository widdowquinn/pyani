import os
import sys

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

# Distribution dictionary of matrix graphics methods
GMETHODS = {"mpl": pyani_graphics.mpl.heatmap, "seaborn": pyani_graphics.sns.heatmap}
# Distribution dictionary of distribution graphics methods
DISTMETHODS = {
    "mpl": pyani_graphics.mpl.distribution,
    "seaborn": pyani_graphics.sns.distribution,
}
# Dictionary of scatter graphics methods
SMETHODS = {"seaborn": pyani_graphics.sns.scatter, "mpl": pyani_graphics.mpl.scatter}

# Dictionary of Bland-ALtman graphics methods
BMETHODS = {
    "seaborn": pyani_graphics.sns.bland_altman,
    "mpl": pyani_graphics.mpl.bland_altman,
}


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

    # Create subdirectories for output
    outsubdir = Path(args.outdir) / "compare"
    os.makedirs(outsubdir, exist_ok=True)

    logger.info("Writing output to: %s", outsubdir)

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

        outsubdir = Path(outsubdir) / f"ref_{ref}_vs_query_{query}"
        os.makedirs(outsubdir, exist_ok=True)
        logger.debug("Outsubdir: %s", outsubdir)

        ref, query = run_dict[str(ref)], run_dict[str(query)]
        logger.info(f"Reference ID: {ref.run_id}, Run ID: {query.run_id}")

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
            scatterstem = (
                Path(outsubdir)
                / f"scatter_{A.name}_run{ref.run_id}_vs_{B.name}_run{query.run_id}"
            )

            plotting_commands.append(
                (
                    get_scatter,
                    [ref.run_id, query.run_id, A, B, scatterstem, outfmts, args],
                )
            )

            # Create Bland-Altman plots
            # info = get_info_text([ref.run_id, query.run_id], run_dict)
            # sys.stderr.write(info)
            blandstem = (
                Path(outsubdir)
                / f"bland-altman_{A.name}_run{ref.run_id}_vs_{B.name}_run{query.run_id}"
            )

            plotting_commands.append(
                (
                    get_bland_altman,
                    [ref.run_id, query.run_id, A, B, blandstem, outfmts, args],
                )
            )

        # Send dataframes for heatmaps, distributions
        for matdata in difference_matrices.values():
            # Write heatmap for each results matrix
            heatstem = (
                Path(outsubdir)
                / f"heatmap_{matdata.name}_run{ref.run_id}_vs_run{query.run_id}"
            )

            plotting_commands.append(
                (
                    get_heatmap,
                    [
                        ref.run_id,
                        query.run_id,
                        matdata,
                        labels,
                        classes,
                        heatstem,
                        outfmts,
                        args,
                    ],
                )
            )
            # Plot distributions of differences to look at normality
            diststem = (
                Path(outsubdir)
                / f"distribution_{matdata.name}_run{ref.run_id}_run{query.run_id}"
            )

            plotting_commands.append(
                (
                    get_distribution,
                    [
                        ref.run_id,
                        query.run_id,
                        matdata,
                        diststem,
                        outfmts,
                        args,
                    ],
                )
            )

        # Record output directories in the logfile
        logger.info(
            "Writing scatter plots to: %s/scatter_*_run%s_run%s.*",
            Path(outsubdir),
            ref.run_id,
            query.run_id,
        )
        logger.info(
            "Writing Bland-Altman plots to: %s/bland-altman_*_run%s_run%s.*",
            Path(outsubdir),
            ref.run_id,
            query.run_id,
        )
        logger.info(
            "Writing heatmaps to: %s/heatmap_*_run%s_run%s.*",
            Path(outsubdir),
            ref.run_id,
            query.run_id,
        )
        logger.info(
            "Writing distribution plots to: %s/distribution_*_run%s_run%s.*",
            Path(outsubdir),
            ref.run_id,
            query.run_id,
        )

        # Run the plotting commands
        for func, argv in plotting_commands:
            result = pool.apply_async(func, argv, {}, callback=logger.debug)
            result.get()

        # Close worker pool
        pool.close()
        pool.join()

        # Generate summary report
        summary_file = (
            Path(outsubdir) / f"summary_run{ref.run_id}_vs_run{query.run_id}.md"
        )
        write_summary(ref, sub_ref, query, sub_query, summary_file)

        # Reset output subdirectory
        # Failure to do this causes nesting
        outsubdir = Path(args.outdir) / "compare"


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
    outfstem: str,
    outfmts: List[str],
    args: Namespace,
) -> None:
    """Write a single heatmap for a comparison between two pyani runs.

    :param run_id:  int, run_id for this run
    :param matdata:  MatrixData object for this heatmap
    :param outfstem:  stem for output graphics files
    :param outfmts:  list of output formats for files
    :param args:  Namespace for command-line arguments
    """
    # Collect logging statements here while in multiprocessing mode
    logs = []
    proname = multiprocessing.current_process().name

    if "absdiff" in matdata.name:
        cmap = ("crest", matdata.data.values.min(), matdata.data.values.max(), None)
    else:
        cmap = ("icefire", matdata.data.values.min(), matdata.data.values.max(), 0)

    for fmt in outfmts:
        outfname = f"{outfstem}.{fmt}"
        logs.append(f"{proname}: Writing graphics to {outfname}")
        params = pyani_graphics.Params(cmap, result_labels, result_classes)

        # Draw heatmap
        GMETHODS[args.method](
            matdata.data,
            outfname,
            title=f"Compare {matdata.name.title().replace('_', ' ')} run {run_a} vs run {run_b}",
            params=params,
        )

        # Be tidy with matplotlib caches
        plt.close("all")

    return "".join(logs)


def get_distribution(
    run_a: int,
    run_b: int,
    matdata: MatrixData,
    outfstem: str,
    outfmts: List[str],
    args: Namespace,
) -> None:
    """Write distribution plots for each matrix type.

    :param run_id:  int, run_id for this run
    :param matdata:  MatrixData object for this distribution plot
    :param outfstem:  stem for output graphics files
    :param outfmts:  list of output formats for files
    :param args:  Namespace for command-line arguments
    """
    # Collect logging statements here while in multiprocessing mode
    logs = []
    proname = multiprocessing.current_process().name

    for fmt in outfmts:
        outfname = f"{outfstem}.{fmt}"
        logs.append(f"{proname}: Writing graphics to {outfname}")

        # Draw distributions
        DISTMETHODS[args.method](
            matdata.data,
            outfname,
            matdata.name,
            title=f"Compare {matdata.name.title().replace('_', ' ')} run {run_a} vs run {run_b}",
        )

        # Be tidy with matplotlib caches
        plt.close("all")

    return "".join(logs)


def get_scatter(
    run_a: int,
    run_b: int,
    matdata1: MatrixData,
    matdata2: MatrixData,
    outfstem: str,
    outfmts: List[str],
    args: Namespace,
) -> None:
    """Write a single scatterplot for a comparison between two pyani runs.

    :param run_a:  int, run_id for the reference
    :param run_b:  int, run_id for the query
    :param matdata1:  MatrixData object for this scatterplot
    :param matdata2:  MatrixData object for this scatterplot
    :param outfstem:  stem for output graphics files
    :param outfmts:  list of output formats for files
    :param args:  Namespace for command-line arguments
    """
    # Collect logging statements here while in multiprocessing mode
    logs = []
    proname = multiprocessing.current_process().name

    extreme = max(abs(matdata1.data.values.min()), abs(matdata1.data.values.max()))
    cmap = ("BuRd", extreme * -1, extreme, None)
    for fmt in outfmts:
        outfname = f"{outfstem}.{fmt}"
        logs.append(f"{proname}: Writing graphics to {outfname}")
        params = pyani_graphics.Params(cmap, {}, {})

        # Draw scatterplot
        SMETHODS[args.method](
            matdata1.data,
            matdata2.data,
            outfname,
            f"{matdata1.name}_{run_a}",
            f"{matdata2.name}_{run_b}",
            title=f"{matdata1.name.title().replace('_', ' ')} run {run_a} vs {matdata2.name.title()} run {run_b}",
            params=params,
        )

        # Be tidy with matplotlib caches
        # plt.close("all")
    return "".join(logs)


def get_bland_altman(
    run_a: int,
    run_b: int,
    matdata1: MatrixData,
    matdata2: MatrixData,
    outfstem: str,
    outfmts: List[str],
    args: Namespace,
) -> None:
    """Write a Bland-Altman plot for a comparison between two pyani runs.

    :param run_a:  int, run_id for the reference
    :param run_b:  int, run_id for the query
    :param matdata1:  MatrixData object for this scatterplot
    :param matdata2:  MatrixData object for this scatterplot
    :param outfstem:  stem for output graphics files
    :param outfmts:   list of output formats for files
    :param args:  Namespace for command-line arguments
    """
    # Collect logging statements here while in multiprocessing mode
    logs = []
    proname = multiprocessing.current_process().name

    extreme = max(abs(matdata1.data.values.min()), abs(matdata1.data.values.max()))
    cmap = ("BuRd", extreme * -1, extreme, None)
    for fmt in outfmts:
        outfname = f"{outfstem}.{fmt}"

        logs.append(f"{proname}: Writing graphics to {outfname}")
        params = pyani_graphics.Params(cmap, {}, {})

        # Draw Bland-Altman
        BMETHODS[args.method](
            matdata1.data,
            matdata2.data,
            outfname,
            # info,
            f"{matdata1.name}",
            f"{matdata2.name}",
            title=f"{matdata1.name.title().replace('_', ' ')} run {run_a} vs {matdata2.name.title()} run {run_b}",
            params=params,
        )

        # Be tidy with matplotlib caches
        plt.close("all")

    return "".join(logs)


def get_info_text(run_ids: List[int], run_dict: Dict) -> str:
    defaults = {
        "anib": {"fragsize": pyani_config.FRAGSIZE},
        "aniblastall": {"fragsize": pyani_config.FRAGSIZE},
        "anim": {
            "maxmatch": False,
            # 'noextend':
        },
        "fastani": {"kmer": 16, "fragLen": 3000, "minFraction": 0.2},
    }
    s = ""
    for id in run_ids:
        id = str(id)
        # ignore = {'dbpath', 'labels', 'classes', 'i', 'indir',  'l', 'log', 'o', 'outdir'}
        keep = {
            "k",
            "kmer",
            "fragLen",
            "minFraction",
            "fragsize",
            "maxmatch",
            "noextend",
        }
        command_list = list(filter(bool, run_dict[id].cmdline.split("-")))
        method = run_dict[id].method.lower()
        arguments = [tuple(x.split()) for x in command_list[1:]]
        legend_info = {k: v for k, v in arguments if k in keep}
        defaults[method].update(legend_info)

        s += f"Run {id}\nMethod: {run_dict[id].method}\n"
        for key in defaults[method]:
            s += f"    {key}: {str(defaults[method][key])}\n"
        s += "\n"
    return s


def write_summary(
    reference: RunData,
    sub_ref: RunMatrices,
    query: RunData,
    sub_query: RunMatrices,
    summary_file: Path,
) -> None:
    """Write out a summary of the comparison.

    :param reference:  RunData object for reference
    :param sub_ref:    RunMatrices object for reference, containing information common with the query
    :param query:      RunData object for query
    :param sub_query:  RunMatrices object for query, containing information common with the reference
    :param summary_file:  Path object for summary file

    """
    output = f"## Run {reference.run_id} (reference) vs Run {query.run_id}  \n"
    output += "### Reference  \n"
    output += f"**ID**: {reference.run_id}  \n"
    output += f"**Method**: {reference.method}  \n"
    output += f"**Command line**: {reference.cmdline}  \n"
    output += f"**Highest ANI**: {sub_ref.identity.data.max().max()}  \n"
    output += f"**Lowest coverage**: {sub_ref.coverage.data.min().min()}  \n"

    output += "### Query  \n"
    output += f"**ID**: {query.run_id}  \n"
    output += f"**Method**: {query.method}  \n"
    output += f"**Command line**: {query.cmdline}  \n"
    output += f"**Highest ANI**: {sub_query.identity.data.max().max()}  \n"
    output += f"**Lowest coverage**: {sub_query.coverage.data.min().min()}"

    with open(summary_file, "w") as ofh:
        ofh.write(output)
