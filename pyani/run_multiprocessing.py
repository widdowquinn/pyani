# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to run a set of command-line jobs using multiprocessing.

For parallelisation on multi-core desktop/laptop systems, etc. we use
Python's multiprocessing module to distribute command-line jobs.
"""

import multiprocessing
import subprocess
import sys

CUMRETVAL = 0


# Run a job dependency graph with multiprocessing
def run_dependency_graph(jobgraph, workers=None, logger=None):
    """Creates and runs pools of jobs based on the passed jobgraph.

    - jobgraph - list of jobs, which may have dependencies.
    - verbose - flag for multiprocessing verbosity
    - logger - a logger module logger (optional)

    The strategy here is to loop over each job in the list of jobs (jobgraph),
    and create/populate a series of Sets of commands, to be run in
    reverse order with multiprocessing_run as asynchronous pools.
    """
    cmdsets = []
    for job in jobgraph:
        cmdsets = populate_cmdsets(job, cmdsets, depth=1)

    # Put command sets in reverse order, and submit to multiprocessing_run
    cmdsets.reverse()
    cumretval = 0
    for cmdset in cmdsets:
        if logger:  # Try to be informative, if the logger module is being used
            logger.info("Command pool now running:")
            for cmd in cmdset:
                logger.info(cmd)
        cumretval += multiprocessing_run(cmdset, workers)
        if logger:  # Try to be informative, if the logger module is being used
            logger.info("Command pool done.")
    return cumretval


def populate_cmdsets(job, cmdsets, depth):
    """Creates a list of sets containing jobs at different depths of the
    dependency tree.

    This is a recursive function (is there something quicker in the itertools
    module?) that descends each 'root' job in turn, populating each
    """
    if len(cmdsets) < depth:
        cmdsets.append(set())
    cmdsets[depth-1].add(job.command)
    if len(job.dependencies) == 0:
        return cmdsets
    for j in job.dependencies:
        cmdsets = populate_cmdsets(j, cmdsets, depth+1)
    return cmdsets


# Run a set of command lines using multiprocessing
def multiprocessing_run(cmdlines, workers=None):
    """Distributes passed command-line jobs using multiprocessing.

    - cmdlines - an iterable of command line strings

    Returns the sum of exit codes from each job that was run. If
    all goes well, this should be 0. Anything else and the calling
    function should act accordingly.
    """
    # Run jobs
    # If workers is None or greater than the number of cores available,
    # it will be set to the maximum number of cores
    pool = multiprocessing.Pool(processes=workers)
    results = [pool.apply_async(subprocess.run, (str(cline), ),
                                {'shell': sys.platform != "win32",
                                 'stdout': subprocess.PIPE,
                                 'stderr': subprocess.PIPE})
               for cline in cmdlines]
    pool.close()
    pool.join()
    return sum([r.get().returncode for r in results])
