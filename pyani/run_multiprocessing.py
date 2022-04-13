# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2013-2019
# (c) University of Strathclyde 2019
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G1 1XQ
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016-2019 The James Hutton Institute
# Copyright (c) 2019 University of Strathclyde
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
"""Code to run a set of command-line jobs using multiprocessing.

For parallelisation on multi-core desktop/laptop systems, etc. we use
Python's multiprocessing module to distribute command-line jobs.
"""

import multiprocessing
import subprocess
import sys

from logging import Logger
from typing import List, Optional

from .pyani_jobs import Job


# Run a job dependency graph with multiprocessing
def run_dependency_graph(
    jobgraph, workers: Optional[int] = None, logger: Optional[Logger] = None
) -> int:
    """Create and run pools of jobs based on the passed jobgraph.

    :param jobgraph:  list of jobs, which may have dependencies.
    :param workers:  int, number of workers to use with multiprocessing
    :param logger: a logger module logger (optional)

    The strategy here is to loop over each job in the list of jobs (jobgraph),
    and create/populate a series of Sets of commands, to be run in
    reverse order with multiprocessing_run as asynchronous pools.
    """
    cmdsets = []  # type: List
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


def populate_cmdsets(job: Job, cmdsets: List, depth: int) -> List:
    """Create list of jobsets at different depths of dependency tree.

    :param job:
    :param cmdsets:
    :param depth:

    This is a recursive function (is there something quicker in the itertools
    module?) that descends each 'root' job in turn, populating each
    """
    if len(cmdsets) < depth:
        cmdsets.append(set())
    cmdsets[depth - 1].add(job.command)

    # Return now if there are no dependencies
    # (I think I can remove this check/return)
    if not job.dependencies:
        return cmdsets

    # There are dependencies, so add these to the command sets
    for j in job.dependencies:
        cmdsets = populate_cmdsets(j, cmdsets, depth + 1)
    return cmdsets


# Run a set of command lines using multiprocessing
def multiprocessing_run(
    cmdlines: List, workers: Optional[int] = None, logger: Optional[Logger] = None
) -> int:
    """Distributes passed command-line jobs using multiprocessing.

    :param cmdlines:  iterable, command line strings
    :param workers:  int, number of workers to use for multiprocessing

    Returns the sum of exit codes from each job that was run. If
    all goes well, this should be 0. Anything else and the calling
    function should act accordingly.

    Sends a warning to the logger if a comparison fails; the warning consists of the specific command line that has failed.
    """
    # Run jobs
    # If workers is None or greater than the number of cores available,
    # it will be set to the maximum number of cores
    pool = multiprocessing.Pool(processes=workers)
    results = [
        pool.apply_async(
            subprocess.run,
            (str(cline),),
            {
                "shell": sys.platform != "win32",
                "stdout": subprocess.PIPE,
                "stderr": subprocess.PIPE,
            },
        )
        for cline in cmdlines
    ]
    pool.close()
    pool.join()

    # Print a warning so individual runs that fail can be identified
    for r in results:
        if r.get().returncode != 0:
            if logger:  # Try to be informative, if the logger module is being used
                logger.warning(f"Comparison failed: {' '.join(r.get().args)}")

    return sum([r.get().returncode for r in results])
