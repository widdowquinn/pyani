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
# Copyright (c) 2013-2019 The James Hutton Institute
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
"""Code to run a set of command-line jobs using SGE/Grid Engine.

For parallelisation on multi-node system, we use some custom code to submit
jobs.
"""

import itertools
import logging
import os

from collections import defaultdict
from pathlib import Path
from typing import Dict, Generator, Iterable, List, Optional, Set

from . import pyani_config
from .pyani_jobs import Job, JobGroup


def split_seq(iterable: Iterable, size: int) -> Generator:
    """Split a passed iterable into chunks of a given size.

    :param iterable:  iterable
    :param size:  int, number of items to retun in each chunk
    """
    elm = iter(iterable)
    item = list(itertools.islice(elm, size))
    while item:
        yield item
        item = list(itertools.islice(elm, size))


# Build a list of SGE jobs from a graph
def build_joblist(jobgraph) -> List:
    """Return a list of jobs, from a passed jobgraph.

    :param jobgraph:
    """
    logger = logging.getLogger(__name__)

    jobset = set()  # type: Set
    for job in jobgraph:
        jobset = populate_jobset(job, jobset, depth=1)

    logger.debug("built jobset: %s", jobset)
    return list(jobset)


# Convert joblist into jobgroups
def compile_jobgroups_from_joblist(
    joblist: List, jgprefix: str, sgegroupsize: int
) -> List:
    """Return list of jobgroups, rather than list of jobs.

    :param joblist:
    :param jgprefix:  str, prefix for SGE jobgroup
    :param sgegroupsize:  int, number of jobs in each SGE jobgroup
    """
    jobcmds = defaultdict(list)  # type: Dict[str, List[str]]
    for job in joblist:
        jobcmds[job.command.split(" ", 1)[0]].append(job.command)
    jobgroups = []  # type: List
    for cmds in list(jobcmds.items()):
        # Break arglist up into batches of sgegroupsize (default: 10,000)
        sublists = split_seq(cmds[1], sgegroupsize)
        count = 0
        for sublist in sublists:
            count += 1
            sge_jobcmdlist = [f'"{jc}"' for jc in sublist]
            jobgroups.append(
                JobGroup(
                    f"{jgprefix}_{count}", "$cmds", arguments={"cmds": sge_jobcmdlist}
                )
            )
    return jobgroups


# Run a job dependency graph, with SGE
def run_dependency_graph(
    jobgraph,
    jgprefix: str = "ANIm_SGE_JG",
    sgegroupsize: int = 10000,
    sgeargs: Optional[str] = None,
) -> None:
    """Create and runs SGE scripts for jobs based on passed jobgraph.

    :param jobgraph: list of jobs, which may have dependencies.
    :param verbose: flag for multiprocessing verbosity
    :param jgprefix: a prefix for the submitted jobs, in the scheduler
    :param sgegroupsize: the maximum size for an array job submission
    :param sgeargs: additional arguments to qsub

    The strategy here is to loop over each job in the dependency graph
    and, because we expect a single main delta-filter (wrapped) job,
    with a single nucmer dependency for each analysis, we can split
    the dependency graph into two lists of corresponding jobs, and
    run the corresponding nucmer jobs before the delta-filter jobs.
    """
    logger = logging.getLogger(__name__)

    logger.debug("Received jobgraph with %d jobs", len(jobgraph))

    jobs_main = []  # Can be run first, before deps
    jobs_deps = []  # Depend on the main jobs

    # Try to be informative by telling the user what jobs will run
    dep_count = 0  # how many dependencies are there
    logger.info("Jobs to run with scheduler")
    for job in jobgraph:
        logger.info("{0}: {1}".format(job.name, job.command))
        jobs_main.append(job)
        if job.dependencies:
            dep_count += len(job.dependencies)
            for dep in job.dependencies:
                logger.info("\t[^ depends on: %s (%s)]", dep.name, dep.command)
                jobs_deps.append(dep)
    logger.info("There are %d job dependencies" % dep_count)
    # Clear dependencies in main group
    for job in jobs_main:
        job.dependencies = []

    # We can use an array (or series of arrays) to schedule our jobs.
    # This cuts down on social problems with long job lists choking up
    # the queue.
    # We split the main and dependent jobs into separate JobGroups.
    # These JobGroups are paired, in order
    logger.info("Compiling main and dependent jobs into separate JobGroups")
    maingroups = compile_jobgroups_from_joblist(
        jobs_main, jgprefix + "_main", sgegroupsize
    )
    depgroups = compile_jobgroups_from_joblist(
        jobs_deps, jgprefix + "_deps", sgegroupsize
    )

    # Assign dependencies to jobgroups
    for mgp, dgp in zip(maingroups, depgroups):
        mgp.add_dependency(dgp)
    jobgroups = maingroups + depgroups

    # Send jobs to scheduler
    logger.info("Running jobs with scheduler...")
    logger.info("Jobs passed to scheduler in order:")
    for job in jobgroups:
        logger.info("\t%s" % job.name)
    build_and_submit_jobs(Path.cwd(), jobgroups, sgeargs)
    logger.info("Waiting for SGE-submitted jobs to finish (polling)")
    for job in jobgroups:
        job.wait()


def populate_jobset(job: Job, jobset: Set, depth: int) -> Set:
    """Create set of jobs reflecting dependency tree.

    :param job:
    :param jobset:
    :param depth:

    The set contains jobs at different depths of the dependency tree,
    retaining dependencies as strings, not Jobs.
    """
    jobset.add(job)
    if not job.dependencies:
        return jobset
    for j in job.dependencies:
        jobset = populate_jobset(j, jobset, depth + 1)
    return jobset


def build_directories(root_dir: Path) -> None:
    """Construct the subdirectories output, stderr, stdout, and jobs.

    :param root_dir:  path of root directory in which to place output

    Subdirectories are created in the passed root directory. These
    subdirectories have the following roles:

        jobs             Stores the scripts for each job
        stderr           Stores the stderr output from SGE
        stdout           Stores the stdout output from SGE
        output           Stores output (if the scripts place the output here)

    - root_dir   Path to the top-level directory for creation of subdirectories
    """
    # If the root directory doesn't exist, create it
    if not root_dir.exists():
        root_dir.mkdir(exist_ok=True)

    # Create subdirectories
    directories = [
        root_dir / subdir for subdir in ("output", "stderr", "stdout", "jobs")
    ]
    for dirname in directories:
        dirname.mkdir(exist_ok=True)


def build_job_scripts(root_dir: Path, jobs: List) -> None:
    """Construct script for each passed Job in the jobs iterable.

    :param root_dir:  Path to output directory
    :param jobs:
    """
    # Loop over the job list, creating each job script in turn, and then adding
    # scriptPath to the Job object
    for job in jobs:
        scriptpath = root_dir / "jobs" / job.name
        with open(scriptpath, "w") as scriptfile:
            scriptfile.write(f"#!/bin/sh\n#$ -S /bin/bash\n{job.script}\n")
        job.scriptpath = scriptpath


def extract_submittable_jobs(waiting: List) -> List:
    """Obtain list of jobs that are able to be submitted from pending list.

    :param waiting:  list of Job objects
    """
    submittable = set()  # Holds jobs that are able to be submitted
    # Loop over each job, and check all the subjobs in that job's dependency
    # list.  If there are any, and all of these have been submitted, then
    # append the job to the list of submittable jobs.
    for job in waiting:
        unsatisfied = sum([(subjob.submitted is False) for subjob in job.dependencies])
        if unsatisfied == 0:
            submittable.add(job)
    return list(submittable)


def submit_safe_jobs(
    root_dir: Path, jobs: Iterable, sgeargs: Optional[str] = None
) -> None:
    """Submit passed list of jobs to SGE server with dir as root for output.

    :param root_dir:  path to output directory
    :param jobs:  iterable of Job objects
    :param sgeargs:  str, additional arguments for qsub
    """
    logger = logging.getLogger(__name__)
    logger.debug("Received %s jobs", len(jobs))

    # Loop over each job, constructing SGE command-line based on job settings
    for job in jobs:
        job.out = root_dir / "stdout"
        job.err = root_dir / "stderr"

        # Add the job name, current working directory, and SGE stdout/stderr
        # directories to the SGE command line
        args = f" -N {job.name} "
        args += " -cwd "
        args += f" -o {job.out} -e {job.err} "

        # If a queue is specified, add this to the SGE command line
        # LP: This has an undeclared variable, not sure why - delete?
        # if job.queue is not None and job.queue in local_queues:
        #    args += local_queues[job.queue]

        # If the job is actually a JobGroup, add the task numbering argument
        if isinstance(job, JobGroup):
            args += f"-t 1:{job.tasks} "

        # If there are dependencies for this job, hold the job until they are
        # complete
        if job.dependencies:
            args += "-hold_jid "
            for dep in job.dependencies:
                args += dep.name + ","
            args = args[:-1]

        # Build the qsub SGE commandline (passing local environment)
        qsubcmd = f"{pyani_config.QSUB_DEFAULT} -V {args} {job.scriptpath}"
        if sgeargs is not None:
            qsubcmd = f"{qsubcmd} {sgeargs}"
        # We've considered Bandit warnings B404,B603 and silence
        # subprocess.call(qsubcmd, shell=False)  # nosec
        os.system(qsubcmd)
        job.submitted = True  # Set the job's submitted flag to True


def submit_jobs(root_dir: Path, jobs: Iterable, sgeargs: Optional[str] = None) -> None:
    """Submit passed jobs to SGE server with passed directory as root.

    :param root_dir:  path to output directory
    :param jobs:  list of Job objects
    :param sgeargs:  str, additional arguments for qsub
    """
    waiting = list(jobs)  # List of jobs still to be done
    # Loop over the list of pending jobs, while there still are any
    while waiting:
        # extract submittable jobs
        submittable = extract_submittable_jobs(waiting)
        # run those jobs
        submit_safe_jobs(root_dir, submittable, sgeargs)
        # remove those from the waiting list
        for job in submittable:
            waiting.remove(job)


def build_and_submit_jobs(
    root_dir: Path, jobs: Iterable, sgeargs: Optional[str] = None
) -> None:
    """Submit passed iterable of Job objects to SGE.

    :param root_dir:  root directory for SGE and job output
    :param jobs:  list of Job objects, describing each job to be submitted
    :param sgeargs:  str, additional arguments to qsub

    This places SGE's output in the passed root directory
    """
    # If the passed set of jobs is not a list, turn it into one. This makes the
    # use of a single JobGroup a little more intutitive
    if not isinstance(jobs, list):
        jobs = [jobs]

    # Build and submit the passed jobs
    build_directories(root_dir)  # build all necessary directories
    build_job_scripts(root_dir, jobs)  # build job scripts
    submit_jobs(root_dir, jobs, sgeargs)  # submit the jobs to SGE
