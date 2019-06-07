# -*- coding: utf-8 -*-
"""Code to run a set of command-line jobs using SGE/Grid Engine.

For parallelisation on multi-node system, we use some custom code to submit
jobs.

(c) The James Hutton Institute 2013-2017
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD2 5DA,
Scotland,
UK

The MIT License

Copyright (c) 2013-2017 The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import itertools
import os
import subprocess  # nosec

from collections import defaultdict

from . import pyani_config
from .pyani_jobs import JobGroup


def split_seq(iterable, size):
    """Split a passed iterable into chunks of a given size."""
    elm = iter(iterable)
    item = list(itertools.islice(elm, size))
    while item:
        yield item
        item = list(itertools.islice(elm, size))


# Build a list of SGE jobs from a graph
def build_joblist(jobgraph):
    """Return a list of jobs, from a passed jobgraph."""
    jobset = set()
    for job in jobgraph:
        jobset = populate_jobset(job, jobset, depth=1)
    return list(jobset)


# Convert joblist into jobgroups
def compile_jobgroups_from_joblist(joblist, jgprefix, sgegroupsize):
    """Return list of jobgroups, rather than list of jobs."""
    jobcmds = defaultdict(list)
    for job in joblist:
        jobcmds[job.command.split(" ", 1)[0]].append(job.command)
    jobgroups = []
    for cmds in list(jobcmds.items()):
        # Break arglist up into batches of sgegroupsize (default: 10,000)
        sublists = split_seq(cmds[1], sgegroupsize)
        count = 0
        for sublist in sublists:
            count += 1
            sge_jobcmdlist = ['"%s"' % jc for jc in sublist]
            jobgroups.append(
                JobGroup(
                    "%s_%d" % (jgprefix, count),
                    "$cmds",
                    arguments={"cmds": sge_jobcmdlist},
                )
            )
    return jobgroups


# Run a job dependency graph, with SGE
def run_dependency_graph(
    jobgraph, logger=None, jgprefix="ANIm_SGE_JG", sgegroupsize=10000, sgeargs=None
):
    """Create and runs SGE scripts for jobs based on passed jobgraph.

    - jobgraph - list of jobs, which may have dependencies.
    - verbose - flag for multiprocessing verbosity
    - logger - a logger module logger (optional)
    - jgprefix - a prefix for the submitted jobs, in the scheduler
    - sgegroupsize - the maximum size for an array job submission
    - sgeargs - additional arguments to qsub

    The strategy here is to loop over each job in the dependency graph
    and, because we expect a single main delta-filter (wrapped) job,
    with a single nucmer dependency for each analysis, we can split
    the dependency graph into two lists of corresponding jobs, and
    run the corresponding nucmer jobs before the delta-filter jobs.
    """
    jobs_main = []  # Can be run first, before deps
    jobs_deps = []  # Depend on the main jobs

    # Try to be informative by telling the user what jobs will run
    dep_count = 0  # how many dependencies are there
    if logger:
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
    build_and_submit_jobs(os.curdir, jobgroups, sgeargs)
    logger.info("Waiting for SGE-submitted jobs to finish (polling)")
    for job in jobgroups:
        job.wait()


def populate_jobset(job, jobset, depth):
    """Create set of jobs reflecting dependency tree.

    The set contains jobs at different depths of the dependency tree,
    retaining dependencies as strings, not Jobs.
    """
    jobset.add(job)
    if not job.dependencies:
        return jobset
    for j in job.dependencies:
        jobset = populate_jobset(j, jobset, depth + 1)
    return jobset


def build_directories(root_dir):
    """Construct the subdirectories output, stderr, stdout, and jobs.

    Subdirectories are created in the passed root directory. These
    subdirectories have the following roles:

        jobs             Stores the scripts for each job
        stderr           Stores the stderr output from SGE
        stdout           Stores the stdout output from SGE
        output           Stores output (if the scripts place the output here)

    - root_dir   Path to the top-level directory for creation of subdirectories
    """
    # If the root directory doesn't exist, create it
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)

    # Create subdirectories
    directories = [
        os.path.join(root_dir, subdir)
        for subdir in ("output", "stderr", "stdout", "jobs")
    ]
    for dirname in directories:
        os.makedirs(dirname, exist_ok=True)


def build_job_scripts(root_dir, jobs):
    """Construct script for each passed Job in the jobs iterable.

    - root_dir      Path to output directory
    """
    # Loop over the job list, creating each job script in turn, and then adding
    # scriptPath to the Job object
    for job in jobs:
        scriptpath = os.path.join(root_dir, "jobs", job.name)
        with open(scriptpath, "w") as scriptfile:
            scriptfile.write("#!/bin/sh\n#$ -S /bin/bash\n%s\n" % job.script)
        job.scriptpath = scriptpath


def extract_submittable_jobs(waiting):
    """Obtain list of jobs that are able to be submitted from pending list.

    - waiting           List of Job objects
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


def submit_safe_jobs(root_dir, jobs, sgeargs=None):
    """Submit passed list of jobs to SGE server with dir as root for output.

    - root_dir      Path to output directory
    - jobs          Iterable of Job objects
    """
    # Loop over each job, constructing SGE command-line based on job settings
    for job in jobs:
        job.out = os.path.join(root_dir, "stdout")
        job.err = os.path.join(root_dir, "stderr")

        # Add the job name, current working directory, and SGE stdout/stderr
        # directories to the SGE command line
        args = " -N %s " % (job.name)
        args += " -cwd "
        args += " -o %s -e %s " % (job.out, job.err)

        # If a queue is specified, add this to the SGE command line
        # LP: This has an undeclared variable, not sure why - delete?
        # if job.queue is not None and job.queue in local_queues:
        #    args += local_queues[job.queue]

        # If the job is actually a JobGroup, add the task numbering argument
        if isinstance(job, JobGroup):
            args += "-t 1:%d " % (job.tasks)

        # If there are dependencies for this job, hold the job until they are
        # complete
        if job.dependencies:
            args += "-hold_jid "
            for dep in job.dependencies:
                args += dep.name + ","
            args = args[:-1]

        # Build the qsub SGE commandline (passing local environment)
        qsubcmd = "%s -V %s %s" % (pyani_config.QSUB_DEFAULT, args, job.scriptpath)
        if sgeargs is not None:
            qsubcmd = "%s %s" % (qsubcmd, sgeargs)
        # We've considered Bandit warnings B404,B603 and silence
        subprocess.call(qsubcmd, shell=False)  # nosec
        job.submitted = True  # Set the job's submitted flag to True


def submit_jobs(root_dir, jobs, sgeargs=None):
    """Submit passed jobs to SGE server with passed directory as root.

    - root_dir       Path to output directory
    - jobs           List of Job objects
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


def build_and_submit_jobs(root_dir, jobs, sgeargs=None):
    """Submit passed iterable of Job objects to SGE.

    This places SGE's output in the passed root directory

    - root_dir   Root directory for SGE and job output
    - jobs       List of Job objects, describing each job to be submitted
    - sgeargs    Additional arguments to qsub
    """
    # If the passed set of jobs is not a list, turn it into one. This makes the
    # use of a single JobGroup a little more intutitive
    if not isinstance(jobs, list):
        jobs = [jobs]

    # Build and submit the passed jobs
    build_directories(root_dir)  # build all necessary directories
    build_job_scripts(root_dir, jobs)  # build job scripts
    submit_jobs(root_dir, jobs, sgeargs)  # submit the jobs to SGE
