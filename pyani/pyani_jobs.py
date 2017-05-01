# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to manage jobs for pyani.

In order to be a little more consistent behind the scenes for schedulers,
and to allow for a fairly hacky approach to scheduing on SGE, a job
dependency graph is used.

Commands to be run are stored in Jobs. A Job's dependency is stored so that
the Job will not be executed until its dependency is executed.

When used in ANI analysis, the way jobs are used depends on the scheduler.

With multiprocessing, we place all root jobs in a single pool; then all
first-level dependencies will go in a second (dependent) pool that is not run
until the first is completed, and so on. It's not very efficient, but should
work equivalently to the original code that handled asynchronous pools
directly.

With SGE, the dependencies can be managed independently, and effectively
interleaved by the scheduler with no need for pools.

This code is essentially a frozen and cut-down version of pysge
(https://github.com/widdowquinn/pysge)
"""

import os
import time

from .pyani_config import SGE_WAIT

###
# CLASSES


# The Job class describes a single command-line job, with dependencies (jobs
# that must be run first.
class Job:
    """Objects in this class represent individual jobs to be run, with a list
    of dependencies (jobs that must be run first).
    """
    def __init__(self, name, command, queue=None):
        """Instantiates a Job object.

        - name           String describing the job (uniquely)
        - command        String, the valid shell command to run the job
        - queue          String, the SGE queue under which the job shall run
        """
        self.name = name                 # Unique name for the job
        self.queue = queue               # The SGE queue to run the job under
        self.command = command           # Command line to run for this job
        self.script = command
        self.scriptPath = None           # Will hold path to the script file
        self.dependencies = []           # List of jobs to be completed first
        self.submitted = False           # Flag: is job submitted?

    def add_dependency(self, job):
        """Add the passed job to the dependency list for this Job.  This
        Job should not execute until all dependent jobs are completed

        - job     Job to be added to the Job's dependency list
        """
        self.dependencies.append(job)

    def remove_dependency(self, job):
        """Remove the passed job from this Job's dependency list

        - job     Job to be removed from the Job's dependency list
        """
        self.dependencies.remove(job)

    def wait(self, interval=SGE_WAIT):
        """Wait until the job finishes, and poll SGE on its status."""
        finished = False
        while not finished:
            time.sleep(interval)
            interval = min(2 * interval, 60)
            finished = os.system("qstat -j %s > /dev/null" % (self.name))


class JobGroup:
    """ Class that stores a group of jobs, permitting parameter sweeps."""
    def __init__(self, name, command, queue=None, arguments=None):
        """ Instantiate a JobGroup object.  JobGroups allow for the use of
        combinatorial parameter sweeps by using the 'command' and 'arguments'
        arguments.

        - name              String, the JobGroup name
        - command           String, the command to be run, with arguments
                                    specified
        - queue             String, the queue for SGE to use
        - arguments         Dictionary, the values for each parameter as
                            lists of strings, keyed by an identifier for
                            the command string

        For example, to use a command 'my_cmd' with the arguments
        '-foo' and '-bar' having values 1, 2, 3, 4 and 'a', 'b', 'c', 'd' in
        all combinations, respectively, you would pass
        command='my_cmd $SGE_TASK_ID -foo $fooargs -bar $barargs'
        arguments='{'fooargs': ['1','2','3','4'],
                    'barargs': ['a','b','c','d']}
        """
        self.name = name               # Set JobQueue name
        self.queue = queue             # Set SGE queue to request
        self.command = command         # Set command string
        self.dependencies = []         # Create empty list for dependencies
        self.submitted = True          # Set submitted Boolean
        if arguments is not None:
            self.arguments = arguments # Dictionary of arguments for command
        else:
            self.arguments = {}
        self.generate_script()         # Make SGE script for sweep/array

    def generate_script(self):
        """Create the SGE script that will run the jobs in the JobGroup, with
        the passed arguments.
        """
        self.script = ""        # Holds the script string
        total = 1               # total number of jobs in this group

        # for now, SGE_TASK_ID becomes TASK_ID, but we base it at zero
        self.script += """let "TASK_ID=$SGE_TASK_ID - 1"\n"""

        # build the array definitions
        for key in self.arguments.keys():
            values = self.arguments[key]
            line = ("%s_ARRAY=( " % (key))
            for value in values:
                line += value
                line += " "
            line += " )\n"
            self.script += line
            total *= len(values)
        self.script += "\n"

        # now, build the decoding logic in the script
        for key in self.arguments.keys():
            count = len(self.arguments[key])
            self.script += """let "%s_INDEX=$TASK_ID %% %d"\n""" % (key, count)
            self.script += """%s=${%s_ARRAY[$%s_INDEX]}\n""" % (key, key, key)
            self.script += """let "TASK_ID=$TASK_ID / %d"\n""" % (count)

        # now, add the command to run the job
        self.script += "\n"
        self.script += self.command
        self.script += "\n"

        # set the number of tasks in this group
        self.tasks = total

    def add_dependency(self, job):
        """Add the passed job to the dependency list for this JobGroup.  This
        JobGroup should not execute until all dependent jobs are completed

        - job      Job, job to be added to the JobGroup's dependency list
        """
        self.dependencies.append(job)

    def remove_dependency(self, job):
        """ Remove the passed job from this JobGroup's dependency list

        - job      Job, job to be removed from the JobGroup's dependency list
        """
        self.dependencies.remove(job)

    def wait(self, interval=SGE_WAIT):
        """Wait for a defined period, then poll SGE for job status."""
        finished = False
        while not finished:
            time.sleep(interval)
            interval = min(2 * interval, 60)
            finished = os.system("qstat -j %s > /dev/null" % (self.name))
