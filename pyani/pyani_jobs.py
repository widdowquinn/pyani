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
from subprocess import Popen, PIPE
from typing import Any, Dict, List, Optional

from .pyani_config import SGE_WAIT

###
# CLASSES


# The Job class describes a single command-line job, with dependencies (jobs
# that must be run first.
class Job(object):

    """Individual job to be run, with list of dependencies."""

    def __init__(self, name: str, command: str, queue: Optional[str] = None) -> None:
        """Instantiate a Job object.

        :param name:           String describing the job (uniquely)
        :param command:        String, the valid shell command to run the job
        :param queue:          String, the SGE queue under which the job shall run
        """
        self.name = name  # Unique name for the job
        self.queue = queue  # The SGE queue to run the job under
        self.command = command  # Command line to run for this job
        self.script = command
        self.scriptPath = None  # type: Optional[Any]
        self.dependencies = []  # type: List[Any]
        self.submitted = False  # type: bool
        self.finished = False  # type: int

    def add_dependency(self, job) -> None:
        """Add passed job to the dependency list for this Job.

        :param job:  Job to be added to the Job's dependency list

        This Job should not execute until all dependent jobs are completed.
        """
        self.dependencies.append(job)

    def remove_dependency(self, job) -> None:
        """Remove passed job from this Job's dependency list.

        :param job:     Job to be removed from the Job's dependency list
        """
        self.dependencies.remove(job)

    def wait(self, interval: float = SGE_WAIT) -> None:
        """Wait until the job finishes, and poll SGE on its status.

        :param interval:  float, number of seconds to wait before polling SGE
        """
        while not self.finished:
            time.sleep(interval)
            interval = min(2.0 * interval, 60)
            self.finished = os.system(f"qstat -j {self.name} > /dev/null")

class JobGroup(object):

    """Class that stores a group of jobs, permitting parameter sweeps."""

    def __init__(
        self,
        name: str,
        command: str,
        scheduler: str,
        queue: Optional[str] = None,
        arguments: Optional[Dict[str, List[Any]]] = None,
        
    ) -> None:
        """Instantiate a JobGroup object.

        JobGroups allow for the use of combinatorial parameter sweeps by
        using the 'command' and 'arguments' arguments.

        :param name:  str, the JobGroup name
        :param command:  str, the command to be run, with arguments specified
        :param queue:   str, , the queue for SGE to use
        :param arguments:  dict, the values for each parameter as lists of strings, keyed by an identifier for the command string

        For example, to use a command 'my_cmd' with the arguments
        '-foo' and '-bar' having values 1, 2, 3, 4 and 'a', 'b', 'c', 'd' in
        all combinations, respectively, you would pass
        command='my_cmd $SLURM_TASK_ID -foo $fooargs -bar $barargs'
        arguments='{'fooargs': ['1','2','3','4'],
                    'barargs': ['a','b','c','d']}
        """
        self.name = name  # Set JobQueue name
        self.queue = queue  # Set SGE queue to request
        self.command = command  # Set command string
        self.dependencies = []  # type: List[Any]
        self.submitted = False  # type: bool
        self.finished = False  # type: int
        self.scheduler = scheduler
        if arguments is not None:
            self.arguments = arguments  # Dictionary of arguments for command
        else:
            self.arguments = {}
        self.generate_script()  # Make SGE script for sweep/array

    def generate_script(self) -> None:
        """Create the SLURM script that will run the jobs in the JobGroup."""
        self.script = ""  # type: str
        total = 1  # total number of jobs in this group

        # for now, SLURM_TASK_ID becomes TASK_ID, but we base it at zero
        #self.script += """let "TASK_ID=$SLURM_TASK_ID - 1"\n"""
        self.script += """let "TASK_ID=$SLURM_ARRAY_TASK_ID - 1"\n"""
        # build the array definitions
        for key in sorted(self.arguments.keys()):
            # The keys are sorted for py3.5 compatibility with tests
            values = self.arguments[key]
            line = "%s_ARRAY=( " % (key)
            for value in values:
                line += value
                line += " "
            line += " )\n"
            self.script += line
            total *= len(values)
        self.script += "\n"

        # now, build the decoding logic in the script
        for key in sorted(self.arguments.keys()):
            # The keys are sorted for py3.5 compatibility with tests
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

    def add_dependency(self, job) -> None:
        """Add the passed job to the dependency list for this JobGroup.

        :param job:  Job, job to be added to the JobGroup's dependency list

        This JobGroup should not execute until all dependent jobs are
        completed
        """
        self.dependencies.append(job)

    def remove_dependency(self, job) -> None:
        """Remove passed job from this JobGroup's dependency list.

        :param job:  Job, job to be removed from the JobGroup's dependency list
        """
        self.dependencies.remove(job)

    def wait(self, interval: float = SGE_WAIT) -> None:
        """Wait for a defined period, then poll SGE for job status.

        :param interval:  int, seconds to wait before polling SGE
        """
        while not self.finished:
            time.sleep(interval)
            interval = min(2 * interval, 60)
            
            if self.scheduler.lower() == "sge" :  # hpc is SGE
                self.finished = os.system("qstat -j %s > /dev/null" % (self.name))

            elif self.scheduler.lower() == "slurm" :  # hpc is SLURM
                print("Scheduler slurm: squeue -n %s" % (self.name), "finished? ", self.finished)
                cmd = "squeue -n %s | tail -n+2 | wc -l" % (self.name)
                count = get_cmd_output(cmd)
            
                if int(count) == 0:
                    self.finished = True
                    print("Finished ", self.finished)

def get_cmd_output(cmd):
    
    """ call subprocess popen to get command stdout """

    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    out, error =  p.communicate()
    return str(out, 'utf-8')
