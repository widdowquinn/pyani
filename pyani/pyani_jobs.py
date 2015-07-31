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
"""

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
    self.script = command            # 
    self.scriptPath = None           # Will hold path to the script file
    self.dependencies = []           # List of jobs that must be submitted
                                     # before this job may be submitted
    self.submitted = True            # Flag indicating whether the job has
                                     # already been submitted

  def add_dependency(self, job):
    """Add the passed job to the dependency list for this JobGroup.  This
    JobGroup should not execute until all dependent jobs are completed

    - job             Job to be added to the JobGroup's dependency list
    """
    self.dependencies.append(job)

  def remove_dependency(self, job):
    """Remove the passed job from this JobGroup's dependency list

    - job             Job to be removed from the JobGroup's dependency list
    """
    self.dependencies.remove(job)

