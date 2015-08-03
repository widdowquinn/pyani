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
    self.script = command            # 
    self.scriptPath = None           # Will hold path to the script file
    self.dependencies = []           # List of jobs that must be submitted
                                     # before this job may be submitted
    self.submitted = True            # Flag indicating whether the job has
                                     # already been submitted

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

  def wait(self, interval=5):
    """Wait for an interval
    """
    finished = False
    while not finished:
        time.sleep(interval)
        interval = min( 2 * interval, 60 )
        finished = os.system( "qstat -j %s > /dev/null" % (self.name) )

###
# FUNCTIONS


def build_directories(root_dir):
  """Constructs the subdirectories output, stderr, stdout, and jobs in the
  passed root directory. These subdirectories have the following roles:

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
  directories = [os.path.join(root_dir, subdir) for subdir in
                 ("output", "stderr", "stdout", "jobs")]
  [os.mkdir(dirname) for dirname in directories if not os.path.exists(dirname)]


def build_job_scripts(root_dir, jobs):
  """Constructs the script for each passed Job in the jobs iterable

  - root_dir      Path to output directory
  """
  # Loop over the job list, creating each job script in turn, and then adding
  # scriptPath to the Job object
  for job in jobs:
      scriptPath = os.path.join(root_dir, "jobs", job.name)
      scriptFile = file(scriptPath, "w")
      scriptFile.write("#!/bin/sh\n#$ -S /bin/bash\n%s\n" % job.script)
      scriptFile.close()
      job.scriptPath = scriptPath


def extract_submittable_jobs(waiting):
  """Obtain a list of jobs that are able to be submitted from the passed
  list of pending jobs

  - waiting           List of Job objects
  """
  submittable = []            # Holds jobs that are able to be submitted
  # Loop over each job, and check all the subjobs in that job's dependency
  # list.  If there are any, and all of these have been submitted, then
  # append the job to the list of submittable jobs.
  for job in waiting:
      unsatisfied = sum([(subjob.submitted is False) for subjob in \
                         job.dependencies])
      if 0 == unsatisfied:
          submittable.append(job)
  return submittable


def submit_safe_jobs(root_dir, jobs):
  """Submit the passed list of jobs to the Grid Engine server, using the passed
  directory as the root for scheduler output.

  - root_dir      Path to output directory
  - jobs          Iterable of Job objects
  """
  # Loop over each job, constructing the SGE command-line based on job settings
  for job in jobs:
      job.out = os.path.join(root_dir, "stdout")
      job.err = os.path.join(root_dir, "stderr")
      # Add the job name, current working directory, and SGE stdout and stderr
      # directories to the SGE command line
      args = " -N %s " % (job.name)
      args += " -cwd "
      args += " -o %s -e %s " % (job.out, job.err)
      # If a queue is specified, add this to the SGE command line
      if job.queue != None and job.queue in local_queues:
          args += local_queues[job.queue]
          #args += "-q %s " % job.queue
      # If there are dependencies for this job, hold the job until they are
      # complete
      if len(job.dependencies) > 0:
          args += "-hold_jid "
          for dep in job.dependencies:
              args += dep.name + ","
          args = args[:-1]
      # Build the qsub SGE commandline
      qsubcmd = ("/usr/nfs/sge_root/bin/lx24-x86/qsub %s %s" % \
                 (args, job.scriptPath)) 
      #print qsubcmd                   # Show the command to the user
      os.system(qsubcmd)               # Run the command
      job.submitted = True             # Set the job's submitted flag to True


def submit_jobs(root_dir, jobs):
  """ Submit each of the passed jobs to the SGE server, using the passed
  directory as root for SGE output.

  - root_dir       Path to output directory
  - jobs           List of Job objects
  """
  waiting = list(jobs)                 # List of jobs still to be done
  # Loop over the list of pending jobs, while there still are any
  while len(waiting) > 0:
      # extract submittable jobs
      submittable = extract_submittable_jobs(waiting)
      # run those jobs
      submit_safe_jobs(directory, submittable)
      # remove those from the waiting list
      map(waiting.remove, submittable)


def build_and_submit_jobs(root_dir, jobs):
  """Submits the passed iterable of Job objects to SGE, placing SGE's output in
  the passed root directory

  - root_dir   Root directory for SGE and job output

  - jobs       List of Job objects, describing each job to be submitted
  """
  # If the passed set of jobs is not a list, turn it into one.  This makes the
  # use of a single JobGroup a little more intutitive
  if type(jobs) != type([1]):
      jobs = [jobs]
    
  # Build and submit the passed jobs
  build_directories(root_dir)       # build all necessary directories
  build_job_scripts(root_dir, jobs) # build job scripts
  submit_jobs(root_dir, jobs)       # submit the jobs to SGE




