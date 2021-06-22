.. _pyani-scheduler:

====================
Use With a Scheduler
====================

^^^^^^^^^^^^^^^^^^^^
Sun/Open Grid Engine
^^^^^^^^^^^^^^^^^^^^

The ``--scheduler SGE`` argument allows one to use ``pyani`` with an an SGE-type scheduler.

In order for this work, one must be able to submit jobs using the ``qsub`` command. By default, this will batch the pairwise comparisons in array jobs of 10,000, in order to avoid clogging the scheduler queue. Each comparison will be run as a single-core task in an array job.

---------------------------
Arguments assigned by Pyani
---------------------------

The following arguments will be automatically set:

.. code-block:: bash

   -N job_name  # this is the value passed to `--name`
    -cwd
    -o ./stdout  # cwd/ + "stdout"
    -e ./stderr  # cwd/ + "stderr"


--------------------
Modifiable arguments
--------------------

The number of pairwise comparisons submitted per chunk can be modified using:

.. code-block:: bash

    --SGEgroupsize *number*


The job prefix to use can be modified using:

.. code-block:: bash

    --jobprefix *prefix*


-------------------------------
Specifying additional arguments
-------------------------------

Additional SGE arguments may be specified with:

.. code-block:: bash

    --SGEargs "<your arguments here>"


Additional arguments must be specified as a string which includes all flags and their values. For instance, to specify a value for the memory resource:

.. code-block:: bash

  "-l h_mem=64G"

Or for both memory and run time:

.. code-block:: bash

  "-l h_mem=64G -l h_rt=02:00:00"


An alternative to listing all desired options on the command line is to pass an ``optionfile``:

.. code-block:: bash

    "-@ optionfile"


This file only contains comments and the flags to be passed (do not include the ``#$`` in front of the arguments):

.. code-block:: bash

   # Memory to assign to the job
    -l h_mem=64G

    # Time to allow for job HH:MM:SS
    -l h_rt=10:20:00

    # Notification email address
    -M email@domain.com

    # Send notifications when job 'b'egins, 'a'borts (or is rescheduled), 'e'nds, or is 's'uspended
    -m baes


For more information on using an SGE/OGE scheduler, see:

- `Open Grid Scheduler <http://gridscheduler.sourceforge.net/>`_
