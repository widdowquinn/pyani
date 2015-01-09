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


# Run a set of command lines using multiprocessing
def multiprocessing_run(cmdlines, verbose=False):
    """Distributes passed command-line jobs using multiprocessing.

    - cmdlines - an iterable of command line strings
    - logger - an optional Python logger module Logger instance
    """
    pool = multiprocessing.Pool()
    completed = []
    if verbose:
        callback_fn = status_callback
    else:
        callback_fn = completed.append
    pool_outputs = [pool.apply_async(subprocess.call,
                                     (str(cline), ),
                                     {'stderr': subprocess.PIPE,
                                      'shell': sys.platform != "win32"},
                                     callback=callback_fn)
                    for cline in cmdlines]
    pool.close()        # Run jobs
    pool.join()         # Collect output


# Callback function with multiprocessing run status
def status_callback(val):
    """Basic callback for multiprocessing.

    - val - return status indicated from multiprocessing
    """
    pass  # We'll do something with this later
