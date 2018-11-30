# python package version
# should match r"^__version__ = '(?P<version>[^']+)'$" for setup.py
"""Module with main code for pyani application/package."""
__version__ = "0.3.0-alpha"

import sys
import traceback


# General exception for scripts
class PyaniException(Exception):
    """General exception for pyani."""

    def __init__(self, msg="Error in pyani package"):
        """Instantiate class."""
        Exception.__init__(self, msg)


# Report last exception as string
def last_exception():
    """Return last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
