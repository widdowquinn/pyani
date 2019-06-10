# python package version
# should match r"^__version__ = '(?P<version>[^']+)'$" for setup.py
"""Module with main code for pyani application/package."""
__version__ = "0.3.0-alpha"


# General exception for scripts
class PyaniException(Exception):
    """General exception for pyani."""
