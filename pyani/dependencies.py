# -*- coding: utf-8 -*-
# (c) University of Strathclyde 2020
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
# Copyright (c) 2020 University of Strathclyde
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
"""Module providing functions useful for reporting on dependencies."""

import pkg_resources

from typing import Generator

from .anib import get_version as get_blast_version
from .aniblastall import get_version as get_blastall_version
from .anim import get_version as get_nucmer_version

# Dict of dependencies that support pyani
REQUIREMENTS = [
    "biopython",
    "matplotlib",
    "namedlist",
    "networkx",
    "numpy",
    "openpyxl",
    "pandas",
    "Pillow",
    "scipy",
    "seaborn",
    "sqlalchemy",
    "tqdm",
    "alembic",
]

DEVELOPMENT = [
    "bandit",
    "black",
    "codecov",
    "coverage",
    "doc8",
    "flake8",
    "jinja2",
    "mypy",
    "pydocstyle",
    "pylint",
    "pytest",
    "pytest-cov",
    "sphinx",
]

PIP = ["pre-commit", "pytest-ordering", "sphinx-rtd-theme"]


def get_versions(deplist: str) -> Generator:
    """Yield package name and version for contents of REQUIREMENTS."""
    depdict = {"REQUIREMENTS": REQUIREMENTS, "DEVELOPMENT": DEVELOPMENT, "PIP": PIP}
    for depname in sorted(depdict[deplist]):
        try:
            version = pkg_resources.get_distribution(depname).version
            loc = pkg_resources.get_distribution(depname).location
        except pkg_resources.DistributionNotFound:
            version = "Not Installed"
            loc = "-"
        yield (depname, version, loc)


def get_requirements() -> Generator:
    """Yield package name and version for Python requirements."""
    return get_versions("REQUIREMENTS")


def get_dev_requirements() -> Generator:
    """Yield package name and version for Python developer requirements."""
    return get_versions("DEVELOPMENT")


def get_pip_requirements() -> Generator:
    """Yield package name and version for Python pip installed requirements."""
    return get_versions("PIP")


def get_tool_versions() -> Generator:
    """Yield tool name and version for third-party installed tools."""
    for name, func in [
        ("blast+", get_blast_version),
        ("nucmer", get_nucmer_version),
        ("blastall", get_blastall_version),
    ]:
        yield (name, func())
