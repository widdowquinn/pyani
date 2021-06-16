# -*- coding: utf-8 -*-
# (c) University of Strathclyde 2020
# Author: Leighton Pritchard
#
# Contact: leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G4 0RE
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
"""Code to implement the ANIblastall average nucleotide identity method."""

import os
import platform
import re
import shutil
import subprocess

from pathlib import Path

from . import pyani_config


def get_version(blast_exe: Path = pyani_config.BLASTALL_DEFAULT) -> str:
    r"""Return BLAST blastall version as a string.

    :param blast_exe:  path to blastall executable

    We expect blastall to return a string as, for example

    .. code-block:: bash

        $ blastall -version
        [blastall 2.2.26] ERROR: Number of database sequences to show \
        one-line descriptions for (V) [ersion] is bad or out of range [? to ?]

    This is concatenated with the OS name.
    """
    blastall_path = Path(shutil.which(blast_exe))  # type:ignore

    if not os.path.isfile(blastall_path):  # no executable
        return f"No blastall at {blastall_path}"

    if not os.access(blastall_path, os.X_OK):  # file exists but not executable
        return f"blastall exists at {blastall_path} but not executable"

    if platform.system() == "Darwin":
        cmdline = [blast_exe, "-version"]
    else:
        cmdline = [blast_exe]

    try:
        result = subprocess.run(
            cmdline,  # type: ignore
            shell=False,
            stdout=subprocess.PIPE,  # type: ignore
            stderr=subprocess.PIPE,
            check=False,  # blastall doesn't return 0
        )
        version = re.search(  # type: ignore
            r"(?<=blastall\s)[0-9\.]*", str(result.stderr, "utf-8")
        ).group()
    except OSError:
        return f"blastall exists at {blastall_path} but could not be executed"

    if 0 == len(version.strip()):
        return f"blastall exists at {blastall_path} but could not retrieve version"

    return f"{platform.system()}_{version} ({blastall_path})"
