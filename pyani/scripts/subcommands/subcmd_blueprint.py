# -*- coding: utf-8 -*-
# (c) The University of Strathclyde 2021–Present
# Author: Bailey Harrington
#
# Contact: bailey.harrington@strath.ac.uk
"""Provides the blueprint subcommand for pyani."""

import logging

from argparse import Namespace


def subcmd_blueprint(args: Namespace) -> int:
    """Does something new.

    :param args:  Namespace, received command-line arguments
    :param logger:  logging object

    A more detailed explanation of what happens.
    """
    logger = logging.getLogger(__name__)

    logger.info("Starting `pyani blueprint`.")
    # Everything else that happens is highly dependent upon the
    # subcommand's purpose

    return 0


# Other functions may be defined down here, as helpers, perhaps.
