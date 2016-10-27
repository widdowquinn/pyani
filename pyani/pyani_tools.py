# Copyright 2016, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to support pyani."""


# Read sequence annotations in from file
def get_labels(filename, logger=None):
    """Returns a dictionary of alternative sequence labels, or None

    - filename - path to file containing tab-separated table of labels

    Input files should be formatted as <key>\t<label>, one pair per line.
    """
    labeldict = {}
    if filename is not None:
        if logger:
            logger.info("Reading labels from %s", filename)
        with open(filename, 'rU') as fh:
            count = 0
            for line in fh.readlines():
                count += 1
                try:
                    key, label = line.strip().split('\t')
                except ValueError:
                    if logger:
                        logger.warning("Problem with class file: %s",
                                       filename)
                        logger.warning("%d: %s", (count, line.strip()))
                        logger.warning("(skipping line)")
                    continue
                else:
                    labeldict[key] = label
    return labeldict
