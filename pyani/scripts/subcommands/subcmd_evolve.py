import logging

from argparse import (
    Namespace,
)

# from pyani import (
#     evolve,
#     pyani_config,
# )

# from pyani.evolve import MutationEvent, MutatableRecord, PyaniEvolveException

from pyani.pyani_tools import termcolor


def subcmd_evolve(args: Namespace) -> None:
    """Produce a set of input files for pyani with a specific topology and
    known mutational relationships.

    :param args:  Namespace, command-line arguments


    """
    logger = logging.getLogger(__name__)

    # announce that we're starting
    logger.info(termcolor("Running evolve", "red"))

    raise NotImplementedError
