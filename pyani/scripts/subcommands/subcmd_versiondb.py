import logging

from argparse import Namespace

from pyani import pyani_orm


def subcmd_versiondb(args: Namespace) -> int:
    """Create an empty pyani database.

    :param args:  Namespace, command-line arguments
    :param logger:  logging object
    """
    # Create logger
    logger = logging.getLogger(__name__)

    # If the database exists, raise an error rather than overwrite
    if not args.dbpath.is_file():
        logger.error("Database %s does not exist (exiting)", args.dbpath)
        raise SystemError(1)

    # If the path to the database doesn't exist, create it
    if args.upgrade:
        logger.info("Upgrading database schema to: %s", args.upgrade)
    elif args.downgrade:
        logger.info("Downgrading database schema to: %s", args.downgrade)

    # Do some stuff in alembic

    return 0
