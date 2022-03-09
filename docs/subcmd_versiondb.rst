. _pyani-subcmd-versiondb:

===============
pyani versiondb
===============

The ``versiondb`` subcommand will migrate a ``pyani`` database between versions. `fastani_merge`_

usage: pyani versiondb [-h] [-l LOGFILE] [-v] [--debug] [--disable_tqdm] [--version] [--citation] [--dbpath DBPATH]
                       (--upgrade [VERSION] | --downgrade [VERSION] | --dry-run DIRECTION START:END) [--alembic_exe ALEMBIC_EXE] [-n NAME] [-c FILE]

One of --upgrade, --downgrade, or --dry-run must be specified.

optional arguments:
  -h, --help            show this help message and exit
  -l LOGFILE, --logfile LOGFILE
                        logfile location (default: None)
  -v, --verbose         report verbose progress to log (default: False)
  --debug               report debug messages to log (default: False)
  --disable_tqdm        Turn off tqdm progress bar (default: False)
  --version
  --citation
  --dbpath DBPATH       path to pyani database (default: .pyani/pyanidb)
  --upgrade [VERSION]   update an existing database to a newer schema; if no argument is given, 'head' will be used (default: None)
  --downgrade [VERSION]
                        revert an existing database to a older schema; if no argument is given, 'base' will be used (default: None)
  --dry-run DIRECTION START:END
                        produce the SQL that would be run in migrations, without altering the database; a direction {upgrade or downgrade} and start and end versions e.g., {head:base,
                        base:head, base:<revision_id>} must be specified (default: None)
  --alembic_exe ALEMBIC_EXE
                        path to alembic executable (default: alembic)
  -n NAME, --name NAME  used to specify an individual database in a multidb setup (default: None)
  -c FILE, --config FILE
                        used to specify a config file for alembic (default: None)


----------
References
----------

.. _SQLite: https://www.sqlite.org/docs.html

.. _fastani_merge: https://github.com/widdowquinn/pyani/pull/299/commits/254346cae24058b745bd9496b4205400da03fb4c
