.. _pyani-run_versiondb:

=================
Running versiondb
=================

``pyani`` allows for migration of databases between versions using `versiondb`_(*versiondb*), backed by `Alembic`_. To run versiondb on an existing database, use the ``pyani versiondb`` subcommand.

In brief, ``versiondb`` can upgrade or downgrade a database; it also has an option to perform a dry run:

1. ``--upgrade`` is used to upgrade a database to a more current version. If the flag is passed, but no value is specified, this defaults to `head` (the most up-to-date version). A backup of the original database is created prior to the upgrade.
2. ``--downgrade`` is used to downgrade a database to an older version. If the flag is passed, but no value is specified, this defaults to `base` (the oldest version). A backup of the original database is created prior to the downgrade.
3. ``--dry-run`` is used to show what changes would be made to the specified database, without actually performing any alterations. An SQL script containing the commands that would be run is output.


.. ATTENTION::
    ``pyani versiondb`` requires that a working copy of `Alembic`_ is available. This should happen automatically with new ``pyani`` installs. If it is missing, please see `pyani-installation`_ for information about installing this package.

For more information about the ``pyani versiondb`` subcommand, please see the `pyani-subcmd-versiondb`_ page, or issue the command ``pyani versiondb -h`` to see the inline help.

----------------------------
Perform versiondb migrations
----------------------------

The basic form of the command is:

.. code-block:: bash

    pyani versiondb [--dbpath DATABASE] {--upgrade [VERSION], --downgrade [vVERSION], --dry-run {up, down} START:END}

One of ``--upgrade``, ``--downgrade``, or ``--dry-run`` must be provided; a migration version can be optionally specified for ``--upgrade`` and ``--downgrade``; if one is not, these will behave as ``--upgrade head`` and ``--downgrade base``, respectively. ``--dry-run`` takes two values; the first, a direction (one of 'up' or 'down'), the second the start and end database versions (given as ``START:END``), specified by their unique hashes, or 'head' and 'base', for the newest/oldest versions, respectively.

All backups and SQL script outputs from these options will be written to the same directory as the database given to ``--dbpath``.

This instructs ``pyani`` to modify the database at ``--dbpath``. For example, the following command for ``versiondb`` upgrades a database ``./test_database`` and writes a backup ``./test_database.upgrade.YYYY-MM-DD_HH-MM-SS.bak`` to the same containing directory:

.. code-block:: bash

    pyani versiondb --dbpath ./test_database --upgrade head

.. ATTENTION::
    To view the result of running ``versiondb``, you will need to use the ``sqlite3 <database>``. Please see `SQLite`_ for more details.

The expectation is that the above command is what most ``pyani`` users will need in terms of migrations, as this will perform the necessary alterations to the database to accommodate changes merged into ``pyani`` with commit `fastani_merge`_, which broke compatibility with previous database schemas.


--------------
Advanced Usage
--------------

However, there may be cases where it is necessary to downgrade a database to a specific, *earlier*, version, or the user wants to perform a dry run. Some users may also want to specify a non-standard location for an Alembic config file, or need to use this in a project with a ``multidb`` setup. Some examples of these advanced cases can be found below.

~~~~~~~~~~~~~~~~~~~
Performing a dry run
~~~~~~~~~~~~~~~~~~~

.. NOTE::

    This option has been implemented, but tests for it are still in the works.

This following command creates an SQL file, ``./test_database.downgrade.YYYY-MM-DD_HH-MM-SS.sql``, (in the same directory as ``./test_database``) containing the raw SQL that would produce the necessary changes to the database to migrate it to the specified version (in this case, downgrading it to ``base``):

.. code-block:: bash

    pyani versiondb --dbpath ./test_database --dry-run down head:base

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Using a different config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, ``Alembic``, and therefore ``pyani versiondb``, assume a config file, ``alembic.ini``, located in a standard location, which in the case of ``pyani``, is among the package files. To use a different file, add the ``-c`` or ``--config`` flag to your ``pyani`` command:

.. code-block:: bash

    pyani versiondb --dbpath ./test_database --upgrade head --config ./config.ini

If you need to specify additional settings for `Alembic`_, or have multiple databases in your ``pyani`` project (especially if not all should be upgraded/downgraded), this is the way you will need to use this option.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Using pyani versiondb in a multidb setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. NOTE::

  For information about how to set up a project with multiple databases managed by `Alembic`_, please see the `Alembic`_ documentation on `working with multiple databases <https://alembic.sqlalchemy.org/en/latest/branches.html#working-with-multiple-bases>`_.

To specify a single database in a ``multidb`` setup, use the ``-n`` or ``--name`` option, along with the designation for the correct database from your ``multidb`` config file:

.. code-block:: bash

    pyani versiondb --dbpath ./test_database --upgrade head --config ./multidb.ini --name database2

----------
References
----------

.. _SQLite: https://www.sqlite.org/docs.html

.. _fastani_merge: https://github.com/widdowquinn/pyani/pull/299/commits/254346cae24058b745bd9496b4205400da03fb4c

.. _Alembic: https://alembic.sqlalchemy.org/en/latest/
