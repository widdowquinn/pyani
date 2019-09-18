.. _pyani-createdb:

===================================
Creating a Local ``pyani`` Database
===================================

``pyani`` stores genome information and analysis results in a persistent local `SQLite3`_ database. This allows for reuse of previous comparisons and reanalysis of datasets without having to rerun the analysis. It also means that the genome comparison results don't have to be stored in full on disk, saving space.

.. NOTE::
    To conduct a ``pyani`` analysis, there needs to be an existing database in-place.

To create a new, empty database you can use the ``pyani createdb`` command.

For more information about the ``pyani createdb`` subcommand, please see the :ref:`pyani-subcmd-createdb` page, or issue the command ``pyani download -h`` to see the inline help.

-------------------------------------
Create a new empty ``pyani`` database
-------------------------------------

The basic form of the command is:

.. code-block:: bash

    pyani createdb

This instructs ``pyani`` to create a new, empty database for analysis at the default location.

.. NOTE::
    The default location for the ``pyani`` database is in a hidden directory: ``.pyani/pyanidb``. All other ``pyani`` subcommands will look in this location for the database, unless told otherwise using the ``--dbpath`` option.

For example:

.. code-block:: bash

    $ls .pyani
    ls: .pyani: No such file or directory
    $ pyani createdb
    $ ls .pyani
    pyanidb

---------------------------------------------------------
Create an empty ``pyani`` database at a specific location
---------------------------------------------------------

.. TIP::
    If you use ``pyani`` for a number of distinct taxa, it can be convenient to create a new database for each project, to avoid performance issues as the database grows in size, filled by data that does not contribute to the analysis.

The following command can be used to specify the location of the newly-created ``pyani`` database:

.. code-block:: bash

    pyani createdb --dbpath <PATH_TO_DATABASE>

where ``<PATH_TO_DATABASE>`` is the intended location of the database. For instance, to create a new database specific for an analysis we'll call ``multitaxa``, we could use the command:

.. code-block:: bash

    $ ls .pyani
    pyanidb
    $ pyani createdb --dbpath .pyani/multitaxadb
    $ ls .pyani
    multitaxadb	pyanidb

The new database can then be specified in other ``pyani`` subcommands, using the ``--dbpath`` option.


.. _SQLite3: https://www.sqlite.org/index.html
