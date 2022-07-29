from logging.config import fileConfig

from sqlalchemy import engine_from_config
from sqlalchemy import pool

from alembic import context

import os

# this is the Alembic Config object, which provides
# access to the values within the .ini file in use.
config = context.config

dbpath = os.environ.get("PYANI_DATABASE")
script_dir = os.environ.get("ALEMBIC_MIGRATIONS_DIR")
url = f"sqlite:///{dbpath}"
config.set_main_option("sqlalchemy.url", url)
config.set_main_option("script_location", "alembic")


# Interpret the config file for Python logging.
# This line sets up loggers basically.
fileConfig(config.config_file_name)

# add your model's MetaData object here
# for 'autogenerate' support
# from myapp import mymodel
# target_metadata = mymodel.Base.metadata
target_metadata = None

# other values from the config, defined by the needs of env.py,
# can be acquired:
# my_important_option = config.get_main_option("my_important_option")
# ... etc.


def run_migrations_offline():
    """Run migrations in 'offline' mode.

    This configures the context with just a URL
    and not an Engine, though an Engine is acceptable
    here as well.  By skipping the Engine creation
    we don't even need a DBAPI to be available.

    Calls to context.execute() here emit the given string to the
    script output.

    """

    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
    )

    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online():
    """Run migrations in 'online' mode.

    In this scenario we need to create an Engine
    and associate a connection with the context.

    """

    connectable = engine_from_config(
        config.get_section(config.config_ini_section),
        prefix="sqlalchemy.",
        poolclass=pool.NullPool,
    )

    with connectable.connect() as connection:
        context.configure(connection=connection, target_metadata=target_metadata)

        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    # Track the current version of the database in a local file
    # this value is used in the --dry-run option
    version_file = os.path.join(
        os.path.dirname(config.config_file_name), "alembic/version.txt"
    )
    if os.path.exists(version_file):
        current_version = open(version_file).read().strip()
    else:
        current_version = None
    context.configure(dialect_name="sqlite", starting_rev=current_version)

    # Perform the dry run
    run_migrations_offline()

    # Write 'new' version to file
    end_version = context.get_revision_argument()
    if end_version and end_version != current_version:
        open(version_file, "w").write(end_version)
    elif end_version is None:
        open(version_file, "w").write("base")
else:
    run_migrations_online()
