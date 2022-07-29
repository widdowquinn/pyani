"""add fastani columns

Revision ID: 92f7f6b1626e
Revises:
Create Date: 2022-02-07 22:57:35.779356

"""
from alembic import op
import sqlalchemy as sa
import sys

from sqlalchemy import (
    Column,
    # Table,
    # MetaData,
    # ForeignKey,
    Integer,
    # String,
    Float,
    # Boolean,
    UniqueConstraint,
)

# revision identifiers, used by Alembic.
revision = "92f7f6b1626e"
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # op.add_column("comparisons", sa.Column("kmersize", sa.Integer))
    # op.add_column("comparisons", sa.Column("minmatch", sa.Float))
    """
    comparisons = Table(
        "comparisons",
        meta,
        Column("comparisons_id", Integer, primary_key=True),
        Column("query_id", Integer, ForeignKey("genomes.genome_id"), nullable=False),
        Column("subject_id", Integer, ForeignKey("genomes.genome_id"), nullable=False),
        Column("aln_length", Integer),
        Column("sim_errs", Integer),
        Column("identity", Float),
        Column("cov_query", Float),
        Column("cov_subject", Float),
        Column("program", String),
        Column("version", String),
        Column("fragsize", Integer),
        Column("maxmatch", Boolean),
        UniqueConstraint(
            "query_id",
            "subject_id",
            "program",
            "version",
            "fragsize",
            "maxmatch",
            name="base_reqs",
        ),
    )
    """
    with op.batch_alter_table("comparisons") as batch_op:
        # batch_op.add_column(sa.Column("kmersize", sa.Integer, default=None))
        # batch_op.add_column(sa.Column("minmatch", sa.Float, default=None))
        batch_op.drop_constraint("base_reqs")
        batch_op.add_column(sa.Column("kmersize", sa.Integer, default=None))
        batch_op.add_column(sa.Column("minmatch", sa.Float, default=None))
        batch_op.create_unique_constraint(
            "fastani_reqs",
            [
                "query_id",
                "subject_id",
                "program",
                "version",
                "fragsize",
                "maxmatch",
                "kmersize",
                "minmatch",
            ],
        )


def downgrade():
    # op.drop_constraint("comparisons", 'kmersize')
    # op.drop_column("comparisons", "kmersize")
    """
    comparisons = Table(
        "comparisons",
        meta,
        Column("comparisons_id", Integer, primary_key=True),
        Column("query_id", Integer, ForeignKey("genomes.genome_id"), nullable=False),
        Column("subject_id", Integer, ForeignKey("genomes.genome_id"), nullable=False),
        Column("aln_length", Integer),
        Column("sim_errs", Integer),
        Column("identity", Float),
        Column("cov_query", Float),
        Column("cov_subject", Float),
        Column("program", String),
        Column("version", String),
        Column("fragsize", Integer),
        Column("maxmatch", Boolean),
        Column("kmersize", Integer),
        Column("minmatch", Float),
        UniqueConstraint(
            "query_id",
            "subject_id",
            "program",
            "version",
            "fragsize",
            "maxmatch",
            "kmersize",
            "minmatch",
            name="fastani_reqs",
        ),
    )
    """
    with op.batch_alter_table("comparisons") as batch_op:
        batch_op.drop_column("kmersize")
        batch_op.drop_column("minmatch")
        batch_op.drop_constraint("fastani_reqs")
        batch_op.create_unique_constraint(
            "base_reqs",
            ["query_id", "subject_id", "program", "version", "fragsize", "maxmatch"],
        )
