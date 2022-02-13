"""add fastani columns

Revision ID: 92f7f6b1626e
Revises:
Create Date: 2022-02-07 22:57:35.779356

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = "92f7f6b1626e"
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("comparisons", sa.Column("kmersize", sa.Integer))
    op.add_column("comparisons", sa.Column("minmatch", sa.Float))
    # with op.batch_alter_table("comparisons") as batch_op:
    # batch_op.add_column(sa.Column('kmersize', sa.Integer))
    # batch_op.add_column(sa.Column('minmatch', sa.Float))


def downgrade():
    op.drop_column("comparisons", "kmersize")
    op.drop_column("comparisons", "minmatch")
    # with op.batch_alter_table("comparisons") as batch_op:
    # batch_op.drop_column('kmersize')
    # batch_op.drop_column('minmatch')
