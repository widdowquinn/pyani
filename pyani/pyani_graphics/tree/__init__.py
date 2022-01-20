import logging
from pyani import pyani_graphics
from scipy.cluster import hierarchy
from ete3 import ClusterTree, Tree, TreeStyle, faces, AttrFace, PhyloTree
from pathlib import Path
import sys
import seaborn as sns

LABEL_DICT = {}


def build_label_dict(fig, axis, params):
    """Label info for tree plots.

    :param fig:  a Seaborn clustermap instance
    :param axis:  one of {'row', 'col'}
    :param params:  plot parameters; this is where the labels come from

    """
    logger = logging.getLogger(__name__)
    if axis == "col":
        for idx, _ in zip(
            fig.dendrogram_col.reordered_ind, fig.ax_heatmap.get_yticklabels()
        ):
            LABEL_DICT[str(idx + 1)] = params.labels.get(_, _.get_text())
    elif axis == "row":
        for idx, _ in zip(
            fig.dendrogram_row.reordered_ind, fig.ax_heatmap.get_xticklabels()
        ):
            LABEL_DICT[str(idx + 1)] = params.labels.get(_, _.get_text())
    logger.debug(f"Label dict: {LABEL_DICT}")
    return LABEL_DICT


def get_newick(node, parentdist, leaf_names, newick=""):
    """Generates a newick formatted file from a tree,
    using recursion to traverse it.

    :param node:  a (portion of a) tree to be traversed
    :param parentdist:  distance from the parent node
    :param leaf_names:  lables that will be attached to the terminal nodes
    :param newick:  the current newick-formatted tree structure

    """
    # logger = logging.getLogger(__name__)
    # logger.debug(f"{type(parentdist)}, {parentdist}")
    # logger.debug(f"{type(node.dist)}, {node.dist}")
    diff = parentdist - node.dist
    if node.is_leaf():
        return f"{leaf_names[node.id]}:{diff:.2f}{newick}"
    else:
        if len(newick) > 0:
            newick = f"):{diff:.2f}{newick}"
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, f",{newick}")
        newick = f"({newick}"
        return newick


def tree(dfr, outfname, title, params, format, args):
    """Generate a newick file and dendrogram plot for the given dataframe.

    :param dfr:  a dataframe
    # :param fig:  a figure produced by sns.clustermap
    :param title:  name of the matrix plot
    :param format:  image file format being used
    :param params:  matrix plot parameters; including labels
    :param args:  Namespace

    """
    logger = logging.getLogger(__name__)

    # Get matrix name and run_id from the plot title
    matname, run_id = title.split("_", 1)[-1].rsplit("_", 1)

    maxfigsize = 120
    calcfigsize = dfr.shape[0] * 1.1
    figsize = min(max(8, calcfigsize), maxfigsize)
    if figsize == maxfigsize:
        scale = maxfigsize / calcfigsize
        sns.set_context("notebook", font_scale=scale)

    # Add a colorbar?
    if params.classes is None:
        col_cb = None
    else:
        col_cb = pyani_graphics.sns.get_colorbar(dfr, params.classes)

    params.colorbar = col_cb
    params.figsize = figsize
    params.linewidths = 0.25

    fig = pyani_graphics.sns.get_clustermap(dfr, params)

    # Dictionary to allow abstraction over axes
    sides = {
        "columns": {
            "axis": fig.dendrogram_col,
            "names": dfr.columns,  # fig.dendrogram_col.reordered_ind,
        },
        "rows": {
            "axis": fig.dendrogram_row,
            "names": dfr.index,  # fig.dendrogram_row.reordered_ind,
        },
    }

    # Create a linkage dendrogram and newick string for both rows and columns
    newicks = {}

    for axis in args.axes:
        # Generate newick format
        tree = hierarchy.to_tree(sides[axis]["axis"].linkage, False)
        logger.debug(f"Names: {sides[axis]['names']}")

        newick = get_newick(tree, tree.dist, sides[axis]["names"], "")
        newicks.update({f"[{axis}_newick_{matname}_{run_id}]": newick})

        # Generate dendrogram
        # if 'dendrogram' in args.tree:
        # if args.tree:
        build_label_dict(fig, axis, params)
        sys.stderr.write(f"Label dict: {LABEL_DICT}\n")
        # figtree = ClusterTree(newick, text_array=matrix)
        figtree = PhyloTree(newick)
        figtree.set_species_naming_function(get_species_name)
        figtree_file = Path(args.outdir) / f"{axis}_tree_{matname}_{run_id}.{format}"
        logger.debug(f"{figtree}")

        # Write the tree to file
        figtree.render(str(figtree_file), layout=tree_layout)

    # Return the newick strings so we can save them in the database (eventually)
    return newicks


def tree_layout(node):

    # Add taxonomy to nodes, and align to right
    if node.is_leaf():
        # if node.name == "F962_00589":
        #     faces.add_face_to_node(
        #         AttrFace("name", fgcolor="white"),
        #         node,
        #         column=0,
        #         position="branch-right",
        #     )
        #     faces.add_face_to_node(
        #         AttrFace("species", fgcolor="white"), node, column=0, position="aligned"
        #     )
        #     node.img_style["bgcolor"] == "darkred"
        # else:

        faces.add_face_to_node(
            AttrFace("name", fgcolor="black"),
            node,
            column=0,
            position="branch-right",
        )
        faces.add_face_to_node(AttrFace("species"), node, column=0, position="aligned")


def get_species_name(node_name_string):
    """Return `Genus species` (where known) for a node."""
    return LABEL_DICT[node_name_string]
