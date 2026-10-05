from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Literal

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

from pavvis._validation import validate_color_by
from pavvis.plots._helpers import _resolve_gene_index

pio.templates.default = "simple_white"

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix

# Distinct colours for metadata annotation strips (colourblind-friendly palette).
_STRIP_COLOURS = [
    "#2a9d8f",  # teal
    "#e9c46a",  # yellow
    "#e76f51",  # orange-red
    "#6a4c93",  # indigo
    "#264653",  # dark teal
    "#f4a261",  # amber
    "#a8dadc",  # light blue
    "#457b9d",  # steel blue
]


def _cluster_order(distance_matrix: np.ndarray) -> list[int]:
    """Return sample indices reordered by average-linkage hierarchical clustering."""
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform

    condensed = squareform(distance_matrix, checks=False)
    Z = linkage(condensed, method="average")
    return leaves_list(Z).tolist()


def _require_jaccard(pm: PavMatrix) -> None:
    if pm._jaccard_matrix is None:
        raise RuntimeError(
            "No Jaccard matrix found. Call pm.compute_jaccard() before using "
            "col_cluster=True or cluster=True."
        )


def pav_heatmap(
    pm: PavMatrix,
    gene_set: Literal["all", "variable", "core", "soft_core", "dispensable", "private"] = "variable",
    col_cluster: bool = False,
    color_by: str | None = None,
    show_gene_labels: bool = False,
) -> go.Figure:
    """Plot the PAV matrix as a binary heatmap (genes × samples).

    Each cell is coloured by absence (0) or presence (1). Genes are rows,
    samples are columns. An optional metadata annotation strip can be added
    above the heatmap to identify sample groups at a glance.

    Args:
        pm: A PavMatrix instance.
        gene_set: Subset of genes to display. Defaults to 'variable' to keep
            the plot readable on large matrices. Accepts 'all', 'variable',
            'core', 'soft_core', 'dispensable', or 'private'.
        col_cluster: If True, reorder sample columns by hierarchical clustering
            on Jaccard distance. Requires pm.compute_jaccard() to have been
            called first. Defaults to False.
        color_by: An optional discrete metadata column. When provided, a
            colour-coded annotation strip is drawn above the heatmap with one
            cell per sample coloured by group. Defaults to None.
        show_gene_labels: If True, display gene identifiers on the y-axis.
            Defaults to False. Only recommended for small gene sets as labels
            become unreadable at scale.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        RuntimeError: If col_cluster=True and pm.compute_jaccard() has not been called.
        ValueError: If color_by is not a discrete metadata column.
    """
    if col_cluster:
        _require_jaccard(pm)
    if color_by is not None:
        validate_color_by(color_by, pm.metadata, expected_type="discrete")

    gene_index = _resolve_gene_index(pm, gene_set)
    n_genes = len(gene_index)
    n_samples = len(pm.sample_ids)
    if n_genes > 5000:
        warnings.warn(
            f"pav_heatmap: {n_genes} genes selected. Rendering may be slow and the plot "
            f"difficult to read. Consider a more specific gene_set (e.g. 'soft_core' or "
            f"'dispensable') to reduce the number of rows.",
            UserWarning, stacklevel=2,
        )
    if n_samples > 200:
        warnings.warn(
            f"pav_heatmap: {n_samples} samples. Rendering may be slow and individual "
            f"columns difficult to distinguish at this scale.",
            UserWarning, stacklevel=2,
        )
    pav_sub = pm.pav.loc[gene_index]  # (genes, samples)

    sample_order = list(pav_sub.columns)
    if col_cluster:
        dist_matrix = 1.0 - pm._jaccard_matrix.values
        order_idx = _cluster_order(dist_matrix)
        sample_order = [pav_sub.columns[i] for i in order_idx]

    pav_sub = pav_sub[sample_order]
    genes = list(pav_sub.index)

    if color_by is not None:
        fig = make_subplots(
            rows=2, cols=1,
            row_heights=[0.06, 0.94],
            shared_xaxes=True,
            vertical_spacing=0.01,
        )

        col_vals = pm.metadata.loc[sample_order, color_by].astype(str)
        categories = sorted(col_vals.dropna().unique())
        cat_to_int = {c: i for i, c in enumerate(categories)}
        strip_z = [[cat_to_int[v] for v in col_vals]]
        colour_scale = [
            [i / max(len(categories) - 1, 1), _STRIP_COLOURS[i % len(_STRIP_COLOURS)]]
            for i in range(len(categories))
        ]

        fig.add_trace(go.Heatmap(
            z=strip_z,
            x=sample_order,
            y=[color_by],
            colorscale=colour_scale,
            showscale=False,
            customdata=[[[v] for v in col_vals]],
            hovertemplate="%{x}<br>" + color_by + ": %{customdata[0]}<extra></extra>",
            zmin=0,
            zmax=len(categories) - 1,
        ), row=1, col=1)

        fig.add_trace(go.Heatmap(
            z=pav_sub.values.tolist(),
            x=sample_order,
            y=genes,
            colorscale=[[0, "#f0f0f0"], [1, "#1d3557"]],
            showscale=False,
            zmin=0, zmax=1,
            hovertemplate="Gene: %{y}<br>Sample: %{x}<br>Present: %{z}<extra></extra>",
        ), row=2, col=1)

        fig.update_yaxes(showticklabels=show_gene_labels, row=2, col=1)
        fig.update_layout(height=600)
    else:
        fig = go.Figure(go.Heatmap(
            z=pav_sub.values.tolist(),
            x=sample_order,
            y=genes,
            colorscale=[[0, "#f0f0f0"], [1, "#1d3557"]],
            showscale=False,
            zmin=0, zmax=1,
            hovertemplate="Gene: %{y}<br>Sample: %{x}<br>Present: %{z}<extra></extra>",
        ))
        fig.update_yaxes(showticklabels=show_gene_labels)

    fig.update_layout(
        title=f"PAV Matrix — {gene_set.replace('_', ' ').title()} genes",
    )
    if color_by is not None:
        fig.update_xaxes(title_text="Samples", row=2, col=1)
    else:
        fig.update_xaxes(title_text="Samples")
    return fig



def sample_similarity_heatmap(
    pm: PavMatrix,
    color_by: str | None = None,
    cluster: bool = True,
) -> go.Figure:
    """Plot a symmetric sample × sample Jaccard similarity heatmap.

    Each cell shows the pairwise Jaccard similarity between two samples computed
    over all genes. Values range from 0 (no shared genes) to 1 (identical PAV
    profiles). An optional annotation strip can be added to both axes to identify
    group membership.

    Args:
        pm: A PavMatrix instance with a cached Jaccard matrix
            (call pm.compute_jaccard() first).
        color_by: An optional discrete metadata column. When provided, a
            colour-coded annotation strip is drawn along both axes.
            Defaults to None.
        cluster: If True, reorder samples by hierarchical clustering on Jaccard
            distance. Defaults to True.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        RuntimeError: If pm.compute_jaccard() has not been called.
        ValueError: If color_by is not a discrete metadata column.
    """
    _require_jaccard(pm)
    if color_by is not None:
        validate_color_by(color_by, pm.metadata, expected_type="discrete")

    n_samples = len(pm.sample_ids)
    if n_samples > 200:
        warnings.warn(
            f"sample_similarity_heatmap: {n_samples} samples. The similarity matrix "
            f"will be {n_samples}×{n_samples} and may be slow to render with individual "
            f"cells difficult to distinguish at this scale.",
            UserWarning, stacklevel=2,
        )

    sample_order = list(pm._jaccard_matrix.index)
    if cluster:
        dist_matrix = 1.0 - pm._jaccard_matrix.values
        order_idx = _cluster_order(dist_matrix)
        sample_order = [pm._jaccard_matrix.index[i] for i in order_idx]

    sim = pm._jaccard_matrix.loc[sample_order, sample_order]

    if color_by is not None:
        col_vals = pm.metadata.loc[sample_order, color_by].astype(str)
        categories = sorted(col_vals.dropna().unique())
        cat_to_int = {c: i for i, c in enumerate(categories)}
        strip_z = [[cat_to_int[v] for v in col_vals]]
        colour_scale = [
            [i / max(len(categories) - 1, 1), _STRIP_COLOURS[i % len(_STRIP_COLOURS)]]
            for i in range(len(categories))
        ]

        fig = make_subplots(
            rows=2, cols=2,
            row_heights=[0.06, 0.94],
            column_widths=[0.06, 0.94],
            shared_xaxes="columns",
            shared_yaxes="rows",
            vertical_spacing=0.01,
            horizontal_spacing=0.01,
        )

        # Top strip (column annotation)
        fig.add_trace(go.Heatmap(
            z=strip_z,
            x=sample_order,
            y=[color_by],
            colorscale=colour_scale,
            showscale=False,
            zmin=0, zmax=len(categories) - 1,
            customdata=[[[v] for v in col_vals]],
            hovertemplate="%{x}<br>" + color_by + ": %{customdata[0]}<extra></extra>",
        ), row=1, col=2)

        # Left strip (row annotation)
        fig.add_trace(go.Heatmap(
            z=[[cat_to_int[v]] for v in col_vals],
            x=[color_by],
            y=sample_order,
            colorscale=colour_scale,
            showscale=False,
            zmin=0, zmax=len(categories) - 1,
            customdata=[[[v]] for v in col_vals],
            hovertemplate="%{y}<br>" + color_by + ": %{customdata[0]}<extra></extra>",
        ), row=2, col=1)

        # Main heatmap
        fig.add_trace(go.Heatmap(
            z=sim.values.tolist(),
            x=sample_order,
            y=sample_order,
            colorscale="Blues",
            zmin=0, zmax=1,
            colorbar=dict(title="Jaccard similarity"),
            hovertemplate="x: %{x}<br>y: %{y}<br>Jaccard: %{z:.3f}<extra></extra>",
        ), row=2, col=2)

        fig.update_xaxes(showticklabels=False, row=2, col=2)
        fig.update_yaxes(showticklabels=False, row=2, col=2)
        fig.update_layout(height=650, width=700)
    else:
        fig = go.Figure(go.Heatmap(
            z=sim.values.tolist(),
            x=sample_order,
            y=sample_order,
            colorscale="Blues",
            zmin=0, zmax=1,
            colorbar=dict(title="Jaccard similarity"),
            hovertemplate="x: %{x}<br>y: %{y}<br>Jaccard: %{z:.3f}<extra></extra>",
        ))
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

    fig.update_layout(title="Sample Jaccard Similarity")
    return fig
