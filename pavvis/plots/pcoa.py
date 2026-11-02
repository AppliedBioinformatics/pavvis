from __future__ import annotations

from typing import TYPE_CHECKING

import plotly.graph_objects as go
import plotly.io as pio

from pavvis._validation import validate_color_by
from pavvis.plots._helpers import _embedding_scatter_traces

pio.templates.default = "simple_white"

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix


def _require_pcoa(pm: PavMatrix) -> None:
    if pm._pcoa_embedding is None:
        raise RuntimeError(
            "No PCoA embedding found. Call pm.compute_pcoa() before plotting."
        )


def _axis_title(pm: PavMatrix, axis: str) -> str:
    return f"{axis} ({pm._pcoa_explained[axis] * 100:.1f}%)"


def scatter(
    pm: PavMatrix,
    color_by: str | None = None,
    components: tuple[int, int] = (1, 2),
) -> go.Figure:
    """Plot two PCoA axes of the samples as a 2D scatter plot.

    Points represent samples. Axis titles show the percentage of variance each axis
    explains. Colouring behaves identically to `pavvis.plots.umap.scatter`: discrete
    columns produce a categorical legend; continuous columns use a Viridis colorscale.

    Args:
        pm: A PavMatrix instance with a cached PCoA (call pm.compute_pcoa() first).
        color_by: An optional metadata column name to colour points by.
            Accepts both discrete and continuous columns. Defaults to None.
        components: The two 1-based PCoA axes to plot on x and y. Defaults to (1, 2).

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        RuntimeError: If pm.compute_pcoa() has not been called.
        ValueError: If color_by is not a valid metadata column, or if a component
            exceeds the number of computed axes.
    """
    _require_pcoa(pm)
    if color_by is not None:
        validate_color_by(color_by, pm.metadata)

    emb = pm._pcoa_embedding
    n_axes = emb.shape[1]
    if any(not 1 <= c <= n_axes for c in components):
        raise ValueError(
            f"components must be between 1 and {n_axes} (the number of computed axes), "
            f"got {components}. Call pm.compute_pcoa(n_components=...) to keep more axes."
        )
    xcol, ycol = (f"PCo{c}" for c in components)

    fig = go.Figure(_embedding_scatter_traces(emb, xcol, ycol, pm.metadata, color_by))
    fig.update_layout(
        title="PCoA (Jaccard distance)",
        xaxis_title=_axis_title(pm, xcol),
        yaxis_title=_axis_title(pm, ycol),
    )
    return fig


def scree(pm: PavMatrix) -> go.Figure:
    """Plot the percentage of variance explained by each computed PCoA axis.

    Percentages are relative to the sum of all positive eigenvalues, so they sum to
    less than 100% when only the leading axes are kept.

    Args:
        pm: A PavMatrix instance with a cached PCoA (call pm.compute_pcoa() first).

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        RuntimeError: If pm.compute_pcoa() has not been called.
    """
    _require_pcoa(pm)
    explained = pm._pcoa_explained * 100

    fig = go.Figure(go.Bar(
        x=explained.index.tolist(),
        y=explained.tolist(),
        marker_color="#2a9d8f",
        hovertemplate="%{x}: %{y:.2f}%<extra></extra>",
    ))
    fig.update_layout(
        title="PCoA Variance Explained",
        xaxis_title="Axis",
        yaxis_title="Variance Explained (%)",
    )
    return fig