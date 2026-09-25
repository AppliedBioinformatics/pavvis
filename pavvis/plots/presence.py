from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import plotly.graph_objects as go

from pavvis._validation import validate_color_by
from pavvis.plots._helpers import _resolve_gene_index, _gene_set_labels

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix


def frequency_histogram(
    pm: PavMatrix,
    n_bins: int | None = None,
) -> go.Figure:
    """Plot the distribution of gene presence frequencies across samples.

    Args:
        pm: A PavMatrix instance.
        n_bins: Number of histogram bins. Defaults to Plotly's automatic binning.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.
    """
    freq = pm.presence_frequency

    fig = go.Figure(
        go.Histogram(
            x=freq.values,
            xbins=dict(start=0, end=1, size=(1 / n_bins) if n_bins is not None else None),
            name="Genes",
        )
    )
    fig.update_layout(
        title="Gene Presence Frequency Distribution",
        xaxis_title="Presence Frequency",
        yaxis_title="Number of Genes (log10)",
        xaxis=dict(range=[0, 1]),
        yaxis=dict(type="log"),
    )

    return fig


def per_sample_boxplot(
    pm: PavMatrix,
    gene_set: Literal["all", "variable", "core"] = "all",
    x: str | None = None,
) -> go.Figure:
    """Plot the distribution of gene counts per sample as a box plot.

    Args:
        pm: A PavMatrix instance.
        gene_set: Which genes to count per sample — 'all', 'variable', or 'core'.
            Defaults to 'all'.
        x: A discrete metadata column to group samples by. Each unique value
            becomes a separate box. Pass None for a single box across all samples.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        ValueError: If x is not None and is not a discrete metadata column.
    """
    if x is not None:
        validate_color_by(x, pm.metadata, expected_type="discrete")

    gene_index = _resolve_gene_index(pm, gene_set)
    counts = pm.pav.loc[gene_index].sum(axis=0)
    title, yaxis_title = _gene_set_labels(gene_set)

    if x is not None:
        categories = sorted(pm.metadata[x].dropna().unique().astype(str))
        traces = []
        for cat in categories:
            samples = pm.metadata.index[pm.metadata[x].astype(str) == cat]
            traces.append(go.Box(
                y=counts.loc[samples].values,
                text=samples.tolist(),
                name=cat,
                boxpoints="outliers",
                hovertemplate="%{text}<br>Count: %{y}<extra></extra>",
            ))
        fig = go.Figure(traces)
    else:
        fig = go.Figure(go.Box(
            y=counts.values,
            text=counts.index.tolist(),
            name="All samples",
            boxpoints="outliers",
            hovertemplate="%{text}<br>Count: %{y}<extra></extra>",
        ))

    fig.update_layout(
        title=title,
        xaxis_title=x if x is not None else "",
        yaxis_title=yaxis_title,
    )

    return fig


def per_sample_scatter(
    pm: PavMatrix,
    x: str,
    gene_set: Literal["all", "variable", "core"] = "all",
    color_by: str | None = None,
) -> go.Figure:
    """Plot gene count per sample against a continuous metadata variable.

    Args:
        pm: A PavMatrix instance.
        x: A continuous metadata column to plot on the x-axis (e.g. 'Longitude').
        gene_set: Which genes to count per sample — 'all', 'variable', or 'core'.
            Defaults to 'all'.
        color_by: An optional discrete metadata column to colour points by.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        ValueError: If x is not a continuous metadata column, or color_by is not
            a discrete metadata column.
    """
    validate_color_by(x, pm.metadata, expected_type="continuous")
    if color_by is not None:
        validate_color_by(color_by, pm.metadata, expected_type="discrete")

    gene_index = _resolve_gene_index(pm, gene_set)
    counts = pm.pav.loc[gene_index].sum(axis=0)
    samples = counts.index.tolist()
    x_values = pm.metadata.loc[samples, x].tolist()
    title, yaxis_title = _gene_set_labels(gene_set)

    if color_by is not None:
        categories = sorted(pm.metadata[color_by].dropna().unique().astype(str))
        traces = []
        for cat in categories:
            cat_samples = pm.metadata.index[pm.metadata[color_by].astype(str) == cat]
            cat_samples = [s for s in samples if s in cat_samples]
            traces.append(go.Scatter(
                x=pm.metadata.loc[cat_samples, x].tolist(),
                y=counts.loc[cat_samples].values,
                text=cat_samples,
                mode="markers",
                name=cat,
                hovertemplate="%{text}<br>" + x + ": %{x}<br>Count: %{y}<extra></extra>",
            ))
        fig = go.Figure(traces)
    else:
        fig = go.Figure(go.Scatter(
            x=x_values,
            y=counts.values,
            text=samples,
            mode="markers",
            name="Samples",
            hovertemplate="%{text}<br>" + x + ": %{x}<br>Count: %{y}<extra></extra>",
        ))

    fig.update_layout(
        title=title,
        xaxis_title=x,
        yaxis_title=yaxis_title,
    )

    return fig
