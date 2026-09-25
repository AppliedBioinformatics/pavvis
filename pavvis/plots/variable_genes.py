from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio


pio.templates.default = "simple_white"

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix

# Category display order, labels, and colors (no red/green)
_CATEGORIES = [
    ("absent",      "Absent",      "#8d99ae"),  # slate grey
    ("private",     "Private",     "#f4a261"),  # amber
    ("dispensable", "Dispensable", "#2a9d8f"),  # teal
    ("soft_core",   "Soft Core",   "#6a4c93"),  # indigo
    ("core",        "Core",        "#1d3557"),  # navy
]


def frequency_histogram(
    pm: PavMatrix,
    n_bins: int = 20,
    exclude_absent: bool = False,
    exclude_core: bool = False,
    log_y: bool = False,
) -> go.Figure:
    """Plot gene presence frequency as a stacked histogram coloured by gene category.

    Each bin is subdivided by category — absent, private, dispensable, soft core,
    and core — so the composition of each frequency range is visible at a glance.
    Categories are determined by the thresholds stored on the PavMatrix instance.

    Args:
        pm: A PavMatrix instance.
        n_bins: Number of equal-width bins across [0, 1]. Defaults to 20.
        exclude_absent: If True, absent genes (frequency = 0) are not shown.
            Defaults to False.
        exclude_core: If True, core genes (frequency = 1) are not shown.
            Defaults to False.
        log_y: If True, display the y-axis on a log10 scale. Defaults to False.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.
    """
    freq = pm.presence_frequency

    # Assign each gene to a category
    def _categorise(f: float) -> str:
        if f == 0.0:
            return "absent"
        if f == 1.0:
            return "core"
        if f >= pm.soft_core_min:
            return "soft_core"
        if f >= pm.dispensable_min:
            return "dispensable"
        return "private"

    categories = freq.map(_categorise)

    # Compute shared bin edges
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_centres = ((bin_edges[:-1] + bin_edges[1:]) / 2).tolist()
    bin_width = bin_edges[1] - bin_edges[0]

    traces = []
    for key, label, color in _CATEGORIES:
        if key == "absent" and exclude_absent:
            continue
        if key == "core" and exclude_core:
            continue

        gene_freqs = freq[categories == key].values
        counts, _ = np.histogram(gene_freqs, bins=bin_edges)

        traces.append(go.Bar(
            x=bin_centres,
            y=counts.tolist(),
            name=label,
            marker_color=color,
            width=bin_width,
            hovertemplate=f"{label}<br>Frequency: %{{x:.2f}}<br>Genes: %{{y}}<extra></extra>",
        ))

    fig = go.Figure(traces)
    fig.update_layout(
        barmode="stack",
        title="Gene Presence Frequency Distribution",
        xaxis_title="Presence Frequency",
        yaxis_title="Number of Genes",
        xaxis=dict(range=[0, 1]),
        yaxis=dict(type="log" if log_y else "linear"),
        legend_title="Category",
    )

    return fig


def exclusive_genes_bar(
    pm: PavMatrix,
    column: str,
    log_y: bool = False,
) -> go.Figure:
    """Plot the number of exclusive variable genes per group in a metadata column.

    An exclusive gene is a variable gene that is present in at least one sample
    within a group and absent from every sample outside that group. Each bar
    represents one group; the height is its exclusive gene count.

    Only discrete metadata columns are supported.

    Args:
        pm: A PavMatrix instance.
        column: A discrete metadata column to group samples by.
        log_y: If True, display the y-axis on a log10 scale. Defaults to False.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        ValueError: If column is not a discrete metadata column.
    """
    counts = pm.exclusive_gene_counts(column)
    categories = list(counts.keys())
    exclusive_counts = list(counts.values())

    fig = go.Figure(go.Bar(
        x=categories,
        y=exclusive_counts,
        marker_color="#2a9d8f",
        hovertemplate="%{x}<br>Exclusive genes: %{y}<extra></extra>",
    ))
    fig.update_layout(
        title=f"Exclusive Variable Genes by {column}",
        xaxis_title=column,
        yaxis_title="Number of Exclusive Genes",
        yaxis=dict(type="log" if log_y else "linear"),
        showlegend=False,
    )

    return fig