from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import plotly.graph_objects as go
import plotly.io as pio

pio.templates.default = "simple_white"

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix


def presence_frequency(
    pm: PavMatrix,
    type: Literal["histogram", "bar"] = "histogram",
    n_bins: int | None = None,
) -> go.Figure:
    """Plot the presence frequency of genes across samples.

    Args:
        pm: A PavMatrix instance.
        type: Chart type — 'histogram' shows the distribution of frequencies
            across all genes; 'bar' shows one bar per gene sorted by frequency.
        n_bins: For 'histogram' type only, number of bins. Defaults to Plotly's
            automatic binning.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.
    """
    freq = pm.presence_frequency

    if type == "histogram":
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

    elif type == "bar":
        freq_sorted = freq.sort_values(ascending=False)
        fig = go.Figure(
            go.Bar(
                x=freq_sorted.index.tolist(),
                y=freq_sorted.values,
                name="Presence Frequency",
            )
        )
        fig.update_layout(
            title="Gene Presence Frequency",
            xaxis_title="Gene",
            yaxis_title="Presence Frequency",
            yaxis=dict(range=[0, 1]),
        )

    else:
        raise ValueError(f"Unknown plot type '{type}'. Expected 'histogram' or 'bar'.")

    return fig