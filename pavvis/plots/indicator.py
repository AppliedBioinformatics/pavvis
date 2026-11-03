from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

from pavvis.plots._helpers import _GROUP_COLOURS

pio.templates.default = "simple_white"

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix


def _significant_indicators(
    pm: PavMatrix, indicators: pd.DataFrame, column: str, q_max: float
) -> tuple[list, np.ndarray, pd.DataFrame]:
    """Validate an indicator_genes result against pm and keep genes with q_value < q_max."""
    groups, X = pm._group_matrix(column)

    unknown_groups = sorted(set(indicators["group"].dropna()) - set(groups))
    if unknown_groups:
        raise ValueError(
            f"indicators contains groups not found in '{column}': {unknown_groups}. "
            f"Was it computed for a different column?"
        )
    missing = [g for g in indicators.index if g not in pm.pav.index]
    if missing:
        raise ValueError(f"Genes in indicators not found in PAV matrix: {missing[:10]}")

    sig = indicators[indicators["q_value"] < q_max]
    if sig.empty:
        raise ValueError(f"No indicator genes with q_value < {q_max}.")
    return groups, X, sig


def _warn_empty_groups(empty: list, q_max: float, shown_as: str) -> None:
    if empty:
        warnings.warn(
            f"no indicator genes with q_value < {q_max} for: {empty}. "
            f"These groups are shown as {shown_as}.",
            UserWarning, stacklevel=3,
        )


def heatmap(
    pm: PavMatrix,
    indicators: pd.DataFrame,
    column: str,
    top_n: int = 10,
    q_max: float = 0.05,
    show_gene_labels: bool = True,
) -> go.Figure:
    """Plot the presence frequency of the top indicator genes in every group.

    Rows are the strongest indicator genes of each group, columns are the groups,
    and each cell is the fraction of that group's samples carrying the gene. Genes
    are grouped by the group they indicate, so a good indicator shows as a dark cell
    in its own group's column and pale cells elsewhere. The cell in a gene's own
    group equals its fidelity B. A colour strip on the left marks each row's
    indicator group, using the same colours as the other group plots.

    For each group, genes with `q_value < q_max` are ranked by `stat` (then by
    `p_value` and `B` to break ties) and the first `top_n` are shown. Groups with
    fewer significant genes show the ones they have.

    Args:
        pm: A PavMatrix instance.
        indicators: The DataFrame returned by
            `pavvis.stats.indicator.indicator_genes(pm, column, ...)`.
        column: The discrete metadata column the indicators were computed for.
        top_n: Maximum number of genes shown per group. Defaults to 10.
        q_max: Genes are shown only if their q-value is below this. Defaults to 0.05.
        show_gene_labels: If True, label rows with gene IDs. Defaults to True.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        ValueError: If column is not a discrete metadata column, if `indicators`
            contains groups not in `column` or genes not in the PAV matrix, or if no
            gene has `q_value < q_max`.

    Warns:
        UserWarning: If some groups have no gene with `q_value < q_max`.
    """
    groups, X, sig = _significant_indicators(pm, indicators, column, q_max)

    blocks, empty = [], []
    for g in groups:
        block = (sig[sig["group"] == g]
                 .sort_values(["stat", "p_value", "B"], ascending=[False, True, False])
                 .head(top_n))
        if block.empty:
            empty.append(g)
        else:
            blocks.append(block)
    _warn_empty_groups(empty, q_max, "columns only")

    rows = pd.concat(blocks)
    genes = rows.index.tolist()
    freq = (pm.pav.loc[genes].values @ X) / X.sum(axis=0)  # (genes, groups)

    group_idx = {g: i for i, g in enumerate(groups)}
    colour_scale = [
        [i / max(len(groups) - 1, 1), _GROUP_COLOURS[i % len(_GROUP_COLOURS)]]
        for i in range(len(groups))
    ]

    fig = make_subplots(rows=1, cols=2, column_widths=[0.03, 0.97], shared_yaxes=True,
                        horizontal_spacing=0.005)
    fig.add_trace(go.Heatmap(
        z=[[group_idx[g]] for g in rows["group"]], x=["indicator"], y=genes,
        colorscale=colour_scale, zmin=0, zmax=len(groups) - 1, showscale=False,
        customdata=[[[g]] for g in rows["group"]],
        hovertemplate="%{y}<br>indicator of " + column + ": %{customdata[0]}<extra></extra>",
        meta={"kind": "strip"},
    ), row=1, col=1)

    info = np.array([[g, f"{s:.3f}", f"{q:.3g}"] for g, s, q in rows[["group", "stat", "q_value"]].values],
                    dtype=object)
    fig.add_trace(go.Heatmap(
        z=freq, x=groups, y=genes,
        colorscale=[[0, "#ffffff"], [1, "#1d3557"]], zmin=0, zmax=1,
        colorbar=dict(title="Presence<br>frequency", thickness=12),
        customdata=np.repeat(info[:, None, :], len(groups), axis=1),
        hovertemplate=(
            "%{y}<br>" + column + ": %{x}<br>frequency: %{z:.2f}<br>"
            "indicator of %{customdata[0]} (stat %{customdata[1]}, q %{customdata[2]})<extra></extra>"
        ),
        meta={"kind": "frequency"},
    ), row=1, col=2)

    # Separate the gene blocks of consecutive groups.
    boundaries = np.cumsum([len(b) for b in blocks])[:-1]
    for b in boundaries:
        fig.add_hline(y=b - 0.5, line=dict(color="#888888", width=0.8), row=1, col=2)

    fig.update_yaxes(autorange="reversed")
    fig.update_yaxes(showticklabels=show_gene_labels, row=1, col=1)
    fig.update_xaxes(showticklabels=False, ticks="", row=1, col=1)
    fig.update_xaxes(title_text=column, row=1, col=2)
    fig.update_layout(
        title=f"Top indicator genes per {column}",
        height=max(400, 16 * len(genes) + 160),
    )
    return fig


def specificity_fidelity(
    pm: PavMatrix,
    indicators: pd.DataFrame,
    column: str,
    q_max: float = 0.05,
    levels: tuple[float, ...] = (0.5, 0.7, 0.9),
    cols: int = 4,
) -> go.Figure:
    """Plot specificity (A) against fidelity (B) of the significant indicator genes.

    One panel per group shows the genes with `q_value < q_max` assigned to that
    group. Genes in the upper right are carried by most of the group's samples and
    are rare elsewhere; high A with low B marks rare but group-restricted genes, and
    high B with lower A marks genes common in the group but shared with others.
    Grey curves mark equal IndVal.g statistic, sqrt(A * B) = level, and are
    labelled in the first panel.

    Args:
        pm: A PavMatrix instance.
        indicators: The DataFrame returned by
            `pavvis.stats.indicator.indicator_genes(pm, column, ...)`.
        column: The discrete metadata column the indicators were computed for.
        q_max: Genes are shown only if their q-value is below this. Defaults to 0.05.
        levels: IndVal.g values drawn as reference curves. Pass `()` for none.
            Defaults to (0.5, 0.7, 0.9).
        cols: Number of panels per row. Defaults to 4.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        ValueError: If column is not a discrete metadata column, if `indicators`
            contains groups not in `column` or genes not in the PAV matrix, or if no
            gene has `q_value < q_max`.

    Warns:
        UserWarning: If some groups have no gene with `q_value < q_max`.
    """
    groups, _, sig = _significant_indicators(pm, indicators, column, q_max)

    cols = min(cols, len(groups))
    rows = -(-len(groups) // cols)
    fig = make_subplots(rows=rows, cols=cols, shared_xaxes=True, shared_yaxes=True,
                        subplot_titles=groups, horizontal_spacing=0.03,
                        vertical_spacing=0.12 if rows > 1 else 0)

    a_grid = np.linspace(0.01, 1, 200)
    empty = []
    for i, g in enumerate(groups):
        row, col = i // cols + 1, i % cols + 1
        for s in levels:
            b_curve = s ** 2 / a_grid
            keep = b_curve <= 1
            fig.add_trace(go.Scatter(
                x=a_grid[keep], y=b_curve[keep], mode="lines", showlegend=False,
                hoverinfo="skip", line=dict(color="#bbbbbb", width=0.8),
                meta={"kind": "level", "level": s},
            ), row=row, col=col)

        genes = sig[sig["group"] == g]
        if genes.empty:
            empty.append(g)
            continue
        fig.add_trace(go.Scatter(
            x=genes["A"], y=genes["B"], text=genes.index, mode="markers", name=g,
            showlegend=False,
            marker=dict(color=_GROUP_COLOURS[i % len(_GROUP_COLOURS)], size=4, opacity=0.5),
            customdata=genes[["stat", "q_value"]].values,
            hovertemplate=("%{text}<br>A: %{x:.3f}<br>B: %{y:.3f}<br>"
                           "stat: %{customdata[0]:.3f}<br>q: %{customdata[1]:.3g}<extra>" + g + "</extra>"),
            meta={"kind": "genes", "group": g},
        ), row=row, col=col)
    _warn_empty_groups(empty, q_max, "empty panels")

    for s in levels:
        fig.add_annotation(x=1, y=s ** 2, xref="x", yref="y", text=f"{s}", showarrow=False,
                           xanchor="left", yanchor="middle", xshift=2,
                           font=dict(color="#888888", size=10))

    fig.update_xaxes(range=[0, 1.02], dtick=0.5)
    fig.update_yaxes(range=[0, 1.02], dtick=0.5)
    fig.update_xaxes(title_text="Specificity (A)", row=rows)
    fig.update_yaxes(title_text="Fidelity (B)", col=1)
    fig.update_layout(
        title="Specificity vs fidelity of indicator genes",
        height=280 * rows + 120,
    )
    return fig
