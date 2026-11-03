from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import plotly.graph_objects as go
import plotly.io as pio

from pavvis._validation import validate_color_by
from pavvis.plots._helpers import _GROUP_COLOURS, _centroid_traces, _embedding_scatter_traces

pio.templates.default = "simple_white"

if TYPE_CHECKING:
    import pandas as pd

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


def _betadisper_for_plot(
    pm: PavMatrix, column: str, betadisper_result: dict | None, type: str | None
) -> dict:
    """Return a supplied betadisper result after checking it matches, or compute one without permutations."""
    if betadisper_result is None:
        from pavvis.stats.permanova import betadisper
        return betadisper(pm, column, type=type or "median", permutations=0)

    if betadisper_result.get("column", column) != column:
        raise ValueError(
            f"betadisper_result was computed for '{betadisper_result['column']}', not '{column}'."
        )
    if type is not None and betadisper_result.get("type", type) != type:
        raise ValueError(
            f"betadisper_result uses type='{betadisper_result['type']}', but type='{type}' was given."
        )
    if list(betadisper_result["distances"].index) != list(pm.sample_ids):
        raise ValueError("betadisper_result samples do not match pm. Was it computed on a different PavMatrix?")
    return betadisper_result


def centroids(
    pm: PavMatrix,
    column: str,
    components: tuple[int, int] = (1, 2),
    hull: bool = True,
    betadisper_result: dict | None = None,
    permanova_result: pd.DataFrame | None = None,
    type: Literal["median", "centroid"] | None = None,
) -> go.Figure:
    """Plot samples and their group centres in PCoA space, like `plot(vegan::betadisper(...))`.

    Each sample is joined to its group's centre by a line ("spider"), so both the
    separation between groups (centre positions) and the spread within groups (line
    lengths) are visible. Sample coordinates and centres come from
    `pavvis.stats.permanova.betadisper`, so the plot shows exactly what the dispersion
    test measured. Centres are computed in the full PCoA space and projected onto the
    two plotted axes, so a centre need not sit in the visual middle of its points.

    The plot does not run any permutation test. Pass the results of
    `pavvis.stats.permanova.betadisper` and `pavvis.stats.permanova.permanova` to
    reuse their ordination and annotate the plot with their statistics, so the
    figure shows exactly the numbers reported elsewhere.

    Clicking a group in the legend toggles its points, lines, centre and hull together.
    Groups are coloured in sorted order; colours repeat after eight groups.

    Args:
        pm: A PavMatrix instance. A PCoA does not need to be computed first.
        column: A discrete metadata column defining the groups.
        components: The two 1-based PCoA axes to plot on x and y. Defaults to (1, 2).
        hull: If True, shade the convex hull of each group with at least 3 samples.
            Defaults to True.
        betadisper_result: The dict returned by `betadisper(pm, column, ...)`. Its
            coordinates and centres are drawn, and its PERMDISP F and p are added to
            the annotation. If None, the ordination is computed without a test.
        permanova_result: The table returned by `permanova(pm, column, ...)`. If
            given, its R², F and p are added to the annotation.
        type: Group centre: 'median' (spatial median, vegan's default) or
            'centroid' (mean). Defaults to the type of `betadisper_result`, or
            'median' if none is given.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        ValueError: If column is not a discrete metadata column with at least 2
            groups, if a component exceeds the number of positive PCoA axes, or if a
            supplied result was computed for a different column, type or PavMatrix.
    """
    bd = _betadisper_for_plot(pm, column, betadisper_result, type)
    type = bd.get("type", type or "median")
    if permanova_result is not None and permanova_result.attrs.get("column", column) != column:
        raise ValueError(
            f"permanova_result was computed for '{permanova_result.attrs['column']}', not '{column}'."
        )
    coords, cents, explained = bd["coordinates"], bd["centroids"], bd["explained"]
    n_axes = coords.shape[1]
    if any(not 1 <= c <= n_axes for c in components):
        raise ValueError(
            f"components must be between 1 and {n_axes} (the number of positive PCoA axes), "
            f"got {components}."
        )
    xcol, ycol = (f"PCo{c}" for c in components)
    labels = pm.metadata.loc[coords.index, column].astype(str)

    fig = go.Figure(_centroid_traces(coords, labels, cents, xcol, ycol, column, type, hull))
    fig.update_layout(
        title=f"Group {type}s in PCoA space (Jaccard distance)",
        xaxis_title=f"{xcol} ({explained[xcol] * 100:.1f}%)",
        yaxis_title=f"{ycol} ({explained[ycol] * 100:.1f}%)",
        legend_title_text=column,
    )
    lines, perms = [], {}
    if permanova_result is not None:
        model = permanova_result.loc["Model"]
        lines.append(f"PERMANOVA: R² = {model['R2']:.3f}, F = {model['F']:.2f}, p = {model['Pr(>F)']:.3g}")
        perms["PERMANOVA"] = permanova_result.attrs.get("permutations")
    if betadisper_result is not None:
        lines.append(f"PERMDISP ({type}): F = {bd['F']:.2f}, p = {bd['p_value']:.3g}")
        perms["PERMDISP"] = bd.get("permutations")
    if lines:
        counts = set(perms.values()) - {None}
        if len(counts) == 1:
            lines.append(f"{counts.pop():,} permutations")
        elif counts:
            lines.append("permutations: " + ", ".join(f"{k} {v:,}" for k, v in perms.items() if v is not None))
        fig.add_annotation(
            xref="paper", yref="paper", x=0.01, y=0.99, xanchor="left", yanchor="top",
            align="left", showarrow=False, bgcolor="rgba(255,255,255,0.8)",
            text="<br>".join(lines),
        )
    return fig


def dispersion_boxplot(
    pm: PavMatrix,
    column: str,
    betadisper_result: dict | None = None,
    type: Literal["median", "centroid"] | None = None,
) -> go.Figure:
    """Plot each sample's distance to its group centre as one box per group.

    The equivalent of `boxplot(vegan::betadisper(...))`. Distances come from
    `pavvis.stats.permanova.betadisper`; run that function for the PERMDISP test
    and pass its result here to reuse it.

    Args:
        pm: A PavMatrix instance.
        column: A discrete metadata column defining the groups.
        betadisper_result: The dict returned by `betadisper(pm, column, ...)`. If
            None, the distances are computed without a test.
        type: Group centre: 'median' (spatial median, vegan's default) or
            'centroid' (mean). Defaults to the type of `betadisper_result`, or
            'median' if none is given.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        ValueError: If column is not a discrete metadata column with at least 2
            groups, or if `betadisper_result` was computed for a different column,
            type or PavMatrix.
    """
    bd = _betadisper_for_plot(pm, column, betadisper_result, type)
    type = bd.get("type", type or "median")
    z = bd["distances"]
    labels = pm.metadata.loc[z.index, column].astype(str)

    traces = []
    for i, group in enumerate(bd["centroids"].index):
        mask = (labels == group).values
        traces.append(go.Box(
            y=z[mask], name=group, text=z.index[mask],
            boxpoints="all", jitter=0.4, pointpos=0,
            marker=dict(color=_GROUP_COLOURS[i % len(_GROUP_COLOURS)], size=4),
            hovertemplate="%{text}: %{y:.3f}<extra></extra>",
        ))
    fig = go.Figure(traces)
    fig.update_layout(
        title=f"Distance to group {type}",
        xaxis_title=column,
        yaxis_title=f"Distance to {type}",
        showlegend=False,
    )
    return fig