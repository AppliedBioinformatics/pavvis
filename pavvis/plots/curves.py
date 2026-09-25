from __future__ import annotations

from typing import TYPE_CHECKING, Sequence

import numpy as np
import plotly.graph_objects as go

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix


def legacy_pangenome_curve(
    pm: PavMatrix,
    permutations: int = 100,
    log_y: bool = False,
) -> go.Figure:
    """Plot pangenome and core genome curves across increasing genome counts.

    For each permutation a random sample order is drawn. At each step k the
    pangenome size (union of all genes seen so far) and core size (intersection
    of all genes seen so far) are recorded. The mean and ±1 std band across
    permutations are plotted for both curves.

    Args:
        pm: A PavMatrix instance.
        permutations: Number of random sample orderings to average over.
            Defaults to 100.
        log_y: If True, display the y-axis on a log10 scale. Defaults to False.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.
    """
    pav = pm.pav.values  # shape: (genes, samples)
    n_samples = pav.shape[1]
    x = list(range(1, n_samples + 1))

    pan_runs = np.empty((permutations, n_samples), dtype=np.int32)
    core_runs = np.empty((permutations, n_samples), dtype=np.int32)

    rng = np.random.default_rng()
    for i in range(permutations):
        order = rng.permutation(n_samples)
        shuffled = pav[:, order]
        pan_runs[i] = (np.cumsum(shuffled, axis=1) > 0).sum(axis=0)
        core_runs[i] = np.cumprod(shuffled, axis=1).sum(axis=0)

    pan_mean = pan_runs.mean(axis=0)
    pan_upper = pan_mean + pan_runs.std(axis=0)
    pan_lower = pan_mean - pan_runs.std(axis=0)

    core_mean = core_runs.mean(axis=0)
    core_upper = core_mean + core_runs.std(axis=0)
    core_lower = core_mean - core_runs.std(axis=0)

    traces = [
        # Pangenome band
        go.Scatter(
            x=x, y=pan_upper.tolist(),
            mode="lines", line=dict(width=0),
            legendgroup="pangenome", showlegend=False, hoverinfo="skip",
        ),
        go.Scatter(
            x=x, y=pan_lower.tolist(),
            mode="lines", line=dict(width=0),
            fill="tonexty", fillcolor="rgba(99,110,250,0.2)",
            legendgroup="pangenome", showlegend=False, hoverinfo="skip",
        ),
        # Core band
        go.Scatter(
            x=x, y=core_upper.tolist(),
            mode="lines", line=dict(width=0),
            legendgroup="core", showlegend=False, hoverinfo="skip",
        ),
        go.Scatter(
            x=x, y=core_lower.tolist(),
            mode="lines", line=dict(width=0),
            fill="tonexty", fillcolor="rgba(239,85,59,0.2)",
            legendgroup="core", showlegend=False, hoverinfo="skip",
        ),
        # Mean lines (these own the legend entries)
        go.Scatter(
            x=x, y=pan_mean.tolist(),
            mode="lines", name="Pangenome",
            line=dict(color="rgb(99,110,250)", width=2),
            legendgroup="pangenome",
        ),
        go.Scatter(
            x=x, y=core_mean.tolist(),
            mode="lines", name="Core genome",
            line=dict(color="rgb(239,85,59)", width=2),
            legendgroup="core",
        ),
    ]

    fig = go.Figure(traces)
    fig.update_layout(
        title=f"Pangenome Curve ({permutations} permutations)",
        xaxis_title="Number of Genomes",
        yaxis_title="Number of Genes",
        yaxis=dict(type="log" if log_y else "linear"),
    )

    return fig


def pangenome_curve(
    pm: PavMatrix,
    permutations: int = 100,
    log_y: bool = False,
) -> go.Figure:
    """Plot core and variable genome curves on separate left/right y-axes.

    For each permutation a random sample order is drawn. At each step k:
      - Core genes: genes present in all k genomes so far (left y-axis, decays).
      - Variable genes: genes present in some but not all k genomes so far
        (right y-axis, rises then plateaus).

    Mean and ±1 std bands across permutations are shown for both curves.
    Clicking a legend entry toggles its mean line and std band together.

    Args:
        pm: A PavMatrix instance.
        permutations: Number of random sample orderings to average over.
            Defaults to 100.
        log_y: If True, display both y-axes on a log10 scale. Defaults to False.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.
    """
    pav = pm.pav.values  # shape: (genes, samples)
    n_genes, n_samples = pav.shape
    x = list(range(1, n_samples + 1))
    axis_type = "log" if log_y else "linear"

    core_runs = np.empty((permutations, n_samples), dtype=np.int32)
    var_runs = np.empty((permutations, n_samples), dtype=np.int32)

    rng = np.random.default_rng()
    for i in range(permutations):
        order = rng.permutation(n_samples)
        shuffled = pav[:, order]
        # Core: gene must be present in every genome seen so far (cumprod stays 1)
        core_counts = np.cumprod(shuffled, axis=1).sum(axis=0)
        # Pangenome: total unique genes seen so far
        pan_counts = (np.cumsum(shuffled, axis=1) > 0).sum(axis=0)
        # Variable: in pangenome but not in core
        core_runs[i] = core_counts
        var_runs[i] = pan_counts - core_counts

    def _band_traces(runs, color_rgb, group, yaxis):
        mean = runs.mean(axis=0)
        std = runs.std(axis=0)
        upper = (mean + std).tolist()
        lower = (mean - std).tolist()
        return [
            go.Scatter(
                x=x, y=upper,
                mode="lines", line=dict(width=0),
                legendgroup=group, showlegend=False, hoverinfo="skip",
                yaxis=yaxis,
            ),
            go.Scatter(
                x=x, y=lower,
                mode="lines", line=dict(width=0),
                fill="tonexty", fillcolor=f"rgba({color_rgb},0.2)",
                legendgroup=group, showlegend=False, hoverinfo="skip",
                yaxis=yaxis,
            ),
        ], mean.tolist()

    core_band_traces, core_mean = _band_traces(core_runs, "239,85,59", "core", "y")
    var_band_traces, var_mean = _band_traces(var_runs, "99,110,250", "variable", "y2")

    traces = [
        *core_band_traces,
        *var_band_traces,
        go.Scatter(
            x=x, y=core_mean,
            mode="lines", name="Core genes",
            line=dict(color="rgb(239,85,59)", width=2),
            legendgroup="core", yaxis="y",
        ),
        go.Scatter(
            x=x, y=var_mean,
            mode="lines", name="Variable genes",
            line=dict(color="rgb(99,110,250)", width=2),
            legendgroup="variable", yaxis="y2",
        ),
    ]

    fig = go.Figure(traces)
    fig.update_layout(
        title=f"Core & Variable Gene Curves ({permutations} permutations)",
        xaxis_title="Number of Genomes",
        yaxis=dict(title="Core Genes", type=axis_type, side="left"),
        yaxis2=dict(title="Variable Genes", type=axis_type, side="right", overlaying="y"),
    )

    return fig


def grouped_variable_curve(
    pm: PavMatrix,
    column: str,
    group_order: Sequence[str],
    permutations: int = 100,
    log_y: bool = False,
    shade_groups: bool = False,
) -> go.Figure:
    """Plot variable gene accumulation grouped by a metadata column.

    Groups are added to the curve in the order specified by `group_order`.
    Within each group the sample order is randomly permuted `permutations`
    times; the group sequence itself is fixed. The result shows how the
    variable gene pool grows as each group is introduced, with a ±1 std
    band reflecting within-group ordering uncertainty. Vertical lines mark
    the boundary where each new group begins.

    Args:
        pm: A PavMatrix instance.
        column: Metadata column used to assign samples to groups.
        group_order: Sequence of group labels defining the order in which
            groups are added to the curve. Every value must be present in
            pm.metadata[column].
        permutations: Number of within-group random orderings to average
            over. Defaults to 100.
        log_y: If True, display the y-axis on a log10 scale. Defaults to False.
        shade_groups: If True, alternating groups are shaded with a light grey
            background to make group boundaries easier to read. Defaults to False.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        ValueError: If `column` is not in pm.metadata, or if any value in
            `group_order` is not found in pm.metadata[column].
    """
    if column not in pm.metadata.columns:
        raise ValueError(
            f"'{column}' is not a metadata column. "
            f"Available columns: {list(pm.metadata.columns)}"
        )

    col = pm.metadata[column].astype(str)
    missing = [g for g in group_order if g not in col.values]
    if missing:
        raise ValueError(
            f"Groups not found in metadata column '{column}': {missing}. "
            f"Available values: {sorted(col.unique())}"
        )

    # Map each group label to the integer column indices in pm.pav
    pav = pm.pav.values  # (genes, samples)
    sample_ids = list(pm.pav.columns)
    group_indices: list[list[int]] = [
        [sample_ids.index(s) for s in pm.metadata.index[col == g]]
        for g in group_order
    ]

    n_total = pav.shape[1]
    var_runs = np.empty((permutations, n_total), dtype=np.int32)

    rng = np.random.default_rng()
    for i in range(permutations):
        # Build a full column order: shuffle within each group, keep group sequence
        order = np.concatenate([
            rng.permutation(idx) for idx in group_indices
        ])
        shuffled = pav[:, order]
        pan_counts = (np.cumsum(shuffled, axis=1) > 0).sum(axis=0)
        core_counts = np.cumprod(shuffled, axis=1).sum(axis=0)
        var_runs[i] = pan_counts - core_counts

    var_mean = var_runs.mean(axis=0)
    var_std = var_runs.std(axis=0)
    var_upper = (var_mean + var_std).tolist()
    var_lower = (var_mean - var_std).tolist()
    x = list(range(1, n_total + 1))

    # Boundary x positions: where each group starts (1-indexed)
    boundaries: list[int] = []
    cumulative = 0
    for idx in group_indices[:-1]:  # no line after the last group
        cumulative += len(idx)
        boundaries.append(cumulative + 1)

    traces = [
        go.Scatter(
            x=x, y=var_upper,
            mode="lines", line=dict(width=0),
            legendgroup="variable", showlegend=False, hoverinfo="skip",
        ),
        go.Scatter(
            x=x, y=var_lower,
            mode="lines", line=dict(width=0),
            fill="tonexty", fillcolor="rgba(99,110,250,0.2)",
            legendgroup="variable", showlegend=False, hoverinfo="skip",
        ),
        go.Scatter(
            x=x, y=var_mean.tolist(),
            mode="lines", name="Variable genes",
            line=dict(color="rgb(99,110,250)", width=2),
            legendgroup="variable",
        ),
    ]

    fig = go.Figure(traces)

    # Group boundary lines, shading, and labels
    cumulative = 0
    for i, (group, idx) in enumerate(zip(group_order, group_indices)):
        x0 = cumulative + 0.5
        x1 = cumulative + len(idx) + 0.5
        mid_x = (x0 + x1) / 2

        if shade_groups and i % 2 == 1:
            fig.add_vrect(
                x0=x0, x1=x1,
                fillcolor="rgba(0,0,0,0.06)", layer="below",
                line_width=0,
            )

        fig.add_annotation(
            x=mid_x, y=1, yref="paper",
            text=group, showarrow=False,
            font=dict(size=11), yanchor="bottom",
        )
        cumulative += len(idx)

    for bx in boundaries:
        fig.add_vline(x=bx - 0.5, line=dict(color="grey", width=1, dash="dash"))

    fig.update_layout(
        title=f"Variable Gene Accumulation by {column} ({permutations} permutations)",
        xaxis_title="Number of Genomes",
        yaxis_title="Number of Variable Genes",
        yaxis=dict(type="log" if log_y else "linear"),
    )

    return fig


def jaccard_similarity_curve(
    pm: PavMatrix,
    column: str,
    group_order: Sequence[str],
    permutations: int = 100,
    shade_groups: bool = False,
) -> go.Figure:
    """Plot mean pairwise Jaccard similarity as genomes are accumulated by group.

    Groups are added in the order specified by `group_order`. Within each group
    the sample order is randomly permuted `permutations` times; the group
    sequence itself is fixed. At each step k the mean Jaccard similarity across
    all C(k,2) pairs of genomes accumulated so far is recorded. A drop when a
    new group is introduced indicates that group is divergent from those already
    seen; stability indicates similarity. The mean and ±1 std band across
    permutations are shown. Clicking the legend entry toggles both.

    Args:
        pm: A PavMatrix instance.
        column: Metadata column used to assign samples to groups.
        group_order: Sequence of group labels defining the order in which
            groups are added. Every value must be present in pm.metadata[column].
        permutations: Number of within-group random orderings to average over.
            Defaults to 100.
        shade_groups: If True, alternating groups are shaded with a light grey
            background. Defaults to False.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        ValueError: If `column` is not in pm.metadata, or if any value in
            `group_order` is not found in pm.metadata[column].
    """
    if column not in pm.metadata.columns:
        raise ValueError(
            f"'{column}' is not a metadata column. "
            f"Available columns: {list(pm.metadata.columns)}"
        )

    col = pm.metadata[column].astype(str)
    missing = [g for g in group_order if g not in col.values]
    if missing:
        raise ValueError(
            f"Groups not found in metadata column '{column}': {missing}. "
            f"Available values: {sorted(col.unique())}"
        )

    pav = pm.pav.values.astype(np.float32)  # (genes, samples)
    sample_ids = list(pm.pav.columns)
    group_indices: list[list[int]] = [
        [sample_ids.index(s) for s in pm.metadata.index[col == g]]
        for g in group_order
    ]

    n_total = pav.shape[1]
    # x starts at 2 — need at least one pair
    x = list(range(2, n_total + 1))
    sim_runs = np.empty((permutations, n_total - 1), dtype=np.float32)

    rng = np.random.default_rng()
    for i in range(permutations):
        order = np.concatenate([
            rng.permutation(idx) for idx in group_indices
        ])
        for k in range(2, n_total + 1):
            sub = pav[:, order[:k]]
            inter = sub.T @ sub
            counts = sub.sum(axis=0)
            union = counts[:, None] + counts[None, :] - inter
            idx_upper = np.triu_indices(k, k=1)
            sim_runs[i, k - 2] = (inter[idx_upper] / union[idx_upper]).mean()

    sim_mean = sim_runs.mean(axis=0)
    sim_std = sim_runs.std(axis=0)

    traces = [
        go.Scatter(
            x=x, y=(sim_mean + sim_std).tolist(),
            mode="lines", line=dict(width=0),
            legendgroup="jaccard", showlegend=False, hoverinfo="skip",
        ),
        go.Scatter(
            x=x, y=(sim_mean - sim_std).tolist(),
            mode="lines", line=dict(width=0),
            fill="tonexty", fillcolor="rgba(0,204,150,0.2)",
            legendgroup="jaccard", showlegend=False, hoverinfo="skip",
        ),
        go.Scatter(
            x=x, y=sim_mean.tolist(),
            mode="lines", name="Mean Jaccard similarity",
            line=dict(color="rgb(0,204,150)", width=2),
            legendgroup="jaccard",
        ),
    ]

    fig = go.Figure(traces)

    # Boundary lines, shading, and group labels
    boundaries: list[float] = []
    cumulative = 0
    for idx in group_indices[:-1]:
        cumulative += len(idx)
        boundaries.append(cumulative + 0.5)

    cumulative = 0
    for i, (group, idx) in enumerate(zip(group_order, group_indices)):
        x0 = cumulative + 0.5
        x1 = cumulative + len(idx) + 0.5
        if shade_groups and i % 2 == 1:
            fig.add_vrect(
                x0=x0, x1=x1,
                fillcolor="rgba(0,0,0,0.06)", layer="below",
                line_width=0,
            )
        fig.add_annotation(
            x=(x0 + x1) / 2, y=1, yref="paper",
            text=group, showarrow=False,
            font=dict(size=11), yanchor="bottom",
        )
        cumulative += len(idx)

    for bx in boundaries:
        fig.add_vline(x=bx, line=dict(color="grey", width=1, dash="dash"))

    fig.update_layout(
        title=f"Jaccard Similarity Curve by {column} ({permutations} permutations)",
        xaxis_title="Number of Genomes",
        yaxis_title="Mean Pairwise Jaccard Similarity",
        yaxis=dict(range=[
            max(0.0, float((sim_mean - sim_std).min()) - 0.05),
            min(1.0, float((sim_mean + sim_std).max()) + 0.05),
        ]),
    )

    return fig
