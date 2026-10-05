from __future__ import annotations

from typing import TYPE_CHECKING

import plotly.graph_objects as go
import plotly.io as pio

from pavvis._validation import validate_color_by, infer_column_type

pio.templates.default = "simple_white"

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix


def _require_embedding(pm: PavMatrix) -> None:
    if pm._umap_embedding is None:
        raise RuntimeError(
            "No UMAP embedding found. Call pm.compute_umap() before plotting."
        )


def scatter(
    pm: PavMatrix,
    color_by: str | None = None,
) -> go.Figure:
    """Plot the UMAP embedding of samples as a 2D scatter plot.

    Uses the first two components of the cached 3D embedding computed by
    pm.compute_umap(). Points represent samples. Color is optional and can
    be driven by any metadata column — discrete columns produce a categorical
    legend; continuous columns use a Viridis colorscale.

    Args:
        pm: A PavMatrix instance with a cached embedding (call pm.compute_umap() first).
        color_by: An optional metadata column name to colour points by.
            Accepts both discrete and continuous columns. Defaults to None.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        RuntimeError: If pm.compute_umap() has not been called.
        ValueError: If color_by is not a valid metadata column.
    """
    _require_embedding(pm)
    if color_by is not None:
        validate_color_by(color_by, pm.metadata)

    emb = pm._umap_embedding
    x = emb["UMAP1"].tolist()
    y = emb["UMAP2"].tolist()
    samples = emb.index.tolist()

    if color_by is None:
        traces = [go.Scatter(
            x=x, y=y, text=samples, mode="markers",
            name="Samples",
            marker=dict(color="#2a9d8f", size=7),
            hovertemplate="%{text}<extra></extra>",
        )]

    elif infer_column_type(pm.metadata[color_by]) == "discrete":
        col = pm.metadata[color_by].astype(str)
        categories = sorted(col.dropna().unique())
        traces = []
        for cat in categories:
            mask = col == cat
            cat_samples = emb.index[mask].tolist()
            traces.append(go.Scatter(
                x=emb.loc[mask, "UMAP1"].tolist(),
                y=emb.loc[mask, "UMAP2"].tolist(),
                text=cat_samples,
                mode="markers", name=cat,
                marker=dict(size=7),
                hovertemplate=f"{color_by}: {cat}<br>%{{text}}<extra></extra>",
            ))

    else:  # continuous
        color_vals = pm.metadata.loc[samples, color_by].tolist()
        traces = [go.Scatter(
            x=x, y=y, text=samples, mode="markers",
            name=color_by,
            marker=dict(
                color=color_vals, colorscale="Viridis", size=7,
                colorbar=dict(title=color_by), showscale=True,
            ),
            hovertemplate=f"{color_by}: %{{marker.color:.3g}}<br>%{{text}}<extra></extra>",
        )]

    fig = go.Figure(traces)
    fig.update_layout(
        title="UMAP Embedding",
        xaxis_title="UMAP1",
        yaxis_title="UMAP2",
    )
    return fig


def _depth_sizes(
    z: list[float],
    size_min: float = 2.0,
    size_max: float = 8.0,
) -> list[float]:
    """Map UMAP3 values linearly onto [size_min, size_max]."""
    import numpy as np
    arr = np.array(z, dtype=float)
    lo, hi = arr.min(), arr.max()
    if hi == lo:
        return [((size_min + size_max) / 2)] * len(z)
    return ((arr - lo) / (hi - lo) * (size_max - size_min) + size_min).tolist()


def scatter3d(
    pm: PavMatrix,
    color_by: str | None = None,
    depth: bool = False,
) -> go.Figure:
    """Plot the UMAP embedding of samples as a 3D scatter plot.

    Uses all three components of the cached embedding computed by
    pm.compute_umap(). Coloring behaves identically to scatter().

    Args:
        pm: A PavMatrix instance with a cached embedding (call pm.compute_umap() first).
        color_by: An optional metadata column name to colour points by.
            Accepts both discrete and continuous columns. Defaults to None.
        depth: If True, scale each point's size by its UMAP3 (z) coordinate so
            that points with a higher z value appear larger, giving a visual
            depth cue. Defaults to False.

    Returns:
        A Plotly Figure. Call .show() to display or .write_html() to export.

    Raises:
        RuntimeError: If pm.compute_umap() has not been called.
        ValueError: If color_by is not a valid metadata column.
    """
    _require_embedding(pm)
    if color_by is not None:
        validate_color_by(color_by, pm.metadata)

    emb = pm._umap_embedding
    x = emb["UMAP1"].tolist()
    y = emb["UMAP2"].tolist()
    z = emb["UMAP3"].tolist()
    samples = emb.index.tolist()

    def _sizes(mask=None) -> list[float] | float:
        if not depth:
            return 4
        z_sub = [z[i] for i, m in enumerate(mask) if m] if mask is not None else z
        return _depth_sizes(z_sub)

    if color_by is None:
        traces = [go.Scatter3d(
            x=x, y=y, z=z, text=samples, mode="markers",
            name="Samples",
            marker=dict(color="#2a9d8f", size=_sizes()),
            hovertemplate="%{text}<extra></extra>",
        )]

    elif infer_column_type(pm.metadata[color_by]) == "discrete":
        col = pm.metadata[color_by].astype(str)
        categories = sorted(col.dropna().unique())
        traces = []
        for cat in categories:
            bool_mask = col == cat
            mask = bool_mask.tolist()
            cat_emb = emb.loc[bool_mask]
            traces.append(go.Scatter3d(
                x=cat_emb["UMAP1"].tolist(),
                y=cat_emb["UMAP2"].tolist(),
                z=cat_emb["UMAP3"].tolist(),
                text=cat_emb.index.tolist(),
                mode="markers", name=cat,
                marker=dict(size=_sizes(mask)),
                hovertemplate=f"{color_by}: {cat}<br>%{{text}}<extra></extra>",
            ))

    else:  # continuous
        color_vals = pm.metadata.loc[samples, color_by].tolist()
        traces = [go.Scatter3d(
            x=x, y=y, z=z, text=samples, mode="markers",
            name=color_by,
            marker=dict(
                color=color_vals, colorscale="Viridis", size=_sizes(),
                colorbar=dict(title=color_by), showscale=True,
            ),
            hovertemplate=f"{color_by}: %{{marker.color:.3g}}<br>%{{text}}<extra></extra>",
        )]

    fig = go.Figure(traces)
    fig.update_layout(
        title="UMAP Embedding (3D)",
        scene=dict(
            xaxis_title="UMAP1",
            yaxis_title="UMAP2",
            zaxis_title="UMAP3",
        ),
    )
    return fig