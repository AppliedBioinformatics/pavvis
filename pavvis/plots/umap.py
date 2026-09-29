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


def scatter3d(
    pm: PavMatrix,
    color_by: str | None = None,
) -> go.Figure:
    """Plot the UMAP embedding of samples as a 3D scatter plot.

    Uses all three components of the cached embedding computed by
    pm.compute_umap(). Coloring behaves identically to scatter().

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
    z = emb["UMAP3"].tolist()
    samples = emb.index.tolist()

    if color_by is None:
        traces = [go.Scatter3d(
            x=x, y=y, z=z, text=samples, mode="markers",
            name="Samples",
            marker=dict(color="#2a9d8f", size=4),
            hovertemplate="%{text}<extra></extra>",
        )]

    elif infer_column_type(pm.metadata[color_by]) == "discrete":
        col = pm.metadata[color_by].astype(str)
        categories = sorted(col.dropna().unique())
        traces = []
        for cat in categories:
            mask = col == cat
            cat_samples = emb.index[mask].tolist()
            traces.append(go.Scatter3d(
                x=emb.loc[mask, "UMAP1"].tolist(),
                y=emb.loc[mask, "UMAP2"].tolist(),
                z=emb.loc[mask, "UMAP3"].tolist(),
                text=cat_samples,
                mode="markers", name=cat,
                marker=dict(size=4),
                hovertemplate=f"{color_by}: {cat}<br>%{{text}}<extra></extra>",
            ))

    else:  # continuous
        color_vals = pm.metadata.loc[samples, color_by].tolist()
        traces = [go.Scatter3d(
            x=x, y=y, z=z, text=samples, mode="markers",
            name=color_by,
            marker=dict(
                color=color_vals, colorscale="Viridis", size=4,
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