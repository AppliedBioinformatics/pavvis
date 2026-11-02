from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from scipy.spatial import ConvexHull, QhullError

from pavvis._validation import infer_column_type

pio.templates.default = "simple_white"

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix

# Distinct colours for metadata groups (colourblind-friendly palette).
_GROUP_COLOURS = [
    "#2a9d8f",  # teal
    "#e9c46a",  # yellow
    "#e76f51",  # orange-red
    "#6a4c93",  # indigo
    "#264653",  # dark teal
    "#f4a261",  # amber
    "#a8dadc",  # light blue
    "#457b9d",  # steel blue
]


def _hex_to_rgba(hex_colour: str, alpha: float) -> str:
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (1, 3, 5))
    return f"rgba({r},{g},{b},{alpha})"


def _centroid_traces(
    coords: pd.DataFrame,
    labels: pd.Series,
    centres: pd.DataFrame,
    xcol: str,
    ycol: str,
    column: str,
    centre_name: str,
    hull: bool,
) -> list[go.Scatter]:
    """Hull, spider, point and centre traces for each group in an ordination.

    Each group gets up to four traces sharing a `legendgroup`, tagged with
    `meta={"kind": ...}` so callers and tests can find them. Groups follow the order
    of `centres.index` and are coloured from `_GROUP_COLOURS`.
    """
    traces = []
    for i, group in enumerate(centres.index):
        colour = _GROUP_COLOURS[i % len(_GROUP_COLOURS)]
        pts = coords.loc[(labels == group).values, [xcol, ycol]]
        cx, cy = centres.loc[group, xcol], centres.loc[group, ycol]

        if hull and len(pts) >= 3:
            try:
                vertices = ConvexHull(pts.values).vertices
            except QhullError:  # collinear or duplicate points
                vertices = None
            if vertices is not None:
                ring = pts.values[np.append(vertices, vertices[0])]
                traces.append(go.Scatter(
                    x=ring[:, 0], y=ring[:, 1], mode="lines", fill="toself",
                    fillcolor=_hex_to_rgba(colour, 0.12),
                    line=dict(color=_hex_to_rgba(colour, 0.5), width=1),
                    legendgroup=group, showlegend=False, hoverinfo="skip", meta={"kind": "hull"},
                ))

        n = len(pts)
        spider_x = np.column_stack([pts[xcol], np.full(n, cx), np.full(n, np.nan)]).ravel()
        spider_y = np.column_stack([pts[ycol], np.full(n, cy), np.full(n, np.nan)]).ravel()
        traces.append(go.Scatter(
            x=spider_x, y=spider_y, mode="lines",
            line=dict(color=_hex_to_rgba(colour, 0.45), width=1),
            legendgroup=group, showlegend=False, hoverinfo="skip", meta={"kind": "spider"},
        ))
        traces.append(go.Scatter(
            x=pts[xcol], y=pts[ycol], mode="markers", name=f"{group} (n={n})",
            text=pts.index, marker=dict(color=colour, size=7, line=dict(color="white", width=0.5)),
            legendgroup=group, meta={"kind": "points"},
            hovertemplate=f"{column}: {group}<br>%{{text}}<extra></extra>",
        ))
        traces.append(go.Scatter(
            x=[cx], y=[cy], mode="markers",
            marker=dict(color=colour, size=15, symbol="diamond", line=dict(color="black", width=1.5)),
            legendgroup=group, showlegend=False, meta={"kind": "centroid"},
            hovertemplate=f"{group} {centre_name}<extra></extra>",
        ))
    return traces


def _resolve_gene_index(
    pm: PavMatrix,
    gene_set: Literal["all", "variable", "core", "soft_core", "dispensable", "private"],
) -> pd.Index | list[str]:
    return {
        "all": pm.pav.index,
        "variable": pm.variable_genes,
        "core": pm.core_genes,
        "soft_core": pm.soft_core_genes,
        "dispensable": pm.dispensable_genes,
        "private": pm.private_genes,
    }[gene_set]


def _gene_set_labels(
    gene_set: Literal["all", "variable", "core", "soft_core", "dispensable", "private"],
) -> tuple[str, str]:
    return {
        "all": ("Gene Count per Sample", "Number of Genes"),
        "variable": ("Variable Gene Count per Sample", "Number of Variable Genes"),
        "core": ("Core Gene Count per Sample", "Number of Core Genes"),
        "soft_core": ("Soft Core Gene Count per Sample", "Number of Soft Core Genes"),
        "dispensable": ("Dispensable Gene Count per Sample", "Number of Dispensable Genes"),
        "private": ("Private Gene Count per Sample", "Number of Private Genes"),
    }[gene_set]


def _embedding_scatter_traces(
    emb: pd.DataFrame,
    xcol: str,
    ycol: str,
    metadata: pd.DataFrame,
    color_by: str | None,
) -> list[go.Scatter]:
    """Build 2D scatter traces for a sample embedding, coloured by an optional metadata column.

    Discrete columns give one trace per category; continuous columns give a single
    trace on a Viridis colorscale. The caller is responsible for validating color_by.
    """
    samples = emb.index.tolist()

    if color_by is None:
        return [go.Scatter(
            x=emb[xcol].tolist(), y=emb[ycol].tolist(), text=samples, mode="markers",
            name="Samples",
            marker=dict(color="#2a9d8f", size=7),
            hovertemplate="%{text}<extra></extra>",
        )]

    if infer_column_type(metadata[color_by]) == "discrete":
        col = metadata.loc[emb.index, color_by].astype(str)
        traces = []
        for cat in sorted(col.dropna().unique()):
            mask = (col == cat).values
            traces.append(go.Scatter(
                x=emb.loc[mask, xcol].tolist(),
                y=emb.loc[mask, ycol].tolist(),
                text=emb.index[mask].tolist(),
                mode="markers", name=cat,
                marker=dict(size=7),
                hovertemplate=f"{color_by}: {cat}<br>%{{text}}<extra></extra>",
            ))
        return traces

    color_vals = metadata.loc[samples, color_by].tolist()
    return [go.Scatter(
        x=emb[xcol].tolist(), y=emb[ycol].tolist(), text=samples, mode="markers",
        name=color_by,
        marker=dict(
            color=color_vals, colorscale="Viridis", size=7,
            colorbar=dict(title=color_by), showscale=True,
        ),
        hovertemplate=f"{color_by}: %{{marker.color:.3g}}<br>%{{text}}<extra></extra>",
    )]
