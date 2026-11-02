from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

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
