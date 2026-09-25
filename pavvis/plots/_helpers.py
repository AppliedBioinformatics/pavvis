from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import pandas as pd
import plotly.io as pio

pio.templates.default = "simple_white"

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix


def _resolve_gene_index(
    pm: PavMatrix,
    gene_set: Literal["all", "variable", "core"],
) -> pd.Index | list[str]:
    return {
        "all": pm.pav.index,
        "variable": pm.variable_genes,
        "core": pm.core_genes,
    }[gene_set]


def _gene_set_labels(
    gene_set: Literal["all", "variable", "core"],
) -> tuple[str, str]:
    return {
        "all": ("Gene Count per Sample", "Number of Genes"),
        "variable": ("Variable Gene Count per Sample", "Number of Variable Genes"),
        "core": ("Core Gene Count per Sample", "Number of Core Genes"),
    }[gene_set]
