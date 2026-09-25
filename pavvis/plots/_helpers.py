from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import pandas as pd
import plotly.io as pio

pio.templates.default = "simple_white"

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix


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
