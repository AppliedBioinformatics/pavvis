from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix


def _indval_components(pav: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Group-equalised specificity (A) and fidelity (B) for every gene and group.

    For binary data the mean of a gene within a group equals its fidelity, so
    A = B / sum(B over groups).

    Args:
        pav: Binary array of shape (genes, samples).
        X: One-hot array of shape (samples, groups).

    Returns:
        Tuple of (A, B), each of shape (genes, groups).
    """
    B = (pav @ X) / X.sum(axis=0)
    total = B.sum(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        A = np.where(total > 0, B / total, 0.0)
    return A, B


def indicator_genes(
    pm: PavMatrix,
    column: str,
    genes: list[str] | None = None,
    permutations: int = 999,
    random_state: int | None = None,
) -> pd.DataFrame:
    """Indicator gene analysis (IndVal.g) for the groups of a discrete metadata column.

    For each gene and group, specificity A is the gene's mean presence in the group
    divided by the sum of its mean presence across all groups (group-size equalised),
    and fidelity B is the fraction of the group's samples carrying the gene. The
    indicator value is sqrt(A·B), on the same [0, 1] scale as
    `indicspecies::multipatt(func="IndVal.g", duleg=TRUE)`. Each gene is assigned to
    the group with the highest value.

    Significance comes from permuting sample labels: p is the fraction of
    permutations where a gene's best-group value is at least as large as observed,
    with the observed value counted once. q-values use the Benjamini–Hochberg
    procedure over the tested genes only.

    Args:
        pm: A PavMatrix instance.
        column: A discrete metadata column defining the groups.
        genes: Gene IDs to test. Defaults to None (all variable genes). Core and
            absent genes cannot be indicators, so including them only weakens the
            multiple-testing correction.
        permutations: Number of label permutations. Defaults to 999.
        random_state: Random seed for reproducibility. Defaults to None.

    Returns:
        DataFrame indexed by gene ID with columns `group`, `A`, `B`, `stat`,
        `p_value` and `q_value`, sorted by `stat` in descending order.

    Raises:
        ValueError: If column is not a discrete metadata column, or if any gene ID
            is not found.
    """
    groups, X = pm._group_matrix(column)
    if genes is None:
        genes = pm.variable_genes
    else:
        unknown = [g for g in genes if g not in pm.pav.index]
        if unknown:
            raise ValueError(f"Genes not found in PAV matrix: {unknown}")

    pav = pm.pav.loc[genes].values.astype(np.float64)
    X = X.astype(np.float64)
    rows = np.arange(len(genes))

    A, B = _indval_components(pav, X)
    stat_all = np.sqrt(A * B)
    best = stat_all.argmax(axis=1)
    stat = stat_all[rows, best]

    rng = np.random.default_rng(random_state)
    n_ge = np.ones(len(genes))  # the observed value counts as one permutation
    for _ in range(permutations):
        A_p, B_p = _indval_components(pav, X[rng.permutation(len(X))])
        n_ge += np.sqrt(A_p * B_p).max(axis=1) >= stat - 1e-12
    p_values = n_ge / (permutations + 1)

    result = pd.DataFrame({
        "group": np.array(groups, dtype=object)[best],
        "A": A[rows, best],
        "B": B[rows, best],
        "stat": stat,
        "p_value": p_values,
        "q_value": false_discovery_control(p_values, method="bh") if len(genes) else p_values,
    }, index=pd.Index(genes, name="gene_id"))
    return result.sort_values("stat", ascending=False, kind="stable")
