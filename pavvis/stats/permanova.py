from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

from pavvis._ordination import pcoa_eig

if TYPE_CHECKING:
    from pavvis.pav_matrix import PavMatrix


def _groups_and_distances(pm: PavMatrix, column: str) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Validated groups, one-hot matrix (samples × groups) and Jaccard distances."""
    groups, X = pm._group_matrix(column)
    if len(groups) < 2:
        raise ValueError(f"'{column}' must have at least 2 groups, found {groups}.")
    pm.compute_jaccard()
    D = 1.0 - pm._jaccard_matrix.loc[pm.sample_ids, pm.sample_ids].values
    return groups, X.astype(np.float64), D


def _within_ss(D2: np.ndarray, X: np.ndarray) -> float:
    """Within-group sum of squares: sum over groups of (sum of squared distances) / (2 n_g)."""
    return float((((D2 @ X) * X).sum(axis=0) / (2 * X.sum(axis=0))).sum())


def permanova(
    pm: PavMatrix,
    column: str,
    permutations: int = 999,
    random_state: int | None = None,
) -> pd.DataFrame:
    """PERMANOVA on Jaccard distance for the groups of a discrete metadata column.

    Equivalent to `vegan::adonis2(dist ~ group, permutations=...)` with Jaccard
    distance (1 − Jaccard similarity, from `pm.compute_jaccard()`). The pseudo-F is
    the ratio of between-group to within-group sums of squared distances, each
    divided by its degrees of freedom. Significance comes from permuting sample
    labels: p is the fraction of permutations with a pseudo-F at least as large as
    observed, with the observed value counted once.

    PERMANOVA is sensitive to differences in dispersion as well as location. Check
    `betadisper` before interpreting a significant result as a shift in centroids.

    Args:
        pm: A PavMatrix instance.
        column: A discrete metadata column defining the groups.
        permutations: Number of label permutations. Defaults to 999.
        random_state: Random seed for reproducibility. Defaults to None.

    Returns:
        An adonis2-style DataFrame with rows `Model`, `Residual` and `Total` and
        columns `Df`, `SumOfSqs`, `R2`, `F` and `Pr(>F)`. `F` and `Pr(>F)` are only
        set on the `Model` row. `attrs` records the `column` and `permutations` used.

    Raises:
        ValueError: If column is not a discrete metadata column or has fewer than
            2 groups.
    """
    groups, X, D = _groups_and_distances(pm, column)
    D2 = D ** 2
    n, k = X.shape
    df_model, df_resid = k - 1, n - k

    ss_total = D2.sum() / (2 * n)
    ss_resid = _within_ss(D2, X)
    ss_model = ss_total - ss_resid
    f_obs = (ss_model / df_model) / (ss_resid / df_resid)

    rng = np.random.default_rng(random_state)
    n_ge = 1  # the observed value counts as one permutation
    for _ in range(permutations):
        ss_w = _within_ss(D2, X[rng.permutation(n)])
        n_ge += ((ss_total - ss_w) / df_model) / (ss_w / df_resid) >= f_obs - 1e-12
    p_value = n_ge / (permutations + 1)

    table = pd.DataFrame(
        {
            "Df": [df_model, df_resid, n - 1],
            "SumOfSqs": [ss_model, ss_resid, ss_total],
            "R2": [ss_model / ss_total, ss_resid / ss_total, 1.0],
            "F": [f_obs, np.nan, np.nan],
            "Pr(>F)": [p_value, np.nan, np.nan],
        },
        index=["Model", "Residual", "Total"],
    )
    table.attrs = {"column": column, "permutations": permutations}
    return table


def _spatial_median(points: np.ndarray, tol: float = 1e-9, max_iter: int = 1000) -> np.ndarray:
    """Geometric median by Weiszfeld iteration (as in vegan's spatialMed)."""
    if points.shape[1] == 0:
        return np.zeros(0)
    y = points.mean(axis=0)
    for _ in range(max_iter):
        d = np.linalg.norm(points - y, axis=1)
        w = 1.0 / np.maximum(d, tol)
        y_new = (points * w[:, None]).sum(axis=0) / w.sum()
        if np.linalg.norm(y_new - y) < tol:
            return y_new
        y = y_new
    return y


def _anova_f(z: np.ndarray, X: np.ndarray) -> float:
    """One-way ANOVA F statistic of z across the groups in one-hot matrix X."""
    n, k = X.shape
    n_g = X.sum(axis=0)
    means = (z @ X) / n_g
    ss_between = (n_g * (means - z.mean()) ** 2).sum()
    ss_within = ((z - X @ means) ** 2).sum()
    return (ss_between / (k - 1)) / (ss_within / (n - k))


def betadisper(
    pm: PavMatrix,
    column: str,
    type: Literal["median", "centroid"] = "median",
    permutations: int = 999,
    random_state: int | None = None,
) -> dict:
    """Multivariate homogeneity of group dispersions (PERMDISP) on Jaccard distance.

    Follows `vegan::betadisper` and `vegan::permutest`. Samples are placed in PCoA
    space using all axes, including those with negative eigenvalues (Jaccard distance
    is not Euclidean). Each group's centre is computed separately in the positive and
    negative axis spaces, and a sample's distance to its group centre is
    z = sqrt(d₊² − d₋²), with negative values of d₊² − d₋² set to zero as in vegan.
    Dispersions are compared with a one-way ANOVA F on z; significance comes from
    refitting the ANOVA to permuted residuals (z minus its group mean), keeping the
    group centres fixed.

    Validated against vegan 2.7.5 on the Amborella PAV data: distances and F match to
    floating-point precision with type='centroid'. With type='median', distances differ
    by ~1e-5 because vegan stops its BFGS search for the spatial median earlier than
    the Weiszfeld iteration used here.

    A significant result means the groups differ in spread. In that case a
    significant PERMANOVA cannot be read purely as a difference in location.

    Args:
        pm: A PavMatrix instance.
        column: A discrete metadata column defining the groups.
        type: Group centre to measure distances from: 'median' (spatial median,
            vegan's default) or 'centroid' (mean). Defaults to 'median'.
        permutations: Number of residual permutations. Defaults to 999.
        random_state: Random seed for reproducibility. Defaults to None.

    Returns:
        Dict with keys:

        - `distances`: pd.Series of distance to group centre, indexed by sample.
        - `centroids`: DataFrame of group centres (groups × `PCo1…`) on the positive
            PCoA axes, matching `pm.compute_pcoa(n_components=None)`.
        - `coordinates`: DataFrame of sample coordinates on the same positive axes.
        - `explained`: pd.Series of the fraction of positive-eigenvalue variation
            each axis explains.
        - `F`: observed F statistic.
        - `p_value`: permutation p-value.
        - `column`, `type`, `permutations`: the arguments used.

    Raises:
        ValueError: If column is not a discrete metadata column or has fewer than
            2 groups, or if type is not 'median' or 'centroid'.
    """
    if type not in ("median", "centroid"):
        raise ValueError(f"type must be 'median' or 'centroid', got '{type}'.")
    groups, X, D = _groups_and_distances(pm, column)

    vals, vecs = pcoa_eig(D)
    coords = vecs * np.sqrt(np.abs(vals))
    pos, neg = coords[:, vals > 0], coords[:, vals < 0]

    centre = _spatial_median if type == "median" else (lambda pts: pts.mean(axis=0))
    membership = X.argmax(axis=1)
    z2 = np.empty(len(D))
    cents = []
    for g in range(len(groups)):
        m = membership == g
        c_pos, c_neg = centre(pos[m]), centre(neg[m])
        z2[m] = ((pos[m] - c_pos) ** 2).sum(axis=1) - ((neg[m] - c_neg) ** 2).sum(axis=1)
        cents.append(c_pos)
    z = np.sqrt(np.maximum(z2, 0.0))  # vegan sets negative squared distances to zero

    # As in vegan::permutest.betadisper, the ANOVA is refitted to permuted residuals.
    f_obs = _anova_f(z, X)
    resid = z - X @ ((z @ X) / X.sum(axis=0))
    rng = np.random.default_rng(random_state)
    n_ge = 1  # the observed value counts as one permutation
    for _ in range(permutations):
        n_ge += _anova_f(resid[rng.permutation(len(z))], X) >= f_obs - 1e-12

    axes = [f"PCo{i + 1}" for i in range(pos.shape[1])]
    return {
        "distances": pd.Series(z, index=pm.sample_ids, name="distance"),
        "centroids": pd.DataFrame(np.array(cents), index=groups, columns=axes),
        "coordinates": pd.DataFrame(pos, index=pm.sample_ids, columns=axes),
        "explained": pd.Series(vals[vals > 0] / vals[vals > 0].sum(), index=axes, name="explained"),
        "F": float(f_obs),
        "p_value": n_ge / (permutations + 1),
        "column": column,
        "type": type,
        "permutations": permutations,
    }