from __future__ import annotations

import numpy as np


def pcoa_eig(D: np.ndarray, tol: float = 1e-10) -> tuple[np.ndarray, np.ndarray]:
    """Eigendecompose the Gower-centred matrix of a distance matrix (classical MDS).

    Non-Euclidean distances such as Jaccard produce negative eigenvalues. These are
    returned alongside the positive ones because betadisper needs both.

    Args:
        D: Symmetric distance matrix of shape (samples, samples).
        tol: Eigenvalues with magnitude below tol × the largest magnitude are dropped
            as numerically zero.

    Returns:
        Tuple of (eigvals, vectors), sorted by eigenvalue in descending order.
        vectors has shape (samples, axes) with unit-length columns.
    """
    n = D.shape[0]
    C = np.eye(n) - 1.0 / n
    B = -0.5 * C @ (D ** 2) @ C
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    keep = np.abs(vals) > tol * np.abs(vals).max()
    return vals[keep], vecs[:, keep]