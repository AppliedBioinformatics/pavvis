from __future__ import annotations

import numpy as np


def accumulation_runs(
    pav: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Accumulate pangenome and core sizes over random sample orderings.

    Args:
        pav: Binary array of shape (genes, samples).
        permutations: Number of random sample orderings to draw.
        rng: NumPy random generator used to draw the orderings.

    Returns:
        Tuple of (pan_runs, core_runs), each of shape (permutations, samples).
        Entry [i, k] is the pangenome (union) or core (intersection) size after
        adding the first k + 1 samples of permutation i.
    """
    n_samples = pav.shape[1]
    pan_runs = np.empty((permutations, n_samples), dtype=np.int32)
    core_runs = np.empty((permutations, n_samples), dtype=np.int32)

    for i in range(permutations):
        order = rng.permutation(n_samples)
        shuffled = pav[:, order]
        pan_runs[i] = (np.cumsum(shuffled, axis=1) > 0).sum(axis=0)
        core_runs[i] = np.cumprod(shuffled, axis=1).sum(axis=0)

    return pan_runs, core_runs
