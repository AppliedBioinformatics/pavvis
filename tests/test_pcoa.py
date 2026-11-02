import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist, squareform

import pavvis.plots.pcoa as pcoa
from pavvis import PavMatrix
from pavvis._ordination import pcoa_eig


def _make_pm(tmp_path, pav: np.ndarray, groups: list[str] | None = None) -> PavMatrix:
    samples = [f"s{i}" for i in range(pav.shape[1])]
    genes = [f"g{i}" for i in range(pav.shape[0])]
    if groups is None:
        groups = ["a"] * len(samples)
    pav_path = tmp_path / "pav.csv"
    meta_path = tmp_path / "metadata.csv"
    pd.DataFrame(pav, index=genes, columns=samples).to_csv(pav_path)
    pd.DataFrame({"group": groups}, index=samples).to_csv(meta_path)
    return PavMatrix(pav_path, meta_path)


def _euclidean_pm(tmp_path, points: np.ndarray) -> PavMatrix:
    """PavMatrix whose cached Jaccard matrix is 1 - D for Euclidean distances D."""
    pm = _make_pm(tmp_path, _random_pav(n_samples=len(points)))
    D = squareform(pdist(points))
    pm._jaccard_matrix = pd.DataFrame(1 - D, index=pm.sample_ids, columns=pm.sample_ids)
    return pm


def _random_pav(n_samples: int = 10, n_genes: int = 50, p: float = 0.6) -> np.ndarray:
    rng = np.random.default_rng(0)
    return (rng.random((n_genes, n_samples)) < p).astype(int)


def _two_group_pav(n_per_group: int = 6) -> np.ndarray:
    """Shared core plus a block of genes present only in each group."""
    rng = np.random.default_rng(1)
    n = 2 * n_per_group
    core = np.ones((20, n), dtype=int)
    block_a = np.zeros((30, n), dtype=int)
    block_a[:, :n_per_group] = 1
    block_b = np.zeros((30, n), dtype=int)
    block_b[:, n_per_group:] = 1
    noise = (rng.random((10, n)) < 0.5).astype(int)
    return np.vstack([core, block_a, block_b, noise])


# pcoa_eig

def test_pcoa_eig_reproduces_euclidean_distances():
    points = np.random.default_rng(0).random((8, 2)) * 0.3
    D = squareform(pdist(points))
    vals, vecs = pcoa_eig(D)
    assert (vals > 0).all()
    coords = vecs * np.sqrt(vals)
    assert np.allclose(squareform(pdist(coords)), D)


def test_pcoa_eig_sorted_descending():
    D = squareform(pdist(np.random.default_rng(0).random((10, 4))))
    vals, _ = pcoa_eig(D)
    assert np.all(np.diff(vals) <= 0)


# compute_pcoa

def test_compute_pcoa_returns_self(tmp_path):
    pm = _make_pm(tmp_path, _random_pav())
    assert pm.compute_pcoa() is pm


def test_compute_pcoa_reproduces_euclidean_distances(tmp_path):
    points = np.random.default_rng(0).random((8, 2)) * 0.3
    pm = _euclidean_pm(tmp_path, points).compute_pcoa()
    emb = pm._pcoa_embedding
    assert list(emb.columns) == ["PCo1", "PCo2"]
    assert list(emb.index) == pm.sample_ids
    assert np.allclose(squareform(pdist(emb.values)), squareform(pdist(points)))


def test_compute_pcoa_explained_variance_bounded_and_non_increasing(tmp_path):
    pm = _make_pm(tmp_path, _random_pav()).compute_pcoa()
    explained = pm._pcoa_explained.values
    assert ((explained >= 0) & (explained <= 1)).all()
    assert np.all(np.diff(explained) <= 0)


def test_compute_pcoa_n_components_truncates(tmp_path):
    pm = _make_pm(tmp_path, _random_pav()).compute_pcoa(n_components=3)
    assert list(pm._pcoa_embedding.columns) == ["PCo1", "PCo2", "PCo3"]
    assert len(pm._pcoa_explained) == 3


def test_compute_pcoa_none_keeps_all_positive_axes(tmp_path):
    pm = _make_pm(tmp_path, _random_pav()).compute_pcoa(n_components=None)
    assert pm._pcoa_explained.sum() == pytest.approx(1.0)


def test_compute_pcoa_cached_for_same_params(tmp_path):
    pm = _make_pm(tmp_path, _random_pav()).compute_pcoa(n_components=3)
    emb = pm._pcoa_embedding
    pm.compute_pcoa(n_components=3)
    assert pm._pcoa_embedding is emb
    pm.compute_pcoa(n_components=2)
    assert pm._pcoa_embedding.shape[1] == 2


def test_compute_pcoa_separates_groups_on_pco1(tmp_path):
    pm = _make_pm(tmp_path, _two_group_pav()).compute_pcoa()
    pco1 = pm._pcoa_embedding["PCo1"].values
    a, b = pco1[:6], pco1[6:]
    assert a.max() < b.min() or b.max() < a.min()


def test_subset_does_not_carry_pcoa_cache(tmp_path):
    pm = _make_pm(tmp_path, _random_pav()).compute_pcoa()
    sub = pm.subset(samples=pm.sample_ids[:5])
    assert sub._pcoa_embedding is None


# plots

def test_scatter_requires_pcoa(tmp_path):
    pm = _make_pm(tmp_path, _random_pav())
    with pytest.raises(RuntimeError, match="compute_pcoa"):
        pcoa.scatter(pm)
    with pytest.raises(RuntimeError, match="compute_pcoa"):
        pcoa.scree(pm)


def test_scatter_axis_titles_show_explained_variance(tmp_path):
    pm = _make_pm(tmp_path, _random_pav()).compute_pcoa()
    fig = pcoa.scatter(pm, components=(1, 3))
    pct1 = pm._pcoa_explained["PCo1"] * 100
    assert fig.layout.xaxis.title.text == f"PCo1 ({pct1:.1f}%)"
    assert fig.layout.yaxis.title.text.startswith("PCo3 (")
    assert list(fig.data[0].y) == pm._pcoa_embedding["PCo3"].tolist()


def test_scatter_discrete_one_trace_per_group(tmp_path):
    pm = _make_pm(tmp_path, _two_group_pav(), groups=["a"] * 6 + ["b"] * 6).compute_pcoa()
    fig = pcoa.scatter(pm, color_by="group")
    assert [t.name for t in fig.data] == ["a", "b"]
    assert [len(t.x) for t in fig.data] == [6, 6]


def test_scatter_invalid_component_raises(tmp_path):
    pm = _make_pm(tmp_path, _random_pav()).compute_pcoa(n_components=2)
    with pytest.raises(ValueError, match="components"):
        pcoa.scatter(pm, components=(1, 3))


def test_umap_scatter_unchanged_after_helper_refactor(tmp_path):
    import pavvis.plots.umap as umap_plots
    pm = _make_pm(tmp_path, _random_pav(n_samples=6), groups=["a", "b"] * 3)
    pm._umap_embedding = pd.DataFrame(
        np.arange(18, dtype=float).reshape(6, 3),
        index=pm.sample_ids, columns=["UMAP1", "UMAP2", "UMAP3"],
    )
    fig = umap_plots.scatter(pm, color_by="group")
    assert [t.name for t in fig.data] == ["a", "b"]
    assert list(fig.data[0].x) == [0.0, 6.0, 12.0]
    assert list(fig.data[0].text) == ["s0", "s2", "s4"]
    assert fig.layout.xaxis.title.text == "UMAP1"


def test_scree_one_bar_per_axis(tmp_path):
    pm = _make_pm(tmp_path, _random_pav()).compute_pcoa(n_components=4)
    fig = pcoa.scree(pm)
    assert list(fig.data[0].x) == ["PCo1", "PCo2", "PCo3", "PCo4"]
    assert np.allclose(fig.data[0].y, pm._pcoa_explained.values * 100)