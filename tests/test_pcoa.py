import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist, squareform

import pavvis.plots.pcoa as pcoa
from pavvis import PavMatrix
from pavvis._ordination import pcoa_eig


def _make_pm(tmp_path, pav: np.ndarray, groups: list[str] | None = None, extra: dict | None = None) -> PavMatrix:
    samples = [f"s{i}" for i in range(pav.shape[1])]
    genes = [f"g{i}" for i in range(pav.shape[0])]
    if groups is None:
        groups = ["a"] * len(samples)
    pav_path = tmp_path / "pav.csv"
    meta_path = tmp_path / "metadata.csv"
    pd.DataFrame(pav, index=genes, columns=samples).to_csv(pav_path)
    pd.DataFrame({"group": groups, **(extra or {})}, index=samples).to_csv(meta_path)
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

# centroids / dispersion_boxplot

GROUPS_6_6 = ["a"] * 6 + ["b"] * 6


def _traces(fig, kind):
    return [t for t in fig.data if (t.meta or {}).get("kind") == kind]


def test_centroids_one_centroid_trace_per_group(tmp_path):
    pm = _make_pm(tmp_path, _two_group_pav(), groups=GROUPS_6_6)
    fig = pcoa.centroids(pm, "group")
    cents = _traces(fig, "centroid")
    assert [t.legendgroup for t in cents] == ["a", "b"]
    assert all(len(t.x) == 1 for t in cents)


def test_centroids_spider_has_three_points_per_sample(tmp_path):
    pm = _make_pm(tmp_path, _two_group_pav(), groups=GROUPS_6_6)
    fig = pcoa.centroids(pm, "group")
    spiders = _traces(fig, "spider")
    assert [len(t.x) for t in spiders] == [18, 18]


def test_centroids_match_betadisper(tmp_path):
    from pavvis.stats.permanova import betadisper
    pm = _make_pm(tmp_path, _two_group_pav(), groups=GROUPS_6_6)
    fig = pcoa.centroids(pm, "group", components=(2, 1), type="centroid")
    bd = betadisper(pm, "group", type="centroid", permutations=0)
    for t in _traces(fig, "centroid"):
        assert t.x[0] == pytest.approx(bd["centroids"].loc[t.legendgroup, "PCo2"])
        assert t.y[0] == pytest.approx(bd["centroids"].loc[t.legendgroup, "PCo1"])
    assert fig.layout.xaxis.title.text.startswith("PCo2 (")


def test_centroids_hull_toggle_and_no_annotation_by_default(tmp_path):
    pm = _make_pm(tmp_path, _two_group_pav(), groups=GROUPS_6_6)
    fig = pcoa.centroids(pm, "group")
    assert len(_traces(fig, "hull")) == 2
    assert not fig.layout.annotations
    fig = pcoa.centroids(pm, "group", hull=False)
    assert not _traces(fig, "hull")


def test_centroids_annotates_supplied_results(tmp_path):
    from pavvis.stats.permanova import betadisper, permanova
    pm = _make_pm(tmp_path, _two_group_pav(), groups=GROUPS_6_6)
    bd = betadisper(pm, "group", type="centroid", permutations=19, random_state=0)
    pa = permanova(pm, "group", permutations=19, random_state=0)
    fig = pcoa.centroids(pm, "group", betadisper_result=bd, permanova_result=pa)
    text = fig.layout.annotations[0].text
    model = pa.loc["Model"]
    assert f"R² = {model['R2']:.3f}" in text and f"p = {model['Pr(>F)']:.3g}" in text
    assert f"PERMDISP (centroid): F = {bd['F']:.2f}, p = {bd['p_value']:.3g}" in text
    assert text.endswith("19 permutations")
    # The supplied ordination is drawn: centres use the result's type.
    for t in _traces(fig, "centroid"):
        assert t.x[0] == pytest.approx(bd["centroids"].loc[t.legendgroup, "PCo1"])

    only_permanova = pcoa.centroids(pm, "group", permanova_result=pa).layout.annotations[0].text
    assert "PERMANOVA" in only_permanova and "PERMDISP" not in only_permanova

    pa99 = permanova(pm, "group", permutations=29, random_state=0)
    mixed = pcoa.centroids(pm, "group", betadisper_result=bd, permanova_result=pa99).layout.annotations[0].text
    assert "PERMANOVA 29" in mixed and "PERMDISP 19" in mixed


def test_centroids_rejects_mismatched_results(tmp_path):
    from pavvis.stats.permanova import betadisper, permanova
    pm = _make_pm(tmp_path, _two_group_pav(), groups=GROUPS_6_6, extra={"other": ["x", "y"] * 6})
    bd = betadisper(pm, "group", permutations=0)
    with pytest.raises(ValueError, match="computed for 'group'"):
        pcoa.centroids(pm, "other", betadisper_result=bd)
    with pytest.raises(ValueError, match="type='median'"):
        pcoa.centroids(pm, "group", betadisper_result=bd, type="centroid")
    with pytest.raises(ValueError, match="computed for 'other'"):
        pcoa.centroids(pm, "group", permanova_result=permanova(pm, "other", permutations=0))
    sub = pm.subset(samples=pm.sample_ids[:10])
    with pytest.raises(ValueError, match="samples do not match"):
        pcoa.centroids(sub, "group", betadisper_result=bd)


def test_centroids_skips_hull_for_small_groups(tmp_path):
    pm = _make_pm(tmp_path, _random_pav(n_samples=8), groups=["a"] * 6 + ["b"] * 2)
    fig = pcoa.centroids(pm, "group")
    assert [t.legendgroup for t in _traces(fig, "hull")] == ["a"]


def test_centroids_invalid_inputs_raise(tmp_path):
    pm = _make_pm(tmp_path, _two_group_pav(), groups=GROUPS_6_6)
    with pytest.raises(ValueError, match="not a metadata column"):
        pcoa.centroids(pm, "nope")
    with pytest.raises(ValueError, match="components"):
        pcoa.centroids(pm, "group", components=(1, 99))


def test_dispersion_boxplot_one_box_per_group(tmp_path):
    from pavvis.stats.permanova import betadisper
    pm = _make_pm(tmp_path, _two_group_pav(), groups=GROUPS_6_6)
    fig = pcoa.dispersion_boxplot(pm, "group")
    assert [t.name for t in fig.data] == ["a", "b"]
    z = betadisper(pm, "group", permutations=0)["distances"]
    assert np.allclose(sorted(fig.data[0].y), sorted(z.iloc[:6]))


def test_dispersion_boxplot_uses_supplied_result(tmp_path):
    from pavvis.stats.permanova import betadisper
    pm = _make_pm(tmp_path, _two_group_pav(), groups=GROUPS_6_6)
    bd = betadisper(pm, "group", type="centroid", permutations=0)
    fig = pcoa.dispersion_boxplot(pm, "group", betadisper_result=bd)
    assert np.allclose(sorted(fig.data[0].y), sorted(bd["distances"].iloc[:6]))
    assert fig.layout.title.text == "Distance to group centroid"


def test_centroids_trace_order_per_group(tmp_path):
    pm = _make_pm(tmp_path, _two_group_pav(), groups=GROUPS_6_6)
    fig = pcoa.centroids(pm, "group")
    assert [t.meta["kind"] for t in fig.data] == ["hull", "spider", "points", "centroid"] * 2
