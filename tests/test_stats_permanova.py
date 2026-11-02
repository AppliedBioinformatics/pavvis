import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f_oneway

from pavvis import PavMatrix
from pavvis.stats.permanova import betadisper, permanova


def _make_pm(tmp_path, pav: np.ndarray, groups: list, extra: dict | None = None) -> PavMatrix:
    samples = [f"s{i}" for i in range(pav.shape[1])]
    genes = [f"g{i}" for i in range(pav.shape[0])]
    pav_path = tmp_path / "pav.csv"
    meta_path = tmp_path / "metadata.csv"
    pd.DataFrame(pav, index=genes, columns=samples).to_csv(pav_path)
    pd.DataFrame({"group": groups, **(extra or {})}, index=samples).to_csv(meta_path)
    return PavMatrix(pav_path, meta_path)


def _random_pav(n_samples: int, n_genes: int = 60, p: float = 0.5, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return (rng.random((n_genes, n_samples)) < p).astype(int)


def _euclidean_pm(tmp_path, points: np.ndarray, groups: list) -> PavMatrix:
    """PavMatrix whose cached Jaccard matrix is 1 - D for Euclidean distances D."""
    pm = _make_pm(tmp_path, _random_pav(len(points)), groups)
    D = squareform(pdist(points))
    pm._jaccard_matrix = pd.DataFrame(1 - D, index=pm.sample_ids, columns=pm.sample_ids)
    return pm


def _two_group_pav(n_per_group: int = 6) -> np.ndarray:
    rng = np.random.default_rng(1)
    n = 2 * n_per_group
    core = np.ones((20, n), dtype=int)
    block_a = np.zeros((30, n), dtype=int)
    block_a[:, :n_per_group] = 1
    block_b = np.zeros((30, n), dtype=int)
    block_b[:, n_per_group:] = 1
    noise = (rng.random((10, n)) < 0.5).astype(int)
    return np.vstack([core, block_a, block_b, noise])


GROUPS_3 = ["a"] * 4 + ["b"] * 5 + ["c"] * 3


@pytest.fixture
def points():
    rng = np.random.default_rng(0)
    offsets = np.repeat([[0.0, 0.0], [0.2, 0.0], [0.0, 0.2]], [4, 5, 3], axis=0)
    return rng.random((12, 2)) * 0.1 + offsets


# permanova

def test_permanova_sums_of_squares_on_euclidean_distances(tmp_path, points):
    pm = _euclidean_pm(tmp_path, points, GROUPS_3)
    tab = permanova(pm, "group", permutations=99, random_state=0)
    labels = np.array(GROUPS_3)
    ss_w = sum(((points[labels == g] - points[labels == g].mean(0)) ** 2).sum() for g in "abc")
    ss_t = ((points - points.mean(0)) ** 2).sum()
    assert tab.loc["Residual", "SumOfSqs"] == pytest.approx(ss_w)
    assert tab.loc["Total", "SumOfSqs"] == pytest.approx(ss_t)
    assert tab.loc["Model", "SumOfSqs"] == pytest.approx(ss_t - ss_w)
    assert tab.loc["Model", "R2"] == pytest.approx((ss_t - ss_w) / ss_t)
    f = ((ss_t - ss_w) / 2) / (ss_w / 9)
    assert tab.loc["Model", "F"] == pytest.approx(f)


def test_permanova_table_layout(tmp_path, points):
    pm = _euclidean_pm(tmp_path, points, GROUPS_3)
    tab = permanova(pm, "group", permutations=9, random_state=0)
    assert list(tab.index) == ["Model", "Residual", "Total"]
    assert list(tab.columns) == ["Df", "SumOfSqs", "R2", "F", "Pr(>F)"]
    assert list(tab["Df"]) == [2, 9, 11]
    assert tab.loc["Total", "R2"] == pytest.approx(1.0)
    assert tab.loc["Model", "R2"] + tab.loc["Residual", "R2"] == pytest.approx(1.0)
    assert tab[["F", "Pr(>F)"]].iloc[1:].isna().all().all()


def test_permanova_separated_groups_significant(tmp_path):
    pm = _make_pm(tmp_path, _two_group_pav(), ["a"] * 6 + ["b"] * 6)
    tab = permanova(pm, "group", permutations=999, random_state=0)
    assert tab.loc["Model", "Pr(>F)"] < 0.05


def test_permanova_random_groups_not_significant(tmp_path):
    pm = _make_pm(tmp_path, _random_pav(12), ["a", "b"] * 6)
    tab = permanova(pm, "group", permutations=999, random_state=0)
    assert tab.loc["Model", "Pr(>F)"] > 0.05


def test_permanova_reproducible(tmp_path):
    pm = _make_pm(tmp_path, _random_pav(12), ["a", "b"] * 6)
    a = permanova(pm, "group", permutations=99, random_state=3)
    b = permanova(pm, "group", permutations=99, random_state=3)
    pd.testing.assert_frame_equal(a, b)


def test_single_group_raises(tmp_path):
    pm = _make_pm(tmp_path, _random_pav(6), ["a"] * 6)
    with pytest.raises(ValueError, match="at least 2 groups"):
        permanova(pm, "group", permutations=9)
    with pytest.raises(ValueError, match="at least 2 groups"):
        betadisper(pm, "group", permutations=9)


def test_invalid_column_raises(tmp_path):
    pm = _make_pm(tmp_path, _random_pav(12), ["a", "b"] * 6, extra={"depth": np.arange(12) * 1.5})
    with pytest.raises(ValueError, match="continuous"):
        permanova(pm, "depth", permutations=9)
    with pytest.raises(ValueError, match="not a metadata column"):
        betadisper(pm, "nope", permutations=9)


# betadisper

def test_betadisper_centroid_matches_euclidean_distance_to_mean(tmp_path, points):
    pm = _euclidean_pm(tmp_path, points, GROUPS_3)
    res = betadisper(pm, "group", type="centroid", permutations=9, random_state=0)
    labels = np.array(GROUPS_3)
    expected = np.empty(len(points))
    for g in "abc":
        m = labels == g
        expected[m] = np.linalg.norm(points[m] - points[m].mean(0), axis=1)
    assert np.allclose(res["distances"].values, expected)
    assert list(res["distances"].index) == pm.sample_ids


def test_betadisper_median_minimises_summed_distance(tmp_path, points):
    pm = _euclidean_pm(tmp_path, points, GROUPS_3)
    med = betadisper(pm, "group", type="median", permutations=9, random_state=0)["distances"]
    cen = betadisper(pm, "group", type="centroid", permutations=9, random_state=0)["distances"]
    labels = pd.Series(GROUPS_3, index=pm.sample_ids)
    assert (med.groupby(labels).sum() <= cen.groupby(labels).sum() + 1e-9).all()


def test_betadisper_f_matches_one_way_anova(tmp_path, points):
    pm = _euclidean_pm(tmp_path, points, GROUPS_3)
    res = betadisper(pm, "group", permutations=99, random_state=0)
    z, labels = res["distances"].values, np.array(GROUPS_3)
    expected = f_oneway(*(z[labels == g] for g in "abc")).statistic
    assert res["F"] == pytest.approx(expected)
    assert 0 < res["p_value"] <= 1


def test_betadisper_centroids_frame(tmp_path, points):
    pm = _euclidean_pm(tmp_path, points, GROUPS_3)
    res = betadisper(pm, "group", type="centroid", permutations=9, random_state=0)
    cents = res["centroids"]
    assert list(cents.index) == ["a", "b", "c"]
    assert list(cents.columns) == ["PCo1", "PCo2"]
    # Centroids are the group means of the sample coordinates from compute_pcoa.
    emb = pm.compute_pcoa(n_components=None)._pcoa_embedding
    assert np.allclose(cents.values, emb.groupby(np.array(GROUPS_3)).mean().values)


def test_betadisper_jaccard_with_negative_eigenvalues(tmp_path):
    pm = _make_pm(tmp_path, _random_pav(12), GROUPS_3)
    res = betadisper(pm, "group", permutations=99, random_state=0)
    z = res["distances"].values
    assert np.isfinite(z).all() and (z >= 0).all()


def test_betadisper_detects_unequal_dispersion(tmp_path):
    rng = np.random.default_rng(0)
    tight = rng.normal(0, 0.01, (10, 2))
    spread = rng.normal(0, 0.15, (10, 2)) + 0.05
    pm = _euclidean_pm(tmp_path, np.vstack([tight, spread]), ["a"] * 10 + ["b"] * 10)
    res = betadisper(pm, "group", permutations=999, random_state=0)
    assert res["p_value"] < 0.05


def test_betadisper_equal_dispersion_not_significant(tmp_path):
    rng = np.random.default_rng(0)
    a = rng.normal(0, 0.05, (10, 2))
    b = rng.normal(0, 0.05, (10, 2)) + 0.3
    pm = _euclidean_pm(tmp_path, np.vstack([a, b]), ["a"] * 10 + ["b"] * 10)
    res = betadisper(pm, "group", permutations=999, random_state=0)
    assert res["p_value"] > 0.05


def test_betadisper_invalid_type_raises(tmp_path, points):
    pm = _euclidean_pm(tmp_path, points, GROUPS_3)
    with pytest.raises(ValueError, match="type"):
        betadisper(pm, "group", type="mean")

def test_betadisper_returns_coordinates_and_explained(tmp_path):
    pm = _make_pm(tmp_path, _random_pav(12), GROUPS_3)
    res = betadisper(pm, "group", permutations=0)
    emb = pm.compute_pcoa(n_components=None)._pcoa_embedding
    assert np.allclose(res["coordinates"].values, emb.values)
    assert list(res["coordinates"].index) == pm.sample_ids
    assert np.allclose(res["explained"].values, pm._pcoa_explained.values)
