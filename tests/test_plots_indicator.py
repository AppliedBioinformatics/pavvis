import numpy as np
import pandas as pd
import pytest

import pavvis.plots.indicator as indicator_plots
from pavvis import PavMatrix
from pavvis.stats.indicator import indicator_genes

GROUPS = ["a"] * 6 + ["b"] * 6 + ["c"] * 6


def _make_pm(tmp_path, pav: np.ndarray, groups: list, extra: dict | None = None) -> PavMatrix:
    samples = [f"s{i}" for i in range(pav.shape[1])]
    genes = [f"g{i}" for i in range(pav.shape[0])]
    pav_path = tmp_path / "pav.csv"
    meta_path = tmp_path / "metadata.csv"
    pd.DataFrame(pav, index=genes, columns=samples).to_csv(pav_path)
    pd.DataFrame({"group": groups, **(extra or {})}, index=samples).to_csv(meta_path)
    return PavMatrix(pav_path, meta_path)


def _pav() -> np.ndarray:
    """Three indicators of a (g0-g2), two of b (g3, g4), none of c, plus shared genes."""
    a, b, c = np.eye(3, dtype=int).repeat(6, axis=1)
    rng = np.random.default_rng(0)
    return np.vstack([
        a, a, np.r_[[1] * 5, [0], [0] * 12],   # g0, g1 perfect for a; g2 in 5/6 of a
        b, b,                                  # g3, g4 perfect for b
        (rng.random((6, 18)) < 0.5).astype(int),  # g5-g10 no association
        np.ones((1, 18), dtype=int),           # g11 core (not tested)
    ])


@pytest.fixture
def setup(tmp_path):
    pm = _make_pm(tmp_path, _pav(), GROUPS, extra={"depth": np.arange(18) * 1.5})
    res = indicator_genes(pm, "group", permutations=999, random_state=0)
    return pm, res


def _trace(fig, kind):
    return next(t for t in fig.data if (t.meta or {}).get("kind") == kind)


def _sig_by_group(res, q_max=0.05):
    sig = res[res.q_value < q_max]
    return sig.groupby("group").size().to_dict()


def test_rows_are_top_n_per_group(setup):
    pm, res = setup
    with pytest.warns(UserWarning, match="c"):
        fig = indicator_plots.heatmap(pm, res, "group", top_n=2)
    expected = sum(min(2, n) for n in _sig_by_group(res).values())
    assert len(_trace(fig, "frequency").y) == expected


def test_rows_grouped_and_sorted_by_stat(setup):
    pm, res = setup
    with pytest.warns(UserWarning):
        fig = indicator_plots.heatmap(pm, res, "group", top_n=10)
    genes = list(_trace(fig, "frequency").y)
    row_groups = res.loc[genes, "group"].tolist()
    assert row_groups == sorted(row_groups)
    for g in set(row_groups):
        stats = res.loc[[x for x in genes if res.loc[x, "group"] == g], "stat"]
        assert stats.is_monotonic_decreasing


def test_cells_are_group_frequencies(setup):
    pm, res = setup
    with pytest.warns(UserWarning):
        fig = indicator_plots.heatmap(pm, res, "group", top_n=10)
    tr = _trace(fig, "frequency")
    assert list(tr.x) == ["a", "b", "c"]
    z = np.array(tr.z, dtype=float)
    labels = np.array(GROUPS)
    for i, gene in enumerate(tr.y):
        row = pm.pav.loc[gene].values
        expected = [row[labels == g].mean() for g in ["a", "b", "c"]]
        assert np.allclose(z[i], expected)
        # The cell for the gene's own group is its fidelity B.
        own = ["a", "b", "c"].index(res.loc[gene, "group"])
        assert z[i, own] == pytest.approx(res.loc[gene, "B"])


def test_strip_colours_rows_by_indicator_group(setup):
    pm, res = setup
    with pytest.warns(UserWarning):
        fig = indicator_plots.heatmap(pm, res, "group", top_n=10)
    strip = _trace(fig, "strip")
    assert list(strip.y) == list(_trace(fig, "frequency").y)
    labels = [row[0][0] for row in strip.customdata]
    assert labels == res.loc[list(strip.y), "group"].tolist()


def test_q_max_filters_rows(setup):
    pm, res = setup
    with pytest.warns(UserWarning):
        loose = indicator_plots.heatmap(pm, res, "group", top_n=10, q_max=1.01)
        strict = indicator_plots.heatmap(pm, res, "group", top_n=10)
    assert len(_trace(loose, "frequency").y) > len(_trace(strict, "frequency").y)


def test_group_without_indicators_warns(setup):
    pm, res = setup
    with pytest.warns(UserWarning, match="no indicator genes.*c"):
        indicator_plots.heatmap(pm, res, "group")


def test_no_significant_genes_raises(setup):
    pm, res = setup
    with pytest.raises(ValueError, match="No indicator genes"):
        indicator_plots.heatmap(pm, res, "group", q_max=0.0)


def test_invalid_inputs_raise(setup):
    pm, res = setup
    with pytest.raises(ValueError, match="not a metadata column"):
        indicator_plots.heatmap(pm, res, "nope")
    with pytest.raises(ValueError, match="continuous"):
        indicator_plots.heatmap(pm, res, "depth")
    bad = res.copy()
    bad.loc[bad.index[0], "group"] = "zzz"
    with pytest.raises(ValueError, match="zzz"):
        indicator_plots.heatmap(pm, bad, "group")
    missing = pd.concat([res, res.iloc[[0]].rename(index={res.index[0]: "not_a_gene"})])
    with pytest.raises(ValueError, match="not_a_gene"):
        indicator_plots.heatmap(pm, missing, "group")


def test_gene_labels_toggle(setup):
    pm, res = setup
    with pytest.warns(UserWarning):
        fig = indicator_plots.heatmap(pm, res, "group", show_gene_labels=False)
    assert fig.layout.yaxis.showticklabels is False


def _gene_traces(fig):
    return {t.meta["group"]: t for t in fig.data if (t.meta or {}).get("kind") == "genes"}


def test_sf_one_panel_per_group(setup):
    pm, res = setup
    with pytest.warns(UserWarning):
        fig = indicator_plots.specificity_fidelity(pm, res, "group")
    titles = [a.text for a in fig.layout.annotations][:3]
    assert titles == ["a", "b", "c"]
    assert set(_gene_traces(fig)) == set(_sig_by_group(res))


def test_sf_points_are_group_a_and_b(setup):
    pm, res = setup
    with pytest.warns(UserWarning):
        fig = indicator_plots.specificity_fidelity(pm, res, "group")
    sig = res[res.q_value < 0.05]
    for g, tr in _gene_traces(fig).items():
        expected = sig[sig.group == g]
        assert sorted(tr.text) == sorted(expected.index)
        got = pd.DataFrame({"A": tr.x, "B": tr.y}, index=list(tr.text)).loc[expected.index]
        assert np.allclose(got.A, expected.A)
        assert np.allclose(got.B, expected.B)


def test_sf_q_max_filters_points(setup):
    pm, res = setup
    with pytest.warns(UserWarning):
        loose = indicator_plots.specificity_fidelity(pm, res, "group", q_max=1.01)
        strict = indicator_plots.specificity_fidelity(pm, res, "group")
    n = lambda f: sum(len(t.x) for t in _gene_traces(f).values())
    assert n(loose) > n(strict)


def test_sf_level_curves(setup):
    pm, res = setup
    with pytest.warns(UserWarning):
        fig = indicator_plots.specificity_fidelity(pm, res, "group", levels=(0.5, 0.8))
    levels = [t for t in fig.data if (t.meta or {}).get("kind") == "level"]
    assert len(levels) == 2 * 3
    for t in levels:
        stat = np.sqrt(np.array(t.x) * np.array(t.y))
        assert np.allclose(stat, t.meta["level"])
    with pytest.warns(UserWarning):
        none = indicator_plots.specificity_fidelity(pm, res, "group", levels=())
    assert not [t for t in none.data if (t.meta or {}).get("kind") == "level"]


def test_sf_warns_and_raises(setup):
    pm, res = setup
    with pytest.warns(UserWarning, match="no indicator genes.*c"):
        indicator_plots.specificity_fidelity(pm, res, "group")
    with pytest.raises(ValueError, match="No indicator genes"):
        indicator_plots.specificity_fidelity(pm, res, "group", q_max=0.0)
    with pytest.raises(ValueError, match="continuous"):
        indicator_plots.specificity_fidelity(pm, res, "depth")
    bad = res.copy()
    bad.loc[bad.index[0], "group"] = "zzz"
    with pytest.raises(ValueError, match="zzz"):
        indicator_plots.specificity_fidelity(pm, bad, "group")
