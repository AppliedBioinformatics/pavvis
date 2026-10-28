import numpy as np
import pandas as pd
import pytest

import pavvis.plots.curves as curves
from pavvis import PavMatrix
from pavvis._accumulation import accumulation_runs


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


def _private_pav(n_samples: int, n_core: int = 50, n_private: int = 5) -> np.ndarray:
    """Core genes plus a block of private genes unique to each sample."""
    core = np.ones((n_core, n_samples), dtype=int)
    private = np.kron(np.eye(n_samples, dtype=int), np.ones((n_private, 1), dtype=int))
    return np.vstack([core, private])


def _random_pav(n_samples: int = 10, n_genes: int = 50, p: float = 0.6) -> np.ndarray:
    rng = np.random.default_rng(0)
    return (rng.random((n_genes, n_samples)) < p).astype(int)


def _mean_traces(fig):
    return [t for t in fig.data if t.showlegend is not False]


# accumulation_runs

def test_accumulation_runs_shapes_and_bounds():
    pav = _private_pav(n_samples=6)
    pan, core = accumulation_runs(pav, permutations=4, rng=np.random.default_rng(1))
    assert pan.shape == (4, 6)
    assert core.shape == (4, 6)
    assert (pan[:, -1] == pav.shape[0]).all()
    assert (core[:, -1] == 50).all()
    assert (np.diff(pan, axis=1) >= 0).all()
    assert (np.diff(core, axis=1) <= 0).all()


def test_accumulation_runs_is_reproducible_with_seed():
    pav = _random_pav()
    a = accumulation_runs(pav, 5, np.random.default_rng(3))
    b = accumulation_runs(pav, 5, np.random.default_rng(3))
    assert all((x == y).all() for x, y in zip(a, b))


# Refactored curves keep their output

def test_legacy_curve_output_unchanged(tmp_path):
    pm = _make_pm(tmp_path, _private_pav(n_samples=8))
    fig = curves.legacy_pangenome_curve(pm, permutations=5, random_state=0)
    assert len(fig.data) == 6
    pan_mean, core_mean = fig.data[-2], fig.data[-1]
    assert list(pan_mean.x) == list(range(1, 9))
    assert list(pan_mean.y) == [50.0 + 5 * k for k in range(1, 9)]
    # A single genome's core is all of its genes; after that only the 50 core genes remain.
    assert list(core_mean.y) == [55.0] + [50.0] * 7


def test_variable_gene_curve_output_unchanged(tmp_path):
    pm = _make_pm(tmp_path, _private_pav(n_samples=8))
    fig = curves.variable_gene_curve(pm, permutations=5, random_state=0)
    assert len(fig.data) == 6
    core_mean, var_mean = fig.data[-2], fig.data[-1]
    assert list(core_mean.y) == [55.0] + [50.0] * 7
    assert list(var_mean.y) == [0.0] + [5.0 * k for k in range(2, 9)]


def test_refactored_curves_are_reproducible_with_seed(tmp_path):
    pm = _make_pm(tmp_path, _random_pav())
    for fn in (curves.legacy_pangenome_curve, curves.variable_gene_curve):
        a = fn(pm, permutations=5, random_state=1)
        b = fn(pm, permutations=5, random_state=1)
        assert list(a.data[-2].y) == list(b.data[-2].y)


# grouped_pangenome_curve

GROUPS = ["x"] * 4 + ["y"] * 3 + ["z"] * 2


def test_grouped_curve_has_one_mean_trace_per_group(tmp_path):
    pm = _make_pm(tmp_path, _private_pav(n_samples=9), GROUPS)
    fig = curves.grouped_pangenome_curve(pm, column="group", permutations=5, random_state=0)
    names = [t.name for t in _mean_traces(fig)]
    assert names == ["x (n=4)", "y (n=3)", "z (n=2)"]


def test_grouped_curve_x_length_matches_group_size(tmp_path):
    pm = _make_pm(tmp_path, _private_pav(n_samples=9), GROUPS)
    fig = curves.grouped_pangenome_curve(pm, column="group", permutations=5, random_state=0)
    for trace, n in zip(_mean_traces(fig), [4, 3, 2]):
        assert list(trace.x) == list(range(1, n + 1))


def test_grouped_curve_values_per_curve_type(tmp_path):
    pm = _make_pm(tmp_path, _private_pav(n_samples=9), GROUPS)
    kwargs = dict(column="group", groups=["x"], permutations=5, random_state=0)
    pan = _mean_traces(curves.grouped_pangenome_curve(pm, curve="pangenome", **kwargs))[0]
    core = _mean_traces(curves.grouped_pangenome_curve(pm, curve="core", **kwargs))[0]
    var = _mean_traces(curves.grouped_pangenome_curve(pm, curve="variable", **kwargs))[0]
    assert list(pan.y) == [55.0, 60.0, 65.0, 70.0]
    assert list(core.y) == [55.0, 50.0, 50.0, 50.0]
    assert list(var.y) == [0.0, 10.0, 15.0, 20.0]


def test_grouped_curve_respects_groups_order(tmp_path):
    pm = _make_pm(tmp_path, _private_pav(n_samples=9), GROUPS)
    fig = curves.grouped_pangenome_curve(
        pm, column="group", groups=["z", "x"], permutations=5, random_state=0,
    )
    assert [t.name for t in _mean_traces(fig)] == ["z (n=2)", "x (n=4)"]


def test_grouped_curve_legend_toggles_band_with_mean(tmp_path):
    pm = _make_pm(tmp_path, _private_pav(n_samples=9), GROUPS)
    fig = curves.grouped_pangenome_curve(pm, column="group", permutations=5, random_state=0)
    assert len(fig.data) == 9
    for i in range(0, 9, 3):
        band_upper, band_lower, mean = fig.data[i:i + 3]
        assert band_upper.legendgroup == band_lower.legendgroup == mean.legendgroup


def test_grouped_curve_is_reproducible_with_seed(tmp_path):
    pm = _make_pm(tmp_path, _random_pav(n_samples=9), GROUPS)
    a = curves.grouped_pangenome_curve(pm, column="group", permutations=5, random_state=2)
    b = curves.grouped_pangenome_curve(pm, column="group", permutations=5, random_state=2)
    assert [list(t.y) for t in a.data] == [list(t.y) for t in b.data]


def test_grouped_curve_skips_single_sample_group_with_warning(tmp_path):
    groups = ["x"] * 4 + ["y"] * 4 + ["solo"]
    pm = _make_pm(tmp_path, _private_pav(n_samples=9), groups)
    with pytest.warns(UserWarning, match="solo"):
        fig = curves.grouped_pangenome_curve(pm, column="group", permutations=5, random_state=0)
    assert [t.name for t in _mean_traces(fig)] == ["x (n=4)", "y (n=4)"]


def test_grouped_curve_unknown_group_raises(tmp_path):
    pm = _make_pm(tmp_path, _private_pav(n_samples=9), GROUPS)
    with pytest.raises(ValueError, match="Groups not found"):
        curves.grouped_pangenome_curve(pm, column="group", groups=["x", "nope"])


def test_grouped_curve_unknown_column_raises(tmp_path):
    pm = _make_pm(tmp_path, _private_pav(n_samples=9), GROUPS)
    with pytest.raises(ValueError, match="not a metadata column"):
        curves.grouped_pangenome_curve(pm, column="missing")


def test_grouped_curve_invalid_curve_raises(tmp_path):
    pm = _make_pm(tmp_path, _private_pav(n_samples=9), GROUPS)
    with pytest.raises(ValueError, match="curve"):
        curves.grouped_pangenome_curve(pm, column="group", curve="other")


def test_grouped_curve_no_plottable_groups_raises(tmp_path):
    pm = _make_pm(tmp_path, _private_pav(n_samples=3), ["a", "b", "c"])
    with pytest.warns(UserWarning), pytest.raises(ValueError, match="at least 2 samples"):
        curves.grouped_pangenome_curve(pm, column="group")
