import numpy as np
import pandas as pd
import pytest

from pavvis import PavMatrix
from pavvis.stats.indicator import indicator_genes


def _make_pm(tmp_path, pav: np.ndarray, groups: list, extra: dict | None = None) -> PavMatrix:
    samples = [f"s{i}" for i in range(pav.shape[1])]
    genes = [f"g{i}" for i in range(pav.shape[0])]
    pav_path = tmp_path / "pav.csv"
    meta_path = tmp_path / "metadata.csv"
    pd.DataFrame(pav, index=genes, columns=samples).to_csv(pav_path)
    pd.DataFrame({"group": groups, **(extra or {})}, index=samples).to_csv(meta_path)
    return PavMatrix(pav_path, meta_path)


GROUPS_6_6 = ["a"] * 6 + ["b"] * 6


def _balanced_pav() -> np.ndarray:
    return np.array([
        [1] * 6 + [0] * 6,          # g0: perfect indicator of a
        [0] * 6 + [1] * 6,          # g1: perfect indicator of b
        [1] * 12,                   # g2: core
        [1, 0] * 6,                 # g3: no association
        [1, 1, 1, 0, 0, 0] * 2,     # g4: no association
    ])


def test_perfect_indicator_has_unit_components(tmp_path):
    pm = _make_pm(tmp_path, _balanced_pav(), GROUPS_6_6)
    res = indicator_genes(pm, "group", permutations=99, random_state=0)
    for gene, group in (("g0", "a"), ("g1", "b")):
        row = res.loc[gene]
        assert row["group"] == group
        assert row["A"] == pytest.approx(1.0)
        assert row["B"] == pytest.approx(1.0)
        assert row["stat"] == pytest.approx(1.0)


def test_core_genes_excluded_by_default(tmp_path):
    pm = _make_pm(tmp_path, _balanced_pav(), GROUPS_6_6)
    res = indicator_genes(pm, "group", permutations=9, random_state=0)
    assert "g2" not in res.index
    assert set(res.index) == set(pm.variable_genes)


def test_genes_argument_restricts_tested_genes(tmp_path):
    pm = _make_pm(tmp_path, _balanced_pav(), GROUPS_6_6)
    res = indicator_genes(pm, "group", genes=["g0", "g3"], permutations=9, random_state=0)
    assert sorted(res.index) == ["g0", "g3"]


def test_unknown_gene_raises(tmp_path):
    pm = _make_pm(tmp_path, _balanced_pav(), GROUPS_6_6)
    with pytest.raises(ValueError, match="not found"):
        indicator_genes(pm, "group", genes=["nope"], permutations=9)


def test_unequal_group_sizes_use_group_equalised_specificity(tmp_path):
    # Present in 2/2 samples of a and 2/6 of b. IndVal.g: A = 1 / (1 + 1/3) = 0.75.
    # The raw-count version would give A = 2/4 = 0.5.
    pav = np.array([[1, 1, 1, 1, 0, 0, 0, 0]])
    pm = _make_pm(tmp_path, pav, ["a"] * 2 + ["b"] * 6)
    row = indicator_genes(pm, "group", permutations=9, random_state=0).loc["g0"]
    assert row["group"] == "a"
    assert row["A"] == pytest.approx(0.75)
    assert row["B"] == pytest.approx(1.0)
    assert row["stat"] == pytest.approx(np.sqrt(0.75))


def test_columns_and_sorting(tmp_path):
    pm = _make_pm(tmp_path, _balanced_pav(), GROUPS_6_6)
    res = indicator_genes(pm, "group", permutations=9, random_state=0)
    assert list(res.columns) == ["group", "A", "B", "stat", "p_value", "q_value"]
    assert res.index.name == "gene_id"
    assert res["stat"].is_monotonic_decreasing


def test_reproducible_with_random_state(tmp_path):
    pm = _make_pm(tmp_path, _balanced_pav(), GROUPS_6_6)
    a = indicator_genes(pm, "group", permutations=99, random_state=1)
    b = indicator_genes(pm, "group", permutations=99, random_state=1)
    pd.testing.assert_frame_equal(a, b)


def test_perfect_indicator_is_significant(tmp_path):
    pm = _make_pm(tmp_path, _balanced_pav(), GROUPS_6_6)
    res = indicator_genes(pm, "group", permutations=999, random_state=0)
    assert res.loc["g0", "p_value"] < 0.05
    assert res.loc["g3", "p_value"] > 0.05
    assert (res["q_value"] >= res["p_value"]).all()


def test_p_value_bounds(tmp_path):
    pm = _make_pm(tmp_path, _balanced_pav(), GROUPS_6_6)
    res = indicator_genes(pm, "group", permutations=99, random_state=0)
    assert ((res["p_value"] >= 1 / 100) & (res["p_value"] <= 1)).all()


def test_invalid_column_raises(tmp_path):
    pm = _make_pm(tmp_path, _balanced_pav(), GROUPS_6_6, extra={"depth": np.arange(12) * 1.5})
    with pytest.raises(ValueError, match="not a metadata column"):
        indicator_genes(pm, "nope", permutations=9)
    with pytest.raises(ValueError, match="continuous"):
        indicator_genes(pm, "depth", permutations=9)
