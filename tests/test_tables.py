import numpy as np
import pandas as pd
import pytest

from pavvis import PavMatrix


# Frequencies with 5 samples: core 1.0, soft 0.8, disp 0.6, priv_a 0.4, single_* 0.2, absent 0.
PAV_CSV = """\
,s1,s2,s3,s4,s5
g_core,1,1,1,1,1
g_soft,1,1,1,1,0
g_disp,1,0,1,0,1
g_priv_a,1,1,0,0,0
g_single_c,0,0,0,0,1
g_single_b,0,0,1,0,0
g_absent,0,0,0,0,0
"""

META_CSV = """\
,group,depth
s1,A,10.0
s2,A,11.0
s3,B,12.0
s4,B,13.0
s5,C,14.0
"""

GENES = ["g_core", "g_soft", "g_disp", "g_priv_a", "g_single_c", "g_single_b", "g_absent"]
SAMPLES = ["s1", "s2", "s3", "s4", "s5"]


def _make_pm(tmp_path, private_min: float = 0.0) -> PavMatrix:
    pav = tmp_path / "pav.csv"
    meta = tmp_path / "metadata.csv"
    pav.write_text(PAV_CSV)
    meta.write_text(META_CSV)
    return PavMatrix(pav, meta, soft_core_min=0.8, dispensable_min=0.5, private_min=private_min)


@pytest.fixture
def pm(tmp_path):
    return _make_pm(tmp_path)


# gene_table

def test_gene_table_index_and_columns(pm):
    gt = pm.gene_table()
    assert list(gt.index) == GENES
    assert gt.index.name == "gene_id"
    assert list(gt.columns) == ["n_present", "frequency", "category"]


def test_gene_table_values(pm):
    gt = pm.gene_table()
    assert gt["n_present"].tolist() == [5, 4, 3, 2, 1, 1, 0]
    assert gt["frequency"].tolist() == pytest.approx([1.0, 0.8, 0.6, 0.4, 0.2, 0.2, 0.0])
    assert gt["category"].tolist() == [
        "core", "soft_core", "dispensable", "private", "private", "private", "absent",
    ]


def test_gene_table_categories_match_gene_counts(pm):
    counts = pm.gene_table()["category"].value_counts().to_dict()
    assert counts == {k: v for k, v in pm.gene_counts.items() if v}


def test_gene_table_categories_match_properties(pm):
    gt = pm.gene_table()
    for category, genes in [
        ("core", pm.core_genes), ("soft_core", pm.soft_core_genes),
        ("dispensable", pm.dispensable_genes), ("private", pm.private_genes),
        ("absent", pm.absent_genes),
    ]:
        assert list(gt.index[gt["category"] == category]) == genes


def test_gene_table_unclassified_below_private_min(tmp_path):
    pm = _make_pm(tmp_path, private_min=0.3)
    gt = pm.gene_table()
    assert gt.loc["g_single_c", "category"] == "unclassified"
    assert gt.loc["g_single_b", "category"] == "unclassified"
    assert gt.loc["g_priv_a", "category"] == "private"


def test_gene_table_with_column_adds_group_counts(pm):
    gt = pm.gene_table(column="group")
    assert list(gt.columns) == [
        "n_present", "frequency", "category", "n_A", "n_B", "n_C", "exclusive_to",
    ]
    assert gt["n_A"].tolist() == [2, 2, 1, 2, 0, 0, 0]
    assert gt["n_B"].tolist() == [2, 2, 1, 0, 0, 1, 0]
    assert gt["n_C"].tolist() == [1, 0, 1, 0, 1, 0, 0]


def test_gene_table_exclusive_to(pm):
    excl = pm.gene_table(column="group")["exclusive_to"]
    assert excl["g_priv_a"] == "A"
    assert excl["g_single_b"] == "B"
    assert excl["g_single_c"] == "C"
    assert excl[["g_core", "g_soft", "g_disp", "g_absent"]].isna().all()


def test_gene_table_exclusive_counts_match_exclusive_gene_counts(pm):
    excl = pm.gene_table(column="group")["exclusive_to"]
    counts = excl.value_counts().to_dict()
    expected = {g: n for g, n in pm.exclusive_gene_counts("group").items() if n}
    assert counts == expected


def test_gene_table_invalid_column_raises(pm):
    with pytest.raises(ValueError, match="not a metadata column"):
        pm.gene_table(column="missing")


def test_gene_table_continuous_column_raises(tmp_path):
    pav = tmp_path / "pav.csv"
    meta = tmp_path / "metadata.csv"
    samples = [f"s{i}" for i in range(12)]
    pd.DataFrame([[1] * 12], index=["g"], columns=samples).to_csv(pav)
    pd.DataFrame({"depth": np.arange(12, dtype=float)}, index=samples).to_csv(meta)
    pm = PavMatrix(pav, meta)
    with pytest.raises(ValueError, match="continuous"):
        pm.gene_table(column="depth")


# exclusive_gene_counts (refactored onto the shared group matrix)

def test_exclusive_gene_counts_unchanged(pm):
    assert pm.exclusive_gene_counts("group") == {"A": 1, "B": 1, "C": 1}


# sample_table

def test_sample_table_index_and_columns(pm):
    st = pm.sample_table()
    assert list(st.index) == SAMPLES
    assert st.index.name == "sample"
    assert list(st.columns) == [
        "n_genes", "n_core", "n_soft_core", "n_dispensable", "n_private", "n_singletons",
    ]


def test_sample_table_values(pm):
    st = pm.sample_table()
    assert st["n_genes"].tolist() == [4, 3, 4, 2, 3]
    assert st["n_core"].tolist() == [1, 1, 1, 1, 1]
    assert st["n_soft_core"].tolist() == [1, 1, 1, 1, 0]
    assert st["n_dispensable"].tolist() == [1, 0, 1, 0, 1]
    assert st["n_private"].tolist() == [1, 1, 1, 0, 1]
    assert st["n_singletons"].tolist() == [0, 0, 1, 0, 1]


def test_sample_table_include_metadata(pm):
    st = pm.sample_table(include_metadata=True)
    assert list(st.columns[-2:]) == ["group", "depth"]
    assert st.loc["s3", "group"] == "B"
    assert st.loc["s3", "n_genes"] == 4
