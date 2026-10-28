import pytest

from pavvis import PavMatrix


PAV_CSV = """\
,sample_a,sample_b,sample_c,sample_d
gene_core,1,1,1,1
gene_ac,1,0,1,0
gene_b,0,1,0,0
gene_absent,0,0,0,0
"""

META_CSV = """\
,species,location,depth
sample_a,human,UK,10.5
sample_b,mouse,US,12.0
sample_c,human,DE,8.25
sample_d,mouse,UK,15.0
"""


@pytest.fixture
def pm(tmp_path):
    pav = tmp_path / "pav.csv"
    meta = tmp_path / "metadata.csv"
    pav.write_text(PAV_CSV)
    meta.write_text(META_CSV)
    return PavMatrix(pav, meta, soft_core_min=0.9, dispensable_min=0.3, private_min=0.0)


def test_subset_returns_new_pav_matrix(pm):
    sub = pm.subset(samples=["sample_a", "sample_b"])
    assert isinstance(sub, PavMatrix)
    assert sub is not pm


def test_subset_by_samples(pm):
    sub = pm.subset(samples=["sample_a", "sample_c"])
    assert sub.sample_ids == ["sample_a", "sample_c"]
    assert list(sub.metadata.index) == ["sample_a", "sample_c"]
    assert sub.gene_ids == pm.gene_ids


def test_subset_keeps_original_order(pm):
    sub = pm.subset(samples=["sample_c", "sample_a"], genes=["gene_b", "gene_core"])
    assert sub.sample_ids == ["sample_a", "sample_c"]
    assert sub.gene_ids == ["gene_core", "gene_b"]


def test_subset_by_genes(pm):
    sub = pm.subset(genes=["gene_core", "gene_ac"])
    assert sub.gene_ids == ["gene_core", "gene_ac"]
    assert sub.sample_ids == pm.sample_ids


def test_subset_by_query(pm):
    sub = pm.subset(query="species == 'human'")
    assert sub.sample_ids == ["sample_a", "sample_c"]


def test_subset_by_numeric_query(pm):
    sub = pm.subset(query="depth > 11")
    assert sub.sample_ids == ["sample_b", "sample_d"]


def test_subset_query_intersects_with_samples(pm):
    sub = pm.subset(samples=["sample_a", "sample_b", "sample_d"], query="location == 'UK'")
    assert sub.sample_ids == ["sample_a", "sample_d"]


def test_subset_samples_genes_and_query_combined(pm):
    sub = pm.subset(
        samples=["sample_a", "sample_b", "sample_c"],
        genes=["gene_ac", "gene_b"],
        query="species == 'human'",
    )
    assert sub.sample_ids == ["sample_a", "sample_c"]
    assert sub.gene_ids == ["gene_ac", "gene_b"]
    assert sub.pav.values.tolist() == [[1, 1], [0, 0]]


def test_subset_recomputes_gene_categories(pm):
    sub = pm.subset(samples=["sample_a", "sample_c"])
    assert sub.core_genes == ["gene_core", "gene_ac"]
    assert sub.absent_genes == ["gene_b", "gene_absent"]


def test_subset_drop_absent_removes_genes_absent_from_subset(pm):
    with pytest.warns(UserWarning, match="dropped 2 of 4 genes absent from all 2 retained samples"):
        sub = pm.subset(samples=["sample_a", "sample_c"], drop_absent=True)
    assert sub.gene_ids == ["gene_core", "gene_ac"]


def test_subset_drop_absent_does_not_warn_when_nothing_dropped(pm, recwarn):
    pm.subset(genes=["gene_core", "gene_ac"], drop_absent=True)
    assert not [w for w in recwarn if issubclass(w.category, UserWarning)]


def test_subset_keeps_absent_genes_by_default(pm):
    sub = pm.subset(samples=["sample_a", "sample_c"])
    assert "gene_absent" in sub.gene_ids
    assert "gene_b" in sub.gene_ids


def test_subset_preserves_thresholds(pm):
    sub = pm.subset(samples=["sample_a", "sample_b"])
    assert sub.thresholds == pm.thresholds


def test_subset_does_not_carry_caches(pm):
    pm.compute_jaccard()
    sub = pm.subset(samples=["sample_a", "sample_b"])
    assert sub._jaccard_matrix is None
    assert sub._umap_embedding is None


def test_subset_does_not_mutate_original(pm):
    sub = pm.subset(samples=["sample_a", "sample_b"], genes=["gene_core"])
    sub.pav.iloc[0, 0] = 0
    sub.metadata.iloc[0, 0] = "changed"
    assert pm.pav.loc["gene_core", "sample_a"] == 1
    assert pm.metadata.loc["sample_a", "species"] == "human"
    assert pm.sample_ids == ["sample_a", "sample_b", "sample_c", "sample_d"]
    assert len(pm.gene_ids) == 4


def test_subset_with_no_arguments_copies_everything(pm):
    sub = pm.subset()
    assert sub.sample_ids == pm.sample_ids
    assert sub.gene_ids == pm.gene_ids


def test_subset_string_samples_raises_type_error(pm):
    with pytest.raises(TypeError, match=r"samples=\['sample_a'\]"):
        pm.subset(samples="sample_a")


def test_subset_string_genes_raises_type_error(pm):
    with pytest.raises(TypeError, match=r"genes=\['gene_core'\]"):
        pm.subset(genes="gene_core")


def test_subset_unknown_sample_raises(pm):
    with pytest.raises(ValueError, match="sample_x"):
        pm.subset(samples=["sample_a", "sample_x"])


def test_subset_unknown_gene_raises(pm):
    with pytest.raises(ValueError, match="gene_x"):
        pm.subset(genes=["gene_core", "gene_x"])


def test_subset_empty_sample_result_raises(pm):
    with pytest.raises(ValueError, match="no samples"):
        pm.subset(query="species == 'cat'")


def test_subset_empty_gene_result_raises(pm):
    with pytest.raises(ValueError, match="no genes"):
        pm.subset(genes=["gene_absent"], drop_absent=True)
