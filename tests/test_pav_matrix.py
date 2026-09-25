import pytest
from pathlib import Path
from pavvis import PavMatrix


def _write_tmp(tmp_path, filename, content):
    p = tmp_path / filename
    p.write_text(content)
    return p


PAV_CSV = """\
,sample_a,sample_b,sample_c
gene_1,1,0,1
gene_2,0,1,1
"""

PAV_CORE_CSV = """\
,sample_a,sample_b,sample_c
gene_core,1,1,1
gene_partial,1,0,1
gene_absent,0,0,0
"""

META_CSV = """\
,species,location
sample_a,human,UK
sample_b,mouse,US
sample_c,human,DE
"""

META_MISMATCH_CSV = """\
,species,location
sample_a,human,UK
sample_x,mouse,US
sample_c,human,DE
"""

PAV_WITH_NAN_CSV = """\
,sample_a,sample_b,sample_c
gene_1,1,,1
gene_2,0,1,1
"""

PAV_WITH_INVALID_CSV = """\
,sample_a,sample_b,sample_c
gene_1,1,2,1
gene_2,0,1,1
"""


def test_init_accepts_path_objects(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.pav.shape == (2, 3)
    assert list(pm.pav.columns) == ["sample_a", "sample_b", "sample_c"]
    assert list(pm.pav.index) == ["gene_1", "gene_2"]
    assert pm.metadata.shape == (3, 2)
    assert list(pm.metadata.index) == ["sample_a", "sample_b", "sample_c"]


def test_init_accepts_strings(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(str(pav_path), str(meta_path))
    assert pm.pav.shape == (2, 3)
    assert pm.metadata.shape == (3, 2)


def test_missing_pav_raises_file_not_found(tmp_path):
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    with pytest.raises(FileNotFoundError, match="pav"):
        PavMatrix(tmp_path / "nonexistent.csv", meta_path)


def test_missing_metadata_raises_file_not_found(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CSV)
    with pytest.raises(FileNotFoundError, match="metadata"):
        PavMatrix(pav_path, tmp_path / "nonexistent.csv")


def test_mismatched_sample_ids_raises_value_error(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_MISMATCH_CSV)
    with pytest.raises(ValueError) as exc_info:
        PavMatrix(pav_path, meta_path)
    msg = str(exc_info.value)
    assert "sample_b" in msg   # in PAV but not metadata
    assert "sample_x" in msg   # in metadata but not PAV


def test_matching_sample_ids_loads_successfully(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert set(pm.pav.columns) == set(pm.metadata.index)


def test_blank_cells_raise_value_error(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_WITH_NAN_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    with pytest.raises(ValueError, match="blank"):
        PavMatrix(pav_path, meta_path)


def test_non_binary_values_raise_value_error(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_WITH_INVALID_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    with pytest.raises(ValueError, match="0 or 1"):
        PavMatrix(pav_path, meta_path)


def test_sample_ids_returns_list_of_strings(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.sample_ids == ["sample_a", "sample_b", "sample_c"]
    assert all(isinstance(s, str) for s in pm.sample_ids)


def test_gene_ids_returns_list_of_strings(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.gene_ids == ["gene_1", "gene_2"]
    assert all(isinstance(s, str) for s in pm.gene_ids)


def test_core_genes_returns_only_genes_present_in_all_samples(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CORE_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.core_genes == ["gene_core"]


def test_core_genes_returns_empty_when_no_gene_is_universal(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.core_genes == []


def test_variable_genes_returns_genes_present_in_some_but_not_all_samples(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CORE_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.variable_genes == ["gene_partial"]


def test_variable_genes_excludes_fully_absent_genes(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CORE_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert "gene_absent" not in pm.variable_genes


def test_variable_genes_excludes_core_genes(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CORE_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert "gene_core" not in pm.variable_genes


def test_absent_genes_returns_genes_with_no_presence_across_samples(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CORE_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.absent_genes == ["gene_absent"]


def test_absent_genes_excludes_core_and_variable_genes(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CORE_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert "gene_core" not in pm.absent_genes
    assert "gene_partial" not in pm.absent_genes


def test_absent_genes_returns_empty_when_all_genes_present_in_at_least_one_sample(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.absent_genes == []


def test_variable_genes_returns_empty_when_all_genes_are_core(tmp_path):
    pav_all_core = """\
,sample_a,sample_b,sample_c
gene_1,1,1,1
gene_2,1,1,1
"""
    pav_path = _write_tmp(tmp_path, "pav.csv", pav_all_core)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.variable_genes == []


# --- Variable gene subcategory fixtures ---
# 4 samples → presence frequencies: 0, 0.25, 0.5, 0.75, 1.0
PAV_FREQ_CSV = """\
,sample_a,sample_b,sample_c,sample_d
gene_absent,0,0,0,0
gene_private,1,0,0,0
gene_dispensable,1,1,0,0
gene_soft_core,1,1,1,0
gene_core,1,1,1,1
"""

META_4_CSV = """\
,species
sample_a,human
sample_b,mouse
sample_c,human
sample_d,mouse
"""


def test_default_thresholds(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_FREQ_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_4_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.thresholds == {"soft_core_min": 0.95, "dispensable_min": 0.15, "private_min": 0.0}


def test_custom_thresholds_stored(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_FREQ_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_4_CSV)
    pm = PavMatrix(pav_path, meta_path, soft_core_min=0.8, dispensable_min=0.3, private_min=0.1)
    assert pm.soft_core_min == 0.8
    assert pm.dispensable_min == 0.3
    assert pm.private_min == 0.1


def test_soft_core_genes(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_FREQ_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_4_CSV)
    # With default thresholds: soft_core = freq >= 0.95 and < 1.0
    # gene_soft_core has freq 0.75 → not soft core with defaults
    # Use custom thresholds to put gene_soft_core (freq=0.75) in soft core
    pm = PavMatrix(pav_path, meta_path, soft_core_min=0.7, dispensable_min=0.3, private_min=0.0)
    assert pm.soft_core_genes == ["gene_soft_core"]


def test_dispensable_genes(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_FREQ_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_4_CSV)
    # gene_dispensable has freq 0.5, gene_soft_core has freq 0.75
    # With thresholds 0.7 / 0.3 / 0.0: dispensable = freq >= 0.3 and < 0.7
    pm = PavMatrix(pav_path, meta_path, soft_core_min=0.7, dispensable_min=0.3, private_min=0.0)
    assert pm.dispensable_genes == ["gene_dispensable"]


def test_private_genes(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_FREQ_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_4_CSV)
    # gene_private has freq 0.25; with thresholds 0.7 / 0.3 / 0.0: private = freq > 0.0 and < 0.3
    pm = PavMatrix(pav_path, meta_path, soft_core_min=0.7, dispensable_min=0.3, private_min=0.0)
    assert pm.private_genes == ["gene_private"]


def test_subcategory_genes_exclude_core_and_absent(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_FREQ_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_4_CSV)
    pm = PavMatrix(pav_path, meta_path, soft_core_min=0.7, dispensable_min=0.3, private_min=0.0)
    for gene_list in [pm.soft_core_genes, pm.dispensable_genes, pm.private_genes]:
        assert "gene_core" not in gene_list
        assert "gene_absent" not in gene_list


def test_subcategory_genes_return_empty_when_none_in_range(tmp_path):
    pav_all_core = """\
,sample_a,sample_b,sample_c,sample_d
gene_1,1,1,1,1
gene_2,1,1,1,1
"""
    pav_path = _write_tmp(tmp_path, "pav.csv", pav_all_core)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_4_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.soft_core_genes == []
    assert pm.dispensable_genes == []
    assert pm.private_genes == []


def test_gene_counts_includes_subcategories(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_FREQ_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_4_CSV)
    pm = PavMatrix(pav_path, meta_path, soft_core_min=0.7, dispensable_min=0.3, private_min=0.0)
    counts = pm.gene_counts
    assert set(counts.keys()) == {"core", "soft_core", "dispensable", "private", "absent"}
    assert counts["core"] == 1
    assert counts["soft_core"] == 1
    assert counts["dispensable"] == 1
    assert counts["private"] == 1
    assert counts["absent"] == 1


def test_invalid_threshold_out_of_range_raises(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_FREQ_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_4_CSV)
    with pytest.raises(ValueError, match="soft_core_min"):
        PavMatrix(pav_path, meta_path, soft_core_min=1.0)


def test_invalid_threshold_ordering_raises(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_FREQ_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_4_CSV)
    with pytest.raises(ValueError, match="private_min < dispensable_min < soft_core_min"):
        PavMatrix(pav_path, meta_path, soft_core_min=0.5, dispensable_min=0.8, private_min=0.0)