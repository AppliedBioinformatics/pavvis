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

META_CSV = """\
,species,location
sample_a,human,UK
sample_b,mouse,US
sample_c,human,DE
"""


def test_init_accepts_path_objects(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(pav_path, meta_path)
    assert pm.pav is not None
    assert pm.metadata is not None


def test_init_accepts_strings(tmp_path):
    pav_path = _write_tmp(tmp_path, "pav.csv", PAV_CSV)
    meta_path = _write_tmp(tmp_path, "meta.csv", META_CSV)
    pm = PavMatrix(str(pav_path), str(meta_path))
    assert pm.pav is not None