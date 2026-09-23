import pytest
import pandas as pd
from pavvis._validation import validate_sample_alignment, validate_pav_values


def _pav(data: dict, index: list[str]) -> pd.DataFrame:
    return pd.DataFrame(data, index=index)


# --- validate_sample_alignment ---

def test_alignment_passes_when_samples_match():
    pav = _pav({"s1": [1, 0], "s2": [0, 1]}, ["g1", "g2"])
    meta = pd.DataFrame({"col": [1, 2]}, index=["s1", "s2"])
    validate_sample_alignment(pav, meta)  # must not raise


def test_alignment_raises_when_pav_has_extra_sample():
    pav = _pav({"s1": [1], "s2": [0], "s3": [1]}, ["g1"])
    meta = pd.DataFrame({"col": [1, 2]}, index=["s1", "s2"])
    with pytest.raises(ValueError) as exc_info:
        validate_sample_alignment(pav, meta)
    assert "s3" in str(exc_info.value)


def test_alignment_raises_when_metadata_has_extra_sample():
    pav = _pav({"s1": [1], "s2": [0]}, ["g1"])
    meta = pd.DataFrame({"col": [1, 2, 3]}, index=["s1", "s2", "s3"])
    with pytest.raises(ValueError) as exc_info:
        validate_sample_alignment(pav, meta)
    assert "s3" in str(exc_info.value)


def test_alignment_raises_when_both_sides_differ():
    pav = _pav({"s1": [1], "s_pav_only": [0]}, ["g1"])
    meta = pd.DataFrame({"col": [1, 2]}, index=["s1", "s_meta_only"])
    with pytest.raises(ValueError) as exc_info:
        validate_sample_alignment(pav, meta)
    msg = str(exc_info.value)
    assert "s_pav_only" in msg
    assert "s_meta_only" in msg


# --- validate_pav_values ---

def test_pav_values_passes_for_binary_matrix():
    pav = _pav({"s1": [1, 0], "s2": [0, 1]}, ["g1", "g2"])
    validate_pav_values(pav)  # must not raise


def test_pav_values_raises_on_blank_cell():
    pav = _pav({"s1": [1, None], "s2": [0, 1]}, ["g1", "g2"])
    with pytest.raises(ValueError, match="blank"):
        validate_pav_values(pav)


def test_pav_values_raises_on_non_binary_value():
    pav = _pav({"s1": [1, 2], "s2": [0, 1]}, ["g1", "g2"])
    with pytest.raises(ValueError, match="0 or 1"):
        validate_pav_values(pav)


def test_pav_values_blank_reported_before_non_binary():
    """Blank check runs first so the error message is about blanks, not 0/1."""
    pav = _pav({"s1": [None, 2], "s2": [0, 1]}, ["g1", "g2"])
    with pytest.raises(ValueError, match="blank"):
        validate_pav_values(pav)