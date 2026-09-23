from __future__ import annotations

import difflib
from typing import Literal

import pandas as pd


def validate_sample_alignment(pav: pd.DataFrame, metadata: pd.DataFrame) -> None:
    pav_samples = set(pav.columns)
    meta_samples = set(metadata.index)

    in_pav_not_meta = sorted(pav_samples - meta_samples)
    in_meta_not_pav = sorted(meta_samples - pav_samples)

    if in_pav_not_meta or in_meta_not_pav:
        raise ValueError(
            f"Sample ID mismatch between PAV matrix and metadata.\n"
            f"  In PAV but not metadata ({len(in_pav_not_meta)}): {in_pav_not_meta}\n"
            f"  In metadata but not PAV ({len(in_meta_not_pav)}): {in_meta_not_pav}"
        )


def validate_pav_values(pav: pd.DataFrame) -> None:
    blank_mask = pav.isna()
    if blank_mask.any().any():
        blank_coords = [
            (gene, sample)
            for gene in pav.index
            for sample in pav.columns
            if blank_mask.loc[gene, sample]
        ]
        raise ValueError(
            f"PAV matrix contains blank cells at {len(blank_coords)} location(s): "
            f"{blank_coords[:10]}{'...' if len(blank_coords) > 10 else ''}"
        )

    invalid_mask = ~pav.isin([0, 1])
    if invalid_mask.any().any():
        invalid_coords = [
            (gene, sample, pav.loc[gene, sample])
            for gene in pav.index
            for sample in pav.columns
            if invalid_mask.loc[gene, sample]
        ]
        raise ValueError(
            f"PAV matrix values must be 0 or 1. Found invalid values at "
            f"{len(invalid_coords)} location(s): "
            f"{invalid_coords[:10]}{'...' if len(invalid_coords) > 10 else ''}"
        )


def validate_color_by(
    color_by: str,
    metadata: pd.DataFrame,
    expected_type: Literal["discrete", "continuous"] | None = None,
) -> None:
    """Raise ValueError if color_by is not a valid metadata column, or is the wrong type.

    Args:
        color_by: The metadata column name to validate.
        metadata: The metadata DataFrame.
        expected_type: If provided, raises if the column does not match
            'discrete' or 'continuous' as determined by infer_column_type.
    """
    if color_by not in metadata.columns:
        close = difflib.get_close_matches(color_by.lower(), metadata.columns.str.lower(), n=1)
        suggestion = f" Did you mean '{close[0]}'?" if close else ""
        raise ValueError(
            f"'{color_by}' is not a metadata column.{suggestion} "
            f"Available columns: {list(metadata.columns)}"
        )
    if expected_type is not None:
        actual_type = infer_column_type(metadata[color_by])
        if actual_type != expected_type:
            raise ValueError(
                f"'{color_by}' is a {actual_type} column but this plot requires "
                f"a {expected_type} column."
            )


def infer_column_type(series: pd.Series) -> Literal["discrete", "continuous"]:
    """Infer whether a metadata column should be treated as discrete or continuous.

    Numeric columns with more than 10 unique values are considered continuous;
    all others (including numeric columns with few unique values) are discrete.

    Args:
        series: A metadata column as a pandas Series.

    Returns:
        'continuous' or 'discrete'.
    """
    if pd.api.types.is_numeric_dtype(series) and series.nunique() > 10:
        return "continuous"
    return "discrete"