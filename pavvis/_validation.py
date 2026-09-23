from __future__ import annotations

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