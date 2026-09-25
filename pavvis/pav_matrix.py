from __future__ import annotations

from pathlib import Path

import pandas as pd

from pavvis._validation import validate_sample_alignment, validate_pav_values, validate_subcategory_thresholds


class PavMatrix:
    def __init__(
        self,
        matrix_path: str | Path,
        metadata_path: str | Path,
        soft_core_min: float = 0.95,
        dispensable_min: float = 0.15,
        private_min: float = 0.0,
    ) -> None:
        """Load and validate a PAV matrix with associated sample metadata.

        Args:
            matrix_path: Path to the PAV matrix CSV. Genes must be rows (unique index),
                samples must be columns (unique headers), values must be 0 or 1.
            metadata_path: Path to the metadata CSV. Samples must be rows (index must
                match PAV matrix column headers exactly).
            soft_core_min: Minimum presence frequency for soft core genes (inclusive).
                Defaults to 0.95.
            dispensable_min: Minimum presence frequency for dispensable genes (inclusive).
                Defaults to 0.15.
            private_min: Minimum presence frequency for private genes (exclusive lower
                bound — genes at exactly this frequency are excluded). Defaults to 0.0.

        Raises:
            FileNotFoundError: If either file path does not exist.
            ValueError: If sample IDs do not match between files, if any PAV value
                is not 0 or 1, if any cell is blank, or if thresholds are invalid.
        """
        matrix_path = Path(matrix_path)
        metadata_path = Path(metadata_path)

        # Check paths exist.
        if not matrix_path.exists():
            raise FileNotFoundError(f"PAV matrix file not found: {matrix_path}")
        if not metadata_path.exists():
            raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

        # Load data from files.
        self.pav: pd.DataFrame = pd.read_csv(matrix_path, index_col=0)
        self.metadata: pd.DataFrame = pd.read_csv(metadata_path, index_col=0)

        # Normalise metadata column names to lowercase with no surrounding whitespace.
        self.metadata.columns = self.metadata.columns.str.strip().str.lower()

        # Cast numeric-looking columns to float64 so dtype-based inference is reliable.
        def _try_numeric(col: pd.Series) -> pd.Series:
            converted = pd.to_numeric(col, errors="coerce")
            return converted if converted.notna().all() else col

        self.metadata = self.metadata.apply(_try_numeric)

        validate_subcategory_thresholds(soft_core_min, dispensable_min, private_min)
        self.soft_core_min = soft_core_min
        self.dispensable_min = dispensable_min
        self.private_min = private_min

        # Checks to affirm files are of the correct structure. (Can add other checks here).
        validate_sample_alignment(self.pav, self.metadata)
        validate_pav_values(self.pav)

    # Built ins.
    def __repr__(self) -> str:
        """Return a summary string of the form 'PavMatrix(N genes × M samples)'."""
        return f"PavMatrix({len(self.gene_ids)} genes × {len(self.sample_ids)} samples)"

    # Properties.
    @property
    def sample_ids(self) -> list[str]:
        """Sample identifiers, taken from the PAV matrix column headers.

        Returns:
            List of sample ID strings.
        """
        return list(self.pav.columns)

    @property
    def gene_ids(self) -> list[str]:
        """Gene identifiers, taken from the PAV matrix row index.

        Returns:
            List of gene ID strings.
        """
        return list(self.pav.index)

    @property
    def metadata_columns(self) -> list[str]:
        """Column names present in the loaded metadata DataFrame.

        Returns:
            List of metadata column name strings.
        """
        return list(self.metadata.columns)

    @property
    def core_genes(self) -> list[str]:
        """Genes present in every sample (all values equal 1).

        Returns:
            List of gene ID strings.
        """
        return list(self.pav.index[self.pav.eq(1).all(axis=1)])

    @property
    def variable_genes(self) -> list[str]:
        """Genes present in at least one but not all samples.

        Returns:
            List of gene ID strings.
        """
        row_sums = self.pav.sum(axis=1)
        mask = (row_sums > 0) & (row_sums < len(self.pav.columns))
        return list(self.pav.index[mask])

    @property
    def absent_genes(self) -> list[str]:
        """Genes absent from every sample (all values equal 0).

        Returns:
            List of gene ID strings.
        """
        return list(self.pav.index[self.pav.eq(0).all(axis=1)])

    @property
    def soft_core_genes(self) -> list[str]:
        """Variable genes present in >= soft_core_min and < 1.0 of samples.

        Returns:
            List of gene ID strings.
        """
        freq = self.presence_frequency
        mask = (freq >= self.soft_core_min) & (freq < 1.0)
        return list(self.pav.index[mask])

    @property
    def dispensable_genes(self) -> list[str]:
        """Variable genes present in >= dispensable_min and < soft_core_min of samples.

        Returns:
            List of gene ID strings.
        """
        freq = self.presence_frequency
        mask = (freq >= self.dispensable_min) & (freq < self.soft_core_min)
        return list(self.pav.index[mask])

    @property
    def private_genes(self) -> list[str]:
        """Variable genes present in > private_min and < dispensable_min of samples.

        Returns:
            List of gene ID strings.
        """
        freq = self.presence_frequency
        mask = (freq > self.private_min) & (freq < self.dispensable_min)
        return list(self.pav.index[mask])

    @property
    def thresholds(self) -> dict[str, float]:
        """Current variable gene subcategory thresholds.

        Returns:
            Dict with keys 'soft_core_min', 'dispensable_min', 'private_min'.
        """
        return {
            "soft_core_min": self.soft_core_min,
            "dispensable_min": self.dispensable_min,
            "private_min": self.private_min,
        }

    @property
    def gene_counts(self) -> dict[str, int]:
        """Count of genes in each category and subcategory.

        Returns:
            Dict with keys 'core', 'soft_core', 'dispensable', 'private', and 'absent'
            mapping to integer counts.
        """
        return {
            "core": len(self.core_genes),
            "soft_core": len(self.soft_core_genes),
            "dispensable": len(self.dispensable_genes),
            "private": len(self.private_genes),
            "absent": len(self.absent_genes),
        }

    @property
    def presence_frequency(self) -> pd.Series:
        """Fraction of samples in which each gene is present.

        Returns:
            pd.Series indexed by gene ID, with values in [0.0, 1.0].
        """
        return self.pav.sum(axis=1) / len(self.sample_ids)
