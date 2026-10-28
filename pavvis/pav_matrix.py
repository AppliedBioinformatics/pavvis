from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
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

        self._umap_embedding: pd.DataFrame | None = None
        self._umap_params: tuple | None = None
        self._jaccard_matrix: pd.DataFrame | None = None

    @classmethod
    def _from_dataframes(
        cls,
        pav: pd.DataFrame,
        metadata: pd.DataFrame,
        soft_core_min: float = 0.95,
        dispensable_min: float = 0.15,
        private_min: float = 0.0,
    ) -> PavMatrix:
        """Construct a PavMatrix directly from DataFrames, skipping file I/O.

        Intended for internal use (e.g. creating subsets for demo plots). The
        caller is responsible for ensuring the DataFrames are already valid.
        """
        obj = cls.__new__(cls)
        obj.pav = pav
        obj.metadata = metadata
        obj.soft_core_min = soft_core_min
        obj.dispensable_min = dispensable_min
        obj.private_min = private_min
        obj._umap_embedding = None
        obj._umap_params = None
        obj._jaccard_matrix = None
        return obj

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

    def subset(
        self,
        samples: list[str] | None = None,
        genes: list[str] | None = None,
        query: str | None = None,
        drop_absent: bool = False,
    ) -> PavMatrix:
        """Return a new PavMatrix restricted to a subset of samples and/or genes.

        Filters are combined: a sample is kept only if it is in `samples` (when
        given) and matches `query` (when given). Samples and genes keep their order
        from the original matrix. Gene categories are recomputed on the subset, so a
        gene that is variable across the full matrix may be core or absent within
        the subset. Thresholds are carried over; cached results (UMAP, Jaccard) are
        not, and the original object is left unchanged.

        Args:
            samples: Sample IDs to keep. Defaults to None (all samples).
            genes: Gene IDs to keep. Defaults to None (all genes).
            query: A pandas query string evaluated against the metadata, e.g.
                ``"clade == 'Modern'"`` or ``"depth > 10"``. Column names are
                lowercase, as normalised on load. Defaults to None.
            drop_absent: If True, remove genes absent from every retained sample.
                A UserWarning reports how many genes were dropped. Defaults to False.

        Returns:
            A new PavMatrix containing copies of the selected data.

        Raises:
            TypeError: If `samples` or `genes` is passed as a single string
                rather than a list.
            ValueError: If any sample or gene ID is not found, or if the subset
                contains no samples or no genes.

        Warns:
            UserWarning: If drop_absent=True removes any genes, stating how many.
        """
        for name, ids in (("samples", samples), ("genes", genes)):
            if isinstance(ids, str):
                raise TypeError(
                    f"{name} must be a list of IDs, not a string. "
                    f"Use {name}=['{ids}'] to select a single ID."
                )

        sample_mask = np.ones(len(self.pav.columns), dtype=bool)
        gene_mask = np.ones(len(self.pav.index), dtype=bool)

        if samples is not None:
            unknown = [s for s in samples if s not in self.pav.columns]
            if unknown:
                raise ValueError(f"Samples not found in PAV matrix: {unknown}")
            sample_mask &= self.pav.columns.isin(samples)

        if query is not None:
            matched = self.metadata.query(query).index
            sample_mask &= self.pav.columns.isin(matched)

        if genes is not None:
            unknown = [g for g in genes if g not in self.pav.index]
            if unknown:
                raise ValueError(f"Genes not found in PAV matrix: {unknown}")
            gene_mask &= self.pav.index.isin(genes)

        if not sample_mask.any():
            raise ValueError("Subset contains no samples.")

        pav = self.pav.loc[gene_mask, sample_mask]
        n_selected = len(pav)
        if drop_absent:
            pav = pav.loc[pav.any(axis=1)]
        if pav.empty:
            raise ValueError("Subset contains no genes.")
        n_dropped = n_selected - len(pav)
        if n_dropped:
            warnings.warn(
                f"subset: dropped {n_dropped:,} of {n_selected:,} genes absent from all "
                f"{pav.shape[1]:,} retained samples (drop_absent=True).",
                UserWarning, stacklevel=2,
            )

        return PavMatrix._from_dataframes(
            pav=pav.copy(),
            metadata=self.metadata.loc[pav.columns].copy(),
            soft_core_min=self.soft_core_min,
            dispensable_min=self.dispensable_min,
            private_min=self.private_min,
        )

    def compute_umap(
        self,
        n_neighbors: int = 15,
        min_dist: float = 0.1,
        metric: str = "jaccard",
        random_state: int | None = None,
    ) -> PavMatrix:
        """Compute a 3D UMAP embedding of the samples and cache it on this object.

        Always computes 3 components so the result can be used for both 2D and 3D
        plots without recomputing. The embedding is computed on the transposed PAV
        matrix (samples × genes). Re-calling with the same parameters returns the
        cached result immediately. Re-calling with different parameters recomputes
        and overwrites the cache.

        Args:
            n_neighbors: UMAP n_neighbors parameter. Controls local vs global
                structure. Defaults to 15.
            min_dist: UMAP min_dist parameter. Controls how tightly points are
                packed. Defaults to 0.1.
            metric: Distance metric passed to UMAP. Defaults to 'jaccard', which
                ignores shared absences and is appropriate for binary PAV data.
            random_state: Random seed for reproducibility. Defaults to None.

        Returns:
            self, so calls can be chained: pm.compute_umap().plot...
        """
        import umap as umap_lib

        params = (n_neighbors, min_dist, metric, random_state)
        if self._umap_params == params and self._umap_embedding is not None:
            return self

        reducer = umap_lib.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            n_components=3,
            random_state=random_state,
        )
        coords = reducer.fit_transform(self.pav.T.values)
        self._umap_embedding = pd.DataFrame(
            coords, index=self.pav.columns, columns=["UMAP1", "UMAP2", "UMAP3"]
        )
        self._umap_params = params
        return self

    def compute_jaccard(self) -> PavMatrix:
        """Compute the pairwise Jaccard similarity matrix for all samples and cache it.

        Jaccard similarity between two samples is the number of genes present in both
        divided by the number of genes present in either, ignoring shared absences.
        The result is an N×N symmetric DataFrame (samples × samples) with values in
        [0, 1]. Re-calling returns the cached result immediately.

        Returns:
            self, so calls can be chained: pm.compute_jaccard().plot...
        """
        if self._jaccard_matrix is not None:
            return self

        pav = self.pav.values.astype(np.float64)  # (genes, samples)
        intersection = pav.T @ pav                  # (samples, samples)
        row_sums = pav.sum(axis=0)                  # (samples,)
        union = row_sums[:, None] + row_sums[None, :] - intersection
        with np.errstate(invalid="ignore"):
            similarity = np.where(union == 0, 0.0, intersection / union)
        self._jaccard_matrix = pd.DataFrame(
            similarity, index=self.pav.columns, columns=self.pav.columns,
        )
        return self

    def exclusive_gene_counts(self, column: str) -> dict[str, int]:
        """Count exclusive variable genes for each group in a discrete metadata column.

        An exclusive gene is a variable gene present in at least one sample within
        a group and absent from every sample outside that group.

        Args:
            column: A discrete metadata column name.

        Returns:
            Dict mapping each group label to its exclusive gene count, sorted
            alphabetically by group label.

        Raises:
            ValueError: If column is not found in metadata.
        """
        groups, group_matrix = self._group_matrix(column)
        pav_var = self.pav.loc[self.variable_genes].values  # (genes, samples)
        group_counts = pav_var @ group_matrix               # (genes, groups)
        # Exclusive: every presence of the gene falls inside the group.
        exclusive = (group_counts > 0) & (group_counts == pav_var.sum(axis=1)[:, None])
        return {g: int(n) for g, n in zip(groups, exclusive.sum(axis=0))}

    def _group_matrix(self, column: str) -> tuple[list[str], np.ndarray]:
        """One-hot encode a discrete metadata column against the PAV sample order.

        Returns:
            Tuple of (groups, matrix), where groups is the sorted list of group
            labels and matrix is an int array of shape (samples, groups) with a 1
            where the sample belongs to the group.
        """
        from pavvis._validation import validate_color_by
        validate_color_by(column, self.metadata, expected_type="discrete")

        col = self.metadata.loc[self.pav.columns, column].astype(str)
        groups = sorted(col.unique())
        labels = col.to_numpy(dtype=object)  # plain array; pandas 3 string arrays don't broadcast
        matrix = (labels[:, None] == np.array(groups, dtype=object)[None, :]).astype(int)
        return groups, matrix

    def _gene_categories(self) -> pd.Series:
        """Category label for every gene, using the same edges as the category properties.

        Genes with a presence frequency in (0, private_min] belong to no category and
        are labelled 'unclassified'.
        """
        freq = self.presence_frequency.values
        labels = np.select(
            [
                freq == 1.0,
                freq == 0.0,
                freq >= self.soft_core_min,
                freq >= self.dispensable_min,
                freq > self.private_min,
            ],
            ["core", "absent", "soft_core", "dispensable", "private"],
            default="unclassified",
        )
        return pd.Series(labels, index=self.pav.index, name="category")

    def gene_table(self, column: str | None = None) -> pd.DataFrame:
        """Summarise every gene's presence and category as a DataFrame.

        The table always has these columns:

        - `n_present`: number of samples the gene is present in.
        - `frequency`: fraction of samples the gene is present in.
        - `category`: 'core', 'soft_core', 'dispensable', 'private', 'absent', or
          'unclassified' for genes with a frequency in (0, private_min].

        When `column` is given, it also has:

        - `n_<group>`: number of samples in each group the gene is present in.
        - `exclusive_to`: the group a variable gene is present in exclusively
          (see `exclusive_gene_counts`), or NaN if it is not exclusive to any group.

        Args:
            column: An optional discrete metadata column used to add the per-group
                columns. Defaults to None.

        Returns:
            DataFrame indexed by gene ID.

        Raises:
            ValueError: If column is not a discrete metadata column.
        """
        n_present = self.pav.sum(axis=1)
        table = pd.DataFrame({
            "n_present": n_present,
            "frequency": n_present / len(self.sample_ids),
            "category": self._gene_categories(),
        })
        table.index.name = "gene_id"

        if column is not None:
            groups, group_matrix = self._group_matrix(column)
            group_counts = self.pav.values @ group_matrix  # (genes, groups)
            for i, g in enumerate(groups):
                table[f"n_{g}"] = group_counts[:, i]

            n_groups_present = (group_counts > 0).sum(axis=1)
            is_variable = (n_present > 0) & (n_present < len(self.sample_ids))
            exclusive = is_variable.values & (n_groups_present == 1)
            table["exclusive_to"] = pd.Series(
                np.array(groups, dtype=object)[group_counts.argmax(axis=1)],
                index=table.index,
            ).where(exclusive)

        return table

    def sample_table(self, include_metadata: bool = False) -> pd.DataFrame:
        """Summarise every sample's gene content as a DataFrame.

        The table has these columns:

        - `n_genes`: number of genes present in the sample.
        - `n_core`, `n_soft_core`, `n_dispensable`, `n_private`: number of genes
          present in the sample from each category.
        - `n_singletons`: number of genes present in this sample and no other.

        Args:
            include_metadata: If True, append the sample's metadata columns.
                Defaults to False.

        Returns:
            DataFrame indexed by sample ID.
        """
        categories = self._gene_categories()
        n_present = self.pav.sum(axis=1)
        table = pd.DataFrame({"n_genes": self.pav.sum(axis=0)})
        for category in ("core", "soft_core", "dispensable", "private"):
            table[f"n_{category}"] = self.pav.loc[categories == category].sum(axis=0)
        table["n_singletons"] = self.pav.loc[n_present == 1].sum(axis=0)
        table.index.name = "sample"

        if include_metadata:
            table = table.join(self.metadata)
        return table

    @property
    def presence_frequency(self) -> pd.Series:
        """Fraction of samples in which each gene is present.

        Returns:
            pd.Series indexed by gene ID, with values in [0.0, 1.0].
        """
        return self.pav.sum(axis=1) / len(self.sample_ids)
