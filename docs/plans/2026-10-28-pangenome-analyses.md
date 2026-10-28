# Pangenome analyses 3–6: implementation plan

## Context
Pavvis currently only plots descriptive views of a PAV matrix. We are adding four analyses
that are standard in pangenome papers:

3. Per-group pangenome curves
4. Subsetting
5. Gene and sample output tables
6. PCoA

The numbering follows the original feature list. Items 1 (Heaps' law openness fit) and
2 (genomic fluidity) were dropped. They are built **one at a time, in this order**. You
review each one before the next starts. I will not run any git operations; you handle
commits.

## Shared groundwork (do with item 3)
- Add `scipy` to `pyproject.toml` dependencies. `heatmap.py` already imports it but it is
  not declared.
- Create a new `pavvis/_accumulation.py` containing
  `accumulation_runs(pav: np.ndarray, permutations: int, rng) -> (pan_runs, core_runs)`.
  It is the permutation loop currently duplicated in `curves.legacy_pangenome_curve`,
  `curves.variable_gene_curve` and `curves.grouped_variable_curve`. Refactor the first two
  to call it; their behaviour must not change.
- Add an optional `random_state: int | None = None` argument to every new function that
  permutes samples (passed to `np.random.default_rng`), so tests are deterministic.
- Move `_STRIP_COLOURS` from `plots/heatmap.py` into `plots/_helpers.py` as a shared
  palette, and import it back into `heatmap.py`.

---

## 3. Per-group pangenome curves
**`curves.grouped_pangenome_curve(pm, column, groups=None, curve="pangenome", permutations=100, log_y=False, random_state=None)`**
- `curve` is one of `"pangenome"`, `"core"` or `"variable"`.
- Run an independent `accumulation_runs` on each group's samples and overlay the curves
  (x = 1..n_g), each with a mean line and ±1 std band. Each legend entry toggles its line and
  band together, using `legendgroup` as the existing curves do.
- Legend labels read `"<group> (n=…)"`. Colours come from the shared palette.
- `groups` restricts which groups are drawn, and in what order; the default is all of them, sorted.
- Groups with fewer than 2 samples are skipped with a `UserWarning`.
- Validate the column as discrete and check that every value in `groups` exists, with the
  same error wording as `grouped_variable_curve`.

**Tests:** the figure has one mean trace per group, each group's x length equals its sample
count, unknown groups raise, and a single-sample group triggers the warning.

## 4. `PavMatrix.subset`
**`subset(samples=None, genes=None, query=None, drop_absent=False) -> PavMatrix`**
- `query` is passed to `self.metadata.query(query)` and intersected with `samples` if both
  are given.
- Unknown sample or gene IDs raise `ValueError`. An empty result also raises.
- `drop_absent=True` removes genes absent from every retained sample.
- Builds the new object with the existing `_from_dataframes`, using `.copy()` of the
  DataFrames and the same thresholds. Caches are not carried over.
- Also add the new cache attributes from item 6 to `__init__` and `_from_dataframes`.

**Tests:** filter by samples, genes and query; check combined filters, errors and
`drop_absent`; check thresholds are preserved; check the original object is not mutated.

## 5. Output tables
- Private `_gene_categories() -> pd.Series`, vectorised with `np.select`. The category edges
  must match the existing properties exactly.
  - If `private_min > 0`, genes with frequency in (0, private_min] are labelled `"unclassified"`.
    They belong to no existing category.
- Private `_group_matrix(column) -> (groups, np.ndarray samples×groups one-hot)`. Refactor
  `exclusive_gene_counts` to use it; its output must not change.

**`gene_table(column=None) -> pd.DataFrame`** (index `gene_id`)
- Columns: `n_present`, `frequency`, `category`.
- With `column`, it also has `n_<group>` (presence count per group, computed as
  `pav @ group_matrix`) and `exclusive_to`. `exclusive_to` is the group name if the gene is
  variable and present only in that group, otherwise NaN.
- Its exclusive counts must equal `exclusive_gene_counts`.

**`sample_table(include_metadata=False) -> pd.DataFrame`** (index sample)
- Columns: `n_genes`, `n_core`, `n_soft_core`, `n_dispensable`, `n_private`, and
  `n_singletons` (genes present only in that sample).
- Optionally join the metadata columns.

**Tests:** hand-built matrix values; category counts equal `gene_counts`;
`exclusive_to` counts equal `exclusive_gene_counts`; singletons are correct.

## 6. PCoA (classical MDS on Jaccard distance)
**`PavMatrix.compute_pcoa(n_components=10) -> PavMatrix`** (returns self for chaining)
- Calls `compute_jaccard()` if needed. D = 1 − J.
- Computes B = −½·C·D²·C (C is the centring matrix), then `np.linalg.eigh`, sorted in
  descending order.
- Keeps only positive eigenvalues, because Jaccard distance is not Euclidean. The docstring
  should say so.
- Coordinates are V·√λ. Explained variance is λ / Σλ₊.
- Caches `_pcoa_embedding` (columns `PCo1…`) and `_pcoa_explained`. Rerunning with the same
  `n_components` is a no-op.

**New module `pavvis/plots/pcoa.py`**
- `scatter(pm, color_by=None, components=(1, 2))`. Axis titles read `"PCo1 (12.3%)"`.
- `scree(pm)`: a bar chart of % explained per axis.
- Pull the discrete/continuous trace-building out of `umap.scatter` into
  `_helpers._embedding_scatter_traces(emb, xcol, ycol, metadata, color_by)`. Use it in both
  `umap.scatter` and `pcoa.scatter`; UMAP output must not change.
- Add `_require_pcoa`, following `_require_embedding`.

**Tests:**
- For Euclidean points (e.g. distances from a 2D configuration fed in through a mocked
  Jaccard matrix), the pairwise distances of the coordinates reproduce D.
- Explained variance lies in [0, 1] and is non-increasing.
- Two clearly separated sample groups split on PCo1.
- Plotting before computing raises `RuntimeError`.

---

## Per-item checklist (each step)
1. Write the tests first, using minimal `tmp_path` CSV fixtures in the style of
   `tests/test_pav_matrix.py`. New plot tests go in `tests/test_plots_<module>.py`.
2. Implement it, using Google-style docstrings that match the existing ones.
3. Docs:
   - Add the function to the relevant `docs/api/plots/*.md` page.
   - For item 6, add `api/plots/pcoa.md` to the `mkdocs.yml` nav.
   - Add a demo plot to `gen_plots.py` (using `tests/data`, with group column `clade`).
4. Run `pytest`, then `python gen_plots.py`, and stop for your review.

## Verification
- `pytest` passes, including all existing tests (this covers the refactors of the curves,
  `exclusive_gene_counts` and `umap.scatter`).
- `python gen_plots.py` produces the new HTML fragments from the about 1,000-sample wheat
  data. Sanity check: PCo1/PCo2 separate the `Modern` and AG clades.
- `mkdocs serve` renders the new API sections.
