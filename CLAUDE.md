# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run all tests
pytest

# Run a single test file
pytest tests/test_pav_matrix.py

# Run a single test by name
pytest tests/test_pav_matrix.py::test_core_genes_returns_only_genes_present_in_all_samples

# Serve documentation locally
mkdocs serve

# Regenerate HTML plot fragments for docs (run before mkdocs build when plots change)
python gen_plots.py
python gen_plots.py --reset   # force-regenerate all

# Install dev dependencies
pip install -e ".[dev]"

# Install docs dependencies
pip install -e ".[docs]"
```

## Architecture

**Pavvis** is a Python library for Presence/Absence Variation (PAV) pangenome analysis. The entry point is `PavMatrix` — all plot functions receive a `PavMatrix` instance.

### Core

- `pavvis/pav_matrix.py` — `PavMatrix` class. Loads two CSVs (genes×samples PAV matrix + sample metadata), validates them, classifies genes into categories (core / soft_core / dispensable / private / absent) based on configurable frequency thresholds, and caches a UMAP embedding. Metadata column names are normalised to lowercase on load; numeric-looking columns are cast to float64.
- `pavvis/_validation.py` — standalone validation helpers called during `__init__`. Also contains `infer_column_type` (numeric + >10 unique values → continuous, else discrete) and `validate_color_by`, which are used by plot functions to select colour scales.

### Plots (`pavvis/plots/`)

Each module is a collection of functions that accept a `PavMatrix` and return a Plotly `Figure`. The default Plotly template is set to `"simple_white"` in `_helpers.py` at import time.

- `presence.py` — `frequency_histogram`, `per_sample_boxplot`, `per_sample_scatter`
- `curves.py` — `legacy_pangenome_curve`, `variable_gene_curve`, `grouped_variable_curve`, `jaccard_similarity_curve`
- `variable_genes.py` — `frequency_histogram`, `exclusive_genes_bar`
- `umap.py` — `scatter`, `scatter3d` (require `pm.compute_umap()` to have been called first, or call it in the chain)
- `_helpers.py` — `_resolve_gene_index` (maps a `gene_set` string to a gene list) and `_gene_set_labels` (axis labels)

### Visualisation (`pavvis/visualisation/`)

Currently a stub (`_base.py` is empty). Reserved for shared rendering utilities.

### Documentation

MkDocs + Material theme with `mkdocstrings` (Google docstring style). Plot pages embed interactive Plotly HTML fragments from `docs/statics/plots/` via `pymdownx.snippets`. Run `python gen_plots.py` before building or serving docs when plot output is stale.

### Test data

`tests/data/` contains `pav.csv` and `metadata.csv` used by `main.py` (manual exploration) and `gen_plots.py`. Tests in `tests/` create their own minimal CSV fixtures via `tmp_path`.