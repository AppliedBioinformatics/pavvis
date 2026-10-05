<p align="center">
  <img src="docs/statics/banner.jpg" alt="pavvis — Python Gene PAV Matrix Analysis & Visualization" width="800"/>
</p>

**pavvis** is a Python package for analysing and visualising gene Presence/Absence Variation (PAV) data from pangenome studies. It is designed to be easy to use and optimised for large-scale genomic datasets.

## Installation

```bash
pip install pavvis
```

For a development install from source:

```bash
git clone https://github.com/<your-username>/Pavvis.git
cd Pavvis
pip install -e ".[dev]"
```

## Quick start

```python
from pavvis import PavMatrix

pm = PavMatrix(
    matrix_path="pav_matrix.csv",
    metadata_path="metadata.csv",
)
print(pm)
# PavMatrix(34295 genes × 1039 samples)
```

## Gene categories

pavvis classifies every gene into categories based on its presence pattern across all samples:

| Category | Definition |
|---|---|
| **Core** | Present in every sample |
| **Soft Core** | Present in ≥ 95% of samples |
| **Dispensable** | Present in ≥ 15% and < 95% of samples |
| **Private** | Present in > 0% and < 15% of samples |
| **Absent** | Absent from every sample |

Thresholds are configurable:

```python
pm = PavMatrix(
    matrix_path="pav_matrix.csv",
    metadata_path="metadata.csv",
    soft_core_min=0.90,
    dispensable_min=0.20,
    private_min=0.05,
)
print(pm.gene_counts)
# {'core': 12, 'soft_core': 84, 'dispensable': 22401, 'private': 6276, 'absent': 5522}
```

## Plots

### Presence & frequency

```python
from pavvis.plots import presence, variable_genes

# Frequency distribution coloured by gene category
variable_genes.frequency_histogram(pm).show()

# Gene counts per sample, grouped by a metadata column
presence.per_sample_boxplot(pm, gene_set="variable", x="species").show()

# Gene counts vs a continuous metadata variable
presence.per_sample_scatter(pm, x="latitude", color_by="species").show()
```

### Pangenome curves

```python
from pavvis.plots import curves

# Core & variable gene accumulation curves
curves.variable_gene_curve(pm, permutations=100).show()

# Variable gene accumulation grouped by metadata, in a defined order
curves.grouped_variable_curve(
    pm,
    column="region",
    group_order=["North", "Central", "South"],
    shade_groups=True,
).show()

# Jaccard similarity curve grouped by metadata
curves.jaccard_similarity_curve(
    pm,
    column="species",
    group_order=["Human", "Mouse", "Rat"],
).show()
```

### Exclusive genes

```python
from pavvis.plots import variable_genes

# Number of genes exclusive to each group
variable_genes.exclusive_genes_bar(pm, column="species").show()

# Or retrieve the raw counts
print(pm.exclusive_gene_counts("species"))
# {'Human': 42, 'Mouse': 17, 'Rat': 3}
```

### Heatmaps

```python
from pavvis.plots import heatmap

# Compute pairwise Jaccard similarity matrix (cached on the object)
pm.compute_jaccard()

# Binary PAV matrix heatmap with sample clustering and metadata annotation
heatmap.pav_heatmap(pm, gene_set="variable", col_cluster=True, color_by="species").show()

# Sample × sample Jaccard similarity heatmap
heatmap.sample_similarity_heatmap(pm, color_by="species", cluster=True).show()

```

### UMAP embedding

```python
from pavvis.plots import umap

# Compute embedding once (cached on the object, Jaccard metric by default)
pm.compute_umap(n_neighbors=15, min_dist=0.1, random_state=42)

umap.scatter(pm, color_by="species").show()      # 2D, discrete colour
umap.scatter(pm, color_by="latitude").show()     # 2D, continuous Viridis
umap.scatter3d(pm, color_by="species").show()    # 3D
umap.scatter3d(pm, color_by="species", depth=True).show()  # 3D with depth scaling
```

## Documentation

Full documentation including concepts, API reference, and getting-started guide is available at
**https://appliedbioinformatics.github.io/pavvis/**

To run the docs locally:

```bash
mkdocs serve
```
