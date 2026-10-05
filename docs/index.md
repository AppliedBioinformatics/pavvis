<img src="statics/banner.jpg" alt="pavvis" class="pavvis-banner" />

**pavvis** is a Python package for analysing and visualising gene Presence/Absence Variation (PAV) data from pangenome studies.

pavvis was designed to be easy to use and optimised for large-scale genomic datasets. To get the best use out of this
Python package, we recommend reading through our [concepts](concepts.md) page, that will help you to understand the data structures
associated with all analysis in this Python package.

## Plotting modules

| Module | Description |
|---|---|
| [Presence](api/plots/presence.md) | Gene presence frequency histograms and per-sample gene count plots |
| [Curves](api/plots/curves.md) | Pangenome accumulation curves and Jaccard similarity curves |
| [Variable Genes](api/plots/variable_genes.md) | Frequency distributions and exclusive gene bar charts |
| [UMAP](api/plots/umap.md) | 2D and 3D UMAP embeddings of sample PAV profiles |
| [Heatmap](api/plots/heatmap.md) | PAV matrix and sample similarity heatmaps |

---

## Quick start

```py title="Getting started with pavvis." linenums="1"
from pavvis import PavMatrix
import pavvis.plots as plots

pm = PavMatrix("pav_matrix.csv", "metadata.csv")
print(pm)
# PavMatrix(34295 genes × 1039 samples)

print(pm.gene_counts)
# {'core': 12, 'variable': 28761, 'absent': 5522}

plots.presence_frequency_histogram(pm).show()
```
