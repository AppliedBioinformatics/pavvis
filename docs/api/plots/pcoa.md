# PCoA plots
This submodule contains functions for visualising the Principal Coordinates Analysis (PCoA) of PAV matrices computed by
`pm.compute_pcoa()`. PCoA (classical multidimensional scaling) places samples in a low-dimensional space so that the
distances between them approximate their pairwise Jaccard distances (1 − Jaccard similarity), ignoring shared absences.

Unlike UMAP, PCoA is a linear, deterministic ordination: each axis has an eigenvalue, so the proportion of variation it
explains can be reported, and distances on the plot relate directly to the distance matrix used by PERMANOVA and
dispersion tests.

!!! note "Prerequisites"
    You must call `pm.compute_pcoa()` before using any function documented on this page. See the
    [PavMatrix API](../pav_matrix.md) for details. Jaccard distance is not Euclidean, so only axes with positive
    eigenvalues are kept, and explained variance is reported relative to the sum of the positive eigenvalues.

---

## PCoA scatter plot

### Description

Plots two PCoA axes as a 2D scatter plot, with the percentage of variation explained in each axis title. Each point
represents a sample. Colouring works the same way as the UMAP scatter plot.

### Function

::: pavvis.plots.pcoa.scatter

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.pcoa as pcoa

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

pm.compute_pcoa()
pcoa.scatter(pm, color_by="clade").show()
pcoa.scatter(pm, color_by="clade", components=(1, 3)).show()
```

--8<-- "docs/statics/plots/pcoa_scatter.html"

---

## Scree plot

### Description

Bar chart of the percentage of variation explained by each computed axis. Use it to judge how much of the structure in
the distance matrix the plotted axes capture.

### Function

::: pavvis.plots.pcoa.scree

### Example

``` py
pm.compute_pcoa(n_components=10)
pcoa.scree(pm).show()
```

--8<-- "docs/statics/plots/pcoa_scree.html"