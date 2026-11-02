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

---

## Group centroids

### Description

Plots samples on two PCoA axes with each group's centre (spatial median or mean) and a line from every sample to its
centre, like `plot(betadisper(...))` in vegan. Shaded convex hulls outline each group. With `stats=True` the plot is
annotated with the PERMANOVA and PERMDISP results, so the figure and the tests it illustrates come from one call.

Centres are computed in the full PCoA space, as in the dispersion test, and then projected onto the plotted axes.

### Function

::: pavvis.plots.pcoa.centroids

### Example

``` py
import pavvis.plots.pcoa as pcoa

pcoa.centroids(pm, "clade", permutations=999, random_state=0).show()
pcoa.centroids(pm, "clade", components=(1, 3), hull=False).show()
```

--8<-- "docs/statics/plots/pcoa_centroids.html"

!!! note "Interpretation"
    Centres far apart with short spider lines indicate groups that differ in location. Long lines in one group and short
    lines in another indicate unequal dispersion, which the PERMDISP result in the annotation tests formally.

---

## Dispersion boxplot

### Description

Distance of each sample to its group centre, one box per group: the equivalent of `boxplot(betadisper(...))`.

### Function

::: pavvis.plots.pcoa.dispersion_boxplot

### Example

``` py
pcoa.dispersion_boxplot(pm, "clade").show()
```

--8<-- "docs/statics/plots/pcoa_dispersion_boxplot.html"
