# UMAP plots

## Description

Functions for visualising the UMAP embedding of samples computed by `pm.compute_umap()`. UMAP (Uniform Manifold Approximation and Projection) reduces the high-dimensional PAV matrix to two or three components, allowing the relationships between samples to be explored visually based on similarity in gene presence profiles.

!!! note "Prerequisites"
    You must call `pm.compute_umap()` before using any function on this page. See the [PavMatrix API](../pav_matrix.md) for details.

---

## scatter

### Description

Plots the first two UMAP components as a 2D scatter plot. Each point represents a sample. Points can be coloured by any metadata column — discrete columns produce a categorical legend, continuous columns use a Viridis colorscale. Clusters of samples in 2D space indicate groups with similar PAV profiles.

### Function

::: pavvis.plots.umap.scatter

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.umap as umap_plots

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

pm.compute_umap()
umap_plots.scatter(pm, color_by="clade").show()
```

--8<-- "docs/statics/plots/umap_scatter.html"

!!! note "Interpretation"
    Tight, well-separated clusters suggest that groups differ substantially in their PAV profiles. Overlapping
    clusters indicate genomic similarity between groups. Points sitting far from any cluster may represent
    samples with unusual gene content worth investigating further.

---

## scatter3d

### Description

Plots all three UMAP components as an interactive 3D scatter plot. Click and drag to rotate. Coloring behaves identically to `scatter()`. The third dimension can reveal structure that is collapsed in the 2D projection, particularly useful for datasets with complex population stratification.

### Function

::: pavvis.plots.umap.scatter3d

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.umap as umap_plots

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

pm.compute_umap()
umap_plots.scatter3d(pm, color_by="clade").show()
```

--8<-- "docs/statics/plots/umap_scatter3d.html"

!!! note "Interpretation"
    Groups that appear merged in the 2D plot may separate clearly when the third component is included. Rotating
    the plot can reveal layered or nested cluster structure that is not visible from the default viewing angle.