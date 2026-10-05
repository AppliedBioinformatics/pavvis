# UMAP plots
This submodule contains functions for visualising the UMAP (Uniform Manifold Approximation and Projection) embedding
of PAV matrices computed by `pm.compute_umap()`. UMAP reduces the high-dimensional PAV matrix — where each gene is a
binary dimension — to two or three components, allowing the relationships between samples to be explored visually
based on similarity in gene presence profiles.

UMAP plots can provide the first sign of relationships in the PAV matrix across sample metadata groups, as each 
point on the umap can be coloured according to both discrete or continuous metadata values. UMAP analysis does not 
provide a statistical framework for making conclusions, but can be very informative in deciding future directions
of research based on the associations between metadata and the PAV matrix. 

!!! note "Prerequisites"
    You must call `pm.compute_umap()` before using any function documented on this page. See the [PavMatrix API](../pav_matrix.md)
    for details on available parameters such as `n_neighbors`, `min_dist`, and `metric`. This is a computationally 
    intensive function for large PAV matrices, but the resultant embeddings are cached to allow for quick rebuilding of
    UMAP plots once the embeddings have been generated. 

---

## 2D UMAP scatter plot

### Description

Plots the first two UMAP components as a 2D scatter plot. Each point represents a sample. Points can be coloured
by any metadata column — discrete columns produce a categorical legend; continuous columns use a Viridis colorscale.
Clusters of samples in 2D space indicate groups with similar PAV profiles.

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
    clusters indicate more similarity between groups. Points sitting far from any cluster may represent
    samples with unusual gene content worth investigating further. In the example plot, we see a small modern cluster
    far away from all other samples. Further analysis on this cluster revealed that these samples all contained an identical
    genomic introgression causing ~400 genes to be flagged as present across these samples, absent in the 
    rest of the sample set. 

---

## 3D UMAP scatter plot

### Description

Plots all three UMAP components as an interactive 3D scatter plot. Click and drag to rotate. Colouring behaves
identically to `scatter()`. The third dimension can reveal additional structure collapsed in the 2D projection,
particularly useful for datasets with complex population stratification and large number of samples.

Passing `depth=True` scales each point's size by its UMAP3 (z) coordinate, so points at higher z values
appear larger. This provides a visual depth cue that makes it easier to read the z-axis when viewing the plot
from a fixed angle, without needing to rotate to a favourable orientation.

!!! tip "Using depth scaling"
    `depth=True` is most effective when the default viewing angle places the z-axis running front-to-back.
    If the depth cue appears inverted, rotate the plot so that larger points are in the foreground.

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

#### Depth scaling example

``` py
umap_plots.scatter3d(pm, color_by="clade", depth=True).show()
```

!!! note "Interpretation"
    Groups that appear merged in the 2D plot may separate clearly when the third component is included. Rotating
    the plot can reveal a layered or nested cluster structure that is not visible from the default viewing angle.
    When `depth=True`, larger points indicate higher UMAP3 values — structure along the z-axis becomes readable
    even without rotating the plot.