# UMAP plots

Functions for visualising the UMAP embedding of samples computed by `pm.compute_umap()`.

!!! note "Prerequisites"
    You must call `pm.compute_umap()` before using any function on this page.

---

## scatter

::: pavvis.plots.umap.scatter

### Example

2D projection coloured by `clade`.

--8<-- "docs/statics/plots/umap_scatter.html"

---

## scatter3d

::: pavvis.plots.umap.scatter3d

### Example

3D projection coloured by `clade`. Click and drag to rotate.

--8<-- "docs/statics/plots/umap_scatter3d.html"
