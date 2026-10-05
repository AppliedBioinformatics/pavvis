# Heatmap plots
This submodule contains functions for visualising the PAV matrix and the relationships between its samples and genes
as heatmaps. Heatmaps provide a direct view of the underlying data that is complementary to association identification
offered by the presence, curve, and variable gene submodules. They are particularly useful for identifying clustering
structure, co-occurrence patterns, and the degree of genomic similarity across the sample set.

---

## PAV matrix heatmap

### Description

Renders the raw PAV matrix as a binary heatmap, with genes as rows and samples as columns. Each cell is
coloured white (absent, 0) or navy (present, 1). Because large PAV matrices can contain tens of thousands of genes,
the `gene_set` parameter (default `'variable'`) should be used to restrict the display to a biologically meaningful
subset.

Heatmaps can be noisy, and we recommend first spending time filtering the PAV matrix to reduce genes and samples to
increase the space available to render genes and samples of interest.

!!! info "Planned functionality"
    We plan to extend pavvis with functions for automatic matrix filtering, such as pruning genes based on their
    degree of presence across the matrix, to make heatmap generation more streamlined for large datasets.

Sample columns can optionally be reordered by hierarchical clustering on Jaccard distance (`col_cluster=True`), which
groups together samples with more similar PAV profiles. A colour-coded annotation strip can be added above the columns to
label each sample by a discrete metadata variable, making it easy to assess whether clustering aligns with known
biological groupings.

!!! note "Prerequisites for clustering"
    `col_cluster=True` requires `pm.compute_jaccard()` to have been called first. See the
    [PavMatrix API](../pav_matrix.md) for details.

### Function

::: pavvis.plots.heatmap.pav_heatmap

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.heatmap as heatmap

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

pm.compute_jaccard()
heatmap.pav_heatmap(pm, gene_set="variable", col_cluster=True, color_by="clade").show()
```

--8<-- "docs/statics/plots/pav_heatmap.html"

!!! note "Interpretation"
    Samples that cluster together share similar variable gene content. When the annotation strip is visible,
    alignment between colour-coded groups and the clustering structure suggests that the metadata variable is a
    good predictor of PAV profile similarity. Isolated samples or samples that do not cluster with their group
    may contain unusual gene content worth investigating further.

---

## Sample similarity heatmap

### Description

Plots a symmetric sample × sample Jaccard similarity matrix as a heatmap. Each cell shows the pairwise Jaccard
similarity between two samples, computed across all genes in the PAV matrix. Jaccard similarity accounts only for
genes that are present in at least one of the two samples, ignoring shared absences.

Samples can be reordered by hierarchical clustering (`cluster=True`, the default) to group together the most
similar samples. An optional annotation strip on both axes labels each sample by a discrete metadata column,
making it straightforward to assess whether biological grouping corresponds to genomic similarity in the PAV matrix.

!!! note "Prerequisites"
    This function requires `pm.compute_jaccard()` to have been called first. See the
    [PavMatrix API](../pav_matrix.md) for details.

### Function

::: pavvis.plots.heatmap.sample_similarity_heatmap

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.heatmap as heatmap

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

pm.compute_jaccard()
heatmap.sample_similarity_heatmap(pm, color_by="clade", cluster=True).show()
```

--8<-- "docs/statics/plots/sample_similarity_heatmap.html"

!!! note "Interpretation"
    Blocks of high similarity along the diagonal (after clustering) indicate groups of samples with closely
    related PAV profiles. A strong correspondence between the annotation strip colours and the similarity blocks
    suggests that the metadata variable captures meaningful biological differentiation. Low similarity between
    all samples indicates a highly open pangenome with a high degree of private gene content across the sample set.
