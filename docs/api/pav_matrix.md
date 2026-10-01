# PavMatrix

::: pavvis.pav_matrix.PavMatrix
    options:
      members: false

---

## Constructor
The `PavMatrix()` class is the core data structure behind Pavvis. By default, all that is required to generate it is
a pav.csv and accompanying metadata.csv file. 

!!! tip "File formats"
    For more information on how these input files should be structured, see the [concepts page](../concepts.md).

::: pavvis.pav_matrix.PavMatrix.__init__

---

## Identity & IDs

::: pavvis.pav_matrix.PavMatrix.__repr__

::: pavvis.pav_matrix.PavMatrix.sample_ids

::: pavvis.pav_matrix.PavMatrix.gene_ids

::: pavvis.pav_matrix.PavMatrix.metadata_columns

---

## Gene categories

::: pavvis.pav_matrix.PavMatrix.gene_counts

::: pavvis.pav_matrix.PavMatrix.core_genes

::: pavvis.pav_matrix.PavMatrix.absent_genes

::: pavvis.pav_matrix.PavMatrix.variable_genes

---

## Variable gene subcategories

::: pavvis.pav_matrix.PavMatrix.thresholds

::: pavvis.pav_matrix.PavMatrix.soft_core_genes

::: pavvis.pav_matrix.PavMatrix.dispensable_genes

::: pavvis.pav_matrix.PavMatrix.private_genes

::: pavvis.pav_matrix.PavMatrix.presence_frequency

---

## Exclusive genes

::: pavvis.pav_matrix.PavMatrix.exclusive_gene_counts

---

## UMAP
Uniform Manifold Approximation Projection (UMAP) is a dimensionality reduction technique that is used to 
visualise high-dimensional data. In the context of PAV, it can help gain insight into the relationships between samples
through similarity in PAV profiles. Generating UMAP embeddings can be done by calling the `compute_umap()` method. This
is not automatically called on initial construction of new `PavMatrix()` objects, as it can be computationally 
intensive, particularly for large datasets.

::: pavvis.pav_matrix.PavMatrix.compute_umap