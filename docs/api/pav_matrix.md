# PavMatrix

::: pavvis.pav_matrix.PavMatrix
    options:
      members: false

---

## Constructor
The `PavMatrix()` class is the core data structure behind pavvis. By default, all that is required to generate it is
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

## Subsetting
`subset()` returns a new `PavMatrix` restricted to selected samples and/or genes, so any analysis or plot can be rerun
on part of the data — for example a single clade or region. Samples can be selected by ID, by a pandas query against
the metadata, or both. The original object is never modified.

::: pavvis.pav_matrix.PavMatrix.subset

```py
modern = pm.subset(query="clade == 'Modern'")
europe_deep = pm.subset(query="continent == 'Europe' and depth > 10")
pair = pm.subset(samples=["Ability", "Agron"], drop_absent=True)

print(modern)  # PavMatrix(34295 genes × 212 samples)
```

!!! note "Gene categories are recomputed"
    Categories are always relative to the samples in the matrix. A gene that is variable across the full population
    may be core, or absent, within a subset. Use `drop_absent=True` to remove genes that no retained sample carries.

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

## Summary tables
`gene_table()` and `sample_table()` collect per-gene and per-sample summaries into pandas DataFrames, ready to filter,
sort, merge with other data (such as functional annotations), or export with `.to_csv()` as supplementary tables.

::: pavvis.pav_matrix.PavMatrix.gene_table

```py
gt = pm.gene_table(column="clade")

gt.loc["gene_42"]                                  # everything about one gene
gt[gt.exclusive_to == "AG3"].index                 # genes found only in AG3
gt[(gt.category == "dispensable") & (gt.n_Modern == 0)]  # dispensable genes absent from Modern lines
gt.to_csv("gene_table.csv")
```

::: pavvis.pav_matrix.PavMatrix.sample_table

```py
st = pm.sample_table(include_metadata=True)

st.sort_values("n_singletons", ascending=False).head(10)  # accessions carrying the most unique genes
st.groupby("clade")["n_private"].mean()                   # private gene content per clade
st["n_genes"].corr(st["depth"])                           # check gene counts against sequencing depth
```

!!! tip "Unclassified genes"
    By default `private_min` is 0, so every gene present in at least one sample gets a category. If you raise it
    (e.g. `private_min=0.01`), genes rarer than that threshold are no longer private, but they aren't absent either.
    `gene_table()` labels these genes `"unclassified"`.

---

## UMAP
Uniform Manifold Approximation Projection (UMAP) is a dimensionality reduction technique that is used to 
visualise high-dimensional data. In the context of PAV, it can help gain insight into the relationships between samples
through similarity in PAV profiles. Generating UMAP embeddings can be done by calling the `compute_umap()` method. This
is not automatically called on initial construction of new `PavMatrix()` objects, as it can be computationally 
intensive, particularly for large datasets.

::: pavvis.pav_matrix.PavMatrix.compute_umap