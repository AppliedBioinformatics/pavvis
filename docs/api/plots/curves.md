# Curve plots
This submodule contains functions for building traditional pangenome curves and other line-plots that can be used t
gain insight into the variable gene profiles of samples within the PAV matrix. 

---

## Pangenome Curves
A traditional pangenome curve is a plot that allows the researcher to assess the degree of completeness for a pangenome.
Pangenomes seek to capture all of the gene content shared across a population (the samples in the PAV matrix). This is 
made up of core gene content present in all samples and variable gene content present in only a subset of samples.

A pangenome curve plots the accumulation of variable gene content as each genome is "added" to the PAV matrix. As more
genomes are added, more of the variable material is captured, leading to a pattern of diminishing returns, where 
eventually, the addition of further genomes will add very little new variable material. 

At each step $k$, a random ordering of samples is drawn. The **pangenome size** is the number of genes seen in at least
one of the first $k$ samples:

$$
\text{Pan}(k) = \left| \{ g : \exists\, i \leq k,\; x_{g,i} = 1 \} \right|
$$

The **core size** is the number of genes present in all $k$ samples so far:

$$
\text{Core}(k) = \left| \{ g : \forall\, i \leq k,\; x_{g,i} = 1 \} \right|
$$

The **variable gene count** is the difference between the two:

$$
\text{Variable}(k) = \text{Pan}(k) - \text{Core}(k)
$$

This is repeated across $P$ random permutations and the mean and ±1 standard deviation band are plotted for both curves.

!!! tip "Performance"
    A high number of permutations can be computationally intensive for large PAV matrices. We recommend starting
    with a lower value (e.g. `permutations=10`) to verify the plot looks as expected, before increasing to a
    higher value for final figures.

### Description

Plots core and variable gene counts on separate left/right y-axes as genomes are accumulated in random order. Each permutation draws a different sample ordering; the mean and ±1 std band across all permutations are shown. A steeply decaying core curve indicates an open pangenome; a flattening variable curve suggests the gene pool is near-saturated.

### Function

::: pavvis.plots.curves.variable_gene_curve

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.curves as curves

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

curves.variable_gene_curve(pm, permutations=50).show()
```

--8<-- "docs/statics/plots/variable_gene_curve.html"

!!! note "Interpretation"
    If the variable gene curve (right axis) is still rising steeply as more genomes are added, the pangenome
    is `open` — adding new samples would likely reveal additional gene content. A variable gene curve that presents a 
    flat gradient indicates that the majority of the variable gene content has been captured for the population by the
    number of samples.
    In this case, the PAV matrix contains no core genes, so we expect the number of core genes to decrease quickly to 
    near-zero as new variable material is added across the plot. 

---

## grouped_variable_curve

### Description

Plots variable gene accumulation with samples added group-by-group in a user-defined order. Within each group the sample order is randomly permuted; the group sequence itself is fixed. Vertical dashed lines mark group boundaries, making it easy to see how much new variable gene content each group contributes.

### Function

::: pavvis.plots.curves.grouped_variable_curve

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.curves as curves

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

GROUP_ORDER = ["Modern", "AG6", "AG5", "AG3", "AG2", "AG7", "AG4", "AG1"]

curves.grouped_variable_curve(
    pm,
    column="clade",
    group_order=GROUP_ORDER,
    permutations=50,
).show()
```

--8<-- "docs/statics/plots/grouped_variable_curve.html"

!!! note "Interpretation"
    A steep rise after a group boundary indicates that the incoming group contributes a large number of genes
    not seen in any previous group. A flat section suggests the new group is genomically similar to those
    already accumulated and adds little novel gene content.

---

## jaccard_similarity_curve

### Description

Plots mean pairwise Jaccard similarity across all accumulated genomes as groups are added in a fixed order. A drop at a group boundary indicates the incoming group is divergent from those already accumulated; a stable curve indicates similarity. The Jaccard metric ignores shared absences, making it well suited to sparse PAV data.

### Function

::: pavvis.plots.curves.jaccard_similarity_curve

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.curves as curves

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

GROUP_ORDER = ["Modern", "AG6", "AG5", "AG3", "AG2", "AG7", "AG4", "AG1"]

curves.jaccard_similarity_curve(
    pm,
    column="clade",
    group_order=GROUP_ORDER,
    permutations=50,
).show()
```

--8<-- "docs/statics/plots/jaccard_similarity_curve.html"

!!! note "Interpretation"
    A sharp drop at a group boundary indicates that the incoming group carries a substantially different gene
    repertoire from those already added. A curve that remains broadly stable across boundaries suggests the
    groups are genomically similar in terms of their PAV profiles.