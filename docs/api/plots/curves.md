# Curve plots
This submodule contains functions for building traditional pangenome curves and other line-plots for gaining insight
into the variable gene profiles of a PAV matrix, including patterns that emerge from the accumulation of variable
genes across groups defined in a `metadata.csv` file.

---

## Pangenome Curves
Pangenome curves allow researchers to assess the degree of 'completeness' for a pangenome. Pangenomes seek to capture
all gene content shared across a population (all samples making up the the PAV matrix). A pangenome is made up of core
gene content present in all samples, and variable gene content present in only a subset of samples.

A pangenome curve plots the accumulation of variable gene content as each genome is "added" to the PAV matrix. As more
genomes are added, more of the variable material is captured, leading to a pattern of diminishing returns, where 
eventually, the addition of further genomes will add very little novel variable material. 

### Traditional pangenome curve.
A Traditional pangenome curve plots the accumulation of genomic content as each genome is added to the PAV matrix, curves
quickly reaching a plateau infer that the pangenome for the species is relatively closed, and that a complete picture
of all shared genomic content can be captured by a small number of samples. A curve that continues to rise without leveling
off indicates that the pangenome is "open" for the population, and that the sample set has not 
captured a complete picture of the shared genomic content for the population. 

### Description

Plots pangenome and core genome size curves on a shared y-axis as genomes are accumulated in random order. For each
permutation a different sample ordering is drawn. At each step $k$, the **pangenome size** is the count of genes with
at least one presence ($x_{g,i} = 1$) across the first $k$ samples; the **core genome size** is the count of genes
present in all $k$ samples. Both quantities are averaged over $P$ random permutations, and the mean and ±1 standard
deviation band are plotted for each curve. A pangenome curve that continues to rise steeply indicates an open pangenome;
a core curve that decays quickly to near-zero confirms that gene content is highly variable across the population.

!!! tip "Performance"
    A high number of permutations can be computationally intensive for large PAV matrices. We recommend starting
    with a lower value (e.g. `permutations=10`) to verify the plot looks as expected, before increasing to a
    higher value for final figures.

### Function

::: pavvis.plots.curves.legacy_pangenome_curve

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.curves as curves

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

curves.legacy_pangenome_curve(pm, permutations=50).show()
```

--8<-- "docs/statics/plots/legacy_pangenome_curve.html"

!!! note "Interpretation"
    A pangenome curve that continues to rise without levelling off as more genomes are added indicates that the pangenome
    is **open** — additional sampling would likely reveal further novel gene content. A curve that plateaus early suggests
    a **closed** pangenome, where the gene repertoire of the population is well-captured by the existing sample set.
    Correspondingly, a core curve that falls steeply toward zero indicates that few genes are universally conserved across
    all samples, which is typical of species with highly variable accessory genomes.

---

### Variable gene accumulation curves.
A variable gene accumulation curve has no fundamental conceptual difference to a traditional pangenome curve,
but seeks to plot the accumulation of variable gene content separate to the core genetic material. As the variable gene
content increases with the addition of new genomes to the pangenome, the core genetic material decreases. This can be 
visualised as two distinct curves across the same plot, represented using 3-axis.

At each step $k$, a random ordering of samples is drawn. The **pangenome size** is the count of genes present in at
least one of the first $k$ samples; the **core size** is the count of genes present in all $k$ samples; and the
**variable gene count** is the difference between the two. This is repeated across $P$ random permutations, and the
mean and ±1 standard deviation band are plotted for both curves.

!!! tip "Performance"
    A high number of permutations can be computationally intensive for large PAV matrices. We recommend starting
    with a lower value (e.g. `permutations=10`) to verify the plot looks as expected, before increasing to a
    higher value for final figures.

### Description

Plots core and variable gene counts on separate left/right y-axes as genomes are accumulated in random order. 
Each permutation draws a different sample ordering; the mean and ±1 std band across all permutations are shown. A 
steeply decaying core curve indicates an open pangenome; a flattening variable curve suggests the gene pool is 
near-saturated.

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
    is open — adding new samples would likely reveal additional gene content. A variable gene curve that presents a 
    flat gradient indicates that the majority of the variable gene content has been captured for the population by the
    number of samples.
    In this case, the PAV matrix contains no core genes, so we expect the number of core genes to decrease quickly to 
    near-zero as new variable material is added across the plot. 

---

## Grouped variable gene accumulation curves.

### Description

Plots variable gene accumulation with samples added group-by-group in a user-defined order. Within each group, the sample
order is randomly permuted `P` times whilst the group sequence itself is fixed. Vertical dashed lines mark group
boundaries, making it easier to assess how much new variable gene content each group contributes to the pangenome.
Adding genomes to the graph one group at a time may allow the researcher to identify groups that provide variable gene
content not already captured in the shared material of previous groups and identify the existence of group-specific
variable gene content.

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

GROUP_ORDER = ["Modern", "AG6", "AG5", "AG3", "AG2", "AG7", "AG4", "AG1"] # Metadata group values for "Clade"

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

## Jaccard Similarity Curves
The Jaccard similarity between two samples is the proportion of genes that are shared between them relative to the total
number of genes present in either. Unlike simple overlap counts, it accounts for the size of each sample's gene set,
making it a useful measure of overall genomic similarity. Because Jaccard similarity is calculated only over genes that
are present in at least one sample, jaccard-distance does not take into account shared absence between samples.

Plotting mean pairwise Jaccard similarity as groups are added one at a time allows the researcher to assess how
genomically distinct each group is from those already accumulated. A drop in similarity when a new group is introduced
suggests that group carries a substantially different gene repertoire; a stable curve indicates the groups are
broadly similar in their PAV profiles.

### Description

Plots mean pairwise Jaccard similarity across all accumulated genomes as groups are added in a fixed order. Within each
group the sample order is randomly permuted `P` times whilst the group sequence itself is fixed. At each step $k$, the
mean Jaccard similarity is calculated across all pairs of the $k$ genomes accumulated so far. The mean and ±1 standard
deviation band across permutations are shown. A drop at a group boundary indicates the incoming group is divergent from
those already accumulated; a stable curve indicates genomic similarity across groups.

!!! tip "Use case"
    This type of curve can be useful for assessing particular groups that may harbour interesting and novel variable
    genes, rare or absent across other groups, that may warrant further investigation.

!!! tip "Performance"
    A high number of permutations can be computationally intensive for large PAV matrices. We recommend starting
    with a lower value (e.g. `permutations=10`) to verify the plot looks as expected, before increasing to a
    higher value for final figures.

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

GROUP_ORDER = ["Modern", "AG6", "AG5", "AG3", "AG2", "AG7", "AG4", "AG1"] # Metadata group values for "Clade"

curves.jaccard_similarity_curve(
    pm,
    column="clade",
    group_order=GROUP_ORDER,
    permutations=50,
).show()
```

--8<-- "docs/statics/plots/jaccard_similarity_curve.html"

!!! note "Interpretation"
    A sharp drop at a group boundary indicates that the incoming group carries a relatively different gene
    repertoire from those already added to the curve. A curve that remains broadly stable across the boundaries of groups
    suggests that the groups are genomically similar in terms of their PAV profiles.