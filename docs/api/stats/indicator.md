# Indicator genes
This submodule finds genes whose presence is characteristic of one group of samples, using the indicator value
(IndVal) of Dufrêne & Legendre (1997) in its group-size-equalised form, IndVal.g. It is the PAV equivalent of
`indicspecies::multipatt(func="IndVal.g", duleg=TRUE)` in R.

For a gene and a group, two components are combined:

- **A (specificity):** how concentrated the gene is in the group. It is the gene's presence frequency in the group
  divided by the sum of its frequencies across all groups, so unequal group sizes do not bias it.
- **B (fidelity):** the fraction of the group's samples that carry the gene.

The statistic is $\sqrt{A \cdot B}$, between 0 and 1. A value of 1 means the gene is present in every sample of the
group and absent from every other group. Each gene is assigned to the group where its value is highest, and a
permutation test asks whether that value is larger than expected if group labels were random.

---

## Indicator gene analysis

### Function

::: pavvis.stats.indicator.indicator_genes

### Example

``` py
from pavvis import PavMatrix
from pavvis.stats.indicator import indicator_genes

pm = PavMatrix("path/to/pav.csv", "path/to/metadata.csv")

res = indicator_genes(pm, "clade", permutations=9999, random_state=42)

sig = res[res.q_value < 0.05]
sig.groupby("group").size()                 # number of indicator genes per clade
sig[sig.group == "Modern"].head(20)         # strongest Modern indicators
res.join(pm.gene_table(column="clade"))     # add per-group presence counts
```

To plot the top indicator genes of each clade, pass `res` to
[`plots.indicator.heatmap`](../plots/indicator.md).

!!! note "Interpretation"
    A high `stat` with high `B` but moderate `A` describes a gene common in the group that also occurs elsewhere. A high
    `A` with low `B` describes a gene that is rare but almost restricted to the group. With 999 permutations the smallest
    possible p-value is 0.001, which limits how many genes can pass a strict FDR threshold in a large matrix. Use more
    permutations when testing tens of thousands of genes.

!!! tip "Which genes are tested"
    By default only variable genes are tested. Core and absent genes cannot be indicators, and including them only weakens
    the multiple-testing correction. Pass `genes=pm.gene_ids` to test every gene, for example to match an R analysis run on
    the full matrix.