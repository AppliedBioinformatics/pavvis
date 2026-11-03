# Indicator gene plots
This submodule visualises the results of [`indicator_genes`](../stats/indicator.md). The plot takes the results table
as input, so the permutation test only has to be run once.

---

## Indicator gene heatmap

### Description

Shows the strongest indicator genes of each group as rows and the groups as columns. Each cell is the fraction of the
group's samples that carry the gene. Rows are grouped by the group each gene indicates, with a colour strip on the
left, so good indicators form a dark diagonal block: common in their own group and rare in the others. The cell for a
gene's own group is its fidelity (B); how pale the rest of the row is reflects its specificity (A).

### Function

::: pavvis.plots.indicator.heatmap

### Example

``` py
from pavvis import PavMatrix
from pavvis.stats.indicator import indicator_genes
import pavvis.plots.indicator as indicator_plots

pm = PavMatrix("path/to/pav.csv", "path/to/metadata.csv")
res = indicator_genes(pm, "clade", permutations=999, random_state=0)

indicator_plots.heatmap(pm, res, "clade", top_n=5).show()
```

--8<-- "docs/statics/plots/indicator_heatmap.html"

!!! note "Interpretation"
    A row that is dark in its own column and white elsewhere is a gene restricted to that group and carried by most of
    its samples. A row that is dark in its own column but also mid-tone in others indicates a gene enriched in, but not
    restricted to, the group. A row that is pale everywhere but darkest in its own column is a rare but group-specific
    gene. Gene IDs can be long; set `show_gene_labels=False` for a compact overview.

---

## Specificity vs fidelity

### Description

Plots the two components of the indicator value against each other for every significant indicator gene, with one
panel per group. Specificity (A) is on the x-axis and fidelity (B) on the y-axis. Grey curves join points of equal
IndVal.g statistic, so genes beyond the outer curve are the strongest indicators. Where the heatmap shows only the top
genes, this plot shows all of them.

### Function

::: pavvis.plots.indicator.specificity_fidelity

### Example

``` py
indicator_plots.specificity_fidelity(pm, res, "clade").show()
```

--8<-- "docs/statics/plots/indicator_specificity_fidelity.html"

!!! note "Interpretation"
    Genes in the upper right are carried by most samples of the group and rarely occur elsewhere. Genes on the right
    but low down are rare yet almost restricted to the group. Genes at the top but further left are common in the group
    and also common in others. A group whose cloud sits far to the right has a distinctive gene content.
