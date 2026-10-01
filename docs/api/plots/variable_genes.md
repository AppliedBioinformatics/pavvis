# Variable gene plots
This submodule contains functions for visualising and exploring the variable gene pool within a PAV matrix. Histograms
provide a good overview of the distribution of variable gene frequencies across the matrix. The degree of skew in the
matix can give insight into how "closed" a pangnome is. For example, a frequency histogram with most genes skewed to 
the left indicates that the majority of variation is more likely to be based on the individual samples rather than
group-based.

---

## frequency_histogram

::: pavvis.plots.variable_genes.frequency_histogram

### Example

Stacked bars show how absent, private, dispensable, soft-core, and core genes are distributed across frequency bins.

--8<-- "docs/statics/plots/vg_frequency_histogram.html"

Code used to generate the plot above: `frequency_histogram(pm, log_y=True)`.
The variable gene frequency histogram above indicates that for this pav matrix, the vast majority of genes have very 
low frequency of presence (only present in one or two samples), atypical compared to most pangenome PAV datasets.

---

## exclusive_genes_bar

::: pavvis.plots.variable_genes.exclusive_genes_bar

### Example

Each bar shows the number of variable genes exclusive to that clade — present in at least one sample of the group and absent from all other groups.

--8<-- "docs/statics/plots/exclusive_genes_bar.html"

Code used to generate the plot above: `pavvis.variable_genes.frequency_histogram(pm, column="clade")`
