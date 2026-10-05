# Variable gene plots
This submodule contains plotting functions for visualising and exploring the variable gene pool within a PAV matrix.
Variable genes — those present in at least one but not all samples — are often the most biologically informative
component of a pangenome. These plots help characterise the distribution and exclusivity of variable gene content
across groups defined in a `metadata.csv` file.

---

## Presence frequency distribution

### Description

The frequency of presence for a single gene can be calculated using the following calculation:

$$
\text{presence frequency} = \frac{\text{number of samples where gene} = 1}{\text{total number of samples}}
$$

The result is a value in the range `[0, 1]`, where `0` indicates a gene is absent in all samples and `1` means it is
present in every sample. After presence frequency has been calculated for each gene, a histogram can be plotted to
summarise the distribution of presence across all genes.

Unlike the [presence frequency histogram](presence.md#presence-frequency-histograms) in the `presence` submodule,
this plot colours each bin by gene category — absent, private, dispensable, soft core, and core — so the composition
of each frequency range is visible at a glance.

This plot can be useful as a first-pass overview of the variable gene pool — a heavily left-skewed distribution
indicates most variable genes are rare across the sample set, which is typical of an open pangenome. A distribution
concentrated at higher frequencies suggests a more closed pangenome where most variable genes are shared broadly
across the sample set.

### Function

::: pavvis.plots.variable_genes.frequency_histogram

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.variable_genes as vg

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

vg.frequency_histogram(pm, log_y=True).show()
```

--8<-- "docs/statics/plots/vg_frequency_histogram.html"

!!! note "Interpretation"
    The stacked colouring reveals the category composition at each frequency bin. A dominant private (low-frequency)
    bar suggests the pangenome is open and highly variable. A dominant soft-core or core bar at high frequencies
    suggests a more closed pangenome where most genes are conserved across samples.

---

## Exclusive genes

### Description

An exclusive gene is a variable gene present in at least one sample within a group and absent from every sample
outside that group. Plotting exclusive gene counts per group allows the researcher to identify which groups harbour
gene content that is not shared with any other group in the dataset.

A high exclusive gene count within a group may indicate lineage-specific content, such as adaptation to a local
environment, unique metabolic capabilities, or a functional repertoire distinct from the rest of the sample set.
Conversely, groups with very few exclusive genes may be more genomically similar to the broader population or may
have fewer samples represented.

### Function

::: pavvis.plots.variable_genes.exclusive_genes_bar

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.variable_genes as vg

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

vg.exclusive_genes_bar(pm, column="clade").show()
```

--8<-- "docs/statics/plots/exclusive_genes_bar.html"

!!! note "Interpretation"
    Tall bars indicate groups with a large number of genes found nowhere else in the dataset, suggesting
    lineage-specific gene content. Groups with very few exclusive genes may be more genomically similar to
    the rest of the population or may simply have fewer samples represented.