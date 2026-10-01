# Variable gene plots

## Description

Functions for visualising and exploring the variable gene pool within a PAV matrix. Variable genes — those present in at least one but not all samples — are often the most biologically informative component of a pangenome. These plots help characterise the distribution and exclusivity of variable gene content across groups.

---

## frequency_histogram

### Description

Plots gene presence frequency as a stacked histogram coloured by gene category (absent, private, dispensable, soft core, core). Each bin is subdivided by category, making it easy to see the composition of each frequency range at a glance. A distribution skewed toward low frequencies indicates the pangenome is open, with most variable genes present in only a small number of samples.

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

## exclusive_genes_bar

### Description

Plots the number of exclusive variable genes per group in a discrete metadata column. An exclusive gene is a variable gene present in at least one sample within a group and absent from every sample outside that group. A high exclusive gene count within a group may indicate lineage-specific content, such as adaptation to a local environment or unique functional repertoires.

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
    the rest of the population, or may simply have fewer samples represented.