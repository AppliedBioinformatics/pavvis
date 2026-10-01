# Gene Presence
This submodule contains plotting functions for visualising gene presence frequency across samples and per-sample gene
counts across the PAV matrix. Assessment of gene presence can help to identify outlier samples and provide insights into
the distribution of gene content across the dataset in relation to sample metadata.

---
## Presence frequency histograms

### Description
The frequency of presence for a single gene can be calculated using the following calculation:

$$
\text{presence frequency} = \frac{\text{number of samples where gene} = 1}{\text{total number of samples}}
$$

The result is a value in the range `[0, 1]`, where `0` indicates a gene is absent in all samples and `1` means it is 
present in every sample.

After `presence frequency` has been calculated for each gene in the PAV matrix, a histogram can be plotted to summarise
the distribution of presence across all genes, the larger the number of histogram bins, the larger the resolution 
available in the final plot. 

This plot can be useful as a first-pass overview of the PAV matrix — a heavily left-skewed distribution indicates 
most genes are rare across the sample set. A heavily right-skewed distribution indicates most genes are 
present in the majority of samples.

### Function

::: pavvis.plots.presence.frequency_histogram

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.presence as presence

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

presence.frequency_histogram(pm).show()
```

--8<-- "docs/statics/plots/frequency_histogram.html"


!!! note "Interpretation"
    In the above plot, we see a right-skewed distribution of presence frequencies. This indicates that most genes in this
    PAV matrix are rare across the sample set — the majority of genes are only present in a few samples. It also confirms that
    there are very few core genes in this PAV matrix, which we could also assess using the class attribute `pm.core_genes`.

---

## Total presence and discrete metadata

### Description

The number of present genes may differ according to sample metadata. A box-plot is a good way to represent the spread
of present genes across all samples for a given metadata column. For each sample, the total number of present genes is
calculated, and box plots are generated across each group in the metadata column. This plot can be useful for assessing
the occurrence of batch effects and identifying groups of samples that may contain relatively diverse gene content.

!!! note "Restriction on column type"
    This plot will only work if the metadata column is a discrete variable. To visualise relationships between total
    gene presence and continuous variables, see the 
    [scatter plot](#total-presence-and-continuous-metadata) function below.

### Function

::: pavvis.plots.presence.per_sample_boxplot

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.presence as presence

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

presence.per_sample_boxplot(pm, x="continent").show()
```

--8<-- "docs/statics/plots/per_sample_boxplot.html"

!!! note "Interpretation"
    In this plot, we can see that several samples from Asia and Europe contain a total number of present genes that are
    outside the `q3` range for the boxplots. These samples may warrant further investigation to assess whether the 
    total gene presence is biologically accurate, or whether this could be an artefact of some preceding data analysis.

---

## Total presence and continuous metadata

### Description

To visualise the relationship between total gene presence for each sample and continuous metadata variables, scatter plots
are a better option. A particularly useful plot involves plotting sampling depth across the x-axis. This allows the researcher
to identify samples that are lacking adequate read depth for calling PAV, as samples at lower depths will have less total 
genes called present. Points can be coloured according to a second metadata column, which may be useful for assessing
relationships against a second variable with the `color_by` parameter.

### Function

::: pavvis.plots.presence.per_sample_scatter

### Example

``` py
from pavvis import PavMatrix
import pavvis.plots.presence as presence

pm = PavMatrix(
    matrix_path="path/to/pav.csv",
    metadata_path="path/to/metadata.csv",
)

presence.per_sample_scatter(pm, x="depth", color_by="clade").show()
```

--8<-- "docs/statics/plots/per_sample_scatter.html"

!!! note "Interpretation"
    In this example, we can see that across the sample set, presence is mostly uniform across depth. We do see some
    evidence of higher presence in very high-depth samples such as that of the right-most point. In this case, before 
    generating our final PAV matrix, we removed all samples with a depth of less than 10x as part of our standard 
    protocol, if these samples were included, we would expect to see these points render at a much lower position along 
    the plots y-axis. 