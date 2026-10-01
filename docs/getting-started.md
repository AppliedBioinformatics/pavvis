# Getting Started

## Installation

=== "pip"
    ``` bash
    pip install pavvis
    ```

=== "uv"
    ``` bash
    uv pip install pavvis
    ```

=== "Development"
    ``` bash
    git clone https://github.com/<your-username>/Pavvis.git
    cd Pavvis
    pip install -e ".[dev]"
    ```

## Loading your data

Pavvis requires two CSV files to generate a `PavMatrix` object. The `PavMatrix` class is central to the
Pavvis workflow.:

- **PAV matrix** — genes as rows, samples as columns, values must be `0` or `1`, i.e a binary matrix.
- **Metadata** — samples as rows, any number of metadata columns

``` py title="Building a PavMatrix object" linenums="1"
from pavvis import PavMatrix

pm = PavMatrix(
    matrix_path="path/to/pav_matrix.csv",
    metadata_path="path/to/metadata.csv",
)

print(pm)
# PavMatrix(34295 genes × 1039 samples)
```

!!! tip "File formats"
    For examples of how these files should be formatted, see the [concepts page](concepts.md).

When building the `PavMatrix()` class, Pavvis will automatically validate that:

- Both files exist on the system
- Sample IDs of the PAV matrix columns match the metadata row indexes
- All PAV values are binary (`0` or `1`) with no blank cells

## Exploring your data
Now you have generated a `PavMatrix` object, you can explore basic summary statistics using the class attributes.

``` py title="Explore a PavMatrix using basic class attributes" linenums="9"
# Get gene category counts
print(pm.gene_counts)
# {'core': 12, 'variable': 28761, 'absent': 5522}

# List gene IDs by category
print(pm.core_genes[:3])
print(pm.variable_genes[:3])
print(pm.absent_genes[:3])

# List exclusive variable genes for a metadata category
print(pm.exclusive_gene_counts("population_group"))
# {'GroupA': 183, 'GroupB': 204}
```

## First plot
Aside from allowing the user to easily group and summarise their data, the Pavvis package also provides functions that
allow the user to generate plots for thier data. All of these functions are contained in the `plots` module.

``` py title="Build histograms to summarise gene presence in your PavMatrix" "linenums="20"
import pavvis.plots as pp

# Histogram of presence frequencies (log y-axis)
pp.presence_frequency_histogram(pm).show()

# Bar chart of top 50 genes by presence frequency
pp.presence_frequency_histogram(pm, type="bar").show()
```
