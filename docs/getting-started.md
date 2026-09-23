# Getting Started

## Installation

```bash
pip install pavvis
```

For a development install from source:

```bash
git clone https://github.com/<your-username>/Pavvis.git
cd Pavvis
pip install -e ".[dev]"
```

## Loading data

Pavvis expects two CSV files:

- **PAV matrix** — genes as rows, samples as columns, values must be `0` or `1`
- **Metadata** — samples as rows, any number of metadata columns

```python
from pavvis import PavMatrix

pm = PavMatrix(
    matrix_path="path/to/pav_matrix.csv",
    metadata_path="path/to/metadata.csv",
)
print(pm)
# PavMatrix(34295 genes × 1039 samples)
```

On loading, Pavvis automatically validates that:

- Both files exist
- Sample IDs in the PAV matrix columns match the metadata row index
- All PAV values are binary (`0` or `1`) with no blank cells

## Exploring your data

```python
# Gene category counts
print(pm.gene_counts)
# {'core': 12, 'variable': 28761, 'absent': 5522}

# Lists of gene IDs by category
print(pm.core_genes[:3])
print(pm.variable_genes[:3])
print(pm.absent_genes[:3])

# Presence frequency (fraction of samples) per gene
print(pm.presence_frequency.describe())
```

## First plot

```python
import pavvis.plots as plots

# Histogram of presence frequencies (log y-axis)
plots.presence_frequency(pm).show()

# Bar chart of top 50 genes by presence frequency
plots.presence_frequency(pm, type="bar").show()
```
