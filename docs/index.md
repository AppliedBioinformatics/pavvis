# Pavvis

**Pavvis** is a Python package for analysing and visualising gene Presence/Absence Variation (PAV) data from pangenome studies.

Pavvis was designed to be easy to use and optimised for large-scale genomic datasets. To get the best use out of this
Python package, we reccommend reading through our [concepts](concepts.md) page, that will help you to understand the data structures
associated with all analysis in this Python package.

## Quick start

```python
from pavvis import PavMatrix
import pavvis.plots as plots

pm = PavMatrix("pav_matrix.csv", "metadata.csv")
print(pm)
# PavMatrix(34295 genes × 1039 samples)

print(pm.gene_counts)
# {'core': 12, 'variable': 28761, 'absent': 5522}

plots.presence_frequency_histogram(pm).show()
```
