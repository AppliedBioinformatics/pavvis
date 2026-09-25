# Concepts

## What is a PAV matrix?

A **Presence/Absence Variation (PAV) matrix** is a binary matrix commonly used in pangenome analysis
to represent whether each gene is present (`1`) or absent (`0`) across a population.

This information can be stored in a number of different ways; however, Pavvis standardises the input data
structures based on the following format:

A gene PAV matrix in Pavvis can be loaded from a CSV file with the following format: 

```
             sample_a  sample_b  sample_c
gene_alpha      1         0         1
gene_beta       1         1         1
gene_gamma      0         0         1
```

- **Rows** — gene identifiers (must be unique). These can be any stable identifier: UniProt
  IDs, structural annotation IDs, nucleotide sequences, or custom labels.
- **Columns** — sample/accession identifiers (must be unique).
- **Values** — strictly `0` (absent) or `1` (present).

## How does Pavvis handle associations with sample metadata?
A separate **metadata CSV** provides sample-level annotations (e.g. species, country of
origin, growth habit). Its row index must match the PAV matrix column headers exactly.

For example, a valid metadata CSV for the matrix above:

```
          population_group  climate  longitude  latitude
sample_a       GroupA        Temperate   -0.27     51.49
sample_b       GroupB        Tropical    -99.13    19.42
sample_c       GroupA        Temperate    4.90     52.37
```

The metadata CSV may contain any number of columns with any mix of discrete or continuous values.

## Gene categories
Pangenomic studies often define gene sets based on their presence or absence across the whole population.
Pavvis classifies every gene into one of three categories based on its presence pattern:

| Category | Definition |
|----------|------------|
| **Core** | Present in every sample (`1` across the entire row) |
| **Variable** | Present in at least one but not all samples |
| **Absent** | Absent from every sample (`0` across the entire row) |

In PAV analysis, variable genes tend to be of the most interest to researchers. Pavvis classifies variable genes in a 
PAV matrix as one of three subcategories to make visualisation easier. These values are based on the frequency of 
gene presence across the entire gene row in the pav matrix. The default values for these subcategories are shown below:

| Variable gene subclass | Default presence frequency |
|------------------------|----------------------------|
| **Soft Core**          | >= 0.95 and < 1            |
| **Dispensable**        | >= 0.15 and < 0.95         |
| **Private**            | > 0 and < 0.15             |

These defaults are based on common thresholds used throughout pangenomic research papers. We realise that depending on
the number of samples in the PAV matrix, these thresholds may not be optimal for the user and the defaults can be
overwritten by following the instructions below when instantiating a new PavMatrix() object.

```python
from pavvis.pav_matrix import PavMatrix

pm = PavMatrix(matrix_path="./pav.csv", 
              metadata_path="./metadata.csv",
              soft_core_min=0.90,
              dispensable_min=0.20,
              private_min=0.05)

print(pm.thresholds)
```

**absent** genes are always defined as all `0` across the PAV matrix.
**core** genes are always defined as all `1` across the PAV matrix.

## Why handle PAV matrices as Python class objects?

Pavvis represents a PAV matrix as a `PavMatrix` Python class rather than a loose
collection of functions. This allows the object to:

- **Validate on load** — bad inputs are caught immediately, not silently propagated
- **Cache derived values** — computationally expensive results (UMAP embeddings,
  pangenome curves) can be stored on the object after the first calculation
- **Carry context** — the PAV data and its metadata travel together, reducing the
  chance of mismatched inputs in multistep workflows
