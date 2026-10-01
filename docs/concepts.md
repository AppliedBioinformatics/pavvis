# Concepts

## What is a PAV matrix?

A **Presence/Absence variation (PAV) matrix** is a binary matrix commonly used in pangenome analysis
to represent whether a gene is present (`1`) or absent (`0`) across a population, often across thousands of genes.

This information can be stored in a number of different ways; however, Pavvis standardises the input data
structures based on the following format:

A gene PAV matrix in Pavvis can be loaded from a CSV file with the following format: 

|            | sample_a | sample_b | sample_c |
|------------|----------|----------|----------|
| gene_alpha |    1     |    0     |    1     |
| gene_beta  |    1     |    1     |    1     |
| gene_gamma |    0     |    0     |    1     |

- **Rows** — gene identifiers (must be unique). These can be any stable identifier: UniProt
  IDs, structural annotation IDs, nucleotide sequences, or custom labels.
- **Columns** — sample/accession identifiers (must be unique).
- **Values** — strictly `0` (absent) or `1` (present).

## PAV and associations with sample metadata
Pavvis is designed to integrate gene PAV analysis with sample level metadata. To do this, Pavvis requires the user
to supply a separate metadata.csv. This enables plots and summary statistics to be generated for both discrete and 
continuous variables. 

In this context, common sample metadata may include columns such as `population_group`, `ancestral_group`, or
`longitude`/`latitude`.

### The metadata.csv file
A separate **metadata CSV** file provides sample-level annotations (e.g. species, country of
origin, growth habit). Its row index must match the PAV matrix column headers exactly.

For example, a valid metadata CSV for the matrix above:

|          | population_group | climate   | longitude | latitude |
|----------|-----------------|-----------|-----------|----------|
| sample_a | GroupA          | Temperate | -0.27     | 51.49    |
| sample_b | GroupB          | Tropical  | -99.13    | 19.42    |
| sample_c | GroupA          | Temperate | 4.90      | 52.37    |

The metadata CSV may contain any number of columns with any mix of discrete or continuous values.

## Gene categories
Pangenomic studies often define gene sets based on their presence or absence across the whole population.
Pavvis classifies every gene into one of three categories based on its presence pattern across all samples in the matrix:

| Category | Definition |
|----------|------------|
| **Core** | Present in every sample (`1` across the entire row) |
| **Variable** | Present in at least one but not all samples |
| **Absent** | Absent from every sample (`0` across the entire row) |

### Variable genes
In PAV analysis, variable genes tend to be of the most interest to researchers. Pavvis classifies variable genes in a 
PAV matrix as one of three subcategories to add more depth to matrix analysis. These values are based on the relative
frequency of each genes presence across the entire gene row in the pav matrix. The default values for these three
subcategories are shown below:

| Variable gene subclass | Default presence frequency |
|------------------------|----------------------------|
| **Soft Core**          | >= 0.95 and < 1            |
| **Dispensable**        | >= 0.15 and < 0.95         |
| **Private**            | > 0 and < 0.15             |

!!! tip "Customising thresholds"
    These defaults are based on common thresholds used throughout pangenomic research papers. Depending on
    the number of samples in the PAV matrix, these thresholds may not be optimal — they can be
    overwritten when instantiating a new `PavMatrix()` object as shown below:

```py title="Building a PavMatrix object with custom variable gene thresholds"
from pavvis.pav_matrix import PavMatrix

pm = PavMatrix(matrix_path="./pav.csv", 
              metadata_path="./metadata.csv",
              soft_core_min=0.90,
              dispensable_min=0.20,
              private_min=0.05)

print(pm.thresholds)
```

### Absent and core genes
In addition to variable genes, Pavvis also defines two other types of gene:

- **absent** genes are always defined as all `0` across the PAV matrix.
- **core** genes are always defined as all `1` across the PAV matrix.

The default value for these gene categories cannot be changed.

### Exclusive genes

An **exclusive gene** is defined as a variable gene that is present in at least one sample belonging to a
specific group whilst bieng completely absent from every sample outside that group. Exclusivity is always
defined relative to a discrete metadata column (e.g. species, population, treatment).

For a gene to be classified as exclusive to group *G*:

1. It must be present (`1`) in **at least one** sample where `metadata[column] == G`.
2. It must be absent (`0`) in **every** sample where `metadata[column] != G`.

Core and absent genes are excluded from this analysis by definition — a core gene is present in
all samples and therefore cannot be exclusive to any group, and an absent gene is present in no
samples.

Exclusive genes are biologically interesting because they represent gene content that is unique
to a particular group — for example, genes found only in a specific species, geographic
population. A high exclusive gene count within a group may indicate lineage-specific gene content relating to processes
such as adaption to local environments.

## Why handle PAV matrices as Python class objects?

Pavvis represents a PAV matrix as a `PavMatrix` Python class rather than a loose
collection of functions. This allows the object to:

- **Validate on load** — bad inputs are caught immediately, not silently propagated.
- **Cache derived values** — computationally expensive results (i.e. UMAP embeddings,
  pangenome curves, calculated jaccard distances) can be stored on the object after the first calculation and utilised 
  by multiple different plotting functions.
- **Carry context** — the PAV data and its metadata always travel together, reducing the
  chance of mismatched inputs in multistep workflows.

!!! note "Example data"
    Throughout this documentation, we provide example data and plots that will allow users to visualise
    the outputs of the functions they are calling. This is particularly useful for the plotting subpackage, allowing
    users to see plots rendered for each function, as well as the input parameter values used. For generating these
    plots, we use a PAV matrix of variable genes that was generated as an outcome of the
    [Watkins Wheat Pangenome](https://pangenome.wheatgenome.info/) in one of our previous research projects.
    The metadata file associated with all Watkins wheat samples is derived from research published by
    [Cheng __et al__. (2024)](https://pubmed.ncbi.nlm.nih.gov/38885696/).