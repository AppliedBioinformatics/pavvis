# PERMANOVA and dispersion
This submodule tests whether groups of samples differ in gene content, using Jaccard distance (1 − Jaccard
similarity, which ignores shared absences). The two functions are the PAV equivalents of `vegan::adonis2` and
`vegan::betadisper` with `vegan::permutest`.

PERMANOVA asks whether group centroids differ more than expected by chance. It is also sensitive to differences in
spread, so a significant PERMANOVA should be read together with the dispersion test: if dispersions differ, the
PERMANOVA result may reflect spread rather than location.

To see both on one figure, use [`pcoa.centroids`](../plots/pcoa.md#group-centroids).

---

## PERMANOVA

### Function

::: pavvis.stats.permanova.permanova

### Example

``` py
from pavvis import PavMatrix
from pavvis.stats.permanova import permanova

pm = PavMatrix("path/to/pav.csv", "path/to/metadata.csv")
permanova(pm, "clade", permutations=9999, random_state=42)
```

```
            Df  SumOfSqs      R2        F  Pr(>F)
Model        7   20.1471  0.0859  13.8202  0.0001
Residual  1030  214.5040  0.9141      NaN     NaN
Total     1037  234.6511  1.0000      NaN     NaN
```

!!! note "Interpretation"
    `R2` is the fraction of the total variation in Jaccard distance explained by the grouping. In the example, clade
    explains 8.6% of the variation in gene content. Small R² values are typical for PAV data, because much of the
    variable gene content differs between individual accessions.

---

## Dispersion (PERMDISP)

### Function

::: pavvis.stats.permanova.betadisper

### Example

``` py
from pavvis.stats.permanova import betadisper

bd = betadisper(pm, "clade", permutations=9999, random_state=42)
bd["F"], bd["p_value"]
bd["distances"].groupby(pm.metadata["clade"]).mean()   # mean distance to centre per clade
```

!!! note "Interpretation"
    A significant result means some groups are more variable in gene content than others. Use
    [`pcoa.dispersion_boxplot`](../plots/pcoa.md#dispersion-boxplot) to see which.