# Curve plots

Functions for accumulation curves that reveal pangenome openness, variable gene growth, and inter-sample similarity.

---

## pangenome_curve

::: pavvis.plots.curves.pangenome_curve

### Example

50 permutations, linear scale. Core genes decay on the left axis; variable genes accumulate on the right.

--8<-- "docs/statics/plots/pangenome_curve.html"

---

## grouped_variable_curve

::: pavvis.plots.curves.grouped_variable_curve

### Example

Groups are added in the order: Modern → AG6 → AG5 → AG3 → AG2 → AG7 → AG4 → AG1. Vertical dashed lines mark each group boundary.

--8<-- "docs/statics/plots/grouped_variable_curve.html"

---

## jaccard_similarity_curve

::: pavvis.plots.curves.jaccard_similarity_curve

### Example

Same group order as above. A drop at a boundary indicates the incoming group is divergent from those already accumulated.

--8<-- "docs/statics/plots/jaccard_similarity_curve.html"
