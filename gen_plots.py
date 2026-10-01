"""Generate demo Plotly HTML fragments for the MkDocs documentation.

Run this script once before building docs, and re-run whenever the demo
data or plotting functions change:

    python gen_plots.py

Output files land in docs/statics/plots/ and are inlined by the docs
pages via the pymdownx.snippets extension.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from pavvis import PavMatrix
import pavvis.plots.presence as presence
import pavvis.plots.curves as curves
import pavvis.plots.variable_genes as vg

DATA_DIR = Path("tests/data")
OUT_DIR = Path("docs/statics/plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

pm = PavMatrix(
    matrix_path=DATA_DIR / "pav.csv",
    metadata_path=DATA_DIR / "metadata.csv",
)

# Group order for grouped curves — must include all values in the column
GROUP_ORDER = ["Modern", "AG6", "AG5", "AG3", "AG2", "AG7", "AG4", "AG1"]


def _save(fig, name: str) -> None:
    html = fig.to_html(
        full_html=False,
        include_plotlyjs="cdn",
        config={"responsive": True},
    )
    out = OUT_DIR / f"{name}.html"
    out.write_text(html, encoding="utf-8")
    print(f"  wrote {out}")


print("Generating presence plots...")
_save(presence.frequency_histogram(pm), "frequency_histogram")
_save(presence.per_sample_boxplot(pm, x="continent"), "per_sample_boxplot")
_save(presence.per_sample_scatter(pm, x="depth", color_by="clade"), "per_sample_scatter")

print("Generating curve plots...")
_save(curves.variable_gene_curve(pm, permutations=5), "pangenome_curve")
_save(
    curves.grouped_variable_curve(pm, column="clade", group_order=GROUP_ORDER, permutations=5, shade_groups=True),
    "grouped_variable_curve",
)
_save(
    curves.jaccard_similarity_curve(pm, column="clade", group_order=GROUP_ORDER, permutations=5),
    "jaccard_similarity_curve",
)

print("Generating variable gene plots...")
_save(vg.frequency_histogram(pm, log_y=True), "vg_frequency_histogram")
_save(vg.exclusive_genes_bar(pm, column="clade"), "exclusive_genes_bar")

print("Generating UMAP plots...")
pm.compute_umap()
import pavvis.plots.umap as umap_plots
_save(umap_plots.scatter(pm, color_by="clade"), "umap_scatter")
_save(umap_plots.scatter3d(pm, color_by="clade"), "umap_scatter3d")

print("Done.")
