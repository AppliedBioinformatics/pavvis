from pathlib import Path

from pavvis import PavMatrix
import pavvis.plots.presence as pp
import pavvis.plots.curves as pc
import pavvis.plots.variable_genes as pvg
import pavvis.plots.umap as pum

DATA_DIR = Path(__file__).parent / "tests" / "data"

pm = PavMatrix(DATA_DIR / "pav.csv", DATA_DIR / "metadata.csv")

print(f"Gene counts: {pm.metadata_columns}")

#pp.frequency_histogram(pm, n_bins=50).show()
#pp.per_sample_boxplot(pm, x="clade", gene_set="variable").show()
#pp.per_sample_scatter(pm, x="depth", gene_set="variable").show()
#pc.legacy_pangenome_curve(pm, permutations=10).show()
#pc.pangenome_curve(pm, permutations=10).show()
#pc.grouped_variable_curve(pm, column="clade",
#                          group_order=["AG1", "AG2", "AG3", "AG4", "AG5", "AG6", "AG7", "Modern"],
#                          permutations=100,
#                          shade_groups=True).show()

#pc.jaccard_similarity_curve(pm,
#                            permutations=3,
#                            column="clade",
#                            group_order=["Modern", "AG1", "AG2", "AG3", "AG4", "AG5", "AG6", "AG7"],
#                            shade_groups=True).show()

#pvg.frequency_histogram(pm, log_y=True).show()
#pvg.exclusive_genes_bar(pm, column="clade").show()
pum.scatter3d(pm.compute_umap(), color_by="continent").show()
pum.scatter(pm, color_by="continent").show()