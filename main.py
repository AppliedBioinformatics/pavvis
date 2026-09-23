from pathlib import Path

from pavvis import PavMatrix
import pavvis.plots as plots

DATA_DIR = Path(__file__).parent / "tests" / "data"

pm = PavMatrix(DATA_DIR / "pav.csv", DATA_DIR / "metadata.csv")

print(f"Gene counts: {pm.metadata_columns}")

plots.presence_frequency_histogram(pm, n_bins=50).show()
plots.presence_per_sample_boxplot(pm, x="clade", gene_set="variable").show()
plots.presence_per_sample_scatter(pm, x="depth", gene_set="variable").show()