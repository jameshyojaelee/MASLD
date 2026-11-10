"""Label check: mean counts per cell of LSEC / stellate / B markers by source refined type (all cells)."""
import sys
import numpy as np, pandas as pd
from scipy import io
SRC = sys.argv[1]; OUT = sys.argv[2]; GENES = sys.argv[3].split(",")
E = f"{SRC}/expression/68f4186f4e6da7a7e1b0edbe"
genes = [l.strip() for l in open(f"{E}/counts_features.tsv") if l.strip()]
bc = [l.strip() for l in open(f"{E}/counts_barcodes.tsv") if l.strip()]
meta = pd.read_csv(f"{SRC}/cluster/coords.txt", sep="\t", skiprows=[1], usecols=["NAME", "refined_cell_type", "condition"]).set_index("NAME").loc[bc]
X = io.mmread(f"{E}/counts_matrix.mtx").T.tocsr()
cols = [genes.index(g) for g in GENES]
sub = X[:, cols].toarray()
df = pd.DataFrame(sub, columns=GENES)
df["type"] = meta.refined_cell_type.values
df["condition"] = meta.condition.values
m = df.groupby(["type"])[GENES].mean().round(3)
m.to_csv(f"{OUT}/label_check_refined_means.tsv", sep="\t")
pd.set_option("display.width", 300); print(m.to_string())
