"""F4abcd blind check: cells per library and per donor behind the scored lineages.
Reads integrated atlas obs, lineage subset obs (h5py, X never read) and cell_scores_all.parquet.
"""
import sys
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

B = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1]); OUT.mkdir(exist_ok=True)
SC = B / "Analysis/SingleCell"


def obs_cols(path, cols):
    out = {}
    with h5py.File(path, "r") as h:
        o = h["obs"]
        for c in cols:
            g = o[c]
            if isinstance(g, h5py.Group):
                cats = g["categories"].asstr()[:]
                codes = g["codes"][:]
                out[c] = np.where(codes >= 0, cats[np.clip(codes, 0, None)], None)
            else:
                out[c] = g.asstr()[:] if g.dtype.kind in "OS" else g[:]
    return pd.DataFrame(out)


# 1. integrated atlas: cells per library x cell type
atlas = obs_cols(SC / "results_gpu_v2/integrated_atlas.h5ad", ["sample", "cell_type", "dataset"])
print("atlas cells", len(atlas))
at = atlas.groupby(["dataset", "sample", "cell_type"]).size().rename("n_atlas").reset_index()
at.to_csv(OUT / "atlas_cells_per_library_celltype.tsv", sep="\t", index=False)

# 2. lineage subsets used for Hotspot scoring
sub = {}
for ct in ["hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells"]:
    p = SC / f"results_gpu_v2/mcp/inputs/atlas_cnmf_{ct}.h5ad"
    d = obs_cols(p, ["sample"])
    sub[ct] = d.groupby("sample").size()
    print(ct, "subset cells", len(d), "libraries", d["sample"].nunique(), "min per library", sub[ct].min())

# 3. scored cells (module 1 rows) per library
pf = pq.ParquetFile(SC / "results_gpu_v2/hotspot_modules/cell_scores_all.parquet")
parts = []
for rg in range(pf.metadata.num_row_groups):
    t = pf.read_row_group(rg, columns=["cell_id", "module", "cell_type"]).to_pandas()
    t = t[t["module"].str.endswith("__1")]
    t["sample"] = t["cell_id"].str.split("_", n=1).str[0]
    parts.append(t.groupby(["cell_type", "sample"]).size())
sc = pd.concat(parts).groupby(level=[0, 1]).sum().rename("n_scored").reset_index()
sc.to_csv(OUT / "scored_cells_per_library.tsv", sep="\t", index=False)
print(sc.groupby("cell_type").n_scored.agg(["min", "median", "count"]))

# 4. atlas libraries with 1-29 cells of a core type: are they scored?
core = {"Hepatocytes": "hepatocytes", "Fibroblasts": "fibroblasts", "Macrophages": "macrophages",
        "Cholangiocytes": "cholangiocytes"}
a4 = at[at.cell_type.isin(core)].copy()
a4["lineage"] = a4.cell_type.map(core)
m = a4.merge(sc.rename(columns={"cell_type": "lineage"}), on=["lineage", "sample"], how="left")
m["n_scored"] = m.n_scored.fillna(0).astype(int)
summ = m.assign(lt30=m.n_atlas < 30).groupby(["lineage", "lt30"]).agg(
    libraries=("sample", "size"), libraries_scored=("n_scored", lambda x: int((x > 0).sum())),
    atlas_cells=("n_atlas", "sum"), scored_cells=("n_scored", "sum")).reset_index()
print(summ.to_string(index=False))
summ.to_csv(OUT / "lt30_library_scoring.tsv", sep="\t", index=False)
m.to_csv(OUT / "atlas_vs_scored_core_lineages.tsv", sep="\t", index=False)
print("T cells in atlas: labels", sorted(set(atlas.cell_type.dropna()) - set(core)))
