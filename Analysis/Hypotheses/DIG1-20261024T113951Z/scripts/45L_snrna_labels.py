"""DIG1 job 45L (input for jobs 45 and 46). GSE256398 snRNA, 17 donors: V3 lineage labels through V0's
build_labelled (QC n_genes >= 200, mt < 5%; Leiden cluster argmax of marker scores; seed 0), so that jobs 45
and 46 use the same labels. V3 labels reproduce exactly only on Intel v4 nodes (RESULT_V0.md); this job
may run elsewhere, so it records how far the hepatocyte pseudobulk and the nuclei-per-lineage counts are
from V3's.
Outputs: hepatocyte nuclei counts (h5ad: X = raw counts, obs donor, var gene_ids), lineage per nucleus,
reproduction tables.
Usage: python -I 45L_snrna_labels.py <new_out_dir>
"""
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
V0 = ROOT / "Analysis/Hypotheses/V0-20261021T122533Z"
V3PB = ROOT / "Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
spec = importlib.util.spec_from_file_location("v0_audit", V0 / "scripts/v0_coverage_audit.py")
v0 = importlib.util.module_from_spec(spec); spec.loader.exec_module(v0)
v3 = v0.load_v3_module()
adata, donors = v0.build_labelled(v3, V0 / "downloads_untrusted")
C = adata.layers["counts"]
lab = adata.obs["lineage"].astype(str).to_numpy(); don = adata.obs["donor"].astype(str).to_numpy()
H = lab == "hepatocyte"
pb = pd.DataFrame({d: np.asarray(C[H & (don == d)].sum(axis=0)).ravel() for d in donors}, index=adata.var["gene_ids"].values)
old = pd.read_csv(V3PB / "pseudobulk_hepatocyte.tsv.gz", sep="\t", index_col=0).drop(columns="symbol")
same = old.shape == pb.shape and list(old.columns) == list(pb.columns)
ndiff = int((old.values != pb.loc[old.index].values).sum()) if same else -1
ct = pd.crosstab(adata.obs["donor"], adata.obs["lineage"])
ct.to_csv(OUT / "nuclei_per_donor_lineage.tsv", sep="\t")
v3ct = pd.read_csv(V3PB / "nuclei_per_donor_lineage.tsv", sep="\t", index_col=0).reindex(index=ct.index, columns=ct.columns).fillna(0)
rep = dict(same_shape=same, n_entries_differing=ndiff, reproduced=same and ndiff == 0,
           hep_nuclei_new=int(H.sum()), hep_nuclei_v3=int(v3ct["hepatocyte"].sum()),
           max_abs_diff_hep_per_donor=int((ct["hepatocyte"] - v3ct["hepatocyte"]).abs().max()),
           total_abs_diff_all_lineages=int((ct - v3ct).abs().to_numpy().sum()))
pd.DataFrame([rep]).to_csv(OUT / "v3_label_reproduction.tsv", sep="\t", index=False)
adata.obs[["donor", "leiden", "lineage"]].to_csv(OUT / "lineage_per_nucleus.tsv.gz", sep="\t")
hep = ad.AnnData(X=sp.csr_matrix(C[H]).astype(np.float32), obs=adata.obs.loc[H, ["donor", "lineage"]].copy(), var=adata.var[["gene_ids"]].copy())
hep.write_h5ad(OUT / "gse256398_hepatocyte_counts.h5ad", compression="gzip")
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
print(pd.DataFrame([rep]).to_string()); print(ct.to_string())
