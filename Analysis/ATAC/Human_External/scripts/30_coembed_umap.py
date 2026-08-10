#!/usr/bin/env python
"""30_coembed_umap.py -- joint UMAP co-embedding of GSE244832 + GSE281367 snATAC.

Loads both processed h5ads (IDENTICAL 500bp tile space, 6,062,095 tiles), concatenates,
select_features -> spectral(LSI) -> harmony(batch=donor_id) -> UMAP, and carries the
existing per-cohort cell_type labels (both assigned by the same LIVER_MARKERS route).

Harmony over DONOR (30 donors; cohort nested within donor) is the standard scATAC batch
correction and the cell-level analog of the pooled DA's `cohort` mean-offset covariate.
Cohort mixing on the resulting UMAP is therefore an EMERGENT outcome, not forced -- an
honest test of whether the two cohorts occupy the same cell-state manifold after the same
class of batch correction the DA relies on. Also emits a quantitative mixing diagnostic
(per-cell-type cohort fraction) so "do they merge" is not judged by eye alone.

Output: Analysis/ATAC/Human_External/results/coembed/umap_export_joint.tsv.gz  (+ .h5ad)
Env: snapatac2.  bigmem (~315K cells x 6M tiles).
"""
import os, logging
import anndata as ad, numpy as np, pandas as pd
import snapatac2 as snap

ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
G244 = f"{ROOT}/Analysis/ATAC/Human_Multiome/results/snapatac2/snapatac2_processed.h5ad"
G281 = f"{ROOT}/Analysis/ATAC/Human_External/snapatac2_fast/snapatac2_processed.h5ad"
OUT  = f"{ROOT}/Analysis/ATAC/Human_External/results/coembed"
os.makedirs(OUT, exist_ok=True)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("coembed")

log.info("loading GSE244832 ..."); a = ad.read_h5ad(G244); a.obs["cohort"] = "GSE244832"
log.info("loading GSE281367 ..."); b = ad.read_h5ad(G281); b.obs["cohort"] = "GSE281367"
# unique donor ids across cohorts; drop per-cohort embeddings so concat is clean
a.obs["donor_id"] = ("244_" + a.obs["donor_id"].astype(str)).values
b.obs["donor_id"] = ("281_" + b.obs["donor_id"].astype(str)).values
for X in (a, b):
    for k in list(X.obsm.keys()): del X.obsm[k]
    for k in list(X.uns.keys()):  del X.uns[k]
keep = ["cell_type", "cohort", "condition", "donor_id"]
a.obs = a.obs[[c for c in keep if c in a.obs]]; b.obs = b.obs[[c for c in keep if c in b.obs]]

log.info("concatenating (%d + %d cells) ...", a.n_obs, b.n_obs)
comb = ad.concat([a, b], join="inner", index_unique="-")
del a, b
log.info("combined: %d cells x %d tiles", comb.n_obs, comb.n_vars)

snap.pp.select_features(comb, n_features=50000); log.info("features selected")
snap.tl.spectral(comb); log.info("spectral (LSI) done")

# harmony over donor (harmonypy, same engine as per-cohort processing)
import harmonypy
emb = comb.obsm["X_spectral"]
log.info("harmonypy on %d comps across %d donors ...", emb.shape[1], comb.obs.donor_id.nunique())
ho = harmonypy.run_harmony(emb, comb.obs, ["donor_id"], max_iter_harmony=20)
comb.obsm["X_spectral_harmony"] = ho.Z_corr.T
log.info("harmony done")

snap.tl.umap(comb, use_rep="X_spectral_harmony"); log.info("umap done")

df = comb.obs[["cell_type", "cohort", "condition", "donor_id"]].copy()
df["UMAP1"] = comb.obsm["X_umap"][:, 0]; df["UMAP2"] = comb.obsm["X_umap"][:, 1]
df.to_csv(f"{OUT}/umap_export_joint.tsv.gz", sep="\t", compression="gzip")
comb.write(f"{OUT}/coembed_processed.h5ad")
log.info("wrote %s/umap_export_joint.tsv.gz (%d cells)", OUT, len(df))

# ---- quantitative mixing diagnostic (not eyeball) ----
ct_coh = pd.crosstab(df.cell_type, df.cohort)
log.info("cell_type x cohort counts:\n%s", ct_coh.to_string())
# per cell type, cohort fraction vs the global cohort fraction (well-mixed => close)
glob = df.cohort.value_counts(normalize=True)
frac = ct_coh.div(ct_coh.sum(1), axis=0)
log.info("per-cell-type GSE281367 fraction (global=%.3f):\n%s",
         glob.get("GSE281367", np.nan), frac.get("GSE281367", pd.Series()).round(3).to_string())
# kBET-style: fraction of each cell's kNN from the same cohort (lower=better mixed).
# cheap proxy: neighbor-cohort purity on a 5k subsample in harmony space.
from sklearn.neighbors import NearestNeighbors
rng = np.random.default_rng(0)
sub = rng.choice(comb.n_obs, size=min(5000, comb.n_obs), replace=False)
Z = comb.obsm["X_spectral_harmony"][sub]; coh = df["cohort"].values[sub]
nn = NearestNeighbors(n_neighbors=30).fit(Z); _, idx = nn.kneighbors(Z)
same = np.mean([(coh[idx[i]] == coh[i]).mean() for i in range(len(sub))])
exp = float((glob**2).sum())  # expected same-cohort fraction if perfectly mixed
log.info("kNN same-cohort purity=%.3f (perfectly-mixed baseline=%.3f; 1.0=fully separated)",
         same, exp)
