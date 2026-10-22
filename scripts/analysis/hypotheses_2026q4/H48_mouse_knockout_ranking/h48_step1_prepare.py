#!/usr/bin/env python3
"""H48 step 1: cell QC, animal verification, pseudobulk, gene universe, signatures.

Writes into OUT/step1/. No knockout score is computed here.
"""
import argparse
import glob
import json
import os

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

SEED = 48
HEP_MARKERS = ["Ttr", "Apoh", "Fgb", "Fga", "Hpx", "Tf"]
FAST_MARKERS = ["Hmgcs2", "Cyp4a14", "Angptl4", "G6pc", "Ppargc1a", "Fgf21"]
N_SIG = 200

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True)
ap.add_argument("--orth", required=True)
ap.add_argument("--sig_sc", required=True)
ap.add_argument("--sig_bulk", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--h5dir", required=True)
a = ap.parse_args()
out = os.path.join(a.out, "step1")
os.makedirs(out, exist_ok=False)
rng = np.random.default_rng(SEED)
log = {}

adata = ad.read_h5ad(a.h5ad)
X = adata.X.tocsr().astype(np.float64)
obs = adata.obs.copy()
var = adata.var.copy()
var["gene_id"] = var["gene_id"].astype(str)
var["gene_name"] = var["gene_name"].astype(str)
log["cells_in_h5ad"] = int(X.shape[0])
log["genes_in_h5ad"] = int(X.shape[1])

# --- cell QC
tot = np.asarray(X.sum(1)).ravel()
mk_idx = [np.where(var["gene_name"].values == g)[0] for g in HEP_MARKERS]
mk_idx = np.concatenate([i for i in mk_idx if len(i)])
mk = np.asarray(X[:, mk_idx].sum(1)).ravel()
frac = np.divide(mk, tot, out=np.zeros_like(mk), where=tot > 0)
# sublib names: concat index i -> i-th sorted GEO h5 file (Amendment 1)
h5_files = sorted(glob.glob(os.path.join(a.h5dir, "GSM8478327_*.h5")))
names = [os.path.basename(f).split("_sample")[0].replace("GSM8478327_", "") for f in h5_files]
assert len(names) == 32
obs["sublib"] = [names[int(i)] for i in obs["sublib"].astype(str)]
bc24 = obs["cell_barcode_24"].astype(str).str[:24]
probe = {"ACTTTAGG": "A", "AACGGGAA": "B", "AGTAGGCT": "C", "ATGTTGAC": "D"}
obs["probe_letter"] = bc24.str[16:].map(probe)
assert (obs["probe_letter"].values == obs["sublib"].str[-1].values).all(), "probe letter != sublib letter"
from collections import Counter
cnt = Counter()
for f in h5_files:
    with h5py.File(f, "r") as h:
        cnt.update(x.decode().split("-")[0] for x in h["matrix/barcodes"][:])
multi = bc24.map(lambda b: cnt[b] > 1).values
no_feature = obs["sublib"].isin(["S7_A", "S7_C", "S8_A", "S8_C"]).values
log["cells_in_sublibs_without_feature_library"] = int(no_feature.sum())
log["cells_with_barcode_in_multiple_h5"] = int(multi.sum())
log["cells_dropped_amendment1"] = int((no_feature | multi).sum())
keep = (obs["tier"].astype(str).values == "confident") & (tot >= 500) & ~no_feature & ~multi
med = np.median(frac[keep])
keep_hep = keep & (frac >= 0.1 * med)
log["cells_confident"] = int((obs["tier"].astype(str).values == "confident").sum())
log["cells_umi_ge500"] = int(keep.sum())
log["hep_marker_median_fraction"] = float(med)
log["cells_after_hepatocyte_filter"] = int(keep_hep.sum())
X = X[keep_hep]
obs = obs.loc[keep_hep].copy()
tot = tot[keep_hep]
obs["sublib"] = obs["sublib"].astype(str)
obs["letter"] = obs["sublib"].str.split("_").str[1]
obs["well"] = obs["sublib"].str.split("_").str[0]
obs["target_gene"] = obs["target_gene"].astype(str)
obs["sgrna_id"] = obs["sgrna_id"].astype(str)

# --- animal verification: pseudobulk per sublib
subs = sorted(obs["sublib"].unique())
G = sp.csr_matrix((np.ones(len(obs)), (pd.Categorical(obs["sublib"], categories=subs).codes, np.arange(len(obs)))),
                  shape=(len(subs), len(obs)))
pb = np.asarray((G @ X).todense())
ncell_sub = np.asarray(G.sum(1)).ravel()
cpm = pb / pb.sum(1, keepdims=True) * 1e6
lc = np.log2(cpm + 1)
expressed = (cpm.mean(0) >= 10)
v = lc[:, expressed].var(0)
top = np.argsort(v)[::-1][:2000]
M = lc[:, expressed][:, top]
C = np.corrcoef(M)
D = 1 - C
np.fill_diagonal(D, 0)
Z = linkage(squareform(D, checks=False), method="average")
cl = fcluster(Z, t=4, criterion="maxclust")
ver = pd.DataFrame({"sublib": subs, "letter": [s.split("_")[1] for s in subs], "n_cells": ncell_sub.astype(int), "cluster": cl})
# fasted marker z
fz = []
for g in FAST_MARKERS:
    j = np.where(var["gene_name"].values == g)[0]
    if len(j):
        x = lc[:, j[0]]
        fz.append((x - x.mean()) / (x.std() + 1e-9))
ver["fasting_z"] = np.mean(fz, 0)
xj = np.where(var["gene_name"].values == "Xist")[0]
ver["Xist_log2cpm"] = lc[:, xj[0]] if len(xj) else np.nan
ver.to_csv(os.path.join(out, "animal_verification_sublibs.tsv"), sep="\t", index=False)
pd.DataFrame(C, index=subs, columns=subs).to_csv(os.path.join(out, "sublib_correlation.tsv"), sep="\t")
# rule: weight by sublib (tiny sublibs with < 200 cells are reported but not counted)
big = ver[ver.n_cells >= 200]
purity = big.groupby("cluster")["letter"].agg(lambda s: s.value_counts().iloc[0] / len(s))
major = big.groupby("cluster")["letter"].agg(lambda s: s.value_counts().index[0])
letters_ok = (purity >= 0.9).all() and major.nunique() == 4 and len(major) == 4
log["animal_rule_cluster_purity"] = purity.to_dict()
log["animal_rule_cluster_major_letter"] = major.to_dict()
log["animal_rule_pass"] = bool(letters_ok)
fl = ver.groupby("letter")["fasting_z"].mean()
log["fasting_z_by_letter"] = fl.to_dict()
log["fasted_letter"] = str(fl.idxmax())
unit_col = "letter" if letters_ok else "sgrna_id"
log["unit"] = "animal (probe-barcode letter)" if letters_ok else "sgRNA (fallback)"

# --- orthologs (one2one, unique both ways)
o = pd.read_csv(a.orth, sep="\t")
o = o[o["ortholog_type"] == "ortholog_one2one"].dropna(subset=["human_gene_ensembl", "mouse_gene_ensembl"])
o = o[~o["mouse_gene_ensembl"].duplicated(keep=False) & ~o["human_gene_ensembl"].duplicated(keep=False)]
m2h = dict(zip(o["mouse_gene_ensembl"], o["human_gene_ensembl"]))
m2hsym = dict(zip(o["mouse_gene_ensembl"], o["human_gene_symbol"]))

# --- gene universe: one2one + CPM >= 5 in pooled NC
is_nc = obs["target_gene"].values == "control"
nc_pool = np.asarray(X[is_nc].sum(0)).ravel()
nc_cpm = nc_pool / nc_pool.sum() * 1e6
has_orth = var["gene_id"].map(m2h).notna().values
univ = has_orth & (nc_cpm >= 5)
log["genes_one2one_in_data"] = int(has_orth.sum())
log["genes_universe"] = int(univ.sum())
log["nc_cells_total"] = int(is_nc.sum())
uv = var.loc[univ, ["gene_id", "gene_name"]].copy()
uv["human_ensembl"] = uv["gene_id"].map(m2h)
uv["human_symbol"] = uv["gene_id"].map(m2hsym)
uv["nc_log_cpm"] = np.log2(nc_cpm[univ])
uv["expr_bin"] = pd.qcut(uv["nc_log_cpm"].rank(method="first"), 20, labels=False)
uv.to_csv(os.path.join(out, "gene_universe.tsv"), sep="\t", index=False)
Xu = X[:, np.where(univ)[0]].tocsr()

# --- samples: (sgRNA, letter) pseudobulks (Amendment 2 fallback)
if letters_ok:
    raise SystemExit("animal rule passed; v2 implements only the prespecified sgRNA fallback")
obs["unit"] = obs["sgrna_id"].astype(str) + "|" + obs["letter"].astype(str)
key = obs["target_gene"] + "|" + obs["unit"]
keys = sorted(key.unique())
codes = pd.Categorical(key, categories=keys).codes
Gk = sp.csr_matrix((np.ones(len(obs)), (codes, np.arange(len(obs)))), shape=(len(keys), len(obs)))
PB = (Gk @ Xu).T.tocsc()  # genes x keys
cells = pd.DataFrame({"key": keys, "target_gene": [k.split("|")[0] for k in keys],
                      "sgrna_id": [k.split("|")[1] for k in keys], "letter": [k.split("|")[2] for k in keys],
                      "n_cells": np.asarray(Gk.sum(1)).ravel().astype(int)})
cells["included"] = cells.n_cells >= 10
scipy.io.mmwrite(os.path.join(out, "pseudobulk_counts.mtx"), sp.csr_matrix(PB).astype(np.int64))
cells.to_csv(os.path.join(out, "pseudobulk_columns.tsv"), sep="\t", index=False)

# NC per-cell counts for null 2 (only cells in included NC samples)
nc_key = key[is_nc]
inc_keys = set(cells.loc[cells.included, "key"])
nc_sel = np.where(is_nc)[0][nc_key.isin(inc_keys).values]
nc_obs = obs.iloc[nc_sel][["sgrna_id", "letter", "sublib"]].copy()
nc_obs["key"] = key.iloc[nc_sel].values
scipy.io.mmwrite(os.path.join(out, "nc_cells_counts.mtx"), sp.csr_matrix(Xu[nc_sel].T).astype(np.int64))
nc_obs.to_csv(os.path.join(out, "nc_cells_meta.tsv"), sep="\t")

# per-knockout eligibility: >= 50 cells total and >= 2 guides with >= 10 cells (over included samples)
rows = []
for g, d in cells[cells.target_gene != "control"].groupby("target_gene"):
    di = d[d.included]
    gcells = di.groupby("sgrna_id").n_cells.sum()
    rows.append(dict(target_gene=g, n_cells_total=int(d.n_cells.sum()), n_cells_in_included_samples=int(di.n_cells.sum()),
                     n_guides_ge10=int((gcells >= 10).sum()), n_samples=int(len(di)),
                     letters=",".join(sorted(di.letter.unique())),
                     cells_per_guide_letter=";".join(f"{s}|{l}:{n}" for s, l, n in zip(d.sgrna_id, d.letter, d.n_cells)),
                     tested=bool(d.n_cells.sum() >= 50 and (gcells >= 10).sum() >= 2)))
elig = pd.DataFrame(rows)
elig.to_csv(os.path.join(out, "knockout_eligibility.tsv"), sep="\t", index=False)
nc = cells[cells.target_gene == "control"]
log["n_targets_in_data"] = int(len(elig))
log["n_targets_tested"] = int(elig.tested.sum())
log["targets_dropped"] = elig.loc[~elig.tested, ["target_gene", "n_cells_total", "n_guides_ge10"]].values.tolist()
log["nc_guides"] = int(nc.sgrna_id.nunique())
log["nc_samples_included"] = int(nc.included.sum())
log["nc_cells_in_included_samples"] = int(nc.loc[nc.included, "n_cells"].sum())
log["cells_per_letter"] = obs.letter.value_counts().to_dict()



# --- signatures
def build_sig(path, label, gene_col, t_col):
    s = pd.read_csv(path)
    s["ens"] = s[gene_col].astype(str).str.replace(r"\.\d+$", "", regex=True)
    s = s.dropna(subset=[t_col])
    s = s[~s["ens"].duplicated(keep=False)]
    n_in = len(s)
    h2row = dict(zip(uv["human_ensembl"], uv["gene_id"]))
    sym2row = dict(zip(uv["human_symbol"], uv["gene_id"]))
    is_ens = s["ens"].str.startswith("ENSG")
    s["mouse_gene_id"] = np.where(is_ens, s["ens"].map(h2row), s["ens"].map(sym2row))
    s.loc[s["mouse_gene_id"].isna() | (s["mouse_gene_id"] == "nan"), "mouse_gene_id"] = np.nan
    s = s[~(s["mouse_gene_id"].notna() & s["mouse_gene_id"].duplicated(keep=False))]
    s_m = s.dropna(subset=["mouse_gene_id"]).copy()
    up = s_m.sort_values(t_col, ascending=False).head(N_SIG)
    up = up[up[t_col] > 0]
    dn = s_m.sort_values(t_col, ascending=True).head(N_SIG)
    dn = dn[dn[t_col] < 0]
    sig = pd.concat([up, dn])
    sig = sig[["ens", "mouse_gene_id", t_col, "logFC", "padj" if "padj" in sig else t_col]].rename(columns={t_col: "human_t"})
    sig["direction"] = np.where(sig.human_t > 0, "up_in_MASLD", "down_in_MASLD")
    sig["mouse_symbol"] = sig["mouse_gene_id"].map(dict(zip(uv.gene_id, uv.gene_name)))
    sig.to_csv(os.path.join(out, f"signature_{label}.tsv"), sep="\t", index=False)
    return dict(n_human_genes=n_in, n_mapped_to_universe=int(len(s_m)), n_sig_up=int(len(up)), n_sig_down=int(len(dn)),
                n_padj_lt_005_in_human_table=int((s["padj"] < 0.05).sum()) if "padj" in s else None,
                min_abs_t_in_signature=float(sig.human_t.abs().min()))


log["signature_S1_sc"] = build_sig(a.sig_sc, "S1_sc_hepatocyte", "gene", "t_stat")
log["signature_S2_bulk_provisional"] = build_sig(a.sig_bulk, "S2_bulk_provisional", "gene", "t")
with open(os.path.join(out, "step1_log.json"), "w") as f:
    json.dump(log, f, indent=1, default=str)
print(json.dumps(log, indent=1, default=str))
