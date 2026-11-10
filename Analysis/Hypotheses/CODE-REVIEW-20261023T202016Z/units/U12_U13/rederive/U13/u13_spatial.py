#!/usr/bin/env python3
"""U13 independent re-derivation of the Figure 5F spatial matched-gene null.

Written from docs/manuscript/METHODS.md, "Spatial transcriptomic analyses":
  * spot score = z-scored log-normalized expression of measured fixed member genes,
    weighted by their L1 weights rescaled to sum to one;
  * residualized on intercept + 16 cell2location q05 abundances + log1p library size
    + log1p detected genes (Figure 5F model), or on hepatocyte q05 only (prespecified model);
  * residual Moran's I on six-neighbour spot graphs within connected tissue islands;
    section value = island-size-weighted mean; donor = mean of its sections; source = mean of donors
    (Vu arrays are their own units);
  * 9,999 matched control-gene sets (expression, detection, biotype, MT/ribo, correlation with
    lineage abundance; bins relaxed by up to two when < 10 candidates); p = (k+1)/(B+1); BH.

The matched-set sampler keeps the producer's draw order (one Generator per source, programs in
registry order, member genes in symbol order, one bounded-integer draw per gene, no repeats inside
a set) so that a same-seed run can be checked against the recorded matched-set SHA256. Pool binning
uses the same float32 sparse operations as the producer for the same reason (tie order at decile
edges depends on float rounding). Everything downstream of the sets is written independently:
residualization by an orthonormal basis of the design applied to every gene once, Moran's I, donor
collapse, p, BH, sign-flip.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
H5 = {
    "GSE192741": P / "Analysis/Spatial/results/cell2location/spatial_model/spatial_deconvolved.h5ad",
    "Vu_et_al_2025": P / "Analysis/Spatial/results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad",
}
PRODUCER_SEED = {"GSE192741": 42, "Vu_et_al_2025": 1042}
HOT = P / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
GENEMETA = P / "data/gencode_v49_gene_metadata.tsv.gz"
R2 = P / ("Analysis/Multimodal_Program_Projection/candidates/"
          "spatial-full-composition-sensitivity-candidate-2026-08-11-r2/full_composition_program_results.tsv")
AN4 = P / "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN4/20261021T124414Z"
AN4_FAMILY_TAG = {"GSE192741": "GSE192741", "Vu_et_al_2025": "Vu"}  # AN4 run-directory names
LINEAGE = {"hepatocytes": "Hepatocytes", "fibroblasts": "Fibroblasts", "macrophages": "Macrophages",
           "cholangiocytes": "Cholangiocytes", "tcells": "T cells"}
STROMAL = "hotspot_hepatocytes_f05c535ae5bbc0b9"
DUCT = "hotspot_hepatocytes_48f39dd4d817a10e"
CONTENT = {STROMAL: "fibroblasts", DUCT: "cholangiocytes"}
K_NEIGHBOURS, EDGE_MULT, MIN_ISLAND = 6, 2.5, 8
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 23), b""):
            h.update(chunk)
    return h.hexdigest()


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p, kind="mergesort")
    adj = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.minimum(adj, 1.0)
    return q


def signflip_one_sided(e):
    e = np.asarray(e, float)
    t = e.sum()
    signs = np.array(list(itertools.product([1.0, -1.0], repeat=len(e))))
    return float(np.mean(signs @ np.abs(e) >= t - 1e-12))


# ---------------------------------------------------------------- data
class Source:
    def __init__(self, dataset):
        import anndata as ad
        log(f"loading {dataset}")
        a = ad.read_h5ad(H5[dataset])
        self.dataset = dataset
        self.genes = pd.Index(a.var_names.astype(str))
        assert self.genes.is_unique
        obs = a.obs.copy()
        obs["sample_id"] = obs["sample_id"].astype(str)
        obs["individual"] = obs["individual"].astype(str)
        self.obs = obs
        x = a.layers["counts"]
        counts = x.tocsr().astype(np.float32) if sparse.issparse(x) else sparse.csr_matrix(np.asarray(x, np.float32))
        factors = [str(f).replace("means_per_cluster_mu_fg_", "") for f in a.uns["mod"]["factor_names"]]
        assert len(factors) == 16 and len(set(factors)) == 16
        self.factors = factors
        self.q05 = np.asarray(a.obsm["q05_cell_abundance_w_sf"], float)
        self.coords = np.asarray(a.obsm["spatial"], float)
        del a
        lib = pd.to_numeric(obs["total_counts"]).to_numpy(float)
        ngenes = pd.to_numeric(obs["n_genes_by_counts"]).to_numpy(float)
        scale = np.divide(1e4, lib, out=np.zeros_like(lib), where=lib > 0).astype(np.float32)
        norm = (sparse.diags(scale) @ counts).tocsr()
        norm.data = np.log1p(norm.data)
        self.norm = norm
        self.detection = np.asarray((counts > 0).mean(axis=0)).ravel()
        del counts
        self.n_spots, self.n_genes = norm.shape
        # pool properties
        gene_mean = np.asarray(norm.mean(axis=0)).ravel()
        meta = pd.read_csv(GENEMETA, sep="\t", usecols=["gene_name", "gene_biotype"]).drop_duplicates("gene_name")
        bt = dict(zip(meta["gene_name"].astype(str), meta["gene_biotype"].astype(str)))
        biot = pd.Series([bt.get(g, "unknown") for g in self.genes])
        self.bt_code = pd.factorize(biot)[0]
        self.mt = np.asarray(self.genes.str.startswith("MT-"), dtype=bool)
        self.ribo = np.asarray(self.genes.str.match(r"^RP[SL]"), dtype=bool)
        self.mean_bin = self.decile(gene_mean)
        self.det_bin = self.decile(self.detection)
        self.lin_bin, self.lin_corr, self.bin_check = {}, {}, {}
        for lin, label in LINEAGE.items():
            ab = self.q05[:, factors.index(label)]
            corr = self.corr_like_producer(ab, gene_mean)
            self.lin_corr[lin] = corr
            self.lin_bin[lin] = self.decile(corr)
        # independent float64 Pearson check of the hepatocyte correlation bins
        hep = self.q05[:, factors.index("Hepatocytes")]
        dense_corr = self.corr_dense(hep)
        b64 = self.decile(dense_corr)
        self.bin_check["hepatocytes_bins_differ_float64_vs_float32"] = int(np.sum(b64 != self.lin_bin["hepatocytes"]))
        self.bin_check["hepatocytes_corr_maxabs_diff"] = float(np.max(np.abs(dense_corr - self.lin_corr["hepatocytes"])))
        # designs
        canon = sorted(factors)
        ab16 = self.q05[:, [factors.index(f) for f in canon]]
        one = np.ones(self.n_spots)
        self.design = {
            "all16": np.column_stack([one, ab16, np.log1p(lib), np.log1p(ngenes)]),
            "hep_only": np.column_stack([one, hep, np.log1p(lib), np.log1p(ngenes)]),
        }
        # graphs
        self.build_graphs()

    @staticmethod
    def decile(values):
        r = pd.Series(np.asarray(values)).rank(method="first")
        return pd.qcut(r, 10, labels=False, duplicates="drop").to_numpy(int)

    def corr_like_producer(self, abundance, gene_mean):
        a = np.asarray(abundance, float)
        a = a - np.nanmean(a)
        a[~np.isfinite(a)] = 0.0
        ass = float(np.sum(a * a))
        sq = np.asarray(self.norm.multiply(self.norm).mean(axis=0)).ravel()
        gss = self.n_spots * np.maximum(sq - gene_mean * gene_mean, 0.0)
        num = np.asarray(self.norm.T @ a).ravel()
        den = np.sqrt(gss * ass)
        return np.divide(num, den, out=np.zeros_like(num), where=den > 0)

    def corr_dense(self, abundance):
        a = abundance - abundance.mean()
        out = np.empty(self.n_genes)
        csc = self.norm.tocsc()
        for s in range(0, self.n_genes, 2000):
            x = csc[:, s:s + 2000].toarray().astype(np.float64)
            x -= x.mean(0)
            den = np.sqrt((x * x).sum(0) * (a * a).sum())
            out[s:s + 2000] = np.divide(x.T @ a, den, out=np.zeros(x.shape[1]), where=den > 0)
        return out

    def build_graphs(self):
        sid = self.obs["sample_id"].to_numpy()
        self.sections, self.islands, audit = [], [], []
        for s in pd.unique(sid):
            idx = np.flatnonzero(sid == s)
            xy = self.coords[idx]
            tree = cKDTree(xy)
            d, nb = tree.query(xy, k=min(K_NEIGHBOURS + 1, len(idx)))
            rows = np.repeat(np.arange(len(idx)), d.shape[1])
            d, nb = d.ravel(), nb.ravel()
            keep = nb != rows
            rows, d, nb = rows[keep], d[keep], nb[keep]
            # keep the k nearest non-self neighbours of each spot
            order = np.lexsort((d, rows))
            rows, d, nb = rows[order], d[order], nb[order]
            rank_in_row = np.arange(len(rows)) - np.searchsorted(rows, rows)
            sel = rank_in_row < K_NEIGHBOURS
            rows, d, nb = rows[sel], d[sel], nb[sel]
            first = np.searchsorted(rows, np.arange(len(idx)))
            nearest = d[first]
            median_nn = float(np.median(nearest[nearest > 0]))
            ok = d <= EDGE_MULT * median_nn
            pairs = np.unique(np.sort(np.column_stack([rows[ok], nb[ok]]), axis=1), axis=0)
            adj = sparse.coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(idx), len(idx)))
            ncomp, lab = connected_components(adj, directed=False)
            sec_islands = []
            for c in range(ncomp):
                mem = np.flatnonzero(lab == c)
                if len(mem) < MIN_ISLAND:
                    continue
                emask = lab[pairs[:, 0]] == c
                if not emask.any():
                    continue
                local = np.full(len(idx), -1)
                local[mem] = np.arange(len(mem))
                sec_islands.append((idx[mem], local[pairs[emask, 0]], local[pairs[emask, 1]]))
            self.sections.append(s)
            self.islands.append(sec_islands)
            audit.append({"sample_id": s, "n_spots": len(idx), "n_islands": len(sec_islands),
                          "n_spots_in_graph": int(sum(len(i[0]) for i in sec_islands)),
                          "n_edges": int(sum(len(i[1]) for i in sec_islands)), "median_nn": median_nn})
        self.graph_audit = pd.DataFrame(audit)
        ind = self.obs[["sample_id", "individual"]].drop_duplicates().set_index("sample_id")["individual"]
        self.units = list(pd.unique(ind.loc[self.sections]))
        D = np.zeros((len(self.sections), len(self.units)))
        for i, s in enumerate(self.sections):
            D[i, self.units.index(ind.loc[s])] = 1.0
        D /= D.sum(axis=0, keepdims=True)
        self.D = D  # section -> donor/array mean
        self.collapse_vec = D.mean(axis=1)  # donor mean then mean of donors

    def residual_matrix(self, design_name):
        X = self.design[design_name]
        U, s, _ = np.linalg.svd(X, full_matrices=False)
        Q = U[:, s > s.max() * max(X.shape) * np.finfo(float).eps]
        RT = np.empty((self.n_genes, self.n_spots), np.float32)
        csc = self.norm.tocsc()
        for st in range(0, self.n_genes, 1000):
            z = csc[:, st:st + 1000].toarray().astype(np.float64)
            sd = z.std(0, ddof=1)
            sd[~(sd > 0)] = 1.0
            z = (z - z.mean(0)) / sd
            z -= Q @ (Q.T @ z)
            RT[st:st + 1000] = z.T
        return RT, Q.shape[1]


# ---------------------------------------------------------------- Moran
def moran_sections(S, islands):
    """S: (b, n_spots). Returns (b, n_sections) island-size-weighted residual Moran's I."""
    out = np.empty((S.shape[0], len(islands)))
    for k, sec in enumerate(islands):
        vals, wts = [], []
        for idx, left, right in sec:
            z = S[:, idx].astype(np.float64)
            z -= z.mean(axis=1, keepdims=True)
            num = np.einsum("ij,ij->i", z[:, left], z[:, right])
            den = np.einsum("ij,ij->i", z, z)
            vals.append(len(idx) * num / (len(left) * den))
            wts.append(len(idx))
        out[:, k] = np.average(np.vstack(vals), axis=0, weights=np.asarray(wts, float))
    return out


G = {}  # worker globals (inherited by fork)


def score_block(RT, sets, weights):
    acc = np.zeros((sets.shape[0], RT.shape[1]), np.float32)
    for j, w in enumerate(weights):
        acc += np.float32(w) * RT[sets[:, j]]
    return acc


def null_task(args):
    key, design, sets, weights, n_sens, top = args
    RT, islands = G["RT"][design], G["islands"]
    B = sets.shape[0]
    secs = np.empty((B, len(islands)))
    for st in range(0, B, 256):
        secs[st:st + 256] = moran_sections(score_block(RT, sets[st:st + 256], weights), islands)
    out = {"sections": secs}
    if n_sens:
        sub = sets[:n_sens]
        eq = np.full(len(weights), 1.0 / len(weights))
        lw = np.delete(weights, top)
        lw = lw / lw.sum()
        out["equal"] = moran_sections(score_block(RT, sub, eq), islands)
        out["leave_top"] = moran_sections(score_block(RT, np.delete(sub, top, axis=1), lw), islands)
    return key, out


# ---------------------------------------------------------------- programs and matching
def load_programs():
    reg = pd.read_csv(HOT / "program_registry_v2.tsv", sep="\t")
    reg["module"] = reg["module"].astype(int)
    reg = reg.sort_values(["cell_type", "module", "program_uid"]).reset_index(drop=True)
    mem = pd.read_csv(HOT / "program_membership_v2.tsv", sep="\t")
    mem = mem[mem["mapped_symbol_status"] == "gencode_v49_unique_symbol_confirmed"].copy()
    mem["gene_symbol"] = mem["mapped_symbol"].astype(str).str.strip()
    return reg, mem


def measured(src, mem, pid):
    m = mem[mem["program_uid"] == pid].groupby("gene_symbol", as_index=False)["original_l1_weight"].sum()
    gi = src.genes.get_indexer(m["gene_symbol"])
    m["gi"] = gi
    m = m[(gi >= 0)].copy()
    m = m[src.detection[m["gi"].to_numpy()] >= 0.01].copy()
    retained = float(m["original_l1_weight"].sum())
    m["w"] = m["original_l1_weight"] / retained if retained > 0 else np.nan
    return m, retained


def candidates(src, gi, lin, excl):
    lb = src.lin_bin[lin]
    base = (src.bt_code == src.bt_code[gi]) & (src.mt == src.mt[gi]) & (src.ribo == src.ribo[gi]) & ~excl
    dm, dd, dl = np.abs(src.mean_bin - src.mean_bin[gi]), np.abs(src.det_bin - src.det_bin[gi]), np.abs(lb - lb[gi])
    c = np.array([], int)
    for lr in (0, 1, 2):
        for er in (0, 1, 2):
            c = np.flatnonzero(base & (dm <= er) & (dd <= er) & (dl <= lr))
            if len(c) >= 10:
                return c, er, lr
    if len(c):
        return c, 2, 2
    c = np.flatnonzero((src.bt_code == src.bt_code[gi]) & (dl <= 2) & ~excl)
    if len(c):
        return c, 3, 2
    c = np.flatnonzero((src.bt_code == src.bt_code[gi]) & ~excl)
    if not len(c):
        raise RuntimeError(f"no controls for gene index {gi}")
    return c, 3, 3


def draw_sets(cands, n_null, rng, n_genes):
    sets = np.empty((n_null, len(cands)), np.int32)
    taken = np.zeros(n_genes, bool)
    for i in range(n_null):
        picked = []
        for j, c in enumerate(cands):
            if picked:
                avail = c[~taken[c]]
                srcset = avail if len(avail) else c
            else:
                srcset = c
            g = int(srcset[int(rng.integers(0, len(srcset)))])
            sets[i, j] = g
            taken[g] = True
            picked.append(g)
        taken[picked] = False
    return sets


# ---------------------------------------------------------------- main
def summarize(src, obs_sec, null_sec):
    """obs_sec: (n_sections,), null_sec: (B, n_sections)."""
    obs_c = float(obs_sec @ src.collapse_vec)
    null_c = null_sec @ src.collapse_vec
    k = int(np.sum(null_c >= obs_c))
    B = len(null_c)
    mu, sd = float(null_c.mean()), float(null_c.std(ddof=1))
    unit_obs = obs_sec @ src.D
    unit_null_mean = (null_sec @ src.D).mean(axis=0)
    return {"observed": obs_c, "null_mean": mu, "null_sd": sd, "excess": obs_c - mu, "k": k, "B": B,
            "p": (k + 1) / (B + 1), "z": (obs_c - mu) / sd if sd > 0 else np.nan,
            "null_q025": float(np.quantile(null_c, 0.025)), "null_q975": float(np.quantile(null_c, 0.975)),
            "unit_excess": unit_obs - unit_null_mean, "unit_observed": unit_obs}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=list(H5))
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--order", choices=["registry", "reversed"], default="registry")
    ap.add_argument("--n-null", type=int, default=9999)
    ap.add_argument("--n-sens", type=int, default=999)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--skip-family", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    meta = {"args": {k: str(v) for k, v in vars(args).items()}, "numpy": np.__version__,
            "inputs_sha256": {"h5ad": sha256_file(H5[args.dataset]),
                              "registry": sha256_file(HOT / "program_registry_v2.tsv"),
                              "membership": sha256_file(HOT / "program_membership_v2.tsv"),
                              "gene_meta": sha256_file(GENEMETA)}}
    src = Source(args.dataset)
    meta["bin_check"] = src.bin_check
    meta["n_spots"], meta["n_genes"] = src.n_spots, src.n_genes
    src.graph_audit.to_csv(args.out / "graph_audit.tsv", sep="\t", index=False)
    log(f"spots {src.n_spots} genes {src.n_genes} sections {len(src.sections)} units {len(src.units)} "
        f"bin_check {src.bin_check}")
    RT = {}
    for d in ("all16", "hep_only"):
        RT[d], rank = src.residual_matrix(d)
        meta[f"design_rank_{d}"] = rank
    log(f"residual matrices built ranks {meta['design_rank_all16']}/{meta['design_rank_hep_only']}")
    G["RT"], G["islands"] = RT, src.islands
    del src.norm
    reg, mem = load_programs()
    r2 = pd.read_csv(R2, sep="\t")
    rec_focal = {(r.program_id, "hepatocytes"): r.matched_set_sha256 for r in r2[r2.dataset == args.dataset].itertuples()}
    fl = pd.read_csv(AN4 / "focal_lineage/focal_lineage_results.tsv", sep="\t")
    for r in fl[(fl.dataset == args.dataset) & (fl.scheme == "content_lineage")].itertuples():
        rec_focal[(r.program_id, "content")] = r.matched_set_sha256
    pool = mp.get_context("fork").Pool(args.workers)

    # ---- focal two-program runs (Figure 5F primary and content-lineage sensitivity)
    focal_ids = [p for p in reg["program_uid"] if p in (STROMAL, DUCT)]
    if args.order == "reversed":
        focal_ids = focal_ids[::-1]
    union2 = set(mem.loc[mem.program_uid.isin([STROMAL, DUCT]), "gene_symbol"])
    excl2 = src.genes.isin(list(union2))
    focal_rows, unit_rows = [], []
    for scheme in ("hepatocytes", "content"):
        rng = np.random.default_rng(args.seed)
        tasks, info = [], {}
        for pid in focal_ids:
            m, retained = measured(src, mem, pid)
            lin = "hepatocytes" if scheme == "hepatocytes" else CONTENT[pid]
            cands = [candidates(src, gi, lin, excl2) for gi in m["gi"]]
            sets = draw_sets([c[0] for c in cands], args.n_null, rng, src.n_genes)
            sha = hashlib.sha256(sets.tobytes()).hexdigest()
            w = m["w"].to_numpy(float)
            top = int(np.argmax(w))
            info[pid] = (m, retained, w, top, sha, sum(c[1] > 0 for c in cands), sum(c[2] > 0 for c in cands))
            designs = ("all16", "hep_only") if scheme == "hepatocytes" else ("all16",)
            for d in designs:
                tasks.append(((pid, d), d, sets, w, args.n_sens if d == "all16" else 0, top))
        res = dict(pool.map(null_task, tasks))
        for (pid, d), out in res.items():
            m, retained, w, top, sha, n_er, n_lr = info[pid]
            idx = m["gi"].to_numpy()
            obs_sec = moran_sections((w @ RT[d][idx])[None, :].astype(np.float32), src.islands)[0]
            s = summarize(src, obs_sec, out["sections"])
            row = {"dataset": args.dataset, "seed": args.seed, "order": args.order, "scheme": scheme, "design": d,
                   "program_id": pid, "n_measured": len(m), "retained_l1": retained,
                   "n_genes_expression_relaxed": n_er, "n_genes_lineage_relaxed": n_lr,
                   **{k: v for k, v in s.items() if not k.startswith("unit_")},
                   "matched_set_sha256": sha, "recorded_sha256": rec_focal.get((pid, scheme), ""),
                   "sha_match": sha == rec_focal.get((pid, scheme), "")}
            if "equal" in out:
                for nm, ww, ii in (("equal", np.full(len(w), 1 / len(w)), idx),
                                   ("leave_top", np.delete(w, top) / np.delete(w, top).sum(), np.delete(idx, top))):
                    o = float(moran_sections((ww @ RT[d][ii])[None, :].astype(np.float32), src.islands)[0] @ src.collapse_vec)
                    row[f"{nm}_observed"] = o
                    row[f"{nm}_null_mean"] = float((out[nm] @ src.collapse_vec).mean())
                row["sensitivity_sign_agree"] = bool(np.sign(s["z"]) == np.sign(row["equal_observed"] - row["equal_null_mean"])
                                                     == np.sign(row["leave_top_observed"] - row["leave_top_null_mean"]))
            focal_rows.append(row)
            for u, e, o in zip(src.units, s["unit_excess"], s["unit_observed"]):
                unit_rows.append({"dataset": args.dataset, "seed": args.seed, "order": args.order, "scheme": scheme,
                                  "design": d, "program_id": pid, "unit": u, "observed": o, "excess": e})
            if d == "all16":
                focal_rows[-1]["n_units_excess_pos"] = int(np.sum(s["unit_excess"] > 0))
                focal_rows[-1]["signflip_p_one_sided"] = signflip_one_sided(s["unit_excess"])
    fr = pd.DataFrame(focal_rows)
    fr["q_bh_two_program"] = np.nan
    for _, g in fr.groupby(["scheme", "design"]):
        fr.loc[g.index, "q_bh_two_program"] = bh(g["p"].to_numpy())
    fr.to_csv(args.out / "focal_results.tsv", sep="\t", index=False)
    pd.DataFrame(unit_rows).to_csv(args.out / "focal_units.tsv", sep="\t", index=False)
    log("focal done\n" + fr[["scheme", "design", "program_id", "observed", "excess", "p", "q_bh_two_program",
                             "sha_match"]].to_string())
    (args.out / "run_meta.json").write_text(json.dumps(meta, indent=2, default=str))
    if args.skip_family:
        return 0

    # ---- 114-program family (descriptive analysis)
    union = set(mem["gene_symbol"])
    excl = src.genes.isin(list(union))
    rng = np.random.default_rng(args.seed)
    order = reg if args.order == "registry" else reg.iloc[::-1]
    rec_family = {}
    hf = AN4 / f"famr2_{AN4_FAMILY_TAG[args.dataset]}_s0/matched_set_hashes_{args.dataset}.tsv"
    if hf.exists():
        rec_family = dict(pd.read_csv(hf, sep="\t")[["program_id", "matched_set_sha256"]].values)
    rows, tasks = [], []
    t_sets = time.time()
    for r in order.itertuples(index=False):
        m, retained = measured(src, mem, r.program_uid)
        testable = len(m) >= 8 and retained >= 0.20
        row = {"dataset": args.dataset, "seed": args.seed, "order": args.order, "program_id": r.program_uid,
               "cell_type": r.cell_type, "module": r.module, "program_name": r.module_name,
               "n_measured": len(m), "retained_l1": retained, "testable": testable}
        if testable:
            cands = [candidates(src, gi, r.cell_type, excl)[0] for gi in m["gi"]]
            sets = draw_sets(cands, args.n_null, rng, src.n_genes)
            sha = hashlib.sha256(sets.tobytes()).hexdigest()
            w = m["w"].to_numpy(float)
            top = int(np.argmax(w))
            row.update({"matched_set_sha256": sha, "recorded_sha256": rec_family.get(r.program_uid, ""),
                        "sha_match": sha == rec_family.get(r.program_uid, ""), "_idx": m["gi"].to_numpy(), "_w": w,
                        "_top": top})
            tasks.append((r.program_uid, "all16", sets, w, args.n_sens, top))
        rows.append(row)
    log(f"family matched sets built in {time.time() - t_sets:.0f}s for {len(tasks)} programs")
    res = {}
    for key, out in pool.imap_unordered(null_task, tasks):
        res[key] = out
    pool.close()
    log("family nulls done")
    out_rows = []
    for row in rows:
        if row["testable"]:
            out = res[row["program_id"]]
            idx, w, top = row.pop("_idx"), row.pop("_w"), row.pop("_top")
            obs_sec = moran_sections((w @ RT["all16"][idx])[None, :].astype(np.float32), src.islands)[0]
            s = summarize(src, obs_sec, out["sections"])
            row.update({k: v for k, v in s.items() if not k.startswith("unit_") and not k.startswith("null_q")})
            eqw = np.full(len(w), 1 / len(w))
            lw = np.delete(w, top) / np.delete(w, top).sum()
            eo = float(moran_sections((eqw @ RT["all16"][idx])[None, :].astype(np.float32), src.islands)[0] @ src.collapse_vec)
            lo = float(moran_sections((lw @ RT["all16"][np.delete(idx, top)])[None, :].astype(np.float32),
                                      src.islands)[0] @ src.collapse_vec)
            row["sensitivity_sign_agree"] = bool(
                np.sign(s["z"]) == np.sign(eo - (out["equal"] @ src.collapse_vec).mean())
                == np.sign(lo - (out["leave_top"] @ src.collapse_vec).mean()))
        out_rows.append(row)
    fam = pd.DataFrame(out_rows)
    t = fam["testable"].astype(bool)
    fam.loc[t, "q_bh_family"] = bh(fam.loc[t, "p"].to_numpy())
    fam.loc[t, "rank_excess"] = fam.loc[t, "excess"].rank(ascending=False, method="min")
    fam.to_csv(args.out / "family_results.tsv", sep="\t", index=False)
    ft = fam[t]
    summ = {"dataset": args.dataset, "seed": args.seed, "order": args.order, "n_testable": int(t.sum()),
            "n_excess_gt0": int((ft.excess > 0).sum()), "n_q_lt_0.05": int((ft.q_bh_family < 0.05).sum()),
            "n_q_lt_0.05_and_sign": int(((ft.q_bh_family < 0.05) & ft.sensitivity_sign_agree.astype(bool)).sum()),
            "n_p_at_floor": int((ft.k == 0).sum()),
            "stromal_rank": int(ft.loc[ft.program_id == STROMAL, "rank_excess"].iloc[0]),
            "ductular_rank": int(ft.loc[ft.program_id == DUCT, "rank_excess"].iloc[0]),
            "n_sha_match": int(ft.sha_match.sum()), "n_sha_recorded": int((ft.recorded_sha256 != "").sum())}
    pd.DataFrame([summ]).to_csv(args.out / "family_summary.tsv", sep="\t", index=False)
    log(f"family summary {summ}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
