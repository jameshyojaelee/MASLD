#!/usr/bin/env python
"""
490c_mouse_null_vectorized.py — Vectorized mouse cross-species conservation null.

Computes a permutation null for each cNMF program's conservation score on
the 562K-cell mouse atlas. For each permutation, a random set of mouse genes
(matched in size to the program's mapped ortholog count) is scored and
aggregated across cells. The observed program score uses the same direct
per-cell mean (log1p-CPM space) so observed and null are directly comparable.

Scoring method: direct per-cell mean of gene-set expression in log1p-CPM space,
rather than the control-gene-binning approach used by sc.tl.score_genes. This
is faster and avoids the dependency on control gene selection. The same
strict 1:1 ortholog chain from 490b is reused.

Precomputation (once):
  1. Load atlas backed; normalize_total → log1p of raw counts.
  2. Build the union of candidate mouse genes that any random HVG draw could
     ever map into = strict-ortho mapping of *every* HVG human gene.  This is
     ≤ ~1.5k mouse genes.  Slice X to (n_cells, n_candidate) → ~2 GB float32
     (one-time densification).
  3. Stack n_perm random gene-set indicators as a sparse matrix
     `M` of shape (n_perm, n_candidate).  Matmul: `cell_scores = X_cand @
     M.T / sz`.  Mean across cells → per-perm null mean (size-matched).
  4. Optional GPU acceleration via cupy when available (`--use-gpu`).

Output (one chunk per array-task):
  validation/mouse/null_chunks/null_chunk_<seed_offset>.tsv
  with columns:  perm_id, program_size, n_strict, score (cell-mean),
                 score_per_donor (median across 92 donors).

The aggregator (490d_mouse_null_aggregate.py) then computes per-program
permutation p (observed quantile vs null), BH-FDR.

Reusable observed scoring is provided via `--mode observed` (computes the
direct-mean observed program scores in the same metric used by the null).
"""
from __future__ import annotations
import argparse, os, gzip, re, time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
ATLAS = ROOT / "Analysis/Deconvolution/reference/reference_mouse.h5ad"
ORTHO = ROOT / "archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"
GENCODE = ROOT / "data/gencode_v49_gene_metadata.tsv.gz"
HVG_FILE = MCP / "cnmf_runs/global/global.overdispersed_genes.txt"
TOPGENES = MCP / "cnmf_annot/global/program_topgenes.k16.tsv"
OUT_DIR = MCP / "validation/mouse"
NULL_CHUNK_DIR = OUT_DIR / "null_chunks"


def _build_chain_maps():
    gm = pd.read_csv(GENCODE, sep="\t")
    sym2ens = dict(zip(gm["gene_name"], gm["ensembl_base"]))
    o = pd.read_csv(ORTHO, sep="\t")
    o = o[o["orthology_type"] == "ortholog_one2one"]
    h2m = dict(zip(o["human_ensembl_gene_id"], o["mouse_ensembl_gene_id"]))

    gtf = ROOT / "../../../../home/jameslee/reference_genome/refdata-gex-GRCm39-2024-A/annotation/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf.gz"
    mouse_gtf = Path(str(gtf))
    m_ens2sym: dict[str, str] = {}
    if mouse_gtf.exists():
        pat = re.compile(r'gene_id "([^"]+)".*gene_name "([^"]+)"')
        with gzip.open(mouse_gtf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if "\tgene\t" not in line:
                    continue
                m = pat.search(line)
                if m:
                    m_ens2sym[m.group(1).split(".")[0]] = m.group(2)
    return sym2ens, h2m, m_ens2sym


def _strict_map_sym_to_mouse(human_syms, sym2ens, h2m, m_ens2sym, mouse_var_set):
    out = []
    for hs in human_syms:
        e = sym2ens.get(hs)
        if not e:
            continue
        mens = h2m.get(e)
        if not mens:
            continue
        msym = m_ens2sym.get(mens)
        if not msym:
            continue
        if msym in mouse_var_set:
            out.append(msym)
    return out


def _load_and_normalize_atlas():
    print(f"[490c] loading atlas {ATLAS} (full-load, then normalize_total + log1p)")
    t0 = time.time()
    a = sc.read_h5ad(ATLAS)
    a.obs_names_make_unique()
    # Atlas X confirmed to be raw counts (max=298, layer 'counts' identical) → normalize.
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    print(f"[490c]   shape={a.shape} loaded+normed in {time.time()-t0:.1f}s")
    return a


def _build_candidate_pool(a, hvg_pool_h, sym2ens, h2m, m_ens2sym):
    """Return:
        cand_mouse_syms: list[str]     mouse symbols reachable from HVG pool
        cand_idx:        np.ndarray    column indices into a.X
        h2cand:          dict[hsym→cand_idx]
    """
    mouse_var_set = set(a.var_names)
    h2cand: dict[str, int] = {}        # human_sym → column in cand_mouse_syms
    msym_to_col: dict[str, int] = {}   # mouse_sym → column (dedup)
    cand_mouse_syms: list[str] = []
    for hs in hvg_pool_h:
        msyms = _strict_map_sym_to_mouse([hs], sym2ens, h2m, m_ens2sym, mouse_var_set)
        if not msyms:
            continue
        msym = msyms[0]
        if msym in msym_to_col:
            h2cand[hs] = msym_to_col[msym]
            continue
        col = len(cand_mouse_syms)
        msym_to_col[msym] = col
        h2cand[hs] = col
        cand_mouse_syms.append(msym)
    var2idx = {v: i for i, v in enumerate(a.var_names)}
    cand_idx = np.array([var2idx[m] for m in cand_mouse_syms], dtype=np.int64)
    return cand_mouse_syms, cand_idx, h2cand


def _slice_dense(a, cand_idx):
    print(f"[490c] slicing X[:, {len(cand_idx)} candidate cols] → dense float32")
    t0 = time.time()
    X = a.X[:, cand_idx]
    if sp.issparse(X):
        Xd = np.asarray(X.todense(), dtype=np.float32)
    else:
        Xd = np.asarray(X, dtype=np.float32)
    print(f"[490c]   X_cand shape={Xd.shape} dtype={Xd.dtype} ({Xd.nbytes/1e9:.2f} GB) in {time.time()-t0:.1f}s")
    return Xd


def _xp(use_gpu: bool):
    if not use_gpu:
        return np, False
    try:
        import cupy as cp
        cp.zeros(1)  # warm-up; raises if no GPU
        return cp, True
    except Exception as e:
        print(f"[490c] cupy unavailable ({e}); falling back to numpy")
        return np, False


def _program_sizes(a, sym2ens, h2m, m_ens2sym):
    topg = pd.read_csv(TOPGENES, sep="\t")
    mouse_var_set = set(a.var_names)
    rows = []
    prog_strict: dict[str, list[str]] = {}
    for p, g in topg.groupby("program"):
        h_syms = g["gene_name"].tolist()[:100]
        strict = _strict_map_sym_to_mouse(h_syms, sym2ens, h2m, m_ens2sym, mouse_var_set)
        prog_strict[str(p)] = strict
        rows.append({"program": str(p), "n_human_top": len(h_syms), "n_strict_1to1": len(strict)})
    return pd.DataFrame(rows), prog_strict


def _cell_to_donor_index(a) -> tuple[np.ndarray, list[str]]:
    """Return (donor_index_per_cell, donor_levels) — donor_index is int."""
    samples = a.obs["sample"].astype(str).values
    levels = sorted(set(samples))
    lev2i = {s: i for i, s in enumerate(levels)}
    idx = np.fromiter((lev2i[s] for s in samples), dtype=np.int32, count=len(samples))
    return idx, levels


def run_observed(args):
    """Compute observed direct-mean scores per program (cell-mean and donor-median).
    Writes a single-row-per-program TSV used as ground truth for null comparison.
    """
    sym2ens, h2m, m_ens2sym = _build_chain_maps()
    a = _load_and_normalize_atlas()

    obs_table, prog_strict = _program_sizes(a, sym2ens, h2m, m_ens2sym)

    # Slice once on the union of all program genes.
    union_syms: list[str] = []
    seen = set()
    for p, gl in prog_strict.items():
        for g in gl:
            if g not in seen:
                seen.add(g)
                union_syms.append(g)
    var2idx = {v: i for i, v in enumerate(a.var_names)}
    union_idx = np.array([var2idx[g] for g in union_syms], dtype=np.int64)
    Xd = _slice_dense(a, union_idx)
    sym2col = {g: i for i, g in enumerate(union_syms)}

    donor_idx, donor_levels = _cell_to_donor_index(a)
    n_donors = len(donor_levels)
    n_cells = Xd.shape[0]

    rows = []
    for p, gl in prog_strict.items():
        sz = len(gl)
        if sz < 5:
            rows.append({"program": p, "n_strict": sz,
                         "obs_score_cell_mean": np.nan,
                         "obs_score_donor_mean": np.nan,
                         "obs_score_donor_median": np.nan})
            continue
        cols = np.array([sym2col[g] for g in gl], dtype=np.int64)
        cell_score = Xd[:, cols].mean(axis=1)         # (n_cells,)
        cell_mean = float(cell_score.mean())
        # donor-aggregated: cell→donor mean
        sums = np.bincount(donor_idx, weights=cell_score, minlength=n_donors)
        counts = np.bincount(donor_idx, minlength=n_donors).astype(np.float64)
        counts[counts == 0] = 1
        donor_means = sums / counts
        rows.append({"program": p, "n_strict": sz,
                     "obs_score_cell_mean": cell_mean,
                     "obs_score_donor_mean": float(donor_means.mean()),
                     "obs_score_donor_median": float(np.median(donor_means))})
    df = pd.DataFrame(rows).sort_values("program", key=lambda s: s.astype(int))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_f = OUT_DIR / f"mouse_observed_directmean_k{args.k}.tsv"
    df.to_csv(out_f, sep="\t", index=False)
    print(f"[490c.observed] wrote {out_f}")
    print(df.to_string(index=False))


def run_null(args):
    sym2ens, h2m, m_ens2sym = _build_chain_maps()

    if not HVG_FILE.exists():
        raise FileNotFoundError(HVG_FILE)
    hvg_pool_h = [ln.strip() for ln in HVG_FILE.read_text().splitlines() if ln.strip()]
    print(f"[490c.null] HVG pool size: {len(hvg_pool_h)}")

    a = _load_and_normalize_atlas()

    cand_mouse, cand_idx, h2cand = _build_candidate_pool(
        a, hvg_pool_h, sym2ens, h2m, m_ens2sym
    )
    print(f"[490c.null] candidate mouse-mappable HVGs: {len(cand_mouse)} / {len(hvg_pool_h)}")

    Xd = _slice_dense(a, cand_idx)              # (n_cells, n_cand) float32
    n_cells, n_cand = Xd.shape

    # Donor index for fast bincount aggregation.
    donor_idx, donor_levels = _cell_to_donor_index(a)
    n_donors = len(donor_levels)
    print(f"[490c.null] donors: {n_donors}")

    # Program sizes (for size-matched null), then unique sizes used.
    obs_table, prog_strict = _program_sizes(a, sym2ens, h2m, m_ens2sym)
    sizes_unique = sorted({len(v) for v in prog_strict.values() if len(v) >= 5})
    print(f"[490c.null] program sizes to fit: {sizes_unique}")

    # GPU?
    xp, on_gpu = _xp(args.use_gpu)
    if on_gpu:
        Xg = xp.asarray(Xd)
        donor_idx_g = xp.asarray(donor_idx)
    else:
        Xg = Xd
        donor_idx_g = donor_idx

    rng = np.random.default_rng(args.seed_offset)
    NULL_CHUNK_DIR.mkdir(parents=True, exist_ok=True)
    out_path = Path(args.output) if args.output else (
        NULL_CHUNK_DIR / f"null_chunk_seed{args.seed_offset:06d}.tsv"
    )

    # We sample human genes from the HVG pool; only those that map into cand_mouse
    # contribute. To match a target size `sz`, oversample then truncate (same logic
    # as 490b but vectorized below).
    h_indices = np.arange(len(hvg_pool_h))
    n_h = len(hvg_pool_h)

    rows = []
    for sz in sizes_unique:
        t0 = time.time()
        # Build a sparse indicator matrix M of shape (n_perm, n_cand).
        # Each row has exactly `sz` ones (or however many we managed to map).
        rows_csr: list[int] = []
        cols_csr: list[int] = []
        sizes_actual: list[int] = []

        oversample = sz * 5 + 50            # generous so most draws yield ≥ sz
        oversample = min(oversample, n_h)
        for it in range(args.n_perms):
            tries = 0
            hits_set: set[int] = set()
            hits_order: list[int] = []
            while True:
                tries += 1
                draw = rng.choice(h_indices, size=oversample, replace=False)
                for hi in draw:
                    hsym = hvg_pool_h[hi]
                    ci = h2cand.get(hsym)
                    if ci is None or ci in hits_set:
                        continue
                    hits_set.add(ci)
                    hits_order.append(ci)
                    if len(hits_order) >= sz:
                        break
                if len(hits_order) >= sz:
                    hits_order = hits_order[:sz]
                    break
                if tries > 5:
                    break
            if len(hits_order) < 5:
                sizes_actual.append(0)
                continue
            sizes_actual.append(len(hits_order))
            for ci in hits_order:
                rows_csr.append(it)
                cols_csr.append(ci)

        # Build CSR.  Zero-rowed perms (size_actual=0) stay all-zero → score=0; we
        # filter them out below by sizes_actual.
        n_eff = args.n_perms
        data = np.ones(len(rows_csr), dtype=np.float32)
        M = sp.csr_matrix(
            (data, (np.asarray(rows_csr, dtype=np.int64),
                    np.asarray(cols_csr, dtype=np.int64))),
            shape=(n_eff, n_cand),
            dtype=np.float32,
        )
        sizes_actual_arr = np.asarray(sizes_actual, dtype=np.int32)

        # cell_scores[perm, cell] = (X_cand[cell, :] · M[perm, :]) / sizes_actual[perm]
        # Equivalently:  S = (X_cand @ M.T) / sizes_actual
        # Shape:          (n_cells, n_perm)
        if on_gpu:
            import cupy as cp
            # M is small (n_perm × n_cand, ≤ 100 × 1.5k = 150k entries) —
            # densify on GPU, then dense @ dense is one cuBLAS call.
            Mt_dense = cp.asarray(M.T.toarray())            # (n_cand, n_perm)
            S = Xg @ Mt_dense                               # (n_cells, n_perm)
            denom = cp.asarray(sizes_actual_arr.clip(min=1).astype(np.float32))
            S = S / denom[None, :]
            cell_mean_per_perm = cp.asnumpy(S.mean(axis=0))
            S_cpu = cp.asnumpy(S)
            donor_idx_cpu = (cp.asnumpy(donor_idx_g)
                             if hasattr(donor_idx_g, 'get') else donor_idx_g)
            del Mt_dense, S
            cp.get_default_memory_pool().free_all_blocks()
        else:
            # CPU path: sparse Mt @ X.T → (n_perm, n_cells), reshape and divide.
            Mt = M.tocsr()                                  # (n_perm, n_cand)
            S_pn = Mt @ Xg.T                                # (n_perm, n_cells)
            S_pn = np.asarray(S_pn)
            denom = sizes_actual_arr.clip(min=1).astype(np.float32)
            S = (S_pn / denom[:, None]).astype(np.float32)  # (n_perm, n_cells)
            S = S.T                                         # (n_cells, n_perm)
            cell_mean_per_perm = S.mean(axis=0)
            S_cpu = S
            donor_idx_cpu = donor_idx

        # Donor-aggregated mean per perm.
        donor_mean_per_perm = np.empty(n_eff, dtype=np.float32)
        donor_median_per_perm = np.empty(n_eff, dtype=np.float32)
        counts_per_donor = np.bincount(donor_idx_cpu, minlength=n_donors).astype(np.float32)
        counts_safe = np.where(counts_per_donor == 0, 1.0, counts_per_donor)
        for it in range(n_eff):
            sums = np.bincount(donor_idx_cpu, weights=S_cpu[:, it].astype(np.float64),
                               minlength=n_donors)
            dmeans = (sums / counts_safe).astype(np.float32)
            donor_mean_per_perm[it] = float(dmeans.mean())
            donor_median_per_perm[it] = float(np.median(dmeans))

        for it in range(n_eff):
            rows.append({
                "perm_id": int(args.seed_offset + it),
                "program_size": int(sz),
                "n_strict": int(sizes_actual_arr[it]),
                "score_cell_mean": float(cell_mean_per_perm[it]) if sizes_actual_arr[it] >= 5 else np.nan,
                "score_donor_mean": float(donor_mean_per_perm[it]) if sizes_actual_arr[it] >= 5 else np.nan,
                "score_donor_median": float(donor_median_per_perm[it]) if sizes_actual_arr[it] >= 5 else np.nan,
            })
        print(f"[490c.null] size={sz}: {n_eff} perms in {time.time()-t0:.1f}s "
              f"(median actual size {int(np.median(sizes_actual_arr))})")

    df = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[490c.null] wrote chunk → {out_path}  ({len(df)} rows)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["observed", "null"], default="null")
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--n-perms", type=int, default=100)
    ap.add_argument("--seed-offset", type=int, default=0)
    ap.add_argument("--output", type=str, default="")
    ap.add_argument("--use-gpu", action="store_true")
    args = ap.parse_args()
    if args.mode == "observed":
        run_observed(args)
    else:
        run_null(args)
