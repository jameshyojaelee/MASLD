#!/usr/bin/env python
"""503b_mofa_pseudobulk_python.py — MOFA+ factorization on pseudobulk data.

Calls mofapy2.run.entry_point directly in Python (mofapy2 0.7.4, Python 3.12).
Used in the multi-method factorization comparison alongside cNMF, plain NMF,
scHPF, and LIGER.

Pipeline
--------
* Inputs: 5 pseudobulk views (donors × top-2500 HVG genes per cell type) at
  Analysis/SingleCell/results_gpu_v2/mcp/pseudobulk_mofa/view_<celltype>.tsv
* Each cell type is treated as an independent view with its own gene set
  (MOFA+'s native multi-view design; Argelaguet 2020 PMID 32393329).
  Genes are NOT intersected across views.
* Donors are aligned across views by union, with NaN padding for missing donors.
* log1p variance-stabilizing transform + per-view scaling
  (set_data_options(scale_views=True)).
* k in {10, 16} matching the cNMF / scHPF / LIGER comparison.

Outputs
-------
* mofa_k{K}.hdf5 — full MOFA+ model
* benchmarks/factorization_showdown_mofa.tsv — long-format showdown rows
* benchmarks/program_jaccard_mofa_vs_cnmf.tsv — per-MOFA-factor best Jaccard
  vs canonical cNMF spectra
* benchmarks/mofa_factors_k{K}.tsv — per-view factor weights, long-format
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
PSEUDO = MCP / "pseudobulk_mofa"
BENCH = MCP / "benchmarks"
BENCH.mkdir(parents=True, exist_ok=True)

CELL_TYPES = ["Hepatocytes", "Endothelial_cells", "Macrophages",
              "Fibroblasts", "Cholangiocytes"]


# -----------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------
def load_views() -> dict[str, pd.DataFrame]:
    """Load 5 pseudobulk views as gene x donor frames; log1p + variance-rank.

    Each view file is genes x donors with header row = donor IDs and column
    1 = gene symbol.
    """
    out: dict[str, pd.DataFrame] = {}
    for ct in CELL_TYPES:
        path = PSEUDO / f"view_{ct}.tsv"
        df = pd.read_csv(path, sep="\t", index_col=0)
        # Drop fully-NaN / fully-zero rows (no signal)
        df = df.loc[~df.isna().all(axis=1)]
        df = df.loc[(df.fillna(0).abs().sum(axis=1) > 0)]
        # Variance-stabilise: log1p of raw counts
        df = np.log1p(df.fillna(0))
        # Top-2500 by per-view variance (matches cNMF HVG count). Inputs are
        # already top-2500 per the existing manifest, but re-rank anyway in
        # case duplicates / zero-variance rows slipped through.
        v = df.var(axis=1)
        keep = v.sort_values(ascending=False).index[: min(2500, len(df))]
        out[ct] = df.loc[keep]
        print(f"[503b] view {ct}: {out[ct].shape} (genes x donors)", flush=True)
    return out


def align_donors(views: dict[str, pd.DataFrame]) -> tuple[dict[str, pd.DataFrame], list[str]]:
    """Align donors across views by union; NaN pad missing entries.

    Returns dict of gene x donor frames with the same ordered donor list.
    """
    donor_union: list[str] = []
    seen: set[str] = set()
    for df in views.values():
        for d in df.columns:
            if d not in seen:
                seen.add(d)
                donor_union.append(d)
    aligned: dict[str, pd.DataFrame] = {}
    for ct, df in views.items():
        missing = [d for d in donor_union if d not in df.columns]
        if missing:
            pad = pd.DataFrame(
                np.full((df.shape[0], len(missing)), np.nan),
                index=df.index, columns=missing,
            )
            df = pd.concat([df, pad], axis=1)
        aligned[ct] = df.loc[:, donor_union]
    print(f"[503b] aligned {len(donor_union)} donors across {len(views)} views",
          flush=True)
    return aligned, donor_union


# -----------------------------------------------------------------------
# MOFA+ runner
# -----------------------------------------------------------------------
def fit_mofa(views: dict[str, pd.DataFrame], donors: list[str], k: int,
             out_hdf5: Path, max_iter: int = 1000,
             conv_mode: str = "medium", seed: int = 42) -> dict:
    """Fit MOFA+ with per-view independent gene sets.

    mofapy2.run.entry_point.set_data_matrix expects a list of lists of numpy
    arrays: data[m][g] -> (samples, features) for view m, group g. Single
    group = nested list with a single inner list of length M.
    """
    from mofapy2.run.entry_point import entry_point  # type: ignore

    ent = entry_point()

    # MOFA+ wants samples x features. Transpose each view's frame.
    # Data layout: data[m][g] -> view m, group g (samples, features).
    # We have a single group ("all"), so wrap each view's matrix as [mat].
    # MOFA+ requires globally-unique feature names across views (chain check),
    # so prefix each gene symbol with the view name. Strip the prefix when
    # writing outputs.
    feature_names: list[list[str]] = []
    view_names = list(views.keys())
    data_nested: list[list[np.ndarray]] = []
    for ct in view_names:
        m = views[ct].T.to_numpy(dtype=np.float64)  # donors x genes
        data_nested.append([m])
        feature_names.append([f"{ct}__{g}" for g in views[ct].index])

    ent.set_data_options(scale_views=True, scale_groups=False)
    ent.set_data_matrix(
        data_nested,
        likelihoods=["gaussian"] * len(view_names),
        views_names=view_names,
        groups_names=["all"],
        samples_names=[donors],
        features_names=feature_names,
    )
    ent.set_model_options(
        factors=k,
        ard_factors=True,
        ard_weights=True,
        spikeslab_factors=False,
        spikeslab_weights=True,
    )
    # MOFA+ ARD prunes low-variance factors at each ELBO check using the
    # dropR2 threshold. For benchmark fairness with cNMF/scHPF/LIGER (which
    # don't prune), set dropR2=None so all k requested factors are retained.
    ent.set_train_options(
        iter=max_iter,
        convergence_mode=conv_mode,
        startELBO=1,
        freqELBO=10,
        dropR2=None,
        gpu_mode=False,
        verbose=False,
        seed=seed,
    )
    ent.build()
    t0 = time.time()
    ent.run()
    runtime = time.time() - t0
    out_hdf5.parent.mkdir(parents=True, exist_ok=True)
    if out_hdf5.exists():
        out_hdf5.unlink()
    ent.save(str(out_hdf5))
    print(f"[503b] MOFA+ k={k} done in {runtime:.1f}s -> {out_hdf5}",
          flush=True)
    return {"runtime_s": runtime}


# -----------------------------------------------------------------------
# Post-fit metric extraction (from HDF5 directly so we don't depend on a
# separate `mofax` package).
# -----------------------------------------------------------------------
def extract_from_hdf5(path: Path) -> dict:
    """Pull per-view weights and variance-explained from the MOFA+ HDF5."""
    out: dict = {"weights": {}, "feature_names": {}, "view_names": []}
    with h5py.File(path, "r") as f:
        view_names = [v.decode() if isinstance(v, bytes) else v
                      for v in f["views/views"][:]]
        out["view_names"] = view_names
        # Per-view weights live at expectations/W/<view>
        for vn in view_names:
            w = f[f"expectations/W/{vn}"][:]
            # Convention in mofapy2: stored as (factors, features) in HDF5
            # (ndim==2). We standardise to (features, factors) in our outputs.
            if w.shape[0] < w.shape[1]:
                w = w.T
            out["weights"][vn] = w
            feats = f[f"features/{vn}"][:]
            decoded = [
                x.decode() if isinstance(x, bytes) else x for x in feats
            ]
            # Strip per-view prefix added in fit_mofa() so downstream Jaccard
            # uses bare HUGO symbols matching cNMF spectra column names.
            prefix = f"{vn}__"
            out["feature_names"][vn] = [
                s[len(prefix):] if s.startswith(prefix) else s
                for s in decoded
            ]
        # variance_explained / r2 per view per factor (single group "all")
        ve_dict: dict[str, np.ndarray] = {}
        if "variance_explained/r2_per_factor" in f:
            grp = f["variance_explained/r2_per_factor"]
            for g_name in grp:
                arr = grp[g_name][:]
                # arr shape: (factors, views) typically
                ve_dict[g_name] = arr
        out["r2_per_factor"] = ve_dict
        # Total ELBO trajectory
        if "training_stats/elbo" in f:
            out["elbo_trajectory"] = f["training_stats/elbo"][:]
        else:
            out["elbo_trajectory"] = np.array([])
    return out


# -----------------------------------------------------------------------
# Jaccard vs cNMF
# -----------------------------------------------------------------------
def load_cnmf_spectra(k: int) -> pd.DataFrame:
    """Return cNMF spectra at canonical dt=0.03 as factors x genes."""
    path = MCP / "cnmf_runs/global" / f"global.gene_spectra_score.k_{k}.dt_0_03.txt"
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df  # rows = factor index, columns = gene names


def top_n_genes(weights: np.ndarray, names: list[str], n: int = 100) -> list[set[str]]:
    """Top-n genes per factor by absolute weight."""
    n_feat, n_factor = weights.shape
    out = []
    for k in range(n_factor):
        order = np.argsort(np.abs(weights[:, k]))[::-1][:n]
        out.append({names[i] for i in order})
    return out


def jaccard(a: set, b: set) -> float:
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def per_factor_jaccard(mofa_extracted: dict, k: int) -> pd.DataFrame:
    """For each MOFA factor, compute its top-100 gene set as the union over
    views (top-100 by absolute weight per view, unioned). Then for each cNMF
    program, compute Jaccard vs the union; report max as best_jaccard_top100.

    Returns long-format frame: method, k, factor, best_jaccard_top100,
    best_match_cnmf_program.
    """
    cnmf = load_cnmf_spectra(k)
    cnmf_top: dict[str, set[str]] = {}
    for idx in cnmf.index:
        row = cnmf.loc[idx]
        top = row.abs().sort_values(ascending=False).index[:100]
        cnmf_top[f"P{idx}"] = set(top)

    # Build MOFA top-100 sets per factor (union across views)
    factor_sets: list[set[str]] = []
    n_factors = next(iter(mofa_extracted["weights"].values())).shape[1]
    for fi in range(n_factors):
        union: set[str] = set()
        for vn, W in mofa_extracted["weights"].items():
            names = mofa_extracted["feature_names"][vn]
            order = np.argsort(np.abs(W[:, fi]))[::-1][:100]
            union |= {names[i] for i in order}
        factor_sets.append(union)

    # Best Jaccard per MOFA factor vs any cNMF program
    rows = []
    for fi, mset in enumerate(factor_sets, start=1):
        best_j = 0.0
        best_p = ""
        for pn, cset in cnmf_top.items():
            j = jaccard(mset, cset)
            if j > best_j:
                best_j = j
                best_p = pn
        rows.append({
            "method": "MOFA",
            "k": k,
            "factor": f"factor_{fi}",
            "best_jaccard_top100": best_j,
            "best_match_cnmf_program": best_p,
        })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------
# Showdown TSV writers
# -----------------------------------------------------------------------
def append_showdown_rows(rows: list[dict], path: Path) -> None:
    df_new = pd.DataFrame(rows)
    if path.exists():
        df_old = pd.read_csv(path, sep="\t")
        df = pd.concat([df_old, df_new], ignore_index=True)
    else:
        df = df_new
    df.to_csv(path, sep="\t", index=False)
    print(f"[503b] wrote {path} ({len(df)} rows)", flush=True)


def write_factor_weights(mofa_extracted: dict, k: int, path: Path) -> None:
    """Long-format: view, factor, gene, weight, abs_weight."""
    out_rows: list[pd.DataFrame] = []
    for vn, W in mofa_extracted["weights"].items():
        names = mofa_extracted["feature_names"][vn]
        n_feat, n_factor = W.shape
        for fi in range(n_factor):
            df = pd.DataFrame({
                "view": vn,
                "factor": f"factor_{fi+1}",
                "gene": names,
                "weight": W[:, fi],
                "abs_weight": np.abs(W[:, fi]),
            })
            out_rows.append(df)
    df_all = pd.concat(out_rows, ignore_index=True)
    df_all.to_csv(path, sep="\t", index=False)
    print(f"[503b] wrote {path} ({len(df_all)} rows)", flush=True)


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------
def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--k", type=int, action="append", required=True,
                   help="Repeatable. e.g. --k 10 --k 16")
    p.add_argument("--smoke", action="store_true",
                   help="Smoke test: 200 iters max, --k must be 2")
    p.add_argument("--out-dir", type=str, default=str(MCP / "mofa_runs"))
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[503b] mode={'SMOKE' if args.smoke else 'FULL'}; k={args.k}",
          flush=True)

    views = load_views()
    aligned, donors = align_donors(views)

    showdown_rows: list[dict] = []

    for k in args.k:
        out_hdf5 = out_dir / f"mofa_k{k}.hdf5"
        max_iter = 200 if args.smoke else 1000
        conv_mode = "fast" if args.smoke else "medium"
        info = fit_mofa(aligned, donors, k, out_hdf5,
                        max_iter=max_iter, conv_mode=conv_mode,
                        seed=args.seed)
        ext = extract_from_hdf5(out_hdf5)

        # Self-validation
        assert len(ext["weights"]) == 5, \
            f"expected 5 views, got {len(ext['weights'])}"
        # MOFA+ may retain fewer factors than requested if ARD prunes them.
        # We log the effective k actually fit. Allow <= k.
        actual_k = next(iter(ext["weights"].values())).shape[1]
        if actual_k != k:
            print(f"[503b]   WARNING: requested k={k}, MOFA+ retained "
                  f"actual_k={actual_k} after ARD pruning", flush=True)
        for vn, W in ext["weights"].items():
            n_feat = len(ext["feature_names"][vn])
            assert W.shape[0] == n_feat, \
                f"view {vn} weight rows {W.shape[0]} != n_feat {n_feat}"
            assert W.shape[1] == actual_k, \
                f"view {vn} weight cols {W.shape[1]} != actual_k {actual_k}"
        elbo = ext["elbo_trajectory"]
        elbo_finite = elbo[np.isfinite(elbo)] if elbo.size else np.array([])
        if elbo_finite.size:
            elbo_final = float(elbo_finite[-1])
            assert np.isfinite(elbo_final), "non-finite ELBO at convergence"
            print(f"[503b]   elbo final={elbo_final:.3e}, "
                  f"n_elbo_pts={elbo_finite.size}, n_iter={elbo.size}",
                  flush=True)
        else:
            elbo_final = float("nan")
            print(f"[503b]   WARNING: no finite ELBO points recorded",
                  flush=True)

        # Variance explained per view per factor
        if ext["r2_per_factor"]:
            grp_name = next(iter(ext["r2_per_factor"]))
            r2 = ext["r2_per_factor"][grp_name]
            # r2 shape: (n_views, n_factors) or (n_factors, n_views) — handle both
            if r2.shape[0] == k:
                r2 = r2.T  # -> (n_views, n_factors)
            view_names = ext["view_names"]
            for vi, vn in enumerate(view_names):
                for fi in range(k):
                    showdown_rows.append({
                        "method": "MOFA",
                        "k": k,
                        "replicate": -1,
                        "metric": f"variance_explained_view_{vn}_factor_{fi+1}",
                        "value": float(r2[vi, fi]),
                        "note": (f"mofapy2 0.7.4; per-view-HVG (top-2500); "
                                 f"max_iter={max_iter}; conv_mode={conv_mode}; "
                                 f"seed={args.seed}"),
                    })
            # Total per-view
            for vi, vn in enumerate(view_names):
                showdown_rows.append({
                    "method": "MOFA",
                    "k": k,
                    "replicate": -1,
                    "metric": f"total_variance_explained_view_{vn}",
                    "value": float(np.nansum(r2[vi, :])),
                    "note": ("mofapy2 0.7.4; sum across factors; "
                             f"seed={args.seed}"),
                })

        # ELBO summary
        if elbo_finite.size:
            showdown_rows.append({
                "method": "MOFA",
                "k": k,
                "replicate": -1,
                "metric": "elbo_final",
                "value": elbo_final,
                "note": (f"mofapy2 0.7.4; n_elbo_pts={elbo_finite.size}; "
                         f"max_iter={max_iter}; seed={args.seed}"),
            })
            showdown_rows.append({
                "method": "MOFA",
                "k": k,
                "replicate": -1,
                "metric": "elbo_n_iter",
                "value": float(elbo.size),
                "note": (f"mofapy2 0.7.4; max_iter={max_iter}; total iters "
                         "incl. NaN-filled non-eval steps"),
            })
        # Runtime
        showdown_rows.append({
            "method": "MOFA",
            "k": k,
            "replicate": -1,
            "metric": "runtime_s",
            "value": float(info["runtime_s"]),
            "note": (f"mofapy2 0.7.4; conv_mode={conv_mode}; "
                     f"max_iter={max_iter}"),
        })

        # Per-factor weights TSV (only for non-smoke)
        if not args.smoke:
            write_factor_weights(
                ext, k,
                BENCH / f"mofa_factors_k{k}.tsv",
            )
            # Jaccard vs cNMF
            jdf = per_factor_jaccard(ext, k)
            jdf_path = BENCH / "program_jaccard_mofa_vs_cnmf.tsv"
            if jdf_path.exists():
                jdf_old = pd.read_csv(jdf_path, sep="\t")
                jdf = pd.concat([jdf_old, jdf], ignore_index=True)
            jdf.to_csv(jdf_path, sep="\t", index=False)
            print(f"[503b] wrote {jdf_path}", flush=True)
            # Showdown summary: max + mean Jaccard
            jac_only = jdf[jdf["k"] == k]["best_jaccard_top100"].astype(float)
            showdown_rows.append({
                "method": "MOFA",
                "k": k,
                "replicate": -1,
                "metric": "jaccard_top100_vs_cnmf",
                "value": float(jac_only.mean()) if len(jac_only) else float("nan"),
                "note": (f"mean over {len(jac_only)} MOFA factors; per-MOFA-factor "
                         "best Jaccard vs canonical cNMF dt=0.03 spectra"),
            })
            showdown_rows.append({
                "method": "MOFA",
                "k": k,
                "replicate": -1,
                "metric": "jaccard_top100_vs_cnmf_max",
                "value": float(jac_only.max()) if len(jac_only) else float("nan"),
                "note": "max over MOFA factors",
            })

    # Write showdown TSV (overwrite, idempotent)
    out_path = BENCH / "factorization_showdown_mofa.tsv"
    if args.smoke:
        # Don't overwrite the canonical file with smoke-test rows; route to a
        # separate file so the caller can inspect it.
        out_path = BENCH / "factorization_showdown_mofa_SMOKE.tsv"
    pd.DataFrame(showdown_rows).to_csv(out_path, sep="\t", index=False)
    print(f"[503b] wrote {out_path} ({len(showdown_rows)} rows)", flush=True)
    print("[503b] DONE.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
