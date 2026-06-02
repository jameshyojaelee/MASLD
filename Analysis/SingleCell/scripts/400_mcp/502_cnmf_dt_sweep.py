#!/usr/bin/env python
"""
502_cnmf_dt_sweep.py — cNMF density_threshold sensitivity scan at k=16.

density_threshold controls which replicate spectra are retained for consensus:
replicates whose pairwise cosine distance to the local density peak exceeds
the threshold are excluded. Lower values are more stringent; higher values
include more replicates but may admit noisy solutions.

Reruns cobj.consensus() at dt ∈ {0.022, 0.03, 0.05} using the existing
factorize output on disk (no refactorize needed).

For each dt:
  - cobj.consensus(k=16, density_threshold=dt) writes consensus spectra,
    usages, and stability metrics to cnmf_tmp/ with the dt-tagged suffix.
  - Cophenetic correlation is read from cobj.k_selection_stats and
    independently recomputed from the retained replicate spectra.
  - Top-100 Jaccard vs the canonical dt=0.03 run checks whether program
    gene content is stable across threshold choices.

Output:
  benchmarks/cnmf_dt_sensitivity_k16.tsv
  benchmarks/factorization_showdown_cnmf.tsv  (appended later by synth step)
"""
from __future__ import annotations

import os
import shutil
import time
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import pdist, squareform

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
CNMF_ROOT = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_runs"
RUN_NAME = "global"
K = 16
DT_LIST = [float(x) for x in os.environ.get("CNMF_DT_LIST", "0.022,0.03,0.05").split(",")]
OUT_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/benchmarks"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _existing_consensus_path(dt):
    """Path to the dt=0.03 canonical consensus spectrum file (already on disk)."""
    dt_s = f"{dt:.2f}".replace(".", "_")
    return CNMF_ROOT / RUN_NAME / f"{RUN_NAME}.spectra.k_{K}.dt_{dt_s}.consensus.txt"


def top_gene_sets(spectra_df, n=100):
    """Given a DataFrame (rows = programs, cols = genes), return list[set[str]]."""
    out = []
    for i in spectra_df.index:
        vals = spectra_df.loc[i].values
        top = np.argsort(vals)[::-1][:n]
        out.append(set(spectra_df.columns[top]))
    return out


def pairwise_jaccard(ref_sets, query_sets):
    n_ref = len(ref_sets)
    n_qry = len(query_sets)
    M = np.zeros((n_qry, n_ref), dtype=float)
    for i, q in enumerate(query_sets):
        for j, r in enumerate(ref_sets):
            denom = len(q | r)
            M[i, j] = (len(q & r) / denom) if denom else 0.0
    return M


def cophenetic_from_iter_spectra(iter_files, genes_in_order):
    """Load all per-iter spectra for k=16, stack rows (programs), compute pairwise
    cosine distance, then cophenetic correlation of average-linkage tree."""
    mats = []
    for f in iter_files:
        arr = np.load(f, allow_pickle=False) if f.suffix == ".npz" else None
        if arr is None:
            continue
        # cNMF .df.npz is a pickled DataFrame wrapped via np.savez — use pandas
    # Load via pandas — the files were saved by cNMF as DataFrames using the .npz trick
    mats = []
    for f in iter_files:
        try:
            # cNMF saves via numpy's .npz format; the DataFrame is reconstructed by cnmf
            # Use the lower-level helper from cnmf module:
            from cnmf import load_df_from_npz  # type: ignore
            df = load_df_from_npz(str(f))
        except Exception:
            # Fallback: try pandas
            df = pd.read_pickle(f)
        # DataFrame expected to have rows = programs, cols = genes
        if not isinstance(df, pd.DataFrame):
            continue
        # reindex to common gene order
        df = df.reindex(columns=genes_in_order).fillna(0.0)
        # unit-row-normalize
        vals = df.values.astype(np.float64)
        norms = np.linalg.norm(vals, axis=1, keepdims=True) + 1e-10
        mats.append(vals / norms)
    if not mats:
        return float("nan"), 0
    B = np.concatenate(mats, axis=0)  # (n_iter*k, n_genes)
    # cosine distance
    D = 1.0 - (B @ B.T)
    np.fill_diagonal(D, 0.0)
    D = np.clip(D, 0.0, 2.0)
    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method="average")
    coph, _ = cophenet(Z, condensed)
    return float(coph), B.shape[0]


def main():
    from cnmf import cNMF
    print(f"[502] cNMF run={RUN_NAME} k={K} dt_list={DT_LIST}")
    cobj = cNMF(output_dir=str(CNMF_ROOT), name=RUN_NAME)

    # Ensure combined spectra exist
    cobj.combine(components=[K], skip_missing_files=True)

    # Pre-load norm_counts in float64 (consensus requires H.dtype==X.dtype)
    norm_counts = sc.read(cobj.paths["normalized_counts"])
    if norm_counts.X.dtype != np.float64:
        norm_counts.X = norm_counts.X.astype(np.float64)

    # Iter-spectra files for this k
    tmp_dir = CNMF_ROOT / RUN_NAME / "cnmf_tmp"
    iter_files = sorted(tmp_dir.glob(f"{RUN_NAME}.spectra.k_{K}.iter_*.df.npz"))
    print(f"[502] found {len(iter_files)} iter spectra files for k={K}")

    # Cophenetic from iteration spectra (ground truth independent of dt)
    genes_order = norm_counts.var_names
    coph_iter, n_stacked = cophenetic_from_iter_spectra(iter_files, genes_order)
    print(f"[502] cophenetic from iter-spectra (k={K}, n_stacked={n_stacked}): {coph_iter:.4f}")

    rows = []
    # Reference top-100 at canonical dt=0.03 (may be overwritten below by new consensus
    # if file missing — we snapshot before touching it)
    ref_consensus_path = _existing_consensus_path(0.03)
    ref_df = None
    if ref_consensus_path.exists():
        ref_df = pd.read_csv(ref_consensus_path, sep="\t", index_col=0)
        print(f"[502] canonical dt=0.03 consensus spectra: {ref_df.shape}")
    else:
        print(f"[502] WARN: canonical consensus file missing: {ref_consensus_path}")

    for dt in DT_LIST:
        t0 = time.time()
        print(f"\n[502] ===== dt = {dt} =====")
        try:
            cobj.consensus(
                k=K,
                density_threshold=dt,
                show_clustering=False,
                close_clustergram_fig=True,
                norm_counts=norm_counts,
            )
            runtime = time.time() - t0
            # Path of freshly written consensus (dt-tagged)
            dt_s = f"{dt:.2f}".replace(".", "_")
            new_path = CNMF_ROOT / RUN_NAME / f"{RUN_NAME}.spectra.k_{K}.dt_{dt_s}.consensus.txt"
            if not new_path.exists():
                # cNMF writes three digits (e.g., 0.022 -> "0_02") — probe alternatives
                alt_candidates = list((CNMF_ROOT / RUN_NAME).glob(f"{RUN_NAME}.spectra.k_{K}.dt_*consensus.txt"))
                print(f"[502] consensus file not at expected path; candidates: {[p.name for p in alt_candidates]}")
                new_path = max(alt_candidates, key=lambda p: p.stat().st_mtime) if alt_candidates else None
            if new_path is None or not new_path.exists():
                raise RuntimeError(f"No consensus file produced for dt={dt}")

            dt_df = pd.read_csv(new_path, sep="\t", index_col=0)
            # Program-level Jaccard vs dt=0.03 reference
            if ref_df is not None:
                common_cols = ref_df.columns.intersection(dt_df.columns)
                ref_top = top_gene_sets(ref_df[common_cols])
                dt_top = top_gene_sets(dt_df[common_cols])
                M = pairwise_jaccard(ref_top, dt_top)  # dt_top rows, ref cols
                # For each dt program, best match to any canonical program
                best_per_dt = M.max(axis=1)
                mean_best = float(np.mean(best_per_dt))
            else:
                best_per_dt = np.array([np.nan] * len(dt_df))
                mean_best = float("nan")

            # Stability stats file (cNMF writes one per k but overwrites across dt)
            stab_path = CNMF_ROOT / RUN_NAME / f"{RUN_NAME}.k_{K}.consensus.stats.df.npz"
            stab = None
            try:
                from cnmf import load_df_from_npz
                if stab_path.exists():
                    stab = load_df_from_npz(str(stab_path))
            except Exception:
                stab = None
            coph_cnmf = float(stab["cophenetic"].iloc[0]) if (stab is not None and "cophenetic" in stab.columns) else float("nan")
            silh_cnmf = float(stab["silhouette"].iloc[0]) if (stab is not None and "silhouette" in stab.columns) else float("nan")

            rows.append({
                "method": "cNMF",
                "k": K,
                "density_threshold": dt,
                "n_programs_retained": len(dt_df),
                "cophenetic_iter_stack": coph_iter,      # dt-independent
                "cophenetic_cnmf_stats": coph_cnmf,      # dt-dependent (from stats file)
                "silhouette_cnmf_stats": silh_cnmf,
                "mean_best_jaccard_vs_dt_0_03": mean_best,
                "runtime_s": runtime,
                "consensus_file": new_path.name,
            })
            # Per-program match details
            per_prog_rows = []
            for i, (prog_idx, jac) in enumerate(zip(dt_df.index, best_per_dt)):
                per_prog_rows.append({
                    "k": K, "density_threshold": dt, "program": int(prog_idx),
                    "best_jaccard_top100_vs_dt_0_03": float(jac),
                })
            pd.DataFrame(per_prog_rows).to_csv(
                OUT_DIR / f"cnmf_dt_program_jaccard_k{K}_dt_{dt_s}.tsv",
                sep="\t", index=False,
            )
        except Exception as e:
            print(f"[502] dt={dt} FAILED: {e}")
            rows.append({
                "method": "cNMF", "k": K, "density_threshold": dt,
                "cophenetic_iter_stack": coph_iter,
                "cophenetic_cnmf_stats": float("nan"),
                "silhouette_cnmf_stats": float("nan"),
                "mean_best_jaccard_vs_dt_0_03": float("nan"),
                "runtime_s": time.time() - t0,
                "error": str(e),
            })

    out = OUT_DIR / f"cnmf_dt_sensitivity_k{K}.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"\n[502] wrote {out}")
    print(pd.DataFrame(rows).to_string(index=False))
    print("[502] DONE.")


if __name__ == "__main__":
    main()
