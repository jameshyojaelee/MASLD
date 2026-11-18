#!/usr/bin/env python3
"""Refit masld-panel-v1.1 on POOL_V2 (557 samples, label-free candidate abundance + biotype filter).

Strict separation:
  - Candidate filter uses ONLY POOL_V2 expression (median CPM >= 5.0, frac_zero <= 0.05) and
    GENCODE v49 annotation (protein_coding, non-chrY, excluding olfactory/testis families).
  - No external cohort (GSE193066, GSE268273, GSE213621, GSE276114) touches selection or fitting.
  - Evaluates GSE268273 (OOD n=109) only AFTER fitting to record expected.json.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
REPO = BENCH.parents[1]
PANELDIR = BENCH / "executions/deployable-panel-20260904T214938Z"
sys.path.insert(0, str(PANELDIR))
sys.path.insert(0, str(BENCH / "executions/standardized-metric-panel-20260902T130923Z"))

import substrate as sb  # noqa: E402
import selection as sel  # noqa: E402
import panel_core as pc  # noqa: E402
from panel_metrics import spearman, tie_ceiling  # noqa: E402

SEED = 20261118
K = 20
CANDIDATE_M = 500
REDUNDANCY = 0.90
MIN_MEDIAN_CPM = 5.0
MAX_FRAC_ZERO = 0.05

EXCLUDED_FAMILIES_RE = re.compile(
    r"^(?:OR\d+[A-Z]|TAAR\d|TAS[12]R\d|VN[12]R|SPANX|GAGE|MAGE|CSAG|CTAG|PRAME|"
    r"SSX\d|DAZ|RBMY|TSPY|PRM\d|TNP\d|SPACA|SPAG|TEX\d|SYCP|HORMAD|PIWIL|DMRT|SOHLH)",
    re.IGNORECASE,
)


def main():
    np.random.seed(SEED)
    W = sb.release_weights()
    AXIS = W["axis"]

    gc = pd.read_csv(REPO / "data/gencode_v49_gene_metadata.tsv.gz", sep="\t", keep_default_na=False)
    gc_map = gc.drop_duplicates("ensembl_base").set_index("ensembl_base")

    V2 = sb.pool_v2(AXIS)
    cpm = V2["cpm"]
    y = V2["y"]
    coh = V2["cohort"]
    cohorts = V2["cohorts"]

    med_cpm = np.median(cpm, axis=0)
    frac_zero = (cpm <= 0).mean(axis=0)

    is_protein_coding = np.array(
        [
            (g in gc_map.index) and (gc_map.at[g, "gene_biotype"] == "protein_coding")
            for g in AXIS
        ],
        dtype=bool,
    )
    not_chry = np.array(
        [
            (g in gc_map.index) and (gc_map.at[g, "chromosome"] != "chrY")
            for g in AXIS
        ],
        dtype=bool,
    )
    not_olf_testis = np.array(
        [
            not bool(EXCLUDED_FAMILIES_RE.match(str(gc_map.at[g, "gene_name"] if g in gc_map.index else "")))
            for g in AXIS
        ],
        dtype=bool,
    )

    eligible = (
        (med_cpm >= MIN_MEDIAN_CPM)
        & (frac_zero <= MAX_FRAC_ZERO)
        & is_protein_coding
        & not_chry
        & not_olf_testis
    )
    print(f"POOL_V2 n={V2['n']}, axis={len(AXIS)}, eligible label-free candidates={int(eligible.sum())}", flush=True)

    Xlog = np.log2(cpm + 1.0)
    MARG, _ = sel.marginal_abs_spearman(Xlog, y, coh, cohorts)
    cand_all = np.argsort(-MARG)
    cand = np.sort(cand_all[eligible[cand_all]][:CANDIDATE_M])

    genes, pairs = sel.tsp_panel(cpm, cand, y, coh, cohorts, k=K, redundancy=REDUNDANCY)
    model = pc.fit_logratio(cpm, genes, pairs, y, coh, cohorts)

    # Reconstruct exact pair_mean (m_p) and b0_raw matching panel_core.fit_logratio
    pos_map = {int(g): i for i, g in enumerate(genes)}
    pl = np.array([[pos_map[int(a)], pos_map[int(b)]] for a, b in pairs], dtype=np.int64)
    Xtr, _ = pc.panel_log(cpm[:, genes])
    D = Xtr[:, pl[:, 0]] - Xtr[:, pl[:, 1]]
    pair_mean = D.mean(axis=0).astype(np.float64)
    coef = np.asarray(model["coef"], dtype=np.float64)
    intercept = float(model["intercept"])
    b0_raw = float(intercept + np.dot(coef, pair_mean))

    pair_a = np.array([AXIS[int(a)] for a, _ in pairs])
    pair_b = np.array([AXIS[int(b)] for _, b in pairs])
    sym_a = np.array([str(gc_map.at[g, "gene_name"] if g in gc_map.index else g) for g in pair_a])
    sym_b = np.array([str(gc_map.at[g, "gene_name"] if g in gc_map.index else g) for g in pair_b])
    med_a = np.array([float(med_cpm[int(a)]) for a, _ in pairs], dtype=np.float64)
    med_b = np.array([float(med_cpm[int(b)]) for _, b in pairs], dtype=np.float64)
    fz_a = np.array([float(frac_zero[int(a)]) for a, _ in pairs], dtype=np.float64)
    fz_b = np.array([float(frac_zero[int(b)]) for _, b in pairs], dtype=np.float64)

    wdir = HERE / "weights"
    wdir.mkdir(parents=True, exist_ok=True)
    npz_path = wdir / "masld_panel_v1_1.npz"
    np.savez_compressed(
        npz_path,
        pair_gene_a=pair_a,
        pair_gene_b=pair_b,
        pair_symbol_a=sym_a,
        pair_symbol_b=sym_b,
        coef=coef,
        intercept=np.float64(intercept),
        pair_mean=pair_mean,
        b0_raw=np.float64(b0_raw),
        floor_cpm=np.float64(model["floor_cpm"]),
        formulation=np.array("logratio"),
        interval_calibration=np.array("NONE"),
        lam=np.float64(model["lam"]),
        loco_spearman_in_pool=np.float64(model["loco_spearman_in_pool"]),
        sigma2_eta=np.float64(model["sigma2_eta"]),
        min_median_cpm_filter=np.float64(MIN_MEDIAN_CPM),
        max_frac_zero_filter=np.float64(MAX_FRAC_ZERO),
        train_median_cpm_a=med_a,
        train_median_cpm_b=med_b,
        train_frac_zero_a=fz_a,
        train_frac_zero_b=fz_b,
    )

    pair_rows = []
    for idx_p in range(len(pairs)):
        pair_rows.append(
            {
                "pair_index": idx_p + 1,
                "gene_a": str(pair_a[idx_p]),
                "symbol_a": str(sym_a[idx_p]),
                "median_cpm_a": float(med_a[idx_p]),
                "frac_zero_a": float(fz_a[idx_p]),
                "gene_b": str(pair_b[idx_p]),
                "symbol_b": str(sym_b[idx_p]),
                "median_cpm_b": float(med_b[idx_p]),
                "frac_zero_b": float(fz_b[idx_p]),
                "coef": float(coef[idx_p]),
                "pair_mean": float(pair_mean[idx_p]),
            }
        )
    pd.DataFrame(pair_rows).to_csv(wdir / "panel_pairs.tsv", sep="\t", index=False)

    # Evaluate frozen release artifact via released score.py on external cohorts
    sys.path.insert(0, str(HERE))
    import score as rel_score  # noqa: E402

    panel_obj = rel_score.load_panel(npz_path)
    EVAL = sb.eval_cohorts(AXIS)
    eval_summary = {}
    for name, E in EVAL.items():
        cohort_genes = AXIS[E["axis_idx"]].tolist()
        s_ext, rep = rel_score.score(E["counts"].T, cohort_genes, panel_obj, is_cq=False)
        r = spearman(s_ext, E["y"])
        c_max = tie_ceiling(E["y"])
        eval_summary[name] = {
            "n": int(len(E["y"])),
            "spearman": float(r),
            "tie_ceiling": float(c_max),
            "frac_of_ceiling": float(r / c_max),
            "pair_coverage_min": float(rep["pair_coverage_min"]),
        }
        print(f"  EVAL {name}: n={len(E['y'])}, rho={r:.4f}, ceil={c_max:.4f}, min_cov={rep['pair_coverage_min']:.3f}", flush=True)

    gse268_rho = eval_summary["GSE268273"]["spearman"]
    expected = {
        "check": "GSE268273 external, model frozen, participants",
        "n": eval_summary["GSE268273"]["n"],
        "spearman_vs_kleiner_fibrosis": round(gse268_rho, 6),
        "tolerance": 0.0005,
        "tie_ceiling": round(eval_summary["GSE268273"]["tie_ceiling"], 6),
        "all_eval_cohorts": eval_summary,
    }
    (HERE / "expected.json").write_text(json.dumps(expected, indent=2) + "\n")

    fit_meta = {
        "seed": SEED,
        "pool": "POOL_V2",
        "n_train": int(V2["n"]),
        "cohorts": cohorts,
        "min_median_cpm": MIN_MEDIAN_CPM,
        "max_frac_zero": MAX_FRAC_ZERO,
        "protein_coding_only": True,
        "excluded_families": "olfactory/testis/chrY",
        "n_eligible_candidates": int(eligible.sum()),
        "n_candidates_considered": int(len(cand)),
        "k_pairs": int(len(pairs)),
        "n_unique_genes": int(len(genes)),
        "min_panel_median_cpm_pool_v2": float(min(med_a.min(), med_b.min())),
        "max_panel_frac_zero_pool_v2": float(max(fz_a.max(), fz_b.max())),
        "lam": float(model["lam"]),
        "loco_spearman_in_pool": float(model["loco_spearman_in_pool"]),
        "intercept": intercept,
        "b0_raw": b0_raw,
        "eval_cohorts": eval_summary,
    }
    (wdir / "fit_metadata.json").write_text(json.dumps(fit_meta, indent=2) + "\n")
    print("Built masld-panel-v1.1 weights and metadata.", flush=True)


if __name__ == "__main__":
    main()
