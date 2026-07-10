#!/usr/bin/env python
"""
28d_score_benchmark.py  --  Pseudobulk benchmark, step 4/4 (scoring).

Loads truth_proportions.tsv and the method estimate tables (rectangle / music /
instaprism), aligns donors + the 16 canonical cell types, and computes per-cell-type
and overall accuracy vs the known simulated proportions:

  RMSE, Pearson r, Spearman rho, MAE, and spillover (mean mass assigned to truly-absent
  cell types, truth == 0).

Primary scope = RAW estimates reindexed to the 16 canonical types (unassigned mass such
as Rectangle's 'Unknown' bucket is dropped and therefore counts as error -- honest). A
secondary 'overall_renorm' scope renormalizes each method's 16 types to sum 1 for an
apples-to-apples simplex comparison.

Writes:
  $BENCH_OUT_DIR/../pseudobulk_benchmark_scores.csv   (long: method, scope, metric, value)
  $BENCH_OUT_DIR/../pseudobulk_benchmark_summary.txt   (ranked human-readable summary)

Runs in the `rectangle` env python.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR = os.environ.get("BENCH_OUT_DIR",
                         os.path.join(ROOT, "Analysis/Deconvolution/rectangle_comparison/pseudobulk"))
SCORE_DIR = os.path.dirname(OUT_DIR.rstrip("/"))          # rectangle_comparison/

METHOD_FILES = {
    "Rectangle":  "rectangle_estimates.tsv",
    "MuSiC":      "music_estimates.tsv",
    "InstaPrism": "instaprism_estimates.tsv",
}


def safe_corr(a, b, kind):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    r = pearsonr(a, b)[0] if kind == "pearson" else spearmanr(a, b)[0]
    return float(r)


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a, float) - np.asarray(b, float)) ** 2)))


def mae(a, b):
    return float(np.mean(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def main():
    truth = pd.read_csv(os.path.join(OUT_DIR, "truth_proportions.tsv"), sep="\t",
                        index_col=0)
    truth.index = truth.index.astype(str)
    canonical = [l.strip() for l in
                 open(os.path.join(OUT_DIR, "canonical_celltypes.txt")) if l.strip()]
    truth = truth.reindex(columns=canonical)

    rows = []
    for method, fname in METHOD_FILES.items():
        path = os.path.join(OUT_DIR, fname)
        if not os.path.exists(path):
            print(f"[28d] SKIP {method}: {fname} not found", flush=True)
            continue
        est = pd.read_csv(path, sep="\t", index_col=0)
        est.index = est.index.astype(str)
        est = est.reindex(index=truth.index, columns=canonical, fill_value=0.0)

        T = truth.values
        E = est.values
        keep = ~np.isnan(T).any(axis=1)                  # donors with valid truth
        T, E = T[keep], E[keep]
        n = T.shape[0]

        # ---- per cell type ----
        for j, ct in enumerate(canonical):
            rows += [
                dict(method=method, scope=ct, metric="rmse",       value=rmse(T[:, j], E[:, j])),
                dict(method=method, scope=ct, metric="pearson",    value=safe_corr(T[:, j], E[:, j], "pearson")),
                dict(method=method, scope=ct, metric="spearman",   value=safe_corr(T[:, j], E[:, j], "spearman")),
                dict(method=method, scope=ct, metric="mae",        value=mae(T[:, j], E[:, j])),
                dict(method=method, scope=ct, metric="mean_truth", value=float(T[:, j].mean())),
                dict(method=method, scope=ct, metric="mean_est",   value=float(E[:, j].mean())),
            ]

        # ---- overall (raw, flattened over donor x type) ----
        tf, ef = T.ravel(), E.ravel()
        absent = T == 0
        rows += [
            dict(method=method, scope="overall", metric="rmse",      value=rmse(tf, ef)),
            dict(method=method, scope="overall", metric="pearson",   value=safe_corr(tf, ef, "pearson")),
            dict(method=method, scope="overall", metric="spearman",  value=safe_corr(tf, ef, "spearman")),
            dict(method=method, scope="overall", metric="mae",       value=mae(tf, ef)),
            dict(method=method, scope="overall", metric="spillover", value=float(E[absent].mean()) if absent.any() else np.nan),
            dict(method=method, scope="overall", metric="assigned_mass_mean", value=float(E.sum(axis=1).mean())),
            dict(method=method, scope="overall", metric="n_donors",  value=float(n)),
        ]

        # ---- overall renormalized (rows -> simplex over 16 types) ----
        Er = E / np.where(E.sum(axis=1, keepdims=True) == 0, np.nan, E.sum(axis=1, keepdims=True))
        Er = np.nan_to_num(Er)
        trf, erf = T.ravel(), Er.ravel()
        rows += [
            dict(method=method, scope="overall_renorm", metric="rmse",     value=rmse(trf, erf)),
            dict(method=method, scope="overall_renorm", metric="pearson",  value=safe_corr(trf, erf, "pearson")),
            dict(method=method, scope="overall_renorm", metric="spearman", value=safe_corr(trf, erf, "spearman")),
            dict(method=method, scope="overall_renorm", metric="mae",      value=mae(trf, erf)),
        ]
        print(f"[28d] scored {method}: n_donors={n}", flush=True)

    scores = pd.DataFrame(rows)
    out_csv = os.path.join(SCORE_DIR, "pseudobulk_benchmark_scores.csv")
    scores.to_csv(out_csv, index=False)
    print(f"[28d] wrote {out_csv} ({len(scores)} rows)", flush=True)

    # -------------------------------------------------------------- summary.txt
    piv = (scores[scores.scope == "overall"]
           .pivot(index="method", columns="metric", values="value"))
    methods = [m for m in METHOD_FILES if m in piv.index]
    piv = piv.loc[methods]
    ranked = piv.sort_values("rmse")

    lines = []
    lines.append("=" * 78)
    lines.append("PSEUDOBULK GROUND-TRUTH DECONVOLUTION BENCHMARK")
    lines.append("=" * 78)
    lines.append(f"output dir : {OUT_DIR}")
    n_donors = int(piv["n_donors"].iloc[0]) if "n_donors" in piv else 0
    lines.append(f"held-out donors scored : {n_donors}")
    lines.append(f"cell types : {len(canonical)}")
    lines.append("")
    lines.append("OVERALL (raw estimates vs known proportions; lower RMSE / higher r = better)")
    lines.append("-" * 78)
    lines.append(f"{'method':<12}{'RMSE':>10}{'Pearson':>10}{'Spearman':>10}"
                 f"{'MAE':>10}{'spillover':>11}")
    for m in ranked.index:
        r = ranked.loc[m]
        lines.append(f"{m:<12}{r.get('rmse', np.nan):>10.4f}{r.get('pearson', np.nan):>10.4f}"
                     f"{r.get('spearman', np.nan):>10.4f}{r.get('mae', np.nan):>10.4f}"
                     f"{r.get('spillover', np.nan):>11.4f}")
    lines.append("")

    pivr = (scores[scores.scope == "overall_renorm"]
            .pivot(index="method", columns="metric", values="value")).reindex(methods)
    lines.append("OVERALL (renormalized to 16-type simplex; apples-to-apples)")
    lines.append("-" * 78)
    lines.append(f"{'method':<12}{'RMSE':>10}{'Pearson':>10}{'Spearman':>10}{'MAE':>10}")
    for m in pivr.sort_values("rmse").index:
        r = pivr.loc[m]
        lines.append(f"{m:<12}{r.get('rmse', np.nan):>10.4f}{r.get('pearson', np.nan):>10.4f}"
                     f"{r.get('spearman', np.nan):>10.4f}{r.get('mae', np.nan):>10.4f}")
    lines.append("")

    winner = ranked.index[0]
    lines.append(f"WINNER (lowest overall raw RMSE): {winner}")
    lines.append("")

    # per-cell-type RMSE table
    ct_rmse = (scores[(scores.metric == "rmse") & (scores.scope.isin(canonical))]
               .pivot(index="scope", columns="method", values="value")
               .reindex(index=canonical, columns=methods))
    lines.append("PER-CELL-TYPE RMSE (lower = better)")
    lines.append("-" * 78)
    hdr = f"{'cell_type':<26}" + "".join(f"{m:>12}" for m in methods)
    lines.append(hdr)
    for ct in canonical:
        row = ct_rmse.loc[ct]
        lines.append(f"{ct:<26}" + "".join(f"{row.get(m, np.nan):>12.4f}" for m in methods))
    lines.append("=" * 78)

    out_txt = os.path.join(SCORE_DIR, "pseudobulk_benchmark_summary.txt")
    with open(out_txt, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"[28d] wrote {out_txt}", flush=True)
    print("[28d] DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
