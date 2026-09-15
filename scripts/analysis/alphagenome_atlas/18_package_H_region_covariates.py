#!/usr/bin/env python3
"""Step 18 (Package H): region-level covariate explanation of the released chromatin predictor's behaviour.

Uses the shipped `release/masld-liver-chromatin-state-v1.2` form (RRR shrinkage floor 1e-7, rank 20, fixed +
cis ridge on its residual) and its saved out-of-fold predictions. No architecture search; no Atlas overlay (the
release's coordinate convention is inferred, not sourced, so only region-level covariates enter).

Endpoints (per region, all 96,460 regions first; the two "reliable" masks are supplementary subsets):
  - out-of-fold skill of the shipped form, per fold (the producer's definition: 1 - SSE / SS about the
    TRAINING-fold mean of that fold's own residualised matrix), reproduced against the release deposit;
  - the local increment: shipped form minus its RRR offset on the same folds; the random-gene increment on the
    same folds is the control (lane A pairing);
  - the external ATAC association (GSE296875 donors, marginal transport) as a separate endpoint.
Covariates: measured signal strength (mean log2 CPM), variability (sd), zero-count fraction, promoter (TSS <= 1 kb),
TSS distance, expressed-gene count within 100 kb, GC. Mappability is not available in this repository.

Uncertainty: per-fold values (five participant folds) give the participant-aware spread; a 1-Mb block bootstrap
over regions gives region-level intervals only. Same region set across forms for any form comparison.

Outputs (tables/): package_H_predictions.json (written before any result), region_covariates.tsv.gz,
region_covariate_associations.tsv, region_covariate_strata.tsv, region_covariate_model.tsv,
region_form_comparison.tsv, reliable_set_covariates.tsv, package_H_results.json
"""

from __future__ import annotations

import json
import math
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pysam
from scipy.stats import spearmanr

import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
BENCH = la.PROJECT / "Analysis/MASLD_Model_Benchmark"
RELEASE = BENCH / "release/masld-liver-chromatin-state-v1.2/weights/chromatin_state_v1_2.npz"
CIS_LANE = BENCH / "executions/chromatin-cis-head-20260908T102817Z/out/cis_head_oof.npz"
STABLE_LANE = BENCH / "executions/chromatin-stable-rrr-forms-20260908T192024Z"
OFFSET_LANE = BENCH / "executions/chromatin-rrr-offset-cis-20260908T173402Z/out/offset_oof.npz"
TRANSFER = STABLE_LANE / "out_finish/stable_transfer_rho.npz"
LANE0 = BENCH / "executions/chromatin-beyond-histology-20260907T161521Z"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
SEED = 20260909
DRAWS = 2000
NPROC = 8
COVARIATES = ["gc", "log10_dist_tss", "n_cis_100k", "signal_mean", "signal_sd", "zero_fraction", "promoter"]
ENDPOINTS = ["skill_shipped_mean", "skill_states_cis_mean", "increment_local_mean", "increment_random_mean",
             "transfer_rho_marginal", "transfer_rho_raw"]

PREDICTIONS = {
    "written_before": "any array beyond shapes was read",
    "H1_locality": "the local increment (shipped minus RRR offset) rises with n_cis_100k (block-bootstrap Spearman lower bound > 0) and falls with TSS distance; the random-gene increment shows |rho| < 0.03 with n_cis_100k",
    "H2_variability": "shipped-form skill rises with region sd of log2 CPM (rho > 0.10) and with mean signal (rho > 0)",
    "H3_gc": "GC explains little: |rho| < 0.05 with shipped-form skill",
    "H4_promoter": "promoter regions have LOWER shipped-form skill than distal regions (median difference < 0), consistent with the release's distal reliable set",
    "H5_transfer": "the external ATAC rho rises with in-cohort skill (rho about +0.20, the release's +0.204 reproduced within 0.02) and shares the sign pattern of H1-H4",
    "H6_model_r2": "the seven covariates together explain between 0.02 and 0.15 of the variance in shipped-form skill",
}


# ---------------------------------------------------------------- helpers (tested)
def parse_region_key(key: str) -> tuple[str, int, int]:
    chrom, rng_ = key.split(":")
    if not chrom.startswith("chr"):
        raise ValueError(f"region key without chr prefix: {key}")
    start, end = rng_.split("-")
    return chrom, int(start), int(end)


def mb_block(chrom: str, start: int, end: int) -> str:
    return f"{chrom}:{int(((start + end) / 2.0) // 1_000_000)}"


def cis_bin(k: int) -> str:
    if k <= 2:
        return str(int(k))
    if k <= 4:
        return "3-4"
    if k <= 9:
        return "5-9"
    return "10+"


def per_fold_region_skill(target: np.ndarray, pred: np.ndarray, fold: np.ndarray, train_mean: dict) -> np.ndarray:
    """Per fold f, per region: 1 - sum_te (target - pred)^2 / sum_te (target - train_mean[f])^2; NaN where the denominator is ~0."""
    folds = sorted(set(int(f) for f in fold))
    out = np.full((len(folds), target.shape[1]), np.nan)
    for i, f in enumerate(folds):
        te = np.where(fold == f)[0]
        num = ((target[te] - pred[te]) ** 2).sum(0)
        den = ((target[te] - np.asarray(train_mean[f])[None, :]) ** 2).sum(0)
        out[i] = 1 - num / np.where(den < 1e-12, np.nan, den)
    return out


def block_bootstrap_spearman(x: np.ndarray, y: np.ndarray, blocks: np.ndarray, draws: int, seed: int) -> dict:
    """Spearman rho with a percentile interval from resampling BLOCKS (1-Mb) with replacement, rows within a block kept together."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, blocks = x[ok], y[ok], np.asarray(blocks)[ok]
    rho = float(spearmanr(x, y).correlation)
    uniq, inv = np.unique(blocks, return_inverse=True)
    members = [np.where(inv == b)[0] for b in range(len(uniq))]
    rng = np.random.default_rng(seed)
    vals = np.empty(draws)
    for d in range(draws):
        pick = rng.integers(0, len(uniq), len(uniq))
        idx = np.concatenate([members[b] for b in pick])
        vals[d] = spearmanr(x[idx], y[idx]).correlation
    return {"rho": rho, "ci_low": float(np.nanpercentile(vals, 2.5)), "ci_high": float(np.nanpercentile(vals, 97.5)),
            "n": int(ok.sum()), "n_blocks": int(len(uniq)), "draws": draws}


def covariate_model(X: np.ndarray, y: np.ndarray, names: list[str]) -> dict:
    """OLS of y on standardised columns; R2 and drop-one partial R2 = (SSE_reduced - SSE_full) / SSE_reduced (0 when the reduced fit is already exact)."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    ok = np.isfinite(y) & np.isfinite(X).all(1)
    X, y = X[ok], y[ok]
    sd = X.std(0); sd[sd == 0] = 1.0
    Z = (X - X.mean(0)) / sd
    A = np.column_stack([np.ones(len(y)), Z])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    sse_full = float(((y - A @ beta) ** 2).sum()); sst = float(((y - y.mean()) ** 2).sum())
    partial, delta_r2 = {}, {}
    for j, name in enumerate(names):
        keep = [c for c in range(A.shape[1]) if c != j + 1]
        b_red, *_ = np.linalg.lstsq(A[:, keep], y, rcond=None)
        sse_red = float(((y - A[:, keep] @ b_red) ** 2).sum())
        partial[name] = 0.0 if sse_red <= 1e-12 * max(sst, 1.0) else (sse_red - sse_full) / sse_red
        delta_r2[name] = (sse_red - sse_full) / sst if sst > 0 else math.nan
    return {"n": int(ok.sum()), "r2": 1 - sse_full / sst if sst > 0 else math.nan, "coef_standardised": dict(zip(names, beta[1:].tolist())),
            "partial_r2": partial, "delta_r2": delta_r2}


def conc_block(raw: np.ndarray) -> np.ndarray:
    """The producer's H3K27ac concentration descriptors (gini, entropy, IQR, top-5% share), copied verbatim from the discrepancy-map lane."""
    lib = raw.sum(1); ncol = raw.shape[1]; Pm = raw / lib[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -np.nansum(np.where(Pm > 0, Pm * np.log(Pm), 0.0), axis=1)
    sa = np.sort(raw, axis=1); gini = ((2 * np.arange(1, ncol + 1) - ncol - 1) * sa).sum(1) / (ncol * sa.sum(1))
    iqr = np.array([float(np.subtract(*np.percentile(np.log2(r[r > 0]), [75, 25]))) for r in raw])
    kt = int(round(0.05 * ncol)); top5 = (-np.sort(-raw, axis=1))[:, :kt].sum(1) / lib
    return np.column_stack([gini, ent, iqr, top5])


def _assoc(args):
    endpoint, covariate, subset, x, y, blocks = args
    r = block_bootstrap_spearman(x, y, blocks, DRAWS, SEED)
    return {"endpoint": endpoint, "covariate": covariate, "region_set": subset, **r}


# ---------------------------------------------------------------- main
def main() -> None:
    pred_path = TABLES / "package_H_predictions.json"
    if pred_path.exists():
        if json.load(pred_path.open()) != PREDICTIONS:
            raise la.ContractError(f"refusing to overwrite a different prespecification: {pred_path}")
    else:
        json.dump(PREDICTIONS, pred_path.open("w"), indent=1)

    rel = np.load(RELEASE, allow_pickle=True)
    cis = np.load(CIS_LANE, allow_pickle=True)
    stable = np.load(STABLE_LANE / "out/stable_oof.npz", allow_pickle=True)
    stable_res = json.load((STABLE_LANE / "out/stable_results.json").open())
    offset = np.load(OFFSET_LANE, allow_pickle=True)
    keys = np.asarray(rel["prof_region_key"]).astype(str)
    n_reg = len(keys)
    if not (np.array_equal(keys, np.asarray(cis["region_key"]).astype(str)) and np.array_equal(keys, np.asarray(rel["chrom_region_axis"]).astype(str))):
        raise la.ContractError("region axes differ between the release and the cis-head deposit")
    fold = np.asarray(cis["fold"], int)
    Cres = np.asarray(cis["Cres_oof"], np.float64)
    for name, z in (("stable", stable), ("offset", offset)):
        if np.abs(np.asarray(z["Cres_oof"], np.float64) - Cres).max() > 1e-6 or not np.array_equal(np.asarray(z["fold"], int), fold):
            raise la.ContractError(f"{name} lane deposit does not share the cis-head residual target or folds")
    la.log(f"deposits aligned: {n_reg} regions, {len(fold)} participants, folds {sorted(set(fold.tolist()))}")

    # Training-fold baselines rebuilt exactly as the producer defines them (fold-specific descriptor residualisation).
    sys.path.insert(0, str(LANE0))
    import common as C0  # noqa: E402
    C0.OUT = ROOT / "logs" / "package_H_substrate"; C0.OUT.mkdir(parents=True, exist_ok=True)
    S = C0.load_substrate()
    if not np.array_equal(np.asarray(S["fold"], int), fold):
        raise la.ContractError("substrate fold vector differs from the deposit")
    h3_raw = np.asarray(S["h3_raw"], np.float64)
    Craw = C0.logcpm(h3_raw)
    CONC = conc_block(h3_raw)
    Zi = np.column_stack([np.ones(len(Craw)), CONC])
    train_mean = {}
    for f in sorted(set(fold.tolist())):
        tr, te = np.where(fold != f)[0], np.where(fold == f)[0]
        b, *_ = np.linalg.lstsq(Zi[tr], Craw[tr], rcond=None)
        Cf = Craw - Zi @ b
        if not np.allclose(Cf[te], Cres[te], atol=1e-6, rtol=0):
            raise la.ContractError(f"fold {f}: recomputed residual != deposited Cres_oof")
        train_mean[f] = Cf[tr].mean(0)

    P_ship = np.asarray(stable["pred_rrr_offset_cis_wide"], np.float64)
    P_rrr = np.asarray(stable["pred_rrr_wide_stable"], np.float64)
    P_rrrA = np.asarray(offset["pred_rrr"], np.float64)
    P_cisA = np.asarray(offset["pred_rrr_offset_cis"], np.float64)
    P_randA = np.asarray(offset["pred_rrr_offset_random"], np.float64)
    TM = np.vstack([train_mean[int(f)] for f in fold])
    whole = 1 - ((Cres - P_ship) ** 2).sum() / ((Cres - TM) ** 2).sum()
    g1 = abs(whole - stable_res["in_cohort"]["skill"]["rrr_offset_cis_wide"]["skill"])
    if g1 > 1e-5:
        raise la.ContractError(f"G1: shipped-form whole-matrix skill not reproduced (|diff| {g1:.2e})")
    sk_ship = per_fold_region_skill(Cres, P_ship, fold, train_mean)
    g2 = float(np.nanmax(np.abs(sk_ship.astype(np.float32) - np.asarray(rel["prof_per_fold_skill_offset_form"], np.float32))))
    if g2 > 1e-5:
        raise la.ContractError(f"G2: per-fold per-region skill of the shipped form not reproduced (max |diff| {g2:.2e})")
    la.log(f"G1 whole-matrix skill {whole:+.4f} reproduced (|diff| {g1:.1e}); G2 per-region per-fold skill reproduced (max |diff| {g2:.1e})")
    sk_rrr = per_fold_region_skill(Cres, P_rrr, fold, train_mean)
    sk_rrrA = per_fold_region_skill(Cres, P_rrrA, fold, train_mean)
    sk_cisA = per_fold_region_skill(Cres, P_cisA, fold, train_mean)
    sk_randA = per_fold_region_skill(Cres, P_randA, fold, train_mean)
    sk_sc = np.asarray(rel["prof_per_fold_skill"], np.float64)
    rel_ship = np.asarray(rel["prof_reliable_offset_form"], bool); rel_sc = np.asarray(rel["prof_reliable"], bool)
    if int(rel_ship.sum()) != int((sk_ship.min(0) > 0.2).sum()):
        raise la.ContractError("shipped-form reliable mask does not equal (per-fold skill > 0.2 in all folds)")

    # Region covariates.
    parsed = [parse_region_key(k) for k in keys]
    chrom = np.array([p[0] for p in parsed]); start = np.array([p[1] for p in parsed]); end = np.array([p[2] for p in parsed])
    fasta = pysam.FastaFile(FASTA)
    gc = np.empty(n_reg)
    for j in range(n_reg):
        seq = fasta.fetch(chrom[j], start[j], end[j]).upper()
        gc[j] = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
    dist = np.asarray(cis["dist_nearest_tss"], float)
    frame = pd.DataFrame({
        "region_key": keys, "chrom": chrom, "start": start, "end": end, "width": end - start,
        "autosomal": ~np.isin(chrom, ["chrX", "chrY", "chrM"]), "mb_block": [mb_block(c, s, e) for c, s, e in parsed],
        "gc": gc, "promoter": np.asarray(cis["promoter"], bool).astype(int), "dist_nearest_tss": dist, "log10_dist_tss": np.log10(dist + 1.0),
        "n_cis_100k": np.asarray(cis["n_cis_100k"], int), "cis_bin": [cis_bin(k) for k in np.asarray(cis["n_cis_100k"], int)],
        "signal_mean": Craw.mean(0), "signal_sd": Craw.std(0, ddof=1), "zero_fraction": (h3_raw == 0).mean(0),
        "skill_shipped_mean": sk_ship.mean(0), "skill_shipped_min": sk_ship.min(0), "skill_states_cis_mean": sk_sc.mean(0),
        "skill_rrr_offset_mean": sk_rrr.mean(0), "increment_local_mean": (sk_ship - sk_rrr).mean(0),
        "increment_local_laneA_mean": (sk_cisA - sk_rrrA).mean(0), "increment_random_mean": (sk_randA - sk_rrrA).mean(0),
        "reliable_shipped_form": rel_ship, "reliable_states_cis": rel_sc,
        "transfer_rho_marginal": np.nan, "transfer_rho_raw": np.nan, "transfer_lane_reliable_flag": False, "in_transfer_set": False,
    })
    for i in range(sk_ship.shape[0]):
        frame[f"skill_shipped_fold{i}"] = sk_ship[i]; frame[f"increment_local_fold{i}"] = sk_ship[i] - sk_rrr[i]
    tr_ = np.load(TRANSFER, allow_pickle=True)
    ridx = np.asarray(tr_["region_index"], int)
    frame.loc[ridx, "transfer_rho_marginal"] = np.asarray(tr_["rho_marginal_rrr_offset_cis"], float)
    frame.loc[ridx, "transfer_rho_raw"] = np.asarray(tr_["rho_raw_rrr_offset_cis"], float)
    frame.loc[ridx, "transfer_lane_reliable_flag"] = np.asarray(tr_["reliable"], bool)
    frame.loc[ridx, "in_transfer_set"] = True
    frame.to_csv(TABLES / "region_covariates.tsv.gz", sep="\t", index=False)
    la.log(f"region table written: {n_reg} regions, {int(frame.in_transfer_set.sum())} in the transfer set, GC median {np.median(gc):.3f}")

    # Associations: block-bootstrap Spearman, all autosomal regions first; reliable subsets supplementary.
    subsets = {"all_autosomal": frame.autosomal.values, "reliable_shipped_form": frame.autosomal.values & rel_ship,
               "reliable_states_cis": frame.autosomal.values & rel_sc, "reliable_union": frame.autosomal.values & (rel_ship | rel_sc)}
    jobs = []
    for subset, mask in subsets.items():
        sub = frame[mask]
        for endpoint in ENDPOINTS:
            for cov in COVARIATES:
                jobs.append((endpoint, cov, subset, sub[endpoint].values, sub[cov].values, sub["mb_block"].values))
        jobs.append(("transfer_rho_marginal", "skill_shipped_mean", subset, sub["transfer_rho_marginal"].values, sub["skill_shipped_mean"].values, sub["mb_block"].values))
        jobs.append(("transfer_rho_raw", "skill_shipped_mean", subset, sub["transfer_rho_raw"].values, sub["skill_shipped_mean"].values, sub["mb_block"].values))
    with Pool(NPROC) as pool:
        assoc = pool.map(_assoc, jobs, chunksize=1)
    # per-fold spread for fold-resolved endpoints (participant-aware)
    auto = frame.autosomal.values
    per_fold = {}
    for endpoint, arr in (("skill_shipped_mean", sk_ship), ("increment_local_mean", sk_ship - sk_rrr), ("increment_random_mean", sk_randA - sk_rrrA)):
        for cov in COVARIATES:
            vals = [float(spearmanr(arr[i][auto], frame[cov].values[auto], nan_policy="omit").correlation) for i in range(arr.shape[0])]
            per_fold[(endpoint, cov)] = vals
    for r in assoc:
        v = per_fold.get((r["endpoint"], r["covariate"])) if r["region_set"] == "all_autosomal" else None
        r["fold_rho_min"] = min(v) if v else ""; r["fold_rho_max"] = max(v) if v else ""; r["fold_rho_values"] = ";".join(f"{x:.4f}" for x in v) if v else ""
        r["uncertainty"] = "block bootstrap over 1-Mb regions (region-level only); fold columns give the participant-fold spread"
    la.write_tsv_once(TABLES / "region_covariate_associations.tsv", assoc, list(assoc[0].keys()))

    # Strata table: endpoint medians by covariate decile / bin.
    strata = []
    sub = frame[auto].copy()
    for cov in COVARIATES:
        if cov == "promoter":
            sub["_stratum"] = np.where(sub.promoter == 1, "promoter(<=1kb)", "distal")
        elif cov == "n_cis_100k":
            sub["_stratum"] = sub.cis_bin
        else:
            sub["_stratum"] = pd.qcut(sub[cov].rank(method="first"), 10, labels=[f"decile{i + 1}" for i in range(10)]).astype(str)
        for stratum, d in sub.groupby("_stratum", sort=False):
            row = {"covariate": cov, "stratum": stratum, "n_regions": len(d), "covariate_median": float(d[cov].median()),
                   "fraction_reliable_shipped": float(d.reliable_shipped_form.mean()), "fraction_reliable_states_cis": float(d.reliable_states_cis.mean())}
            for endpoint in ENDPOINTS:
                q = d[endpoint].quantile([0.25, 0.5, 0.75])
                row[f"{endpoint}_median"] = float(q[0.5]); row[f"{endpoint}_iqr"] = f"{q[0.25]:.4f};{q[0.75]:.4f}"; row[f"{endpoint}_n"] = int(d[endpoint].notna().sum())
            strata.append(row)
    la.write_tsv_once(TABLES / "region_covariate_strata.tsv", strata, list(strata[0].keys()))

    # Covariate model (OLS on standardised covariates), per endpoint; per fold for fold-resolved endpoints.
    model_rows = []
    X = frame.loc[auto, COVARIATES].values
    for endpoint in ENDPOINTS:
        m = covariate_model(X, frame.loc[auto, endpoint].values, COVARIATES)
        model_rows.append({"endpoint": endpoint, "fold": "mean_over_folds" if "fold" not in endpoint else "", "n": m["n"], "r2": m["r2"],
                           **{f"coef_{k}": v for k, v in m["coef_standardised"].items()}, **{f"partial_r2_{k}": v for k, v in m["partial_r2"].items()}})
    for i in range(sk_ship.shape[0]):
        for endpoint, arr in (("skill_shipped", sk_ship), ("increment_local", sk_ship - sk_rrr)):
            m = covariate_model(X, arr[i][auto], COVARIATES)
            model_rows.append({"endpoint": endpoint, "fold": str(i), "n": m["n"], "r2": m["r2"],
                               **{f"coef_{k}": v for k, v in m["coef_standardised"].items()}, **{f"partial_r2_{k}": v for k, v in m["partial_r2"].items()}})
    la.write_tsv_once(TABLES / "region_covariate_model.tsv", model_rows, list(model_rows[0].keys()))

    # Same-region-set form comparison: shipped vs states + cis, paired per region, by stratum.
    comp = []
    diff = sk_ship - sk_sc
    for set_name, mask in (("all_autosomal", auto), ("reliable_union_fixed_subset", auto & (rel_ship | rel_sc))):
        d = frame[mask].assign(_diff=diff.mean(0)[mask])
        groups = [("all", d)] + [(f"promoter={k}", g) for k, g in d.groupby("promoter")] + [(f"cis_bin={k}", g) for k, g in d.groupby("cis_bin")] \
            + [(f"signal_sd_quartile={k}", g) for k, g in d.groupby(pd.qcut(d.signal_sd.rank(method="first"), 4, labels=["q1", "q2", "q3", "q4"]))]
        for label, g in groups:
            fold_means = [float(np.nanmean(diff[i][g.index.values])) for i in range(diff.shape[0])]
            comp.append({"region_set": set_name, "stratum": label, "n_regions": len(g), "mean_skill_shipped": float(g.skill_shipped_mean.mean()),
                         "mean_skill_states_cis": float(g.skill_states_cis_mean.mean()), "mean_paired_difference": float(g._diff.mean()),
                         "fraction_shipped_higher": float((g._diff > 0).mean()), "fold_means": ";".join(f"{x:+.4f}" for x in fold_means),
                         "fold_min": min(fold_means), "fold_max": max(fold_means), "note": "same regions for both forms; paired per region and per fold"})
    la.write_tsv_once(TABLES / "region_form_comparison.tsv", comp, list(comp[0].keys()))

    # Reliable-set covariate profile (supplementary).
    prof = []
    for set_name, mask in (("all", np.ones(n_reg, bool)), ("all_autosomal", auto), ("reliable_shipped_form", rel_ship), ("reliable_states_cis", rel_sc), ("reliable_both", rel_ship & rel_sc)):
        d = frame[mask]
        prof.append({"region_set": set_name, "n_regions": len(d), "promoter_fraction": float(d.promoter.mean()), "median_dist_tss": float(d.dist_nearest_tss.median()),
                     "median_n_cis_100k": float(d.n_cis_100k.median()), "median_gc": float(d.gc.median()), "median_signal_mean": float(d.signal_mean.median()),
                     "median_signal_sd": float(d.signal_sd.median()), "median_zero_fraction": float(d.zero_fraction.median()),
                     "mean_skill_shipped": float(d.skill_shipped_mean.mean()), "mean_increment_local": float(d.increment_local_mean.mean()),
                     "mean_transfer_rho_marginal": float(d.transfer_rho_marginal.mean()), "n_in_transfer_set": int(d.in_transfer_set.sum())})
    la.write_tsv_once(TABLES / "reliable_set_covariates.tsv", prof, list(prof[0].keys()))

    # Score the written predictions and record the model-card number.
    A = {(r["endpoint"], r["covariate"]): r for r in assoc if r["region_set"] == "all_autosomal"}
    prom_med = frame[auto].groupby("promoter").skill_shipped_mean.median()
    res = {
        "guards": {"G1_whole_matrix_skill_abs_diff": g1, "G2_per_region_skill_max_abs_diff": g2, "whole_matrix_skill": float(whole)},
        "n_regions": n_reg, "n_autosomal": int(auto.sum()), "n_reliable_shipped_form": int(rel_ship.sum()), "n_reliable_states_cis": int(rel_sc.sum()),
        "n_reliable_both": int((rel_ship & rel_sc).sum()),
        "model_card_number": {"headline": 4847, "shipped_form_mask": int(rel_ship.sum()), "note": "README/MODEL_CARD headline the 4,847 states+cis set; the shipped form's own mask is the smaller number"},
        "predictions_scored": {
            "H1_locality": bool(A[("increment_local_mean", "n_cis_100k")]["ci_low"] > 0 and A[("increment_local_mean", "log10_dist_tss")]["rho"] < 0 and abs(A[("increment_random_mean", "n_cis_100k")]["rho"]) < 0.03),
            "H2_variability": bool(A[("skill_shipped_mean", "signal_sd")]["rho"] > 0.10 and A[("skill_shipped_mean", "signal_mean")]["rho"] > 0),
            "H3_gc": bool(abs(A[("skill_shipped_mean", "gc")]["rho"]) < 0.05),
            "H4_promoter": bool(prom_med.get(1, math.nan) < prom_med.get(0, math.nan)),
            "H5_transfer": bool(abs(A[("transfer_rho_marginal", "skill_shipped_mean")]["rho"] - 0.204) < 0.02),
            "H6_model_r2": bool(0.02 <= model_rows[0]["r2"] <= 0.15),
        },
        "key_values": {f"{k[0]}~{k[1]}": {"rho": A[k]["rho"], "ci": [A[k]["ci_low"], A[k]["ci_high"]]} for k in A} | {"promoter_median_skill": {str(k): float(v) for k, v in prom_med.items()}, "model_r2_shipped": model_rows[0]["r2"]},
        "limits": ["a static region annotation cannot explain donor differences", "an RNA-to-chromatin model cannot be inverted into an enhancer-to-gene model",
                   "region resampling gives region-level intervals only; the five fold values are the participant-aware spread", "mappability is not available in this repository",
                   "the release's region coordinate convention is inferred, so GC is computed on the keyed interval as written (a 1-bp offset is immaterial at this width)"],
    }
    json.dump(res, (TABLES / "package_H_results.json").open("w"), indent=1)
    la.log(f"Package H: predictions scored {res['predictions_scored']}; model R2 {model_rows[0]['r2']:.4f}")


if __name__ == "__main__":
    main()
