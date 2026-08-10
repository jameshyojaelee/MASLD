#!/usr/bin/env python3
"""Independent Currin calibration with 1-Mb block bootstrap uncertainty."""

from __future__ import annotations

import argparse
import bisect
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
BASE = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2"

LABEL = SCORE = BETA = BLOCK_INDEX = None


def peak_index(path: Path):
    raw = {}
    with path.open() as handle:
        for line in handle:
            chrom, start, end = line.split("\t")[:3]
            raw.setdefault(chrom, []).append((int(start), int(end)))
    out = {}
    for chrom, intervals in raw.items():
        merged = []
        for start, end in sorted(intervals):
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        out[chrom] = ([x[0] for x in merged], [x[1] for x in merged])
    return out


def in_peak(index, chrom, pos1):
    if chrom not in index:
        return 0
    starts, ends = index[chrom]
    point = int(pos1) - 1
    i = bisect.bisect_right(starts, point) - 1
    return int(i >= 0 and point < ends[i])


def auc_ap(y, score):
    if len(np.unique(y)) < 2:
        return math.nan, math.nan
    return float(roc_auc_score(y, score)), float(average_precision_score(y, score))


def one_bootstrap(seed):
    rng = np.random.default_rng(seed)
    keys = list(BLOCK_INDEX)
    sampled = rng.choice(len(keys), size=len(keys), replace=True)
    idx = np.concatenate([BLOCK_INDEX[keys[i]] for i in sampled])
    y = LABEL[idx]
    result = {}
    for name, values in SCORE.items():
        if name == "signed_logfc":
            continue
        auc, ap = auc_ap(y, values[idx])
        result[f"{name}_auroc"] = auc
        result[f"{name}_auprc"] = ap
    result["ips_minus_peak_overlap_auroc"] = result["ips_auroc"] - result["peak_overlap_auroc"]
    oriented = idx[(y == 1) & np.isfinite(BETA[idx])]
    if len(oriented) > 2:
        result["direction_spearman"] = float(spearmanr(SCORE["signed_logfc"][oriented], BETA[oriented]).statistic)
        result["direction_sign_concordance"] = float(np.mean(np.sign(SCORE["signed_logfc"][oriented]) == np.sign(BETA[oriented])))
    else:
        result["direction_spearman"] = math.nan
        result["direction_sign_concordance"] = math.nan
    return result


def ci(values):
    z = np.asarray([x for x in values if np.isfinite(x)])
    return {"lower_95": float(np.quantile(z, 0.025)), "upper_95": float(np.quantile(z, 0.975)), "n_bootstrap_valid": len(z)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstraps", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "32")))
    ap.add_argument("--base", type=Path, default=BASE)
    ap.add_argument("--calibration-dir", type=Path)
    ap.add_argument("--scores-dir", type=Path)
    ap.add_argument("--peak-file", type=Path)
    args = ap.parse_args()
    base = args.base.resolve()
    out = args.calibration_dir.resolve() if args.calibration_dir else base / "currin_calibration"
    scores_dir = args.scores_dir.resolve() if args.scores_dir else out / "scores"
    peak_file = args.peak_file.resolve() if args.peak_file else base / "hepatocyte.idr_pooled_summit.narrowPeak"
    meta = pd.read_csv(out / "currin_scoring_metadata.tsv.gz", sep="\t")
    score_frames = []
    for fold in range(5):
        path = scores_dir / f"fold{fold}.variant_scores.tsv"
        frame = pd.read_csv(path, sep="\t")
        frame["heldout_fold_scored"] = fold
        score_frames.append(frame)
    scores = pd.concat(score_frames, ignore_index=True)
    if scores["variant_id"].duplicated().any():
        raise SystemExit("Currin score output contains duplicate variants")
    data = meta.merge(scores, left_on="variant_id_hg38", right_on="variant_id", how="left", validate="one_to_one")
    if data["logfc"].isna().any() or len(data) != len(meta):
        raise SystemExit(f"Currin score completeness failure: {data['logfc'].isna().sum()} missing")
    if not (data["heldout_fold"] == data["heldout_fold_scored"]).all():
        raise SystemExit("one or more Currin variants was scored by a non-held-out chromosome model")
    peaks = peak_index(peak_file)
    data["model_peak_overlap"] = [in_peak(peaks, c, p) for c, p in zip(data["chr_x"], data["pos_hg38"])]
    data.to_csv(out / "currin_scored_primary.tsv.gz", sep="\t", index=False, compression="gzip")

    y = data["label"].to_numpy(dtype=int)
    metric_arrays = {
        "ips": data["abs_logfc_x_jsd_x_active_allele_quantile"].to_numpy(float),
        "ies": data["abs_logfc_x_jsd"].to_numpy(float),
        "abs_logfc": data["abs_logfc"].to_numpy(float),
        "jsd": data["jsd"].to_numpy(float),
        "peak_overlap": data["model_peak_overlap"].to_numpy(float),
        "signed_logfc": data["logfc"].to_numpy(float),
    }
    points = {}
    for name in ("ips", "ies", "abs_logfc", "jsd", "peak_overlap"):
        auc, aprc = auc_ap(y, metric_arrays[name])
        points[name] = {"auroc": auc, "auprc": aprc}
    beta = pd.to_numeric(data["beta_alt"], errors="coerce").to_numpy(float)
    oriented = (y == 1) & np.isfinite(beta)
    points["direction"] = {
        "n_oriented_positive": int(oriented.sum()),
        "spearman": float(spearmanr(metric_arrays["signed_logfc"][oriented], beta[oriented]).statistic),
        "sign_concordance": float(np.mean(np.sign(metric_arrays["signed_logfc"][oriented]) == np.sign(beta[oriented]))),
    }
    points["ips_minus_peak_overlap_auroc"] = points["ips"]["auroc"] - points["peak_overlap"]["auroc"]

    global LABEL, SCORE, BETA, BLOCK_INDEX
    LABEL, SCORE, BETA = y, metric_arrays, beta
    BLOCK_INDEX = {str(block): idx.to_numpy() for block, idx in data.groupby("block_1mb", sort=True).groups.items()}
    seeds = np.random.SeedSequence(42).generate_state(args.bootstraps).tolist()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        boot = list(pool.map(one_bootstrap, seeds, chunksize=max(1, args.bootstraps // (args.workers * 4))))
    intervals = {}
    for key in boot[0]:
        intervals[key] = ci([row[key] for row in boot])

    # Freeze the published ChromBPNet IPS as primary discrimination metric.
    # Direction is evaluated independently using signed logFC (ALT vs REF).
    discrimination_pass = (
        intervals["ips_auroc"]["lower_95"] > 0.50
        and intervals["ips_minus_peak_overlap_auroc"]["lower_95"] > 0.0
    )
    direction_pass = (
        intervals["direction_spearman"]["lower_95"] > 0.0
        and intervals["direction_sign_concordance"]["lower_95"] > 0.50
    )
    verdict = {
        "currin_external_calibration_pass": bool(discrimination_pass and direction_pass),
        "discrimination_pass": bool(discrimination_pass),
        "direction_pass": bool(direction_pass),
        "primary_metric": "absolute IPS = abs(log2 ALT/REF count effect) * profile JSD * active-allele accessibility quantile",
        "simple_baseline": "binary overlap with an IDR-filtered adult-hepatocyte training peak",
        "gate_rule": "IPS AUROC lower95>0.5 AND paired block-bootstrap lower95(IPS AUROC - peak-overlap AUROC)>0 AND directional Spearman lower95>0 AND sign-concordance lower95>0.5",
        "n_variants": len(data),
        "n_positive": int(y.sum()),
        "n_background": int((y == 0).sum()),
        "n_blocks": len(BLOCK_INDEX),
        "n_bootstraps": args.bootstraps,
        "bootstrap_unit": "source 1-Mb genomic block",
        "point_estimates": points,
        "bootstrap_intervals": intervals,
        "effect_orientation": "ChromBPNet logFC and Currin beta_alt are both ALT versus REF",
        "scores_dir": str(scores_dir),
        "peak_file": str(peak_file),
    }
    with (out / "currin_external_calibration.json").open("w") as handle:
        json.dump(verdict, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(verdict, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
