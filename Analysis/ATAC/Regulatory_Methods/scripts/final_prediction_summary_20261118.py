#!/usr/bin/env python3
"""Summarize 20261118 hepatocyte refits for BOTH total OOF MSE and
donor-mean-residualized within-donor MSE across baseline, pairing_shuffle,
and shuffle_cross_donor scenarios.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


A = Path(__file__).resolve().parents[1]
OUT_ROOT = A / "results/20261118_prediction_refits"
METHODS = ["promoter_only", "distance", "cicero"]
KEYS = ["outer_fold", "gene_id", "donor_id"]
SCENARIOS = [
    "baseline",
    "pairing_shuffle_43",
    "pairing_shuffle_44",
    "pairing_shuffle_45",
    "shuffle_cross_donor_43",
    "shuffle_cross_donor_44",
    "shuffle_cross_donor_45",
]


def read_predictions(path: Path) -> pd.DataFrame:
    parts = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=100000, dtype={"donor_id": str}):
        chunk = chunk[chunk.method.isin(METHODS)]
        parts.append(chunk)
    df = pd.concat(parts, ignore_index=True)
    assert not df.duplicated(KEYS + ["method"]).any()
    return df


def matched_metric(frame: pd.DataFrame, metric_col: str) -> pd.DataFrame:
    wide = frame.pivot(index=KEYS, columns="method", values=metric_col)
    return wide.dropna(subset=METHODS)


def estimate(donor_df: pd.DataFrame, scenario: str, baseline: str, metric_name: str, keys_count: int) -> dict:
    values = (donor_df[baseline] - donor_df["cicero"]).to_numpy()
    salt = int(hashlib.sha256((scenario + baseline + metric_name).encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng((42 + salt) % 2**32)
    means = np.array([rng.choice(values, len(values), replace=True).mean() for _ in range(2000)])
    lo, hi = np.quantile(means, [0.025, 0.975])
    return dict(
        metric=metric_name,
        scenario=scenario,
        baseline=baseline,
        n_donors=len(values),
        n_gene_donor_fold_pairs=keys_count,
        mean_baseline=float(donor_df[baseline].mean()),
        mean_cicero=float(donor_df["cicero"].mean()),
        mean_gain=float(values.mean()),
        ci_low=float(lo),
        ci_high=float(hi),
        relative_error_reduction_percent=float(100.0 * values.mean() / donor_df[baseline].mean()),
        n_positive_donors=int((values > 0).sum()),
    )


def summarize_metric(frames: dict[str, pd.DataFrame], metric_col: str, metric_label: str):
    bwide = matched_metric(frames["baseline"], metric_col)
    bdonor = bwide.groupby(level="donor_id")[METHODS].mean()
    rows = [estimate(bdonor, "baseline", b, metric_label, len(bwide)) for b in ["promoter_only", "distance"]]
    donor_rows = [bdonor.reset_index().assign(metric=metric_label, scenario="baseline")]
    changes = []
    for name in SCENARIOS:
        if name == "baseline":
            continue
        wide = matched_metric(frames[name], metric_col)
        common = wide.index.intersection(bwide.index)
        a = bwide.loc[common].groupby(level="donor_id")[METHODS].mean()
        b = wide.loc[common].groupby(level="donor_id")[METHODS].mean()
        rows.extend(estimate(b, name, comp, metric_label, len(common)) for comp in ["promoter_only", "distance"])
        donor_rows.append(b.reset_index().assign(metric=metric_label, scenario=name))
        for comp in ["promoter_only", "distance"]:
            diffs = (a[comp] - a["cicero"]) - (b[comp] - b["cicero"])
            for d_id, v in diffs.items():
                changes.append(dict(metric=metric_label, scenario=name, baseline=comp, donor_id=d_id, baseline_gain_minus_scenario_gain=float(v)))

    change_df = pd.DataFrame(changes)
    csum = []
    for (sc, comp), grp in change_df.groupby(["scenario", "baseline"], sort=True):
        vals = grp.baseline_gain_minus_scenario_gain.to_numpy()
        salt = int(hashlib.sha256((sc + comp + metric_label).encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng((42 + salt) % 2**32)
        samp = np.array([rng.choice(vals, len(vals), replace=True).mean() for _ in range(2000)])
        lo, hi = np.quantile(samp, [0.025, 0.975])
        csum.append(
            dict(
                metric=metric_label,
                scenario=sc,
                baseline=comp,
                n_donors=len(vals),
                baseline_gain_minus_scenario_gain=float(vals.mean()),
                ci_low=float(lo),
                ci_high=float(hi),
            )
        )
    return pd.DataFrame(rows), pd.concat(donor_rows, ignore_index=True), pd.DataFrame(csum)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=OUT_ROOT)
    args = p.parse_args()
    assert os.environ.get("SLURM_JOB_ID"), "Must run on compute node via SLURM"

    summary_dir = args.root / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    frames = {}
    for name in SCENARIOS:
        parts = []
        for fold in range(5):
            fdir = args.root / name / f"fold{fold}"
            comp = json.loads((fdir / "completion.json").read_text())
            assert comp["status"] == "completed" and not comp.get("technical_pilot", False)
            parts.append(read_predictions(fdir / "predictions.tsv"))
        frames[name] = pd.concat(parts, ignore_index=True)

    res_tot, don_tot, ch_tot = summarize_metric(frames, "standardized_squared_error", "total_oof_standardized_mse")
    res_wd, don_wd, ch_wd = summarize_metric(frames, "within_donor_standardized_squared_error", "within_donor_standardized_mse")
    res_wraw, _, ch_wraw = summarize_metric(frames, "within_donor_squared_error", "within_donor_raw_mse")

    res_tot.to_csv(summary_dir / "comparisons_total_oof.tsv", sep="\t", index=False)
    res_wd.to_csv(summary_dir / "comparisons_within_donor.tsv", sep="\t", index=False)
    res_wraw.to_csv(summary_dir / "comparisons_within_donor_raw.tsv", sep="\t", index=False)
    pd.concat([don_tot, don_wd], ignore_index=True).to_csv(summary_dir / "donor_scores.tsv", sep="\t", index=False)
    pd.concat([ch_tot, ch_wd, ch_wraw], ignore_index=True).to_csv(summary_dir / "paired_scenario_change_summary.tsv", sep="\t", index=False)

    lines = ["Hepatocyte Cicero 20261118 (Total OOF vs Within-Donor Residualized MSE + Cross-Donor Shuffle)", ""]
    for df_r, label in [(res_tot, "TOTAL OOF"), (res_wd, "WITHIN-DONOR RESIDUALIZED")]:
        lines.append(f"--- {label} ---")
        for r in df_r.to_dict("records"):
            lines.append(
                f"{r['scenario']} vs {r['baseline']}: gain={r['mean_gain']:.6g} [{r['ci_low']:.6g}, {r['ci_high']:.6g}], "
                f"rer={r['relative_error_reduction_percent']:.3f}%, pos={r['n_positive_donors']}/{r['n_donors']}"
            )
        lines.append("")
    (summary_dir / "REPORT.txt").write_text("\n".join(lines) + "\n")
    (summary_dir / "completion.json").write_text(json.dumps({"status": "completed", "scenarios": SCENARIOS}, indent=2) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
