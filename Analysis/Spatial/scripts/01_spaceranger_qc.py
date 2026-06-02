#!/usr/bin/env python3
"""
01_spaceranger_qc.py — Automated QC assessment of SpaceRanger outputs.

Parses metrics_summary.csv from each sample and applies go/no-go thresholds.
Generates a summary CSV and flags samples that fail QC.
"""

import pathlib
import sys
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    SPATIAL_ROOT, RESULTS_DIR, load_config, load_dataset_config,
    parse_spaceranger_metrics, assess_spaceranger_qc,
    print_header, print_step, save_csv,
)


def main():
    print_header("01: SpaceRanger QC Assessment")

    config = load_config()
    datasets = load_dataset_config()
    thresholds = config["qc"]["spaceranger"]

    all_metrics = []

    for dataset_name, ds_config in datasets["datasets"].items():
        sr_base = RESULTS_DIR / "spaceranger" / dataset_name
        if not sr_base.exists():
            print(f"  Skipping {dataset_name}: no SpaceRanger output directory")
            continue

        print(f"\n  Dataset: {dataset_name}")
        sample_dirs = sorted([d for d in sr_base.iterdir() if d.is_dir()])

        for sample_dir in sample_dirs:
            outs_dir = sample_dir / "outs"
            sample_id = sample_dir.name
            print_step(f"Checking {sample_id}...")

            metrics = parse_spaceranger_metrics(outs_dir)
            if metrics.get("status") == "MISSING":
                print(f"    FAIL: metrics_summary.csv not found")
                metrics["sample_id"] = sample_id
                metrics["dataset"] = dataset_name
                metrics["pass_qc"] = False
                metrics["failures"] = "metrics_summary.csv missing"
                all_metrics.append(metrics)
                continue

            pass_qc, failures = assess_spaceranger_qc(metrics, thresholds)
            metrics["sample_id"] = sample_id
            metrics["dataset"] = dataset_name
            metrics["pass_qc"] = pass_qc
            metrics["failures"] = "; ".join(failures) if failures else ""

            status = "PASS" if pass_qc else "FAIL"
            print(f"    {status}: {sample_id}")
            if failures:
                for f in failures:
                    print(f"      - {f}")

            all_metrics.append(metrics)

    # Summary
    df = pd.DataFrame(all_metrics)
    save_csv(df, "spaceranger_qc_summary.csv", subdir="qc")

    n_pass = df["pass_qc"].sum()
    n_total = len(df)
    print(f"\n  QC Summary: {n_pass}/{n_total} samples pass ({n_pass/max(n_total,1)*100:.0f}%)")

    # Per-dataset summary
    for dataset in df["dataset"].unique():
        ds_df = df[df["dataset"] == dataset]
        ds_pass = ds_df["pass_qc"].sum()
        print(f"    {dataset}: {ds_pass}/{len(ds_df)} pass")

    print_header("01: QC Complete")


if __name__ == "__main__":
    main()
