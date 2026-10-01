#!/usr/bin/env python3
"""Input-only check of actual variance versus SD thresholds in donor RNA masks."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915"
FIXTURE = BASE / "model/donor_fixture_21773589"
SEED = 20260915


def stable(value):
    return int(hashlib.sha256(value.encode()).hexdigest()[:15], 16)


def thresholds(values, train):
    variance = values[train].var(axis=0)
    admitted = values[train].std(axis=0) > 1e-10
    removable = variance > 1e-10
    return variance, admitted, removable


def main(args):
    assert os.environ.get("SLURM_JOB_ID"), "Compute allocation required"
    args.out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, args.out / "executed_review_donor_mask_precheck.py")
    synthetic = np.array([[0., 0., 0.], [0., 1e-6, 1.], [0., 2e-6, 2.], [0., 3e-6, 3.]])
    variance, admitted, removable = thresholds(synthetic, np.arange(4))
    np.testing.assert_array_equal(admitted & ~removable, [False, True, False])
    # Deliberately do not access target_log2cpm_full_library or any fitted predictions.
    with np.load(FIXTURE / "paired.npz", allow_pickle=False) as stored:
        x = stored["rna_log2cpm_full_library"]
        fold = stored["donor_fold"]
        genes = stored["gene_ids"]
        chrom = stored["gene_chrom"]
    regions = pd.read_csv(FIXTURE / "regions.tsv", sep="\t")
    assert x.shape == (99, 42163) and len(set(genes)) == 42163 and np.isfinite(x).all()
    known = chrom != "unknown"
    assert known.sum() == 42126 and (~known).sum() == 37
    report, affected = [], []
    for held in range(5):
        for stage in ("inner", "outer"):
            train = np.flatnonzero((fold != held) & ((fold != (held + 1) % 5) if stage == "inner" else True))
            variance, admitted, removable = thresholds(x, train)
            for mode in ("trained_regions", "held_regions"):
                evaluated = regions if mode == "trained_regions" else regions.loc[regions.region_role.eq("held")]
                for target_chrom in sorted(evaluated.chrom.unique()):
                    explicit = (chrom == target_chrom) & removable
                    keep = known & ~explicit & admitted
                    leaked = keep & (chrom == target_chrom)
                    random_remove = np.zeros(len(genes), bool)
                    rng = np.random.default_rng(SEED + stable("random_chromosome|" + target_chrom))
                    random_remove[rng.choice(np.flatnonzero(known & removable), int(explicit.sum()), replace=False)] = True
                    random_keep = known & ~random_remove & admitted
                    assert keep.sum() == random_keep.sum(), "Effective random exclusion dimensions differ"
                    report.append({"mode": mode, "held_fold": held, "stage": stage, "target_chromosome": target_chrom,
                        "training_donors": len(train), "unrestricted_admitted_features": int(admitted.sum()),
                        "annotated_admitted_features": int((known & admitted).sum()),
                        "explicit_chromosome_removals": int(explicit.sum()), "effective_chromosome_admitted_features": int(keep.sum()),
                        "effective_random_admitted_features": int(random_keep.sum()),
                        "surviving_target_chromosome_features": int(leaked.sum()),
                        "minimum_surviving_training_variance": float(variance[leaked].min()) if leaked.any() else np.nan,
                        "maximum_surviving_training_variance": float(variance[leaked].max()) if leaked.any() else np.nan})
                    for index in np.flatnonzero(leaked):
                        affected.append({"mode": mode, "held_fold": held, "stage": stage, "target_chromosome": target_chrom,
                                         "gene_id": str(genes[index]), "training_variance": float(variance[index])})
    frame = pd.DataFrame(report)
    frame.to_csv(args.out / "actual_training_mask_thresholds.tsv", sep="\t", index=False)
    pd.DataFrame(affected, columns=["mode", "held_fold", "stage", "target_chromosome", "gene_id", "training_variance"]).to_csv(
        args.out / "surviving_target_chromosome_features.tsv", sep="\t", index=False)
    receipt = {"status": "actual_mask_threshold_mismatch_present" if affected else "no_actual_mask_threshold_mismatch",
        "job_id": os.environ["SLURM_JOB_ID"], "synthetic_threshold_invariant_passed": True,
        "partitions_checked": len(frame), "affected_partitions": int(frame.surviving_target_chromosome_features.gt(0).sum()),
        "distinct_affected_gene_ids": len({row["gene_id"] for row in affected}),
        "maximum_surviving_features_per_partition": int(frame.surviving_target_chromosome_features.max()),
        "all_effective_random_dimensions_match": True, "unrestricted_features": 42163, "annotated_universe": 42126,
        "thresholds": {"explicit_removal_training_variance_greater_than": 1e-10, "projection_admission_training_SD_greater_than": 1e-10},
        "interpretation": "Input-only exclusion validity check; nonzero survivors limit chromosome-exclusion claims, without invalidating unrestricted fits or implying a measured performance effect.",
        "molecular_effect_outcomes_or_model_predictions_read": False, "participant_values_exported": False,
        "source_files_edited": False, "new_training_or_selection": False}
    (args.out / "checks.json").write_text(json.dumps(receipt, indent=2) + "\n")
    assert sum(path.stat().st_size for path in args.out.iterdir()) < 10000000


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
