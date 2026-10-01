#!/usr/bin/env python3
"""Full existing Currin lead population: native scoring input and matched development results.

This population is significance-selected and historically deduplicated by q-value.
It is NOT an outcome-independent all-variant/peak population. No new selection
on molecular outcomes occurs here. Native scores retain log2 track-sum units;
their RMSE against caQTL beta is inapplicable until training-only calibration.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJ / "scripts/analysis/alphagenome_program"))
import i1_common as C


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(out):
    out.mkdir(parents=True, exist_ok=True)
    target = out / "currin_variants.tsv"
    if target.exists():
        raise FileExistsError(target)
    labels = pd.read_csv(C.C2_LABELS, sep="\t")
    labels["key"] = [C.key_of(*r) for r in labels[["chr", "pos_hg38", "ref", "alt"]].itertuples(index=False, name=None)]
    if labels.key.duplicated().any() or len(labels) != 32322:
        raise ValueError("Existing Currin comparison population changed")
    if labels.groupby("chr").heldout_fold.nunique().max() != 1:
        raise ValueError("Chromosome split permits overlapping windows")
    labels.to_csv(target, sep="\t", index=False)
    prior_path = PROJ / "GWAS/finemapping/results/alphagenome_program/ad1-gate-20260915T114615Z/raw/gate_1048576.tsv"
    prior = pd.read_csv(prior_path, sep="\t")
    meta = json.loads(prior_path.with_name("gate_1048576_run.json").read_text())
    if meta["checkpoint"] != str(C.CHECKPOINT) or meta["fasta"] != C.FASTA_PATH:
        raise ValueError("Archived native input definitions changed")
    if not ((prior.length_bp == 1048576) & (prior.card_tag == "l40s") &
            (prior.mask_width_bp == 501) & (prior.aggregation == "DIFF_LOG2_SUM")).all():
        raise ValueError("Archived scores differ from this native recipe")
    prior = prior.loc[prior.key.isin(labels.key)].copy()
    if prior.key.duplicated().any():
        raise ValueError("Duplicate archived scores")
    check = labels.set_index("key").loc[prior.key]
    if not np.allclose(check.beta_alt, prior.beta_alt, rtol=0, atol=1e-12):
        raise ValueError("Reused score labels do not match")
    prior.to_csv(out / "native_1048576.tsv", sep="\t", index=False)
    # The extraction table deliberately contains no outcome columns.
    cols = ["key", "chr", "pos_hg38", "ref", "alt", "peak_id", "peak_start_hg38", "peak_stop_hg38", "heldout_fold", "block_1mb"]
    labels[cols].to_csv(out / "sequence_manifest.tsv", sep="\t", index=False)
    receipt = {"eligible_existing_leads": len(labels), "reused_native_rows": len(prior),
               "remaining_scoring_attempts": len(labels)-len(prior),
               "population": "existing_significance_selected_deduplicated_leads",
               "broader_outcome_independent_population": False,
               "folds": "existing_whole_chromosome_holdouts",
               "label_sha256": digest(C.C2_LABELS), "reused_score_sha256": digest(prior_path),
               "native_score_units": "ALT_minus_REF_log2_1_plus_501bp_predicted_track_sum",
               "measured_units": "source_FastQTL_ALT_dosage_beta",
               "protected_outcomes_read": False}
    (out / "population.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt), flush=True)


def evaluate(out):
    labels = pd.read_csv(out / "currin_variants.tsv", sep="\t")
    scores = pd.read_csv(out / "native_1048576.tsv", sep="\t")
    if scores.key.duplicated().any():
        raise ValueError("Duplicate native predictions")
    reject_path = out / "native_1048576_rejected.tsv"
    rejected = pd.read_csv(reject_path, sep="\t") if reject_path.exists() else pd.DataFrame(columns=["key"])
    unresolved = set(labels.key) - set(scores.key) - set(rejected.key)
    if unresolved:
        raise ValueError(f"{len(unresolved)} variants remain unattempted; full-population comparison not ready")
    if set(scores.key) & set(rejected.key):
        raise ValueError("A variant is both scored and rejected")
    frame = labels.merge(scores[["key", "local_atac_liver", "local_dnase_liver"]], on="key", validate="one_to_one")
    old = pd.read_csv(C.C2_OOF, sep="\t")
    old["key"] = old.variant_id.str.replace("chr", "", regex=False)
    old_arms = [x for x in old if x.endswith("__ensemble")]
    frame = frame.merge(old[["key", *old_arms]], on="key", validate="one_to_one")
    if len(frame) != len(scores):
        raise ValueError("Comparator coverage differs")
    y = frame.beta_alt.to_numpy()
    folds = frame.heldout_fold.to_numpy()
    calibration = []
    for native in ("local_atac_liver", "local_dnase_liver"):
        x = frame[native].to_numpy()
        predicted = np.full(len(frame), np.nan)
        for fold in range(5):
            train, test = folds != fold, folds == fold
            # No intercept: zero and swap symmetry are preserved by calibration.
            slope = float(np.dot(x[train], y[train]) / max(np.dot(x[train], x[train]), 1e-20))
            predicted[test] = slope * x[test]
            calibration.append({"arm": native, "heldout_fold": fold, "training_slope": slope,
                                "training_n": int(train.sum()), "intercept": 0.0})
        frame[native + "__train_calibrated"] = predicted
    all_arms = ["local_atac_liver", "local_dnase_liver", "local_atac_liver__train_calibrated", "local_dnase_liver__train_calibrated", *old_arms]
    rows = []
    for arm in all_arms:
        pred = frame[arm].to_numpy()
        native = arm in ("local_atac_liver", "local_dnase_liver")
        rows.append({"arm": arm, "n": len(frame), "eligible_n": len(labels),
                     "coverage": len(frame)/len(labels),
                     "blocks": frame.block_1mb.nunique(), "macro_spearman": C.macro_over_folds(y, pred, folds),
                     "pooled_spearman": C.fast_spearman(y, pred),
                     "caQTL_beta_RMSE": np.nan if native else float(np.sqrt(np.mean((pred-y)**2))),
                     "caQTL_beta_MAE": np.nan if native else float(np.mean(np.abs(pred-y))),
                     "error_status": "inapplicable_unmatched_units" if native else "native_label_units",
                     "evidence_status": "development_chromosome_heldout"})
    pd.DataFrame(rows).to_csv(out / "matched_performance.tsv", sep="\t", index=False)
    pd.DataFrame(calibration).to_csv(out / "training_calibration.tsv", sep="\t", index=False)
    frame.to_csv(out / "matched_predictions.tsv.gz", sep="\t", index=False)
    (out / "completion.json").write_text(json.dumps({"eligible": len(labels), "scored": len(scores),
        "rejected": rejected.key.nunique(), "unresolved": 0, "protected_evaluation": False,
        "uncertainty_status": "paired_block_bootstrap_not_yet_run", "specialists": "separate_eligible_subset_pending"}, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["prepare", "evaluate"])
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    (prepare if args.mode == "prepare" else evaluate)(args.out)
