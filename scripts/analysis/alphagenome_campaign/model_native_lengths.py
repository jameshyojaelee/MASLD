#!/usr/bin/env python3
"""Match native readouts across context lengths before interpreting adaptation.

The 501-bp track readout, tracks, checkpoint and GPU remain fixed. All native
calibration uses the same eligible training rows in folds 2--4. Fold 1 is an
already inspected development population; fold 0 is not evaluated here.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts/analysis/alphagenome_program"))
import i1_common as C
from model_final import describe_both

LENGTHS = (2048, 16384)
SEED = 20260916


def emit(path, value):
    path.write_text(json.dumps(value, indent=2, default=str) + "\n")


def calibrate_scores(train, valid, columns):
    if set(train.heldout_fold) != {2, 3, 4} or set(valid.heldout_fold) != {1}:
        raise ValueError("Native calibration must use folds 2--4, evaluation fold 1")
    if set(train.key) & set(valid.key):
        raise ValueError("Training/evaluation identities overlap")
    result = valid.copy()
    receipts = []
    for column in columns:
        x = train[column].to_numpy(dtype=float)
        y = train.beta_alt.to_numpy(dtype=float)
        if not np.isfinite(x).all() or not np.isfinite(y).all() or np.dot(x, x) <= 0:
            raise ValueError("Unidentifiable native calibration")
        slope = float(x @ y / (x @ x))
        result[column + "_calibrated"] = slope * valid[column]
        receipts.append(dict(arm=column, slope=slope, intercept=0,
                             training_n=len(train), training_folds="2,3,4"))
    return result, receipts


def check_scores(path, labels, length):
    scores = pd.read_csv(path, sep="\t")
    if scores.key.duplicated().any() or not set(scores.key).issubset(set(labels.key)):
        raise ValueError("Duplicate or unexpected native identity")
    expected = labels.set_index("key").loc[scores.key]
    for column in ("chr", "ref", "alt"):
        if not np.array_equal(scores[column].to_numpy(), expected[column].to_numpy()):
            raise ValueError(f"Mismatched native {column}")
    for column in ("pos_hg38", "heldout_fold", "beta_alt"):
        np.testing.assert_allclose(scores[column], expected[column], rtol=0, atol=1e-12)
    if not ((scores.length_bp == length) & (scores.card_tag == "l40s") &
            (scores.mask_width_bp == 501) & (scores.aggregation == "DIFF_LOG2_SUM")).all():
        raise ValueError("Native hardware/readout differs")
    if not np.isfinite(scores[["local_atac_liver", "local_dnase_liver"]]).all().all():
        raise ValueError("Nonfinite native scores")
    return scores


def run_scorer(manifest, length, output, scoring_python, resume=False):
    command = [str(scoring_python), str(ROOT / "scripts/analysis/alphagenome_program/i1_01_score_local.py"),
               "--variants", str(manifest), "--length", str(length), "--out", str(output),
               "--require-card", "l40s", "--flush-every", "100"]
    if resume:
        command.append("--resume")
    subprocess.run(command, check=True)
    metadata = json.loads(output.with_name(output.stem + "_run.json").read_text())
    if (metadata["checkpoint"] != str(C.CHECKPOINT) or metadata["fasta"] != C.FASTA_PATH or
            metadata["atac_liver_track_names"] != C.AG_ATAC_TRACKS or
            metadata["dnase_liver_track_names"] != C.AG_DNASE_TRACKS):
        raise ValueError("Scoring source or track identity changed")
    return metadata


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Scientific workload requires a compute node")
    started = time.monotonic()
    args.out.mkdir(parents=True, exist_ok=False)
    emit(args.out / "design.json", dict(
        question="Does context length explain the native versus 2kb-adapter difference?",
        lengths=[2048, 16384, 1048576], readout="same 501bp DIFF_LOG2_SUM liver ATAC/DNase tracks",
        training_folds=[2, 3, 4], development_fold=1, fold0_evaluated=False,
        inference_only=True, protected_outcomes_read=False,
        population="historically selected Currin leads; common native coverage primary",
        family="all reported native-length and deposited-adaptation contrasts against native ATAC at 2kb and 1Mb",
        correction="BH separately within complete historical-bin and chromosome exploratory families",
        uncertainty="conditional on fitted models and source effects; no retraining or selection uncertainty",
        stopping_rule="100 outcome-independent scoring rows per length; projection plus25percent and20min analysis must fit remaining allocation",
        no_new_hyperparameter_selection=True, seed=SEED))
    labels = pd.read_csv(args.native / "currin_variants.tsv", sep="\t")
    labels = labels.loc[labels.heldout_fold.isin([1, 2, 3, 4])].copy()
    if labels.key.duplicated().any() or labels.groupby("chr").heldout_fold.nunique().max() != 1:
        raise ValueError("Canonical identities or chromosome folds disagree")
    manifest = args.out / "permitted_variants.tsv"
    labels.to_csv(manifest, sep="\t", index=False)
    # Sequence position, not effect magnitude, chooses representative pilot rows.
    pilot = labels.iloc[np.linspace(0, len(labels)-1, 100, dtype=int)].copy()
    pilot_manifest = args.out / "pilot_variants.tsv"
    pilot.to_csv(pilot_manifest, sep="\t", index=False)
    pilot_meta = {}
    for length in LENGTHS:
        destination = args.out / f"pilot_{length}.tsv"
        pilot_meta[length] = run_scorer(pilot_manifest, length, destination, args.scoring_python)
        scored = check_scores(destination, pilot, length)
        if len(scored) < 90:
            raise ValueError("Pilot has too few scored rows for a runtime projection")
    projected = sum(
        (len(labels)-m["variants_scored"])*m["warm_call_seconds_mean"] +
        m["checkpoint_load_seconds"] + m["first_variant_seconds"]
        for m in pilot_meta.values()) * 1.25 + 1200
    remaining = args.max_hours*3600 - (time.monotonic()-started)
    emit(args.out / "projection.json", dict(projected_remaining_seconds=projected,
        remaining_seconds=remaining, admitted=projected < remaining-300, pilot=pilot_meta))
    if projected >= remaining-300:
        emit(args.out / "DEFERRED.json", {"reason": "measured projection exceeds allocation"})
        return
    tables = {}
    coverage = []
    for length in LENGTHS:
        destination = args.out / f"native_{length}.tsv"
        shutil.copyfile(args.out / f"pilot_{length}.tsv", destination)
        # Pilot rejections are retried and recorded once by the full scorer.
        run_scorer(manifest, length, destination, args.scoring_python, resume=True)
        scores = check_scores(destination, labels, length)
        rejected_path = destination.with_name(destination.stem + "_rejected.tsv")
        rejected = pd.read_csv(rejected_path, sep="\t") if rejected_path.exists() else pd.DataFrame(columns=["key"])
        if rejected.key.duplicated().any() or set(scores.key) & set(rejected.key):
            raise ValueError("Scored/rejected identities are not disjoint")
        if set(scores.key) | set(rejected.key) != set(labels.key):
            raise ValueError("Incomplete or unexpected native attempt accounting")
        tables[length] = scores
        coverage.append(dict(length=length, eligible=len(labels), scored=len(scores), rejected=len(rejected)))
    old = pd.read_csv(args.native / "native_1048576.tsv", sep="\t")
    old = old.loc[old.key.isin(labels.key)]
    # The prior full-population scorer has the same pinned definitions.
    temporary = args.out / "existing_1048576_permitted.tsv"
    old.to_csv(temporary, sep="\t", index=False)
    tables[1048576] = check_scores(temporary, labels, 1048576)
    coverage.append(dict(length=1048576, eligible=len(labels), scored=len(old),
                         rejected=len(labels)-len(old)))
    pd.DataFrame(coverage).to_csv(args.out / "coverage.tsv", sep="\t", index=False)
    common = labels.copy()
    raw_columns = []
    for length, scores in tables.items():
        names = {f"local_{channel}_liver": f"native_{channel}_{length}" for channel in ("atac", "dnase")}
        common = common.merge(scores[["key", *names]].rename(columns=names), on="key", validate="one_to_one")
        raw_columns.extend(names.values())
    train = common.loc[common.heldout_fold.isin([2, 3, 4])].copy()
    valid = common.loc[common.heldout_fold == 1].copy()
    result, receipts = calibrate_scores(train, valid, raw_columns)
    train[["key", "peak_id", "heldout_fold", "beta_alt", *raw_columns]].to_csv(
        args.out / "calibration_training.tsv.gz", sep="\t", index=False)
    pd.DataFrame(receipts).to_csv(args.out / "training_calibration.tsv", sep="\t", index=False)
    adapt = pd.read_csv(args.adaptation, sep="\t", index_col=0).rename_axis("variant_id").reset_index()
    adapt["key"] = adapt.variant_id.str.removeprefix("chr")
    arms = [c for c in adapt if c.startswith(("adapter_", "partial_", "frozen_"))]
    if len(arms) != 12 or adapt.key.duplicated().any():
        raise ValueError("Expected all twelve deposited adaptation/control arms")
    result = result.merge(adapt[["key", "observed_beta", *arms]], on="key", validate="one_to_one")
    if len(result) != len(valid):
        raise ValueError("Adapter validation coverage differs from common native population")
    np.testing.assert_allclose(result.observed_beta, result.beta_alt, rtol=1e-7, atol=1e-7)
    result = result.drop(columns="observed_beta")
    native_arms = [c + "_calibrated" for c in raw_columns]
    describe_both(result, native_arms + arms,
                  ["native_atac_2048_calibrated", "native_atac_1048576_calibrated"],
                  "matched_native_lengths", args.out, 1000, len(adapt))
    emit(args.out / "complete.json", dict(status="development_comparison_complete",
        native_attempts=coverage, common_training_n=len(train), common_validation_n=len(result),
        original_adapter_validation_n=len(adapt), all_arms=18, original_adapter_training_unchanged=True,
        inference_seconds=time.monotonic()-started, independent_reconstruction="pending",
        main_limit="input length is matched for 2kb native/adapted comparison; native track readout and scalar head remain distinct",
        protected_outcomes_read=False, formal_five_seed_finalist=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--adaptation", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--scoring-python", type=Path, required=True)
    parser.add_argument("--max-hours", type=float, default=6)
    main(parser.parse_args())
