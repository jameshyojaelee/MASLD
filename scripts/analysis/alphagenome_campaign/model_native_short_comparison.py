#!/usr/bin/env python3
"""Compare 2-kb native and scalar heads without imposing 1-Mb coverage.

Declared while native-length inference was running, before its performance
tables existed. Complements, rather than replaces, the common-length population.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from model_native_lengths import calibrate_scores, check_scores
from model_final import describe_both


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Compute-node execution required")
    if not (args.lengths / "complete.json").is_file():
        raise ValueError("Full context-length scoring/comparison must finish first")
    args.out.mkdir(parents=True, exist_ok=False)
    labels = pd.read_csv(args.native / "currin_variants.tsv", sep="\t")
    labels = labels.loc[labels.heldout_fold.isin([1, 2, 3, 4])].copy()
    scores = check_scores(args.lengths / "native_2048.tsv", labels, 2048)
    names = {"local_atac_liver": "native_atac_2048", "local_dnase_liver": "native_dnase_2048"}
    table = labels.merge(scores[["key", *names]].rename(columns=names), on="key", validate="one_to_one")
    train = table.loc[table.heldout_fold.isin([2, 3, 4])].copy()
    valid = table.loc[table.heldout_fold.eq(1)].copy()
    calibrated, slopes = calibrate_scores(train, valid, list(names.values()))
    adapt = pd.read_csv(args.adaptation, sep="\t", index_col=0).rename_axis("variant_id").reset_index()
    adapt["key"] = adapt.variant_id.str.removeprefix("chr")
    arms = [c for c in adapt if c.startswith(("adapter_", "partial_", "frozen_"))]
    if len(arms) != 12 or adapt.key.duplicated().any():
        raise ValueError("Expected complete twelve-arm adaptation population")
    expected = labels.set_index("key").loc[adapt.key]
    if not expected.heldout_fold.eq(1).all():
        raise ValueError("Adaptation rows cross the permitted validation fold")
    np.testing.assert_allclose(adapt.observed_beta, expected.beta_alt, rtol=1e-7, atol=1e-7)
    matched = calibrated.merge(adapt[["key", "observed_beta", *arms]], on="key", validate="one_to_one")
    np.testing.assert_allclose(matched.observed_beta, matched.beta_alt, rtol=1e-7, atol=1e-7)
    matched = matched.drop(columns="observed_beta")
    omitted = adapt.loc[~adapt.key.isin(matched.key), ["key"]].copy()
    omitted["reason"] = "native_2kb_window_not_eligible"
    omitted.to_csv(args.out / "excluded_adapter_validation_rows.tsv", sep="\t", index=False)
    train[["key", "peak_id", "heldout_fold", "beta_alt", *names.values()]].to_csv(
        args.out / "calibration_training.tsv.gz", sep="\t", index=False)
    pd.DataFrame(slopes).to_csv(args.out / "training_calibration.tsv", sep="\t", index=False)
    native = [x + "_calibrated" for x in names.values()]
    describe_both(matched, [*native, *arms], ["native_atac_2048_calibrated"],
                  "same_length_2kb_eligible_population", args.out, 1000, len(adapt))
    # Reconcile metadata concerns identified independently during inference.
    pilots = [json.loads((args.lengths / f"pilot_{length}_run.json").read_text()) for length in (2048, 16384)]
    original = json.loads((args.lengths / "projection.json").read_text())
    conservative = sum(m["total_seconds"] / m["variants_attempted"] *
                       (len(labels)-m["variants_scored"]) for m in pilots) * 1.25 + 1500
    review = dict(actual_bootstrap_seed=20260915, initial_design_seed=20260916,
        seed_disposition="initial design seed was misrecorded; unchanged producer uses20260915; no seed chosen after results",
        initial_projection_scope="prediction-call time; excluded some host processing",
        conservative_host_inclusive_projection_seconds=conservative,
        remaining_allocation_at_original_admission_seconds=original["remaining_seconds"],
        conservative_projection_fits_original_allocation=conservative < original["remaining_seconds"]-300,
        calculation="pilot total loop time including first compilation per attempted row; 25percent buffer plus25min startup/analysis",
        full_inference_completed=(args.lengths / "complete.json").is_file())
    (args.out / "length_run_metadata_review.json").write_text(json.dumps(review, indent=2)+"\n")
    complete = dict(status="same_length_development_comparison_complete",
        training_n=len(train), validation_n=len(matched), original_adapter_validation_n=len(adapt),
        calibration_training_folds=[2, 3, 4], evaluation_fold=1, fold0_evaluated=False,
        source_selection="historically selected Currin caQTL leads", input_length_bp=2048,
        native_readout="501bp ATAC/DNase track sums; distinct from shared scalar head",
        one_megabase_coverage_required=False, adapter_fits_unchanged=True,
        original_common_length_comparison_preserved=True,
        independent_numerical_reconstruction="pending", source_beta_units="FastQTL ALT-dosage coefficient",
        no_external_generalization=True, uncertainty="fixed source effects/predictions; no fitting/selection uncertainty",
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.out / "complete.json").write_text(json.dumps(complete, indent=2)+"\n")
    print(json.dumps(complete, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--lengths", type=Path, required=True)
    parser.add_argument("--adaptation", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
