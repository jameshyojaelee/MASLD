#!/usr/bin/env python3
"""Independent native-length and unrestricted 2-kb reconstruction.

Only this reviewer's checked numerical helpers are reused; no producer
inferential function is imported. Fold-0 effect strings are discarded before
numeric conversion. Source identities and training calibration precede every
result comparison. All MSE families and three key ranking contrasts are checked.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd

import review_native_final as R

ATAC = ["UBERON:0001114 ATAC-seq", "UBERON:0001115 ATAC-seq", "UBERON:0002107 ATAC-seq"]
DNASE = ["UBERON:0001114 DNase-seq", "UBERON:0001115 DNase-seq"]
IDENTITY = ["key", "chr", "pos_hg38", "ref", "alt", "peak_id", "peak_start_hg38",
            "peak_stop_hg38", "heldout_fold", "block_1mb"]
BASE = R.ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model"
CHECKPOINT = R.ROOT / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-weights-20260915T103442Z/checkpoints"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"


def permitted_labels(path, variant_column):
    fields = [variant_column, *[c for c in IDENTITY if c != "key"], "beta_alt"]
    frame = pd.read_csv(path, sep="\t", usecols=fields, dtype={"beta_alt": str})
    total = len(frame)
    frame = frame.loc[frame.heldout_fold.isin([1, 2, 3, 4])].copy()
    frame["key"] = frame[variant_column].str.removeprefix("chr")
    frame["beta_alt"] = pd.to_numeric(frame.beta_alt)
    return R.indexed(frame), total


def metadata(audit, path, length):
    meta = json.loads(path.read_text())
    for field, expected in (("checkpoint", str(CHECKPOINT)), ("fasta", FASTA), ("length_bp", length),
                            ("atac_liver_track_names", ATAC), ("dnase_liver_track_names", DNASE)):
        audit.equal(path.name + ":" + field, np.asarray(meta[field]), np.asarray(expected))
    return meta


def score_table(audit, path, length, labels):
    frame = pd.read_csv(path, sep="\t", dtype={"beta_alt": str})
    if length == 1048576:
        frame = frame.loc[frame.heldout_fold.isin([1, 2, 3, 4])].copy()
    elif not frame.heldout_fold.isin([1, 2, 3, 4]).all():
        raise ValueError("New scoring included a forbidden fold")
    frame["beta_alt"] = pd.to_numeric(frame.beta_alt)
    frame = R.indexed(frame)
    if not set(frame.key).issubset(labels.index):
        raise ValueError("Unexpected native variant")
    R.check_source_rows(audit, f"raw{length}", frame, labels)
    audit.equal(f"raw{length}:source_chromosome", frame.chr, labels.loc[frame.key, "chr"])
    for field, expected in (("length_bp", length), ("mask_width_bp", 501),
                            ("aggregation", "DIFF_LOG2_SUM"), ("card_tag", "l40s"), ("n_non_acgt", 0)):
        audit.equal(f"raw{length}:{field}", frame[field], np.repeat(expected, len(frame)))
    audit.equal(f"raw{length}:window_length", frame.window_end-frame.window_start0, np.full(len(frame), length))
    audit.equal(f"raw{length}:mask_span", frame.mask_end0_offset-frame.mask_start0_offset, np.full(len(frame), 501))
    audit.equal(f"raw{length}:genomic501bp", frame.window_start0+frame.mask_start0_offset, frame.pos_hg38-1-250)
    for channel, tracks in (("atac", ATAC), ("dnase", DNASE)):
        columns = ["local_" + channel + "__" + name.split()[0] for name in tracks]
        values = frame[columns].to_numpy(float)
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite native component tracks")
        audit.close(f"raw{length}:{channel}:mean", frame["local_"+channel+"_liver"], values.mean(axis=1), atol=1e-12)
    return frame


def prepare(audit, args):
    labels, total = permitted_labels(R.LABELS, "lead_variant_id")
    archive, archive_total = permitted_labels(args.native / "currin_variants.tsv", "key")
    audit.equal("original_source_counts", [total, archive_total], [32322, 32322])
    audit.equal("permitted_identity_order", archive[IDENTITY].to_numpy(), labels[IDENTITY].to_numpy())
    audit.close("permitted_labels", archive.beta_alt, labels.beta_alt, atol=1e-12)
    full_meta = json.loads((args.lengths / "complete.json").read_text())
    short_meta = json.loads((args.short / "complete.json").read_text())
    if full_meta["status"] != "development_comparison_complete" or short_meta["status"] != "same_length_development_comparison_complete":
        raise ValueError("Native-length or short analysis incomplete; not evaluated")
    if short_meta["fold0_evaluated"] or short_meta["one_megabase_coverage_required"]:
        raise ValueError("Short comparison eligibility/exposure declaration differs")
    manifest = R.indexed(R.read(args.lengths / "permitted_variants.tsv"))
    audit.equal("scorer_identity_order", manifest[IDENTITY].to_numpy(), labels[IDENTITY].to_numpy())
    audit.close("scorer_input_labels", manifest.beta_alt, labels.beta_alt, atol=1e-12)
    tables, coverage = {}, []
    for length in (2048, 16384, 1048576):
        folder = args.native if length == 1048576 else args.lengths
        tables[length] = score_table(audit, folder / f"native_{length}.tsv", length, labels)
        metadata(audit, folder / f"native_{length}_run.json", length)
        rejected_path = folder / f"native_{length}_rejected.tsv"
        rejected = R.read(rejected_path) if rejected_path.exists() else pd.DataFrame(columns=["key"])
        if length == 1048576:
            rejected = rejected.loc[rejected.key.isin(labels.index)].copy()
        if rejected.key.duplicated().any() or set(rejected.key) & set(tables[length].key):
            raise ValueError("Scored/rejected identities conflict")
        audit.equal(f"attempts{length}", sorted(set(rejected.key) | set(tables[length].key)), sorted(labels.key))
        if len(rejected):
            R.check_source_rows(audit, f"rejected{length}", rejected, labels)
        coverage.append(dict(length=length, eligible=len(labels), scored=len(tables[length]), rejected=len(rejected),
            training_scored=int(tables[length].heldout_fold.isin([2, 3, 4]).sum()),
            validation_scored=int(tables[length].heldout_fold.eq(1).sum())))
    reuse = R.ROOT / "GWAS/finemapping/results/alphagenome_program/ad1-gate-20260915T114615Z/raw/gate_1048576_run.json"
    metadata(audit, reuse, 1048576)
    reported = R.read(args.lengths / "coverage.tsv").set_index("length")
    for row in coverage:
        for col in ("eligible", "scored", "rejected"):
            audit.equal(f"coverage:{row['length']}:{col}", [reported.loc[row["length"], col]], [row[col]])
    adapt = pd.read_csv(args.adaptation, sep="\t", index_col=0).rename_axis("variant_id").reset_index()
    adapt["key"] = adapt.variant_id.str.removeprefix("chr")
    adapt = R.indexed(adapt)
    audit.equal("original_adapter_validation", sorted(adapt.key), sorted(labels.loc[labels.heldout_fold.eq(1), "key"]))
    audit.close("original_adapter_labels", adapt.observed_beta, labels.loc[adapt.key, "beta_alt"].to_numpy(np.float32).astype(float), atol=1e-12)
    runs = json.loads((args.sweep / "runs.json").read_text())
    expected = [r["name"] for r in runs if r.get("status") == "completed_fixed_exposure" and r.get("steps_requested") == 5000]
    arms = [c for c in adapt if c.startswith(("adapter_", "partial_", "frozen_"))]
    audit.equal("original_model_family", sorted(arms), sorted(expected))
    audit.equal("original_model_count", [len(arms)], [12])
    for run in runs:
        if run["name"] not in arms:
            continue
        pred = R.read(Path(run["output"]) / "validation_predictions.tsv")
        pred["key"] = pred.variant_id.str.removeprefix("chr")
        pred = R.indexed(pred)
        audit.equal(run["name"] + ":original_population", sorted(pred.key), sorted(adapt.key))
        audit.close(run["name"] + ":original_predictions", adapt[run["name"]], pred.loc[adapt.key, "predicted_beta"], atol=1e-12)
    return labels, tables, adapt, arms, coverage, full_meta, short_meta


def calibration(audit, folder, name, lengths, labels, tables, adapt, adapter_arms):
    common = set(labels.key)
    for length in lengths:
        common &= set(tables[length].key)
    train = labels.loc[labels.key.isin(common) & labels.heldout_fold.isin([2, 3, 4])]
    valid = labels.loc[labels.key.isin(common) & labels.heldout_fold.eq(1)]
    frame = R.indexed(R.read(folder / (name + "_matched_predictions.tsv.gz")))
    audit.equal(name + ":exact_validation", sorted(frame.key), sorted(valid.key))
    R.check_source_rows(audit, name, frame, labels)
    for field in ("peak_start_hg38", "peak_stop_hg38"):
        audit.equal(name + ":" + field, frame[field], labels.loc[frame.key, field])
    audit.equal(name + ":fold1", frame.heldout_fold, np.ones(len(frame), dtype=int))
    audit.equal(name + ":chromosomes", frame.uncertainty_chromosome.astype(str), frame.key.str.split(":").str[0])
    saved = R.indexed(R.read(folder / "calibration_training.tsv.gz"))
    audit.equal(name + ":calibration_population", sorted(saved.key), sorted(train.key))
    R.check_source_rows(audit, name + ":training", saved, labels)
    coefficients = R.read(folder / "training_calibration.tsv").set_index("arm")
    native_arms, rows = [], []
    for length in lengths:
        for channel in ("atac", "dnase"):
            raw, source = f"native_{channel}_{length}", f"local_{channel}_liver"
            value = R.slope(tables[length].loc[train.key, source].to_numpy(float), train.beta_alt.to_numpy(float))
            audit.close(name + ":training_raw:" + raw, saved[raw], tables[length].loc[saved.key, source], atol=1e-12)
            audit.close(name + ":validation_raw:" + raw, frame[raw], tables[length].loc[frame.key, source], atol=1e-12)
            audit.close(name + ":slope:" + raw, [coefficients.loc[raw, "slope"]], [value])
            audit.close(name + ":intercept:" + raw, [coefficients.loc[raw, "intercept"]], [0], atol=0)
            audit.equal(name + ":training_n:" + raw, [coefficients.loc[raw, "training_n"]], [len(train)])
            audit.equal(name + ":training_folds:" + raw, [coefficients.loc[raw, "training_folds"]], ["2,3,4"])
            audit.close(name + ":calibrated:" + raw, frame[raw + "_calibrated"], value * frame[raw])
            native_arms.append(raw + "_calibrated")
            rows.append(dict(population=name, arm=raw, training_n=len(train), independent_slope=value))
    audit.equal(name + ":all_calibrations", sorted(coefficients.index), sorted(a.removesuffix("_calibrated") for a in native_arms))
    for arm in adapter_arms:
        audit.close(name + ":unchanged_adapter:" + arm, frame[arm], adapt.loc[frame.key, arm], atol=1e-12)
    arms = native_arms + adapter_arms
    if not np.isfinite(frame[["beta_alt", *arms]].to_numpy()).all():
        raise ValueError("Incomplete family; no successful-arm subset accepted")
    return frame, arms, rows, len(train)


def reconstruct(audit, folder, name, frame, arms, anchors, rank_pairs, original_n):
    y, p = frame.beta_alt.to_numpy(float), frame[arms].to_numpy(float)
    folds = frame.heldout_fold.to_numpy(int)
    mse = np.mean((p-y[:, None])**2, axis=0)
    rho = np.array([R.spearman(y, a) for a in p.T])
    pairs = [(a, b) for a in arms for b in anchors if a != b]
    results, points, rank_results = [], [], []
    for scheme, column, suffix in (("historical_1Mb", "block_1mb", ""),
                                   ("chromosome", "uncertainty_chromosome", "_chromosome_sensitivity")):
        prefix = name + suffix
        if json.loads((folder / (prefix + "_incomplete_arms.json")).read_text()):
            raise ValueError("Incomplete arm cannot disappear from matched family")
        metrics = R.indexed(R.read(folder / (prefix + "_performance.tsv")), "arm")
        audit.equal(prefix + ":models", sorted(metrics.arm), sorted(arms))
        metrics = metrics.loc[arms]
        code, units, weight = R.sampling_weights(frame[column], folds, 1000, 20260915)
        draws = R.error_samples(y, p, code, weight)
        errors = np.quantile(np.sqrt(draws), [.025, .975], axis=0)
        for field, expected in (("RMSE_beta_units", np.sqrt(mse)), ("MAE_beta_units", np.mean(abs(p-y[:, None]), axis=0)),
            ("pooled_signed_spearman", rho), ("macro_spearman", np.full(len(arms), np.nan)),
            ("RMSE_low95", errors[0]), ("RMSE_high95", errors[1])):
            audit.close(prefix + ":" + field, metrics[field], expected)
        for field, expected in (("n", len(y)), ("source_eligible_n", original_n), ("blocks", len(units)), ("folds", 1)):
            audit.equal(prefix + ":" + field, metrics[field], np.full(len(arms), expected))
        audit.close(prefix + ":coverage", metrics.coverage_fraction, np.full(len(arms), len(y)/original_n))
        table = R.read(folder / (prefix + "_paired_contrasts.tsv"))
        if table.duplicated(["arm", "comparator"]).any():
            raise ValueError("Repeated contrast")
        audit.equal(prefix + ":complete_family", sorted(zip(table.arm, table.comparator)), sorted(pairs))
        table = table.set_index(["arm", "comparator"]).loc[pairs]
        effects, cis, probabilities = [], [], []
        for arm, comparator in pairs:
            i, j = arms.index(arm), arms.index(comparator)
            benefit = mse[j]-mse[i]
            samples = draws[:, j]-draws[:, i]
            effects.append(benefit)
            cis.append(np.quantile(samples, [.025, .975]))
            probabilities.append((1+np.count_nonzero(np.abs(samples-benefit) >= abs(benefit)))/1001)
        cis, q = np.asarray(cis), R.adjust_bh(probabilities)
        for field, expected in (("MSE_improvement", effects), ("MSE_improvement_low95", cis[:, 0]),
            ("MSE_improvement_high95", cis[:, 1]), ("p_nominal_centered_block_bootstrap_MSE", probabilities),
            ("BH_q_complete_exploratory_MSE_family", q),
            ("rho_improvement", [rho[arms.index(a)]-rho[arms.index(b)] for a, b in pairs])):
            audit.close(prefix + ":" + field, table[field], expected)
        for field, expected in (("n", len(y)), ("planned_family_contrasts", len(pairs)),
                               ("estimable_family_contrasts", len(pairs)), ("resampling_units", len(units))):
            audit.equal(prefix + ":" + field, table[field], np.full(len(pairs), expected))
        selected = list(dict.fromkeys(a for pair in rank_pairs for a in pair))
        rank_draws, _ = R.rank_samples(y, p[:, [arms.index(a) for a in selected]], folds, code, weight)
        for j, arm in enumerate(selected):
            audit.close(prefix + ":rho_interval:" + arm, metrics.loc[arm, ["rho_low95", "rho_high95"]].to_numpy(float), R.interval(rank_draws[:, j]))
        for arm, comparator in rank_pairs:
            ci = R.interval(rank_draws[:, selected.index(arm)]-rank_draws[:, selected.index(comparator)])
            audit.close(prefix + ":rank_gain:" + arm + ":" + comparator,
                table.loc[(arm, comparator), ["rho_improvement_low95", "rho_improvement_high95"]].to_numpy(float), ci)
            rank_results.append(dict(population=name, scheme=scheme, arm=arm, comparator=comparator,
                rho_gain=rho[arms.index(arm)]-rho[arms.index(comparator)], CI_low=ci[0], CI_high=ci[1]))
        for i, arm in enumerate(arms):
            points.append(dict(population=name, scheme=scheme, arm=arm, n=len(y), source_eligible_n=original_n,
                RMSE_beta_units=np.sqrt(mse[i]), MAE_beta_units=np.mean(abs(p[:, i]-y)), signed_spearman=rho[i],
                RMSE_low95=errors[0, i], RMSE_high95=errors[1, i]))
        for i, (arm, comparator) in enumerate(pairs):
            results.append(dict(population=name, scheme=scheme, arm=arm, comparator=comparator, n=len(y),
                resampling_units=len(units), planned_family=len(pairs), MSE_improvement=effects[i],
                CI_low=cis[i, 0], CI_high=cis[i, 1], nominal_p=probabilities[i], BH_q=q[i]))
        print(f"Checked {prefix}: {len(y)} variants, {len(arms)} models, {len(pairs)} complete MSE contrasts", flush=True)
    return points, results, rank_results


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Use a compute node for scientific review")
    args.out.mkdir(parents=True, exist_ok=False)
    audit = R.Audit()
    started = time.monotonic()
    try:
        R.scientific_tests(audit)
        labels, tables, adapt, arms, coverage, full_meta, short_meta = prepare(audit, args)
        designs = [
            (args.lengths, "matched_native_lengths", [2048, 16384, 1048576],
             ["native_atac_2048_calibrated", "native_atac_1048576_calibrated"],
             [("adapter_r16_last3_symmetric", "native_atac_2048_calibrated"),
              ("native_atac_1048576_calibrated", "native_atac_2048_calibrated")]),
            (args.short, "same_length_2kb_eligible_population", [2048], ["native_atac_2048_calibrated"],
             [("adapter_r16_last3_symmetric", "native_atac_2048_calibrated")]),
        ]
        ready, slopes = [], []
        for folder, name, lengths, anchors, rank_pairs in designs:
            frame, all_arms, receipts, ntrain = calibration(audit, folder, name, lengths, labels, tables, adapt, arms)
            ready.append((folder, name, frame, all_arms, anchors, rank_pairs))
            slopes.extend(receipts)
            meta = full_meta if len(lengths) == 3 else short_meta
            audit.equal(name + ":receipt_counts", [ntrain, len(frame), len(adapt)],
                [meta["common_training_n"] if len(lengths) == 3 else meta["training_n"],
                 meta["common_validation_n"] if len(lengths) == 3 else meta["validation_n"], meta["original_adapter_validation_n"]])
        excluded = R.read(args.short / "excluded_adapter_validation_rows.tsv")
        audit.equal("short:excluded_validation", sorted(excluded.key), sorted(set(adapt.key)-set(ready[1][2].key)))
        report = json.loads((args.short / "length_run_metadata_review.json").read_text())
        audit.equal("actual_bootstrap_seed", [report["actual_bootstrap_seed"]], [20260915])
        pilots = [metadata(audit, args.lengths / f"pilot_{length}_run.json", length) for length in (2048, 16384)]
        conservative = sum(m["total_seconds"]/m["variants_attempted"]*(len(labels)-m["variants_scored"]) for m in pilots)*1.25+1500
        projection = json.loads((args.lengths / "projection.json").read_text())
        audit.close("host_projection_receipt", [report["conservative_host_inclusive_projection_seconds"]], [conservative], atol=1e-8)
        audit.equal("host_projection_disposition", [report["conservative_projection_fits_original_allocation"]],
                    [conservative < projection["remaining_seconds"]-300])
        all_points, all_contrasts, all_ranks = [], [], []
        for folder, name, frame, all_arms, anchors, rank_pairs in ready:
            points, contrasts, ranks = reconstruct(audit, folder, name, frame, all_arms, anchors, rank_pairs, len(adapt))
            all_points.extend(points); all_contrasts.extend(contrasts); all_ranks.extend(ranks)
        for filename, rows in (("performance.tsv", all_points), ("MSE_contrasts.tsv", all_contrasts),
                               ("key_rank_intervals.tsv", all_ranks), ("source_coverage.tsv", coverage),
                               ("training_calibration.tsv", slopes)):
            pd.DataFrame(rows).to_csv(args.out / filename, sep="\t", index=False)
        result = dict(status="pass", populations=2, model_point_rows=len(all_points), MSE_contrasts=len(all_contrasts),
            key_ranking_intervals=len(all_ranks), bootstrap_draws=1000, bootstrap_seed=20260915,
            eligible_validation_variants=len(adapt), eligible_training_variants=int(labels.heldout_fold.isin([2, 3, 4]).sum()),
            units="FastQTL_ALT_dosage_beta; MSE_beta_squared; native_raw_log2_track_sum_training_calibrated",
            common_population="shared_validation_and_calibration_coverage_all3lengths",
            unrestricted2kb_population="all_available2kb_validation_and_training_folds2_3_4",
            algorithm="prior_reviewer's_independent_lstsq_group_multiplicity_weighted_rank_helpers; no_producer_inferential_imports",
            rank_scope="all_point_correlations; three_key_contrast_intervals_both_units; remaining_rank_intervals_not_reconstructed",
            biological_n="138_original_cohort; per_variant_effective_n_unknown; variant_counts_are_not_participant_n",
            limits=["historically_significance_selected_Currin_leads", "one_previously_inspected_development_fold_and_one_adapter_seed",
                "common_population_adapter_training_larger_than_native_common_calibration",
                "no_refitting_or_model_selection_uncertainty", "1Mb_bins_not_verified_LD_units; only5validation_chromosomes",
                "chromosome_holdout_does_not_establish_pretraining_or_source_cohort_independence",
                "same2kb_input_still_distinct_native_track_readout_and_learned_scalar_head"],
            protected_external_outcomes_read=False, fold0_effects_numerically_loaded=False, fits_repeated=False, adoption=False,
            elapsed_seconds=time.monotonic()-started, versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__),
            source_hashes={str(p): R.digest(p) for p in [R.LABELS, args.adaptation, args.lengths / "complete.json", args.short / "complete.json",
                args.lengths / "native_2048.tsv", args.lengths / "native_16384.tsv", args.native / "native_1048576.tsv",
                Path(__file__), Path(R.__file__)]})
        (args.out / "checks.json").write_text(json.dumps(result, indent=2)+"\n")
        print(json.dumps(result, indent=2), flush=True)
    except Exception as exc:
        (args.out / "checks.json").write_text(json.dumps(dict(status="fail", error=repr(exc),
            completed_checks=len(audit.rows), disposition="stop_no_definition_changes_or_success_only_population"), indent=2)+"\n")
        raise
    finally:
        pd.DataFrame(audit.rows).to_csv(args.out / "check_differences.tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, default=BASE / "native")
    parser.add_argument("--lengths", type=Path, required=True)
    parser.add_argument("--short", type=Path, required=True)
    parser.add_argument("--adaptation", type=Path, default=BASE / "sweep_summary_21773124/matched_predictions.tsv.gz")
    parser.add_argument("--sweep", type=Path, default=BASE / "sweep_21772726")
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
