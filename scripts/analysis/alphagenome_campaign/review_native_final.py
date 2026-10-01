#!/usr/bin/env python3
"""Independent reconstruction of the four completed Currin development tables.

No producer inferential helper is imported. Source identities and training-only
zero-intercept calibrations are checked before numerical comparisons. MSE draws
use bootstrap multiplicities and grouped sufficient statistics; selected rank
intervals use frequency-weighted ranks without expanding sampled rows. These
check the saved, fixed predictions, not model fitting or selection uncertainty.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[3]
C2 = ROOT / "GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z"
LABELS = C2 / "inputs/currin_lead_labels.tsv.gz"
LEGACY = C2 / "tables/oof_predictions.tsv.gz"
SPECIALISTS = ROOT / "GWAS/finemapping/results/alphagenome_program/c1-endpoint2-v2-20260915T084559Z/tables/endpoint1_matched_tierA4.tsv.gz"
POPULATIONS = ("native_full_existing_leads", "fixed_representation_common_rows",
               "adaptation_single_split", "specialist_eligible_common_rows")
KEY_RANK_PAIRS = {
    POPULATIONS[0]: ("local_atac_liver__train_calibrated", "hyenadna__delta_ridge__ensemble"),
    POPULATIONS[1]: ("frozen_16384_target_shared_mlp64", "local_atac_liver__training_calibrated"),
    POPULATIONS[2]: ("adapter_r16_last3_symmetric", "native_atac_calibrated"),
    POPULATIONS[3]: ("chrombpnet_adult_hep_gse281367__training_calibrated", "local_atac_liver__training_calibrated"),
}


def read(path):
    return pd.read_csv(path, sep="\t")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def indexed(frame, key="key"):
    if frame[key].isna().any() or frame[key].duplicated().any():
        raise ValueError(f"Missing/duplicate {key}")
    return frame.set_index(key, drop=False)


class Audit:
    def __init__(self):
        self.rows = []

    def equal(self, label, observed, expected):
        a, b = np.asarray(observed), np.asarray(expected)
        if a.shape != b.shape or not np.array_equal(a, b):
            raise ValueError(f"Identity/coverage mismatch: {label}")
        self.rows.append(dict(check=label, n=a.size, max_absolute_difference=0.0, tolerance=0.0))

    def close(self, label, observed, expected, atol=2e-11, rtol=0.0):
        a, b = np.asarray(observed, dtype=float), np.asarray(expected, dtype=float)
        if a.shape != b.shape or not np.allclose(a, b, atol=atol, rtol=rtol, equal_nan=True):
            raise ValueError(f"Numerical mismatch: {label}")
        valid = np.isfinite(a) & np.isfinite(b)
        difference = float(np.max(np.abs(a[valid] - b[valid]))) if valid.any() else 0.0
        self.rows.append(dict(check=label, n=a.size, max_absolute_difference=difference,
                              tolerance=atol, relative_tolerance=rtol))


def spearman(y, p):
    if len(y) < 3 or np.ptp(y) == 0 or np.ptp(p) == 0:
        return np.nan
    return float(np.corrcoef(rankdata(y), rankdata(p))[0, 1])


def macro(y, p, folds):
    if set(folds) != set(range(5)):
        return np.nan
    rs = np.array([spearman(y[folds == f], p[folds == f]) for f in range(5)])
    return float(np.tanh(np.mean(np.arctanh(np.clip(rs, -.999999, .999999)))))


def adjust_bh(p):
    """BH via sorted reverse running minimum, independently implemented."""
    order = sorted(range(len(p)), key=lambda i: (p[i], i))
    adjusted = np.ones(len(p))
    running = 1.0
    for rank in range(len(order), 0, -1):
        i = order[rank - 1]
        running = min(running, float(p[i]) * len(p) / rank)
        adjusted[i] = running
    return adjusted


def slope(x, y):
    if not np.isfinite(x).all() or not np.isfinite(y).all() or float(x @ x) <= 1e-20:
        raise ValueError("Nonfinite or degenerate calibration training data")
    return float(np.linalg.lstsq(np.asarray(x)[:, None], y, rcond=None)[0][0])


def check_source_rows(audit, label, frame, labels, float32_labels=False):
    source = labels.loc[frame.key]
    audit.equal(label + ":canonical_variant", frame.key, source.key)
    for column in ("peak_id", "heldout_fold", "block_1mb", "ref", "alt", "pos_hg38"):
        if column in frame:
            audit.equal(label + ":" + column, frame[column], source[column])
    if "beta_alt" in frame:
        if float32_labels:
            # Producer serializes float32 targets; preserve this precision for
            # reconstruction, while checking the actual source-native labels.
            audit.close(label + ":source_beta_float32", frame.beta_alt, source.beta_alt,
                        atol=1.2e-7, rtol=1.2e-7)
        else:
            audit.close(label + ":source_beta", frame.beta_alt, source.beta_alt, atol=1e-12)


def check_calibration(audit, frame, raw_columns, suffix, receipt, label):
    """Reconstruct five training-only slopes with a separate least-squares solver."""
    rows = []
    for raw in raw_columns:
        for f in range(5):
            train = frame.heldout_fold.ne(f).to_numpy()
            test = ~train
            x = frame[raw].to_numpy(dtype=float)
            y = frame.beta_alt.to_numpy(dtype=float)
            value = slope(x[train], y[train])
            record = receipt.loc[receipt.raw.eq(raw) & receipt.heldout_fold.eq(f)]
            if len(record) != 1:
                raise ValueError(f"Missing calibration record: {label}/{raw}/{f}")
            audit.equal(f"{label}:{raw}:fold{f}:train_n", record.training_n, [int(train.sum())])
            audit.close(f"{label}:{raw}:fold{f}:slope", record.slope, [value])
            audit.close(f"{label}:{raw}:fold{f}:intercept", record.intercept, [0.0], atol=0)
            audit.close(f"{label}:{raw}:fold{f}:held_predictions", frame.loc[test, raw + suffix],
                        x[test] * value)
            rows.append(dict(population=label, raw=raw, heldout_fold=f, training_n=int(train.sum()),
                             independent_training_slope=value, intercept=0))
    return rows


def verify_inputs(args, audit):
    labels = read(LABELS)
    labels["key"] = labels.lead_variant_id.str.removeprefix("chr")
    labels = indexed(labels)
    canonical = (labels.chr.astype(str).str.removeprefix("chr") + ":" +
                 labels.pos_hg38.astype(str) + ":" + labels.ref + ":" + labels.alt)
    audit.equal("source:canonical_alleles", labels.key, canonical)
    audit.equal("source:population", [len(labels)], [32322])
    if labels.groupby("chr").heldout_fold.nunique().max() != 1:
        raise ValueError("A source chromosome crosses folds")
    native = args.model / "native"
    completion = json.loads((native / "completion.json").read_text())
    population = json.loads((native / "population.json").read_text())
    if (population["label_sha256"] != digest(LABELS) or
            population["measured_units"] != "source_FastQTL_ALT_dosage_beta" or
            population["native_score_units"] != "ALT_minus_REF_log2_1_plus_501bp_predicted_track_sum" or
            population["broader_outcome_independent_population"] or population["protected_outcomes_read"]):
        raise ValueError("Native population source identity/units/exposure declaration changed")
    scores = indexed(read(native / "native_1048576.tsv"))
    rejected = indexed(read(native / "native_1048576_rejected.tsv"))
    if completion["unresolved"] or set(scores.key) & set(rejected.key):
        raise ValueError("Unresolved/contradictory native attempt accounting")
    audit.equal("native:all_attempts", sorted(set(scores.key) | set(rejected.key)), sorted(labels.key))
    audit.equal("native:completion_counts", [len(labels), len(scores), len(rejected)],
                [completion["eligible"], completion["scored"], completion["rejected"]])
    check_source_rows(audit, "native_raw", scores, labels)
    check_source_rows(audit, "native_rejected", rejected, labels)
    for col, value in (("length_bp", 1048576), ("mask_width_bp", 501),
                       ("aggregation", "DIFF_LOG2_SUM"), ("card_tag", "l40s"), ("n_non_acgt", 0)):
        audit.equal("native_recipe:" + col, scores[col], np.repeat(value, len(scores)))
    for col in ("local_atac_liver", "local_dnase_liver"):
        if not np.isfinite(scores[col]).all():
            raise ValueError("Nonfinite native score")
    legacy = read(LEGACY)
    legacy["key"] = legacy.variant_id.str.removeprefix("chr")
    legacy = indexed(legacy)
    check_source_rows(audit, "legacy", legacy, labels)
    old = [c for c in legacy if c.endswith("__ensemble")]
    frozen = indexed(read(args.model / "frozen/comparisons/predictions.tsv.gz"))
    check_source_rows(audit, "fixed_original", frozen, labels, float32_labels=True)
    matched = np.ones(len(labels), dtype=bool)
    for length in (2048, 16384):
        folder = args.model / "frozen" / str(length)
        manifest = indexed(read(folder / "manifest.tsv"))
        audit.equal(f"manifest{length}:population_order", manifest.key, labels.key)
        check_source_rows(audit, f"manifest{length}", manifest, labels)
        meta = json.loads((folder / "complete.json").read_text())
        if not meta["complete"] or meta["card_tag"] != "l40s":
            raise ValueError("Fixed representation extraction incomplete/hardware inconsistent")
        with np.load(folder / "coverage.npz", allow_pickle=False) as masks:
            for pool in ("variant", "symmetric", "target"):
                matched &= masks[pool]
    audit.equal("fixed:six_way_coverage", sorted(frozen.key), sorted(labels.loc[matched, "key"]))
    fa = [c for c in frozen if c.startswith("frozen_")]
    if len(fa) != 12:
        raise ValueError("Expected all twelve fixed recipes")
    special = read(SPECIALISTS)
    special = indexed(special.loc[special.model_peak_overlap.eq(1)].copy())
    present = special.key.isin(labels.index)
    same_target = pd.Series(False, index=special.index)
    same_target.loc[present] = (special.loc[present, "peak_id"].to_numpy() ==
                                labels.loc[special.loc[present, "key"], "peak_id"].to_numpy())
    admitted = special.loc[present & same_target].copy()
    # TierA4 and C2 use the same floor(pos1 / 1e6) grouping, serialized
    # respectively as chrN~bin and chrN:bin. Check source semantics before
    # converting the delimiter for the canonical comparison; never regroup.
    specialist_blocks = ("chr" + admitted.key.str.split(":").str[0] + "~" +
                         (admitted.pos_hg38 // 1_000_000).astype(str))
    audit.equal("specialist:source_block_formula", admitted.block_1mb, specialist_blocks)
    admitted["block_1mb"] = admitted.block_1mb.str.replace("~", ":", regex=False)
    check_source_rows(audit, "specialist_exact_target", admitted, labels)
    raw_special = ["chrombpnet_adult_hep", "chrombpnet_adult_hep_gse281367", "borzoi_atac", "borzoi_dnase"]
    special_common = admitted.loc[admitted.key.isin(scores.index)].copy()
    special_common = special_common.loc[np.isfinite(special_common[raw_special]).all(axis=1)]
    exclusions = read(args.final / "specialist_other_target_or_population_exclusions.tsv")
    audit.equal("specialist:target_exclusions", sorted(exclusions.key), sorted(special.loc[~(present & same_target), "key"]))
    if len(read(args.final / "specialist_source_disagreements.tsv")):
        raise ValueError("Producer recorded contradictory exact-target effects/folds")

    runs = json.loads((args.sweep / "runs.json").read_text())
    completed = [r for r in runs if r.get("status") == "completed_fixed_exposure" and r.get("steps_requested") == 5000]
    if len(completed) != 12 or len({r["name"] for r in completed}) != 12:
        raise ValueError("Incomplete fixed-exposure sweep")
    original_validation = labels.loc[labels.heldout_fold.eq(1)]
    original_training = labels.loc[labels.heldout_fold.isin([2, 3, 4])]
    run_tables = {}
    for run in completed:
        path = Path(run["output"])
        pred = read(path / "validation_predictions.tsv")
        pred["key"] = pred.variant_id.str.removeprefix("chr")
        pred = indexed(pred)
        audit.equal(run["name"] + ":validation_population", sorted(pred.key), sorted(original_validation.key))
        audit.close(run["name"] + ":source_labels", pred.observed_beta,
                    labels.loc[pred.key, "beta_alt"].to_numpy(dtype=np.float32).astype(float), atol=1e-12)
        split = json.loads((path / "split.json").read_text())
        cost = json.loads((path / "feasibility.json").read_text())
        cfg = json.loads((path / "config.json").read_text())
        if (split["held_fold"] != 0 or split["validation_fold"] != 1 or split["held_outcomes_evaluated"] or
                split["training_n"] != len(original_training) or split["validation_n"] != len(original_validation) or
                split["training_native_slope"] != 0 or split["seed"] != 1103 or cfg["length"] != 2048 or
                cost["completed_steps"] != 5000 or cost["validation_examples"] != len(pred) or
                not all(cost[k] for k in ("frozen_parameters_unchanged", "running_statistics_unchanged", "save_reload_agreement"))):
            raise ValueError("Sweep exposure/state declaration inconsistent")
        if not np.isfinite(pred.predicted_beta).all():
            raise ValueError("Incomplete sweep predictions")
        run_tables[run["name"]] = pred

    frame_dict, arm_dict, anchor_dict, calibration_rows, coverage = {}, {}, {}, [], []
    for name in POPULATIONS:
        frame = indexed(read(args.final / (name + "_matched_predictions.tsv.gz")))
        check_source_rows(audit, name, frame, labels, float32_labels=(name == POPULATIONS[1]))
        audit.equal(name + ":chromosome", frame.uncertainty_chromosome.astype(str), frame.key.str.split(":").str[0])
        if name == POPULATIONS[0]:
            expected = set(scores.key)
            anchors = ["local_atac_liver__train_calibrated", "local_dnase_liver__train_calibrated"]
            arms = anchors + old
            receipt = read(native / "training_calibration.tsv").rename(columns={"arm": "raw", "training_slope": "slope"})
            calibration_rows += check_calibration(audit, frame, ["local_atac_liver", "local_dnase_liver"], "__train_calibrated", receipt, name)
        elif name == POPULATIONS[1]:
            expected = set(frozen.key) & set(scores.key)
            anchors = ["local_atac_liver__training_calibrated", "local_dnase_liver__training_calibrated"]
            arms = fa + anchors + old
            for arm in fa:
                audit.close(name + ":source_predictions:" + arm, frame[arm], frozen.loc[frame.key, arm], atol=1e-12)
            audit.close(name + ":serialized_targets", frame.beta_alt, frozen.loc[frame.key, "beta_alt"], atol=0)
            calibration_rows += check_calibration(audit, frame, ["local_atac_liver", "local_dnase_liver"], "__training_calibrated",
                read(args.final / "fixed_representation_native_calibration.tsv"), name)
        elif name == POPULATIONS[2]:
            expected = set(original_validation.key) & set(scores.key)
            anchors = ["native_atac_calibrated", "native_dnase_calibrated"]
            arms = list(run_tables) + ["allele_identity", "zero_effect"] + anchors
            audit.equal(name + ":fold1_only", frame.heldout_fold, np.ones(len(frame), dtype=int))
            for arm, pred in run_tables.items():
                audit.close(name + ":source_predictions:" + arm, frame[arm], pred.loc[frame.key, "predicted_beta"], atol=1e-12)
            nt = original_training.loc[original_training.key.isin(scores.index)]
            receipt = read(args.final / "sweep_reconstruction/native_training_calibration.tsv").set_index("output")
            for output in ("atac", "dnase"):
                raw = "local_" + output + "_liver"
                value = slope(scores.loc[nt.key, raw].to_numpy(), nt.beta_alt.to_numpy())
                audit.close(name + ":slope:" + output, [receipt.loc[output, "slope"]], [value])
                audit.equal(name + ":calibration_n:" + output, [receipt.loc[output, "training_n"]], [len(nt)])
                audit.close(name + ":intercept:" + output, [receipt.loc[output, "intercept"]], [0], atol=0)
                audit.close(name + ":calibrated:" + output, frame["native_" + output + "_calibrated"],
                            value * scores.loc[frame.key, raw])
                calibration_rows.append(dict(population=name, raw=raw, heldout_fold=1,
                    training_n=len(nt), independent_training_slope=value, intercept=0))
            eye = np.eye(4)
            lookup = dict(zip("ACGT", range(4)))
            def identity(rows):
                return np.stack([eye[lookup[a]] - eye[lookup[r]] for r, a in zip(rows.ref, rows.alt)])
            x = identity(original_training)
            weights = np.linalg.solve(x.T @ x + np.eye(4) * 1000, x.T @ original_training.beta_alt.to_numpy())
            audit.close(name + ":identity_ridge", frame.allele_identity, identity(labels.loc[frame.key]) @ weights)
            audit.close(name + ":zero", frame.zero_effect, np.zeros(len(frame)), atol=0)
        else:
            expected = set(special_common.key)
            raw = raw_special + ["local_atac_liver", "local_dnase_liver"]
            arms = [r + "__training_calibrated" for r in raw]
            anchors = arms[-2:]
            for col in raw_special:
                audit.close(name + ":source_scores:" + col, frame[col], special_common.loc[frame.key, col], atol=1e-12)
            calibration_rows += check_calibration(audit, frame, raw, "__training_calibrated",
                read(args.final / "specialist_training_calibration.tsv"), name)
        audit.equal(name + ":exact_population", sorted(frame.key), sorted(expected))
        if name != POPULATIONS[2]:
            for raw in ("local_atac_liver", "local_dnase_liver"):
                audit.close(name + ":raw_native:" + raw, frame[raw], scores.loc[frame.key, raw], atol=1e-12)
        if name in POPULATIONS[:2]:
            for arm in old:
                audit.close(name + ":legacy_predictions:" + arm, frame[arm], legacy.loc[frame.key, arm], atol=1e-12)
        if not np.isfinite(frame[["beta_alt", *arms]].to_numpy()).all():
            raise ValueError("Nonfinite predictions/labels: " + name)
        eligible = len(original_validation) if name == POPULATIONS[2] else len(labels)
        conditional = len(expected) if name == POPULATIONS[2] else len(labels)
        coverage.append(dict(population=name, matched_rows=len(frame), original_endpoint_eligible_rows=eligible,
            fraction_of_original_endpoint=len(frame)/eligible, printed_source_eligible_n=conditional,
            printed_denominator_scope="native_available_validation_only" if name == POPULATIONS[2] else "all32322existing_leads",
            fixed_short_length_common_rows=len(frozen) if name == POPULATIONS[1] else None,
            specialist_original_overlap_rows=len(special) if name == POPULATIONS[3] else None,
            specialist_exact_target_rows=len(admitted) if name == POPULATIONS[3] else None,
            arms=len(arms), planned_contrasts=sum(a != b for a in arms for b in anchors)))
        frame_dict[name], arm_dict[name], anchor_dict[name] = frame, arms, anchors
    return frame_dict, arm_dict, anchor_dict, calibration_rows, coverage


def sampling_weights(blocks, folds, draws, seed):
    codes, unique = pd.factorize(np.asarray(blocks, dtype=str), sort=True)
    unit_folds = np.array([np.unique(folds[codes == i]) for i in range(len(unique))], dtype=object)
    if any(len(v) != 1 for v in unit_folds):
        raise ValueError("Uncertainty unit crosses heldout folds")
    groups = [np.flatnonzero(np.array([int(v[0]) for v in unit_folds]) == f) for f in sorted(set(folds))]
    rng = np.random.default_rng(seed)
    weights = np.zeros((draws, len(unique)), dtype=np.int32)
    for d in range(draws):
        for group in groups:
            # integers indexes sorted units, exactly the declared sampling law;
            # grouped counts replace producer's repeated unit summation.
            sampled = group[rng.integers(0, len(group), len(group))]
            weights[d] += np.bincount(sampled, minlength=len(unique))
    return codes, unique, weights


def error_samples(y, predictions, codes, weights):
    squared = (predictions - y[:, None]) ** 2
    # pandas grouped sufficient statistics provide a distinct reduction path.
    sums = pd.DataFrame(squared).groupby(codes, sort=True).sum().to_numpy()
    counts = np.bincount(codes, minlength=weights.shape[1])
    denominator = weights @ counts
    return (weights @ sums) / denominator[:, None]


class FrequencyRanks:
    def __init__(self, values):
        self.unique, self.codes = np.unique(values, return_inverse=True)

    def ranks(self, weights):
        frequency = np.bincount(self.codes, weights=weights, minlength=len(self.unique))
        midrank = np.cumsum(frequency) - (frequency - 1) / 2
        return midrank[self.codes]


def frequency_correlation(y_ranks, p_ranks, weights):
    total = weights.sum()
    if total < 3:
        return np.nan
    ry, rp = y_ranks.ranks(weights), p_ranks.ranks(weights)
    ry -= np.sum(weights * ry) / total
    rp -= np.sum(weights * rp) / total
    variance = np.sum(weights * ry ** 2) * np.sum(weights * rp ** 2)
    if variance <= 0:
        return np.nan
    return float(np.sum(weights * ry * rp) / np.sqrt(variance))


def rank_samples(y, predictions, folds, codes, weights):
    y_ranks = FrequencyRanks(y)
    p_ranks = [FrequencyRanks(p) for p in predictions.T]
    pooled = np.full((len(weights), predictions.shape[1]), np.nan)
    macros = pooled.copy()
    for d, unit_weights in enumerate(weights):
        row_weights = unit_weights[codes]
        for j, ranked in enumerate(p_ranks):
            pooled[d, j] = frequency_correlation(y_ranks, ranked, row_weights)
            if set(folds) == set(range(5)):
                per_fold = [frequency_correlation(y_ranks, ranked, row_weights * (folds == f)) for f in range(5)]
                macros[d, j] = np.tanh(np.mean(np.arctanh(np.clip(per_fold, -.999999, .999999))))
    return pooled, macros


def interval(draws):
    finite = draws[np.isfinite(draws)]
    if len(finite) < int(.95 * len(draws)):
        return np.array([np.nan, np.nan])
    return np.quantile(finite, [.025, .975])


def scientific_tests(audit):
    y = np.array([2., 2., 5., 9., 0., 5.])
    p = np.array([7., 7., 1., 6., 2., 1.])
    w = np.array([3, 0, 1, 2, 4, 1])
    idx = np.repeat(np.arange(len(y)), w)
    audit.close("test:frequency_ranks_equal_explicit_replicates", [frequency_correlation(FrequencyRanks(y), FrequencyRanks(p), w)],
                [spearman(y[idx], p[idx])], atol=1e-14)
    x = np.arange(1., 7.)
    target = x * .7
    training = np.arange(4)
    first = slope(x[training], target[training])
    target[4:] = 1e9
    audit.close("test:held_outcome_calibration_invariance", [slope(x[training], target[training])], [first], atol=0)
    predictions = np.column_stack([y, y + 2, y + 2])
    code, _, weights = sampling_weights(["a", "a", "b", "b", "c", "c"], np.zeros(6, dtype=int), 20, 37)
    draws = error_samples(y, predictions, code, weights)
    audit.close("test:constant_error_bootstrap", draws[:, 1], np.full(20, 4.), atol=0)
    audit.close("test:identical_model_contrasts", draws[:, 1] - draws[:, 2], np.zeros(20), atol=0)
    arbitrary = np.column_stack([p, p * .3])
    grouped = error_samples(y, arbitrary, code, weights)
    expanded = []
    for draw in weights:
        rows = np.repeat(np.arange(len(y)), draw[code])
        expanded.append(np.mean((arbitrary[rows] - y[rows, None]) ** 2, axis=0))
    audit.close("test:grouped_errors_equal_expanded_rows", grouped, expanded, atol=1e-13)
    audit.close("test:BH_complete_family", adjust_bh([.01, .04, .03, 1]), [.04, .0533333333333333, .0533333333333333, 1], atol=1e-14)


def reconstruct(args, audit, frames, arms_by_population, anchors_by_population, coverage):
    analysis = json.loads((args.final / "analysis.json").read_text())
    if analysis["bootstrap_draws"] != 1000 or analysis["seed"] != 20260915:
        raise ValueError("Unexpected registered draw count/seed")
    points, contrasts_out, rank_out = [], [], []
    for name, frame in frames.items():
        arms, anchors = arms_by_population[name], anchors_by_population[name]
        y, p = frame.beta_alt.to_numpy(float), frame[arms].to_numpy(float)
        folds = frame.heldout_fold.to_numpy(int)
        mse = np.mean((p - y[:, None]) ** 2, axis=0)
        rho = np.array([spearman(y, x) for x in p.T])
        macros = np.array([macro(y, x, folds) for x in p.T])
        original_coverage = next(r for r in coverage if r["population"] == name)
        pairs = [(a, b) for a in arms for b in anchors if a != b]
        for scheme, column, suffix in (("historical_1Mb", "block_1mb", ""),
                                       ("chromosome", "uncertainty_chromosome", "_chromosome_sensitivity")):
            prefix = name + suffix
            label = prefix + ":"
            incomplete = json.loads((args.final / (prefix + "_incomplete_arms.json")).read_text())
            if incomplete:
                raise ValueError("Current complete comparison contains incomplete arms; no success-only reconstruction")
            reported = indexed(read(args.final / (prefix + "_performance.tsv")), "arm")
            audit.equal(label + "model_family", sorted(reported.arm), sorted(arms))
            reported = reported.loc[arms]
            code, unique, weight = sampling_weights(frame[column], folds, 1000, 20260915)
            draws = error_samples(y, p, code, weight)
            for col, values in (("RMSE_beta_units", np.sqrt(mse)), ("MAE_beta_units", np.mean(np.abs(p-y[:, None]), axis=0)),
                                ("pooled_signed_spearman", rho), ("macro_spearman", macros),
                                ("RMSE_low95", np.quantile(np.sqrt(draws), .025, axis=0)),
                                ("RMSE_high95", np.quantile(np.sqrt(draws), .975, axis=0))):
                audit.close(label + col, reported[col], values)
            for col, value in (("n", len(y)), ("blocks", len(unique)), ("folds", len(set(folds))),
                               ("source_eligible_n", original_coverage["printed_source_eligible_n"])):
                audit.equal(label + col, reported[col], np.full(len(arms), value))
            audit.close(label + "printed_coverage", reported.coverage_fraction,
                        np.full(len(arms), len(y)/original_coverage["printed_source_eligible_n"]))
            table = read(args.final / (prefix + "_paired_contrasts.tsv"))
            if table.duplicated(["arm", "comparator"]).any():
                raise ValueError("Duplicated exploratory contrast")
            audit.equal(label + "complete_pairs", sorted(zip(table.arm, table.comparator)), sorted(pairs))
            table = table.set_index(["arm", "comparator"]).loc[pairs]
            point, cis, nominal = [], [], []
            for a, b in pairs:
                i, j = arms.index(a), arms.index(b)
                difference = draws[:, j] - draws[:, i]
                effect = mse[j] - mse[i]
                ci = np.quantile(difference, [.025, .975])
                prob = (1 + np.count_nonzero(np.abs(difference - effect) >= abs(effect))) / 1001
                point.append(effect); cis.append(ci); nominal.append(prob)
            cis = np.asarray(cis)
            q = adjust_bh(nominal)
            for col, values in (("MSE_improvement", point), ("MSE_improvement_low95", cis[:, 0]),
                                ("MSE_improvement_high95", cis[:, 1]), ("p_nominal_centered_block_bootstrap_MSE", nominal),
                                ("BH_q_complete_exploratory_MSE_family", q),
                                ("rho_improvement", [rho[arms.index(a)]-rho[arms.index(b)] for a, b in pairs])):
                audit.close(label + col, table[col], values)
            for col, value in (("n", len(y)), ("resampling_units", len(unique)),
                               ("planned_family_contrasts", len(pairs)), ("estimable_family_contrasts", len(pairs))):
                audit.equal(label + col, table[col], np.full(len(pairs), value))
            key_pair = KEY_RANK_PAIRS[name]
            rank_indices = [arms.index(a) for a in key_pair]
            rd, md = rank_samples(y, p[:, rank_indices], folds, code, weight)
            for j, arm in enumerate(key_pair):
                audit.close(label + arm + ":rho_interval", reported.loc[arm, ["rho_low95", "rho_high95"]].to_numpy(float), interval(rd[:, j]))
                audit.close(label + arm + ":macro_interval", reported.loc[arm, ["macro_low95", "macro_high95"]].to_numpy(float), interval(md[:, j]))
            if key_pair in pairs:
                rank_comparison, sign = key_pair, 1
            else:
                rank_comparison, sign = key_pair[::-1], -1
            rank_ci = interval(sign * (rd[:, 0] - rd[:, 1]))
            audit.close(label + "key_rho_difference_interval",
                table.loc[rank_comparison, ["rho_improvement_low95", "rho_improvement_high95"]].to_numpy(float), rank_ci)
            rank_out.append(dict(population=name, scheme=scheme, arm=rank_comparison[0], comparator=rank_comparison[1],
                rho_improvement_low95=rank_ci[0], rho_improvement_high95=rank_ci[1],
                valid_draws=int(np.isfinite(rd).all(axis=1).sum()), interval_scope="post_result_key_contrast_reconstruction_not_selection"))
            for i, arm in enumerate(arms):
                points.append(dict(population=name, scheme=scheme, arm=arm, n=len(y), RMSE_beta_units=float(np.sqrt(mse[i])),
                    MAE_beta_units=float(np.mean(np.abs(p[:, i]-y))), pooled_signed_spearman=rho[i], macro_spearman=macros[i],
                    RMSE_low95=float(np.quantile(np.sqrt(draws[:, i]), .025)), RMSE_high95=float(np.quantile(np.sqrt(draws[:, i]), .975))))
            for i, (a, b) in enumerate(pairs):
                contrasts_out.append(dict(population=name, scheme=scheme, arm=a, comparator=b, n=len(y),
                    resampling_units=len(unique), planned_family=len(pairs), MSE_improvement=point[i],
                    MSE_improvement_low95=cis[i, 0], MSE_improvement_high95=cis[i, 1], p_nominal=nominal[i], BH_q=q[i]))
            print(f"Checked {prefix}: {len(y)} variants, {len(arms)} models, {len(pairs)} MSE contrasts", flush=True)
    return points, contrasts_out, rank_out


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Execute scientific reconstruction on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    audit = Audit()
    try:
        scientific_tests(audit)
        frames, arms, anchors, calibration, coverage = verify_inputs(args, audit)
        # Identity/calibration failures above abort before any result uncertainty.
        pd.DataFrame(calibration).to_csv(args.out / "independent_training_calibration.tsv", sep="\t", index=False)
        pd.DataFrame(coverage).to_csv(args.out / "coverage_and_denominators.tsv", sep="\t", index=False)
        points, contrasts, ranks = reconstruct(args, audit, frames, arms, anchors, coverage)
        pd.DataFrame(points).to_csv(args.out / "reconstructed_performance.tsv", sep="\t", index=False)
        pd.DataFrame(contrasts).to_csv(args.out / "reconstructed_MSE_contrasts.tsv", sep="\t", index=False)
        pd.DataFrame(ranks).to_csv(args.out / "reconstructed_key_rank_intervals.tsv", sep="\t", index=False)
        result = dict(status="pass", original_source_labels=32322, populations=4, uncertainty_schemes=2,
            model_point_rows=len(points), MSE_contrasts=len(contrasts), key_rank_contrasts=len(ranks),
            all_model_point_correlations_checked=True, all_model_RMSE_intervals_checked=True,
            rank_interval_scope="two_preidentified_models_and_one_pair_per_population_both_units; other_rank_intervals_not_reconstructed",
            algorithm="independent_lstsq_calibration; pandas_group_sums_and_bootstrap_multiplicities; weighted_rank_ties",
            producer_inferential_helpers_imported=False, bootstrap_seed=20260915, bootstrap_draws=1000,
            units="source_FastQTL_ALT_dosage_beta; MSE_in_beta_squared; native_raw_track_log2_sum_scores_training_calibrated",
            source_population="significance_selected_historical_Currin_leads_not_outcome_independent_all_variant_peak_population",
            participant_n="138_source_cohort_participants; per_variant_effective_association_n_unknown; reported_n_is_variants",
            exposure_limits=["adaptation_length2048_vs_native1048576; practical_comparison_not_same_length_inferiority",
                "adaptation_one_seed_fold1_validation_training_folds2_3_4_fold0_closed; fixed_recipes_fivefold_development",
                "legacy_DNA_five_seed_ensembles_and_larger_original_training_exposure",
                "frozen_heads_training28600_rows_native_common27138; training_coverage_asymmetry_remains",
                "uncertainty_conditions_on_saved_fits_no_refitting_or_model_selection_uncertainty",
                "historical1Mb_bins_not_verified_LD_independence; whole_chr_sensitivity_few_units",
                "zero_output_has_defined_MSE_but_undefined_ranking_retained",
                "adaptation_printed_100percent_coverage_conditional6441_of6834_original_validation_rows"],
            protected_external_outcomes_read=False, model_fitting_reperformed=False, adopted=False,
            elapsed_seconds=time.monotonic()-start,
            versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
            sources={str(p): digest(p) for p in [LABELS, LEGACY, SPECIALISTS, args.final / "analysis.json",
                args.model / "native/native_1048576.tsv", args.model / "frozen/comparisons/predictions.tsv.gz",
                args.sweep / "runs.json", Path(__file__)]})
        (args.out / "checks.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
        print(json.dumps(result, indent=2), flush=True)
    except Exception as exc:
        (args.out / "checks.json").write_text(json.dumps(dict(status="fail", error=repr(exc),
            disposition="stop_identity_calibration_or_numerical_mismatch_no_definitions_changed",
            passed_checks=len(audit.rows)), indent=2)+"\n")
        raise
    finally:
        pd.DataFrame(audit.rows).to_csv(args.out / "check_differences.tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--sweep", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
