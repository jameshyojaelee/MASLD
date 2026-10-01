#!/usr/bin/env python3
"""Held-participant prediction under fixed missing-gene and batch scenarios."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from chromatin_acquisition import (FIX, WEIGHTS, SEED as BASE_SEED, panel,
                                   logcpm, concentration, rrr_fit_predict,
                                   local_fit_predict)

SEED = 20260923


def query_transform(counts, present, train_mean, train_sd, marginal, expected_mass=None):
    """Observed-gene CPM or its training-reference mass correction; absent-gene mean."""
    if expected_mass is not None and marginal:
        raise ValueError("Reference mass correction and batch marginal mapping are separate routes")
    if expected_mass is not None and (not np.isfinite(expected_mass) or not 0 < expected_mass <= 1):
        raise ValueError("Expected observed RNA mass must be positive and at most one")
    if expected_mass is None or expected_mass == 1:
        values = logcpm(counts[:, present])
    else:
        observed = counts[:, present]
        library = observed.sum(1)
        if np.any(library <= 0):
            raise ValueError("Reference mass correction requires positive observed RNA totals")
        values = np.log2(observed/library[:, None]*1e6*expected_mass+1)
    out = np.tile(train_mean, (len(counts), 1))
    if marginal:
        sd = values.std(0)
        values = train_mean[present] + (values-values.mean(0))/np.where(sd < 1e-12, 1, sd)*train_sd[present]
    out[:, present] = values
    return out


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    args.out.mkdir(parents=True, exist_ok=False)
    weights = np.load(WEIGHTS, allow_pickle=True)
    selected, keys, cis = panel(weights)
    ids = pd.read_csv(FIX/"molecular/participant_axis.tsv", sep="\t").participant_id.astype(str).to_numpy()
    fold = pd.read_csv(FIX/"folds/participant_outer_folds.tsv", sep="\t").set_index("participant_id").loc[ids, "outer_fold"].to_numpy(int)
    raw = np.asarray(np.load(FIX/"molecular/rna_values.npy", mmap_mode="r"), float)
    h3 = np.asarray(np.load(FIX/"molecular/h3k27ac_counts.npy", mmap_mode="r"), float)
    x, yraw, descriptors = logcpm(raw), logcpm(h3), concentration(h3)
    genes = raw.shape[1]
    all_present = np.ones(genes, bool)
    rows, checks, mass_coverage, masks, mask_identities = [], [], [], [], []
    routes = ("raw", "marginal", "reference_mass") if args.reference_mass_correction else ("raw", "marginal")
    for outer in sorted(set(fold)):
        tr, te = np.flatnonzero(fold != outer), np.flatnonzero(fold == outer)
        coefficient = np.linalg.lstsq(np.column_stack([np.ones(len(tr)), descriptors[tr]]), yraw[tr], rcond=None)[0]
        y = yraw-np.column_stack([np.ones(len(yraw)), descriptors])@coefficient
        mean, sd = x[tr].mean(0), x[tr].std(0)
        # Each training participant contributes equally, independently of library size.
        expected_gene_mass = (raw[tr]/raw[tr].sum(1)[:, None]).mean(0)
        rng = np.random.default_rng(SEED+int(outer))
        patterns = {"complete": all_present}
        for fraction in (.9, .75, .6):
            present = np.zeros(genes, bool)
            present[rng.choice(genes, round(fraction*genes), replace=False)] = True
            patterns[f"random_retained_{fraction:.2f}"] = present
        abundance = raw[tr].mean(0)
        size = round(.1*genes)
        # Many low-expression genes tie at zero; make their order independent
        # of platform-specific quicksort tie ordering in future runs.
        abundance_order = np.lexsort((np.arange(genes), abundance))
        for label, missing in (("training_top10pct_absent", abundance_order[-size:]),
                               ("training_bottom10pct_absent", abundance_order[:size])):
            present = all_present.copy()
            present[missing] = False
            patterns[label] = present
        blocks, meta = [], []
        pattern_mass = {}
        for scenario, present in patterns.items():
            masks.append(present.copy())
            mask_identities.append(f"fold{outer}|{scenario}")
            expected_mass = 1.0 if present.all() else float(expected_gene_mass[present].sum())
            pattern_mass[scenario] = expected_mass
            mass_coverage.append({"outer_fold": int(outer), "scenario": scenario,
                                  "n_present_genes": int(present.sum()), "n_genes": genes,
                                  "gene_coverage": float(present.mean()),
                                  "expected_observed_RNA_mass_fraction": expected_mass})
            for route in routes:
                blocks.append(query_transform(raw[te], present, mean, sd, route == "marginal",
                                               expected_mass if route == "reference_mass" else None))
                meta.append((scenario, route, te, present.mean()))
        probe = np.array(sorted(te, key=lambda j: ids[j]))[:10]
        for route in routes:
            blocks.append(query_transform(raw[probe], all_present, mean, sd, route == "marginal",
                                           1.0 if route == "reference_mass" else None))
            meta.append(("complete_batch10", route, probe, 1.0))
        for scale in (.5, 2.0):
            for route in routes:
                blocks.append(query_transform(raw[te]*scale, all_present, mean, sd, route == "marginal",
                                               1.0 if route == "reference_mass" else None))
                meta.append((f"library_scale_{scale:g}", route, te, 1.0))
        expanded = np.concatenate([x]+blocks)
        query = np.arange(len(x), len(expanded))
        # The same baseline is fitted once per outer fold for all perturbations.
        parts = np.array_split(np.random.default_rng(BASE_SEED+len(tr)).permutation(tr), 3)
        oof = np.zeros((len(tr), len(selected)))
        locations = {int(p): j for j, p in enumerate(tr)}
        for part in parts:
            _, pred = rrr_fit_predict(x, y, np.setdiff1d(tr, part), part, selected)
            oof[[locations[int(p)] for p in part]] = pred
        _, pred_global = rrr_fit_predict(expanded, y, tr, query, selected)
        local = local_fit_predict(expanded, y[tr][:, selected]-oof, tr, query, cis)
        prediction_blocks = {}
        offset = 0
        for scenario, route, people, coverage in meta:
            take = slice(offset, offset+len(people))
            global_block = pred_global[take].copy()
            if route == "marginal":
                scale = global_block.std(0)
                # Existing marginal route uses training cross-fitted moments,
                # avoiding the larger variance of in-sample global predictions.
                global_block = oof.mean(0)+(global_block-global_block.mean(0))/np.where(scale < 1e-12, 1, scale)*oof.std(0)
            prediction_blocks[scenario, route] = (people, global_block, global_block+local[take])
            offset += len(people)
        for (scenario, route), (people, global_block, corrected) in prediction_blocks.items():
            target = y[people][:, selected]
            ref = np.sum((target-y[tr][:, selected].mean(0))**2, axis=1)
            baseline = prediction_blocks["complete", route]
            lookup = {int(p): j for j, p in enumerate(baseline[0])}
            position = np.array([lookup[int(p)] for p in people])
            for model, prediction, base in (("global_RRR", global_block, baseline[1][position]),
                                             ("RRR_local", corrected, baseline[2][position])):
                loss = np.sum((target-prediction)**2, axis=1)
                shift = np.sum((prediction-base)**2, axis=1)
                for j, person in enumerate(people):
                    rows.append({"outer_fold": int(outer), "participant_id": ids[person],
                                 "scenario": scenario, "route": route, "model": model,
                                 "coverage": float(patterns.get(scenario, all_present).mean()),
                                 "expected_observed_RNA_mass_fraction": pattern_mass.get(scenario, 1.0),
                                 "reference_SSE": float(ref[j]), "SSE": float(loss[j]),
                                 "prediction_shift_SSE": float(shift[j]), "n_regions": len(selected)})
                if scenario.startswith("library_scale_") or (scenario == "complete_batch10" and route in ("raw", "reference_mass")):
                    maximum = float(np.max(np.abs(prediction-base)))
                    if maximum > 1e-8:
                        raise ValueError(f"Numerical invariance failed: {scenario}/{route}/{model}")
                    checks.append({"outer_fold": int(outer), "scenario": scenario, "route": route,
                                   "model": model, "max_absolute_difference": maximum})
                if scenario == "complete" and route == "reference_mass":
                    raw_prediction = prediction_blocks["complete", "raw"][1 if model == "global_RRR" else 2]
                    maximum = float(np.max(np.abs(prediction-raw_prediction)))
                    if maximum > 1e-8:
                        raise ValueError(f"Complete reference-mass input differs from raw: {model}")
                    checks.append({"outer_fold": int(outer), "scenario": scenario, "route": route,
                                   "model": model, "comparison": "complete_reference_mass_vs_raw",
                                   "max_absolute_difference": maximum})
        print(json.dumps({"outer_fold": int(outer), "scenarios": len(meta), "participants": len(te)}), flush=True)
    errors = pd.DataFrame(rows)
    errors.to_csv(args.out/"held_participant_errors.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(mass_coverage).to_csv(args.out/"reference_mass_coverage.tsv", sep="\t", index=False)
    gene_ids = pd.read_csv(FIX/"molecular/rna_feature_axis.tsv", sep="\t").stable_gene_id.astype(str).to_numpy()
    np.savez_compressed(args.out/"missing_gene_masks.npz", present=np.stack(masks),
                        fold_scenario=np.array(mask_identities), gene_id=gene_ids)
    rng = np.random.default_rng(SEED)
    summaries = []
    for (scenario, route, model), frame in errors.groupby(["scenario", "route", "model"]):
        frame = frame.sort_values("participant_id")
        if frame.participant_id.duplicated().any():
            raise ValueError("Duplicate held participant")
        draw = rng.integers(0, len(frame), size=(5000, len(frame)))
        loss, ref = frame.SSE.to_numpy(), frame.reference_SSE.to_numpy()
        boot = 1-loss[draw].sum(1)/ref[draw].sum(1)
        summaries.append({"scenario": scenario, "route": route, "model": model,
                          "n_participants": len(frame), "coverage": float(frame.coverage.mean()),
                          "expected_observed_RNA_mass_fraction": float(frame.expected_observed_RNA_mass_fraction.mean()),
                          "skill": float(1-loss.sum()/ref.sum()),
                          "ci_low95": float(np.quantile(boot, .025)), "ci_high95": float(np.quantile(boot, .975)),
                          "profile_RMS_shift_vs_full_fold": float(np.sqrt(frame.prediction_shift_SSE.sum()/(len(frame)*len(selected))))})
    pd.DataFrame(summaries).to_csv(args.out/"summary.tsv", sep="\t", index=False)
    paired = []
    for (scenario, route, model), frame in errors.groupby(["scenario", "route", "model"]):
        if scenario == "complete":
            continue
        base = errors.loc[errors.scenario.eq("complete") & errors.route.eq(route) & errors.model.eq(model)]
        matched = frame.merge(base, on="participant_id", suffixes=("", "_complete"), validate="one_to_one")
        if len(matched) != len(frame) or not np.allclose(matched.reference_SSE, matched.reference_SSE_complete):
            raise ValueError("Perturbation contrast differs in evaluation population or target")
        gain = (matched.SSE_complete-matched.SSE).to_numpy()
        ref = matched.reference_SSE.to_numpy()
        draw = rng.integers(0, len(matched), size=(5000, len(matched)))
        boot = gain[draw].sum(1)/ref[draw].sum(1)
        paired.append({"scenario": scenario, "contrast": route+"_perturbed_minus_complete", "model": model,
                       "n_paired_participants": len(matched), "delta_skill": float(gain.sum()/ref.sum()),
                       "ci_low95": float(np.quantile(boot, .025)), "ci_high95": float(np.quantile(boot, .975))})
    for (scenario, model), frame in errors.groupby(["scenario", "model"]):
        comparisons = [("marginal", "raw")]
        if args.reference_mass_correction:
            comparisons += [("reference_mass", "raw"), ("reference_mass", "marginal")]
        for candidate, comparator in comparisons:
            first = frame.loc[frame.route.eq(comparator)]
            second = frame.loc[frame.route.eq(candidate)]
            matched = first.merge(second, on="participant_id", suffixes=("_comparator", "_candidate"), validate="one_to_one")
            if len(matched) != len(first) or not np.allclose(matched.reference_SSE_comparator, matched.reference_SSE_candidate):
                raise ValueError("Transformation contrast differs in evaluation population or target")
            gain = (matched.SSE_comparator-matched.SSE_candidate).to_numpy()
            ref = matched.reference_SSE_comparator.to_numpy()
            draw = rng.integers(0, len(matched), size=(5000, len(matched)))
            boot = gain[draw].sum(1)/ref[draw].sum(1)
            paired.append({"scenario": scenario, "contrast": candidate+"_minus_"+comparator, "model": model,
                           "n_paired_participants": len(matched), "delta_skill": float(gain.sum()/ref.sum()),
                           "ci_low95": float(np.quantile(boot, .025)), "ci_high95": float(np.quantile(boot, .975))})
    pd.DataFrame(paired).to_csv(args.out/"paired_sensitivity_contrasts.tsv", sep="\t", index=False)
    (args.out/"numerical_invariants.json").write_text(json.dumps(checks, indent=2)+"\n")
    design_name = "reference_mass_robustness_addendum.json" if args.reference_mass_correction else "held_input_robustness_addendum.json"
    (args.out/"scope.json").write_text(json.dumps({"seed": SEED, "n_participants": len(ids),
        "global_response_regions": y.shape[1], "local_and_evaluation_regions": len(selected),
        "design_sha256": hashlib.sha256(Path(__file__).with_name(design_name).read_bytes()).hexdigest(),
        "reference_mass_correction": args.reference_mass_correction,
        "reference_mass_assumption": "observed-gene relative composition compatible with the training mean of per-participant RNA proportions; same missing genes can carry different RNA mass in a new person",
        "evidence_state": "adaptive single-cohort development sensitivity after observed missing-gene renormalization failure" if args.reference_mass_correction else "prespecified within-development input sensitivity",
        "uncertainty": "paired participant sampling conditional on fitted models; descriptive intervals without multiplicity correction; no retraining uncertainty or external transport",
        "model_scope": "fixed rank12 bounded recipe, not exact released outer-fold tuned recipe",
        "depth_limit": "fractional inputs preserved; library scaling is not read sampling",
        "slurm_job_id": os.environ["SLURM_JOB_ID"]}, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reference-mass-correction", action="store_true",
                        help="Add the adaptive training-reference RNA-mass transformation comparison")
    main(parser.parse_args())
