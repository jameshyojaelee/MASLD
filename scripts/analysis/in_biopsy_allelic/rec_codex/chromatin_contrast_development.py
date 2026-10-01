#!/usr/bin/env python3
"""Nested participant development of a fixed within-histology loss comparison."""
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from weighted_chromatin import donor_weights, fit_predict

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "scripts/analysis/masld_two_models"))
from chromatin_acquisition import FIX, JOIN, WEIGHTS, concentration, logcpm, panel

DESIGN = Path(__file__).with_name("chromatin_contrast_design.json")
EXPECTED_PAIRS = ROOT / "GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/within_histology_21849320/pairs_fixed_from_metadata.tsv"
EXPECTED_PANEL = ROOT / "GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/chromatin_acquisition_21849197/fixed_region_panel.tsv"
FIELDS = ("steatosis", "ballooning", "lobular_inflammation", "fibrosis", "sex")
ARMS = {"unweighted_ridge": None, "unweighted_rank12": 12,
        "weighted_ridge": None, "weighted_rank12": 12}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def paired_indices(indices, strata, ids):
    groups = {}
    for i in indices:
        groups.setdefault(strata[i], []).append(int(i))
    pairs = []
    for key in sorted(groups):
        ordered = sorted(groups[key], key=lambda i: hashlib.sha256(f"20260923|{ids[i]}".encode()).digest())
        pairs.extend(zip(ordered[::2], ordered[1::2]))
    return np.asarray(pairs, int).reshape(-1, 2)


def residual_target(raw, descriptors, train):
    z = np.column_stack([np.ones(len(raw)), descriptors])
    coefficient = np.linalg.lstsq(z[train], raw[train], rcond=None)[0]
    return raw - z @ coefficient


def losses(prediction, target, query, pairs):
    if not np.isfinite(prediction).all() or not np.isfinite(target).all():
        raise ValueError("Nonfinite prediction or residual target")
    locations = {int(i): j for j, i in enumerate(query)}
    a = np.array([locations[int(i)] for i in pairs[:, 0]])
    b = np.array([locations[int(i)] for i in pairs[:, 1]])
    pair_error = (prediction[a] - prediction[b]) - (target[pairs[:, 0]] - target[pairs[:, 1]])
    result = float(np.sum(pair_error ** 2)), float(np.sum((prediction - target[query]) ** 2))
    if not all(np.isfinite(value) for value in result):
        raise ValueError("Nonfinite candidate loss")
    return result


def choose(candidates, profile_limit=None):
    eligible = [case for case in candidates if profile_limit is None or case["profile_sse"] <= profile_limit]
    if not eligible:
        raise ValueError("No candidate satisfies the prespecified profile safeguard")
    best = min(case["pair_sse"] for case in eligible)
    tied = [case for case in eligible if case["pair_sse"] <= best + 1e-10 * max(1.0, abs(best))]
    return min(tied, key=lambda case: (case["strength"], -case["penalty"]))


def main():
    job = os.environ.get("SLURM_JOB_ID")
    if not job:
        raise RuntimeError("CPU SLURM allocation required")
    design = json.loads(DESIGN.read_text())
    out = Path(os.environ["CODEX_REC_OUTPUT"]) / f"chromatin_contrast_development_{job}"
    out.mkdir(exist_ok=False)
    (out / "fixed_design.json").write_text(json.dumps(design, indent=2) + "\n")
    axis = pd.read_csv(FIX / "molecular/participant_axis.tsv", sep="\t")
    ids = axis.participant_id.astype(str).to_numpy()
    if len(ids) != 99 or len(set(ids)) != 99:
        raise ValueError("Expected 99 unique development donors")
    fold_table = pd.read_csv(FIX / "folds/participant_outer_folds.tsv", sep="\t").set_index("participant_id")
    if not fold_table.index.is_unique or set(fold_table.index.astype(str)) != set(ids):
        raise ValueError("Fold donor axis differs")
    fold = fold_table.loc[ids, "outer_fold"].to_numpy(int)
    join_table = pd.read_csv(JOIN, sep="\t").set_index("participant_id")
    if not join_table.index.is_unique or set(join_table.index.astype(str)) != set(ids) or set(fold) != set(range(5)):
        raise ValueError("Authoritative donor axis or fold universe differs")
    join = join_table.loc[ids]
    if not np.array_equal(fold, join.outer_fold.to_numpy(int)) or join[list(FIELDS)].isna().any().any():
        raise ValueError("Fold or exact matching metadata mismatch")
    strata = [tuple(str(value) for value in row) for row in join[list(FIELDS)].to_numpy()]
    expected_pairs = pd.read_csv(EXPECTED_PAIRS, sep="\t")
    planned_pairs = np.concatenate([paired_indices(np.flatnonzero(fold == f), strata, ids) for f in range(5)])
    planned_ids = {(ids[a], ids[b]) for a, b in planned_pairs}
    if len(planned_pairs) != 22 or len(set(planned_pairs.ravel())) != 44 or planned_ids != set(zip(expected_pairs.participant_a, expected_pairs.participant_b)):
        raise ValueError("Exact existing evaluation pair identities differ")
    with np.load(WEIGHTS, allow_pickle=True) as weights:
        selected, keys, _ = panel(weights)
    expected_panel = pd.read_csv(EXPECTED_PANEL, sep="\t")
    if not np.array_equal(selected, expected_panel.region_index.to_numpy(int)) or not np.array_equal(keys, expected_panel.region_key.astype(str).to_numpy()):
        raise ValueError("Exact existing panel differs")
    pd.DataFrame({"region_index": selected, "region_key": keys}).to_csv(out / "fixed_regions.tsv", sep="\t", index=False)
    input_paths = (JOIN, WEIGHTS, EXPECTED_PAIRS, EXPECTED_PANEL, FIX / "molecular/participant_axis.tsv",
                   FIX / "folds/participant_outer_folds.tsv", FIX / "molecular/rna_values.npy",
                   FIX / "molecular/h3k27ac_counts.npy", ROOT / "scripts/analysis/masld_two_models/chromatin_acquisition.py")
    receipt = {"inputs": {str(path): sha256(path) for path in input_paths},
               "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
               "design_sha256": sha256(DESIGN), "panel_sha256": sha256(out / "fixed_regions.tsv"),
               "slurm_job_id": job}
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    rna = np.asarray(np.load(FIX / "molecular/rna_values.npy", mmap_mode="r"), float)
    h3 = np.asarray(np.load(FIX / "molecular/h3k27ac_counts.npy", mmap_mode="r"), float)
    if rna.shape[0] != 99 or h3.shape != (99, 96460):
        raise ValueError("Molecular donor/region axes differ")
    if not np.isfinite(rna).all() or not np.isfinite(h3).all() or np.any(rna < 0) or np.any(h3 < 0):
        raise ValueError("Nonfinite or negative molecular values")
    if np.any(rna.sum(1) <= 0) or np.any(h3.sum(1) <= 0):
        raise ValueError("Empty molecular library")
    x, raw_y, descriptors = logcpm(rna), logcpm(h3)[:, selected], concentration(h3)
    predictions = {arm: np.zeros_like(raw_y) for arm in (*ARMS, "best_unweighted")}
    truth = np.zeros_like(raw_y)
    choices, candidate_records, all_pairs = [], [], []
    for outer in range(5):
        train, test = np.flatnonzero(fold != outer), np.flatnonzero(fold == outer)
        candidates = {arm: [] for arm in ARMS}
        for arm, rank in ARMS.items():
            strengths = design["contrast_strengths"] if arm.startswith("weighted_") else [0.0]
            for strength in strengths:
                for penalty in design["penalty_fractions"]:
                    pair_sse, profile_sse = 0.0, 0.0
                    for inner in range(5):
                        if inner == outer:
                            continue
                        internal = np.flatnonzero((fold != outer) & (fold != inner))
                        valid = np.flatnonzero((fold != outer) & (fold == inner))
                        pairs = paired_indices(valid, strata, ids)
                        if not len(pairs):
                            raise ValueError("Empty inner contrast evaluation")
                        y = residual_target(raw_y, descriptors, internal)
                        w = donor_weights([strata[i] for i in internal], strength)
                        pred = fit_predict(x[internal], y[internal], x[valid], w, penalty, rank)
                        pair_loss, profile_loss = losses(pred, y, valid, pairs)
                        pair_sse += pair_loss
                        profile_sse += profile_loss
                    case = {"outer_fold": outer, "arm": arm, "strength": strength, "penalty": penalty,
                            "pair_sse": pair_sse, "profile_sse": profile_sse}
                    candidates[arm].append(case)
                    candidate_records.append(case)
        unweighted = [choose(candidates[arm]) for arm in ("unweighted_ridge", "unweighted_rank12")]
        baseline = min(unweighted, key=lambda case: (case["pair_sse"], case["arm"]))
        profile_reference = min(case["profile_sse"] for arm in ("unweighted_ridge", "unweighted_rank12") for case in candidates[arm])
        (out / f"inner_candidates_outer_{outer}.json").write_text(json.dumps(candidates, indent=2, allow_nan=False) + "\n")
        kept = {case["arm"]: case for case in unweighted}
        for arm in ("weighted_ridge", "weighted_rank12"):
            try:
                kept[arm] = choose(candidates[arm], 1.01 * profile_reference)
            except ValueError:
                refusal = {"status": "prespecified_profile_safeguard_refusal_incomplete_evaluation",
                           "outer_fold": outer, "refused_arm": arm,
                           "minimum_arm_profile_sse": min(case["profile_sse"] for case in candidates[arm]),
                           "best_unweighted_profile_sse": profile_reference,
                           "profile_sse_limit": 1.01 * profile_reference,
                           "complete_outer_folds": outer, "checks": "No threshold, arm, split or candidate changed",
                           "receipt": receipt}
                (out / "summary.json").write_text(json.dumps(refusal, indent=2, allow_nan=False) + "\n")
                print(json.dumps(refusal, allow_nan=False), flush=True)
                return
        y = residual_target(raw_y, descriptors, train)
        truth[test] = y[test]
        for arm, case in kept.items():
            w = donor_weights([strata[i] for i in train], case["strength"])
            predictions[arm][test] = fit_predict(x[train], y[train], x[test], w, case["penalty"], ARMS[arm])
            choices.append(case)
        predictions["best_unweighted"][test] = predictions[baseline["arm"]][test]
        all_pairs.extend(paired_indices(test, strata, ids).tolist())
        (out / f"outer_{outer}.json").write_text(json.dumps({"choices": kept, "baseline": baseline,
                                                           "inner_profile_reference": profile_reference}, indent=2) + "\n")
        print(f"outer_{outer} saved", flush=True)
    pairs = np.asarray(all_pairs, int)
    if len(pairs) != 22 or len(set(pairs.ravel())) != 44:
        raise ValueError("Existing disjoint 22-pair census differs")
    pd.DataFrame({"participant_a": ids[pairs[:, 0]], "participant_b": ids[pairs[:, 1]]}).to_csv(out / "fixed_pairs.tsv", sep="\t", index=False)
    difference = truth[pairs[:, 0]] - truth[pairs[:, 1]]
    reference = np.sum(difference ** 2, axis=1)
    if not np.isfinite(reference).all() or np.any(reference <= 0):
        raise ValueError("Nonpositive or nonfinite contrast reference denominator")
    rng = np.random.default_rng(design["seed"])
    draws = rng.integers(0, 22, (4000, 22))
    pair_losses = {}
    for arm, pred in predictions.items():
        if not np.isfinite(pred).all():
            raise ValueError("Nonfinite held prediction")
        residual = (pred[pairs[:, 0]] - pred[pairs[:, 1]]) - difference
        pair_losses[arm] = np.sum(residual ** 2, axis=1)
    comparisons = {}
    for arm in ("weighted_ridge", "weighted_rank12"):
        delta = pair_losses["best_unweighted"] - pair_losses[arm]
        point = float(delta.sum() / reference.sum())
        samples = delta[draws].sum(1) / reference[draws].sum(1)
        comparisons[arm] = {"delta_pair_skill": point, "conditional_ci95": np.quantile(samples, [0.025, 0.975]).tolist(),
                            "planning_point_margin_met": point >= 0.01}
    result = {"status": "inspected_source_development_no_independent_confirmation", "participants": 99, "pairs": 22,
              "regions": len(selected), "comparisons": comparisons,
              "pair_skill": {arm: float(1 - loss.sum() / reference.sum()) for arm, loss in pair_losses.items()},
              "profile_sse": {arm: float(np.sum((pred - truth) ** 2)) for arm, pred in predictions.items()},
              "choices": choices, "numpy": np.__version__, "seed": design["seed"],
              "design_sha256": hashlib.sha256(DESIGN.read_bytes()).hexdigest(),
              "code_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                              for name in ("chromatin_contrast_development.py", "weighted_chromatin.py")}}
    np.savez_compressed(out / "development_predictions.npz", truth=truth, participant_ids=ids, pairs=pairs, **predictions)
    pd.DataFrame(candidate_records).to_csv(out / "inner_candidate_losses.tsv", sep="\t", index=False)
    (out / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(out), "comparisons": comparisons}), flush=True)


if __name__ == "__main__":
    main()
