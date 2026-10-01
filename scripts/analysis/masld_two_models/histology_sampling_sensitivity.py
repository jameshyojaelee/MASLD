#!/usr/bin/env python3
"""Proportional histology-stratified random acquisition on the original splits.

The original equal-count stratum balancing is retained as a separately named
policy. This sensitivity selects donors using recorded histology and fixed
random seeds only, then reuses the original target transform and refit method.
"""
import argparse
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import chromatin_acquisition as original

HERE = Path(__file__).resolve().parent
ADDENDUM = HERE / "histology_sampling_addendum.json"
REFERENCE = (original.ROOT / "GWAS/finemapping/results/alphagenome_campaign/"
             "two-models-20260922T203900EDT/chromatin_acquisition_21849197")
POLICY = "histology_proportional_random"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def proportional_selections(hist, pool, budgets, seed, outer):
    """Largest-remainder proportional quotas with randomized donor selection.

    Returns positions within pool. Neither RNA nor H3K27ac is an input.
    Sorting strata and using one permutation per stratum make random draws
    reproducible. Integer proportional quotas can be nonnested over budgets.
    """
    if not np.isfinite(hist[pool]).all():
        raise ValueError("Missing histology requires a prespecified stratum rule")
    _, stratum = np.unique(hist[pool], axis=0, return_inverse=True)
    counts = np.bincount(stratum)
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), int(outer)]))
    tie_priority = rng.random(len(counts))
    permutations = [rng.permutation(np.flatnonzero(stratum == j))
                    for j in range(len(counts))]
    output = {}
    for k in budgets:
        if not 1 <= k <= len(pool):
            raise ValueError("Acquisition budget is outside selectable pool")
        expected = k * counts.astype(float) / len(pool)
        quotas = np.floor(expected).astype(int)
        remaining = k - int(quotas.sum())
        order = np.lexsort((tie_priority, -(expected - quotas)))
        quotas[order[:remaining]] += 1
        chosen = np.concatenate([positions[:quota]
                                 for positions, quota in zip(permutations, quotas)])
        if len(chosen) != k or len(np.unique(chosen)) != k:
            raise AssertionError("Stratified selection is not without replacement")
        if np.any(quotas > counts) or not np.all(np.abs(quotas - expected) < 1):
            raise AssertionError("Integer quotas are not proportional")
        output[int(k)] = {"positions": chosen, "stratum": stratum,
                          "stratum_counts": counts, "stratum_quotas": quotas}
    return output


def split(ids, folds, outer):
    evaluation = np.flatnonzero(folds == outer)
    available = np.flatnonzero(folds != outer)
    available = np.array(sorted(available, key=lambda j: hashlib.sha256(
        (str(original.SEED) + "|acquisition|" + ids[j]).encode()).digest()))
    base, pool = available[:45], available[45:]
    if set(base) & set(pool) or set(available) & set(evaluation):
        raise AssertionError("Acquisition donor split overlaps")
    return base, pool, evaluation


def self_check():
    hist = np.array([[0, 0, 0, 0]] * 7 + [[1, 0, 0, 1]] * 5
                    + [[2, 1, 1, 2]] * 3, dtype=float)
    pool = np.arange(len(hist))
    budgets = [3, 8, 13]
    a = proportional_selections(hist, pool, budgets, 2026092201, 2)
    b = proportional_selections(hist, pool, budgets, 2026092201, 2)
    c = proportional_selections(hist, pool, budgets, 2026092202, 2)
    for k in budgets:
        assert np.array_equal(a[k]["positions"], b[k]["positions"])
        assert len(a[k]["positions"]) == k
        assert np.all(np.bincount(a[k]["stratum"][a[k]["positions"]], minlength=3)
                      == a[k]["stratum_quotas"])
    assert any(not np.array_equal(a[k]["positions"], c[k]["positions"])
               for k in budgets)
    # The helper cannot read outcomes; changing an unrelated outcome matrix
    # leaves the same fixed-seed histology-only selection unchanged.
    ids = np.array([f"person{i:02d}" for i in range(99)])
    folds = np.arange(99) % 5
    for outer in range(5):
        base, candidates, evaluation = split(ids, folds, outer)
        assert len(base) == 45 and len(set(base) | set(candidates) | set(evaluation)) == 99
    print("PASS: proportional quotas, no replacement, fixed-seed reproducibility, donor separation")


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    if args.self_check:
        self_check()
        return
    if args.out is None:
        raise ValueError("--out is required for acquisition analysis")
    if args.out.exists():
        raise FileExistsError(args.out)
    spec = json.loads(ADDENDUM.read_text())
    if original.SEED != spec["original_split_seed"] or original.N_REGIONS != 1000:
        raise ValueError("Original acquisition design differs")
    reference = args.reference
    reference_fold = pd.read_csv(reference / "acquisition_per_fold.tsv", sep="\t")
    reference_people = pd.read_csv(reference / "acquisition_per_participant.tsv.gz", sep="\t")
    reference_summary = pd.read_csv(reference / "acquisition_summary.tsv", sep="\t")
    axis = pd.read_csv(original.FIX / "molecular/participant_axis.tsv", sep="\t")
    ids = axis.participant_id.astype(str).to_numpy()
    if len(ids) != 99 or len(set(ids)) != 99:
        raise ValueError("Expected 99 distinct participants")
    folds = pd.read_csv(original.FIX / "folds/participant_outer_folds.tsv", sep="\t").set_index(
        "participant_id").loc[ids, "outer_fold"].to_numpy(int)
    join = pd.read_csv(original.JOIN, sep="\t").set_index("participant_id").loc[ids]
    hist = join[spec["histology_fields"]].to_numpy(float)
    with np.load(original.WEIGHTS, allow_pickle=True) as weights:
        selected, keys, cis = original.panel(weights)
    reference_panel = pd.read_csv(reference / "fixed_region_panel.tsv", sep="\t")
    if not np.array_equal(selected, reference_panel.region_index.to_numpy()):
        raise ValueError("Annotation-fixed region indices differ from reference")
    if not np.array_equal(keys, reference_panel.region_key.to_numpy()):
        raise ValueError("Annotation-fixed region identities differ from reference")
    args.out.mkdir(parents=True)
    rna = np.asarray(np.load(original.FIX / "molecular/rna_values.npy", mmap_mode="r"), float)
    h3 = np.asarray(np.load(original.FIX / "molecular/h3k27ac_counts.npy", mmap_mode="r"), float)
    x, yraw = original.logcpm(rna), original.logcpm(h3)
    concentration = original.concentration(h3)
    rows, people, selections, split_rows = [], [], [], []
    for outer in sorted(set(folds)):
        base, pool, evaluation = split(ids, folds, outer)
        expected = reference_fold.loc[reference_fold.outer_fold.eq(outer),
            ["n_base", "n_pool", "n_evaluation"]].drop_duplicates().to_numpy()
        if expected.shape != (1, 3) or not np.array_equal(expected[0],
                [len(base), len(pool), len(evaluation)]):
            raise ValueError("Original fold census differs")
        for role, donors in (("base", base), ("pool", pool), ("evaluation", evaluation)):
            split_rows.extend({"outer_fold": int(outer), "participant_id": ids[j], "role": role}
                              for j in donors)
        Zbase = np.column_stack([np.ones(len(base)), concentration[base]])
        coef = np.linalg.lstsq(Zbase, yraw[base], rcond=None)[0]
        y = yraw - np.column_stack([np.ones(len(yraw)), concentration]) @ coef
        _, base_all = original.fit_profile(x, y, base, np.r_[pool, evaluation], cis, selected)
        error_pool = np.sum((y[pool][:, selected] - base_all[:len(pool)]) ** 2, axis=1)
        truth = y[evaluation][:, selected]
        base_loss = np.sum((truth - base_all[len(pool):]) ** 2, axis=1)
        previous = reference_people.loc[reference_people.outer_fold.eq(outer),
            ["participant_id", "base_SSE"]].drop_duplicates().set_index("participant_id")
        if len(previous) != len(evaluation) or not np.allclose(
                base_loss, previous.loc[ids[evaluation], "base_SSE"], rtol=1e-8, atol=1e-8):
            raise ValueError("Original participant-level baseline loss was not reproduced")
        for repeat, seed in enumerate(spec["repeat_seeds"]):
            draws = proportional_selections(hist, pool, spec["n_additional_measurements"], seed, outer)
            for fraction, k in zip(spec["budget_fractions_of_all_99"], spec["n_additional_measurements"]):
                draw = draws[k]
                pos = draw["positions"]
                chosen = pool[pos]
                training = np.r_[base, chosen]
                _, refit = original.fit_profile(x, y, training, evaluation, cis, selected)
                loss = np.sum((truth - refit) ** 2, axis=1)
                rows.append({"outer_fold": int(outer), "policy": POLICY, "random_repeat": repeat,
                    "repeat_seed": seed, "budget_fraction": fraction, "n_acquired": k,
                    "n_base": len(base), "n_pool": len(pool), "n_evaluation": len(evaluation),
                    "captured_pool_error_fraction": float(error_pool[pos].sum() / error_pool.sum()),
                    "evaluation_base_MSE": float(base_loss.mean() / len(selected)),
                    "evaluation_refit_MSE": float(loss.mean() / len(selected)),
                    "evaluation_MSE_reduction": float((base_loss.mean() - loss.mean()) / len(selected))})
                people.extend({"outer_fold": int(outer), "participant_id": ids[donor], "policy": POLICY,
                    "random_repeat": repeat, "repeat_seed": seed, "budget_fraction": fraction,
                    "base_SSE": float(base_loss[j]), "refit_SSE": float(loss[j])}
                    for j, donor in enumerate(evaluation))
                chosen_positions = set(pos)
                selections.extend({"outer_fold": int(outer), "participant_id": ids[donor],
                    "random_repeat": repeat, "repeat_seed": seed, "budget_fraction": fraction,
                    "stratum": int(draw["stratum"][j]),
                    "stratum_histology": "|".join(f"{v:g}" for v in hist[donor]),
                    "stratum_pool_n": int(draw["stratum_counts"][draw["stratum"][j]]),
                    "stratum_selected_n": int(draw["stratum_quotas"][draw["stratum"][j]]),
                    "selected": j in chosen_positions} for j, donor in enumerate(pool))
            print(json.dumps({"completed_fold": int(outer), "repeat": repeat,
                              "seed": seed, "refits": 3}), flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(args.out / "acquisition_per_fold.tsv", sep="\t", index=False)
    pd.DataFrame(people).to_csv(args.out / "acquisition_per_participant.tsv.gz", sep="\t", index=False)
    pd.DataFrame(selections).to_csv(args.out / "selected_pool_participants.tsv", sep="\t", index=False)
    pd.DataFrame(split_rows).to_csv(args.out / "participant_splits.tsv", sep="\t", index=False)
    summary = table.groupby(["budget_fraction", "policy"]).agg(
        n_folds=("outer_fold", "nunique"), n_fit_comparisons=("outer_fold", "size"),
        mean_error_capture=("captured_pool_error_fraction", "mean"),
        mean_refit_MSE_reduction=("evaluation_MSE_reduction", "mean"),
        mean_base_MSE=("evaluation_base_MSE", "mean")).reset_index()
    summary.to_csv(args.out / "acquisition_summary.tsv", sep="\t", index=False)
    reference_summary["policy"] = reference_summary.policy.replace(
        {"histology": "histology_equal_count_balancing"})
    pd.concat([reference_summary, summary], ignore_index=True).to_csv(
        args.out / "acquisition_summary_with_reference.tsv", sep="\t", index=False)
    averaged = table.groupby(["outer_fold", "budget_fraction"]).mean(numeric_only=True)
    contrasts = []
    for name in ("random", "histology"):
        control = reference_fold.loc[reference_fold.policy.eq(name)].groupby(
            ["outer_fold", "budget_fraction"]).mean(numeric_only=True)
        for index in averaged.index:
            contrasts.append({"outer_fold": int(index[0]), "budget_fraction": index[1],
                "reference_policy": "histology_equal_count_balancing" if name == "histology" else name,
                "error_capture_difference": float(averaged.loc[index, "captured_pool_error_fraction"]
                    - control.loc[index, "captured_pool_error_fraction"]),
                "MSE_reduction_difference": float(averaged.loc[index, "evaluation_MSE_reduction"]
                    - control.loc[index, "evaluation_MSE_reduction"])})
    pd.DataFrame(contrasts).to_csv(args.out / "paired_fold_differences.tsv", sep="\t", index=False)
    report = {"status": spec["status"], "participants": 99, "regions": len(selected),
        "policy": POLICY, "refits": len(table), "specification": spec,
        "baseline_reconstruction": "all five folds reproduce participant-level original baseline SSE",
        "target_unit": "original base-training-residual H3K27ac log2 CPM; MSE in squared units",
        "uncertainty": "descriptive paired fold differences only; no new independent confirmation or superiority claim",
        "source_sha256": {str(p): sha256(p) for p in [ADDENDUM, Path(__file__),
            Path(original.__file__), reference / "fixed_region_panel.tsv",
            reference / "acquisition_per_fold.tsv", reference / "acquisition_per_participant.tsv.gz"]},
        "environment": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
            "platform": platform.platform(), "slurm_job_id": os.environ["SLURM_JOB_ID"]}}
    (args.out / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path)
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--self-check", action="store_true")
    main(parser.parse_args())
