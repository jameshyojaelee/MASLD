"""Synthetic algebra and frozen-weight metadata checks; never reads RNA values."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys


ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
MEMBERS = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
COVERAGE = REC / "program_measurement_admission_21995500/program_coverage.tsv"
HASHES = {
    MEMBERS: "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b",
    COVERAGE: "e571982d2d522d543492fcc5e3455305cf4ca3a14c4b03e012bf66e59fbf721b",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def ratio(counts, weights, observed):
    """Separate nonnegative count-weighted estimand; no pseudocount."""
    require(len(counts) == len(weights) == len(observed) > 0, "Shape mismatch.")
    require(all(math.isfinite(x) and x >= 0 for x in counts), "Invalid counts.")
    require(all(math.isfinite(x) and x >= 0 for x in weights), "Signed/invalid weights refused.")
    require(all(type(x) is bool for x in observed), "Mask must be Boolean.")
    numerator = math.fsum(c * w for c, w, o in zip(counts, weights, observed) if not o)
    denominator = math.fsum(c * w for c, w, o in zip(counts, weights, observed) if o)
    require(denominator > 0, "No positive observed weighted contribution: abstain.")
    result = numerator / denominator
    require(math.isfinite(result), "Ratio exceeds finite numerical range: abstain.")
    return {"missing": numerator, "observed": denominator, "ratio": result,
            "log1p_ratio": math.log1p(result), "total_over_observed": 1 + result}


def rejected(call):
    try:
        call()
    except ValueError:
        return True
    return False


def main():
    require(os.environ.get("SLURM_JOB_ID"), "Run synthetic checks on a compute node.")
    for path, expected in HASHES.items():
        require(sha(path) == expected, "Frozen metadata hash changed.")
    with MEMBERS.open(newline="") as handle:
        members = list(csv.DictReader(handle, delimiter="\t"))
    with COVERAGE.open(newline="") as handle:
        coverage = list(csv.DictReader(handle, delimiter="\t"))
    require(len(members) == 7093 and len(coverage) == 117, "Membership shape changed.")
    require(all(float(m["canonical_weight_text"]) >= 0 and
                float(m["source_weight"]) >= 0 and
                float(m["original_l1_weight"]) >= 0 for m in members),
            "Actual canonical/source/original L1 weights include a negative weight.")
    groups = {}
    for row in members:
        groups.setdefault(row["program_uid"], []).append(row)
    require(len(groups) == 117, "Program count changed.")
    for uid, rows in groups.items():
        original = [float(m["original_l1_weight"]) for m in rows]
        canonical = [float(m["canonical_weight_text"]) for m in rows]
        total = math.fsum(canonical)
        require(total > 0 and math.isclose(math.fsum(original), 1, abs_tol=1e-10), uid)
        require(all(math.isclose(w, c / total, rel_tol=1e-9, abs_tol=1e-12)
                    for w, c in zip(original, canonical)), "L1 mapping altered canonical weights.")
    complete = [r for r in coverage if int(r["complete_fixture_target"]) == 1]
    require(len(complete) == 8, "Do not relax complete-target admission.")
    structural = [r["program_uid"] for r in complete if int(r["admitted_panel_members"]) == 0]
    require(structural == ["hotspot_tcells_05f8b3c103f336ef"], "Structural ratio roster changed.")

    weights, mask, counts = [0.2, 0.3, 0.5], [True, True, False], [4.0, 2.0, 10.0]
    base = ratio(counts, weights, mask)
    require(math.isclose(base["ratio"], 25 / 7), "Hand-derived ratio wrong.")
    errors = []
    for scale in [1e-6, 0.1, 1, 100, 1e6]:
        changed = ratio([scale * x for x in counts], weights, mask)
        errors.append(abs(changed["ratio"] - base["ratio"]))
        require(math.isclose(changed["ratio"], base["ratio"], rel_tol=1e-12), "Common scale invariance failed.")
        require(math.isclose(changed["log1p_ratio"], base["log1p_ratio"], rel_tol=1e-12), "Log1p ratio failed.")
    weight_scaled = ratio(counts, [7 * w for w in weights], mask)
    require(math.isclose(weight_scaled["ratio"], base["ratio"], rel_tol=1e-12), "Common weight scale failed.")
    require(math.isclose(base["total_over_observed"] * base["observed"],
                         base["observed"] + base["missing"]), "Count-unit correction identity failed.")
    zero_missing = ratio([4, 2, 0], weights, mask)
    require(zero_missing["ratio"] == zero_missing["log1p_ratio"] == 0, "Zero missing case failed.")
    failures = {
        "zero_observed": rejected(lambda: ratio([0, 0, 10], weights, mask)),
        "no_observed_members": rejected(lambda: ratio(counts, weights, [False] * 3)),
        "negative_weight": rejected(lambda: ratio(counts, [0.2, -0.3, 0.5], mask)),
        "negative_count": rejected(lambda: ratio([-1, 2, 10], weights, mask)),
        "nonfinite_count": rejected(lambda: ratio([math.nan, 2, 10], weights, mask)),
        "shape_mismatch": rejected(lambda: ratio(counts, weights, [True])),
    }
    require(all(failures.values()), "An invalid case was accepted.")

    # Additive count offsets depend on depth; log1p(M/O) itself does not.
    offset_before = (base["missing"] + 1) / (base["observed"] + 1)
    offset_after = (100 * base["missing"] + 1) / (100 * base["observed"] + 1)
    require(not math.isclose(offset_before, offset_after), "Offset counterexample failed.")
    # Unobserved background changes logCPM ratios despite unchanged program counts.
    logcpm_ratios = [math.log2(1 + 1e6 * 8 / total) /
                     math.log2(1 + 1e6 * 2 / total) for total in [100, 1000]]
    require(not math.isclose(*logcpm_ratios), "LogCPM denominator counterexample failed.")
    # Gene-specific capture changes observed composition and cannot cancel globally.
    bias = [10, 1, 1]
    biased = ratio([x * b for x, b in zip(counts, bias)], weights, mask)
    require(not math.isclose(biased["ratio"], base["ratio"]), "Capture bias counterexample failed.")
    panel_original = [counts[0] / sum(counts[:2]), counts[1] / sum(counts[:2])]
    panel_biased = [bias[i] * counts[i] / sum(bias[j] * counts[j] for j in range(2)) for i in range(2)]
    require(panel_original != panel_biased, "Panel capture composition counterexample failed.")
    same_input_other_missing = ratio([4, 2, 100], weights, mask)
    require(same_input_other_missing["observed"] == base["observed"] and
            math.isclose(same_input_other_missing["ratio"], 250 / 7),
            "Hand-derived identical-input counterexample failed.")
    require(same_input_other_missing["ratio"] != base["ratio"],
            "Observed expression alone must not identify missing expression algebraically.")
    out = REC / ("program_ratio_invariants_" + os.environ["SLURM_JOB_ID"])
    out.mkdir(exist_ok=False)
    summary = {
        "scope": "synthetic_algebra_and_frozen_membership_metadata_only",
        "job_id": os.environ["SLURM_JOB_ID"], "script_sha256": sha(Path(__file__)),
        "python": sys.version, "platform": platform.platform(),
        "input_sha256": {str(p): sha(p) for p in HASHES},
        "canonical_weights_all_nonnegative": True, "membership_rows": len(members),
        "programs": len(groups), "complete_target_roster": len(complete),
        "structurally_defined_ratio_programs": len(complete) - len(structural),
        "structural_abstentions": structural, "synthetic_hand_check": base,
        "maximum_common_scale_ratio_difference": max(errors),
        "invalid_case_refusals": failures,
        "count_offset_counterexample": [offset_before, offset_after],
        "logcpm_ratio_background_counterexample": logcpm_ratios,
        "gene_capture_bias_counterexample": {"original": base["ratio"], "biased": biased["ratio"]},
        "same_observed_input_distinct_missing_ratio": [base["ratio"], same_input_other_missing["ratio"]],
        "checks_passed": True, "expression_arrays_read": False, "outcomes_read": False,
        "fits_performed": False,
        "limit": "Algebra establishes a separate count-weighted ratio, not atlas score recovery or cross-assay validation.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(out), "checks_passed": True}))


if __name__ == "__main__":
    main()
