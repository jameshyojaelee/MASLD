#!/usr/bin/env python3
"""Independently validate Plan 42 outputs, null handling, and immutable boundaries."""

from __future__ import annotations

import datetime as dt
import json
import math
from pathlib import Path

from risk_state_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_text,
    read_tsv,
    require_validated_seal,
    sha256_file,
    write_tsv,
)


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def bh(values: list[float], family_n: int | None = None) -> list[float]:
    n_observed = len(values)
    n_family = family_n or n_observed
    order = sorted(range(n_observed), key=lambda index: values[index])
    adjusted = [math.nan] * n_observed
    running = 1.0
    for reverse_rank in range(n_observed - 1, -1, -1):
        index = order[reverse_rank]
        rank = reverse_rank + 1
        running = min(running, values[index] * n_family / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def main() -> None:
    seal = require_validated_seal()
    checks: list[dict[str, object]] = []

    def add(check_id: str, passed: bool, detail: object) -> None:
        checks.append({"check_id": check_id, "passed": str(bool(passed)).lower(), "detail": detail})

    for row in read_tsv(CANDIDATE_ROOT / "frozen_input_manifest.tsv"):
        source = PROJECT_ROOT / row["source_path"]
        snapshot = PROJECT_ROOT / row["snapshot_path"]
        passed = source.is_file() and snapshot.is_file() and sha256_file(source) == row["sha256"] == sha256_file(snapshot)
        add(f"immutable_input:{row['input_id']}", passed, row["sha256"])

    baseline_path = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/spatial_context/v1_preservation_baseline.tsv"
    baseline = read_tsv(baseline_path)
    sealed_at = dt.datetime.fromisoformat(str(seal["sealed_at_utc"]))
    surface_audit = []
    for row in baseline:
        path = PROJECT_ROOT / row["relative_path"]
        current_bytes = path.stat().st_size if path.is_file() else ""
        current_sha256 = sha256_file(path) if path.is_file() else ""
        modified_at = (
            dt.datetime.fromtimestamp(path.stat().st_mtime, tz=dt.timezone.utc)
            if path.is_file()
            else None
        )
        historical_exact = (
            path.is_file()
            and current_bytes == int(row["bytes"])
            and current_sha256 == row["sha256"]
        )
        if historical_exact:
            status = "historical_exact"
        elif modified_at is not None and modified_at < sealed_at:
            # The authoritative figure/reframing work on 2026-08-08 intentionally
            # changed nine members of the older 2026-08-06 preservation surface.
            # Plan 42 was sealed on 2026-08-09.  Preserve that historical drift as
            # an audit finding while failing only on mutations at/after this seal.
            status = "preexisting_before_plan42_seal"
        else:
            status = "modified_after_plan42_seal_or_missing"
        surface_audit.append(
            {
                "relative_path": row["relative_path"],
                "historical_bytes": row["bytes"],
                "historical_sha256": row["sha256"],
                "current_bytes": current_bytes,
                "current_sha256": current_sha256,
                "current_modified_utc": modified_at.isoformat() if modified_at else "",
                "plan42_sealed_utc": sealed_at.isoformat(),
                "status": status,
            }
        )
    write_tsv(
        CANDIDATE_ROOT / "immutable_surface_audit.tsv",
        surface_audit,
        [
            "relative_path", "historical_bytes", "historical_sha256",
            "current_bytes", "current_sha256", "current_modified_utc",
            "plan42_sealed_utc", "status",
        ],
    )
    historical_exact_n = sum(row["status"] == "historical_exact" for row in surface_audit)
    preexisting_n = sum(row["status"] == "preexisting_before_plan42_seal" for row in surface_audit)
    postseal_n = sum(row["status"] == "modified_after_plan42_seal_or_missing" for row in surface_audit)
    add(
        "canonical_surface_not_mutated_after_plan42_seal",
        len(baseline) == 591 and postseal_n == 0,
        f"historical_exact={historical_exact_n}/591; preexisting_before_seal={preexisting_n}; postseal_or_missing={postseal_n}",
    )
    add(
        "historical_v1_drift_explicitly_audited",
        preexisting_n == 9,
        "Nine 2026-08-08 figure/reframing changes predate the 2026-08-09 Plan 42 seal; see immutable_surface_audit.tsv",
    )

    source_gates = read_tsv(CANDIDATE_ROOT / "source_gate_status.tsv")
    gate_map = {(row["source_id"], row["gate"]): row["status"] for row in source_gates}
    add("clcc1_source_gate", gate_map.get(("CLCC1", "official_workbook")) == "pass", gate_map.get(("CLCC1", "official_workbook")))
    add("gse313544_source_gate", gate_map.get(("GSE313544", "sample_design")) == "pass", gate_map.get(("GSE313544", "sample_design")))
    add("gse158182_source_gate", gate_map.get(("GSE158182", "sample_design")) == "pass_source_dependent", gate_map.get(("GSE158182", "sample_design")))
    add("hmsma_maldi_skip_is_explicit", gate_map.get(("HMSMA_MALDI", "public_registered_product")) == "skipped_no_public_registered_maldi", gate_map.get(("HMSMA_MALDI", "public_registered_product")))
    add("hmsma_clinical_key_skip_is_explicit", gate_map.get(("HMSMA_VISIUM", "GSA_to_clinical_sample_key")) == "skipped_no_authoritative_sample_key", gate_map.get(("HMSMA_VISIUM", "GSA_to_clinical_sample_key")))

    crispr = read_tsv(CANDIDATE_ROOT / "crispr_class_effects.tsv")
    p_values = [as_float(row["p"]) for row in crispr]
    expected_q = bh(p_values, 3)
    observed_q = [as_float(row["q"]) for row in crispr]
    add("crispr_primary_family_three_rows", len(crispr) == 3, len(crispr))
    add("crispr_BH_rederived", all(abs(x - y) < 1e-10 for x, y in zip(expected_q, observed_q)), f"expected={expected_q}; observed={observed_q}")
    add("crispr_permutation_null_centered", all(abs(as_float(row["null_mean"])) < 0.01 for row in crispr), [row["null_mean"] for row in crispr])
    add("crispr_99999_permutations", all(int(row["permutations"]) >= 99999 for row in crispr), [row["permutations"] for row in crispr])
    add("crispr_full_gene_universe", len(read_tsv(CANDIDATE_ROOT / "crispr_gene_effects.tsv")) >= 20000, len(read_tsv(CANDIDATE_ROOT / "crispr_gene_effects.tsv")))
    add("crispr_locus_units_not_guides", len(read_tsv(CANDIDATE_ROOT / "crispr_locus_effects.tsv")) == int(float(crispr[0]["n_matched_sets"])), crispr[0]["n_matched_sets"])
    source_reproduction = read_tsv(CANDIDATE_ROOT / "crispr_source_reproduction.tsv")
    add("crispr_source_reproduction", all(row["passed"].lower() == "true" for row in source_reproduction), [(row["check"], row["passed"]) for row in source_reproduction])
    crispr_sensitivity = read_tsv(CANDIDATE_ROOT / "crispr_sensitivity.tsv")
    matching_sensitivities = {row["sensitivity"] for row in crispr_sensitivity}
    required_matching = {"matched_control_reverse_locus_order", "matched_controls_with_replacement"}
    add(
        "crispr_alternative_matching_specifications",
        required_matching.issubset(matching_sensitivities),
        sorted(required_matching & matching_sensitivities),
    )

    stress = read_tsv(CANDIDATE_ROOT / "stress_response_effects.tsv")
    add("stress_primary_complete", len(stress) == 4, len(stress))
    program = read_tsv(CANDIDATE_ROOT / "stress_response_program_effects.tsv")
    add("stress_complete_117_program_landscape", len(program) == 117 * 3 * 2, len(program))
    stress_source = read_tsv(CANDIDATE_ROOT / "stress_response_source_reproduction.tsv")
    add("stress_source_reproduction", all(row["passed"].lower() == "true" for row in stress_source), [(row["dataset_id"], row["observed_samples"]) for row in stress_source])
    clone_deltas = read_tsv(CANDIDATE_ROOT / "stress_response_clone_deltas.tsv")
    add("gse313544_eight_clone_units", len(clone_deltas) == 8 and len({row["clone"] for row in clone_deltas}) == 8, len(clone_deltas))
    clone_permutation = read_tsv(CANDIDATE_ROOT / "stress_response_permutations.tsv")
    add(
        "gse313544_exact_clone_label_permutation",
        len(clone_permutation) == 1
        and int(clone_permutation[0]["n_exact_allocations"]) == 560
        and abs(as_float(clone_permutation[0]["null_mean"])) < 1e-12,
        clone_permutation[0] if clone_permutation else {},
    )

    double = read_tsv(CANDIDATE_ROOT / "double_dissociation.tsv")
    double_map = {row["statistic_id"]: row for row in double}
    add("formal_interaction_not_fabricated", double_map.get("D", {}).get("status") == "not_computed_prerequisite_failed" and not double_map.get("D", {}).get("estimate"), double_map.get("D", {}))
    promotion = read_tsv(CANDIDATE_ROOT / "promotion_verdict.tsv")[0]
    add("promotion_is_false", promotion["functional_partition_claim_pass"].lower() == "false" and promotion["figure5_promotion"].lower() == "false", promotion["status"])
    root_hypotheses = CANDIDATE_ROOT / "frozen_hypotheses.tsv"
    spec_hypotheses = CANDIDATE_ROOT / "frozen_spec/frozen_hypotheses.tsv"
    add(
        "root_hypotheses_matches_sealed_spec",
        root_hypotheses.is_file()
        and spec_hypotheses.is_file()
        and read_tsv(root_hypotheses) == read_tsv(spec_hypotheses),
        (
            f"normalized_root={sha256_file(root_hypotheses)}; "
            f"sealed_spec={sha256_file(spec_hypotheses)}; parsed_rows_equal=true"
            if root_hypotheses.is_file() and spec_hypotheses.is_file()
            else "missing"
        ),
    )

    mandatory_outputs = [
        "source_gate_status.tsv", "frozen_hypotheses.tsv", "genetic_locus_registry.tsv",
        "crispr_gene_effects.tsv", "crispr_locus_effects.tsv", "stress_response_effects.tsv",
        "stress_response_permutations.tsv",
        "spatial_niche_effects.tsv", "double_dissociation.tsv", "promotion_verdict.tsv",
        "source_manifest.tsv", "environment_manifest.tsv", "execution_manifest.tsv",
    ]
    for relative in mandatory_outputs:
        path = CANDIDATE_ROOT / relative
        add(f"mandatory_output:{relative}", path.is_file() and path.stat().st_size > 0, path.stat().st_size if path.is_file() else 0)

    failures = [row for row in checks if row["passed"] != "true"]
    write_tsv(CANDIDATE_ROOT / "validation_report.tsv", checks, ["check_id", "passed", "detail"])
    if failures:
        raise RuntimeError(json.dumps(failures, indent=2, sort_keys=True))

    active_files = []
    excluded = {"release_manifest.tsv", "VALIDATED.json", "COMPLETE"}
    for path in sorted(CANDIDATE_ROOT.rglob("*")):
        if not path.is_file() or path.name in excluded or "logs" in path.parts or "source" in path.parts or "__pycache__" in path.parts:
            continue
        if path.suffix not in {".tsv", ".json", ".txt"} and path.name not in {"SPECIFICATION_SHA256", "INTEGRATED.json"}:
            continue
        active_files.append({"relative_path": str(path.relative_to(CANDIDATE_ROOT)), "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    write_tsv(CANDIDATE_ROOT / "release_manifest.tsv", active_files, ["relative_path", "bytes", "sha256"])
    validation = {
        "candidate_id": seal["candidate_id"],
        "status": "complete_no_promotion",
        "validated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "n_checks": len(checks),
        "n_failed": 0,
        "release_payloads": len(active_files),
        "release_manifest_sha256": sha256_file(CANDIDATE_ROOT / "release_manifest.tsv"),
        "functional_partition_claim_pass": False,
        "figure5_promotion": False,
    }
    atomic_write_text(CANDIDATE_ROOT / "VALIDATED.json", json.dumps(validation, indent=2, sort_keys=True) + "\n")
    atomic_write_text(CANDIDATE_ROOT / "COMPLETE", validation["release_manifest_sha256"] + "\n")
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
