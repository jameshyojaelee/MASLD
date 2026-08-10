#!/usr/bin/env python3
"""Validate and freeze the terminal Plan 44 no-promotion release."""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
CANDIDATE = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates/"
    "multicellular-assembly-response-2026-08-09"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    temporary.replace(path)


def atomic_text(path: Path, value: str) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value)
    temporary.replace(path)


def truth(value: str) -> bool:
    return str(value).strip().lower() == "true"


def main() -> None:
    checks: list[dict[str, object]] = []

    def add(check_id: str, passed: bool, detail: object) -> None:
        checks.append(
            {
                "check_id": check_id,
                "passed": str(bool(passed)).lower(),
                "detail": str(detail),
            }
        )

    preflight = json.loads((CANDIDATE / "PREFLIGHT_VALIDATED.json").read_text())
    add(
        "preflight_manifest_unchanged",
        sha256(CANDIDATE / "release_manifest.tsv")
        == preflight["release_manifest_sha256"]
        == "e4927dff7deaef91373ecd68237ec9434a848af7965f2b1557d18eb6edce27e9",
        preflight["release_manifest_sha256"],
    )

    upstream_ok = True
    for row in read_tsv(CANDIDATE / "upstream_immutability_baseline.tsv"):
        path = ROOT / row["path"]
        upstream_ok &= (
            path.is_file()
            and path.stat().st_size == int(row["bytes"])
            and sha256(path) == row["sha256"]
        )
    add("predecessor_releases_byte_identical", upstream_ok, "Plan 20/43 frozen inputs")

    source_ok = True
    source_rows = read_tsv(CANDIDATE / "terminal_source_manifest.tsv")
    for row in source_rows:
        path = CANDIDATE / row["relative_path"]
        source_ok &= (
            path.is_file()
            and path.stat().st_size == int(row["bytes"])
            and sha256(path) == row["sha256"]
        )
    add("public_source_checksums", source_ok and len(source_rows) == 7, len(source_rows))

    geometry = read_tsv(CANDIDATE / "frozen_lineage_geometry.tsv")
    lineages = {}
    for row in geometry:
        lineages[row["cell_type"]] = lineages.get(row["cell_type"], 0) + 1
    add(
        "complete_frozen_117_program_geometry",
        len(geometry) == 117
        and lineages
        == {
            "cholangiocytes": 28,
            "fibroblasts": 29,
            "hepatocytes": 30,
            "macrophages": 19,
            "tcells": 11,
        },
        lineages,
    )

    statistics = read_tsv(CANDIDATE / "statistical_rederivation.tsv")
    add(
        "independent_statistical_rederivation",
        len(statistics) >= 13 and all(truth(row["passed"]) for row in statistics),
        f"checks={len(statistics)}",
    )

    mps_gate = read_tsv(CANDIDATE / "mps/mps_gate_status.tsv")[0]
    human_gate = read_tsv(CANDIDATE / "human_response/human_gate_status.tsv")[0]
    add(
        "mps_gate_fail_closed",
        mps_gate["status"] == "complete_nonconfirmatory"
        and not truth(mps_gate["all_primary_direction"])
        and float(mps_gate["conjunction_p"]) == 1.0,
        mps_gate,
    )
    add(
        "human_gate_fail_closed",
        human_gate["status"] == "complete_nonconfirmatory"
        and not truth(human_gate["primary_direction_met"])
        and int(human_gate["n_participants"]) == 19,
        human_gate,
    )

    verdict = {
        row["component_id"]: row
        for row in read_tsv(CANDIDATE / "promotion_verdict.tsv")
    }
    required_verdicts = {
        "controlled_multicellular_assembly",
        "paired_human_remodeling_disassembly",
        "supplementary_perturbseq_observability",
        "complete_lineage_assembly_disassembly",
        "public_data_escalation_program",
    }
    add("complete_verdict_family", required_verdicts == set(verdict), sorted(verdict))
    add(
        "promotion_mechanically_false",
        all(not truth(row["gate_pass"]) for row in verdict.values())
        and verdict["complete_lineage_assembly_disassembly"]["status"]
        == "complete_no_promotion"
        and verdict["public_data_escalation_program"]["status"] == "terminal_stop",
        "all required gates false",
    )

    perturbseq = read_tsv(CANDIDATE / "perturbseq_observability.tsv")[0]
    add(
        "nonblocking_perturbseq_stop_explicit",
        perturbseq["status"] == "not_run_nonblocking_after_primary_failure"
        and perturbseq["estimate"] == "",
        perturbseq["status"],
    )
    sensitivity = read_tsv(CANDIDATE / "sensitivity.tsv")
    stopped = {
        row["sensitivity"]
        for row in sensitivity
        if row["status"] == "not_run_after_primary_direction_falsification"
    }
    add(
        "nondecisive_sensitivity_stop_explicit",
        stopped
        == {
            "sample_level_voom_clustered_bootstrap",
            "leave_one_condition_family_out",
            "source_batch_adjustment",
            "family_overlap_reduction",
        },
        sorted(stopped),
    )

    mandatory = [
        "ADJUDICATED.json",
        "mps_sample_manifest.tsv",
        "mps_condition_design.tsv",
        "mps_lineage_effects.tsv",
        "human_pair_manifest.tsv",
        "human_response_topology.tsv",
        "perturbseq_observability.tsv",
        "sensitivity.tsv",
        "promotion_verdict.tsv",
        "terminal_source_gate_status.tsv",
        "terminal_source_manifest.tsv",
        "environment_manifest.tsv",
        "execution_manifest.tsv",
        "statistical_rederivation.tsv",
    ]
    for relative in mandatory:
        path = CANDIDATE / relative
        add(
            f"mandatory_output:{relative}",
            path.is_file() and path.stat().st_size > 0,
            path.stat().st_size if path.is_file() else 0,
        )

    failures = [row for row in checks if row["passed"] != "true"]
    write_tsv(CANDIDATE / "validation_report.tsv", checks, ["check_id", "passed", "detail"])
    if failures:
        raise RuntimeError(json.dumps(failures, indent=2, sort_keys=True))

    excluded_parts = {"logs", "archive", "source_expression", "__pycache__"}
    excluded_names = {
        "release_manifest_terminal.tsv",
        "VALIDATED.json",
        "COMPLETE",
    }
    payload = []
    for path in sorted(CANDIDATE.rglob("*")):
        relative = path.relative_to(CANDIDATE)
        if (
            not path.is_file()
            or excluded_parts.intersection(relative.parts)
            or path.name in excluded_names
        ):
            continue
        payload.append(
            {
                "relative_path": str(relative),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    write_tsv(
        CANDIDATE / "release_manifest_terminal.tsv",
        payload,
        ["relative_path", "bytes", "sha256"],
    )
    validation = {
        "release_id": "multicellular-assembly-response-2026-08-09",
        "status": "complete_no_promotion",
        "validated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "checks_passed": len(checks),
        "checks_failed": 0,
        "release_payloads": len(payload),
        "release_manifest_sha256": sha256(
            CANDIDATE / "release_manifest_terminal.tsv"
        ),
        "figure5_promotion": False,
        "canonical_promotion_authorized": False,
        "public_data_escalation_ended": True,
    }
    atomic_text(
        CANDIDATE / "VALIDATED.json",
        json.dumps(validation, indent=2, sort_keys=True) + "\n",
    )
    atomic_text(CANDIDATE / "COMPLETE", validation["release_manifest_sha256"] + "\n")
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
