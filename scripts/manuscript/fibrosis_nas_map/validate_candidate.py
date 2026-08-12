#!/usr/bin/env python3
"""Independently validate the fibrosis-NAS candidate workstream."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import subprocess
from pathlib import Path


ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
WS = ROOT / "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v8"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def rows(path: Path) -> list[dict[str, str]]:
    require(path.is_file(), f"Missing table: {path}")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_manifest(root: Path) -> None:
    manifest = rows(root / "artifact_manifest.tsv")
    for item in manifest:
        path = root / item["artifact"]
        require(path.is_file(), f"Manifest artifact missing: {path}")
        require(sha256(path) == item["sha256"], f"Artifact hash drift: {path}")
        require(path.stat().st_size == int(item["size_bytes"]), f"Artifact size drift: {path}")


def main() -> None:
    contract = WS / "contract/analysis_contract.json"
    discovery = WS / "discovery"
    holdout = WS / "holdout"
    figures = WS / "figures"
    for path in (contract, discovery / "DISCOVERY_READY.json", holdout / "HOLDOUT_READY.json"):
        require(path.is_file(), f"Required READY artifact missing: {path}")

    contract_data = json.loads(contract.read_text())
    require(contract_data["features"]["outcome_blind_testable_programs"] == 113,
            "Contracted testable-program count drift")
    for item in contract_data["input_manifest"].values():
        path = ROOT / item["path"]
        require(path.is_file(), f"Contract input missing: {path}")
        require(sha256(path) == item["sha256"], f"Contract input hash drift: {path}")
    for relative, expected in contract_data["code_manifest"].items():
        path = ROOT / relative
        require(path.is_file(), f"Contract code missing: {path}")
        require(sha256(path) == expected, f"Contract code hash drift: {path}")

    d_ready = json.loads((discovery / "DISCOVERY_READY.json").read_text())
    h_ready = json.loads((holdout / "HOLDOUT_READY.json").read_text())
    require(d_ready["state"] == "DISCOVERY_COMPLETE_HOLDOUT_UNOPENED", "Bad discovery state")
    require(h_ready["state"] == "HOLDOUT_COMPLETE_NO_REFIT", "Bad holdout state")
    require(d_ready["contract_sha256"] == sha256(contract), "Discovery contract hash mismatch")
    require(h_ready["contract_sha256"] == sha256(contract), "Holdout contract hash mismatch")
    require(d_ready["n_discovery"] == 469, "Discovery n drift")
    require(d_ready["n_programs"] == 117, "Program registry n drift")
    require(d_ready["n_testable_programs"] == 113, "Testable program n drift")
    require(d_ready["holdout_accessed"] is False, "Discovery claims holdout access")
    require(h_ready["discovery_refit"] is False, "Holdout stage claims a discovery refit")
    require((h_ready["n_pairs"], h_ready["n_changed_pairs"], h_ready["n_unchanged_pairs"], h_ready["n_discordant_pairs"]) == (54, 46, 8, 11), "Paired census drift")

    sample_rows = rows(discovery / "sample_manifest.tsv")
    require(all(item["dataset"] != "GSE193066" for item in sample_rows), "HOLDOUT LEAK in discovery manifest")
    require(sum(item["role"] == "joint_fibrosis_nas_discovery" for item in sample_rows) == 469, "Discovery role count drift")

    program_map = rows(discovery / "program_map.tsv")
    require(len(program_map) == 117, "Program map must retain 117 programs")
    require(sum(item["testable"] == "TRUE" for item in program_map) == 113, "Program-map testability drift")
    require(sum(item["robust_display"] == "TRUE" for item in program_map) == 2, "Robust-display count drift")
    require(all(item["evidence_state"] == "untestable" for item in program_map if item["testable"] == "FALSE"), "Untestable program mislabeled")

    feature_scores = rows(discovery / "feature_scores.tsv.gz")
    required_score_fields = {
        "sample", "participant", "cohort", "biopsy", "view", "feature_id",
        "score_definition", "score", "testable",
    }
    require(required_score_fields.issubset(feature_scores[0]), "Feature-score schema is incomplete")
    require(all(item["cohort"] != "GSE193066" for item in feature_scores),
            "HOLDOUT LEAK in discovery feature scores")
    composition_features = {
        item["feature_id"] for item in feature_scores if item["view"] == "composition"
    }
    require(len(composition_features) == 22, "Composition score registry must retain 22 lineages")

    meta = rows(discovery / "meta_effects.tsv")
    primary = [item for item in meta if item["view"] == "hotspot_program" and item["score_definition"] == "weighted_mean_z" and item["model_kind"] == "primary" and item["estimable"] == "TRUE"]
    require(len(primary) == 226, f"Primary multiplicity family must contain 226 tests, found {len(primary)}")
    require(all(item["q_value"] not in {"", "NA"} for item in primary), "Primary family has missing BH q-values")
    require(all(item["family_size"] == "226" for item in primary), "Primary BH denominator drift")

    effect = rows(discovery / "effect_map.tsv")
    primary_effect = [
        item for item in effect
        if item["view"] == "hotspot_program"
        and item["score_definition"] == "weighted_mean_z"
        and item["model_kind"] == "primary"
    ]
    require(len(primary_effect) == 1170,
            f"Complete cohort/meta primary effect map must have 1170 rows, found {len(primary_effect)}")
    require(sum(item["effect_level"] == "meta" and item["estimable"] == "TRUE"
                for item in primary_effect) == 226, "Primary meta-effect count drift")
    require(sum(item["effect_level"] == "cohort" and item["estimable"] == "TRUE"
                for item in primary_effect) == 904, "Primary cohort-effect count drift")

    linearity = rows(discovery / "linearity_contrasts.tsv")
    require(linearity, "Categorical adjacent-contrast audit is empty")
    require(all(item["multiplicity_family"] ==
                "categorical_adjacent_contrasts_all_programs_axes_cohorts"
                for item in linearity), "Linearity multiplicity-family drift")
    require(all(int(item["family_size"]) == len(linearity) for item in linearity),
            "Linearity BH denominator drift")
    require(all(item["q_value"] not in {"", "NA"} for item in linearity
                if item["estimable"] == "TRUE"),
            "Estimable linearity contrast lacks a BH q-value")

    loco = rows(discovery / "validation_summary.tsv")
    require(sum(item["split"] == "leave_one_cohort_out" for item in loco) == 8, "LOCO fold count drift")
    require(sum(item["split"] == "leave_one_cohort_out_global" for item in loco) == 2, "LOCO global count drift")

    paired = rows(holdout / "paired_validation_summary.tsv")
    primary_pair = [item for item in paired if item["view"] == "hotspot_program" and item["score_definition"] == "weighted_mean_z"]
    require(len(primary_pair) == 1, "Missing unique primary paired result")
    require(int(primary_pair[0]["n_participants"]) == 46, "Primary paired test must use 46 changed participants")
    paired_programs = rows(holdout / "paired_program_validation.tsv")
    require(len(paired_programs) == 113, "Paired program family must contain every testable program")

    transport = rows(holdout / "transport_validation.tsv")
    require(len(transport) == 20, f"Expected 20 view-by-native-contrast transports, found {len(transport)}")
    program_transport = {
        (item["cohort"], item["axis"]): int(item["n_participants"])
        for item in transport
        if item["view"] == "hotspot_program"
    }
    require(program_transport == {
        ("GSE240729", "fibrosis_marginal"): 64,
        ("GSE126848", "disease_control"): 53,
        ("GSE126848", "nash_nafl"): 31,
        ("GSE167523", "nash_nafl"): 96,
        ("GSE213621", "disease_control"): 361,
    }, f"Native transport biological n drift: {program_transport}")

    combined_validation = rows(holdout / "validation_summary.tsv")
    required_validation_fields = {
        "split", "cohort", "view", "axis", "metric", "estimate", "ci_lower",
        "ci_upper", "empirical_p", "multiplicity_adjustment", "biological_n",
    }
    require(required_validation_fields.issubset(combined_validation[0]),
            "Combined validation-summary schema is incomplete")

    validate_manifest(discovery)
    validate_manifest(holdout)
    validate_manifest(figures)

    pdfs = sorted(figures.glob("*.pdf"))
    require(len(pdfs) == 11, f"Expected 11 candidate PDFs, found {len(pdfs)}")
    if subprocess.call(["which", "pdfinfo"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0:
        for pdf in pdfs:
            output = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
            pages = [line for line in output.splitlines() if line.startswith("Pages:")]
            require(pages and pages[0].split(":", 1)[1].strip() == "1", f"PDF is not one page: {pdf}")

    report = {
        "state": "VALIDATED_CANDIDATE_NOT_PROMOTED",
        "contract_sha256": sha256(contract),
        "discovery_ready_sha256": sha256(discovery / "DISCOVERY_READY.json"),
        "holdout_ready_sha256": sha256(holdout / "HOLDOUT_READY.json"),
        "programs": 117,
        "testable_programs": 113,
        "primary_tests": 226,
        "paired_participants": 54,
        "candidate_pdfs": len(pdfs),
        "composition_status": rows(discovery / "composition_status.tsv")[0],
        "promotion_decision": rows(figures / "promotion_decision.tsv"),
    }
    destination = WS / "VALIDATION_REPORT.json"
    if destination.exists():
        raise SystemExit(f"Refusing to overwrite validation report: {destination}")
    destination.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    os.chmod(destination, 0o440)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
