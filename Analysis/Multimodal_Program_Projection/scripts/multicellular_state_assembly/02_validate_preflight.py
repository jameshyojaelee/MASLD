#!/usr/bin/env python3
"""Independently validate the outcome-blind Plan 44 preflight seal."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
CANDIDATE = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates/"
    "multicellular-assembly-response-2026-08-09"
)
EXPECTED_LINEAGES = Counter(
    {
        "hepatocytes": 30,
        "fibroblasts": 29,
        "macrophages": 19,
        "cholangiocytes": 28,
        "tcells": 11,
    }
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def soft_samples(path: Path) -> list[dict[str, list[str]]]:
    samples: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                current = {"accession": [line.split(" = ", 1)[1]]}
                samples.append(current)
            elif current is not None and line.startswith("!") and " = " in line:
                key, value = line.split(" = ", 1)
                current.setdefault(key, []).append(value)
    return samples


def only(sample: dict[str, list[str]], key: str) -> str:
    values = sample.get(key, [])
    assert len(values) == 1, (key, values)
    return values[0]


def characteristic(sample: dict[str, list[str]], field: str) -> str:
    values = [
        value[len(field) + 2 :]
        for value in sample.get("!Sample_characteristics_ch1", [])
        if value.startswith(field + ": ")
    ]
    assert len(values) == 1, (field, values)
    return values[0]


def main() -> None:
    checks: list[str] = []

    release = rows(CANDIDATE / "release_manifest.tsv")
    assert len(release) == 11, len(release)
    for row in release:
        path = CANDIDATE / row["relative_path"]
        assert path.is_file(), path
        assert path.stat().st_size == int(row["bytes"]), path
        assert sha256(path) == row["sha256"], path
    checks.append("release_manifest_hashes")

    ready = json.loads((CANDIDATE / "PREFLIGHT_READY.json").read_text())
    assert ready["status"] == "PREFLIGHT_READY_OUTCOMES_UNOPENED"
    assert ready["new_expression_outcomes_accessed"] is False
    assert ready["release_manifest_sha256"] == sha256(
        CANDIDATE / "release_manifest.tsv"
    )
    checks.append("ready_attestation")

    geometry = rows(CANDIDATE / "frozen_lineage_geometry.tsv")
    assert len(geometry) == 117
    assert Counter(row["cell_type"] for row in geometry) == EXPECTED_LINEAGES
    for lineage in EXPECTED_LINEAGES:
        loading = sum(
            abs(float(row["within_lineage_loading"]))
            for row in geometry
            if row["cell_type"] == lineage
        )
        assert abs(loading - 1.0) <= 1e-12, (lineage, loading)
    assert all(
        row["external_outcomes_read_at_freeze"].lower() in {"false", "0", "no", ""}
        for row in geometry
    )
    checks.append("lineage_geometry")

    hypotheses = rows(CANDIDATE / "frozen_hypotheses.tsv")
    assert [row["hypothesis_id"] for row in hypotheses] == [
        "MPS_FAT_HEPATOCYTE_SPECIFICITY",
        "MPS_TGFB_REMODELING_SPECIFICITY",
        "MPS_NPC_REMODELING_SPECIFICITY",
        "HUMAN_REMODELING_DISASSEMBLY",
        "HEPARG_OBSERVABILITY_BOUNDARY",
    ]
    assert [row["expected_direction"] for row in hypotheses] == [
        "positive",
        "positive",
        "positive",
        "negative",
        "positive",
    ]
    checks.append("frozen_hypotheses")

    source_manifest = rows(CANDIDATE / "source_manifest.tsv")
    assert len(source_manifest) == 3
    assert all(
        row["content_class"] == "metadata_only_no_expression_outcome"
        for row in source_manifest
    )
    for row in source_manifest:
        path = CANDIDATE / row["relative_path"]
        assert path.stat().st_size == int(row["bytes"])
        assert sha256(path) == row["sha256"]
    checks.append("metadata_only_sources")

    gse168 = soft_samples(CANDIDATE / "source_metadata/GSE168285_family.soft.gz")
    assert len(gse168) == 179
    design: Counter[tuple[str, str]] = Counter()
    for sample in gse168:
        experiment = only(sample, "!Sample_title").split("_", 1)[0]
        treatment = characteristic(sample, "treatment")
        design[(experiment, treatment)] += 1
    assert len(design) == 60
    assert Counter(design.values()) == Counter({3: 59, 2: 1})
    assert {key[0] for key in design} == {
        "NAFT18012",
        "NAFT18013",
        "NAFT18014",
        "NAFT18015",
    }
    checks.append("gse168_design")

    gse175 = soft_samples(CANDIDATE / "source_metadata/GSE175448_family.soft.gz")
    assert len(gse175) == 38
    participants: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for sample in gse175:
        description = only(sample, "!Sample_description")
        participant, source_arm, *_ = description.split("_")
        arm = "cenicriviroc" if source_arm == "drugA" else "placebo"
        participants[participant].append(
            (
                arm,
                characteristic(sample, "timing of biopsy"),
                characteristic(sample, "fibrosis improvement"),
            )
        )
    assert len(participants) == 19
    arm_response: Counter[tuple[str, str]] = Counter()
    for participant_rows in participants.values():
        assert len(participant_rows) == 2
        assert {row[1] for row in participant_rows} == {"pre-treatment", "post-treatment"}
        post = next(row for row in participant_rows if row[1] == "post-treatment")
        arm_response[(post[0], post[2])] += 1
    assert arm_response == Counter(
        {
            ("cenicriviroc", "improved"): 4,
            ("cenicriviroc", "not_improved"): 5,
            ("placebo", "improved"): 3,
            ("placebo", "not_improved"): 7,
        }
    )
    checks.append("gse175_pairs")

    gse238 = soft_samples(CANDIDATE / "source_metadata/GSE238219_family.soft.gz")
    assert len(gse238) == 5
    targets = rows(CANDIDATE / "gse238219_target_registry.tsv")
    assert len(targets) == 11
    genetic_only = {
        row["canonical_symbol"]
        for row in targets
        if row["evidence_class"] == "genetic_only"
    }
    assert genetic_only == {"ILRUN", "GPAM", "PPP1R3B", "TNKS"}
    assert next(row for row in targets if row["source_symbol"] == "C6orf106")[
        "canonical_symbol"
    ] == "ILRUN"
    assert {
        row["canonical_symbol"]
        for row in targets
        if row["any_direct_tier1_support"] == "true"
    } == {"GPAM"}
    checks.append("gse238_alias_and_overlap")

    gates = rows(CANDIDATE / "source_gate_status.tsv")
    assert len(gates) == 3
    assert all(row["outcomes_accessed"] == "false" for row in gates)
    assert {row["status"] for row in gates} == {
        "pass_model_specific_no_donor_key",
        "pass_authoritative_pairs",
        "pass_supplementary_only",
    }
    checks.append("source_gate_semantics")

    design_gate = rows(CANDIDATE / "mps_design_gate.tsv")
    assert len(design_gate) == 1
    design_row = design_gate[0]
    assert design_row["status"] == "pass_full_rank_estimable"
    assert int(design_row["n_source_samples"]) == 179
    assert int(design_row["n_condition_means"]) == 60
    assert int(design_row["n_long_rows"]) == 180
    assert int(design_row["n_model_columns"]) == 20
    assert int(design_row["design_rank"]) == 20
    assert float(design_row["primary_vif_max"]) < 10
    assert design_row["preflight_manifest_sha256"] == sha256(
        CANDIDATE / "release_manifest.tsv"
    )
    assert design_row["new_expression_outcomes_accessed"] == "false"
    checks.append("mps_design_gate")

    upstream = rows(CANDIDATE / "upstream_immutability_baseline.tsv")
    assert len(upstream) >= 4
    for row in upstream:
        path = ROOT / row["path"]
        assert path.is_file(), path
        assert path.stat().st_size == int(row["bytes"]), path
        assert sha256(path) == row["sha256"], path
    checks.append("upstream_byte_identity")

    allowed = {
        "PREFLIGHT_READY.json",
        "PREFLIGHT_VALIDATED.json",
        "frozen_hypotheses.tsv",
        "frozen_lineage_geometry.tsv",
        "gse238219_target_registry.tsv",
        "mps_design_gate.tsv",
        "preflight_validation.tsv",
        "release_manifest.tsv",
        "source_gate_status.tsv",
        "source_known_findings.tsv",
        "source_manifest.tsv",
        "upstream_immutability_baseline.tsv",
        "source_metadata/GSE168285_family.soft.gz",
        "source_metadata/GSE175448_family.soft.gz",
        "source_metadata/GSE238219_family.soft.gz",
    }
    observed = {
        path.relative_to(CANDIDATE).as_posix()
        for path in CANDIDATE.rglob("*")
        if path.is_file()
    }
    assert observed <= allowed, sorted(observed - allowed)
    checks.append("no_unsealed_payloads")

    validation = {
        "status": "PREFLIGHT_VALIDATED_OUTCOMES_UNOPENED",
        "checks_passed": len(checks),
        "checks": checks,
        "release_manifest_sha256": sha256(CANDIDATE / "release_manifest.tsv"),
        "validator_sha256": sha256(Path(__file__)),
        "new_expression_outcomes_accessed": False,
    }
    (CANDIDATE / "PREFLIGHT_VALIDATED.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "PLAN44_PREFLIGHT_VALIDATED\t"
        f"checks={len(checks)}\tmanifest={validation['release_manifest_sha256']}"
    )


if __name__ == "__main__":
    main()
