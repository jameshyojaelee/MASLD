#!/usr/bin/env python3
"""Seal program, evidence-class, source, contrast, and promotion rules outcome-blindly."""

from __future__ import annotations

import csv
import datetime as dt
import json
import shutil
import subprocess
from pathlib import Path

from public_functional_common import (
    CANDIDATE_ID,
    CANDIDATE_ROOT,
    CONFIG_ROOT,
    PROJECT_ROOT,
    SCRIPT_ROOT,
    atomic_write_text,
    read_tsv,
    require_within,
    sha256_file,
    stable_json_sha256,
    write_tsv,
)


ALLOWED_CLASSES = {
    "disease_state_only",
    "genetic_only",
    "convergent",
    "neither",
    "indeterminate_not_jointly_testable",
}
PRIMARY_MODULES = {("hepatocytes", "8"), ("hepatocytes", "20")}


def main() -> None:
    require_within(CANDIDATE_ROOT).mkdir(parents=True, exist_ok=True)
    seal_path = CANDIDATE_ROOT / "SEALED.json"
    if seal_path.exists():
        raise RuntimeError(
            f"Candidate already sealed; validate rather than overwrite: {seal_path}"
        )

    frozen_inputs = read_tsv(CONFIG_ROOT / "frozen_inputs.tsv")
    manifest: list[dict[str, object]] = []
    for row in frozen_inputs:
        source = PROJECT_ROOT / row["relative_path"]
        if not source.is_file():
            raise FileNotFoundError(source)
        observed = sha256_file(source)
        if observed != row["expected_sha256"]:
            raise RuntimeError(
                f"Frozen input drift for {row['input_id']}: {observed} != {row['expected_sha256']}"
            )
        destination = CANDIDATE_ROOT / "frozen_inputs" / source.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        if sha256_file(destination) != observed:
            raise RuntimeError(f"Copy verification failed: {destination}")
        manifest.append(
            {
                "input_id": row["input_id"],
                "source_path": row["relative_path"],
                "snapshot_path": str(destination.relative_to(PROJECT_ROOT)),
                "sha256": observed,
                "size_bytes": source.stat().st_size,
                "role": row["role"],
            }
        )

    registry = read_tsv(CANDIDATE_ROOT / "frozen_inputs/program_registry_v2.tsv")
    memberships = read_tsv(CANDIDATE_ROOT / "frozen_inputs/program_membership_v2.tsv")
    external = read_tsv(CANDIDATE_ROOT / "frozen_inputs/external_test_programs.tsv")
    classes = read_tsv(CANDIDATE_ROOT / "frozen_inputs/frozen_evidence_classes.tsv")
    if len(registry) != 117:
        raise RuntimeError(f"Expected 117 registry rows, found {len(registry)}")
    if {(r["cell_type"], r["module"]) for r in external} != PRIMARY_MODULES:
        raise RuntimeError("Externally frozen program family is not hepatocyte modules 8 and 20")
    observed_classes = {r["static_class"] for r in classes}
    if not observed_classes.issubset(ALLOWED_CLASSES):
        raise RuntimeError(f"Unexpected evidence classes: {sorted(observed_classes - ALLOWED_CLASSES)}")

    primary_uids = {r["program_uid"] for r in external}
    primary_memberships = [r for r in memberships if r["program_uid"] in primary_uids]
    if len(primary_memberships) != 71:
        raise RuntimeError(f"Expected 71 frozen primary memberships, found {len(primary_memberships)}")
    write_tsv(
        CANDIDATE_ROOT / "primary_program_membership.tsv",
        primary_memberships,
        list(primary_memberships[0]),
    )

    config_manifest = []
    for config in sorted(CONFIG_ROOT.glob("*.tsv")):
        destination = CANDIDATE_ROOT / "frozen_spec" / config.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(config, destination)
        config_manifest.append(
            {
                "config": str(config.relative_to(PROJECT_ROOT)),
                "snapshot": str(destination.relative_to(PROJECT_ROOT)),
                "sha256": sha256_file(config),
                "size_bytes": config.stat().st_size,
            }
        )

    gate_rows = [
        {"gate_order": 1, "gate_id": "PCLS_INDUCTION", "required": "true", "rule": "q<0.05 in frozen two-program family; expected direction; >=5/6 donors when six pairs exist; >=4 pairs required"},
        {"gate_order": 2, "gate_id": "PCLS_ROBUSTNESS", "required": "true", "rule": "equal-weight, leave-top-gene, leave-one-donor-out, and culture-time controls pass"},
        {"gate_order": 3, "gate_id": "SCHLO_LOCALIZATION", "required": "true", "rule": "directionally consistent in both source replicates using replicate-by-cell-type pseudobulk"},
        {"gate_order": 4, "gate_id": "HUMAN_REVERSAL", "required": "true", "rule": "lifestyle response-vs-nonresponse q<0.05 and leave-one-participant-out stable"},
        {"gate_order": 5, "gate_id": "ORTHOGONAL_REVERSAL", "required": "true", "rule": "RYGB direction agrees and ACMSD or prespecified mechanistic sensitivity agrees"},
        {"gate_order": 6, "gate_id": "PROVENANCE", "required": "true", "rule": "no source-headline gene, ambiguous mapping, inferred donor, or outcome-selected subset is load-bearing"},
    ]
    write_tsv(
        CANDIDATE_ROOT / "main_figure_gate_spec.tsv",
        gate_rows,
        ["gate_order", "gate_id", "required", "rule"],
    )

    known_findings = [
        {"dataset_id": "GSE200418", "finding": "source reports donor-paired nutrient-condition effects in human PCLS", "use": "source reproduction only; not atlas support"},
        {"dataset_id": "GSE207889", "finding": "source reports injury-specific multicellular HLO responses", "use": "source reproduction only; no cell-level inference"},
        {"dataset_id": "GSE253380", "finding": "source reports TLC-065/ACMSD-inhibitor effects in HLO models", "use": "source reproduction only; line-generalization source-gated"},
        {"dataset_id": "GSE106737", "finding": "source defines RYGB responders, lifestyle responders, lifestyle nonresponders, and baseline-only participants", "use": "source reproduction only"},
        {"dataset_id": "CLCC1", "finding": "source reports genome-wide and selected-library neutral-lipid CRISPR screens", "use": "source reproduction only; universes never mixed"},
        {"dataset_id": "MYOJIN_HLF", "finding": "sealed atlas-class/program analysis is nonconfirmatory for palmitate survival", "use": "reuse exact terminal artifacts; no re-mining"},
    ]
    write_tsv(
        CANDIDATE_ROOT / "known_source_findings.tsv",
        known_findings,
        ["dataset_id", "finding", "use"],
    )

    git_commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True
    ).strip()
    contract = {
        "candidate_id": CANDIDATE_ID,
        "sealed_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "git_commit": git_commit,
        "registry_rows": len(registry),
        "membership_rows": len(memberships),
        "primary_membership_rows": len(primary_memberships),
        "evidence_class_rows": len(classes),
        "primary_program_uids": sorted(primary_uids),
        "input_manifest": manifest,
        "config_manifest": config_manifest,
        "prohibitions": [
            "no program discovery, renaming, reweighting, or outcome-selected subset",
            "no cell, slice, guide, well, AOI, or repeated biopsy treated as a donor",
            "no heterogeneous effect-unit universal score",
            "no alternative dataset or contrast search after unsealing",
            "no canonical output writes or figure promotion from this workstream",
        ],
    }
    contract["specification_sha256"] = stable_json_sha256(contract)
    write_tsv(
        CANDIDATE_ROOT / "frozen_input_manifest.tsv",
        manifest,
        ["input_id", "source_path", "snapshot_path", "sha256", "size_bytes", "role"],
    )
    write_tsv(
        CANDIDATE_ROOT / "config_manifest.tsv",
        config_manifest,
        ["config", "snapshot", "sha256", "size_bytes"],
    )
    atomic_write_text(
        CANDIDATE_ROOT / "SEALED.json", json.dumps(contract, indent=2, sort_keys=True) + "\n"
    )
    atomic_write_text(
        CANDIDATE_ROOT / "SPECIFICATION_SHA256", contract["specification_sha256"] + "\n"
    )
    print(json.dumps({"status": "sealed", **contract}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
