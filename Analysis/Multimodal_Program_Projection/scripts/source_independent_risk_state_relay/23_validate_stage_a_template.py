#!/usr/bin/env python3
"""Validate the target-independent Plan 45 Stage A template."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

from relay_common import PROJECT_ROOT, sha256_file


CANDIDATE_ID = os.environ.get(
    "PLAN45_STAGE_A_TEMPLATE_ID",
    "source-independent-risk-state-relay-stage-a-template-v2-2026-08-10",
).strip()
ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates" / CANDIDATE_ID
ELIGIBLE = {"Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Macrophages"}


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def main() -> None:
    seal_path = ROOT / "STAGE_A_TEMPLATE_SEALED.json"
    names = [
        "stage_a_recipient_axis_registry", "stage_a_required_metadata_schema",
        "stage_a_condition_template", "stage_a_estimand_registry",
        "stage_a_score_contract", "stage_a_randomization_blinding_contract",
        "stage_a_template_status", "stage_a_template_input_manifest",
    ]
    paths = {name: ROOT / f"{name}.tsv" for name in names}
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("target_freeze_status") != "prohibited":
        raise RuntimeError("Stage A template improperly freezes a target")
    if seal.get("experimental_outcomes_inspected") is not False:
        raise RuntimeError("Stage A template indicates outcome access")
    if seal.get("invalid_source_monoculture_recipient_comparison_prohibited") is not True:
        raise RuntimeError("Invalid recipient/source-monoculture contrast is not prohibited")
    for name, path in paths.items():
        if not path.is_file() or sha256_file(path) != seal["output_sha256"][name]:
            raise RuntimeError(f"Stage A template artifact/hash failure: {name}")

    axes = list(rows(paths["stage_a_recipient_axis_registry"]))
    passing = {
        row["recipient_lineage"]
        for row in axes
        if row["primary_recipient_eligible"].lower() == "true"
    }
    if passing != ELIGIBLE or len(axes) != 5:
        raise RuntimeError("Stage A recipient-axis universe drift")
    hepatocyte = next(row for row in axes if row["recipient_lineage"] == "Hepatocytes")
    if hepatocyte["primary_use"] != "secondary_nonconfirmatory_only":
        raise RuntimeError("Hepatocyte improperly promoted as recipient endpoint")

    fields = {row["field"] for row in rows(paths["stage_a_required_metadata_schema"])}
    required_fields = {
        "biological_unit_id", "background_id", "differentiation_id", "well_id",
        "target_locus_uid", "target_gene", "guide_id", "risk_direction",
        "source_lineage", "edited_lineage", "recipient_lineage", "culture_context",
        "time_role", "challenge_id", "n_cells", "viability_pass",
        "recipient_substate_fractions_path", "blinded_label", "unblinding_status",
    }
    if not required_fields.issubset(fields):
        raise RuntimeError("Stage A metadata schema lacks a load-bearing field")

    conditions = list(rows(paths["stage_a_condition_template"]))
    contexts = {row["culture_context"] for row in conditions}
    required_contexts = {
        "source_only", "direct_mosaic", "transwell_recipient",
        "conditioned_medium_recipient", "reciprocal_recipient_edit",
    }
    if contexts != required_contexts:
        raise RuntimeError("Stage A condition context drift")
    if any(
        row["culture_context"] == "source_only"
        and row["recipient_score_permitted"].lower() == "true"
        for row in conditions
    ):
        raise RuntimeError("Recipient score improperly permitted in source-only culture")

    estimands = {row["estimand_id"]: row for row in rows(paths["stage_a_estimand_registry"])}
    if set(estimands) != {"CIS01", "RELAY01", "ORIGIN01", "ROUTE01", "COMP01", "PHENO01"}:
        raise RuntimeError("Stage A estimand registry drift")
    relay = estimands["RELAY01"]
    if "within direct_mosaic" not in relay["contrast"] or "four eligible recipient axes" not in relay["family"]:
        raise RuntimeError("Primary relay estimand is not the corrected within-mosaic contrast")
    if "source_only" in relay["contrast"]:
        raise RuntimeError("Primary relay estimand revives invalid source-only comparison")

    scoring = {row["rule"]: row["value"] for row in rows(paths["stage_a_score_contract"])}
    if scoring.get("coverage_gate") != "at least 500 retained genes and at least 80% original axis L1 loading":
        raise RuntimeError("Stage A score coverage gate drift")
    if set(scoring.get("recipient_axes", "").split(";")) != ELIGIBLE:
        raise RuntimeError("Stage A score includes an unauthorized recipient axis")
    if "cells/wells/lanes are technical" not in scoring.get("biological_unit", ""):
        raise RuntimeError("Stage A biological-unit firewall absent")

    forbidden_gene_names = {
        "ACSL5", "CENPQ", "ERCC2", "IL18R1", "NECAB2",
        "PNPLA6", "RANBP17", "SHMT1", "ZBTB41",
    }
    for path in paths.values():
        text = path.read_text(encoding="utf-8")
        if any(gene in text for gene in forbidden_gene_names):
            raise RuntimeError("Stage A template contains a legacy target name")
    print(
        "Plan 45 Stage A template validation passed: axes=4; contexts=5; "
        "estimands=6; metadata fields=" + str(len(fields)) + "; target freeze prohibited"
    )


if __name__ == "__main__":
    main()
