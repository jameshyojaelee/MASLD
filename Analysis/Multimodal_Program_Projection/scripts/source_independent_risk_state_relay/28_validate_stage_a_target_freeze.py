#!/usr/bin/env python3
"""Independently rederive the Plan 45 Stage-A target-freeze verdict."""

from __future__ import annotations

import json
from collections import defaultdict

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


DIRECT_STRATA = {
    "direct_masld_mash_diagnosis",
    "mri_pdff_or_histologic_steatosis",
}
PRIMARY_RECIPIENTS = {
    "Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Macrophages"
}


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def rank_key(row: dict[str, str]) -> tuple[object, ...]:
    return (
        -int(row["n_concordant_gwas_evidence_families"]),
        -float(row["best_orientation_consensus"]),
        -float(row["best_susie_pp4"]),
        -float(row["selected_exact_shared_posterior"]),
        row["coarse_locus_uid"],
        row["ensembl_id"],
    )


def independent_selection(rows: list[dict[str, str]]) -> tuple[str, set[str]]:
    eligible = [row for row in rows if yes(row["candidate_gate_pass"])]
    by_locus: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in eligible:
        by_locus[row["coarse_locus_uid"]].append(row)
    pool = sorted((sorted(values, key=rank_key)[0] for values in by_locus.values()), key=rank_key)
    if not pool:
        return "no_target_freeze", set()
    directions = {row["oriented_risk_effect"] for row in pool}
    strata = {
        value for row in pool for value in row["phenotype_strata"].split(";") if value
    }
    if len(pool) < 4 or len(directions) < 2 or not DIRECT_STRATA <= strata:
        return "single_locus_mechanism_only", {pool[0]["target_uid"]}

    selected: list[dict[str, str]] = []
    remaining = list(pool)
    uncovered_directions = set(directions)
    uncovered_strata = set(DIRECT_STRATA)
    while remaining and len(selected) < min(6, len(pool)):
        ordered = sorted(
            remaining,
            key=lambda row: (
                -(
                    int(row["oriented_risk_effect"] in uncovered_directions)
                    + len(set(row["phenotype_strata"].split(";")) & uncovered_strata)
                ),
                *rank_key(row),
            ),
        )
        chosen = ordered[0]
        selected.append(chosen)
        remaining.remove(chosen)
        uncovered_directions.discard(chosen["oriented_risk_effect"])
        uncovered_strata -= set(chosen["phenotype_strata"].split(";"))
    return "balanced_four_to_six_locus_screen", {row["target_uid"] for row in selected}


def main() -> None:
    seal_path = CANDIDATE_ROOT / "FROZEN_STAGE_A_TARGETS.json"
    if not seal_path.is_file():
        raise RuntimeError("Missing Stage-A target-freeze seal")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    allowed_statuses = {
        "stage_a_screen_targets_frozen",
        "single_locus_stage_a_target_frozen_architecture_generalization_prohibited",
        "no_stage_a_target_frozen",
    }
    if seal.get("status") not in allowed_statuses:
        raise RuntimeError("Invalid Stage-A target-freeze status")
    if seal.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Stage-A target freeze reports scientific outcome access")
    if seal.get("stage_b_design_frozen") is not False:
        raise RuntimeError("Stage-A adjudication improperly froze Stage B")
    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Stage-A target-freeze hash mismatch: {name}")

    manifests = read_tsv(CANDIDATE_ROOT / "target_freeze_input_manifest.tsv")
    for row in manifests:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Target-freeze source drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Target-freeze source hash drift: {path}")
    template_rows = [row for row in manifests if row["role"] == "stage_a_template"]
    if len(template_rows) != 1:
        raise RuntimeError("Target-freeze manifest lacks one Stage-A template seal")
    template_seal_path = PROJECT_ROOT / template_rows[0]["source_path"]
    template_root = template_seal_path.parent
    template_seal = json.loads(template_seal_path.read_text(encoding="utf-8"))
    for stem, expected in template_seal["output_sha256"].items():
        path = template_root / f"{stem}.tsv"
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Stage-A template source drift: {path}")
    condition_rows = read_tsv(template_root / "stage_a_condition_template.tsv")
    expected_conditions = {row["arm_id"]: row for row in condition_rows}

    registry = read_tsv(CANDIDATE_ROOT / "frozen_target_registry.tsv")
    guides = read_tsv(CANDIDATE_ROOT / "frozen_stage_a_guides.tsv")
    exact = read_tsv(CANDIDATE_ROOT / "stage_b_exact_editability_registry.tsv")
    design = read_tsv(CANDIDATE_ROOT / "stage_a_design.tsv")
    status = read_tsv(CANDIDATE_ROOT / "target_freeze_status.tsv")
    if len(status) != 1:
        raise RuntimeError("Target-freeze status table must have one row")
    mode, expected_uids = independent_selection(registry)
    observed_uids = {
        row["target_uid"] for row in registry
        if row["selection_status"] == "selected_for_stage_a"
    }
    if mode != seal["selection_mode"] or expected_uids != observed_uids:
        raise RuntimeError("Independent target selection does not reproduce the seal")
    if observed_uids != set(seal["frozen_target_uids"]):
        raise RuntimeError("Frozen target universe differs from the seal")
    if len(observed_uids) != int(seal["n_stage_a_targets_frozen"]):
        raise RuntimeError("Frozen target count differs from the seal")
    if bool(observed_uids) != seal["experimental_targets_frozen"]:
        raise RuntimeError("Experimental target-freeze boolean mismatch")

    selected_registry = [row for row in registry if row["target_uid"] in observed_uids]
    if any(not yes(row["candidate_gate_pass"]) for row in selected_registry):
        raise RuntimeError("An ineligible target was selected")
    if any(row["dominant_celltype"] != "Hepatocytes" for row in selected_registry):
        raise RuntimeError("A selected target violates the hepatocyte-source gate")
    if any(
        row["accessibility_evidence_class"] != "two_cohort_hepatocyte_accessibility"
        for row in selected_registry
    ):
        raise RuntimeError("A selected target lacks independent hepatocyte accessibility")
    if len({row["coarse_locus_uid"] for row in selected_registry}) != len(selected_registry):
        raise RuntimeError("More than one selected target came from one physical locus")
    if mode == "balanced_four_to_six_locus_screen":
        if not 4 <= len(selected_registry) <= 6:
            raise RuntimeError("Balanced screen is outside the frozen 4-6 locus range")
        if len({row["oriented_risk_effect"] for row in selected_registry}) < 2:
            raise RuntimeError("Balanced screen lacks both risk-expression directions")
        observed_strata = {
            value for row in selected_registry
            for value in row["phenotype_strata"].split(";") if value
        }
        if not DIRECT_STRATA <= observed_strata:
            raise RuntimeError("Balanced screen lacks both direct phenotype strata")
    elif mode == "single_locus_mechanism_only" and len(selected_registry) != 1:
        raise RuntimeError("Single-locus fallback did not freeze exactly one target")

    guides_by_target: dict[str, list[dict[str, str]]] = defaultdict(list)
    exact_by_target: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in guides:
        guides_by_target[row["target_uid"]].append(row)
    for row in exact:
        exact_by_target[row["target_uid"]].append(row)
    if set(guides_by_target) != observed_uids or set(exact_by_target) != observed_uids:
        raise RuntimeError("Guide/editability tables do not cover exactly the frozen targets")
    registry_by_uid = {row["target_uid"]: row for row in selected_registry}
    for target_uid in observed_uids:
        target_guides = guides_by_target[target_uid]
        target_exact = exact_by_target[target_uid]
        if len(target_guides) != 2 or len({row["guide_sequence"] for row in target_guides}) != 2:
            raise RuntimeError(f"Stage-A target lacks two nonidentical guides: {target_uid}")
        if len(target_exact) != 2 or len({row["pegrna_id"] for row in target_exact}) != 2:
            raise RuntimeError(f"Stage-B editability lacks two pegRNAs: {target_uid}")
        expected_mode = {
            "risk_increases_expression": "CRISPRa",
            "risk_decreases_expression": "CRISPRi",
        }[registry_by_uid[target_uid]["oriented_risk_effect"]]
        if any(row["perturbation_mode"] != expected_mode for row in target_guides):
            raise RuntimeError(f"Guide mode disagrees with genetic orientation: {target_uid}")
        if len({row["design_uid"] for row in target_exact}) != 1:
            raise RuntimeError(f"PegRNAs do not edit the same exact variant: {target_uid}")

    design_by_target_guide_arm: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in design:
        design_by_target_guide_arm[(row["target_uid"], row["guide_id"], row["arm_id"])].append(row)
    for target_uid in observed_uids:
        for guide in guides_by_target[target_uid]:
            target_rows = [
                row for row in design
                if row["target_uid"] == target_uid and row["guide_id"] == guide["guide_id"]
            ]
            if not target_rows:
                raise RuntimeError(f"Frozen guide lacks Stage-A design rows: {guide['guide_id']}")
            observed_arms = {row["arm_id"] for row in target_rows}
            if observed_arms != set(expected_conditions):
                raise RuntimeError(
                    f"Stage-A arm universe mismatch: {target_uid}/{guide['guide_id']}"
                )
            for arm_id in observed_arms:
                arm_rows = design_by_target_guide_arm[(target_uid, guide["guide_id"], arm_id)]
                permitted = yes(expected_conditions[arm_id]["recipient_score_permitted"])
                if any(
                    row["recipient_score_permitted"]
                    != expected_conditions[arm_id]["recipient_score_permitted"]
                    for row in arm_rows
                ):
                    raise RuntimeError(
                        f"Recipient-permission drift: {target_uid}/{guide['guide_id']}/{arm_id}"
                    )
                observed_recipients = {row["recipient_lineage"] for row in arm_rows}
                expected_recipients = PRIMARY_RECIPIENTS if permitted else {"not_applicable"}
                if observed_recipients != expected_recipients:
                    raise RuntimeError(
                        f"Recipient assignment mismatch: {target_uid}/{guide['guide_id']}/{arm_id}"
                    )

    print(
        "STAGE_A_TARGET_FREEZE_VALIDATION_PASS "
        f"mode={mode} targets={len(observed_uids)} guides={len(guides)} "
        f"stage_b_design_frozen=false outcomes_opened=false"
    )


if __name__ == "__main__":
    main()
