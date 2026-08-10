#!/usr/bin/env python3
"""Independently validate Stage-A unblinding chronology and post-QC pairs."""

from __future__ import annotations

import json
import datetime as dt
from collections import defaultdict

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def utc(value: str) -> dt.datetime:
    parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise RuntimeError(f"Timestamp lacks timezone: {value!r}")
    return parsed.astimezone(dt.timezone.utc)


def main() -> None:
    seal_path = CANDIDATE_ROOT / "STAGE_A_UNBLINDING_AUTHORIZED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "stage_a_unblinding_authorized_for_prespecified_post_qc_pairs":
        raise RuntimeError("Invalid Stage-A unblinding status")
    if seal.get("sample_exclusions_frozen_before_unblinding") is not True:
        raise RuntimeError("Sample exclusions were not frozen before unblinding")
    if seal.get("scientific_condition_labels_opened") is not True or seal.get("unblinding_authorized") is not True:
        raise RuntimeError("Condition-key authorization state drift")
    if seal.get("scientific_outcomes_inspected") is not False or seal.get("stage_b_design_frozen") is not False:
        raise RuntimeError("Unblinding release improperly opens outcomes or Stage B")
    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Unblinding output hash mismatch: {name}")

    manifest = read_tsv(CANDIDATE_ROOT / "unblinding_input_manifest.tsv")
    expected_roles = {
        "stage_a_execution_seal", "blinded_qc_seal", "qc_unblinding_contract_seal",
        "unblinding_key", "unblinding_signoff",
    }
    if {row["role"] for row in manifest} != expected_roles:
        raise RuntimeError("Unblinding input-manifest universe drift")
    for row in manifest:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Unblinding source drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Unblinding source hash drift: {path}")

    execution_row = next(row for row in manifest if row["role"] == "stage_a_execution_seal")
    execution_root = (PROJECT_ROOT / execution_row["source_path"]).parent
    randomization = read_tsv(execution_root / "frozen_randomization_manifest.tsv")
    analysis = read_tsv(CANDIDATE_ROOT / "unblinded_analysis_sample_manifest.tsv")
    identity_fields = [
        "sample_id", "blinded_label", "biological_unit_id", "background_id",
        "differentiation_id", "target_uid", "guide_id", "arm_id", "time_role",
        "challenge_role", "challenge_id", "assay_id",
    ]
    expected_identity = {
        row["sample_id"]: tuple(row[field] for field in identity_fields)
        for row in randomization
    }
    observed_identity = {
        row["sample_id"]: tuple(row[field] for field in identity_fields)
        for row in analysis
    }
    if len(analysis) != len(expected_identity) or observed_identity != expected_identity:
        raise RuntimeError("Independent unblinding-key bijection failed")

    qc_row = next(row for row in manifest if row["role"] == "blinded_qc_seal")
    qc_root = (PROJECT_ROOT / qc_row["source_path"]).parent
    qc_seal = json.loads((PROJECT_ROOT / qc_row["source_path"]).read_text(encoding="utf-8"))
    qc_status = {
        row["sample_id"]: row for row in read_tsv(qc_root / "blinded_sample_qc_status.tsv")
    }
    for row in analysis:
        qc = qc_status[row["sample_id"]]
        expected_status = (
            "eligible_after_preunblinding_qc"
            if yes(qc["technical_qc_pass"])
            else "excluded_before_condition_key_opened"
        )
        if row["analysis_status"] != expected_status:
            raise RuntimeError("Unblinded analysis status differs from frozen blinded QC")

    signoff_row = next(row for row in manifest if row["role"] == "unblinding_signoff")
    signoffs = read_tsv(PROJECT_ROOT / signoff_row["source_path"])
    if len(signoffs) != 1:
        raise RuntimeError("Independent unblinding signoff count drift")
    signoff = signoffs[0]
    if (
        not yes(signoff["key_opened_after_qc_seal"])
        or yes(signoff["scientific_outcomes_inspected"])
        or signoff["data_manager"] == signoff["analysis_lead"]
        or utc(signoff["signed_utc"]) <= utc(str(qc_seal["created_utc"]))
    ):
        raise RuntimeError("Independent unblinding chronology/signoff validation failed")

    pair_rows = read_tsv(CANDIDATE_ROOT / "post_qc_pair_completeness.tsv")
    analysis_by_key: dict[tuple[str, ...], dict[str, list[dict[str, str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    pairs = {
        "source": ("SRC_CTRL", "SRC_RISK"),
        "mosaic": ("MOS_CTRL", "MOS_RISK"),
        "reciprocal": ("RECIP_CTRL", "RECIP_RISK"),
    }
    arms = {row["arm_id"] for row in analysis}
    pairs["transfer"] = (
        ("TRW_CTRL", "TRW_RISK")
        if {"TRW_CTRL", "TRW_RISK"} <= arms
        else ("CM_CTRL", "CM_RISK")
    )
    arm_pair = {arm: pair_id for pair_id, pair in pairs.items() for arm in pair}
    for row in analysis:
        pair_id = arm_pair[row["arm_id"]]
        key = (
            row["target_uid"], row["guide_id"], row["biological_unit_id"],
            row["background_id"], row["differentiation_id"], row["assay_id"],
            row["time_role"], row["challenge_role"], pair_id,
        )
        analysis_by_key[key][row["arm_id"]].append(row)
    rederived = {}
    for key, by_arm in analysis_by_key.items():
        pair = pairs[key[-1]]
        rederived[key] = all(
            any(row["analysis_status"] == "eligible_after_preunblinding_qc" for row in by_arm.get(arm, []))
            for arm in pair
        )
    observed_pairs = {
        (
            row["target_uid"], row["guide_id"], row["biological_unit_id"],
            row["background_id"], row["differentiation_id"], row["assay_id"],
            row["time_role"], row["challenge_role"], row["pair_id"],
        ): yes(row["post_qc_pair_complete"])
        for row in pair_rows
    }
    if observed_pairs != rederived:
        raise RuntimeError("Independent post-QC pair rederivation failed")

    contract_row = next(row for row in manifest if row["role"] == "qc_unblinding_contract_seal")
    contract = json.loads((PROJECT_ROOT / contract_row["source_path"]).read_text(encoding="utf-8"))
    minimum = int(contract["minimum_backgrounds_per_claim_bearing_comparison"])
    comparison_rows = read_tsv(CANDIDATE_ROOT / "post_qc_comparison_eligibility.tsv")
    grouped: dict[tuple[str, ...], set[str]] = defaultdict(set)
    for row in pair_rows:
        if yes(row["post_qc_pair_complete"]):
            key = (
                row["target_uid"], row["guide_id"], row["assay_id"],
                row["time_role"], row["challenge_role"], row["pair_id"],
            )
            grouped[key].add(row["background_id"])
    observed_eligibility = {
        (
            row["target_uid"], row["guide_id"], row["assay_id"], row["time_role"],
            row["challenge_role"], row["pair_id"],
        ): yes(row["claim_bearing_analysis_eligible"])
        for row in comparison_rows
    }
    expected_eligibility = {
        key: len(backgrounds) >= minimum
        for key, backgrounds in grouped.items()
    }
    for key in observed_eligibility:
        expected_eligibility.setdefault(key, False)
    if observed_eligibility != expected_eligibility:
        raise RuntimeError("Independent claim-bearing eligibility rederivation failed")
    if sum(observed_eligibility.values()) != int(seal["n_claim_bearing_eligible_comparisons"]):
        raise RuntimeError("Claim-bearing comparison count drift")
    print(
        "STAGE_A_UNBLINDING_VALIDATION_PASS "
        f"samples={len(analysis)} pairs={len(pair_rows)} eligible={sum(observed_eligibility.values())} "
        "outcomes_opened=false stage_b_frozen=false"
    )


if __name__ == "__main__":
    main()
