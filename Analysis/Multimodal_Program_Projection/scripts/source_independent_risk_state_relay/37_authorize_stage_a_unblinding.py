#!/usr/bin/env python3
"""Authorize Stage-A condition unblinding only after blinded QC is immutable."""

from __future__ import annotations

import csv
import datetime as dt
import json
import os
from collections import defaultdict
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT, PROJECT_ROOT, atomic_write_json, read_tsv, sha256_file,
    write_tsv,
)


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def utc(value: str) -> dt.datetime:
    parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise RuntimeError(f"Timestamp lacks timezone: {value!r}")
    return parsed.astimezone(dt.timezone.utc)


def candidate_source(variable: str) -> Path:
    raw = os.environ.get(variable, "").strip()
    if not raw:
        raise RuntimeError(f"{variable} is required")
    path = Path(raw)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    path = path.resolve()
    allowed = (
        PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates"
    ).resolve()
    if allowed not in path.parents:
        raise RuntimeError(f"{variable} escapes the candidate root: {path}")
    return path


def header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def exact_rows(input_root: Path, contract_root: Path, basename: str) -> list[dict[str, str]]:
    path = input_root / basename
    template = contract_root / basename.replace(".tsv", "_template.tsv")
    if not path.is_file() or not template.is_file() or header(path) != header(template):
        raise RuntimeError(f"Unblinding schema/source failure: {basename}")
    rows = read_tsv(path)
    if not rows:
        raise RuntimeError(f"Unblinding input is empty: {basename}")
    return rows


def pair_map(arms: set[str]) -> dict[str, tuple[str, str]]:
    transfer = (
        ("TRW_CTRL", "TRW_RISK")
        if {"TRW_CTRL", "TRW_RISK"} <= arms
        else ("CM_CTRL", "CM_RISK")
    )
    result = {
        "source": ("SRC_CTRL", "SRC_RISK"),
        "mosaic": ("MOS_CTRL", "MOS_RISK"),
        "reciprocal": ("RECIP_CTRL", "RECIP_RISK"),
        "transfer": transfer,
    }
    if not set().union(*(set(pair) for pair in result.values())) <= arms:
        raise RuntimeError("Unblinding key lacks a complete prespecified arm-pair universe")
    return result


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite Stage-A unblinding authorization: {CANDIDATE_ROOT}")
    execution_root = candidate_source("PLAN45_STAGE_A_EXECUTION_ROOT")
    qc_root = candidate_source("PLAN45_BLINDED_QC_ROOT")
    contract_root = candidate_source("PLAN45_STAGE_A_QC_CONTRACT_ROOT")
    input_root = candidate_source("PLAN45_UNBLINDING_INPUT_ROOT")

    execution_seal_path = execution_root / "FROZEN_STAGE_A_EXECUTION.json"
    qc_seal_path = qc_root / "BLINDED_STAGE_A_QC_SEALED.json"
    contract_seal_path = contract_root / "STAGE_A_QC_UNBLINDING_CONTRACT_SEALED.json"
    execution_seal = json.loads(execution_seal_path.read_text(encoding="utf-8"))
    qc_seal = json.loads(qc_seal_path.read_text(encoding="utf-8"))
    contract_seal = json.loads(contract_seal_path.read_text(encoding="utf-8"))
    if execution_seal.get("status") != "stage_a_execution_frozen_ready_for_blinded_outcome_generation":
        raise RuntimeError("Unblinding requires a validated execution freeze")
    if qc_seal.get("status") not in {
        "blinded_stage_a_qc_sealed", "blinded_stage_a_qc_sealed_no_samples_pass"
    }:
        raise RuntimeError("Unblinding requires a sealed blinded-QC release")
    if not qc_seal.get("sample_exclusions_frozen") or qc_seal.get("scientific_condition_labels_opened") is not False:
        raise RuntimeError("Blinded-QC release did not freeze exclusions before labels")
    if contract_seal.get("status") != "sealed_outcome_blind_stage_a_qc_unblinding_contract":
        raise RuntimeError("Invalid QC/unblinding contract")
    for root, seal in [(execution_root, execution_seal), (qc_root, qc_seal), (contract_root, contract_seal)]:
        for name, expected in seal["output_sha256"].items():
            path = root / name
            if not path.is_file() or sha256_file(path) != expected:
                raise RuntimeError(f"Unblinding dependency drift: {path}")

    qc_manifest = read_tsv(qc_root / "blinded_qc_input_manifest.tsv")
    qc_execution = next(row for row in qc_manifest if row["role"] == "stage_a_execution_seal")
    qc_contract = next(row for row in qc_manifest if row["role"] == "qc_unblinding_contract_seal")
    if (
        qc_execution["source_path"] != str(execution_seal_path.relative_to(PROJECT_ROOT))
        or qc_execution["sha256"] != sha256_file(execution_seal_path)
        or qc_contract["source_path"] != str(contract_seal_path.relative_to(PROJECT_ROOT))
        or qc_contract["sha256"] != sha256_file(contract_seal_path)
    ):
        raise RuntimeError("Blinded-QC release is not bound to the supplied execution/contract")

    key_rows = exact_rows(input_root, contract_root, "unblinding_key.tsv")
    signoff_rows = exact_rows(input_root, contract_root, "unblinding_signoff.tsv")
    if len(signoff_rows) != 1:
        raise RuntimeError("Unblinding signoff must contain exactly one row")
    signoff = signoff_rows[0]
    key_path = input_root / "unblinding_key.tsv"
    if (
        Path(signoff["qc_seal_path"]).as_posix()
        != qc_seal_path.relative_to(PROJECT_ROOT).as_posix()
        or signoff["qc_seal_sha256"] != sha256_file(qc_seal_path)
        or Path(signoff["unblinding_key_path"]).as_posix()
        != key_path.relative_to(PROJECT_ROOT).as_posix()
        or signoff["unblinding_key_sha256"] != sha256_file(key_path)
    ):
        raise RuntimeError("Unblinding signoff is not bound to the exact QC seal/key")
    if not yes(signoff["key_opened_after_qc_seal"]) or yes(signoff["scientific_outcomes_inspected"]):
        raise RuntimeError("Unblinding order/outcome firewall failed")
    if not signoff["data_manager"] or not signoff["analysis_lead"] or signoff["data_manager"] == signoff["analysis_lead"] or not signoff["signed_utc"]:
        raise RuntimeError("Unblinding requires separate data-manager and analysis identities")
    if utc(signoff["signed_utc"]) <= utc(str(qc_seal["created_utc"])):
        raise RuntimeError("Unblinding signoff predates or equals the blinded-QC seal")

    randomization = read_tsv(execution_root / "frozen_randomization_manifest.tsv")
    identity_fields = [
        "sample_id", "blinded_label", "biological_unit_id", "background_id",
        "differentiation_id", "target_uid", "guide_id", "arm_id", "time_role",
        "challenge_role", "challenge_id", "assay_id",
    ]
    expected = {
        row["sample_id"]: tuple(row[field] for field in identity_fields)
        for row in randomization
    }
    observed = {
        row["sample_id"]: tuple(row[field] for field in identity_fields)
        for row in key_rows
    }
    if len(observed) != len(key_rows) or observed != expected:
        raise RuntimeError("Unblinding key is not an exact bijection to the frozen randomization")
    for row in key_rows:
        if not yes(row["key_generated_before_qc"]) or not yes(row["key_opened_after_qc_seal"]):
            raise RuntimeError("Unblinding key chronology is invalid")
        if not yes(row["review_concordant"]) or row["reviewer_1"] == row["reviewer_2"]:
            raise RuntimeError("Unblinding key lacks independent review")

    qc_status = {
        row["sample_id"]: row
        for row in read_tsv(qc_root / "blinded_sample_qc_status.tsv")
    }
    if set(qc_status) != set(expected):
        raise RuntimeError("Blinded-QC sample universe differs from the unblinding key")
    analysis_rows = []
    for row in key_rows:
        qc = qc_status[row["sample_id"]]
        analysis_rows.append({
            **{field: row[field] for field in identity_fields},
            "technical_qc_pass": qc["technical_qc_pass"],
            "preunblinding_exclusion_status": qc["preunblinding_exclusion_status"],
            "preunblinding_exclusion_reason": qc["preunblinding_exclusion_reason"],
            "analysis_status": (
                "eligible_after_preunblinding_qc"
                if yes(qc["technical_qc_pass"])
                else "excluded_before_condition_key_opened"
            ),
        })

    arms = {row["arm_id"] for row in key_rows}
    pairs = pair_map(arms)
    arm_to_pair = {arm: pair_id for pair_id, pair in pairs.items() for arm in pair}
    grouped: dict[tuple[str, ...], dict[str, list[dict[str, str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in analysis_rows:
        pair_id = arm_to_pair[row["arm_id"]]
        group_key = (
            row["target_uid"], row["guide_id"], row["biological_unit_id"],
            row["background_id"], row["differentiation_id"], row["assay_id"],
            row["time_role"], row["challenge_role"], pair_id,
        )
        grouped[group_key][row["arm_id"]].append(row)

    pair_rows = []
    for key in sorted(grouped):
        target_uid, guide_id, unit_id, background_id, differentiation_id, assay_id, time_role, challenge_role, pair_id = key
        ctrl, risk = pairs[pair_id]
        by_arm = grouped[key]
        ctrl_pass = any(row["analysis_status"] == "eligible_after_preunblinding_qc" for row in by_arm.get(ctrl, []))
        risk_pass = any(row["analysis_status"] == "eligible_after_preunblinding_qc" for row in by_arm.get(risk, []))
        pair_rows.append({
            "target_uid": target_uid,
            "guide_id": guide_id,
            "biological_unit_id": unit_id,
            "background_id": background_id,
            "differentiation_id": differentiation_id,
            "assay_id": assay_id,
            "time_role": time_role,
            "challenge_role": challenge_role,
            "pair_id": pair_id,
            "control_arm": ctrl,
            "risk_arm": risk,
            "control_qc_pass": str(ctrl_pass).lower(),
            "risk_qc_pass": str(risk_pass).lower(),
            "post_qc_pair_complete": str(ctrl_pass and risk_pass).lower(),
        })

    comparison_groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in pair_rows:
        key = (
            row["target_uid"], row["guide_id"], row["assay_id"], row["time_role"],
            row["challenge_role"], row["pair_id"],
        )
        comparison_groups[key].append(row)
    minimum_backgrounds = int(contract_seal["minimum_backgrounds_per_claim_bearing_comparison"])
    comparison_rows = []
    for key in sorted(comparison_groups):
        rows = comparison_groups[key]
        complete = [row for row in rows if yes(row["post_qc_pair_complete"])]
        n_backgrounds = len({row["background_id"] for row in complete})
        comparison_rows.append({
            "target_uid": key[0], "guide_id": key[1], "assay_id": key[2],
            "time_role": key[3], "challenge_role": key[4], "pair_id": key[5],
            "n_complete_biological_units": len(complete),
            "n_complete_backgrounds": n_backgrounds,
            "claim_bearing_analysis_eligible": str(n_backgrounds >= minimum_backgrounds).lower(),
            "eligibility_reason": (
                "pass_minimum_backgrounds"
                if n_backgrounds >= minimum_backgrounds
                else "insufficient_post_qc_backgrounds"
            ),
        })

    input_manifest = []
    for role, path in [
        ("stage_a_execution_seal", execution_seal_path),
        ("blinded_qc_seal", qc_seal_path),
        ("qc_unblinding_contract_seal", contract_seal_path),
        ("unblinding_key", key_path),
        ("unblinding_signoff", input_root / "unblinding_signoff.tsv"),
    ]:
        input_manifest.append({
            "role": role,
            "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })

    outputs: list[Path] = []
    analysis_path = CANDIDATE_ROOT / "unblinded_analysis_sample_manifest.tsv"
    write_tsv(analysis_path, analysis_rows, list(analysis_rows[0]))
    outputs.append(analysis_path)
    pair_path = CANDIDATE_ROOT / "post_qc_pair_completeness.tsv"
    write_tsv(pair_path, pair_rows, list(pair_rows[0]))
    outputs.append(pair_path)
    comparison_path = CANDIDATE_ROOT / "post_qc_comparison_eligibility.tsv"
    write_tsv(comparison_path, comparison_rows, list(comparison_rows[0]))
    outputs.append(comparison_path)
    manifest_path = CANDIDATE_ROOT / "unblinding_input_manifest.tsv"
    write_tsv(manifest_path, input_manifest, ["role", "source_path", "size_bytes", "sha256"])
    outputs.append(manifest_path)
    signoff_path = CANDIDATE_ROOT / "frozen_unblinding_signoff.tsv"
    write_tsv(signoff_path, signoff_rows, header(input_root / "unblinding_signoff.tsv"))
    outputs.append(signoff_path)

    payload = {
        "status": "stage_a_unblinding_authorized_for_prespecified_post_qc_pairs",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "n_samples": len(analysis_rows),
        "n_samples_included_after_blinded_qc": sum(
            row["analysis_status"] == "eligible_after_preunblinding_qc"
            for row in analysis_rows
        ),
        "n_post_qc_pairs": len(pair_rows),
        "n_post_qc_complete_pairs": sum(yes(row["post_qc_pair_complete"]) for row in pair_rows),
        "n_claim_bearing_eligible_comparisons": sum(
            yes(row["claim_bearing_analysis_eligible"]) for row in comparison_rows
        ),
        "sample_exclusions_frozen_before_unblinding": True,
        "scientific_condition_labels_opened": True,
        "scientific_outcomes_inspected": False,
        "unblinding_authorized": True,
        "stage_b_design_frozen": False,
        "next_gate": "analyze only the frozen post-QC complete Stage-A comparison universe",
        "output_sha256": {path.name: sha256_file(path) for path in outputs},
    }
    atomic_write_json(CANDIDATE_ROOT / "STAGE_A_UNBLINDING_AUTHORIZED.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
