#!/usr/bin/env python3
"""Prepare checksum-bound coordinator selection and manual PASS-06 attestations.

This helper records an attestation; it does not substitute for scientific or UI
review.  Every manual-review confirmation must be supplied explicitly on the
command line, and the exact MASLD Gene Catalog analysis release ID is frozen.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import generate_evidence_passports as contract


SCRIPT_PATH = Path(__file__).resolve()


def _require_text(value: str, label: str) -> str:
    value = value.strip()
    if not value:
        raise contract.PassportContractError("ATTESTATION_FIELD", f"{label} is required")
    return value


def sign_selection(args: argparse.Namespace) -> Path:
    selection = args.selection.resolve()
    if not selection.is_file():
        raise contract.PassportContractError("SELECTION_MISSING", str(selection))
    with selection.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != contract.SELECTION_REQUIRED_COLUMNS:
            raise contract.PassportContractError(
                "SELECTION_SCHEMA",
                f"expected exact columns {contract.SELECTION_REQUIRED_COLUMNS}",
            )
        rows = list(reader)
    if not rows or any(row["fixture_only"].lower() != "false" for row in rows):
        raise contract.PassportContractError(
            "SELECTION_FIXTURE_FLAG", "production selection rows must set fixture_only=false"
        )
    input_ids = {row["input_id"] for row in rows}
    signed_ids = {
        args.gen_terminal_closure_input_id,
        args.gen_terminal_ready_input_id,
        args.plan13_terminal_ready_input_id,
    }
    if not signed_ids.issubset(input_ids):
        raise contract.PassportContractError(
            "SELECTION_SIGNED_GATE_REFERENCE", str(sorted(signed_ids - input_ids))
        )
    signature_path = selection.with_name("passport_input_selection.signature.json")
    contract.atomic_write_json(
        signature_path,
        {
            "attestation_version": contract.SELECTION_ATTESTATION_VERSION,
            "selection_sha256": contract.sha256_file(selection),
            "signed_by": _require_text(args.signed_by, "signed_by"),
            "signed_at_utc": _require_text(args.signed_at_utc, "signed_at_utc"),
            "decision_register_id": _require_text(
                args.decision_register_id, "decision_register_id"
            ),
            "authority_document": _require_text(
                args.authority_document, "authority_document"
            ),
            "fixture_only": False,
            "coordinator_attested": True,
            "analysis_release_id": contract.PRODUCTION_ANALYSIS_RELEASE_ID,
            "gen_terminal_closure_input_id": args.gen_terminal_closure_input_id,
            "gen_terminal_ready_input_id": args.gen_terminal_ready_input_id,
            "plan13_terminal_ready_input_id": args.plan13_terminal_ready_input_id,
        },
    )
    return signature_path


def prepare_manual(args: argparse.Namespace) -> tuple[Path, Path]:
    import validate_evidence_passports as validator

    confirmations = {
        "call_states_reviewed": args.confirm_call_states,
        "provenance_states_reviewed": args.confirm_provenance_states,
        "visible_boundary_pass": args.confirm_visible_boundary,
        "visible_testability_pass": args.confirm_visible_testability,
        "visible_source_dependence_pass": args.confirm_visible_source_dependence,
        "visible_falsifier_pass": args.confirm_visible_falsifier,
    }
    missing = [name for name, passed in confirmations.items() if not passed]
    if missing:
        raise contract.PassportContractError(
            "MANUAL_ACCEPTANCE_CONFIRMATION",
            f"manual reviewer must explicitly confirm: {missing}",
        )
    selection = args.selection.resolve()
    signature_path = selection.with_name("passport_input_selection.signature.json")
    if not selection.is_file() or not signature_path.is_file():
        raise contract.PassportContractError(
            "SIGNED_SELECTION_MISSING", f"selection={selection};signature={signature_path}"
        )
    contract.validate_signed_selection(selection, allow_fixture=False)
    signature = json.loads(signature_path.read_text(encoding="utf-8"))
    if signature.get("analysis_release_id") != contract.PRODUCTION_ANALYSIS_RELEASE_ID:
        raise contract.PassportContractError(
            "MANUAL_ACCEPTANCE_RELEASE", str(signature.get("analysis_release_id"))
        )
    bundle = args.bundle.resolve()
    if not bundle.is_dir() or bundle.is_symlink():
        raise contract.PassportContractError(
            "MANUAL_ACCEPTANCE_BUNDLE", f"unsafe bundle={bundle}"
        )
    if selection != (bundle / "passport_input_selection.tsv").resolve():
        raise contract.PassportContractError(
            "MANUAL_ACCEPTANCE_SELECTION",
            "the reviewed bundle must contain the exact signed selection",
        )
    validator.check_candidate_ui(bundle, fixture_mode=False)
    ui_paths = {
        "ui_index_sha256": bundle / "portal_candidate/index.html",
        "ui_contract_sha256": bundle / "portal_candidate/ui_contract.json",
        "ui_source_manifest_sha256": (
            bundle / "portal_candidate/ui_source_manifest.tsv"
        ),
        "ui_review_manifest_sha256": (
            bundle / "portal_candidate/review/review_manifest.tsv"
        ),
    }
    output = args.output.resolve()
    try:
        output.relative_to(bundle)
    except ValueError:
        pass
    else:
        raise contract.PassportContractError(
            "MANUAL_ACCEPTANCE_OUTPUT",
            "manual acceptance must be staged outside the reviewed bundle",
        )
    manual_signature = output.with_name("passport_manual_acceptance.signature.json")
    if (
        output.exists()
        or output.is_symlink()
        or manual_signature.exists()
        or manual_signature.is_symlink()
    ):
        raise contract.PassportContractError(
            "MANUAL_ACCEPTANCE_EXISTS",
            "refusing to overwrite a staged manual decision record",
        )
    ui_hashes = {
        field: contract.sha256_file(path) for field, path in ui_paths.items()
    }
    row = {
        "review_id": _require_text(args.review_id, "review_id"),
        "analysis_release_id": contract.PRODUCTION_ANALYSIS_RELEASE_ID,
        "selection_sha256": contract.sha256_file(selection),
        "reviewer": _require_text(args.reviewer, "reviewer"),
        "reviewed_at_utc": _require_text(args.reviewed_at_utc, "reviewed_at_utc"),
        **ui_hashes,
        "call_states_reviewed": "true",
        "provenance_states_reviewed": "true",
        "hero_genes_reviewed": "THRB;HKDC1;GLP1R;MTARC1",
        "visible_boundary_pass": "true",
        "visible_testability_pass": "true",
        "visible_source_dependence_pass": "true",
        "visible_falsifier_pass": "true",
        "decision": "accepted",
    }
    contract.atomic_write_tsv(output, [row], contract.MANUAL_ACCEPTANCE_COLUMNS)
    contract.atomic_write_json(
        manual_signature,
        {
            "attestation_version": contract.MANUAL_ACCEPTANCE_ATTESTATION_VERSION,
            "acceptance_sha256": contract.sha256_file(output),
            **ui_hashes,
            "signed_by": _require_text(args.reviewer, "reviewer"),
            "signed_at_utc": _require_text(args.reviewed_at_utc, "reviewed_at_utc"),
            "decision_register_id": _require_text(
                args.decision_register_id, "decision_register_id"
            ),
            "fixture_only": False,
        },
    )
    return output, manual_signature


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    selection = subparsers.add_parser("selection", help="attest a completed selection TSV")
    selection.add_argument("--selection", type=Path, required=True)
    selection.add_argument("--signed-by", required=True)
    selection.add_argument("--signed-at-utc", required=True)
    selection.add_argument("--decision-register-id", required=True)
    selection.add_argument(
        "--authority-document",
        default="docs/PAPER.md",
    )
    selection.add_argument("--gen-terminal-closure-input-id", required=True)
    selection.add_argument("--gen-terminal-ready-input-id", required=True)
    selection.add_argument("--plan13-terminal-ready-input-id", required=True)

    manual = subparsers.add_parser("manual", help="record completed manual portal review")
    manual.add_argument("--selection", type=Path, required=True)
    manual.add_argument("--bundle", type=Path, required=True)
    manual.add_argument("--output", type=Path, required=True)
    manual.add_argument("--review-id", required=True)
    manual.add_argument("--reviewer", required=True)
    manual.add_argument("--reviewed-at-utc", required=True)
    manual.add_argument("--decision-register-id", required=True)
    manual.add_argument("--confirm-call-states", action="store_true")
    manual.add_argument("--confirm-provenance-states", action="store_true")
    manual.add_argument("--confirm-visible-boundary", action="store_true")
    manual.add_argument("--confirm-visible-testability", action="store_true")
    manual.add_argument("--confirm-visible-source-dependence", action="store_true")
    manual.add_argument("--confirm-visible-falsifier", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "selection":
        print(f"PASS: selection attestation={sign_selection(args)}")
    else:
        output, signature = prepare_manual(args)
        print(f"PASS: manual acceptance={output};signature={signature}")


if __name__ == "__main__":
    try:
        main()
    except contract.PassportContractError as exc:
        raise SystemExit(str(exc)) from exc
