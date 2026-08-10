#!/usr/bin/env python3
"""Create an editable response scaffold without mutating the sealed handoff."""

from __future__ import annotations

import datetime as dt
import json
import os
from pathlib import Path

from relay_common import (
    CANDIDATE_ID, CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_json,
    atomic_write_text, read_tsv, sha256_file, write_tsv,
)


RESPONSE_FILES = [
    "platform_feasibility_questions.tsv",
    "role_and_blinding_firewall.tsv",
    "collaborator_deliverable_register.tsv",
]


def candidate_source(variable: str) -> Path:
    raw = os.environ.get(variable, "").strip()
    if not raw:
        raise RuntimeError(f"{variable} is required")
    path = Path(raw)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    path = path.resolve()
    allowed = (PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates").resolve()
    if allowed not in path.parents:
        raise RuntimeError(f"{variable} escapes the candidate root: {path}")
    return path


def validate_handoff(root: Path) -> tuple[Path, dict[str, object]]:
    seal_path = root / "STAGE_A_COLLABORATOR_HANDOFF_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if (
        seal.get("status") != "sealed_target_independent_stage_a_collaborator_handoff"
        or seal.get("bundle_read_only") is not True
        or seal.get("experimental_targets_frozen") is not False
        or seal.get("scientific_outcomes_inspected") is not False
    ):
        raise RuntimeError("Invalid or outcome-open collaborator handoff")
    for relative, expected in seal["output_sha256"].items():  # type: ignore[union-attr]
        path = root / str(relative)
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Collaborator handoff drift: {relative}")
    return seal_path, seal


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite response scaffold: {CANDIDATE_ROOT}")
    handoff_root = candidate_source("PLAN45_COLLABORATOR_HANDOFF_ROOT")
    handoff_seal_path, handoff_seal = validate_handoff(handoff_root)
    CANDIDATE_ROOT.mkdir(parents=True)

    outputs: list[Path] = []
    for name in RESPONSE_FILES:
        source = handoff_root / name
        destination = CANDIDATE_ROOT / name
        atomic_write_text(destination, source.read_text(encoding="utf-8"))
        outputs.append(destination)

    identity_path = CANDIDATE_ROOT / "response_identity.tsv"
    write_tsv(
        identity_path,
        [{
            "response_package_id": CANDIDATE_ID,
            "parent_candidate_id": handoff_root.name,
            "parent_seal_sha256": sha256_file(handoff_seal_path),
            "organization": "",
            "prepared_by": "",
            "created_utc": "",
            "response_status": "editable_draft",
        }],
        [
            "response_package_id", "parent_candidate_id", "parent_seal_sha256",
            "organization", "prepared_by", "created_utc", "response_status",
        ],
    )
    outputs.append(identity_path)
    evidence_path = CANDIDATE_ROOT / "evidence_file_register.tsv"
    write_tsv(
        evidence_path,
        [],
        ["evidence_id", "relative_path", "description", "sha256", "size_bytes"],
    )
    outputs.append(evidence_path)
    instructions_path = CANDIDATE_ROOT / "EDITABLE_RESPONSE_README.md"
    atomic_write_text(
        instructions_path,
        """# Editable collaborator response

This directory is an editable response workspace, not a sealed scientific
release. Do not edit the parent handoff candidate.

1. Fill every answer/evidence/owner/timestamp field in the feasibility table.
2. Assign and timestamp every role while preserving the frozen firewall.
3. Set deliverable D01 to `complete`; other statuses may remain `not_started`.
4. Copy supporting files under `evidence/` and register every file with its
   relative path, SHA256, byte size, and a unique evidence ID.
5. Complete the organization, preparer, and UTC timestamp in
   `response_identity.tsv`, then change `response_status` to `submitted`.
6. Run the separate adjudication command with a new immutable candidate ID.

Do not add target identities, guides, allocations, or scientific outcomes.
""",
    )
    outputs.append(instructions_path)
    marker = {
        "status": "editable_collaborator_response_scaffold",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "response_package_id": CANDIDATE_ID,
        "parent_candidate_id": handoff_root.name,
        "parent_seal_sha256": sha256_file(handoff_seal_path),
        "parent_bundle_read_only": handoff_seal.get("bundle_read_only"),
        "mutable_response_workspace": True,
        "target_identity_present": False,
        "scientific_outcomes_present": False,
        "target_disclosure_permitted": False,
        "producer_sha256": sha256_file(SCRIPT_ROOT / "45_build_collaborator_response_scaffold.py"),
        "initial_template_sha256": {
            str(path.relative_to(CANDIDATE_ROOT)): sha256_file(path) for path in outputs
        },
        "next_gate": "completed response adjudication under a distinct immutable candidate ID",
    }
    atomic_write_json(CANDIDATE_ROOT / "DRAFT_EDITABLE.json", marker)
    print(json.dumps(marker, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
