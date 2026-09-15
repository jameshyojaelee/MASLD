#!/usr/bin/env python3
"""Correct the scGPT checkpoint bundle to match what is on disk.

The bundle was written 2026-08-22 and states that the weight bodies were not
downloaded.  That was true when written.  The bodies were acquired the next day
by job 21015743 and the bundle was never updated, so it went on asserting a
falsehood that read as licence compliance.

This makes three factual corrections and nothing else: the download flag, and
the two placeholder digests.  Every value is checked against both its expected
stale content and the live bytes before anything is written.  The declared
licence is deliberately left UNDECLARED, because correcting the record must not
quietly imply a licence that nobody has granted.

The bundle is pinned by one file, its own native_task_disposition, which is
updated in the same run so no binding is left stranded.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


BUNDLE = "config/artifacts/models/scgpt/checkpoints.json"
PINNER = "config/artifacts/models/scgpt/native_task_disposition_20260825.json"

# bundle id -> (output-file path, expected bytes, expected digest of the body)
WEIGHT_BODIES = {
    "scgpt_continual": (
        "executions/scgpt-acquisition-21015743/scgpt_continual/best_model.pt",
        207861754,
        "ad0252a1971e0cd619b7116dbab3177432236c4537225d54280a2aa7e5fe402a",
    ),
    "scgpt_whole_human": (
        "executions/scgpt-acquisition-21015743/scgpt_whole_human/best_model.pt",
        205385258,
        "6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922",
    ),
}

STALE_EVIDENCE = (
    "Weight bodies were not downloaded; their SHA-256 values remain unresolved."
)
CORRECTED_EVIDENCE = (
    "Weight bodies were acquired on 2026-08-23 by job 21015743 and are held at "
    "executions/scgpt-acquisition-21015743; their SHA-256 values were re-derived "
    "from the bodies on 2026-08-25 and are recorded above. The weight licence "
    "remains UNDECLARED: acquisition does not grant terms, and derivative-weight "
    "redistribution stays prohibited pending a declared licence."
)

STALE_BLOCKER_TEXT = "weight_content_downloaded is false."
CORRECTED_BLOCKER_TEXT = (
    "weight_content_downloaded is true as of the 2026-08-25 correction; the "
    "weight licence is UNDECLARED and derivative-weight redistribution is "
    "prohibited pending declared terms."
)


class ScgptCorrectionError(RuntimeError):
    """Raised when the tree does not match the audited pre-correction state."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(64 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify_bodies(root: Path) -> dict[str, str]:
    """Re-derive each weight digest from the body itself."""

    resolved: dict[str, str] = {}
    for bundle_id, (relative, size, expected) in WEIGHT_BODIES.items():
        body = root / relative
        if not body.is_file():
            raise ScgptCorrectionError(f"weight body missing: {body}")
        if body.stat().st_size != size:
            raise ScgptCorrectionError(
                f"{relative} is {body.stat().st_size} bytes, expected {size}"
            )
        actual = digest(body)
        if actual != expected:
            raise ScgptCorrectionError(
                f"{relative} digest {actual} differs from the audited {expected}"
            )
        resolved[bundle_id] = actual
    return resolved


def correct_bundle(document: dict, resolved: dict[str, str]) -> list[dict]:
    changes: list[dict] = []
    if document["weight_content_downloaded"] is not False:
        raise ScgptCorrectionError(
            "weight_content_downloaded is not the audited false; re-audit first"
        )
    document["weight_content_downloaded"] = True
    changes.append({"field": "weight_content_downloaded", "from": False, "to": True})

    for bundle_id, resolved_digest in resolved.items():
        artifacts = {
            item["path"]: item for item in document["bundles"][bundle_id]["artifacts"]
        }
        body = artifacts["best_model.pt"]
        if body["sha256"] != "UNRESOLVED":
            raise ScgptCorrectionError(
                f"{bundle_id} best_model.pt is already {body['sha256']}"
            )
        body["sha256"] = resolved_digest
        changes.append(
            {
                "field": f"bundles.{bundle_id}.artifacts[best_model.pt].sha256",
                "from": "UNRESOLVED",
                "to": resolved_digest,
            }
        )

    if STALE_EVIDENCE not in document["source_evidence"]:
        raise ScgptCorrectionError("source_evidence does not carry the audited claim")
    document["source_evidence"] = document["source_evidence"].replace(
        STALE_EVIDENCE, CORRECTED_EVIDENCE
    )
    changes.append({"field": "source_evidence", "from": STALE_EVIDENCE, "to": CORRECTED_EVIDENCE})

    # The licence is untouched on purpose.
    if document["weight_license"] != "UNDECLARED":
        raise ScgptCorrectionError("weight_license changed unexpectedly")
    return changes


def correct_pinner(document: dict, bundle_digest: str) -> list[dict]:
    changes: list[dict] = []
    for entry in document["bound_evidence"]:
        if entry["path"] == BUNDLE:
            changes.append(
                {"field": "bound_evidence[checkpoints.json].sha256",
                 "from": entry["sha256"], "to": bundle_digest}
            )
            entry["sha256"] = bundle_digest
    detail = document["family_specific_blockers"]["detail"]
    if STALE_BLOCKER_TEXT in detail:
        document["family_specific_blockers"]["detail"] = detail.replace(
            STALE_BLOCKER_TEXT, CORRECTED_BLOCKER_TEXT
        )
        changes.append({"field": "family_specific_blockers.detail",
                        "from": STALE_BLOCKER_TEXT, "to": CORRECTED_BLOCKER_TEXT})
    if not changes:
        raise ScgptCorrectionError("pinner did not reference the bundle")
    return changes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        raise ScgptCorrectionError("correction receipt output exists")

    resolved = verify_bodies(args.root)
    bundle_path = args.root / BUNDLE
    pinner_path = args.root / PINNER
    document = json.loads(bundle_path.read_text(encoding="utf-8"))
    before = digest(bundle_path)
    bundle_changes = correct_bundle(document, resolved)

    receipt = {
        "schema_version": "masld-bench-scgpt-weight-state-correction-v1",
        "applied": bool(args.apply),
        "bundle": BUNDLE,
        "bundle_sha256_before": before,
        "bundle_changes": bundle_changes,
        "weight_license_unchanged": "UNDECLARED",
        "weights_deleted": False,
        "checkpoint_opened": False,
    }
    if args.apply:
        bundle_path.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        after = digest(bundle_path)
        pinner = json.loads(pinner_path.read_text(encoding="utf-8"))
        receipt["pinner_changes"] = correct_pinner(pinner, after)
        pinner_path.write_text(
            json.dumps(pinner, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        receipt["bundle_sha256_after"] = after
        receipt["pinner"] = PINNER

    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
