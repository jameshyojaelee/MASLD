#!/usr/bin/env python3
"""Correct the Geneformer checkpoint bundle and normalise its licence keys.

Two separate things are wrong, and only one of them is a rights matter.

The factual defect is ``weight_content_downloaded: false``.  All three
model.safetensors bodies are on disk at the frozen sizes with digests that
re-derive, and the bundle's own overlay recorded this staleness earlier today
and deliberately deferred it.  This closes that deferred item.

The other is a schema trap.  Geneformer is the only one of 104 bundles carrying
a singular ``license`` key instead of the ``weight_license`` / ``code_license``
pair.  Grepping for ``weight_license`` therefore found nothing and the absence
was reported as a rights gap, when ``license: "Apache-2.0"`` was present the
whole time.  Normalising the key removes the trap rather than just this
instance of it.

The licence VALUE is preserved exactly.  Apache-2.0 is a declared, permissive
grant: this correction must not deny terms that were granted, just as the scGPT
correction must not imply terms that were not.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


BUNDLE = "config/artifacts/models/geneformer/checkpoints.json"
PINNER = "config/artifacts/models/geneformer/native_task_disposition_20260825.json"

DECLARED_LICENCE = "Apache-2.0"

# Bodies whose presence makes the recorded false stale.
WEIGHT_BODIES = {
    "Geneformer-V1-10M/model.safetensors": 41183536,
    "Geneformer-V2-104M/model.safetensors": 417571156,
    "Geneformer-V2-316M/model.safetensors": 1265455076,
}
ACQUISITION = "executions/geneformer-acquisition-21014827/source"


class GeneformerCorrectionError(RuntimeError):
    """Raised when the tree does not match the audited pre-correction state."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(64 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify_bodies(root: Path, document: dict) -> int:
    """Confirm every recorded weight body exists and re-derives."""

    recorded = {a["path"]: a for a in document["artifacts"]}
    total = 0
    for relative, size in WEIGHT_BODIES.items():
        body = root / ACQUISITION / relative
        if not body.is_file():
            raise GeneformerCorrectionError(f"weight body missing: {body}")
        if body.stat().st_size != size:
            raise GeneformerCorrectionError(
                f"{relative} is {body.stat().st_size} bytes, expected {size}"
            )
        if digest(body) != recorded[relative]["sha256"]:
            raise GeneformerCorrectionError(f"{relative} digest differs from the bundle")
        total += size
    return total


def normalise_licence(document: dict) -> list[dict]:
    """Move the singular ``license`` onto the standard pair, value preserved.

    The code licence half is carried across too: the same Apache-2.0 grant
    covers the ``code_source_artifacts`` in the bundle, and dropping it would
    lose information the singular key was holding.
    """

    if "weight_license" in document:
        raise GeneformerCorrectionError("weight_license already present")
    value = document.get("license")
    if value != DECLARED_LICENCE:
        raise GeneformerCorrectionError(
            f"license is {value!r}, not the audited {DECLARED_LICENCE!r}"
        )
    document["weight_license"] = value
    document["code_license"] = value
    del document["license"]
    return [
        {"field": "license", "from": value, "to": "<removed, split into the standard pair>"},
        {"field": "weight_license", "from": "<key absent>", "to": value},
        {"field": "code_license", "from": "<key absent>", "to": value},
    ]


def correct_bundle(document: dict) -> list[dict]:
    if document["weight_content_downloaded"] is not False:
        raise GeneformerCorrectionError(
            "weight_content_downloaded is not the audited false; re-audit first"
        )
    document["weight_content_downloaded"] = True
    changes = [{"field": "weight_content_downloaded", "from": False, "to": True}]
    changes.extend(normalise_licence(document))

    # Guard the direction of travel. The scGPT correction must never upgrade a
    # licence; this one must never downgrade a declared, permissive one.
    if document["weight_license"] != DECLARED_LICENCE:
        raise GeneformerCorrectionError("weight_license was downgraded")
    return changes


def correct_pinner(document: dict, bundle_digest: str) -> list[dict]:
    changes: list[dict] = []
    for entry in document["bound_evidence"]:
        if entry["path"] == BUNDLE:
            changes.append({"field": "bound_evidence[checkpoints.json].sha256",
                            "from": entry["sha256"], "to": bundle_digest})
            entry["sha256"] = bundle_digest
    if not changes:
        raise GeneformerCorrectionError("pinner did not reference the bundle")
    for record in document.get("recorded_inconsistencies", []):
        if record.get("field", "").endswith("weight_content_downloaded"):
            record["action"] = "corrected_2026_08_25"
            record["correction"] = (
                "weight_content_downloaded set to true after re-deriving all "
                "three bodies from the bytes. This overlay recorded the "
                "staleness and deferred the fix; the deferred item is now closed."
            )
            changes.append({"field": "recorded_inconsistencies.action",
                            "from": "recorded_only_not_corrected_by_this_overlay",
                            "to": "corrected_2026_08_25"})
    return changes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        raise GeneformerCorrectionError("correction receipt output exists")

    bundle_path = args.root / BUNDLE
    pinner_path = args.root / PINNER
    document = json.loads(bundle_path.read_text(encoding="utf-8"))
    total_bytes = verify_bodies(args.root, document)
    before = digest(bundle_path)
    changes = correct_bundle(document)

    receipt = {
        "schema_version": "masld-bench-geneformer-weight-state-correction-v1",
        "applied": bool(args.apply),
        "bundle": BUNDLE,
        "bundle_sha256_before": before,
        "bundle_changes": changes,
        "weight_license_value_preserved": DECLARED_LICENCE,
        "rights_gap_found": False,
        "redistribution_prohibited": False,
        "champion_rights_blocked": False,
        "weight_bodies_verified": len(WEIGHT_BODIES),
        "weight_bytes_unique": total_bytes,
        "weights_deleted": False,
        "digests_corrected": 0,
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
