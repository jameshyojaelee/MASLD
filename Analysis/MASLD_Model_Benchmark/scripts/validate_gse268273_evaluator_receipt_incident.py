#!/usr/bin/env python3
"""Freeze the content-free GSE268273 evaluator-receipt incident record."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from masld_bench.artifacts import freeze_tree, write_json_exclusive


EXPECTED = {
    "path_opened": (
        "executions/model-data-074-21080832/fixture/evaluator_only/receipt.json"
    ),
    "timestamp": "2026-08-24T23:25:55-04:00",
    "turn": "/root/gse268273_ood_transfer:post_geneformer_queue_followup",
    "outcomes.tsv_not_opened": True,
    "labels_not_joined": True,
    "disposition": "evaluator_receipt_exposed",
}


class IncidentFirewallError(ValueError):
    """Raised when the content-free incident receipt differs."""


def run(*, source: Path, output: Path) -> None:
    if output.exists():
        raise IncidentFirewallError("refusing to overwrite incident artifact")
    value = json.loads(source.read_text(encoding="utf-8"))
    if value != EXPECTED or set(value) != set(EXPECTED):
        raise IncidentFirewallError("incident receipt fields or values differ")
    output.mkdir(mode=0o750)
    write_json_exclusive(output / "incident_receipt.json", value)
    freeze_tree(
        output,
        {
            "artifact_class": "gse268273_evaluator_receipt_exposure_incident",
            "disposition": "evaluator_receipt_exposed",
            "clean_or_sealed_champion_eligible": False,
            "independent_unexposed_rederivation_or_audit_required": True,
            "status": "contained_fail_closed",
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    run(source=arguments.source, output=arguments.output)
    print(json.dumps({"output": arguments.output.as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
