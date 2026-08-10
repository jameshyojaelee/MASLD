#!/usr/bin/env python3
"""Build a route-specific editable, outcome-free Plan 46 response workspace."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import shutil
from pathlib import Path

from clinical_response_common import (
    CENSUS_METRICS, DEFAULT_OUTREACH_ROOT, ROUTE_ROLES, SCRIPT_ROOT,
    candidate_root, sha256_file, validate_outreach, write_json, write_tsv,
)


ROUTE_TABLES = (
    "clinical_cohort_feasibility_questions.tsv",
    "assay_availability_matrix.tsv",
    "clinical_role_firewall.tsv",
    "minimum_data_dictionary.tsv",
    "route_specific_questions.tsv",
    "public_inventory_audit.tsv",
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--route-id", required=True, choices=sorted(ROUTE_ROLES))
    parser.add_argument("--candidate-id", required=True)
    parser.add_argument("--outreach-root", default=str(DEFAULT_OUTREACH_ROOT))
    return parser.parse_args()


def main() -> None:
    args = arguments()
    root = candidate_root(args.candidate_id)
    if root.exists():
        raise RuntimeError(f"Refusing to overwrite clinical response scaffold: {root}")
    outreach = Path(args.outreach_root).resolve()
    seal_path, _, routes = validate_outreach(outreach)
    route = routes[args.route_id]
    source_root = outreach / "routes" / args.route_id
    root.mkdir(parents=True)
    outputs: list[Path] = []
    for name in ROUTE_TABLES:
        destination = root / name
        shutil.copyfile(source_root / name, destination)
        outputs.append(destination)

    identity = root / "clinical_response_identity.tsv"
    write_tsv(identity, [{
        "response_package_id": root.name,
        "parent_outreach_candidate_id": outreach.name,
        "parent_outreach_seal_sha256": sha256_file(seal_path),
        "route_id": args.route_id,
        "candidate_cohort_role": route["candidate_cohort_role"],
        "independence_group": route["independence_group"],
        "organization": route["organization"],
        "cohort_uid": "",
        "institution_uid": "",
        "trial_uid": "",
        "named_custodian": "",
        "governance_route": "",
        "credible_completion_utc": "",
        "prepared_by": "",
        "submitted_utc": "",
        "single_missing_gate_recoverable": "no",
        "missing_gate_id": "",
        "recovery_due_utc": "",
        "response_status": "editable_draft",
    }], [
        "response_package_id", "parent_outreach_candidate_id",
        "parent_outreach_seal_sha256", "route_id", "candidate_cohort_role",
        "independence_group", "organization", "cohort_uid", "institution_uid",
        "trial_uid", "named_custodian", "governance_route",
        "credible_completion_utc", "prepared_by", "submitted_utc",
        "single_missing_gate_recoverable", "missing_gate_id",
        "recovery_due_utc", "response_status",
    ])
    outputs.append(identity)
    census = root / "blinded_cohort_census.tsv"
    write_tsv(census, [
        {"metric": metric, "value": "", "evidence_path": ""}
        for metric in CENSUS_METRICS
    ], ["metric", "value", "evidence_path"])
    outputs.append(census)
    evidence = root / "evidence_file_register.tsv"
    write_tsv(evidence, [], ["evidence_id", "relative_path", "description", "sha256", "size_bytes"])
    outputs.append(evidence)
    readme = root / "EDITABLE_CLINICAL_RESPONSE_README.md"
    readme.write_text(
        """# Editable metadata-only clinical response

This workspace is not a scientific release. Fill every CF and route-specific
answer with `yes`, `no`, or `unresolved`; provide evidence, owner, and UTC
attestation; complete the aggregate blinded census and assay matrix; name the
custodian and governance route; and register every evidence file under
`evidence/`. Use no participant rows and include no expression, protein,
spatial, treatment-linked, responder-linked, or other molecular outcome.

Set `response_status` to `submitted` only when complete. A submitted response
still authorizes no participant-data request, outcome access, analysis, or
paper promotion. It is adjudicated into a separate immutable candidate.
""",
        encoding="utf-8",
    )
    outputs.append(readme)
    marker = {
        "status": "editable_metadata_only_clinical_response",
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "route_id": args.route_id,
        "candidate_cohort_role": route["candidate_cohort_role"],
        "independence_group": route["independence_group"],
        "parent_outreach_candidate_id": outreach.name,
        "parent_outreach_seal_sha256": sha256_file(seal_path),
        "mutable_response_workspace": True,
        "participant_rows_present": False,
        "clinical_outcomes_present": False,
        "molecular_outcomes_present": False,
        "permission_assumed": False,
        "analysis_authorized": False,
        "paper_promotion_authorized": False,
        "producer_sha256": sha256_file(SCRIPT_ROOT / "57_build_clinical_response_scaffold.py"),
        "initial_template_sha256": {path.name: sha256_file(path) for path in outputs},
    }
    write_json(root / "DRAFT_EDITABLE.json", marker)
    print(json.dumps(marker, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

