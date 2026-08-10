#!/usr/bin/env python3
"""Classify the frozen active consumer inventory without running downstream work."""

from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
import json
from pathlib import Path


ALLOWED_STATUS = {"rebuilt", "verified_invariant", "retired", "out_of_scope"}
REQUIRED_CANONICAL_CHAIN = (
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/07_consensus_degs.R",
    "RNA-seq/27a_assemble_evidence_atlas.R",
    "RNA-seq/45a_integrate_all_sources.R",
    "RNA-seq/75_integrate_causal_overhaul.R",
    "RNA-seq/217_stratified_causal_atlas.R",
    "RNA-seq/46d_convergence_evidence.R",
    "RNA-seq/46e_fisher_combined_test.R",
    "RNA-seq/27b_benchmark_presets.R",
    "scripts/portal/generate_convergence_evidence_json.py",
)
IMMEDIATE_POST_ACCEPTANCE = {
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/07_consensus_degs.R",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_frozen_project_root(project: Path, root: Path) -> None:
    expected = (root / "source_snapshot").resolve(strict=True)
    if project != expected:
        raise SystemExit("--project-root must be the candidate's frozen source_snapshot")


def classify_verdict_state(verdict: dict) -> tuple[str, str, bool, bool, bool, bool]:
    overall_status = verdict.get("overall_status")
    completion_status = verdict.get("completion_status")
    escalation_required = verdict.get("escalation_required")
    accepted_stable = (
        overall_status == "STABLE"
        and completion_status == "PRIMARY_GATES_COMPLETE"
        and escalation_required is False
    )
    material = (
        overall_status == "MATERIAL"
        and completion_status == "MATERIAL_TARGETED_ESCALATION_REQUIRED"
        and escalation_required is True
    )
    pending = (
        overall_status == "PENDING_TARGETED_ESCALATION"
        and completion_status == "PENDING_TARGETED_ESCALATION"
        and escalation_required is True
    )
    if sum((accepted_stable, material, pending)) != 1:
        raise SystemExit("Verdict acceptance/completion/escalation state is inconsistent")
    return (
        overall_status, completion_status, escalation_required,
        accepted_stable, material, pending,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args()
    project = args.project_root.resolve(strict=True)
    root = args.run_root.resolve(strict=True)
    if not (root / ".bg001_candidate_root").is_file():
        raise SystemExit("Missing candidate sentinel")
    require_frozen_project_root(project, root)
    verdict_path = root / "comparisons/verdict.json"
    verdict = json.loads(verdict_path.read_text())
    (
        overall_status, completion_status, escalation_required,
        accepted_stable, material, pending,
    ) = classify_verdict_state(verdict)

    inventory_path = root / "contract/dependency_source_inventory.tsv"
    with inventory_path.open(newline="") as handle:
        inventory = list(csv.DictReader(handle, delimiter="\t"))
    retirement_path = (
        project
        / "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/dependency_retirement_registry.tsv"
    )
    with retirement_path.open(newline="") as handle:
        retirement = list(csv.DictReader(handle, delimiter="\t"))
    if not retirement or any(not row.get("path_glob") or not row.get("reason") for row in retirement):
        raise SystemExit("Explicit dependency retirement registry is empty or malformed")
    rows = []
    for item in inventory:
        snapshot = project / item["snapshot_relative_path"]
        observed_hash = sha256(snapshot)
        if observed_hash != item["source_sha256"]:
            raise SystemExit(f"Frozen consumer source hash drift: {item['consumer_path']}")
        roles = set(item["consumed_roles"].split("|"))
        retirement_hits = [
            row for row in retirement
            if fnmatch.fnmatch(item["consumer_path"], row["path_glob"])
        ]
        retired = bool(retirement_hits)
        if retired:
            status = "retired"
            proposed = "none"
            reason = "|".join(sorted({row["reason"] for row in retirement_hits}))
        else:
            status = "out_of_scope"
            direct = bool(roles & {
                "raw_counts", "qc", "merged_raw", "dge", "canonical_deg", "consensus",
                "multi_evidence_atlas", "convergence_evidence",
            })
            if direct or material:
                proposed = "rebuild"
            elif accepted_stable:
                proposed = "verify_invariant"
            else:
                proposed = "hold_pending_targeted_escalation"
            reason = "candidate_execution_stops_before_approved_downstream_rebuild"
        rows.append(
            {
                "record_type": "consumer",
                "consumer_path": item["consumer_path"],
                "snapshot_relative_path": item["snapshot_relative_path"],
                "consumed_roles": item["consumed_roles"],
                "direct_roles": item.get("direct_roles", ""),
                "direct_access_modes": item.get("direct_access_modes", ""),
                "effective_access_modes": item.get("effective_access_modes", ""),
                "wrapper_targets": item.get("wrapper_targets", ""),
                "evidence_lines": item.get("evidence_lines", ""),
                "source_kind": item.get("source_kind", ""),
                "source_sha256": observed_hash,
                "status": status,
                "proposed_disposition": proposed,
                "reason": reason,
            }
        )

    consumer_rows = {row["consumer_path"]: row for row in rows if row["record_type"] == "consumer"}
    missing_chain = [path for path in REQUIRED_CANONICAL_CHAIN if path not in consumer_rows]
    if missing_chain:
        raise SystemExit(f"Frozen dependency DAG omitted required canonical chain scripts: {missing_chain}")
    for ordinal, path in enumerate(REQUIRED_CANONICAL_CHAIN, start=1):
        row = consumer_rows[path]
        if path in IMMEDIATE_POST_ACCEPTANCE:
            row["proposed_disposition"] = "rebuild"
            row["reason"] = f"required_immediate_post_acceptance_step_{ordinal:02d}"
        else:
            row["proposed_disposition"] = "defer_for_bg031"
            row["reason"] = f"required_bg031_synchronized_chain_step_{ordinal:02d}"

    # Explicit chain nodes prevent a text-search inventory from obscuring the
    # required post-acceptance topology.
    chain_rows = (
        ("candidate::corrected_counts", "raw_counts", "rebuilt", "retain_candidate", "five validated fragment matrices"),
        ("candidate::qc_and_dge", "qc|merged_raw|dge", "rebuilt", "retain_candidate", "four isolated arm substrates"),
        ("candidate::05h_deg", "canonical_deg", "rebuilt", "retain_candidate", "F_five joint candidate; not canonical"),
        ("chain::07_consensus", "canonical_deg|consensus", "out_of_scope", "rebuild", "always rebuild after candidate acceptance"),
        ("chain::atlas_bulk_fields", "canonical_deg|consensus", "out_of_scope", "defer_for_bg031", "mandatory rebuild in the BG-031-synchronized combined-atlas chain"),
        ("chain::continuous_logfc_rankings", "canonical_deg", "out_of_scope", "rebuild", "coefficient-sensitive consumers"),
        ("chain::combined_atlas_release", "canonical_deg|consensus", "out_of_scope", "defer_for_bg031", "requires accepted frozen BG-031 genetic snapshot"),
    )
    for path, roles, status, proposed, reason in chain_rows:
        rows.append(
            {
                "record_type": "chain_node",
                "consumer_path": path,
                "snapshot_relative_path": "",
                "consumed_roles": roles,
                "direct_roles": roles,
                "direct_access_modes": "",
                "effective_access_modes": "contract",
                "wrapper_targets": "",
                "evidence_lines": "explicit required BG-001 post-acceptance chain node",
                "source_kind": "contract",
                "source_sha256": "",
                "status": status,
                "proposed_disposition": proposed,
                "reason": reason,
            }
        )
    if any(row["status"] not in ALLOWED_STATUS for row in rows):
        raise SystemExit("Dependency classifier emitted an invalid status")
    rows.sort(key=lambda row: (row["record_type"], row["consumer_path"]))

    output = root / "comparisons/dependency_manifest.tsv"
    with output.open("w", newline="") as handle:
        fields = (
            "record_type", "consumer_path", "snapshot_relative_path", "consumed_roles",
            "direct_roles", "direct_access_modes", "effective_access_modes", "wrapper_targets",
            "evidence_lines", "source_kind",
            "source_sha256", "status", "proposed_disposition", "reason",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    counts = {status: sum(row["status"] == status for row in rows) for status in sorted(ALLOWED_STATUS)}
    disposition = [
        "# BG-001 proposed downstream disposition",
        "",
        f"Candidate verdict: **{overall_status}** ({completion_status}).",
        "",
        "No canonical or downstream artifact was promoted or rebuilt in this run.",
        "Candidate acceptance is eligible only when completion_status is PRIMARY_GATES_COMPLETE and escalation_required is false.",
        "After explicit candidate acceptance, rebuild 07 consensus and independent continuous-logFC consumers whose actual fields changed.",
        "Rebuild atlas bulk-expression fields as part of the combined atlas/evidence-class/convergence/release chain only after BG-031 supplies an accepted frozen genetic snapshot.",
        "If the candidate is MATERIAL, expand rebuilding to every active transitive consumer and perform claim-level review.",
        "",
        "Status counts: " + ", ".join(f"{key}={value}" for key, value in counts.items()),
    ]
    (root / "comparisons/downstream_disposition.md").write_text("\n".join(disposition) + "\n")

    verdict["dependency_manifest_sha256"] = sha256(output)
    verdict_path.write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {output} with {len(rows)} frozen consumers/chain nodes; no downstream promotion performed")


if __name__ == "__main__":
    main()
