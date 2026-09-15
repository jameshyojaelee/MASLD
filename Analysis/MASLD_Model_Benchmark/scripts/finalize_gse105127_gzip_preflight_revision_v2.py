#!/usr/bin/env python3
"""Freeze terminal evidence that every revised GSE105127 gzip stage passed."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


class GSE105127RevisionV2FinalizationError(ValueError):
    """Raised when the terminal v2 requirements are incomplete."""


def load_receipt(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    manifest = verify_frozen_tree(root)
    receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
    return manifest, receipt


def finalize(
    *,
    control_root: Path,
    activation_root: Path,
    stage_links_root: Path,
    independent_audit_root: Path,
    output: Path,
) -> dict[str, object]:
    if output.exists():
        raise GSE105127RevisionV2FinalizationError("refusing to overwrite terminal revision receipt")
    control_manifest = verify_frozen_tree(control_root)
    control = json.loads((control_root / "control.json").read_text(encoding="utf-8"))
    audit_manifest = verify_frozen_tree(independent_audit_root)
    _, activation = load_receipt(activation_root)
    if (
        control.get("status") != "source_audited_revision_locked_before_submission"
        or control.get("independent_source_audit_root") != str(independent_audit_root.resolve(strict=True))
        or control.get("independent_source_audit_artifacts_sha256")
        != sha256_file(independent_audit_root / "ARTIFACTS.json")
        or audit_manifest.get("metadata", {}).get("status")
        != "passed_audit_source_mechanics_verified_promotion_blocked"
        or activation.get("status") != "assays_ready_no_fit_evaluator_implementation_pending"
        or activation.get("labels_accessed") is not False
        or activation.get("model_fit_run") is not False
        or activation.get("model_score_run") is not False
    ):
        raise GSE105127RevisionV2FinalizationError("control, audit, or activation contract differs")
    expected = [("rsem_reference", None), ("cpg_crosswalk", None)]
    expected.extend(("rrbs_collapse", bundle) for bundle in range(8))
    expected.extend(("rna_quantification", bundle) for bundle in range(8))
    links: dict[str, dict[str, object]] = {}
    file_totals: dict[str, int] = {}
    job_ids: set[str] = set()
    for stage, bundle_id in expected:
        name = stage if bundle_id is None else f"{stage}_bundle_{bundle_id:02d}"
        _, link = load_receipt(stage_links_root / name)
        if (
            link.get("status") != "job_local_preflight_bound_to_completed_stage"
            or link.get("stage") != stage
            or link.get("bundle_id") != bundle_id
            or link.get("slurm_job_id") in (None, "")
            or link.get("molecular_values_opened_by_linker") is not False
            or link.get("labels_accessed") is not False
            or link.get("outcomes_accessed") is not False
            or link.get("fit_or_score_performed") is not False
        ):
            raise GSE105127RevisionV2FinalizationError("stage-link contract differs")
        links[name] = {
            "artifacts_sha256": sha256_file(stage_links_root / name / "ARTIFACTS.json"),
            "slurm_job_id": link["slurm_job_id"],
            "preflight_file_count": link["preflight_file_count"],
        }
        file_totals[stage] = file_totals.get(stage, 0) + int(link["preflight_file_count"])
        job_ids.add(str(link["slurm_job_id"]))
    if file_totals != {
        "rsem_reference": 2,
        "cpg_crosswalk": 59,
        "rrbs_collapse": 57,
        "rna_quantification": 57,
    }:
        raise GSE105127RevisionV2FinalizationError("preflight file census differs")
    receipt = {
        "schema_version": "masld-bench-gse105127-gzip-preflight-revision-v2-terminal",
        "status": "assays_ready_no_fit_gzip_contract_repaired",
        "revision_id": control["revision_id"],
        "participants": 19,
        "participant_zone_rows": 57,
        "pairing_topology": "adjacent_section",
        "job_local_gzip_preflight_stages": [
            "rsem_reference",
            "rrbs_collapse",
            "rna_quantification",
            "cpg_crosswalk",
        ],
        "preflight_file_totals": file_totals,
        "stage_link_count": len(links),
        "distinct_stage_job_count": len(job_ids),
        "stage_links": links,
        "control_artifacts_sha256": sha256_file(control_root / "ARTIFACTS.json"),
        "independent_source_audit_artifacts_sha256": sha256_file(
            independent_audit_root / "ARTIFACTS.json"
        ),
        "activation_artifacts_sha256": sha256_file(activation_root / "ARTIFACTS.json"),
        "prior_campaign_and_results_remain_development_only": True,
        "new_campaign_is_assay_ready_development_evidence_only": True,
        "clean_or_sealed_champion_eligible_by_this_receipt": False,
        "molecular_values_opened_by_terminal_validator": False,
        "labels_accessed": False,
        "outcomes_accessed": False,
        "fit_or_score_performed": False,
    }
    output.mkdir(mode=0o750)
    write_json_exclusive(output / "receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse105127_gzip_preflight_revision_v2_terminal",
            "revision_id": control["revision_id"],
            "status": receipt["status"],
            "molecular_values_opened": False,
            "labels_accessed": False,
            "fit_or_score_performed": False,
        },
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control-root", required=True, type=Path)
    parser.add_argument("--activation-root", required=True, type=Path)
    parser.add_argument("--stage-links-root", required=True, type=Path)
    parser.add_argument("--independent-audit-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(finalize(**vars(args)), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
