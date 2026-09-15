#!/usr/bin/env python3
"""Freeze the read-only GSE105127 gzip-preflight-v2 campaign control."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re

from masld_bench.artifacts import (
    canonical_hash,
    freeze_tree,
    publish_directory_noreplace,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)
from masld_bench.hashing import sha256_file


PLAN_ARTIFACTS_SHA256 = "c4f91d2a3fedc7b3b1802a0d6cf3d1b1803ce1bea70cc343f6e18f744783b973"
REFERENCE_ARTIFACTS_SHA256 = "9575aa24dc39f855b961eecf53e29243b4a45fca6770c7337e05901bc7d1ec89"
RUNTIME_ARTIFACTS_SHA256 = "f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb"

LOCKED_PROJECT_FILES = (
    "config/campaigns/gse105127_gzip_preflight_revision_v2.toml",
    "scripts/build_gse105127_cpg_crosswalk.py",
    "scripts/build_gse105127_rsem_reference.py",
    "scripts/collapse_gse105127_rrbs_bundle.py",
    "scripts/consolidate_gse105127_rna_quantification.py",
    "scripts/finalize_gse105127_gzip_preflight_revision_v2.py",
    "scripts/finalize_gse105127_production_activation.py",
    "scripts/freeze_gse105127_gzip_stage_link_v2.py",
    "scripts/preflight_gse105127_gzip_inputs_v2.py",
    "scripts/prepare_gse105127_gzip_preflight_revision_v2.py",
    "scripts/quantify_gse105127_rna_bundle_v2.py",
    "src/masld_bench/artifacts.py",
    "src/masld_bench/hashing.py",
    "slurm/gse105127_gzip_v2_common.sh",
    "slurm/gse105127_v2_consolidate_cpu.sbatch",
    "slurm/gse105127_v2_crosswalk_io.sbatch",
    "slurm/gse105127_v2_finalize_cpu.sbatch",
    "slurm/gse105127_v2_prepare_cpu.sbatch",
    "slurm/gse105127_v2_rna_quant_cpu.sbatch",
    "slurm/gse105127_v2_rrbs_collapse_io.sbatch",
    "slurm/gse105127_v2_rsem_reference_cpu.sbatch",
    "tests/unit/test_preflight_gse105127_gzip_inputs_v2.py",
)
VALIDATED_PROJECT_FILES = frozenset(
    {
        "config/campaigns/gse105127_gzip_preflight_revision_v2.toml",
        "scripts/finalize_gse105127_gzip_preflight_revision_v2.py",
        "scripts/freeze_gse105127_gzip_stage_link_v2.py",
        "scripts/preflight_gse105127_gzip_inputs_v2.py",
        "scripts/prepare_gse105127_gzip_preflight_revision_v2.py",
        "scripts/quantify_gse105127_rna_bundle_v2.py",
        "src/masld_bench/artifacts.py",
        "src/masld_bench/hashing.py",
        "slurm/gse105127_gzip_v2_common.sh",
        "slurm/gse105127_v2_consolidate_cpu.sbatch",
        "slurm/gse105127_v2_crosswalk_io.sbatch",
        "slurm/gse105127_v2_finalize_cpu.sbatch",
        "slurm/gse105127_v2_prepare_cpu.sbatch",
        "slurm/gse105127_v2_rna_quant_cpu.sbatch",
        "slurm/gse105127_v2_rrbs_collapse_io.sbatch",
        "slurm/gse105127_v2_rsem_reference_cpu.sbatch",
        "tests/unit/test_preflight_gse105127_gzip_inputs_v2.py",
    }
)
AUDITED_REUSED_IMPLEMENTATION_SHA256 = {
    "scripts/build_gse105127_cpg_crosswalk.py": "2676dbeb4db530a90f0816d34df480649d66b50755b0978922d95f418d37736e",
    "scripts/build_gse105127_rsem_reference.py": "d360d7b4204bf9603c3796f9569b609ba12a3b1bf9d1091ca6fb1fe714234938",
    "scripts/collapse_gse105127_rrbs_bundle.py": "ee830a4b51d01df1b1a8beced6e7484ec92448e0d4f5304198776b6911e7d7e6",
    "scripts/consolidate_gse105127_rna_quantification.py": "cdb77fcdaded3f8ac73b93a7c98047f98c4a5a41abd11ff961e4e97c66e87d0f",
    "scripts/finalize_gse105127_production_activation.py": "fb891b529c545972c14f120c341256c993e3d54cbf13f8032499036d7b9d9ac3",
}


class GSE105127RevisionV2PreparationError(ValueError):
    """Raised when the revised campaign control is not reproducible."""


def require_artifacts(root: Path, expected_sha256: str, label: str) -> dict[str, object]:
    if sha256_file(root / "ARTIFACTS.json") != expected_sha256:
        raise GSE105127RevisionV2PreparationError(f"{label} artifact identity differs")
    return verify_frozen_tree(root)


def require_runtime_manifest(root: Path) -> dict[str, object]:
    """Admit the legacy runtime by its original checksum-manifest requirements."""

    artifacts = root / "ARTIFACTS.json"
    if sha256_file(artifacts) != RUNTIME_ARTIFACTS_SHA256:
        raise GSE105127RevisionV2PreparationError("runtime artifact identity differs")
    manifest = json.loads(artifacts.read_text(encoding="utf-8"))
    expected_prefix = str(root.resolve(strict=True) / "env")
    if (
        manifest.get("schema_version") != "masld-bench-transcriptformer-runtime-v1"
        or manifest.get("status") != "pass_cpu_runtime"
        or manifest.get("environment_prefix") != expected_prefix
        or manifest.get("python_no_user_site_required") is not True
        or not (root / "env/bin/python").is_file()
    ):
        raise GSE105127RevisionV2PreparationError("runtime manifest contract differs")
    declared = {
        "conda-explicit.txt": manifest.get("conda_explicit_sha256"),
        "pip-freeze.txt": manifest.get("pip_freeze_sha256"),
        "validation/import-probe.json": manifest.get("import_probe_sha256"),
        "wheels/transcriptformer-0.6.1-py3-none-any.whl": manifest.get("wheel_sha256"),
    }
    if any(
        not isinstance(expected, str) or sha256_file(root / relative) != expected
        for relative, expected in declared.items()
    ):
        raise GSE105127RevisionV2PreparationError("runtime declared member identity differs")
    checksums: dict[str, str] = {}
    for line in (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        fields = line.split(None, 1)
        if len(fields) != 2 or re.fullmatch(r"[0-9a-f]{64}", fields[0]) is None:
            raise GSE105127RevisionV2PreparationError("runtime SHA256SUMS is malformed")
        checksums[fields[1].lstrip("*")] = fields[0]
    if not checksums or any(sha256_file(root / relative) != expected for relative, expected in checksums.items()):
        raise GSE105127RevisionV2PreparationError("runtime SHA256SUMS member differs")
    import_probe = json.loads(
        (root / "validation/import-probe.json").read_text(encoding="utf-8")
    )
    if (
        import_probe.get("status") != "pass_cpu_import"
        or import_probe.get("versions", {}).get("transcriptformer") != "0.6.1"
        or import_probe.get("versions", {}).get("torch") != "2.5.1"
        or (root / "validation/pip-check.txt").read_text(encoding="utf-8").strip()
        != "No broken requirements found."
    ):
        raise GSE105127RevisionV2PreparationError("runtime CPU import contract differs")
    return manifest


def source_lock(project_root: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for relative in LOCKED_PROJECT_FILES:
        path = project_root / relative
        if not path.is_file():
            raise GSE105127RevisionV2PreparationError(f"locked project file is absent: {relative}")
        result[relative] = sha256_file(path)
    return result


def validation_source_lock(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split(None, 1)
        if len(fields) != 2 or re.fullmatch(r"[0-9a-f]{64}", fields[0]) is None:
            raise GSE105127RevisionV2PreparationError("validation source lock is malformed")
        result[fields[1].lstrip("*")] = fields[0]
    return result


def build_spec(
    project_root: Path,
    source_campaign_root: Path,
    independent_audit_root: Path,
    independent_audit_sha256: str,
    validation_root: Path,
    validation_sha256: str,
) -> dict[str, object]:
    require_artifacts(source_campaign_root / "plan", PLAN_ARTIFACTS_SHA256, "source plan")
    require_artifacts(source_campaign_root / "reference", REFERENCE_ARTIFACTS_SHA256, "legacy reference")
    require_runtime_manifest(
        project_root / "executions/environments/transcriptformer-torch2.5.1-21070738"
    )
    audit_manifest = require_artifacts(
        independent_audit_root, independent_audit_sha256, "independent source audit"
    )
    audit_receipt = json.loads(
        (independent_audit_root / "audit/independent_source_dag_audit.json").read_text(
            encoding="utf-8"
        )
    )
    metadata = audit_manifest.get("metadata", {})
    if (
        metadata.get("status") != "passed_audit_source_mechanics_verified_promotion_blocked"
        or metadata.get("compressed_legacy_source_usable_by_sequence_jobs") is not False
        or metadata.get("decompressed_legacy_reference_eligible_for_lift_mechanics") is not True
        or metadata.get("molecular_values_opened") is not False
        or metadata.get("modeling_or_scoring_authorized") is not False
        or metadata.get("strict_every_sequence_job_begins_with_gzip_t_contract_passed") is not False
        or metadata.get("current_dag_promotion_eligible") is not False
        or audit_receipt.get("rna_or_methylation_values_parsed") is not False
        or audit_receipt.get("labels_read") is not False
        or audit_receipt.get("outcomes_read") is not False
        or audit_receipt.get("predictions_read") is not False
        or audit_receipt.get("fit_or_score_performed") is not False
        or audit_receipt.get("clean_or_sealed_champion_eligible_by_this_receipt") is not False
    ):
        raise GSE105127RevisionV2PreparationError("independent source audit contract differs")
    validation_manifest = require_artifacts(
        validation_root, validation_sha256, "gzip-preflight-v2 validation"
    )
    validation_receipt = json.loads(
        (validation_root / "receipt.json").read_text(encoding="utf-8")
    )
    if (
        validation_manifest.get("metadata", {}).get("status") != "passed"
        or validation_manifest.get("metadata", {}).get("molecular_values_opened") is not False
        or validation_receipt.get("status") != "passed"
        or validation_receipt.get("synthetic_tests") != 9
        or validation_receipt.get("truncated_gzip_rejected") is not True
        or validation_receipt.get("trailing_garbage_rejected") is not True
        or validation_receipt.get("additional_member_rejected") is not True
        or validation_receipt.get("wrong_digest_rejected") is not True
        or validation_receipt.get("same_job_wrapper_order_verified") is not True
        or validation_receipt.get("quantifier_v2_receipt_hash_regression_passed") is not True
        or validation_receipt.get("molecular_values_opened") is not False
        or validation_receipt.get("labels_accessed") is not False
        or validation_receipt.get("outcomes_accessed") is not False
        or validation_receipt.get("fit_or_score_performed") is not False
    ):
        raise GSE105127RevisionV2PreparationError("gzip-preflight-v2 validation contract differs")
    locks = source_lock(project_root)
    validated_locks = validation_source_lock(validation_root / "source.sha256")
    if any(validated_locks.get(relative) != locks[relative] for relative in VALIDATED_PROJECT_FILES):
        raise GSE105127RevisionV2PreparationError("validated v2 project source differs")
    if any(
        locks.get(relative) != expected
        for relative, expected in AUDITED_REUSED_IMPLEMENTATION_SHA256.items()
    ):
        raise GSE105127RevisionV2PreparationError("independently audited reused implementation differs")
    return {
        "schema_version": "masld-bench-gse105127-gzip-preflight-revision-v2-control",
        "status": "source_audited_revision_locked_before_submission",
        "source_campaign_root": str(source_campaign_root.resolve(strict=True)),
        "source_plan_artifacts_sha256": PLAN_ARTIFACTS_SHA256,
        "source_reference_artifacts_sha256": REFERENCE_ARTIFACTS_SHA256,
        "independent_source_audit_root": str(independent_audit_root.resolve(strict=True)),
        "independent_source_audit_artifacts_sha256": independent_audit_sha256,
        "validation_root": str(validation_root.resolve(strict=True)),
        "validation_artifacts_sha256": validation_sha256,
        "runtime_artifacts_sha256": RUNTIME_ARTIFACTS_SHA256,
        "runtime_admission": "legacy_checksum_manifest_nonfrozen_environment",
        "project_file_sha256": locks,
        "participants": 19,
        "participant_zone_rows": 57,
        "bundles": 8,
        "pairing_topology": "adjacent_section",
        "reused_artifacts": [
            "frozen_label_free_plan",
            "independently_audited_raw_RNA_FASTQ",
            "source_plan_and_prior_collapse_receipt_audited_raw_RRBS_BED",
            "independently_audited_exact_decompressed_GRCh37_reference",
            "independently_audited_lift_chains",
        ],
        "development_only_artifacts_not_reused": [
            "prior_RSEM_reference",
            "prior_RRBS_collapses",
            "prior_RNA_quantifications",
            "prior_CpG_crosswalk",
            "prior_consolidation",
            "prior_activation",
        ],
        "required_job_local_preflights": [
            "rsem_reference",
            "rrbs_collapse",
            "rna_quantification",
            "cpg_crosswalk",
        ],
        "molecular_values_opened_during_control_preparation": False,
        "labels_accessed": False,
        "outcomes_accessed": False,
        "fit_or_score_performed": False,
    }


def revision_id(spec: dict[str, object]) -> str:
    return canonical_hash(spec)[:16]


def prepare(
    *,
    project_root: Path,
    source_campaign_root: Path,
    independent_audit_root: Path,
    independent_audit_sha256: str,
    validation_root: Path,
    validation_sha256: str,
    output_base: Path,
) -> dict[str, object]:
    if any(
        re.fullmatch(r"[0-9a-f]{64}", value) is None
        for value in (independent_audit_sha256, validation_sha256)
    ):
        raise GSE105127RevisionV2PreparationError("audit or validation SHA-256 is invalid")
    spec = build_spec(
        project_root.resolve(strict=True),
        source_campaign_root.resolve(strict=True),
        independent_audit_root.resolve(strict=True),
        independent_audit_sha256,
        validation_root.resolve(strict=True),
        validation_sha256,
    )
    identifier = revision_id(spec)
    spec["revision_id"] = identifier
    campaign = output_base.resolve(strict=True) / f"gse105127-gzip-preflight-v2-{identifier}"
    if campaign.exists():
        manifest = verify_frozen_tree(campaign / "control")
        observed = json.loads((campaign / "control/control.json").read_text(encoding="utf-8"))
        if observed != spec or manifest.get("metadata", {}).get("revision_id") != identifier:
            raise GSE105127RevisionV2PreparationError("existing revision control differs")
        return {"campaign_root": str(campaign), "revision_id": identifier, "status": "verified_replay"}
    job = os.environ.get("SLURM_JOB_ID", "nojob")
    stage = output_base / f".{campaign.name}.control-{job}.staging"
    if stage.exists():
        raise GSE105127RevisionV2PreparationError("control staging path already exists")
    control = stage / "control"
    control.mkdir(parents=True, mode=0o750)
    write_json_exclusive(control / "control.json", spec)
    lines = [f"{digest}  {project_root / relative}" for relative, digest in sorted(spec["project_file_sha256"].items())]
    lines.extend(
        (
            f"{independent_audit_sha256}  {independent_audit_root / 'ARTIFACTS.json'}",
            f"{validation_sha256}  {validation_root / 'ARTIFACTS.json'}",
            f"{PLAN_ARTIFACTS_SHA256}  {source_campaign_root / 'plan/ARTIFACTS.json'}",
            f"{REFERENCE_ARTIFACTS_SHA256}  {source_campaign_root / 'reference/ARTIFACTS.json'}",
            f"{RUNTIME_ARTIFACTS_SHA256}  {project_root / 'executions/environments/transcriptformer-torch2.5.1-21070738/ARTIFACTS.json'}",
        )
    )
    write_text_exclusive(control / "source-lock.sha256", "\n".join(lines) + "\n")
    freeze_tree(
        control,
        {
            "artifact_class": "gse105127_gzip_preflight_revision_v2_control",
            "revision_id": identifier,
            "status": "source_audited_revision_locked_before_submission",
            "labels_accessed": False,
            "fit_or_score_performed": False,
        },
    )
    verify_frozen_tree(control)
    publish_directory_noreplace(stage, campaign)
    verify_frozen_tree(campaign / "control")
    return {"campaign_root": str(campaign), "revision_id": identifier, "status": "created"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--source-campaign-root", required=True, type=Path)
    parser.add_argument("--independent-audit-root", required=True, type=Path)
    parser.add_argument("--independent-audit-sha256", required=True)
    parser.add_argument("--validation-root", required=True, type=Path)
    parser.add_argument("--validation-sha256", required=True)
    parser.add_argument("--output-base", required=True, type=Path)
    parser.add_argument("--revision-id-only", action="store_true")
    args = parser.parse_args()
    if args.revision_id_only:
        spec = build_spec(
            args.project_root.resolve(strict=True),
            args.source_campaign_root.resolve(strict=True),
            args.independent_audit_root.resolve(strict=True),
            args.independent_audit_sha256,
            args.validation_root.resolve(strict=True),
            args.validation_sha256,
        )
        print(revision_id(spec))
    else:
        print(
            json.dumps(
                prepare(
                    project_root=args.project_root,
                    source_campaign_root=args.source_campaign_root,
                    independent_audit_root=args.independent_audit_root,
                    independent_audit_sha256=args.independent_audit_sha256,
                    validation_root=args.validation_root,
                    validation_sha256=args.validation_sha256,
                    output_base=args.output_base,
                ),
                sort_keys=True,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
