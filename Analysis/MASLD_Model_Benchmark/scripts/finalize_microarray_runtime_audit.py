#!/usr/bin/env python3
"""Combine read-only reader, annotation, and inductive-separation audits."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import tomllib


class MicroarrayRuntimeAuditError(RuntimeError):
    """Raised when the runtime audit cannot support its disposition."""


RUNTIME_PACKAGES = (
    "affy", "affyio", "affxparser", "oligo", "frma", "SCAN.UPC",
    "hgu133plus2cdf", "hgu133plus2.db", "hgu133plus2frmavecs",
    "pd.hugene.2.0.st", "hugene20sttranscriptcluster.db",
    "AnnotationDbi", "org.Hs.eg.db",
)
RUNTIME_CAPABILITIES = (
    "Calvin_header_reader_both_platforms", "GPL570_exact_frma", "GPL570_SCAN",
    "GPL570_affy_summary", "GPL16686_SCAN", "GPL16686_oligo_core_summary",
    "GPL16686_transcript_cluster_annotation", "GPL570_expression_rank_baseline",
    "GPL16686_expression_rank_baseline",
)
RUNTIME_FIXED_KEYS = (
    "schema_version", "status", "runtime_id", "R_version", "R_home",
    "library_paths", "reader.package", "reader.function_name",
    "reader.function_present", "reader.formals",
    *(f"package.{package}.installed" for package in RUNTIME_PACKAGES),
    *(f"package.{package}.version" for package in RUNTIME_PACKAGES),
    "reader.GPL570.status", "reader.GPL16686.status",
    "reader.GPL570.error", "reader.GPL16686.error",
    *(f"capability.{capability}" for capability in RUNTIME_CAPABILITIES),
    "raw_probe_grid_rank_as_gene_expression_allowed",
    "ordinary_all_sample_RMA_run", "quantile_normalization_run",
    "CEL_expression_values_exported", "model_training_activated", "labels_read",
)


BASELINE_RUNTIME_IDS = frozenset(
    {"module_R_4.4.3", "module_R_4.4.1", "micromamba_rnaseq"}
)


def require_runtime_inventory(runtime_ids: set[str]) -> None:
    """Require every baseline runtime; admit newly built runtimes alongside.

    The inventory is a superset check rather than an equality check so that an
    activated runtime can be added without dropping the three existing runtimes
    that establish what the pre-existing installations could and could not do.
    """

    if not BASELINE_RUNTIME_IDS.issubset(runtime_ids):
        raise MicroarrayRuntimeAuditError("R runtime inventory is incomplete")


def read_single_array_capability(path: Path | None) -> dict[str, object] | str:
    """Read the frozen single-array capability probe, when one is supplied."""

    if path is None:
        return "not_provided"
    probe = json.loads(path.read_text(encoding="utf-8"))
    capabilities = probe.get("capabilities", {})
    if (
        probe.get("status") != "pass_single_array_runtime_capability"
        or probe.get("all_sample_RMA_run") is not False
        or probe.get("across_array_quantile_normalization_run") is not False
        or probe.get("GEO_series_matrix_read") is not False
        or probe.get("expression_values_exported") is not False
        or probe.get("labels_read") is not False
        or probe.get("model_training_activated") is not False
        or probe.get("sealed_outcomes_read") is not False
        or capabilities.get("calvin_header_parse_both_platforms") is not True
        or capabilities.get("GPL570_single_array_summarization") is not True
        or capabilities.get("GPL16686_single_array_summarization") is not True
        or capabilities.get("probe_to_entrez_to_ensembl_both_platforms") is not True
    ):
        raise MicroarrayRuntimeAuditError("single-array capability probe differs")
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "status": probe["status"],
        "arrays_read_per_platform": probe["arrays_read_per_platform"],
        "capabilities": capabilities,
    }


def read_runtime_receipt(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("key", "value"):
            raise MicroarrayRuntimeAuditError("R runtime receipt fields differ")
        rows = [dict(row) for row in reader]
    result = {row["key"]: row["value"] for row in rows}
    if len(result) != len(rows):
        raise MicroarrayRuntimeAuditError("R runtime receipt keys are duplicated")
    if tuple(result) != RUNTIME_FIXED_KEYS:
        raise MicroarrayRuntimeAuditError("R runtime receipt fixed schema differs")
    if any(value == "" for value in result.values()):
        raise MicroarrayRuntimeAuditError(
            "R runtime receipt missing states must be explicit, not blank"
        )
    return result


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation", type=Path, required=True)
    parser.add_argument("--inductive-fixture", type=Path, required=True)
    parser.add_argument("--runtime-receipt", type=Path, action="append", required=True)
    parser.add_argument("--preprocessing-contract", type=Path, required=True)
    parser.add_argument("--label-contract", type=Path, required=True)
    parser.add_argument("--single-array-capability", type=Path, default=None)
    parser.add_argument("--dependency-install-performed", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    preparation = json.loads(
        (arguments.preparation / "runtime_preparation_audit.json").read_text(
            encoding="utf-8"
        )
    )
    inductive = json.loads(
        (arguments.inductive_fixture / "inductive_fixture_audit.json").read_text(
            encoding="utf-8"
        )
    )
    runtimes = [read_runtime_receipt(path) for path in arguments.runtime_receipt]
    require_runtime_inventory({value["runtime_id"] for value in runtimes})
    with arguments.preprocessing_contract.open("rb") as handle:
        preprocessing = tomllib.load(handle)
    with arguments.label_contract.open("rb") as handle:
        labels = tomllib.load(handle)
    if (
        preprocessing.get("ordinary_all_sample_RMA_allowed_for_scored_transfer") is not False
        or preprocessing.get("GEO_series_matrix_allowed_for_scored_transfer") is not False
        or inductive.get("status") != "pass_synthetic_held_array_isolation"
        or labels.get("schema_version") != "masld-bench-microarray-label-semantics-v1"
    ):
        raise MicroarrayRuntimeAuditError("inductive or label firewall differs")

    executable_readers = [
        runtime["runtime_id"]
        for runtime in runtimes
        if runtime.get("capability.Calvin_header_reader_both_platforms") == "true"
    ]
    gpl570_summarizers = [
        runtime["runtime_id"]
        for runtime in runtimes
        if runtime.get("capability.GPL570_expression_rank_baseline") == "true"
    ]
    gpl16686_summarizers = [
        runtime["runtime_id"]
        for runtime in runtimes
        if runtime.get("capability.GPL16686_expression_rank_baseline") == "true"
    ]
    capability = read_single_array_capability(arguments.single_array_capability)
    activated = isinstance(capability, dict)
    if gpl570_summarizers and activated:
        gpl570_status = "active_single_array_summarizer_with_frozen_platform_fixture"
    elif gpl570_summarizers:
        gpl570_status = "candidate_existing_runtime_requires_single_array_fixture"
    else:
        gpl570_status = "blocked_no_executable_single_array_summarizer"
    if gpl16686_summarizers and activated:
        gpl16686_status = "active_single_array_transcript_cluster_summarizer"
    elif gpl16686_summarizers:
        gpl16686_status = "candidate_existing_runtime_requires_single_array_fixture"
    else:
        gpl16686_status = "blocked_no_executable_single_array_transcript_cluster_summarizer"
    if preparation.get("gpl16686_to_gencode_v49_crosswalk_complete") is not False:
        raise MicroarrayRuntimeAuditError("GPL16686 crosswalk state differs")

    summary = {
        "schema_version": "masld-bench-microarray-runtime-audit-v1",
        "status": (
            "pass_runtime_inventory_with_single_array_summarizer_activated"
            if activated
            else "pass_runtime_inventory_with_biological_activation_blocked"
        ),
        "preparation_audit_sha256": sha256_file(
            arguments.preparation / "runtime_preparation_audit.json"
        ),
        "inductive_fixture_sha256": sha256_file(
            arguments.inductive_fixture / "inductive_fixture_audit.json"
        ),
        "preprocessing_contract_sha256": sha256_file(arguments.preprocessing_contract),
        "label_contract_sha256": sha256_file(arguments.label_contract),
        "runtime_receipts": {
            runtime["runtime_id"]: {
                "path": str(path),
                "sha256": sha256_file(path),
                "status": runtime["status"],
            }
            for runtime, path in zip(runtimes, arguments.runtime_receipt, strict=True)
        },
        "Calvin_header_reader_runtimes": sorted(executable_readers),
        "GPL570_expression_summarizer_runtimes": sorted(gpl570_summarizers),
        "GPL16686_expression_summarizer_runtimes": sorted(gpl16686_summarizers),
        "GPL570_biological_preprocessing_status": gpl570_status,
        "GPL16686_biological_preprocessing_status": gpl16686_status,
        "single_array_capability_probe": capability,
        "GPL16686_annotation_status": preparation["gpl16686_authoritative_mapping_state"],
        "GPL16686_to_GENCODE_v49_crosswalk_complete": False,
        "synthetic_training_only_quantile_fixture_passed": True,
        "synthetic_per_array_rank_fixture_passed": True,
        "synthetic_fixture_is_biological_expression_authority": False,
        "raw_probe_grid_rank_as_gene_expression_allowed": False,
        "GEO_series_matrix_allowed_for_overlap_audit_only": True,
        "GEO_series_matrix_allowed_for_model_input": False,
        "ordinary_all_sample_RMA_allowed": False,
        "dependency_install_performed": bool(arguments.dependency_install_performed),
        "normalization_run_on_biological_arrays": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "next_step": (
            "bind_the_activated_runtime_hash_into_the_microarray_TaskSpecs_and_fit_the_"
            "outer_training_only_preprocessing_object_before_any_held_array_is_scored"
            if activated
            else "obtain_and_hash_exact_authoritative_GPL16686_transcript_cluster_"
            "annotation_and_activate_only_an_existing_or_separately_approved_single_"
            "array_summarizer"
        ),
    }
    arguments.output.mkdir(parents=True, exist_ok=False)
    (arguments.output / "microarray_runtime_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
