#!/usr/bin/env python3
"""Verify every read-only source copied into the Resource-paper candidate."""

from __future__ import annotations

import csv
import json
from pathlib import Path, PurePosixPath

from accepted_run_provenance import (
    EXPECTED_GRAPH_BYTES,
    EXPECTED_GRAPH_FILES,
    verify_accepted_graph,
)
from resource_contract import (
    ACCEPTED_RUN_INPUT,
    BG_RUN_ID,
    CANONICAL_REFERENCE,
    CANONICAL_REFERENCE_SHA256,
    CANDIDATE_ID,
    CANDIDATE_ROOT,
    GENCODE_GTF,
    GENCODE_GTF_BYTES,
    GENCODE_GTF_SHA256,
    GSE193066_METADATA_SHA256,
    REVIEWED_CODE_FILES,
    ContractError,
    require,
    require_regular_file,
    sha256,
)
from snapshot_io import parse_sha256_manifest


EXPECTED_CENSUS = [
    {
        "layer": "matched_raw_nine_cohort",
        "n_samples": "1281",
        "n_genes": "86369",
        "role": "resource_source",
    },
    {
        "layer": "pass_technical_nine_cohort",
        "n_samples": "1257",
        "n_genes": "24196",
        "role": "stage_context_dge",
    },
    {
        "layer": "pooled_five_cohort",
        "n_samples": "844",
        "n_genes": "23370",
        "role": "disease_control_model",
    },
]
EXPECTED_OBJECT_CENSUS = {
    "raw_n_genes": 86369,
    "raw_n_samples": 1281,
    "legacy_n_genes": 24196,
    "legacy_n_samples": 1257,
    "legacy_n_cohorts": 9,
    "five_n_genes": 23370,
    "five_n_samples": 844,
    "five_n_cohorts": 5,
}


def load_json(path: Path) -> dict[str, object]:
    require_regular_file(path)
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def load_tsv(
    path: Path,
    expected_fields: list[str] | None = None,
) -> tuple[list[str], list[dict[str, str]]]:
    require_regular_file(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if expected_fields is not None:
        require(fields == expected_fields, f"TSV schema drift: {path}")
    require(
        all(None not in row and all(value is not None for value in row.values()) for row in rows),
        f"malformed TSV row width: {path}",
    )
    return fields, rows


def _safe_relative(raw: str, label: str) -> PurePosixPath:
    path = PurePosixPath(raw)
    require(
        not path.is_absolute() and ".." not in path.parts and path.as_posix() not in {"", "."},
        f"unsafe {label}: {raw}",
    )
    return path


def verify() -> dict[str, object]:
    require(CANDIDATE_ROOT.name == CANDIDATE_ID, "candidate identity drift")
    require(CANDIDATE_ROOT.is_dir(), "candidate root is missing")
    require(not CANDIDATE_ROOT.is_symlink(), "candidate root may not be a symlink")
    sentinel = CANDIDATE_ROOT / ".resource_candidate_root"
    require_regular_file(sentinel)
    require(
        sentinel.read_text(encoding="utf-8").strip() == CANDIDATE_ID,
        "sentinel drift",
    )
    for path in CANDIDATE_ROOT.rglob("*"):
        require(not path.is_symlink(), f"candidate contains a symlink: {path}")

    source_manifest = CANDIDATE_ROOT / "manifests/source_selection.tsv"
    _, rows = load_tsv(
        source_manifest,
        [
            "source_scope",
            "source_path",
            "candidate_path",
            "role",
            "size_bytes",
            "sha256",
        ],
    )
    require(rows, "source-selection manifest is empty")
    observed: dict[str, dict[str, str]] = {}
    root_resolved = CANDIDATE_ROOT.resolve(strict=True)
    for row in rows:
        relative = _safe_relative(row["candidate_path"], "candidate path")
        relative_text = relative.as_posix()
        require(relative_text not in observed, f"duplicate candidate path: {relative_text}")
        scope = row["source_scope"]
        source_text = row["source_path"]
        if scope == "project_relative":
            _safe_relative(source_text, "project source path")
        elif scope == "external_absolute":
            source_path = Path(source_text)
            require(
                source_path.is_absolute() and ".." not in source_path.parts,
                f"unsafe external source path: {source_text}",
            )
        else:
            raise ContractError(f"unknown source scope: {scope}")
        candidate_path = CANDIDATE_ROOT.joinpath(*relative.parts)
        require_regular_file(candidate_path)
        require(
            candidate_path.resolve(strict=True).is_relative_to(root_resolved),
            f"candidate path escapes root: {relative_text}",
        )
        require(
            candidate_path.stat().st_size == int(row["size_bytes"]),
            f"byte drift: {relative_text}",
        )
        require(
            sha256(candidate_path) == row["sha256"],
            f"hash drift: {relative_text}",
        )
        role_tokens = row["role"].split(";")
        require(
            role_tokens == sorted(set(role_tokens)) and all(role_tokens),
            f"malformed source-selection roles: {relative_text}",
        )
        observed[relative_text] = row

    graph = verify_accepted_graph(ACCEPTED_RUN_INPUT)
    require(len(graph) == EXPECTED_GRAPH_FILES, "mirrored graph file-count drift")
    require(
        sum(size for _, size, _ in graph.values()) == EXPECTED_GRAPH_BYTES,
        "mirrored graph byte-count drift",
    )
    for relative, (role, size, digest) in graph.items():
        candidate_relative = f"inputs/BG001-DECISION/{relative}"
        require(candidate_relative in observed, f"unselected BG graph path: {relative}")
        row = observed[candidate_relative]
        require(
            int(row["size_bytes"]) == size and row["sha256"] == digest,
            f"BG graph selection identity drift: {relative}",
        )
        require(
            set(role.split(";")).issubset(set(row["role"].split(";"))),
            f"BG graph selection role drift: {relative}",
        )
        require(
            row["source_scope"] == "project_relative"
            and row["source_path"]
            == f"results/remediation/bg001/{BG_RUN_ID}/{relative}",
            f"BG graph source provenance drift: {relative}",
        )

    reviewed_manifest = CANDIDATE_ROOT / "code/reviewed_code.sha256"
    reviewed_hash = sha256(reviewed_manifest)
    reviewed_rows = parse_sha256_manifest(
        reviewed_manifest,
        separator="  ",
        expected_rows=len(REVIEWED_CODE_FILES),
    )
    expected_reviewed_paths = {
        f"scripts/manuscript/resource_f_five/{name}"
        for name in REVIEWED_CODE_FILES
    }
    require(
        {raw for _, raw in reviewed_rows} == expected_reviewed_paths,
        "reviewed-code manifest file-set drift",
    )
    for digest, raw in reviewed_rows:
        name = PurePosixPath(raw).name
        candidate_relative = f"code/{name}"
        require(candidate_relative in observed, f"reviewed code not selected: {name}")
        require(
            observed[candidate_relative]["sha256"] == digest
            and sha256(CANDIDATE_ROOT / candidate_relative) == digest,
            f"reviewed code identity drift: {name}",
        )
    require(
        observed["code/reviewed_code.sha256"]["sha256"] == reviewed_hash,
        "reviewed manifest selection drift",
    )
    expected_code_files = {
        *(f"code/{name}" for name in REVIEWED_CODE_FILES),
        "code/reviewed_code.sha256",
        "code/frozen_bg001_runtime/analysis_runtime_contract.py",
        "code/frozen_bg001_runtime/safe_io.py",
    }
    actual_code_files = {
        path.relative_to(CANDIDATE_ROOT).as_posix()
        for path in (CANDIDATE_ROOT / "code").rglob("*")
        if path.is_file()
    }
    require(actual_code_files == expected_code_files, "candidate code file-set drift")

    gtf_relative = (
        "inputs/BG001-DECISION/external_reference/"
        "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
    )
    require(gtf_relative in observed, "external GTF snapshot missing")
    gtf_row = observed[gtf_relative]
    require(
        gtf_row["source_scope"] == "external_absolute"
        and gtf_row["source_path"] == GENCODE_GTF.as_posix()
        and int(gtf_row["size_bytes"]) == GENCODE_GTF_BYTES
        and gtf_row["sha256"] == GENCODE_GTF_SHA256,
        "external GTF provenance drift",
    )

    source_manifest_hash = sha256(source_manifest)
    adoption = load_json(CANDIDATE_ROOT / "manifests/adoption_contract.json")
    expected_adoption = {
        "schema": "masld-resource-f-five-candidate-v1",
        "candidate_id": CANDIDATE_ID,
        "bg001_run_id": BG_RUN_ID,
        "bg001_graph_file_count": EXPECTED_GRAPH_FILES,
        "bg001_graph_size_bytes": EXPECTED_GRAPH_BYTES,
        "approved_pooled_layer": "human_bulk_disease_control_fivecohort_candidate",
        "resource_source_layer": "human_bulk_corrected_ninecohort_candidate",
        "canonical_promotion_authorized": False,
        "figure_main_write_authorized": False,
        "portal_deployment_authorized": False,
        "corrected_coloc_status": "pending",
        "pooled_reproduction_required": True,
        "source_manifest_sha256": source_manifest_hash,
        "reviewed_code_manifest_sha256": reviewed_hash,
        "gse193066_metadata_sha256": GSE193066_METADATA_SHA256,
        "external_gtf_sha256": GENCODE_GTF_SHA256,
        "external_gtf_size_bytes": GENCODE_GTF_BYTES,
    }
    for key, value in expected_adoption.items():
        require(adoption.get(key) == value, f"adoption contract drift: {key}")

    object_census_path = CANDIDATE_ROOT / "manifests/bulk_object_census.tsv"
    _, object_rows = load_tsv(object_census_path, ["metric", "value"])
    require(len(object_rows) == len(EXPECTED_OBJECT_CENSUS), "object-census row-count drift")
    require(
        {row["metric"]: int(row["value"]) for row in object_rows}
        == EXPECTED_OBJECT_CENSUS,
        "bulk object census drift",
    )
    require(
        adoption.get("bulk_object_census_sha256") == sha256(object_census_path),
        "adoption object-census binding drift",
    )

    census_path = CANDIDATE_ROOT / "manifests/sample_census.tsv"
    _, census_rows = load_tsv(
        census_path,
        ["layer", "n_samples", "n_genes", "role"],
    )
    require(census_rows == EXPECTED_CENSUS, "sample census drift")

    completion = load_json(CANDIDATE_ROOT / "BASE_SNAPSHOT_COMPLETE.json")
    expected_completion = {
        "status": "BASE_SNAPSHOT_COMPLETE",
        "candidate_id": CANDIDATE_ID,
        "source_file_count": len(rows),
        "source_manifest_sha256": source_manifest_hash,
        "sample_census_sha256": sha256(census_path),
        "bulk_object_census_sha256": sha256(object_census_path),
        "reviewed_code_manifest_sha256": reviewed_hash,
        "adoption_contract_sha256": sha256(
            CANDIDATE_ROOT / "manifests/adoption_contract.json"
        ),
        "sentinel_sha256": sha256(sentinel),
        "canonical_promotion_authorized": False,
    }
    for key, value in expected_completion.items():
        require(completion.get(key) == value, f"base completion drift: {key}")

    base_manifest = CANDIDATE_ROOT / "manifests/base_artifact_manifest.tsv"
    _, base_rows = load_tsv(
        base_manifest,
        ["candidate_path", "size_bytes", "sha256"],
    )
    base_observed: set[str] = set()
    for row in base_rows:
        relative = _safe_relative(row["candidate_path"], "base path")
        relative_text = relative.as_posix()
        require(relative_text not in base_observed, f"duplicate base path: {relative_text}")
        base_observed.add(relative_text)
        path = CANDIDATE_ROOT.joinpath(*relative.parts)
        require_regular_file(path)
        require(
            path.stat().st_size == int(row["size_bytes"]),
            f"base byte drift: {relative_text}",
        )
        require(sha256(path) == row["sha256"], f"base hash drift: {relative_text}")
    excluded_roots = {"workstreams", "logs"}
    actual_base = set()
    for path in CANDIDATE_ROOT.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(CANDIDATE_ROOT)
        relative_text = relative.as_posix()
        if (
            relative.parts[0] in excluded_roots
            or relative_text.startswith("manifests/jobs/")
            or relative_text
            in {
                "BASE_SNAPSHOT_COMPLETE.json",
                "manifests/base_artifact_manifest.tsv",
            }
        ):
            continue
        actual_base.add(relative_text)
    require(base_observed == actual_base, "base-artifact manifest file-set mismatch")
    require(
        completion.get("base_file_count") == len(base_rows),
        "base file-count drift",
    )
    require(
        completion.get("base_artifact_manifest_sha256") == sha256(base_manifest),
        "base seal drift",
    )

    copied_canonical = (
        CANDIDATE_ROOT / "inputs/CANONICAL-REFERENCE/canonical_deg_results.csv"
    )
    copied_gse193066 = (
        CANDIDATE_ROOT
        / "inputs/BULK-NINE-COHORT/source_metadata/GSE193066_metadata.tsv"
    )
    require_regular_file(copied_canonical)
    require(
        sha256(copied_canonical) == CANONICAL_REFERENCE_SHA256,
        "copied canonical reference drift",
    )
    require_regular_file(CANONICAL_REFERENCE)
    require(
        sha256(CANONICAL_REFERENCE) == CANONICAL_REFERENCE_SHA256,
        "live canonical reference drift",
    )
    require_regular_file(copied_gse193066)
    require(
        sha256(copied_gse193066) == GSE193066_METADATA_SHA256,
        "copied GSE193066 metadata drift",
    )

    return {
        "status": "CANDIDATE_SNAPSHOT_VERIFIED",
        "candidate_id": CANDIDATE_ID,
        "source_file_count": len(rows),
        "source_manifest_sha256": source_manifest_hash,
        "bg001_graph_file_count": len(graph),
        "bg001_graph_size_bytes": EXPECTED_GRAPH_BYTES,
        "canonical_promotion_authorized": False,
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), indent=2, sort_keys=True))
    except ContractError as error:
        raise SystemExit(f"SNAPSHOT VERIFICATION ERROR: {error}") from error
