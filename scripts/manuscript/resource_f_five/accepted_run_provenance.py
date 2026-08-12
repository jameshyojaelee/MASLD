#!/usr/bin/env python3
"""Validate and enumerate the complete sealed BG-001 provenance graph."""

from __future__ import annotations

import csv
import json
from pathlib import Path, PurePosixPath

from resource_contract import (
    BASELINE_COMPATIBILITY_SHA256,
    BASELINE_FILES_MANIFEST_SHA256,
    BASELINE_FROZEN_SHA256,
    BASELINE_MANIFESTS_MANIFEST_SHA256,
    BG_ROOT,
    EXPECTED_HASHES,
    GENCODE_GTF,
    GENCODE_GTF_SHA256,
    RECOUNT_ARTIFACTS_MANIFEST_SHA256,
    RECOUNT_COMPLETE_SHA256,
    STRUCTURAL_VALIDATION_SHA256,
    require,
    require_regular_file,
    sha256,
)
from snapshot_io import parse_sha256_manifest


Selection = tuple[Path, str, int, str]
GraphEntry = tuple[str, int, str]

EXPECTED_GRAPH_FILES = 3330
EXPECTED_GRAPH_BYTES = 3_097_918_532
EXPECTED_SOURCE_FILES = 2450
KNOWN_UNSEALED_SOURCE_EXTRA = (
    "source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/"
    "bg001_remediation/__pycache__/safe_io.cpython-310.pyc"
)
CONTRACT_FILES = {
    "analysis_environment.conda.json",
    "analysis_environment.explicit.txt",
    "analysis_native_dependencies.tsv",
    "analysis_r_dependency_files.tsv",
    "analysis_runtime_contract.json",
    "dependency_source_inventory.tsv",
    "dirty_worktree.patch",
    "rendered_figure_label_check.txt",
    "run_contract.json",
    "source_manifest.tsv",
    "source_regression_check.txt",
    "untracked_dependency_sources.tsv",
}
RECOUNT_TRANSACTION_MANIFESTS = {
    *{f"bam_hash_shards/task_{index}.tsv" for index in range(8)},
    "bam_manifest.draft.tsv",
    "bam_manifest.freeze.json",
    "bam_manifest.tsv",
    "bam_manifest_finalize.postflight.json",
    "bam_manifest_finalize.preflight.json",
    "count_manifest.tsv",
    "count_overrides.tsv",
}
MANIFEST_FILES = {
    "bam_hash_shards/task_0.tsv",
    "bam_hash_shards/task_1.tsv",
    "bam_hash_shards/task_2.tsv",
    "bam_hash_shards/task_3.tsv",
    "bam_hash_shards/task_4.tsv",
    "bam_hash_shards/task_5.tsv",
    "bam_hash_shards/task_6.tsv",
    "bam_hash_shards/task_7.tsv",
    "bam_manifest.draft.tsv",
    "bam_manifest.freeze.json",
    "bam_manifest.tsv",
    "bam_manifest_finalize.postflight.json",
    "bam_manifest_finalize.preflight.json",
    "baseline_frozen_files.sha256",
    "baseline_manifests.sha256",
    "count_manifest.tsv",
    "count_overrides.tsv",
    "protected_source_manifest.tsv",
    "read_count_overrides.tsv",
    "read_count_sources.tsv",
    "recount_artifacts.sha256",
}
BASELINE_BOUND_INPUTS = {
    "manifests/read_count_sources.tsv",
    "manifests/read_count_overrides.tsv",
    "manifests/protected_source_manifest.tsv",
    "manifests/bam_manifest_finalize.preflight.json",
    "manifests/bam_manifest_finalize.postflight.json",
    "contract/run_contract.json",
    "contract/source_manifest.tsv",
    "contract/analysis_runtime_contract.json",
    "contract/analysis_environment.explicit.txt",
    "contract/analysis_environment.conda.json",
    "contract/analysis_r_dependency_files.tsv",
    "contract/analysis_native_dependencies.tsv",
}
ANALYSIS_BOUND_INPUTS = {
    "BASELINE_FROZEN.json",
    "RECOUNT_COMPLETE",
    *{
        f"arms/{arm}/{name}"
        for arm in ("F_five", "F_legacy", "F_locked", "R0")
        for name in ("ARM_COMPLETE.json", "provenance/artifact_manifest.tsv")
    },
    *{
        f"contract/{name}"
        for name in (
            "analysis_environment.conda.json",
            "analysis_environment.explicit.txt",
            "analysis_native_dependencies.tsv",
            "analysis_r_dependency_files.tsv",
            "analysis_runtime_contract.json",
            "run_contract.json",
            "source_manifest.tsv",
        )
    },
}


def _load_json(path: Path) -> dict[str, object]:
    require_regular_file(path)
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def _strict_tsv(path: Path, fields: list[str]) -> list[dict[str, str]]:
    require_regular_file(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(list(reader.fieldnames or []) == fields, f"TSV schema drift: {path}")
        rows = list(reader)
    require(
        all(None not in row and all(value is not None for value in row.values()) for row in rows),
        f"malformed TSV row width: {path}",
    )
    return rows


def _safe_relative(raw: str) -> str:
    path = PurePosixPath(raw)
    require(not path.is_absolute() and ".." not in path.parts, f"unsafe relative path: {raw}")
    normalized = path.as_posix()
    require(normalized not in {"", "."}, f"empty relative path: {raw}")
    return normalized


def _historical_absolute_to_relative(raw: str) -> str:
    path = Path(raw)
    require(path.is_absolute() and ".." not in path.parts, f"unsafe historical path: {raw}")
    try:
        return path.relative_to(BG_ROOT).as_posix()
    except ValueError as error:
        raise RuntimeError(f"historical path escapes accepted BG root: {raw}") from error


def _add(
    graph: dict[str, GraphEntry],
    root: Path,
    relative: str,
    role: str,
    *,
    expected_hash: str | None = None,
    expected_size: int | None = None,
) -> None:
    relative = _safe_relative(relative)
    path = root / relative
    require_regular_file(path)
    size = path.stat().st_size
    digest = sha256(path)
    if expected_size is not None:
        require(size == expected_size, f"sealed size drift: {relative}")
    if expected_hash is not None:
        require(digest == expected_hash, f"sealed hash drift: {relative}")
    entry = (role, size, digest)
    if relative in graph:
        old_role, old_size, old_digest = graph[relative]
        require(
            old_size == size and old_digest == digest,
            f"conflicting sealed identity: {relative}",
        )
        roles = sorted(set(old_role.split(";")) | set(role.split(";")))
        graph[relative] = (";".join(roles), size, digest)
        return
    graph[relative] = entry


def _manifest_payloads(
    graph: dict[str, GraphEntry],
    root: Path,
    manifest_relative: str,
    *,
    manifest_hash: str,
    role: str,
) -> set[str]:
    _add(
        graph,
        root,
        manifest_relative,
        f"{role}_manifest",
        expected_hash=manifest_hash,
    )
    rows = _strict_tsv(
        root / manifest_relative,
        ["relative_path", "size_bytes", "sha256"],
    )
    require(rows, f"empty artifact manifest: {manifest_relative}")
    payloads: set[str] = set()
    manifest_parent = PurePosixPath(manifest_relative).parent
    artifact_root = (
        manifest_parent.parent
        if manifest_parent.name == "provenance"
        else manifest_parent
    )
    for row in rows:
        child = _safe_relative(row["relative_path"])
        relative = (artifact_root / child).as_posix()
        require(relative not in payloads, f"duplicate artifact path: {relative}")
        payloads.add(relative)
        _add(
            graph,
            root,
            relative,
            role,
            expected_hash=row["sha256"],
            expected_size=int(row["size_bytes"]),
        )
    return payloads


def _baseline(graph: dict[str, GraphEntry], root: Path) -> None:
    identities = {
        "BASELINE_FROZEN.json": BASELINE_FROZEN_SHA256,
        "BASELINE_FROZEN": BASELINE_COMPATIBILITY_SHA256,
        "manifests/baseline_frozen_files.sha256": BASELINE_FILES_MANIFEST_SHA256,
        "manifests/baseline_manifests.sha256": BASELINE_MANIFESTS_MANIFEST_SHA256,
    }
    for relative, digest in identities.items():
        _add(graph, root, relative, "bg001_baseline_transaction", expected_hash=digest)

    baseline = _load_json(root / "BASELINE_FROZEN.json")
    require(
        set(baseline)
        == {
            "baseline_frozen_files_manifest",
            "baseline_frozen_files_manifest_sha256",
            "baseline_manifests_manifest",
            "baseline_manifests_manifest_sha256",
            "bound_input_count",
            "compatibility_marker",
            "compatibility_marker_contract",
            "frozen_file_count",
            "run_contract_sha256",
            "run_id",
            "schema",
            "source_manifest_sha256",
        },
        "baseline seal schema drift",
    )
    require(baseline["schema"] == "bg001-baseline-seal-v1", "baseline schema drift")
    require(baseline["run_id"] == BG_ROOT.name, "baseline run ID drift")
    require(baseline["frozen_file_count"] == 624, "baseline frozen count drift")
    require(baseline["bound_input_count"] == 12, "baseline bound-input count drift")
    require(
        baseline["baseline_frozen_files_manifest_sha256"]
        == BASELINE_FILES_MANIFEST_SHA256
        and baseline["baseline_manifests_manifest_sha256"]
        == BASELINE_MANIFESTS_MANIFEST_SHA256,
        "baseline manifest binding drift",
    )
    require(
        (root / "BASELINE_FROZEN").read_text(encoding="ascii")
        == BASELINE_FROZEN_SHA256 + "\n",
        "baseline compatibility marker drift",
    )

    frozen_rows = parse_sha256_manifest(
        root / "manifests/baseline_frozen_files.sha256",
        separator="  ",
        expected_rows=624,
    )
    frozen_relatives: set[str] = set()
    for digest, raw in frozen_rows:
        relative = _historical_absolute_to_relative(raw)
        require(relative.startswith("frozen_sets/"), f"non-frozen baseline payload: {relative}")
        require(relative not in frozen_relatives, f"duplicate baseline payload: {relative}")
        frozen_relatives.add(relative)
        _add(
            graph,
            root,
            relative,
            "bg001_baseline_transaction",
            expected_hash=digest,
        )
    actual_frozen = {
        path.relative_to(root).as_posix()
        for path in (root / "frozen_sets").rglob("*")
        if path.is_file()
    }
    require(actual_frozen == frozen_relatives, "baseline frozen-set inventory drift")

    input_rows = parse_sha256_manifest(
        root / "manifests/baseline_manifests.sha256",
        separator="  ",
        expected_rows=12,
    )
    input_relatives: set[str] = set()
    for digest, raw in input_rows:
        relative = _historical_absolute_to_relative(raw)
        input_relatives.add(relative)
        _add(
            graph,
            root,
            relative,
            "bg001_baseline_transaction",
            expected_hash=digest,
        )
    require(input_relatives == BASELINE_BOUND_INPUTS, "baseline bound-input inventory drift")
    require(
        baseline["run_contract_sha256"] == sha256(root / "contract/run_contract.json")
        and baseline["source_manifest_sha256"] == sha256(root / "contract/source_manifest.tsv"),
        "baseline contract binding drift",
    )
    baseline_relatives = {
        *identities,
        *frozen_relatives,
        *input_relatives,
    }
    require(len(baseline_relatives) == 640, "baseline closure cardinality drift")
    require(
        sum((root / relative).stat().st_size for relative in baseline_relatives)
        == 1_261_932_287,
        "baseline closure byte-count drift",
    )


def _source_and_run_contract(graph: dict[str, GraphEntry], root: Path) -> None:
    sentinel_hash = "1005091a7bf8db81cde82977ea7908b9cf846b4a278ae817203c2347cb252174"
    _add(
        graph,
        root,
        ".bg001_candidate_root",
        "bg001_run_identity",
        expected_hash=sentinel_hash,
    )
    require(
        (root / ".bg001_candidate_root").read_text(encoding="utf-8").strip()
        == BG_ROOT.name,
        "accepted run sentinel drift",
    )
    contract = _load_json(root / "contract/run_contract.json")
    require(contract.get("schema_version") == "1.2", "run-contract schema drift")
    require(contract.get("run_id") == BG_ROOT.name, "run-contract ID drift")
    require(contract.get("candidate_root") == str(BG_ROOT), "historical candidate-root drift")
    require(contract.get("project_root") == str(BG_ROOT.parents[3]), "historical project-root drift")
    featurecounts = contract.get("featurecounts")
    require(
        featurecounts
        == {
            "arguments": ["-p", "--countReadPairs", "-B", "-s", "2"],
            "gtf": str(GENCODE_GTF),
            "gtf_sha256": GENCODE_GTF_SHA256,
            "version": "2.1.1",
        },
        "featureCounts/GTF run contract drift",
    )

    support_hashes = {
        "contract/dependency_source_inventory.tsv": contract["dependency_source_inventory_sha256"],
        "contract/dirty_worktree.patch": contract["dirty_patch_sha256"],
        "contract/rendered_figure_label_check.txt": contract["rendered_figure_label_check_sha256"],
        "contract/source_regression_check.txt": contract["source_regression_check_sha256"],
        "contract/untracked_dependency_sources.tsv": contract["untracked_dependency_sources_sha256"],
    }
    for relative, digest in support_hashes.items():
        _add(graph, root, relative, "bg001_run_contract_support", expected_hash=str(digest))

    source_manifest = root / "contract/source_manifest.tsv"
    require(
        sha256(source_manifest) == contract["source_manifest_sha256"],
        "source-manifest/run-contract binding drift",
    )
    rows = _strict_tsv(source_manifest, ["relative_path", "size_bytes", "sha256"])
    require(len(rows) == EXPECTED_SOURCE_FILES, "source-snapshot row-count drift")
    expected: set[str] = set()
    for row in rows:
        child = _safe_relative(row["relative_path"])
        relative = f"source_snapshot/{child}"
        require(relative not in expected, f"duplicate source-snapshot path: {child}")
        expected.add(relative)
        _add(
            graph,
            root,
            relative,
            "bg001_frozen_source_snapshot",
            expected_hash=row["sha256"],
            expected_size=int(row["size_bytes"]),
        )
    actual = {
        path.relative_to(root).as_posix()
        for path in (root / "source_snapshot").rglob("*")
        if path.is_file()
    }
    if root.resolve(strict=True) == BG_ROOT.resolve(strict=True):
        require(
            actual == expected | {KNOWN_UNSEALED_SOURCE_EXTRA},
            "accepted source-snapshot inventory drift",
        )
    else:
        require(actual == expected, "mirrored source-snapshot inventory drift")


def _bam_and_recount(graph: dict[str, GraphEntry], root: Path) -> None:
    run_contract = _load_json(root / "contract/run_contract.json")
    draft_hash = str(run_contract["bam_manifest_draft_sha256"])
    _add(
        graph,
        root,
        "manifests/bam_manifest.draft.tsv",
        "bg001_bam_transaction",
        expected_hash=draft_hash,
    )
    postflight = _load_json(root / "manifests/bam_manifest_finalize.postflight.json")
    preflight = root / "manifests/bam_manifest_finalize.preflight.json"
    require(postflight.get("status") == "PASS", "BAM postflight status drift")
    require(postflight.get("run_id") == BG_ROOT.name, "BAM postflight run ID drift")
    require(postflight.get("included_bams") == 820, "BAM count drift")
    require(postflight.get("preflight_sha256") == sha256(preflight), "BAM preflight binding drift")
    inputs = postflight.get("input_sha256")
    require(isinstance(inputs, dict), "BAM postflight input schema drift")
    require(
        set(inputs)
        == {
            *{f"manifests/bam_hash_shards/task_{index}.tsv" for index in range(8)},
            "manifests/bam_manifest.draft.tsv",
        },
        "BAM postflight input inventory drift",
    )
    for relative, digest in inputs.items():
        _add(
            graph,
            root,
            str(relative),
            "bg001_bam_transaction",
            expected_hash=str(digest),
        )
    transaction = postflight.get("transaction")
    require(isinstance(transaction, dict), "BAM transaction schema drift")
    transaction_hashes = {
        "manifests/bam_manifest.freeze.json": transaction["freeze_json_sha256"],
        "manifests/bam_manifest.tsv": transaction["manifest_sha256"],
        "BAM_MANIFEST_FROZEN": transaction["marker_sha256"],
    }
    for relative, digest in transaction_hashes.items():
        _add(
            graph,
            root,
            relative,
            "bg001_bam_transaction",
            expected_hash=str(digest),
        )
    require(
        (root / "BAM_MANIFEST_FROZEN").read_text(encoding="ascii")
        == str(transaction["manifest_sha256"]) + "\n",
        "BAM compatibility marker drift",
    )
    freeze = _load_json(root / "manifests/bam_manifest.freeze.json")
    require(freeze.get("run_id") == BG_ROOT.name, "BAM freeze run ID drift")
    require(freeze.get("included_bams") == 820, "BAM freeze count drift")
    require(freeze.get("draft_sha256") == draft_hash, "BAM freeze draft binding drift")
    require(
        freeze.get("manifest_sha256") == transaction["manifest_sha256"],
        "BAM freeze manifest binding drift",
    )
    require(
        freeze.get("shard_sha256")
        == {
            Path(relative).name: digest
            for relative, digest in inputs.items()
            if relative.startswith("manifests/bam_hash_shards/")
        },
        "BAM freeze shard binding drift",
    )

    completion_path = root / "RECOUNT_COMPLETE"
    manifest_path = root / "manifests/recount_artifacts.sha256"
    structural_path = root / "comparisons/structural_validation.json"
    for relative, digest in {
        "RECOUNT_COMPLETE": RECOUNT_COMPLETE_SHA256,
        "manifests/recount_artifacts.sha256": RECOUNT_ARTIFACTS_MANIFEST_SHA256,
        "comparisons/structural_validation.json": STRUCTURAL_VALIDATION_SHA256,
    }.items():
        _add(graph, root, relative, "bg001_recount_transaction", expected_hash=digest)
    completion = _load_json(completion_path)
    require(
        completion
        == {
            "recount_artifacts_sha256": RECOUNT_ARTIFACTS_MANIFEST_SHA256,
            "structural_validation_sha256": STRUCTURAL_VALIDATION_SHA256,
        },
        "recount completion binding drift",
    )
    rows = parse_sha256_manifest(manifest_path, separator="\t", expected_rows=87)
    payloads: set[str] = set()
    for digest, raw in rows:
        relative = _safe_relative(raw)
        require(relative not in payloads, f"duplicate recount artifact: {relative}")
        payloads.add(relative)
        _add(
            graph,
            root,
            relative,
            "bg001_recount_transaction",
            expected_hash=digest,
        )
    expected_payloads = (
        {"BAM_MANIFEST_FROZEN"}
        | {
            path.relative_to(root).as_posix()
            for path in (root / "counts").rglob("*")
            if path.is_file()
        }
        | {f"manifests/{relative}" for relative in RECOUNT_TRANSACTION_MANIFESTS}
    )
    require(payloads == expected_payloads, "recount artifact inventory drift")

    structural = _load_json(structural_path)
    require(structural.get("status") == "PASS", "structural validation status drift")
    require(
        structural.get("recount_artifacts_sha256") == RECOUNT_ARTIFACTS_MANIFEST_SHA256,
        "structural/recount binding drift",
    )
    for relative, key in (
        ("manifests/count_manifest.tsv", "count_manifest_sha256"),
        ("manifests/count_overrides.tsv", "count_overrides_sha256"),
    ):
        _add(
            graph,
            root,
            relative,
            "bg001_recount_support",
            expected_hash=str(structural[key]),
        )


def _analysis(graph: dict[str, GraphEntry], root: Path) -> None:
    analysis_hash = EXPECTED_HASHES["ANALYSIS_COMPLETE.json"]
    _add(
        graph,
        root,
        "ANALYSIS_COMPLETE.json",
        "bg001_analysis_transaction",
        expected_hash=analysis_hash,
    )
    analysis = _load_json(root / "ANALYSIS_COMPLETE.json")
    require(analysis.get("schema") == "bg001-analysis-artifact-seal-v1", "analysis schema drift")
    require(analysis.get("run_id") == BG_ROOT.name, "analysis run ID drift")
    require(analysis.get("manifest") == "comparisons/artifact_manifest.tsv", "analysis manifest path drift")
    require(analysis.get("manifest_sha256") == EXPECTED_HASHES["comparisons/artifact_manifest.tsv"],
            "analysis manifest hash drift")
    bound = analysis.get("bound_inputs")
    require(isinstance(bound, dict) and set(bound) == ANALYSIS_BOUND_INPUTS,
            "analysis bound-input inventory drift")
    for relative, digest in bound.items():
        _add(
            graph,
            root,
            str(relative),
            "bg001_analysis_bound_input",
            expected_hash=str(digest),
        )

    for arm in ("F_five", "F_legacy", "F_locked", "R0"):
        complete_relative = f"arms/{arm}/ARM_COMPLETE.json"
        manifest_relative = f"arms/{arm}/provenance/artifact_manifest.tsv"
        complete_hash = str(bound[complete_relative])
        manifest_hash = str(bound[manifest_relative])
        _add(
            graph,
            root,
            complete_relative,
            f"bg001_{arm.lower()}_seal",
            expected_hash=complete_hash,
        )
        complete = _load_json(root / complete_relative)
        require(complete.get("schema") == "bg001-analysis-artifact-seal-v1", f"{arm} schema drift")
        require(complete.get("run_id") == BG_ROOT.name, f"{arm} run ID drift")
        require(complete.get("label") == arm, f"{arm} label drift")
        require(complete.get("manifest") == manifest_relative, f"{arm} manifest path drift")
        require(complete.get("manifest_sha256") == manifest_hash, f"{arm} manifest hash drift")
        payloads = _manifest_payloads(
            graph,
            root,
            manifest_relative,
            manifest_hash=manifest_hash,
            role=f"bg001_{arm.lower()}_artifact",
        )
        require(complete.get("artifact_count") == len(payloads), f"{arm} artifact count drift")
        actual = {
            path.relative_to(root).as_posix()
            for path in (root / f"arms/{arm}").rglob("*")
            if path.is_file()
        }
        require(
            actual == payloads | {complete_relative, manifest_relative},
            f"{arm} sealed file-set drift",
        )

    comparison_manifest = "comparisons/artifact_manifest.tsv"
    comparison_payloads = _manifest_payloads(
        graph,
        root,
        comparison_manifest,
        manifest_hash=str(analysis["manifest_sha256"]),
        role="bg001_comparison_artifact",
    )
    require(analysis.get("artifact_count") == len(comparison_payloads), "comparison artifact count drift")
    actual_comparisons = {
        path.relative_to(root).as_posix()
        for path in (root / "comparisons").rglob("*")
        if path.is_file()
    }
    require(
        actual_comparisons == comparison_payloads | {comparison_manifest},
        "comparison sealed file-set drift",
    )
    _add(
        graph,
        root,
        "promotion/CANDIDATE_ACCEPTANCE.md",
        "bg001_acceptance_decision",
        expected_hash=EXPECTED_HASHES["promotion/CANDIDATE_ACCEPTANCE.md"],
    )


def verify_accepted_graph(
    root: Path,
) -> dict[str, GraphEntry]:
    """Verify the full 3,330-file sealed graph at a source or mirrored root."""

    root = root.resolve(strict=True)
    require(root.is_dir() and not root.is_symlink(), f"invalid accepted graph root: {root}")
    for path in root.rglob("*"):
        require(not path.is_symlink(), f"symlink in accepted graph: {path}")

    graph: dict[str, GraphEntry] = {}
    _baseline(graph, root)
    _source_and_run_contract(graph, root)
    _bam_and_recount(graph, root)
    _analysis(graph, root)

    actual_contract = {
        path.name
        for path in (root / "contract").iterdir()
        if path.is_file()
    }
    if root == BG_ROOT.resolve(strict=True):
        require(
            actual_contract == CONTRACT_FILES | {"git_status_porcelain_v2.txt"},
            "accepted contract directory drift",
        )
    else:
        require(actual_contract == CONTRACT_FILES, "mirrored contract directory drift")

    actual_manifests = {
        path.relative_to(root / "manifests").as_posix()
        for path in (root / "manifests").rglob("*")
        if path.is_file()
    }
    require(actual_manifests == MANIFEST_FILES, "manifest directory closure drift")
    for relative in CONTRACT_FILES:
        require(f"contract/{relative}" in graph, f"unbound contract file: {relative}")
    for relative in MANIFEST_FILES:
        require(f"manifests/{relative}" in graph, f"unbound manifest file: {relative}")

    require(len(graph) == EXPECTED_GRAPH_FILES, "accepted graph file-count drift")
    require(
        sum(size for _, size, _ in graph.values()) == EXPECTED_GRAPH_BYTES,
        "accepted graph byte-count drift",
    )
    return graph


def collect_accepted_run_selections() -> dict[str, Selection]:
    """Return the verified, path-preserving sealed accepted-run graph."""

    graph = verify_accepted_graph(BG_ROOT)
    selections: dict[str, Selection] = {}
    for relative, (role, size, digest) in graph.items():
        destination = f"inputs/BG001-DECISION/{relative}"
        selections[destination] = (BG_ROOT / relative, role, size, digest)
    return selections
