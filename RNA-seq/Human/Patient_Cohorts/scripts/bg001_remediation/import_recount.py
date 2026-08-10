#!/usr/bin/env python3
"""Fail-closed, hash-verified import of one sealed BG-001 recount into a fresh run.

The count-producing bytes remain attributed to the source execution snapshot.
Only path/run-ID metadata that cannot be valid in another candidate is regenerated;
the exact originals are archived in the destination import lineage.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


sys.dont_write_bytecode = True

RUN_RE = re.compile(r"^bg001-fragment-v211-gencode49-[0-9]{8}T[0-9]{6}Z$")
SHA_RE = re.compile(r"^[0-9a-f]{64}$")
REMEDIATION_REL = Path("RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation")
EXPECTED_GTF_SHA256 = "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
EXPECTED_FEATURECOUNTS = {
    "version": "2.1.1",
    "arguments": ["-p", "--countReadPairs", "-B", "-s", "2"],
    "gtf": "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz",
    "gtf_sha256": EXPECTED_GTF_SHA256,
}
EXPECTED_SCHEDULER_ATTESTATION_SCHEMA = "bg001-source-slurm-attestation-v1"

# These files determined the BAM identity, featureCounts invocation, validation,
# chunk membership, keyed merge, and publication of the count bytes.  Analysis,
# protected-set, and import-aware post-count metadata code may intentionally differ.
PRODUCER_FILENAMES = (
    "bam_manifest_hash_array.sbatch",
    "bam_manifest_finalize.sbatch",
    "build_bam_manifest.py",
    "finalize_bam_manifest.py",
    "hash_bam_manifest_task.sh",
    "recount_array.sbatch",
    "run_featurecounts_contract.sh",
    "validate_featurecounts.py",
    "merge_validate.sbatch",
    "merge_gse213621.R",
    "publish_gse213621.py",
    "verify_run_contract.py",
)
PRODUCER_PATHS = tuple(str(REMEDIATION_REL / name) for name in PRODUCER_FILENAMES)
DESTINATION_GUARD_MITIGATION_PATHS = tuple(str(REMEDIATION_REL / name) for name in (
    "bam_manifest_guarded_hash_array.sbatch",
    "bam_manifest_guarded_finalize.sbatch",
    "bam_manifest_finalization_guard.py",
    "safe_io.py",
))
REGRESSION_CHECKER_REL = str(REMEDIATION_REL / "validate_featurecounts_sources.py")
REGRESSION_CHECKER_DIGESTS = {
    "legacy_055119Z": {
        "sha256": "0bcbec61fda2cb0bfab01519e57506446978b5f07e2d6635b840cd36c9d47ea9",
        "size": 11445,
    },
    "slurm_directive_aware": {
        "sha256": "1615db6646ea15b022b82ab7708a8d777dd21e3bfd664c8b536beb4ede2c6454",
        "size": 12415,
    },
    "guarded_routes_v1": {
        "sha256": "543d5623627a5787cdfb90b855c09f682e280f6f1ef1442ceec6eb1a5d0b370e",
        "size": 18289,
    },
}
LEGACY_REGRESSION_CHECKER_EXPECTED_OUTPUT = (
    "PASS: 19 fixed human, 1 dynamic shared, and 2 mouse paired producers "
    "satisfy the fragment-counting contract\n"
)
GUARDED_REGRESSION_CHECKER_EXPECTED_OUTPUT = (
    "PASS: 14 fixed human, 5 guarded custom, 1 guarded dynamic shared, and 2 mouse "
    "paired producers satisfy the fragment-counting contract\n"
)
GUARDED_CHECKER_ON_LEGACY_EXPECTED = {
    "returncode": 1,
    "stdout_sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    "stderr_sha256": "d9e4e1f226f491d96d9a016e57cf1d871310efb11512280e9829ce8b474aa3ba",
    "stdout_bytes": 0,
    "stderr_bytes": 1211,
}
POSTCOUNT_METADATA_FILENAMES = (
    "finalize_recount.py",
    "make_count_overrides.py",
    "verify_recount_artifacts.py",
)
PROVENANCE_PAYLOADS = (
    "bam_checksums.tsv",
    "environment.txt",
    "gene_counts.txt",
    "gene_counts.txt.summary",
    "featureCounts.log",
    "validation.json",
)
MERGE_PAYLOADS = (
    "gene_counts.txt",
    "gene_counts.txt.summary",
    "featureCounts.log",
    "environment.txt",
    "validation.json",
)
MERGE_PUBLISH_ORDER = (
    "environment.txt",
    "featureCounts.log",
    "gene_counts.txt.summary",
    "validation.json",
    "gene_counts.txt",
)
PATH_BOUND_COUNT_FILES = {
    *(f"counts/{dataset}/provenance.sha256" for dataset in (
        "GSE130970", "GSE135251", "GSE174478", "GSE240729"
    )),
    *(f"counts/GSE213621/chunks/chunk{index}/provenance.sha256" for index in range(4)),
    "counts/GSE213621/merge_publication.json",
    "counts/GSE213621/MERGE_COMPLETE",
}
PROTECTED_GENE_REL = "frozen_sets/protected_genes.tsv"
NOT_TESTABLE_REL = "frozen_sets/named_genes_not_baseline_testable.tsv"
RENDERED_VALIDATION_REL = "frozen_sets/rendered_figure_label_validation.txt"
SUBSTRATE_CENSUS_REL = "frozen_sets/integration_substrate_census.tsv"
PROTECTED_SOURCE_MANIFEST_REL = "manifests/protected_source_manifest.tsv"
PROTECTED_SOURCE_PREFIX = "frozen_sets/protected_sources/"
RUN_BOUND_BASELINE_MANIFESTS = {
    "manifests/read_count_sources.tsv",
    "manifests/read_count_overrides.tsv",
}
RUNTIME_BASELINE_INPUTS = {
    "contract/analysis_runtime_contract.json",
    "contract/analysis_environment.explicit.txt",
    "contract/analysis_environment.conda.json",
    "contract/analysis_r_dependency_files.tsv",
    "contract/analysis_native_dependencies.tsv",
}
BAM_GUARD_BASELINE_INPUTS = {
    "manifests/bam_manifest_finalize.preflight.json",
    "manifests/bam_manifest_finalize.postflight.json",
}
FRESH_CONTROL_BASELINE_INPUTS = {
    "contract/run_contract.json",
    "contract/source_manifest.tsv",
    *BAM_GUARD_BASELINE_INPUTS,
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fsync_file(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fsync_tree_directories(root: Path) -> None:
    directories = [path for path in root.rglob("*") if path.is_dir() and not path.is_symlink()]
    for path in sorted(directories, key=lambda value: len(value.parts), reverse=True):
        fsync_directory(path)
    fsync_directory(root)


def atomic_write_text(path: Path, text: str, *, replace_existing: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not replace_existing and (path.exists() or path.is_symlink()):
        raise SystemExit(f"Refusing to overwrite existing candidate artifact: {path}")
    temporary = path.with_name(f".{path.name}.tmp.{uuid.uuid4().hex}")
    with temporary.open("x") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    if not replace_existing and (path.exists() or path.is_symlink()):
        temporary.unlink()
        raise SystemExit(f"Candidate artifact appeared during atomic write: {path}")
    os.replace(temporary, path)
    fsync_directory(path.parent)


def atomic_write_json(
    path: Path,
    payload: dict[str, object],
    *,
    replace_existing: bool = False,
) -> None:
    atomic_write_text(
        path,
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        replace_existing=replace_existing,
    )


def contained(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def require_regular(path: Path, label: str, *, nonempty: bool = True) -> None:
    if not path.is_file() or path.is_symlink():
        raise SystemExit(f"{label} is missing, nonregular, or symlinked: {path}")
    if nonempty and path.stat().st_size == 0:
        raise SystemExit(f"{label} is empty: {path}")


def copy_independent_regular(
    source: Path,
    target: Path,
    expected_sha256: str,
    label: str,
) -> None:
    """Copy one file and prove that the destination is not a hardlink or symlink."""
    # COMPLETE is intentionally a zero-byte task marker; the expected digest
    # still authenticates it, so independence must not require nonempty bytes.
    require_regular(source, label, nonempty=False)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        raise SystemExit(f"Refusing existing copy destination: {target}")
    shutil.copy2(source, target, follow_symlinks=False)
    if not target.is_file() or target.is_symlink():
        raise SystemExit(f"Copied artifact is nonregular or symlinked: {target}")
    source_stat = source.stat()
    target_stat = target.stat()
    if (
        (source_stat.st_dev, source_stat.st_ino) == (target_stat.st_dev, target_stat.st_ino)
        or target_stat.st_nlink != 1
    ):
        raise SystemExit(f"Copied artifact is hard-linked instead of independent: {target}")
    fsync_file(target)
    if sha256(target) != expected_sha256:
        raise SystemExit(f"Copied artifact digest differs: {target}")


def load_json(path: Path, label: str) -> dict[str, object]:
    require_regular(path, label)
    try:
        payload = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        raise SystemExit(f"Malformed {label}: {path}") from exc
    if not isinstance(payload, dict):
        raise SystemExit(f"{label} is not a JSON object: {path}")
    return payload


def validate_candidate_root(raw_root: Path, role: str) -> tuple[Path, dict[str, object]]:
    if raw_root.is_symlink():
        raise SystemExit(f"{role} run root may not be a symlink: {raw_root}")
    root = raw_root.resolve(strict=True)
    if not root.is_dir() or not RUN_RE.fullmatch(root.name):
        raise SystemExit(f"Invalid {role} BG-001 run root: {root}")
    sentinel = root / ".bg001_candidate_root"
    require_regular(sentinel, f"{role} candidate sentinel")
    if sentinel.read_text().strip() != root.name:
        raise SystemExit(f"{role} candidate sentinel does not name its run")
    if (root / "REJECTED").exists() or (root / "REJECTED").is_symlink():
        raise SystemExit(f"{role} run is permanently rejected: {root}")
    contract = load_json(root / "contract/run_contract.json", f"{role} run contract")
    project = Path(str(contract.get("project_root", ""))).resolve(strict=True)
    expected_parent = (project / "results/remediation/bg001").resolve(strict=True)
    if root.parent != expected_parent:
        raise SystemExit(f"{role} run escapes the project candidate parent: {root}")
    if (
        contract.get("run_id") != root.name
        or Path(str(contract.get("candidate_root", ""))).resolve(strict=True) != root
        or Path(str(contract.get("reference_root", ""))).resolve(strict=True)
        != (root / "source_snapshot").resolve(strict=True)
    ):
        raise SystemExit(f"{role} run contract paths/ID do not describe the candidate")
    return root, contract


def frozen_script(root: Path, filename: str) -> Path:
    path = root / "source_snapshot" / REMEDIATION_REL / filename
    require_regular(path, f"frozen helper {filename}")
    return path


def require_frozen_importer_invocation(root: Path) -> None:
    expected = frozen_script(root, "import_recount.py").resolve(strict=True)
    observed = Path(__file__).resolve(strict=True)
    if observed != expected:
        raise SystemExit(
            "Import must execute the destination's frozen snapshot: "
            f"{sys.executable} {expected} ..."
        )


def run_checked(command: list[str], *, label: str) -> str:
    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    if result.returncode != 0:
        detail = (result.stderr or result.stdout).strip()
        raise SystemExit(f"{label} failed ({result.returncode}): {detail[:4000]}")
    return result.stdout + result.stderr


def run_recorded(command: list[str]) -> dict[str, object]:
    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    return {
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "stdout_sha256": hashlib.sha256(result.stdout.encode()).hexdigest(),
        "stderr_sha256": hashlib.sha256(result.stderr.encode()).hexdigest(),
        "stdout_bytes": len(result.stdout.encode()),
        "stderr_bytes": len(result.stderr.encode()),
    }


def verify_frozen_contract(root: Path) -> str:
    return run_checked(
        [
            sys.executable,
            str(frozen_script(root, "verify_run_contract.py")),
            "--run-root",
            str(root),
            "--require-bam-frozen",
        ],
        label=f"frozen run-contract verification for {root.name}",
    )


def verify_source_recount(root: Path) -> str:
    return run_checked(
        [
            sys.executable,
            str(frozen_script(root, "verify_recount_artifacts.py")),
            "--run-root",
            str(root),
        ],
        label=f"sealed recount verification for {root.name}",
    )


def verify_structured_analysis_baseline(root: Path) -> str:
    """Run the candidate's own analysis-contract verifier for fresh seals."""
    if not (root / "BASELINE_FROZEN.json").exists():
        return "legacy baseline compatibility accepted by importer only\n"
    return run_checked(
        [
            sys.executable,
            str(frozen_script(root, "verify_analysis_contract.py")),
            "--run-root",
            str(root),
            "--require-baseline-frozen",
        ],
        label=f"structured analysis baseline verification for {root.name}",
    )


def load_source_manifest(root: Path, contract: dict[str, object]) -> dict[str, dict[str, object]]:
    path = root / "contract/source_manifest.tsv"
    require_regular(path, "source manifest")
    if sha256(path) != contract.get("source_manifest_sha256"):
        raise SystemExit(f"Source-manifest digest differs from run contract: {root}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("relative_path", "size_bytes", "sha256"):
            raise SystemExit(f"Unexpected source-manifest schema: {path}")
        rows = list(reader)
    result: dict[str, dict[str, object]] = {}
    for row in rows:
        relative = row["relative_path"]
        rel = Path(relative)
        if rel.is_absolute() or ".." in rel.parts or relative in result:
            raise SystemExit(f"Unsafe or duplicate source-manifest path: {relative}")
        try:
            size = int(row["size_bytes"])
        except ValueError as exc:
            raise SystemExit(f"Invalid source-manifest size: {relative}") from exc
        if size < 0 or not SHA_RE.fullmatch(row["sha256"]):
            raise SystemExit(f"Invalid source-manifest record: {relative}")
        result[relative] = {"size": size, "sha256": row["sha256"]}
    return result


def source_manifest_differences(
    source: dict[str, dict[str, object]],
    destination: dict[str, dict[str, object]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for relative in sorted(set(source) | set(destination)):
        left, right = source.get(relative), destination.get(relative)
        if left == right:
            continue
        status = "changed" if left and right else ("source_only" if left else "destination_only")
        rows.append(
            {
                "relative_path": relative,
                "status": status,
                "execution_source_sha256": str(left["sha256"]) if left else "",
                "destination_analysis_sha256": str(right["sha256"]) if right else "",
                "execution_source_size": str(left["size"]) if left else "",
                "destination_analysis_size": str(right["size"]) if right else "",
            }
        )
    return rows


def validate_destination_guard_mitigations(
    destination_root: Path,
    destination_manifest: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for relative in DESTINATION_GUARD_MITIGATION_PATHS:
        record = destination_manifest.get(relative)
        path = destination_root / "source_snapshot" / relative
        if record is None:
            raise SystemExit(f"Destination source manifest lacks BAM guard mitigation: {relative}")
        require_regular(path, f"destination BAM guard mitigation {relative}")
        size = path.stat().st_size
        digest = sha256(path)
        if size != int(record.get("size", -1)) or digest != record.get("sha256"):
            raise SystemExit(f"Destination BAM guard mitigation differs from source manifest: {relative}")
        rows.append({"relative_path": relative, "size_bytes": size, "sha256": digest})
    return rows


def compare_producer_manifests(
    source_root: Path,
    destination_root: Path,
    source_manifest: dict[str, dict[str, object]],
    destination_manifest: dict[str, dict[str, object]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for relative in PRODUCER_PATHS:
        source_record = source_manifest.get(relative)
        destination_record = destination_manifest.get(relative)
        if source_record is None or destination_record is None:
            raise SystemExit(f"Recount producer is absent from one source snapshot: {relative}")
        if source_record != destination_record:
            raise SystemExit(f"Recount producer contract differs across runs: {relative}")
        source_path = source_root / "source_snapshot" / relative
        destination_path = destination_root / "source_snapshot" / relative
        for path, label in ((source_path, "execution"), (destination_path, "analysis")):
            require_regular(path, f"{label} producer")
            if path.stat().st_size != source_record["size"] or sha256(path) != source_record["sha256"]:
                raise SystemExit(f"{label} producer bytes differ from its frozen manifest: {relative}")
        rows.append({"relative_path": relative, "sha256": str(source_record["sha256"])})
    return rows


def compare_regression_checker(
    source_root: Path,
    destination_root: Path,
    source_manifest: dict[str, dict[str, object]],
    destination_manifest: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Validate the one whitelisted non-counting checker transition.

    This checker never constructs a BAM manifest, invokes featureCounts, merges
    counts, or validates count artifacts. Its only accepted byte transition is
    the reviewed 055119Z preamble scanner to the guarded-route checker.
    All actual recount producer bytes remain exact in compare_producer_manifests.
    """
    source_record = source_manifest.get(REGRESSION_CHECKER_REL)
    destination_record = destination_manifest.get(REGRESSION_CHECKER_REL)
    if source_record is None or destination_record is None:
        raise SystemExit("Source regression checker is absent from one frozen snapshot")
    by_digest = {
        value["sha256"]: (name, value["size"])
        for name, value in REGRESSION_CHECKER_DIGESTS.items()
    }
    source_info = by_digest.get(str(source_record.get("sha256")))
    destination_info = by_digest.get(str(destination_record.get("sha256")))
    if source_info is None or destination_info is None:
        raise SystemExit("Source regression checker digest is not in the narrow transition whitelist")
    if source_record.get("size") != source_info[1] or destination_record.get("size") != destination_info[1]:
        raise SystemExit("Source regression checker size differs from its whitelisted digest")
    if destination_info[0] != "guarded_routes_v1":
        raise SystemExit("Destination must use the reviewed guarded-route regression checker")
    if source_info[0] not in ("legacy_055119Z", "guarded_routes_v1"):
        raise SystemExit("Unsupported source regression-checker generation")

    source_path = source_root / "source_snapshot" / REGRESSION_CHECKER_REL
    destination_path = destination_root / "source_snapshot" / REGRESSION_CHECKER_REL
    for path, record, label in (
        (source_path, source_record, "execution-source"),
        (destination_path, destination_record, "destination"),
    ):
        require_regular(path, f"{label} source regression checker")
        if path.stat().st_size != record["size"] or sha256(path) != record["sha256"]:
            raise SystemExit(f"{label} source regression checker differs from its frozen manifest")

    source_snapshot = source_root / "source_snapshot"
    destination_snapshot = destination_root / "source_snapshot"
    source_checker_on_source = run_checked(
        [sys.executable, str(source_path), "--project-root", str(source_snapshot)],
        label="execution checker on execution snapshot",
    )
    destination_checker_on_source = run_recorded(
        [sys.executable, str(destination_path), "--project-root", str(source_snapshot)],
    )
    destination_checker_on_destination = run_checked(
        [sys.executable, str(destination_path), "--project-root", str(destination_snapshot)],
        label="destination checker on destination snapshot",
    )
    expected_source_output = (
        LEGACY_REGRESSION_CHECKER_EXPECTED_OUTPUT
        if source_info[0] == "legacy_055119Z"
        else GUARDED_REGRESSION_CHECKER_EXPECTED_OUTPUT
    )
    if source_checker_on_source != expected_source_output:
        raise SystemExit("Execution checker/execution snapshot did not reproduce its exact producer census")
    if destination_checker_on_destination != GUARDED_REGRESSION_CHECKER_EXPECTED_OUTPUT:
        raise SystemExit("Destination checker/destination snapshot did not reproduce exact 14/5/1/2 census")
    if source_info[0] == "legacy_055119Z":
        observed_cross = {
            key: destination_checker_on_source[key]
            for key in GUARDED_CHECKER_ON_LEGACY_EXPECTED
        }
        if observed_cross != GUARDED_CHECKER_ON_LEGACY_EXPECTED:
            raise SystemExit(
                "Guarded destination checker on legacy snapshot differs from the reviewed exact failure"
            )
    elif (
        destination_checker_on_source["returncode"] != 0
        or destination_checker_on_source["stdout"] + destination_checker_on_source["stderr"]
        != GUARDED_REGRESSION_CHECKER_EXPECTED_OUTPUT
    ):
        raise SystemExit("Guarded checker did not reproduce the guarded source snapshot")
    diff = "".join(difflib.unified_diff(
        source_path.read_text().splitlines(keepends=True),
        destination_path.read_text().splitlines(keepends=True),
        fromfile=f"{source_root.name}/{REGRESSION_CHECKER_REL}",
        tofile=f"{destination_root.name}/{REGRESSION_CHECKER_REL}",
    ))
    return {
        "relative_path": REGRESSION_CHECKER_REL,
        "execution_generation": source_info[0],
        "destination_generation": destination_info[0],
        "execution_sha256": source_record["sha256"],
        "destination_sha256": destination_record["sha256"],
        "execution_size": source_record["size"],
        "destination_size": destination_record["size"],
        "execution_checker_on_execution_snapshot": source_checker_on_source,
        "destination_checker_on_execution_snapshot": destination_checker_on_source,
        "destination_checker_on_destination_snapshot": destination_checker_on_destination,
        "unified_diff": diff,
    }


def parse_digest_manifest(root: Path, path: Path) -> dict[str, str]:
    require_regular(path, "digest manifest")
    records: dict[str, str] = {}
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        fields = line.split("\t", 1)
        if len(fields) != 2 or not SHA_RE.fullmatch(fields[0]):
            raise SystemExit(f"Malformed digest manifest row {number}: {path}")
        relative = fields[1]
        rel = Path(relative)
        if rel.is_absolute() or ".." in rel.parts or relative in records:
            raise SystemExit(f"Escaping or duplicate digest path: {relative}")
        artifact = (root / rel).resolve(strict=True)
        if not contained(artifact, root) or artifact.is_symlink() or not artifact.is_file():
            raise SystemExit(f"Digest artifact escapes or is nonregular: {artifact}")
        if sha256(artifact) != fields[0]:
            raise SystemExit(f"Digest artifact changed: {artifact}")
        records[relative] = fields[0]
    if not records:
        raise SystemExit(f"Empty digest manifest: {path}")
    return records


def read_digest_records(path: Path) -> dict[str, str]:
    """Parse a relative digest manifest without assuming its artifacts are colocated."""
    require_regular(path, "archived digest manifest")
    records: dict[str, str] = {}
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        fields = line.split("\t", 1)
        if len(fields) != 2 or not SHA_RE.fullmatch(fields[0]):
            raise SystemExit(f"Malformed archived digest row {number}: {path}")
        relative = fields[1]
        rel = Path(relative)
        if rel.is_absolute() or ".." in rel.parts or relative in records:
            raise SystemExit(f"Escaping or duplicate archived digest path: {relative}")
        records[relative] = fields[0]
    if not records:
        raise SystemExit(f"Empty archived digest manifest: {path}")
    return records


def parse_absolute_hash_manifest(root: Path, relative: str) -> dict[str, str]:
    path = root / relative
    require_regular(path, relative)
    records: dict[str, str] = {}
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        match = re.fullmatch(r"([0-9a-f]{64})\s+(.+)", line)
        if not match:
            raise SystemExit(f"Malformed baseline digest row {number}: {path}")
        artifact = Path(match.group(2)).resolve(strict=True)
        if not contained(artifact, root) or artifact.is_symlink() or not artifact.is_file():
            raise SystemExit(f"Baseline digest path escapes or is nonregular: {artifact}")
        rel = str(artifact.relative_to(root))
        if rel in records or sha256(artifact) != match.group(1):
            raise SystemExit(f"Duplicate or changed baseline artifact: {artifact}")
        records[rel] = match.group(1)
    if not records:
        raise SystemExit(f"Empty baseline digest manifest: {path}")
    return records


def baseline_generation(root: Path) -> str:
    """Validate and classify the legacy or structured baseline transaction.

    A zero-byte BASELINE_FROZEN marker is accepted only here for a historical
    execution source. Fresh candidates must have the structured JSON seal plus
    a nonempty compatibility marker that is exactly the JSON SHA-256.
    """
    legacy = root / "BASELINE_FROZEN"
    structured = root / "BASELINE_FROZEN.json"
    if structured.exists() or structured.is_symlink():
        require_regular(structured, "BASELINE_FROZEN.json")
        require_regular(legacy, "BASELINE_FROZEN compatibility marker")
        expected = sha256(structured) + "\n"
        try:
            observed = legacy.read_text(encoding="ascii")
        except UnicodeError as exc:
            raise SystemExit("BASELINE_FROZEN compatibility marker is not ASCII") from exc
        if observed != expected:
            raise SystemExit("Fresh BASELINE_FROZEN marker does not bind BASELINE_FROZEN.json")
        marker = load_json(structured, "structured baseline marker")
        if marker.get("schema") != "bg001-baseline-seal-v1" or marker.get("run_id") != root.name:
            raise SystemExit("Structured baseline marker schema/run ID differs")
        return "structured_v1"
    if structured.is_symlink():
        raise SystemExit("Legacy source has a symlinked structured baseline marker")
    require_regular(legacy, "legacy BASELINE_FROZEN", nonempty=False)
    if legacy.stat().st_size != 0:
        raise SystemExit("Legacy BASELINE_FROZEN must be the historical zero-byte marker")
    unexpected_runtime = [
        relative for relative in sorted(RUNTIME_BASELINE_INPUTS)
        if (root / relative).exists() or (root / relative).is_symlink()
    ]
    if unexpected_runtime:
        raise SystemExit(
            "Legacy baseline contains a partial/unsealed runtime record: "
            + ",".join(unexpected_runtime)
        )
    unexpected_guard = [
        relative for relative in sorted(BAM_GUARD_BASELINE_INPUTS)
        if (root / relative).exists() or (root / relative).is_symlink()
    ]
    if unexpected_guard:
        raise SystemExit(
            "Legacy execution source contains unsealed/synthesized BAM guard artifacts: "
            + ",".join(unexpected_guard)
        )
    return "legacy_zero_marker"


def verify_baseline(root: Path) -> tuple[dict[str, str], str]:
    generation = baseline_generation(root)
    records: dict[str, str] = {}
    for relative in (
        "manifests/baseline_frozen_files.sha256",
        "manifests/baseline_manifests.sha256",
    ):
        for artifact, digest in parse_absolute_hash_manifest(root, relative).items():
            if artifact in records:
                raise SystemExit(f"Baseline artifact is sealed twice: {artifact}")
            records[artifact] = digest
    for required in (
        PROTECTED_GENE_REL,
        PROTECTED_SOURCE_MANIFEST_REL,
        "frozen_sets/locked_merged_dge.rds",
        "frozen_sets/locked_meta_matched.rds",
        "frozen_sets/locked_sample_qc_report.csv",
        "frozen_sets/canonical_deg_results.csv",
        "frozen_sets/canonical_samples.tsv",
        "frozen_sets/canonical_genes.tsv",
        "frozen_sets/finite_susie_447.tsv",
        "manifests/read_count_overrides.tsv",
    ):
        if required not in records:
            raise SystemExit(f"Baseline seal lacks required artifact: {required}")
    if generation == "structured_v1":
        required_fresh = (
            RUNTIME_BASELINE_INPUTS
            | FRESH_CONTROL_BASELINE_INPUTS
            | {NOT_TESTABLE_REL, RENDERED_VALIDATION_REL, SUBSTRATE_CENSUS_REL}
        )
        missing = required_fresh - set(records)
        if missing:
            raise SystemExit(f"Structured baseline seal lacks fresh inputs: {sorted(missing)}")
        validation = root / RENDERED_VALIDATION_REL
        if not validation.read_text(errors="strict").startswith("exit_status=0\n"):
            raise SystemExit("Frozen rendered-label validation did not record exit_status=0")
        with (root / SUBSTRATE_CENSUS_REL).open(newline="") as handle:
            census_rows = list(csv.DictReader(handle, delimiter="\t"))
        expected_census = {
            "matched_samples": "1281",
            "pass_technical_samples": "1260",
            "dge_samples": "1260",
            "canonical_model_samples": "846",
        }
        if len(census_rows) != 1 or census_rows[0] != expected_census:
            raise SystemExit("Structured baseline has an unexpected integration-substrate census")
        guard_script = root / "source_snapshot" / REMEDIATION_REL / "bam_manifest_finalization_guard.py"
        run_checked(
            [sys.executable, str(guard_script), "verify", "--run-root", str(root)],
            label="structured baseline BAM-manifest guard verification",
        )
    elif set(records) & (RUNTIME_BASELINE_INPUTS | FRESH_CONTROL_BASELINE_INPUTS):
        raise SystemExit("Legacy baseline unexpectedly seals fresh runtime/control inputs")
    return records, generation


def normalized_protected_sources(root: Path) -> dict[str, str]:
    project = Path(str(load_json(root / "contract/run_contract.json", "run contract")["project_root"])).resolve()
    source_snapshot = (root / "source_snapshot").resolve(strict=True)
    frozen_root = (root / "frozen_sets/protected_sources").resolve(strict=True)
    path = root / PROTECTED_SOURCE_MANIFEST_REL
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("source_path", "snapshot_path", "frozen_path"):
            raise SystemExit(f"Unexpected protected-source manifest schema: {path}")
        rows = list(reader)
    normalized: dict[str, str] = {}
    for row in rows:
        # source_path is provenance, not a live dependency.  It may legitimately
        # disappear after the source run, so normalize it lexically while the
        # frozen snapshot/copy below remain strict, hash-checked dependencies.
        source = Path(row["source_path"])
        snapshot = Path(row["snapshot_path"]).resolve(strict=True)
        frozen = Path(row["frozen_path"]).resolve(strict=True)
        if (
            not source.is_absolute()
            or ".." in source.parts
            or not contained(source, project)
            or not contained(snapshot, source_snapshot)
            or not contained(frozen, frozen_root)
        ):
            raise SystemExit("Protected-source manifest path escapes its declared root")
        source_rel = str(source.relative_to(project))
        if str(snapshot.relative_to(source_snapshot)) != source_rel or str(frozen.relative_to(frozen_root)) != source_rel:
            raise SystemExit(f"Protected-source path mapping differs: {source_rel}")
        if source_rel in normalized or any(path.is_symlink() for path in (snapshot, frozen)):
            raise SystemExit(f"Duplicate or symlinked protected source: {source_rel}")
        snapshot_sha = sha256(snapshot)
        if sha256(frozen) != snapshot_sha:
            raise SystemExit(f"Frozen protected source differs from its snapshot: {source_rel}")
        normalized[source_rel] = snapshot_sha
    observed = {
        str(path.relative_to(frozen_root))
        for path in frozen_root.rglob("*")
        if path.is_file() and not path.is_symlink()
    }
    if observed != set(normalized):
        raise SystemExit("Protected-source manifest and frozen protected-source tree differ")
    return normalized


def protected_gene_rows(root: Path) -> dict[str, dict[str, object]]:
    path = root / PROTECTED_GENE_REL
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("gene", "symbol", "protected_reason", "protected_source"):
            raise SystemExit(f"Unexpected protected-gene schema: {path}")
        rows = list(reader)
    result: dict[str, dict[str, object]] = {}
    for row in rows:
        gene = row["gene"]
        reasons = {value for value in row["protected_reason"].split("|") if value}
        sources = {value for value in row["protected_source"].split("|") if value}
        if not gene or gene in result or not row["symbol"] or not reasons or not sources:
            raise SystemExit(f"Malformed or duplicate protected-gene row: {gene!r}")
        result[gene] = {"symbol": row["symbol"], "reasons": reasons, "sources": sources}
    return result


def not_testable_rows(root: Path) -> set[tuple[str, str, str]]:
    path = root / NOT_TESTABLE_REL
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = ("symbol", "protected_reason", "protected_source")
        if tuple(reader.fieldnames or ()) != expected:
            raise SystemExit(f"Unexpected not-testable protected-gene schema: {path}")
        rows = list(reader)
    result: set[tuple[str, str, str]] = set()
    for row in rows:
        record = tuple(row[field] for field in expected)
        if any(not value for value in record) or record in result:
            raise SystemExit(f"Malformed or duplicate not-testable protected-gene row: {record}")
        result.add(record)
    return result


def normalized_read_count_manifests(root: Path) -> dict[str, list[tuple[str, ...]]]:
    sources_path = root / "manifests/read_count_sources.tsv"
    with sources_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_fields = (
            "dataset", "source_count_path", "frozen_count_path",
            "source_summary_path", "frozen_summary_path",
        )
        if tuple(reader.fieldnames or ()) != expected_fields:
            raise SystemExit("Unexpected read_count_sources.tsv schema")
        source_rows = list(reader)
    normalized_sources: list[tuple[str, ...]] = []
    for row in source_rows:
        frozen_count = Path(row["frozen_count_path"]).resolve(strict=True)
        frozen_summary = Path(row["frozen_summary_path"]).resolve(strict=True)
        if not contained(frozen_count, root) or not contained(frozen_summary, root):
            raise SystemExit("Frozen read-count source path escapes its run")
        normalized_sources.append(
            (
                row["dataset"],
                row["source_count_path"],
                str(frozen_count.relative_to(root)),
                row["source_summary_path"],
                str(frozen_summary.relative_to(root)),
            )
        )

    overrides_path = root / "manifests/read_count_overrides.tsv"
    with overrides_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("dataset", "count_path"):
            raise SystemExit("Unexpected read_count_overrides.tsv schema")
        override_rows = list(reader)
    normalized_overrides: list[tuple[str, ...]] = []
    for row in override_rows:
        count = Path(row["count_path"]).resolve(strict=True)
        if not contained(count, root):
            raise SystemExit("Read-count override path escapes its run")
        normalized_overrides.append((row["dataset"], str(count.relative_to(root))))
    if (
        len(normalized_sources) != 9
        or len({row[0] for row in normalized_sources}) != 9
        or len(normalized_overrides) != 9
        or len({row[0] for row in normalized_overrides}) != 9
    ):
        raise SystemExit("Read-count baseline manifests do not contain exactly nine datasets")
    return {
        "read_count_sources": sorted(normalized_sources),
        "read_count_overrides": sorted(normalized_overrides),
    }


def require_monotonic_protected_genes(
    source_genes: dict[str, dict[str, object]],
    destination_genes: dict[str, dict[str, object]],
) -> None:
    removed_genes = set(source_genes) - set(destination_genes)
    if removed_genes:
        raise SystemExit(f"Destination removed protected genes: {sorted(removed_genes)[:5]}")
    for gene, old in source_genes.items():
        new = destination_genes[gene]
        if (
            old["symbol"] != new["symbol"]
            or not old["reasons"].issubset(new["reasons"])
            or not old["sources"].issubset(new["sources"])
        ):
            raise SystemExit(f"Destination protection is not monotonic for {gene}")


def compare_baselines(source_root: Path, destination_root: Path) -> list[dict[str, str]]:
    source, source_generation = verify_baseline(source_root)
    destination, destination_generation = verify_baseline(destination_root)
    if destination_generation != "structured_v1":
        raise SystemExit("Recount import destination must use the structured fresh baseline seal")

    source_runtime = set(source) & RUNTIME_BASELINE_INPUTS
    destination_runtime = set(destination) & RUNTIME_BASELINE_INPUTS
    if source_runtime not in (set(), RUNTIME_BASELINE_INPUTS):
        raise SystemExit("Execution-source baseline has a partial analysis runtime sidecar set")
    if destination_runtime != RUNTIME_BASELINE_INPUTS:
        raise SystemExit("Destination baseline lacks the exact analysis runtime sidecar set")
    if source_runtime == RUNTIME_BASELINE_INPUTS:
        changed_runtime = [
            relative for relative in sorted(RUNTIME_BASELINE_INPUTS)
            if source[relative] != destination[relative]
        ]
        if changed_runtime:
            raise SystemExit(
                "Fresh-to-fresh analysis runtime sidecars are not byte-identical: "
                + ",".join(changed_runtime)
            )

    allowed = {
        PROTECTED_GENE_REL,
        NOT_TESTABLE_REL,
        RENDERED_VALIDATION_REL,
        SUBSTRATE_CENSUS_REL,
        PROTECTED_SOURCE_MANIFEST_REL,
    }
    differences: list[dict[str, str]] = []
    for relative in sorted(set(source) | set(destination)):
        left, right = source.get(relative), destination.get(relative)
        if left == right:
            continue
        if relative in RUNTIME_BASELINE_INPUTS:
            if source_runtime:
                # Any fresh-to-fresh difference was rejected above.
                raise SystemExit(f"Analysis runtime sidecar changed: {relative}")
            if left is not None or right is None:
                raise SystemExit(f"Invalid legacy-to-fresh runtime transition: {relative}")
            change_class = "analysis_runtime_added"
        elif relative in FRESH_CONTROL_BASELINE_INPUTS:
            if right is None:
                raise SystemExit(f"Destination removed a structured baseline control input: {relative}")
            if relative in BAM_GUARD_BASELINE_INPUTS and left is None:
                change_class = "legacy_source_no_guard_to_guarded_destination"
            elif relative in BAM_GUARD_BASELINE_INPUTS:
                change_class = "guarded_source_to_guarded_destination_rebind"
            else:
                change_class = "structured_control_rebind"
        else:
            protected = relative in allowed or relative.startswith(PROTECTED_SOURCE_PREFIX)
            run_bound = relative in RUN_BOUND_BASELINE_MANIFESTS
            if not protected and not run_bound:
                raise SystemExit(f"Locked non-protected baseline artifact differs: {relative}")
            change_class = (
                "protected_output_expansion_or_refresh" if protected else "run_path_rebase"
            )
        differences.append(
            {
                "relative_path": relative,
                "execution_source_sha256": left or "",
                "destination_analysis_sha256": right or "",
                "change_class": change_class,
            }
        )

    if normalized_read_count_manifests(source_root) != normalized_read_count_manifests(destination_root):
        raise SystemExit("Run-path-normalized read-count baseline manifests differ")

    source_protected_sources = normalized_protected_sources(source_root)
    destination_protected_sources = normalized_protected_sources(destination_root)
    removed_sources = set(source_protected_sources) - set(destination_protected_sources)
    if removed_sources:
        raise SystemExit(f"Destination removed protected-source inputs: {sorted(removed_sources)[:5]}")

    source_genes = protected_gene_rows(source_root)
    destination_genes = protected_gene_rows(destination_root)
    require_monotonic_protected_genes(source_genes, destination_genes)
    source_not_testable = not_testable_rows(source_root)
    destination_not_testable = not_testable_rows(destination_root)
    removed_not_testable = source_not_testable - destination_not_testable
    if removed_not_testable:
        raise SystemExit(
            "Destination removed not-testable protected-gene evidence: "
            f"{sorted(removed_not_testable)[:5]}"
        )
    return differences


def compare_bam_freezes(source_root: Path, destination_root: Path) -> dict[str, object]:
    source_manifest = source_root / "manifests/bam_manifest.tsv"
    destination_manifest = destination_root / "manifests/bam_manifest.tsv"
    source_sha, destination_sha = sha256(source_manifest), sha256(destination_manifest)
    if source_sha != destination_sha:
        raise SystemExit("Destination frozen BAM manifest is not byte-identical to the execution manifest")
    if (
        (source_root / "BAM_MANIFEST_FROZEN").read_text().strip() != source_sha
        or (destination_root / "BAM_MANIFEST_FROZEN").read_text().strip() != source_sha
    ):
        raise SystemExit("One BAM_MANIFEST_FROZEN marker differs from the shared manifest")
    source_freeze = load_json(source_root / "manifests/bam_manifest.freeze.json", "source BAM freeze")
    destination_freeze = load_json(destination_root / "manifests/bam_manifest.freeze.json", "destination BAM freeze")
    for key in (
        "schema_version", "draft_sha256", "manifest_sha256", "included_bams",
        "documented_unavailable", "shard_sha256",
    ):
        if source_freeze.get(key) != destination_freeze.get(key):
            raise SystemExit(f"Destination BAM-freeze contract differs: {key}")
    if source_freeze.get("run_id") != source_root.name or destination_freeze.get("run_id") != destination_root.name:
        raise SystemExit("BAM-freeze metadata does not name its candidate")
    for relative in ("manifests/bam_manifest.draft.tsv",) + tuple(
        f"manifests/bam_hash_shards/task_{index}.tsv" for index in range(8)
    ):
        if sha256(source_root / relative) != sha256(destination_root / relative):
            raise SystemExit(f"Destination independently frozen BAM input differs: {relative}")
    return {
        "bam_manifest_sha256": source_sha,
        "draft_sha256": source_freeze["draft_sha256"],
        "shard_sha256": source_freeze["shard_sha256"],
        "included_bams": source_freeze["included_bams"],
        "documented_unavailable": source_freeze["documented_unavailable"],
    }


def environment_values(path: Path, key: str) -> set[str]:
    result: set[str] = set()
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split("\t", 1)
        if len(fields) == 2 and fields[0] == key:
            result.add(fields[1])
    return result


def load_tool_pins(source_root: Path) -> dict[str, str]:
    path = source_root / "manifests/count_manifest.tsv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 5:
        raise SystemExit("Source count manifest does not have exactly five cohorts")
    fields = (
        "pinned_environment_sha256",
        "pinned_featurecounts_executable_sha256",
        "pinned_samtools_executable_sha256",
        "pinned_samtools_version",
        "source_manifest_sha256",
    )
    pins: dict[str, str] = {}
    for field in fields:
        values = {row.get(field, "") for row in rows}
        if len(values) != 1 or not next(iter(values)):
            raise SystemExit(f"Source count manifest does not pin one {field}")
        pins[field] = next(iter(values))
    for field in fields[:3] + ("source_manifest_sha256",):
        if not SHA_RE.fullmatch(pins[field]):
            raise SystemExit(f"Invalid source count-manifest digest: {field}")

    environment = source_root / "counts/GSE130970/environment.txt"
    for key in (
        "featureCounts_version", "featureCounts_executable", "featureCounts_executable_sha256",
        "samtools_version", "samtools_executable", "samtools_executable_sha256",
        "gtf_path", "gtf_sha256",
    ):
        values = environment_values(environment, key)
        if len(values) != 1:
            raise SystemExit(f"Source environment does not contain exactly one {key}")
        pins[key] = next(iter(values))
    if (
        pins["featureCounts_version"] != "featureCounts v2.1.1"
        or pins["gtf_sha256"] != EXPECTED_GTF_SHA256
        or pins["pinned_featurecounts_executable_sha256"] != pins["featureCounts_executable_sha256"]
        or pins["pinned_samtools_executable_sha256"] != pins["samtools_executable_sha256"]
        or pins["pinned_samtools_version"] != pins["samtools_version"]
        or pins["pinned_environment_sha256"] != sha256(environment)
    ):
        raise SystemExit("Source count-manifest tool pins differ from task environment")
    for executable_key, digest_key in (
        ("featureCounts_executable", "featureCounts_executable_sha256"),
        ("samtools_executable", "samtools_executable_sha256"),
    ):
        executable = Path(pins[executable_key]).resolve(strict=True)
        require_regular(executable, executable_key)
        if sha256(executable) != pins[digest_key]:
            raise SystemExit(f"Pinned tool binary changed after source recount: {executable}")
    gtf = Path(pins["gtf_path"]).resolve(strict=True)
    require_regular(gtf, "GENCODE v49 GTF")
    if sha256(gtf) != EXPECTED_GTF_SHA256:
        raise SystemExit("GENCODE v49 GTF bytes differ from the recount contract")
    return pins


def preflight_destination_empty(root: Path) -> None:
    counts = root / "counts"
    if not counts.is_dir() or counts.is_symlink() or any(counts.iterdir()):
        raise SystemExit("Destination counts directory is missing, symlinked, or nonempty")
    forbidden = (
        "RECOUNT_COMPLETE", "RECOUNT_IMPORTED", "RECOUNT_IMPORT_IN_PROGRESS",
        "RECOUNT_IMPORT_FAILED.json", "ANALYSIS_COMPLETE",
        "manifests/count_manifest.tsv", "manifests/count_overrides.tsv",
        "manifests/recount_artifacts.sha256", "manifests/recount_import",
        "comparisons/structural_validation.json",
    )
    for relative in forbidden:
        path = root / relative
        if path.exists() or path.is_symlink():
            raise SystemExit(f"Destination already contains recount/import state: {relative}")
    for arm in ("R0", "F_locked", "F_legacy", "F_five"):
        arm_root = root / "arms" / arm
        if not arm_root.is_dir() or arm_root.is_symlink():
            raise SystemExit(f"Destination arm is missing or symlinked: {arm}")
        files = [path for path in arm_root.rglob("*") if path.is_file() and path.name != ".bg001_candidate_root"]
        if files:
            raise SystemExit(f"Destination analysis arm was already attempted: {files[0]}")


def source_state(root: Path, records: dict[str, str]) -> dict[str, tuple[int, int, int, int, str]]:
    result: dict[str, tuple[int, int, int, int, str]] = {}
    for relative, digest in records.items():
        path = root / relative
        stat = path.stat()
        result[relative] = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, digest)
    return result


def write_tsv(path: Path, fieldnames: tuple[str, ...], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        raise SystemExit(f"Refusing to overwrite existing candidate TSV: {path}")
    temporary = path.with_name(f".{path.name}.tmp.{uuid.uuid4().hex}")
    with temporary.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    if path.exists() or path.is_symlink():
        temporary.unlink()
        raise SystemExit(f"Candidate TSV appeared during atomic write: {path}")
    os.replace(temporary, path)
    fsync_directory(path.parent)


def write_rebased_provenance(stage_task: Path, final_task: Path) -> None:
    lines = []
    for name in PROVENANCE_PAYLOADS:
        path = stage_task / name
        require_regular(path, f"copied task payload {name}")
        lines.append(f"{sha256(path)}  {final_task / name}\n")
    atomic_write_text(
        stage_task / "provenance.sha256",
        "".join(lines),
        replace_existing=True,
    )


def write_rebased_merge_publication(stage_parent: Path, destination_run_id: str) -> None:
    hashes = {}
    for name in MERGE_PAYLOADS:
        path = stage_parent / name
        require_regular(path, f"copied merged payload {name}")
        hashes[name] = sha256(path)
    payload = {
        "status": "VALIDATED_STAGING",
        "run_id": destination_run_id,
        "publish_order": list(MERGE_PUBLISH_ORDER),
        "hashes": hashes,
    }
    publication = stage_parent / "merge_publication.json"
    atomic_write_json(publication, payload, replace_existing=True)
    atomic_write_text(
        stage_parent / "MERGE_COMPLETE",
        sha256(publication) + "\n",
        replace_existing=True,
    )


def copy_and_rebase_counts(
    source_root: Path,
    destination_root: Path,
    sealed: dict[str, str],
    staging: Path,
) -> tuple[Path, Path, list[dict[str, str]]]:
    stage_counts = staging / "counts"
    stage_lineage = staging / "lineage"
    originals = stage_lineage / "source_originals"
    stage_counts.mkdir(parents=True)
    originals.mkdir(parents=True)
    inventory: list[dict[str, str]] = []

    for relative, digest in sorted(sealed.items()):
        source = source_root / relative
        if relative.startswith("counts/"):
            target = staging / relative
            disposition = "copied_exact"
        else:
            target = originals / relative
            disposition = "archived_original"
        copy_independent_regular(source, target, digest, f"source recount artifact {relative}")
        inventory.append(
            {
                "source_relative_path": relative,
                "source_sha256": digest,
                "disposition": disposition,
                "preserved_destination_relative_path": relative
                if relative.startswith("counts/")
                else str(
                    Path("manifests/recount_import/source_originals")
                    / target.relative_to(originals)
                ),
                "effective_destination_relative_path": relative if relative.startswith("counts/") else "",
                "effective_sha256": digest if relative.startswith("counts/") else "",
            }
        )

    inventory_by_path = {row["source_relative_path"]: row for row in inventory}
    if set(PATH_BOUND_COUNT_FILES) - set(inventory_by_path):
        raise SystemExit("Source seal lacks path-bound recount metadata")
    for relative in sorted(PATH_BOUND_COUNT_FILES):
        source = source_root / relative
        archive = originals / relative
        copy_independent_regular(
            source,
            archive,
            sealed[relative],
            f"source path-bound artifact {relative}",
        )
        row = inventory_by_path[relative]
        row["disposition"] = "archived_original_and_rebased"
        row["preserved_destination_relative_path"] = str(
            Path("manifests/recount_import/source_originals") / relative
        )

    task_relatives = [
        f"counts/{dataset}" for dataset in ("GSE130970", "GSE135251", "GSE174478", "GSE240729")
    ] + [f"counts/GSE213621/chunks/chunk{index}" for index in range(4)]
    for relative in task_relatives:
        write_rebased_provenance(staging / relative, destination_root / relative)
        inventory_by_path[f"{relative}/provenance.sha256"]["effective_sha256"] = sha256(
            staging / relative / "provenance.sha256"
        )
    merged = stage_counts / "GSE213621"
    write_rebased_merge_publication(merged, destination_root.name)
    for name in ("merge_publication.json", "MERGE_COMPLETE"):
        inventory_by_path[f"counts/GSE213621/{name}"]["effective_sha256"] = sha256(merged / name)

    observed = {
        str(path.relative_to(staging))
        for path in stage_counts.rglob("*")
        if path.is_file() and not path.is_symlink()
    }
    expected = {relative for relative in sealed if relative.startswith("counts/")}
    if observed != expected:
        raise SystemExit("Staged count artifact set differs from the source sealed count set")
    for relative in expected - PATH_BOUND_COUNT_FILES:
        if sha256(staging / relative) != sealed[relative]:
            raise SystemExit(f"Staged count artifact changed: {relative}")
    return stage_counts, stage_lineage, inventory


def archive_extra_source_lineage(
    source_root: Path,
    stage_lineage: Path,
) -> None:
    extras = (
        "manifests/recount_artifacts.sha256",
        "comparisons/structural_validation.json",
        "RECOUNT_COMPLETE",
        "contract/run_contract.json",
        "contract/source_manifest.tsv",
        "contract/source_regression_check.txt",
    )
    for relative in extras:
        source = source_root / relative
        require_regular(source, f"source lineage artifact {relative}", nonempty=relative != "RECOUNT_COMPLETE")
        target = stage_lineage / "source_originals" / relative
        if target.exists() or target.is_symlink():
            if sha256(target) != sha256(source):
                raise SystemExit(f"Conflicting archived lineage artifact: {relative}")
            continue
        expected = sha256(source)
        copy_independent_regular(source, target, expected, f"source lineage artifact {relative}")
    for name in PRODUCER_FILENAMES + POSTCOUNT_METADATA_FILENAMES:
        relative = str(Path("source_snapshot") / REMEDIATION_REL / name)
        source = source_root / relative
        require_regular(source, f"source lineage code {name}")
        target = stage_lineage / "source_originals" / relative
        expected = sha256(source)
        copy_independent_regular(source, target, expected, f"source lineage code {name}")


def archive_scheduler_attestation(
    source_root: Path,
    stage_lineage: Path,
    scheduler_attestation: dict[str, object],
) -> dict[str, object]:
    raw_sacct = scheduler_attestation.get("raw_sacct_psv")
    validation = scheduler_attestation.get("validation")
    if not isinstance(raw_sacct, str) or not raw_sacct:
        raise SystemExit("Source scheduler attestation lacks raw sacct evidence")
    if not isinstance(validation, dict) or validation.get("status") != "PASS":
        raise SystemExit("Source scheduler attestation is not PASS")
    ledger = source_root / "contract/slurm_submission.tsv"
    expected_ledger_sha256 = validation.get("ledger_sha256")
    expected_raw_sha256 = validation.get("raw_sacct_sha256")
    if sha256(ledger) != expected_ledger_sha256:
        raise SystemExit("Source SLURM submission ledger changed after attestation")
    if sha256_bytes(raw_sacct.encode()) != expected_raw_sha256:
        raise SystemExit("Raw sacct evidence differs from its attestation digest")

    scheduler_root = stage_lineage / "slurm"
    scheduler_root.mkdir()
    archived_ledger = scheduler_root / "source_submission_ledger.tsv"
    archived_raw = scheduler_root / "sacct.psv"
    archived_validation = scheduler_root / "validation.json"
    copy_independent_regular(
        ledger,
        archived_ledger,
        str(expected_ledger_sha256),
        "source SLURM submission ledger",
    )
    atomic_write_text(archived_raw, raw_sacct)
    atomic_write_json(archived_validation, validation)
    if sha256(archived_raw) != expected_raw_sha256:
        raise SystemExit("Archived raw sacct evidence differs from its captured digest")
    return {
        "source_submission_ledger_sha256": sha256(archived_ledger),
        "raw_sacct_sha256": sha256(archived_raw),
        "validation_sha256": sha256(archived_validation),
        "logical_jobs": validation.get("logical_jobs"),
        "completed_jobs": validation.get("completed_jobs"),
        "cancelled_unstarted_jobs": validation.get("cancelled_unstarted_jobs"),
        "observed_array_max_concurrency": validation.get(
            "observed_array_max_concurrency"
        ),
        "schema_version": validation.get("schema_version"),
    }


def archive_regression_checker_transition(
    source_root: Path,
    destination_root: Path,
    stage_lineage: Path,
    transition: dict[str, object],
) -> None:
    directory = stage_lineage / "regression_checker_transition"
    directory.mkdir()
    source = source_root / "source_snapshot" / REGRESSION_CHECKER_REL
    destination = destination_root / "source_snapshot" / REGRESSION_CHECKER_REL
    copy_independent_regular(
        source,
        directory / "execution_checker.py",
        str(transition["execution_sha256"]),
        "execution-source regression checker",
    )
    copy_independent_regular(
        destination,
        directory / "destination_checker.py",
        str(transition["destination_sha256"]),
        "destination regression checker",
    )
    for key, filename in (
        ("execution_checker_on_execution_snapshot", "execution_checker_on_execution_snapshot.txt"),
        ("destination_checker_on_destination_snapshot", "destination_checker_on_destination_snapshot.txt"),
        ("unified_diff", "checker_transition.diff"),
    ):
        atomic_write_text(directory / filename, str(transition[key]))
    atomic_write_json(
        directory / "destination_checker_on_execution_snapshot.json",
        transition["destination_checker_on_execution_snapshot"],
    )
    metadata = {
        key: transition[key]
        for key in (
            "relative_path", "execution_generation", "destination_generation",
            "execution_sha256", "destination_sha256", "execution_size", "destination_size",
        )
    }
    metadata["validated_producer_census"] = {
        "execution_source": {
            "fixed_human": 19 if transition["execution_generation"] == "legacy_055119Z" else 14,
            "guarded_custom": 0 if transition["execution_generation"] == "legacy_055119Z" else 5,
            "guarded_dynamic_shared": 0 if transition["execution_generation"] == "legacy_055119Z" else 1,
            "dynamic_shared": 1 if transition["execution_generation"] == "legacy_055119Z" else 0,
            "mouse": 2,
        },
        "destination": {
            "fixed_human": 14,
            "guarded_custom": 5,
            "guarded_dynamic_shared": 1,
            "mouse": 2,
        },
    }
    atomic_write_json(directory / "transition.json", metadata)


def publish_staging(
    destination_root: Path,
    staging: Path,
    stage_counts: Path,
    stage_lineage: Path,
) -> Path:
    destination_counts = destination_root / "counts"
    lineage = destination_root / "manifests/recount_import"
    if any(destination_counts.iterdir()) or lineage.exists() or lineage.is_symlink():
        raise SystemExit("Destination changed after import preflight")
    destination_counts.rmdir()
    try:
        os.replace(stage_counts, destination_counts)
    except BaseException:
        destination_counts.mkdir(exist_ok=True)
        raise
    os.replace(stage_lineage, lineage)
    staging.rmdir()
    fsync_directory(destination_root)
    fsync_directory(destination_root / "manifests")
    return lineage


def load_scheduler_attestation_module():
    path = Path(__file__).resolve().with_name("slurm_attestation.py")
    require_regular(path, "frozen destination SLURM attestation helper")
    spec = importlib.util.spec_from_file_location(
        f"bg001_slurm_attestation_{uuid.uuid4().hex}", path
    )
    if spec is None or spec.loader is None:
        raise SystemExit("Cannot load frozen destination SLURM attestation helper")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def capture_scheduler_attestation(source_root: Path) -> dict[str, object]:
    module = load_scheduler_attestation_module()
    captured = module.capture(source_root)
    if not isinstance(captured, dict):
        raise SystemExit("Frozen SLURM attestation helper returned malformed evidence")
    return captured


def load_verifier_module(destination_root: Path):
    path = frozen_script(destination_root, "verify_recount_artifacts.py")
    spec = importlib.util.spec_from_file_location(f"bg001_verify_{uuid.uuid4().hex}", path)
    if spec is None or spec.loader is None:
        raise SystemExit("Cannot load destination recount verifier")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_recount_seal_without_completion(destination_root: Path) -> tuple[dict[str, object], str, str]:
    digest_manifest = destination_root / "manifests/recount_artifacts.sha256"
    structural = destination_root / "comparisons/structural_validation.json"
    for path in (digest_manifest, structural, destination_root / "RECOUNT_COMPLETE"):
        if path.exists() or path.is_symlink():
            raise SystemExit(f"Refusing existing destination recount seal: {path}")
    verifier = load_verifier_module(destination_root)
    files, details = verifier.validate_semantics(destination_root, allow_pending_import=True)
    digest_text = "".join(
        f"{sha256(path)}\t{path.relative_to(destination_root)}\n" for path in sorted(files)
    )
    atomic_write_text(digest_manifest, digest_text)
    details["recount_artifacts_sha256"] = sha256(digest_manifest)
    atomic_write_json(structural, details)
    return details, sha256(digest_manifest), sha256(structural)


def write_lineage_manifest(lineage: Path) -> str:
    manifest = lineage / "import_lineage.sha256"
    if manifest.exists() or manifest.is_symlink():
        raise SystemExit("Refusing existing import lineage digest")
    paths = sorted(path for path in lineage.rglob("*") if path.is_file() and not path.is_symlink())
    atomic_write_text(
        manifest,
        "".join(f"{sha256(path)}\t{path.relative_to(lineage)}\n" for path in paths),
    )
    return sha256(manifest)


def verify_archived_scheduler_attestation(
    lineage: Path,
    lineage_json: dict[str, object],
) -> None:
    scheduler = lineage / "slurm"
    ledger = scheduler / "source_submission_ledger.tsv"
    raw = scheduler / "sacct.psv"
    validation_path = scheduler / "validation.json"
    for path, label in (
        (ledger, "archived source SLURM submission ledger"),
        (raw, "archived raw sacct evidence"),
        (validation_path, "archived scheduler validation"),
    ):
        require_regular(path, label)
    validation = load_json(validation_path, "archived scheduler validation")
    normalized = validation.get("normalized_records")
    concurrency = validation.get("observed_array_max_concurrency")
    module = load_scheduler_attestation_module()
    try:
        parsed_raw = module.parse_sacct(raw.read_text())
    except SystemExit as exc:
        raise SystemExit("Archived raw sacct evidence is malformed") from exc
    normalized_without_stage = []
    if isinstance(normalized, list):
        normalized_without_stage = [
            {key: value for key, value in record.items() if key != "stage"}
            for record in normalized
            if isinstance(record, dict)
        ]
    parsed_by_job = {record.get("JobID"): record for record in parsed_raw}
    normalized_by_job = {
        record.get("JobID"): record for record in normalized_without_stage
    }
    if (
        validation.get("status") != "PASS"
        or validation.get("schema_version") != EXPECTED_SCHEDULER_ATTESTATION_SCHEMA
        or validation.get("source_run_id") != lineage_json.get("execution_source_run_id")
        or validation.get("ledger_sha256") != sha256(ledger)
        or validation.get("raw_sacct_sha256") != sha256(raw)
        or validation.get("logical_jobs") != 20
        or validation.get("completed_jobs") != 19
        or validation.get("cancelled_unstarted_jobs") != 1
        or not isinstance(normalized, list)
        or len(normalized) != 20
        or len(parsed_raw) != 20
        or len(parsed_by_job) != 20
        or len(parsed_by_job) != len(parsed_raw)
        or parsed_by_job != normalized_by_job
        or validation.get("sacct_fields") != list(module.SACCT_FIELDS)
        or concurrency != {"bam_manifest_hash": 4, "fragment_recount": 4}
    ):
        raise SystemExit("Archived source scheduler attestation is internally inconsistent")
    summary = lineage_json.get("scheduler_attestation")
    if not isinstance(summary, dict):
        raise SystemExit("Import lineage lacks a scheduler-attestation summary")
    expected_summary = {
        "source_submission_ledger_sha256": sha256(ledger),
        "raw_sacct_sha256": sha256(raw),
        "validation_sha256": sha256(validation_path),
        "logical_jobs": 20,
        "completed_jobs": 19,
        "cancelled_unstarted_jobs": 1,
        "observed_array_max_concurrency": concurrency,
        "schema_version": EXPECTED_SCHEDULER_ATTESTATION_SCHEMA,
    }
    if summary != expected_summary:
        raise SystemExit("Import-lineage scheduler summary differs from archived evidence")
    compatibility_expected = {
        "source_scheduler_attestation_schema": EXPECTED_SCHEDULER_ATTESTATION_SCHEMA,
        "source_scheduler_ledger_sha256": sha256(ledger),
        "source_scheduler_raw_sacct_sha256": sha256(raw),
        "source_scheduler_logical_jobs": 20,
        "source_scheduler_completed_jobs": 19,
        "source_scheduler_cancelled_unstarted_jobs": 1,
        "source_scheduler_observed_array_max_concurrency": concurrency,
    }
    if any(lineage_json.get(key) != value for key, value in compatibility_expected.items()):
        raise SystemExit("Import-lineage scheduler fields differ from archived evidence")


def verify_lineage_digest(destination_root: Path, marker: dict[str, object]) -> dict[str, object]:
    lineage = destination_root / "manifests/recount_import"
    if not lineage.is_dir() or lineage.is_symlink():
        raise SystemExit("Missing or symlinked recount import lineage")
    manifest = lineage / "import_lineage.sha256"
    if sha256(manifest) != marker.get("import_lineage_sha256"):
        raise SystemExit("Import-lineage manifest digest differs")
    records = parse_digest_manifest(lineage, manifest)
    symlinks = [path for path in lineage.rglob("*") if path.is_symlink()]
    if symlinks:
        raise SystemExit(f"Import lineage contains a symlink: {symlinks[0]}")
    observed = {
        str(path.relative_to(lineage))
        for path in lineage.rglob("*")
        if path.is_file() and not path.is_symlink() and path != manifest
    }
    if set(records) != observed:
        raise SystemExit("Import-lineage digest coverage differs")
    destination_seal = read_digest_records(
        destination_root / "manifests/recount_artifacts.sha256"
    )
    lineage_json = load_json(lineage / "lineage.json", "import lineage")
    if sha256(lineage / "lineage.json") != marker.get("lineage_json_sha256"):
        raise SystemExit("Import lineage JSON digest differs")
    verify_archived_scheduler_attestation(lineage, lineage_json)
    source_digest_path = lineage / "source_originals/manifests/recount_artifacts.sha256"
    source_digest = read_digest_records(source_digest_path)
    if sha256(source_digest_path) != lineage_json.get("source_recount_artifacts_sha256"):
        raise SystemExit("Archived source recount seal differs from lineage")
    source_structural_path = (
        lineage / "source_originals/comparisons/structural_validation.json"
    )
    source_complete_path = lineage / "source_originals/RECOUNT_COMPLETE"
    source_structural = load_json(source_structural_path, "archived source structural validation")
    source_complete = load_json(source_complete_path, "archived source RECOUNT_COMPLETE")
    if (
        sha256(source_structural_path)
        != lineage_json.get("source_structural_validation_sha256")
        or sha256(source_complete_path) != lineage_json.get("source_recount_complete_sha256")
        or source_structural.get("status") != "PASS"
        or source_structural.get("run_id") != lineage_json.get("execution_source_run_id")
        or source_complete.get("recount_artifacts_sha256") != sha256(source_digest_path)
        or source_complete.get("structural_validation_sha256") != sha256(source_structural_path)
    ):
        raise SystemExit("Archived source structural/complete seal is not internally consistent")
    source_contract_path = lineage / "source_originals/contract/run_contract.json"
    source_manifest_path = lineage / "source_originals/contract/source_manifest.tsv"
    source_contract = load_json(source_contract_path, "archived execution-source run contract")
    if (
        source_contract.get("run_id") != lineage_json.get("execution_source_run_id")
        or source_contract.get("source_manifest_sha256")
        != lineage_json.get("execution_source_manifest_sha256")
        or sha256(source_manifest_path) != source_contract.get("source_manifest_sha256")
    ):
        raise SystemExit("Archived execution-source contract/manifest lineage differs")
    inventory_path = lineage / "import_inventory.tsv"
    with inventory_path.open(newline="") as handle:
        inventory = list(csv.DictReader(handle, delimiter="\t"))
    by_source = {row["source_relative_path"]: row for row in inventory}
    if len(by_source) != len(inventory) or set(by_source) != set(source_digest):
        raise SystemExit("Import inventory does not cover the exact source sealed artifact set")
    for relative, expected in source_digest.items():
        row = by_source[relative]
        if row["source_sha256"] != expected:
            raise SystemExit(f"Import inventory source digest differs: {relative}")
        preserved_relative = row["preserved_destination_relative_path"]
        preserved = (destination_root / preserved_relative).resolve(strict=True)
        if not contained(preserved, destination_root) or preserved.is_symlink():
            raise SystemExit(f"Preserved source artifact differs: {relative}")
        lineage_prefix = "manifests/recount_import/"
        if preserved_relative.startswith(lineage_prefix):
            lineage_relative = preserved_relative[len(lineage_prefix):]
            preserved_digest = records.get(lineage_relative)
        else:
            preserved_digest = destination_seal.get(preserved_relative)
        if preserved_digest != expected:
            raise SystemExit(f"Preserved source artifact digest differs: {relative}")
        effective_relative = row["effective_destination_relative_path"]
        if effective_relative:
            effective = (destination_root / effective_relative).resolve(strict=True)
            if (
                not contained(effective, destination_root)
                or effective.is_symlink()
                or destination_seal.get(effective_relative) != row["effective_sha256"]
            ):
                raise SystemExit(f"Effective imported artifact differs: {relative}")
    return lineage_json


def verify_imported_destination(
    destination_raw: Path,
    *,
    allow_in_progress: bool = False,
) -> dict[str, object]:
    root, contract = validate_candidate_root(destination_raw, "destination")
    require_frozen_importer_invocation(root)
    if (
        ((root / "RECOUNT_IMPORT_IN_PROGRESS").exists() and not allow_in_progress)
        or (root / "RECOUNT_IMPORT_FAILED.json").exists()
    ):
        raise SystemExit("Destination import is incomplete or failed")
    marker = load_json(root / "RECOUNT_IMPORTED", "RECOUNT_IMPORTED")
    complete = load_json(root / "RECOUNT_COMPLETE", "RECOUNT_COMPLETE")
    if (
        marker.get("status") != "PASS"
        or marker.get("destination_run_id") != root.name
        or marker.get("destination_analysis_source_manifest_sha256") != contract["source_manifest_sha256"]
        or complete.get("recount_imported_sha256") != sha256(root / "RECOUNT_IMPORTED")
        or complete.get("recount_artifacts_sha256") != marker.get("destination_recount_artifacts_sha256")
        or complete.get("structural_validation_sha256") != marker.get("destination_structural_validation_sha256")
    ):
        raise SystemExit("Imported recount markers do not form one cryptographic lineage")
    lineage = verify_lineage_digest(root, marker)
    if (
        lineage.get("status") != "PASS"
        or lineage.get("destination_run_id") != root.name
        or lineage.get("execution_source_run_id") != marker.get("execution_source_run_id")
        or lineage.get("execution_source_manifest_sha256")
        != marker.get("execution_source_manifest_sha256")
    ):
        raise SystemExit("Import lineage JSON differs from RECOUNT_IMPORTED")
    verify_frozen_contract(root)
    verify_structured_analysis_baseline(root)
    verify_baseline(root)
    run_checked(
        [sys.executable, str(frozen_script(root, "verify_recount_artifacts.py")), "--run-root", str(root)],
        label="destination imported recount verification",
    )
    return marker


def prepare_import(
    source_raw: Path,
    destination_raw: Path,
) -> tuple[
    Path, Path, dict[str, object], dict[str, object], dict[str, str],
    list[dict[str, str]], dict[str, object], list[dict[str, str]],
    list[dict[str, str]], list[dict[str, str]], dict[str, object],
    dict[str, object], dict[str, str], dict[str, str],
]:
    source_root, source_contract = validate_candidate_root(source_raw, "source")
    destination_root, destination_contract = validate_candidate_root(destination_raw, "destination")
    require_frozen_importer_invocation(destination_root)
    if source_root == destination_root:
        raise SystemExit("Source and destination runs must differ")
    preflight_destination_empty(destination_root)
    source_contract_log = verify_frozen_contract(source_root)
    destination_contract_log = verify_frozen_contract(destination_root)
    source_analysis_log = verify_structured_analysis_baseline(source_root)
    destination_analysis_log = verify_structured_analysis_baseline(destination_root)
    verify_source_recount(source_root)
    scheduler_attestation = capture_scheduler_attestation(source_root)
    structural = load_json(source_root / "comparisons/structural_validation.json", "source structural validation")
    if structural.get("status") != "PASS" or structural.get("run_id") != source_root.name:
        raise SystemExit("Source structural validation is not PASS for the source run")
    if source_contract.get("featurecounts") != EXPECTED_FEATURECOUNTS:
        raise SystemExit("Source featureCounts contract differs from the BG-001 import contract")
    if destination_contract.get("featurecounts") != source_contract.get("featurecounts"):
        raise SystemExit("Destination featureCounts/GTF contract differs from source execution")
    if (
        destination_contract.get("active_affected_cohorts") != source_contract.get("active_affected_cohorts")
        or destination_contract.get("canonical_affected_cohorts") != source_contract.get("canonical_affected_cohorts")
    ):
        raise SystemExit("Destination affected-cohort roles differ from source execution")

    source_manifest = load_source_manifest(source_root, source_contract)
    destination_manifest = load_source_manifest(destination_root, destination_contract)
    guard_mitigations = validate_destination_guard_mitigations(
        destination_root, destination_manifest
    )
    producers = compare_producer_manifests(
        source_root, destination_root, source_manifest, destination_manifest
    )
    regression_checker = compare_regression_checker(
        source_root, destination_root, source_manifest, destination_manifest
    )
    source_differences = source_manifest_differences(source_manifest, destination_manifest)
    bam = compare_bam_freezes(source_root, destination_root)
    baseline_differences = compare_baselines(source_root, destination_root)
    tools = load_tool_pins(source_root)
    if tools["source_manifest_sha256"] != source_contract["source_manifest_sha256"]:
        raise SystemExit("Source count manifest does not name the source execution snapshot")
    sealed = parse_digest_manifest(
        source_root, source_root / "manifests/recount_artifacts.sha256"
    )
    if not any(relative.startswith("counts/") for relative in sealed):
        raise SystemExit("Source recount seal contains no count artifacts")
    logs = {
        "source_contract": source_contract_log,
        "destination_contract": destination_contract_log,
        "source_analysis_baseline": source_analysis_log,
        "destination_analysis_baseline": destination_analysis_log,
    }
    return (
        source_root, destination_root, source_contract, destination_contract, sealed,
        producers, regression_checker, guard_mitigations, source_differences,
        baseline_differences, bam, scheduler_attestation, tools, logs,
    )


def perform_import(source_raw: Path, destination_raw: Path, *, check_only: bool) -> dict[str, object]:
    (
        source_root, destination_root, source_contract, destination_contract, sealed,
        producers, regression_checker, guard_mitigations, source_differences,
        baseline_differences, bam, scheduler_attestation, tools, preflight_logs,
    ) = prepare_import(source_raw, destination_raw)
    scheduler_validation = scheduler_attestation.get("validation")
    if not isinstance(scheduler_validation, dict):
        raise SystemExit("Source scheduler attestation validation is malformed")
    compatibility = {
        "status": "PASS",
        "execution_source_run_id": source_root.name,
        "destination_run_id": destination_root.name,
        "execution_source_manifest_sha256": source_contract["source_manifest_sha256"],
        "destination_analysis_source_manifest_sha256": destination_contract["source_manifest_sha256"],
        "bam_manifest_sha256": bam["bam_manifest_sha256"],
        "producer_files": len(producers),
        "destination_guard_mitigation_files": len(guard_mitigations),
        "destination_guard_mitigation_sha256": {
            str(row["relative_path"]): str(row["sha256"])
            for row in guard_mitigations
        },
        "regression_checker_transition": (
            f"{regression_checker['execution_generation']}->"
            f"{regression_checker['destination_generation']}"
        ),
        "bam_guard_transition": (
            "legacy_source_no_guard->guarded_destination"
            if baseline_generation(source_root) == "legacy_zero_marker"
            else "guarded_source->guarded_destination"
        ),
        "source_manifest_differences": len(source_differences),
        "protected_baseline_differences": len(baseline_differences),
        "source_scheduler_attestation_schema": scheduler_validation.get("schema_version"),
        "source_scheduler_ledger_sha256": scheduler_validation.get("ledger_sha256"),
        "source_scheduler_raw_sacct_sha256": scheduler_validation.get("raw_sacct_sha256"),
        "source_scheduler_logical_jobs": scheduler_validation.get("logical_jobs"),
        "source_scheduler_completed_jobs": scheduler_validation.get("completed_jobs"),
        "source_scheduler_cancelled_unstarted_jobs": scheduler_validation.get(
            "cancelled_unstarted_jobs"
        ),
        "source_scheduler_observed_array_max_concurrency": scheduler_validation.get(
            "observed_array_max_concurrency"
        ),
        "checked_utc": utc_now(),
    }
    if check_only:
        return compatibility

    in_progress = destination_root / "RECOUNT_IMPORT_IN_PROGRESS"
    staging = destination_root / f".recount_import_staging.{uuid.uuid4().hex}"
    if any(destination_root.glob(".recount_import_staging.*")):
        raise SystemExit("Destination contains a prior recount-import staging directory")
    atomic_write_json(
        in_progress,
        {
            "status": "IN_PROGRESS",
            "execution_source_run_id": source_root.name,
            "destination_run_id": destination_root.name,
            "staging_directory": staging.name,
            "started_utc": utc_now(),
        },
    )
    staging.mkdir(mode=0o700)
    source_before = source_state(source_root, sealed)
    stage_counts, stage_lineage, inventory = copy_and_rebase_counts(
        source_root, destination_root, sealed, staging
    )
    archive_extra_source_lineage(source_root, stage_lineage)
    scheduler_lineage = archive_scheduler_attestation(
        source_root, stage_lineage, scheduler_attestation
    )
    archive_regression_checker_transition(
        source_root, destination_root, stage_lineage, regression_checker
    )
    write_tsv(
        stage_lineage / "producer_contract.tsv",
        ("relative_path", "sha256"),
        producers,
    )
    write_tsv(
        stage_lineage / "destination_guard_mitigation_bundle.tsv",
        ("relative_path", "size_bytes", "sha256"),
        guard_mitigations,
    )
    write_tsv(
        stage_lineage / "source_manifest_differences.tsv",
        (
            "relative_path", "status", "execution_source_sha256", "destination_analysis_sha256",
            "execution_source_size", "destination_analysis_size",
        ),
        source_differences,
    )
    write_tsv(
        stage_lineage / "baseline_differences.tsv",
        ("relative_path", "execution_source_sha256", "destination_analysis_sha256", "change_class"),
        baseline_differences,
    )
    write_tsv(
        stage_lineage / "import_inventory.tsv",
        (
            "source_relative_path", "source_sha256", "disposition",
            "preserved_destination_relative_path", "effective_destination_relative_path",
            "effective_sha256",
        ),
        inventory,
    )
    for label, output in preflight_logs.items():
        atomic_write_text(stage_lineage / f"{label}_verification.txt", output)
    atomic_write_json(stage_lineage / "bam_compatibility.json", bam)
    atomic_write_json(stage_lineage / "tool_provenance.json", tools)

    verify_frozen_contract(source_root)
    verify_source_recount(source_root)
    if source_state(source_root, sealed) != source_before:
        raise SystemExit("Source recount artifacts changed during import copy")
    fsync_tree_directories(stage_counts)
    fsync_tree_directories(stage_lineage)
    lineage = publish_staging(destination_root, staging, stage_counts, stage_lineage)

    # Reverify the destination snapshot immediately before executing its
    # post-count helpers, then again after they return.
    verify_frozen_contract(destination_root)
    finalize_output = run_checked(
        [
            sys.executable,
            str(frozen_script(destination_root, "finalize_recount.py")),
            "--run-root",
            str(destination_root),
            "--artifact-origin",
            "imported",
            "--execution-source-run-id",
            source_root.name,
            "--execution-source-manifest-sha256",
            str(source_contract["source_manifest_sha256"]),
        ],
        label="destination imported recount finalization",
    )
    atomic_write_text(lineage / "destination_finalize_recount.txt", finalize_output)
    override_output = run_checked(
        [
            sys.executable,
            str(frozen_script(destination_root, "make_count_overrides.py")),
            "--run-root",
            str(destination_root),
            "--count-manifest",
            str(destination_root / "manifests/count_manifest.tsv"),
            "--baseline-overrides",
            str(destination_root / "manifests/read_count_overrides.tsv"),
            "--output",
            str(destination_root / "manifests/count_overrides.tsv"),
        ],
        label="destination imported count-override construction",
    )
    atomic_write_text(lineage / "destination_make_count_overrides.txt", override_output)
    verify_frozen_contract(destination_root)
    verify_frozen_contract(source_root)
    verify_source_recount(source_root)
    if source_state(source_root, sealed) != source_before:
        raise SystemExit("Source recount artifacts changed during destination finalization")

    details, recount_digest, structural_digest = write_recount_seal_without_completion(destination_root)
    lineage_payload: dict[str, object] = {
        **compatibility,
        "status": "PASS",
        "schema_version": "1.0",
        "execution_source_root": str(source_root),
        "destination_root": str(destination_root),
        "source_recount_artifacts_sha256": sha256(
            source_root / "manifests/recount_artifacts.sha256"
        ),
        "source_structural_validation_sha256": sha256(
            source_root / "comparisons/structural_validation.json"
        ),
        "source_recount_complete_sha256": sha256(source_root / "RECOUNT_COMPLETE"),
        "destination_recount_artifacts_sha256": recount_digest,
        "destination_structural_validation_sha256": structural_digest,
        "destination_count_manifest_sha256": sha256(
            destination_root / "manifests/count_manifest.tsv"
        ),
        "destination_count_overrides_sha256": sha256(
            destination_root / "manifests/count_overrides.tsv"
        ),
        "source_sealed_artifacts": len(sealed),
        "source_count_artifacts_copied": sum(
            relative.startswith("counts/") for relative in sealed
        ),
        "path_bound_artifacts_rebased": sorted(PATH_BOUND_COUNT_FILES),
        "protected_baseline_policy": (
            "all non-protected baseline hashes exact; protected source set and gene reasons/sources monotonic"
        ),
        "tool_provenance": tools,
        "scheduler_attestation": scheduler_lineage,
        "structural_details": details,
        "completed_utc": utc_now(),
    }
    atomic_write_json(lineage / "lineage.json", lineage_payload)
    lineage_manifest_digest = write_lineage_manifest(lineage)
    imported_payload = {
        "status": "PASS",
        "schema_version": "1.0",
        "execution_source_run_id": source_root.name,
        "execution_source_manifest_sha256": source_contract["source_manifest_sha256"],
        "destination_run_id": destination_root.name,
        "destination_analysis_source_manifest_sha256": destination_contract["source_manifest_sha256"],
        "destination_recount_artifacts_sha256": recount_digest,
        "destination_structural_validation_sha256": structural_digest,
        "import_lineage_sha256": lineage_manifest_digest,
        "lineage_json_sha256": sha256(lineage / "lineage.json"),
    }
    atomic_write_json(destination_root / "RECOUNT_IMPORTED", imported_payload)
    atomic_write_json(
        destination_root / "RECOUNT_COMPLETE",
        {
            "recount_artifacts_sha256": recount_digest,
            "structural_validation_sha256": structural_digest,
            "recount_imported_sha256": sha256(destination_root / "RECOUNT_IMPORTED"),
        },
    )
    try:
        verify_imported_destination(destination_root, allow_in_progress=True)
    except BaseException:
        atomic_write_json(
            destination_root / "RECOUNT_IMPORT_FAILED.json",
            {
                "status": "FAILED_POST_SEAL_VERIFICATION",
                "execution_source_run_id": source_root.name,
                "destination_run_id": destination_root.name,
                "failed_utc": utc_now(),
            },
        )
        raise
    in_progress.unlink()
    fsync_directory(destination_root)
    marker = verify_imported_destination(destination_root)
    return marker


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--source-run", type=Path)
    mode.add_argument("--verify-import", action="store_true")
    parser.add_argument("--destination-run", required=True, type=Path)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.verify_import:
        if args.check_only:
            raise SystemExit("--check-only is valid only with --source-run")
        result = verify_imported_destination(args.destination_run)
    else:
        assert args.source_run is not None
        result = perform_import(args.source_run, args.destination_run, check_only=args.check_only)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt as exc:
        raise SystemExit("Interrupted; destination remains fail-closed") from exc
