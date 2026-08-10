#!/usr/bin/env python3
"""Pre/postflight the legacy BG-001 BAM-manifest publication transaction.

The frozen producer intentionally remains byte-identical to the producer used by
the source recount.  That producer uses three predictable temporary names.  This
guard makes that retained behavior safe against other users by requiring an
owner-only publication directory, proves that none of those names existed just
before execution, and independently reconstructs and verifies every published
byte immediately afterwards.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import stat
import subprocess
import sys
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from typing import Any

# A frozen source snapshot must not gain __pycache__ entries when this guard
# imports its colocated helper or invokes the frozen contract verifier.
sys.dont_write_bytecode = True

from safe_io import SafeIOError, publish_new_bytes  # noqa: E402


EXPECTED = {
    "GSE130970": 78,
    "GSE135251": 216,
    "GSE174478": 93,
    "GSE213621": 367,
    "GSE240729": 66,
}
DOCUMENTED_UNAVAILABLE = {
    ("GSE174478", "SRR14551000"): "documented_failed_run",
    ("GSE240729", "SRR25630203"): "documented_failed_run",
}
MANIFEST_FIELDS = (
    "dataset",
    "sample_id",
    "bam_path",
    "layout",
    "strandedness",
    "expected_status",
    "exclusion_reason",
    "bam_size",
    "bam_mtime",
    "bam_sha256",
)
SHARD_FIELDS = (
    "task",
    "dataset",
    "sample_id",
    "bam_path",
    "bam_size",
    "bam_mtime",
    "bam_sha256",
)
TASK_DATASETS = ("GSE130970", "GSE135251", "GSE174478", "GSE240729")
TASK_CARDINALITIES = (78, 216, 93, 66, 92, 92, 92, 91)
SHA256_RE = re.compile(r"[0-9a-f]{64}")

PRE_ARTIFACT = Path("manifests/bam_manifest_finalize.preflight.json")
POST_ARTIFACT = Path("manifests/bam_manifest_finalize.postflight.json")
FINAL_OUTPUTS = (
    Path("manifests/bam_manifest.tsv"),
    Path("manifests/bam_manifest.freeze.json"),
    Path("BAM_MANIFEST_FROZEN"),
)
# These names exactly match Path.with_suffix() in finalize_bam_manifest.py.
PREDICTABLE_TEMPORARIES = (
    Path("manifests/bam_manifest.tsv.tmp"),
    Path("manifests/bam_manifest.freeze.json.tmp"),
    Path("BAM_MANIFEST_FROZEN.tmp"),
)


class GuardError(RuntimeError):
    """A fail-closed BAM-manifest transaction check failed."""


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def absent_even_if_dangling(path: Path) -> bool:
    try:
        os.lstat(path)
    except FileNotFoundError:
        return True
    return False


def _identity(path: Path, status: os.stat_result) -> dict[str, int | str]:
    return {
        "path": str(path),
        "device": status.st_dev,
        "inode": status.st_ino,
        "uid": status.st_uid,
        "mode": stat.S_IMODE(status.st_mode),
        "nlink": status.st_nlink,
    }


def require_owner_only_directory(path: Path, label: str) -> dict[str, int | str]:
    try:
        status = os.lstat(path)
    except FileNotFoundError as exc:
        raise GuardError(f"Missing {label}: {path}") from exc
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
        raise GuardError(f"{label} must be a non-symlink directory: {path}")
    if status.st_uid != os.geteuid():
        raise GuardError(f"{label} is not owned by effective UID {os.geteuid()}: {path}")
    permissions = stat.S_IMODE(status.st_mode)
    if permissions & 0o700 != 0o700 or permissions & 0o077:
        raise GuardError(
            f"{label} must be owner-only rwx (group/other mode bits zero): "
            f"{path} mode={permissions:#05o}"
        )
    return _identity(path, status)


def require_private_parent(path: Path, label: str) -> None:
    try:
        status = os.lstat(path)
    except FileNotFoundError as exc:
        raise GuardError(f"Missing {label}: {path}") from exc
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
        raise GuardError(f"{label} must be a non-symlink directory: {path}")
    if status.st_uid != os.geteuid() or stat.S_IMODE(status.st_mode) & 0o022:
        raise GuardError(
            f"{label} must be owned by effective UID and not group/other writable: "
            f"{path} mode={stat.S_IMODE(status.st_mode):#05o}"
        )


def require_safe_regular(
    path: Path,
    label: str,
    *,
    owner_only: bool,
    single_link: bool,
) -> dict[str, int | str]:
    try:
        status = os.lstat(path)
    except FileNotFoundError as exc:
        raise GuardError(f"Missing {label}: {path}") from exc
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISREG(status.st_mode):
        raise GuardError(f"{label} must be a non-symlink regular file: {path}")
    if status.st_uid != os.geteuid():
        raise GuardError(f"{label} is not owned by effective UID {os.geteuid()}: {path}")
    permissions = stat.S_IMODE(status.st_mode)
    if permissions & 0o111:
        raise GuardError(f"{label} must not be executable: {path} mode={permissions:#05o}")
    if owner_only and permissions & 0o077:
        raise GuardError(f"{label} must be owner-only: {path} mode={permissions:#05o}")
    if not permissions & 0o400:
        raise GuardError(f"{label} is not owner-readable: {path} mode={permissions:#05o}")
    if single_link and status.st_nlink != 1:
        raise GuardError(f"{label} must have exactly one hard link: {path}")
    return _identity(path, status)


def canonical_run_root(raw_root: Path) -> Path:
    root = Path(os.path.abspath(raw_root))
    require_private_parent(root.parent, "candidate parent directory")
    require_owner_only_directory(root, "run root")
    resolved = root.resolve(strict=True)
    if resolved != root:
        raise GuardError(f"Run-root path traverses a symlink or is non-canonical: {root}")
    manifests = root / "manifests"
    require_owner_only_directory(manifests, "manifests directory")
    require_owner_only_directory(manifests / "bam_hash_shards", "BAM hash-shard directory")
    return root


def read_regular_bytes(
    path: Path,
    label: str,
    *,
    owner_only: bool = True,
    single_link: bool = True,
) -> tuple[bytes, dict[str, int | str]]:
    before = require_safe_regular(
        path, label, owner_only=owner_only, single_link=single_link
    )
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    try:
        opened = os.fstat(descriptor)
        if opened.st_dev != before["device"] or opened.st_ino != before["inode"]:
            raise GuardError(f"{label} changed while it was opened: {path}")
        chunks: list[bytes] = []
        while True:
            block = os.read(descriptor, 8 * 1024 * 1024)
            if not block:
                break
            chunks.append(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if (
        after.st_dev != opened.st_dev
        or after.st_ino != opened.st_ino
        or after.st_size != opened.st_size
        or after.st_mtime_ns != opened.st_mtime_ns
        or after.st_ctime_ns != opened.st_ctime_ns
    ):
        raise GuardError(f"{label} changed while it was read: {path}")
    return b"".join(chunks), before


def parse_tsv_bytes(
    payload: bytes, expected_fields: tuple[str, ...], label: str
) -> list[dict[str, str]]:
    if b"\r" in payload:
        raise GuardError(f"{label} must use LF-only line endings")
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise GuardError(f"{label} is not UTF-8") from exc
    reader = csv.DictReader(StringIO(text, newline=""), delimiter="\t")
    if tuple(reader.fieldnames or ()) != expected_fields:
        raise GuardError(f"{label} schema/order differs from the frozen contract")
    try:
        rows = list(reader)
    except csv.Error as exc:
        raise GuardError(f"Malformed TSV in {label}") from exc
    if any(None in row for row in rows):
        raise GuardError(f"{label} contains a row with extra fields")
    return rows


def validate_draft(rows: list[dict[str, str]]) -> None:
    if len(rows) != sum(EXPECTED.values()) + len(DOCUMENTED_UNAVAILABLE):
        raise GuardError("Draft BAM manifest does not contain exactly 820+2 rows")
    included = [row for row in rows if row["expected_status"] == "included"]
    unavailable = [
        row for row in rows if row["expected_status"] == "documented_unavailable"
    ]
    if len(included) != 820 or len(unavailable) != 2:
        raise GuardError("Draft BAM manifest status counts differ from 820 included + 2 unavailable")
    if any(
        row["expected_status"] not in {"included", "documented_unavailable"}
        for row in rows
    ):
        raise GuardError("Draft BAM manifest contains an unexpected status")
    keys = [(row["dataset"], row["sample_id"]) for row in included]
    paths = [row["bam_path"] for row in included]
    if len(keys) != len(set(keys)) or len(paths) != len(set(paths)):
        raise GuardError("Draft BAM manifest contains duplicate included samples or paths")
    for dataset, expected in EXPECTED.items():
        observed = sum(row["dataset"] == dataset for row in included)
        if observed != expected:
            raise GuardError(f"{dataset} draft count differs: expected {expected}, found {observed}")
    for row in included:
        if (
            row["layout"] != "paired"
            or row["strandedness"] != "2"
            or row["exclusion_reason"]
            or row["bam_sha256"]
        ):
            raise GuardError(f"Invalid included draft semantics: {row['dataset']}/{row['sample_id']}")
        if (
            not row["sample_id"]
            or not row["bam_path"]
            or not Path(row["bam_path"]).is_absolute()
            or not row["bam_size"].isdigit()
            or not row["bam_mtime"].isdigit()
        ):
            raise GuardError(f"Incomplete included draft provenance: {row['dataset']}/{row['sample_id']}")
    unavailable_by_key = {
        (row["dataset"], row["sample_id"]): row for row in unavailable
    }
    if set(unavailable_by_key) != set(DOCUMENTED_UNAVAILABLE):
        raise GuardError("Documented-unavailable sample identities differ from the contract")
    for key, reason in DOCUMENTED_UNAVAILABLE.items():
        row = unavailable_by_key[key]
        if (
            row["layout"] != "paired"
            or row["strandedness"] != "2"
            or row["exclusion_reason"] != reason
            or any(row[field] for field in ("bam_path", "bam_size", "bam_mtime", "bam_sha256"))
        ):
            raise GuardError(f"Invalid documented-unavailable row: {key}")


def expected_task_rows(rows: list[dict[str, str]], task: int) -> list[dict[str, str]]:
    included = [row for row in rows if row["expected_status"] == "included"]
    if task < 4:
        return [row for row in included if row["dataset"] == TASK_DATASETS[task]]
    gse213621 = [row for row in included if row["dataset"] == "GSE213621"]
    start = (task - 4) * 92
    stop = min(start + 92, len(gse213621))
    return gse213621[start:stop]


def input_paths(root: Path) -> list[Path]:
    return [root / "manifests/bam_manifest.draft.tsv"] + [
        root / f"manifests/bam_hash_shards/task_{task}.tsv" for task in range(8)
    ]


def input_records(root: Path) -> dict[str, dict[str, int | str]]:
    expected_shards = {f"task_{task}.tsv" for task in range(8)}
    shard_dir = root / "manifests/bam_hash_shards"
    observed = {path.name for path in shard_dir.iterdir()}
    if observed != expected_shards:
        raise GuardError(
            "BAM hash-shard directory differs from exactly task_0.tsv..task_7.tsv: "
            f"missing={sorted(expected_shards-observed)}, extra={sorted(observed-expected_shards)}"
        )
    records: dict[str, dict[str, int | str]] = {}
    for path in input_paths(root):
        payload, identity = read_regular_bytes(path, "BAM-manifest input")
        relative = str(path.relative_to(root))
        records[relative] = {
            **identity,
            "size": len(payload),
            "sha256": sha256_bytes(payload),
        }
    return records


def require_absent(paths: tuple[Path, ...], root: Path, label: str) -> None:
    present = [str(path) for relative in paths if not absent_even_if_dangling(path := root / relative)]
    if present:
        raise GuardError(f"{label} must be absent: {present}")


def run_contract_verifier(root: Path, require_frozen: bool) -> dict[str, str]:
    verifier = root / (
        "source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/"
        "bg001_remediation/verify_run_contract.py"
    )
    verifier_bytes, _ = read_regular_bytes(verifier, "frozen run-contract verifier")
    command = [sys.executable, str(verifier), "--run-root", str(root)]
    if require_frozen:
        command.append("--require-bam-frozen")
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise GuardError(f"Frozen verify_run_contract.py failed: {detail}")
    return {
        "command": " ".join(command),
        "script_sha256": sha256_bytes(verifier_bytes),
        "stdout": completed.stdout.strip(),
    }


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def publish_json(path: Path, record: dict[str, Any], root: Path) -> None:
    payload = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
    try:
        publish_new_bytes(path, payload, root, mode=0o600)
    except SafeIOError as exc:
        raise GuardError(str(exc)) from exc
    require_safe_regular(path, "guard artifact", owner_only=True, single_link=True)


def preflight(root: Path, *, invoke_verifier: bool = True) -> dict[str, Any]:
    root = canonical_run_root(root)
    require_absent(FINAL_OUTPUTS, root, "BAM-manifest final outputs")
    require_absent(PREDICTABLE_TEMPORARIES, root, "predictable BAM-manifest temporaries")
    require_absent((PRE_ARTIFACT, POST_ARTIFACT), root, "guard transaction artifacts")
    inputs = input_records(root)
    draft_payload, _ = read_regular_bytes(
        root / "manifests/bam_manifest.draft.tsv", "draft BAM manifest"
    )
    validate_draft(parse_tsv_bytes(draft_payload, MANIFEST_FIELDS, "draft BAM manifest"))
    verifier_record = run_contract_verifier(root, require_frozen=False) if invoke_verifier else None
    record: dict[str, Any] = {
        "schema_version": "1.0",
        "phase": "preflight",
        "status": "PASS",
        "run_id": root.name,
        "run_root": require_owner_only_directory(root, "run root"),
        "manifests_directory": require_owner_only_directory(root / "manifests", "manifests directory"),
        "inputs": inputs,
        "final_outputs_absent": [str(path) for path in FINAL_OUTPUTS],
        "predictable_temporaries_absent": [str(path) for path in PREDICTABLE_TEMPORARIES],
        "security_contract": (
            "The frozen producer's predictable temporary names were absent inside "
            "owner-only directories immediately before producer execution."
        ),
        "contract_verifier": verifier_record,
        "created_utc": utc_now(),
    }
    publish_json(root / PRE_ARTIFACT, record, root)
    return record


def load_preflight(root: Path) -> tuple[dict[str, Any], str]:
    payload, _ = read_regular_bytes(root / PRE_ARTIFACT, "BAM-manifest preflight artifact")
    try:
        record = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise GuardError("BAM-manifest preflight artifact is invalid JSON") from exc
    if (
        record.get("schema_version") != "1.0"
        or record.get("phase") != "preflight"
        or record.get("status") != "PASS"
        or record.get("run_id") != root.name
    ):
        raise GuardError("BAM-manifest preflight artifact identity differs")
    return record, sha256_bytes(payload)


def _without_device(value: Any) -> Any:
    """Recursively drop st_dev from recorded stat identities.

    st_dev is NOT stable across a remount of a clustered filesystem. BG-001 stages
    are separated by days (the recount alone runs ~40 h), so comparing st_dev
    across stages raises a false tamper alarm whenever GPFS is remounted between
    the guarded finalize and a later stage — observed 2026-08-09, when every path
    in the project moved from device 61 to 63 with all inodes unchanged.

    Every other recorded field — path, inode, uid, mode, nlink — still pins the
    object, and st_dev cannot be altered without remounting the filesystem itself,
    so excluding it here does not weaken the guard. The intra-operation TOCTOU
    checks that compare st_dev between open() and fstat() are deliberately left
    untouched: a remount cannot occur inside a single open/fstat pair.
    """
    if isinstance(value, dict):
        return {k: _without_device(v) for k, v in value.items() if k != "device"}
    if isinstance(value, list):
        return [_without_device(v) for v in value]
    return value


def validate_preflight_continuity(root: Path, preflight_record: dict[str, Any]) -> None:
    current_root = require_owner_only_directory(root, "run root")
    current_manifests = require_owner_only_directory(root / "manifests", "manifests directory")
    if _without_device(preflight_record.get("run_root")) != _without_device(current_root):
        raise GuardError("Run-root identity or mode changed after preflight")
    if _without_device(preflight_record.get("manifests_directory")) != _without_device(current_manifests):
        raise GuardError("Manifests-directory identity or mode changed after preflight")
    current_inputs = input_records(root)
    if _without_device(preflight_record.get("inputs")) != _without_device(current_inputs):
        raise GuardError("Draft BAM manifest or hash shards changed after preflight")


def reconstruct_manifest(
    root: Path,
) -> tuple[bytes, dict[str, str], list[int]]:
    draft_payload, _ = read_regular_bytes(
        root / "manifests/bam_manifest.draft.tsv", "draft BAM manifest"
    )
    rows = parse_tsv_bytes(draft_payload, MANIFEST_FIELDS, "draft BAM manifest")
    validate_draft(rows)
    hash_by_key: dict[tuple[str, str], str] = {}
    shard_hashes: dict[str, str] = {}
    cardinalities: list[int] = []
    for task in range(8):
        path = root / f"manifests/bam_hash_shards/task_{task}.tsv"
        payload, _ = read_regular_bytes(path, f"BAM hash shard task {task}")
        observed = parse_tsv_bytes(payload, SHARD_FIELDS, f"BAM hash shard task {task}")
        expected_rows = expected_task_rows(rows, task)
        cardinalities.append(len(observed))
        if len(observed) != TASK_CARDINALITIES[task] or len(observed) != len(expected_rows):
            raise GuardError(
                f"Task {task} cardinality differs: expected {TASK_CARDINALITIES[task]}, "
                f"found {len(observed)}"
            )
        for expected_row, observed_row in zip(expected_rows, observed):
            key = (expected_row["dataset"], expected_row["sample_id"])
            if key in hash_by_key:
                raise GuardError(f"Duplicate BAM hash attribution: {key}")
            for field in ("dataset", "sample_id", "bam_path", "bam_size", "bam_mtime"):
                if observed_row[field] != expected_row[field]:
                    raise GuardError(f"Task {task} provenance differs for {key}: {field}")
            digest = observed_row["bam_sha256"]
            if observed_row["task"] != str(task) or not SHA256_RE.fullmatch(digest):
                raise GuardError(f"Task {task} has invalid task or SHA-256 for {key}")
            bam = Path(expected_row["bam_path"])
            try:
                bam_status = os.lstat(bam)
            except FileNotFoundError as exc:
                raise GuardError(f"BAM disappeared after hashing: {bam}") from exc
            if stat.S_ISLNK(bam_status.st_mode) or not stat.S_ISREG(bam_status.st_mode):
                raise GuardError(f"BAM is symlinked or non-regular after hashing: {bam}")
            if (
                str(bam_status.st_size) != expected_row["bam_size"]
                or str(int(bam_status.st_mtime)) != expected_row["bam_mtime"]
            ):
                raise GuardError(f"BAM size/mtime changed after hashing: {bam}")
            hash_by_key[key] = digest
        shard_hashes[path.name] = sha256_bytes(payload)
    if cardinalities[4:] != [92, 92, 92, 91] or len(hash_by_key) != 820:
        raise GuardError("Exact 92/92/92/91 GSE213621 chunk contract was not reconstructed")
    final_rows = [dict(row) for row in rows]
    for row in final_rows:
        if row["expected_status"] == "included":
            row["bam_sha256"] = hash_by_key[(row["dataset"], row["sample_id"])]
    output = StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=MANIFEST_FIELDS,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(final_rows)
    return output.getvalue().encode("utf-8"), shard_hashes, cardinalities


def validate_freeze_transaction(
    root: Path,
    reconstructed: bytes,
    shard_hashes: dict[str, str],
) -> dict[str, Any]:
    final_payload, final_identity = read_regular_bytes(
        root / FINAL_OUTPUTS[0], "final BAM manifest"
    )
    freeze_payload, freeze_identity = read_regular_bytes(
        root / FINAL_OUTPUTS[1], "BAM-manifest freeze JSON"
    )
    marker_payload, marker_identity = read_regular_bytes(
        root / FINAL_OUTPUTS[2], "BAM-manifest frozen marker"
    )
    if final_payload != reconstructed:
        raise GuardError("Final BAM manifest is not byte-identical to independent reconstruction")
    final_sha = sha256_bytes(final_payload)
    expected_marker = (final_sha + "\n").encode()
    if marker_payload != expected_marker:
        raise GuardError("BAM_MANIFEST_FROZEN is not the exact final-manifest SHA-256 marker")
    try:
        freeze_record = json.loads(freeze_payload)
    except json.JSONDecodeError as exc:
        raise GuardError("BAM-manifest freeze record is invalid JSON") from exc
    expected_keys = {
        "schema_version",
        "run_id",
        "draft_sha256",
        "manifest_sha256",
        "included_bams",
        "documented_unavailable",
        "shard_sha256",
        "finalized_utc",
    }
    if set(freeze_record) != expected_keys:
        raise GuardError("BAM-manifest freeze record fields differ from the frozen producer contract")
    finalized = freeze_record.get("finalized_utc")
    try:
        finalized_time = datetime.fromisoformat(finalized)
    except (TypeError, ValueError) as exc:
        raise GuardError("BAM-manifest freeze record has invalid finalized_utc") from exc
    if finalized_time.tzinfo is None or finalized_time.utcoffset() != timezone.utc.utcoffset(None):
        raise GuardError("BAM-manifest freeze record finalized_utc is not UTC")
    draft_sha = sha256(root / "manifests/bam_manifest.draft.tsv")
    if (
        freeze_record["schema_version"] != "1.0"
        or freeze_record["run_id"] != root.name
        or freeze_record["draft_sha256"] != draft_sha
        or freeze_record["manifest_sha256"] != final_sha
        or freeze_record["included_bams"] != 820
        or freeze_record["documented_unavailable"] != 2
        or freeze_record["shard_sha256"] != shard_hashes
    ):
        raise GuardError("BAM-manifest freeze record hashes or counts differ")
    return {
        "manifest_sha256": final_sha,
        "freeze_json_sha256": sha256_bytes(freeze_payload),
        "marker_sha256": sha256_bytes(marker_payload),
        "outputs": {
            str(FINAL_OUTPUTS[0]): final_identity,
            str(FINAL_OUTPUTS[1]): freeze_identity,
            str(FINAL_OUTPUTS[2]): marker_identity,
        },
    }


def postflight(root: Path, *, invoke_verifier: bool = True) -> dict[str, Any]:
    root = canonical_run_root(root)
    require_absent(PREDICTABLE_TEMPORARIES, root, "predictable BAM-manifest temporaries")
    require_absent((POST_ARTIFACT,), root, "BAM-manifest postflight artifact")
    preflight_record, preflight_sha = load_preflight(root)
    validate_preflight_continuity(root, preflight_record)
    reconstructed, shard_hashes, cardinalities = reconstruct_manifest(root)
    transaction = validate_freeze_transaction(root, reconstructed, shard_hashes)
    verifier_record = run_contract_verifier(root, require_frozen=True) if invoke_verifier else None
    record: dict[str, Any] = {
        "schema_version": "1.0",
        "phase": "postflight",
        "status": "PASS",
        "run_id": root.name,
        "preflight_sha256": preflight_sha,
        "input_sha256": {
            relative: details["sha256"]
            for relative, details in sorted(preflight_record["inputs"].items())
        },
        "included_bams": 820,
        "documented_unavailable": 2,
        "task_cardinalities": cardinalities,
        "gse213621_chunk_cardinalities": cardinalities[4:],
        "transaction": transaction,
        "predictable_temporaries_absent": [str(path) for path in PREDICTABLE_TEMPORARIES],
        "contract_verifier": verifier_record,
        "created_utc": utc_now(),
    }
    publish_json(root / POST_ARTIFACT, record, root)
    return record


def _require_utc_timestamp(value: object, label: str) -> None:
    try:
        parsed = datetime.fromisoformat(value)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise GuardError(f"{label} is not an ISO-8601 timestamp") from exc
    if parsed.tzinfo is None or parsed.utcoffset() != timezone.utc.utcoffset(None):
        raise GuardError(f"{label} is not UTC")


def verify_artifacts(root: Path, *, invoke_verifier: bool = True) -> dict[str, Any]:
    """Reconstruct and reverify an already-published guarded transaction."""

    root = canonical_run_root(root)
    require_absent(PREDICTABLE_TEMPORARIES, root, "predictable BAM-manifest temporaries")
    preflight_record, preflight_sha = load_preflight(root)
    validate_preflight_continuity(root, preflight_record)
    reconstructed, shard_hashes, cardinalities = reconstruct_manifest(root)
    transaction = validate_freeze_transaction(root, reconstructed, shard_hashes)

    postflight_payload, _ = read_regular_bytes(
        root / POST_ARTIFACT, "BAM-manifest postflight artifact"
    )
    try:
        postflight_record = json.loads(postflight_payload)
    except json.JSONDecodeError as exc:
        raise GuardError("BAM-manifest postflight artifact is invalid JSON") from exc
    expected_keys = {
        "schema_version", "phase", "status", "run_id", "preflight_sha256",
        "input_sha256", "included_bams", "documented_unavailable",
        "task_cardinalities", "gse213621_chunk_cardinalities", "transaction",
        "predictable_temporaries_absent", "contract_verifier", "created_utc",
    }
    if set(postflight_record) != expected_keys:
        raise GuardError("BAM-manifest postflight fields differ from the guarded contract")
    input_sha256 = {
        relative: details["sha256"]
        for relative, details in sorted(preflight_record["inputs"].items())
    }
    if (
        postflight_record["schema_version"] != "1.0"
        or postflight_record["phase"] != "postflight"
        or postflight_record["status"] != "PASS"
        or postflight_record["run_id"] != root.name
        or postflight_record["preflight_sha256"] != preflight_sha
        or postflight_record["input_sha256"] != input_sha256
        or postflight_record["included_bams"] != 820
        or postflight_record["documented_unavailable"] != 2
        or postflight_record["task_cardinalities"] != cardinalities
        or postflight_record["gse213621_chunk_cardinalities"] != cardinalities[4:]
        # transaction.outputs carries per-file stat identities, so it needs the
        # same remount-tolerant comparison as validate_preflight_continuity.
        or _without_device(postflight_record["transaction"]) != _without_device(transaction)
        or postflight_record["predictable_temporaries_absent"]
        != [str(path) for path in PREDICTABLE_TEMPORARIES]
    ):
        raise GuardError("BAM-manifest postflight record differs from independent reconstruction")
    _require_utc_timestamp(preflight_record.get("created_utc"), "preflight created_utc")
    _require_utc_timestamp(postflight_record.get("created_utc"), "postflight created_utc")

    if invoke_verifier:
        observed_post = run_contract_verifier(root, require_frozen=True)
        verifier_path = str(root / (
            "source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/"
            "bg001_remediation/verify_run_contract.py"
        ))

        def validate_verifier_record(
            record: object,
            *,
            require_frozen: bool,
            expected_stdout: str,
        ) -> None:
            if not isinstance(record, dict) or set(record) != {"command", "script_sha256", "stdout"}:
                raise GuardError("Contract-verifier provenance fields differ")
            command = shlex.split(str(record["command"]))
            expected_tail = [verifier_path, "--run-root", str(root)]
            if require_frozen:
                expected_tail.append("--require-bam-frozen")
            if (
                len(command) != len(expected_tail) + 1
                or not Path(command[0]).name.startswith("python")
                or command[1:] != expected_tail
                or record["script_sha256"] != observed_post["script_sha256"]
                or record["stdout"] != expected_stdout
            ):
                raise GuardError("Contract-verifier provenance differs on revalidation")

        # The historical preflight necessarily observed the draft state.
        # Re-running after publication auto-detects the frozen marker, so derive
        # the exact expected prior status line from the verified frozen status.
        expected_pre_stdout = re.sub(
            r"BAM manifest frozen$", "BAM manifest draft", observed_post["stdout"]
        )
        validate_verifier_record(
            preflight_record.get("contract_verifier"),
            require_frozen=False,
            expected_stdout=expected_pre_stdout,
        )
        validate_verifier_record(
            postflight_record.get("contract_verifier"),
            require_frozen=True,
            expected_stdout=observed_post["stdout"],
        )
    return postflight_record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fail-closed guard around the frozen BAM-manifest finalizer"
    )
    subparsers = parser.add_subparsers(dest="phase", required=True)
    for phase in ("pre", "post", "verify"):
        child = subparsers.add_parser(phase)
        child.add_argument("--run-root", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        if args.phase == "pre":
            record = preflight(args.run_root)
        elif args.phase == "post":
            record = postflight(args.run_root)
        else:
            record = verify_artifacts(args.run_root)
    except GuardError as exc:
        raise SystemExit(f"FAIL: {exc}") from exc
    if args.phase == "verify":
        print(f"PASS: BAM-manifest guard artifacts verified; run={record['run_id']}")
    else:
        artifact = PRE_ARTIFACT if args.phase == "pre" else POST_ARTIFACT
        print(f"PASS: BAM-manifest {args.phase}flight; wrote {artifact}; run={record['run_id']}")


if __name__ == "__main__":
    main()
