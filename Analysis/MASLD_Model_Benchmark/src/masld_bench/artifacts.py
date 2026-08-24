"""Immutable artifact helpers used by every benchmark component.

The control plane never overwrites a named artifact.  Callers create a new
campaign or attempt identifier when content changes.
"""

from __future__ import annotations

import ctypes
from dataclasses import dataclass
import errno
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Any, Iterable, Mapping

from .hashing import canonical_json, canonical_sha256, sha256_bytes, sha256_file


class ArtifactError(RuntimeError):
    """Raised when an immutable artifact cannot be created or verified."""


_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_AT_FDCWD = -100
_RENAME_NOREPLACE = 1


def reject_symlink_components(path: Path, *, label: str) -> Path:
    """Return an absolute lexical path after rejecting every existing symlink.

    ``Path.resolve`` must not run first because it erases the evidence that a
    configured output root or one of its parents was a symlink.
    """

    absolute = path if path.is_absolute() else Path.cwd() / path
    current = Path(absolute.anchor)
    for part in absolute.parts[1:]:
        current /= part
        if current.is_symlink():
            raise ArtifactError(f"{label} may not traverse a symlink: {current}")
    return absolute


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _renameat2_noreplace(source: Path, target: Path) -> bool:
    """Try Linux RENAME_NOREPLACE and report an unsupported filesystem."""

    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        return False
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    result = renameat2(
        _AT_FDCWD,
        os.fsencode(source),
        _AT_FDCWD,
        os.fsencode(target),
        _RENAME_NOREPLACE,
    )
    if result == 0:
        return True
    error_number = ctypes.get_errno()
    if error_number in {errno.EEXIST, errno.ENOTEMPTY}:
        raise ArtifactError(
            f"refusing to replace existing publication target: {target}"
        )
    unsupported = {
        errno.ENOSYS,
        errno.EINVAL,
        getattr(errno, "ENOTSUP", errno.EINVAL),
        getattr(errno, "EOPNOTSUPP", errno.EINVAL),
    }
    if error_number in unsupported:
        return False
    raise ArtifactError(
        f"atomic directory publication failed for {target}: "
        f"{os.strerror(error_number)}"
    )


def _reserved_directory_publication(source: Path, target: Path) -> None:
    """Publish on filesystems lacking RENAME_NOREPLACE without replacement.

    ``mkdir`` reserves the final name exclusively. Frozen trees become usable
    only when ``COMPLETE`` is linked as the final directory entry, so an
    interrupted publication is visible only as an invalid, fail-closed tree.
    This is the GPFS-compatible fallback for the artifact protocol.
    """

    source_stat = source.stat(follow_symlinks=False)
    parent_stat = target.parent.stat(follow_symlinks=False)
    if source_stat.st_dev != parent_stat.st_dev:
        raise ArtifactError(
            "directory publication requires source and target on one filesystem"
        )
    try:
        os.mkdir(target, mode=0o700)
    except FileExistsError as error:
        raise ArtifactError(
            f"refusing to replace existing publication target: {target}"
        ) from error

    moved: list[str] = []
    try:
        os.chmod(source, (source_stat.st_mode & 0o777) | 0o700)
        completion = source / "COMPLETE"
        has_completion = completion.is_file() and not completion.is_symlink()
        entries = sorted(
            (entry for entry in source.iterdir() if entry.name != "COMPLETE"),
            key=lambda entry: (entry.name == "ARTIFACTS.json", entry.name),
        )
        for entry in entries:
            os.rename(entry, target / entry.name)
            moved.append(entry.name)
        _fsync_directory(target)
        if has_completion:
            os.link(completion, target / "COMPLETE", follow_symlinks=False)
            _fsync_directory(target)
            completion.unlink()
        os.chmod(target, source_stat.st_mode & 0o777)
        source.rmdir()
        _fsync_directory(target.parent)
        if source.parent != target.parent:
            _fsync_directory(source.parent)
    except Exception as error:
        if not (target / "COMPLETE").exists():
            for name in reversed(moved):
                published = target / name
                if os.path.lexists(published) and not os.path.lexists(source / name):
                    try:
                        os.rename(published, source / name)
                    except OSError:
                        break
            try:
                target.rmdir()
            except OSError:
                pass
            if source.exists():
                try:
                    os.chmod(source, source_stat.st_mode & 0o777)
                except OSError:
                    pass
        if isinstance(error, ArtifactError):
            raise
        raise ArtifactError(
            f"reserved directory publication failed for {target}: {error}"
        ) from error


def publish_directory_noreplace(source: str | Path, target: str | Path) -> None:
    """Publish one directory without ever replacing an existing target."""

    source_path = reject_symlink_components(
        Path(source), label="publication source"
    )
    target_path = reject_symlink_components(
        Path(target), label="publication target"
    )
    reject_symlink_components(target_path.parent, label="publication parent")
    if not source_path.is_dir() or not target_path.parent.is_dir():
        raise ArtifactError(
            "directory publication requires an existing source and target parent"
        )
    if _renameat2_noreplace(source_path, target_path):
        _fsync_directory(target_path.parent)
        return
    _reserved_directory_publication(source_path, target_path)


def canonical_hash(value: Any) -> str:
    """Use the package-wide canonical identity algorithm."""

    return canonical_sha256(value)


def _atomic_link(payload: bytes, target: Path, mode: int) -> None:
    """Write fully, then atomically link into place without replacement."""

    target = reject_symlink_components(target, label="immutable artifact path")
    target.parent.mkdir(parents=True, exist_ok=True)
    reject_symlink_components(target, label="immutable artifact path")
    if not target.parent.is_dir():
        raise ArtifactError(f"immutable artifact parent is not a directory: {target.parent}")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, mode)
        try:
            os.link(temporary, target)
        except FileExistsError as error:
            raise ArtifactError(f"refusing to overwrite immutable artifact: {target}") from error
        directory_descriptor = os.open(target.parent, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    finally:
        temporary.unlink(missing_ok=True)


def write_bytes_exclusive(
    path: str | Path, payload: bytes, *, mode: int = 0o644
) -> "ArtifactRecord":
    target = Path(path)
    _atomic_link(payload, target, mode)
    return ArtifactRecord.from_path(target)


def write_text_exclusive(
    path: str | Path, text: str, *, mode: int = 0o644
) -> "ArtifactRecord":
    return write_bytes_exclusive(path, text.encode("utf-8"), mode=mode)


def write_json_exclusive(
    path: str | Path, value: Any, *, mode: int = 0o644
) -> "ArtifactRecord":
    return write_bytes_exclusive(path, (canonical_json(value) + "\n").encode("utf-8"), mode=mode)


@dataclass(frozen=True, slots=True)
class ArtifactRecord:
    path: str
    sha256: str
    size_bytes: int

    @classmethod
    def from_path(cls, path: str | Path, *, relative_to: str | Path | None = None) -> "ArtifactRecord":
        configured = reject_symlink_components(Path(path), label="artifact path")
        try:
            resolved = configured.resolve(strict=True)
        except OSError as error:
            raise ArtifactError(f"artifact is missing or unreadable: {configured}") from error
        if not resolved.is_file():
            raise ArtifactError(f"artifact is not a regular file: {resolved}")
        if relative_to is not None:
            base_configured = reject_symlink_components(
                Path(relative_to), label="artifact base path"
            )
            try:
                base = base_configured.resolve(strict=True)
                display_path = resolved.relative_to(base).as_posix()
            except (OSError, ValueError) as error:
                raise ArtifactError(
                    f"artifact does not reside below its declared base: {configured}"
                ) from error
        else:
            display_path = resolved.as_posix()
        before = resolved.stat(follow_symlinks=False)
        digest = sha256_file(resolved)
        after = resolved.stat(follow_symlinks=False)
        identity_before = (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
        )
        identity_after = (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
        )
        if identity_before != identity_after or resolved.is_symlink():
            raise ArtifactError(f"artifact changed while hashing: {resolved}")
        return cls(
            path=display_path,
            sha256=digest,
            size_bytes=after.st_size,
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "sha256": self.sha256,
            "size_bytes": self.size_bytes,
        }


def inventory_tree(
    root: str | Path,
    *,
    excluded_names: Iterable[str] = ("ARTIFACTS.json", "COMPLETE"),
) -> tuple[ArtifactRecord, ...]:
    configured = reject_symlink_components(
        Path(root), label="immutable-tree root"
    )
    try:
        base = configured.resolve(strict=True)
    except OSError as error:
        raise ArtifactError(f"artifact root does not exist: {configured}") from error
    exclusions = frozenset(excluded_names)
    if not base.is_dir():
        raise ArtifactError(f"artifact root does not exist: {base}")
    records = []
    for path in sorted(base.rglob("*")):
        if path.is_symlink():
            raise ArtifactError(f"symlinks are forbidden in immutable trees: {path}")
        relative = path.relative_to(base).as_posix()
        if path.is_file() and relative not in exclusions:
            try:
                path.resolve(strict=True).relative_to(base)
            except (OSError, ValueError) as error:
                raise ArtifactError(f"artifact escapes immutable tree: {path}") from error
            records.append(ArtifactRecord.from_path(path, relative_to=base))
    return tuple(records)


def freeze_tree(root: str | Path, metadata: Mapping[str, Any] | None = None) -> str:
    """Inventory a directory and add immutable manifest and completion marker."""

    configured = reject_symlink_components(
        Path(root), label="immutable-tree root"
    )
    try:
        base = configured.resolve(strict=True)
    except OSError as error:
        raise ArtifactError(f"artifact root does not exist: {configured}") from error
    records = inventory_tree(base)
    manifest = {
        "schema_version": "masld-bench-artifacts-v1",
        "metadata": dict(metadata or {}),
        "artifacts": [record.as_dict() for record in records],
    }
    manifest_record = write_json_exclusive(base / "ARTIFACTS.json", manifest)
    completion = {
        "schema_version": "masld-bench-complete-v1",
        "manifest_sha256": manifest_record.sha256,
        "artifact_count": len(records),
    }
    write_json_exclusive(base / "COMPLETE", completion)
    return manifest_record.sha256


def verify_frozen_tree(root: str | Path) -> dict[str, Any]:
    configured = reject_symlink_components(
        Path(root), label="immutable-tree root"
    )
    try:
        base = configured.resolve(strict=True)
    except OSError as error:
        raise ArtifactError(f"immutable-tree root does not exist: {configured}") from error
    manifest_path = base / "ARTIFACTS.json"
    completion_path = base / "COMPLETE"
    if manifest_path.is_symlink() or completion_path.is_symlink():
        raise ArtifactError("frozen-tree control files may not be symlinks")
    if not manifest_path.is_file() or not completion_path.is_file():
        raise ArtifactError(f"incomplete frozen artifact tree: {base}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ArtifactError(f"invalid frozen-tree control file: {error}") from error
    if not isinstance(manifest, Mapping) or not isinstance(completion, Mapping):
        raise ArtifactError("frozen-tree control files must contain JSON objects")
    if set(manifest) != {"schema_version", "metadata", "artifacts"}:
        raise ArtifactError("artifact manifest must contain exactly its registered fields")
    if manifest.get("schema_version") != "masld-bench-artifacts-v1":
        raise ArtifactError("unsupported artifact-manifest schema")
    if not isinstance(manifest.get("metadata"), Mapping):
        raise ArtifactError("artifact-manifest metadata must be an object")
    if set(completion) != {
        "schema_version",
        "manifest_sha256",
        "artifact_count",
    }:
        raise ArtifactError("completion marker must contain exactly its registered fields")
    if completion.get("schema_version") != "masld-bench-complete-v1":
        raise ArtifactError("unsupported completion-marker schema")
    raw_artifacts = manifest.get("artifacts")
    if not isinstance(raw_artifacts, list):
        raise ArtifactError("artifact manifest must contain an artifacts array")
    normalized_items: list[dict[str, Any]] = []
    seen_paths: set[str] = set()
    for raw in raw_artifacts:
        if not isinstance(raw, Mapping):
            raise ArtifactError("artifact manifest entries must be objects")
        if set(raw) != {"path", "sha256", "size_bytes"}:
            raise ArtifactError(
                "artifact manifest entries must contain exactly path, sha256, and size_bytes"
            )
        relative = Path(str(raw.get("path", "")))
        if relative.is_absolute() or not relative.parts or ".." in relative.parts:
            raise ArtifactError(f"unsafe artifact-manifest path: {relative}")
        relative_text = relative.as_posix()
        if relative_text in seen_paths:
            raise ArtifactError(f"duplicate artifact-manifest path: {relative_text}")
        seen_paths.add(relative_text)
        checksum = str(raw.get("sha256", ""))
        size = raw.get("size_bytes")
        if not _SHA256.fullmatch(checksum):
            raise ArtifactError(f"invalid artifact-manifest SHA-256: {relative_text}")
        if isinstance(size, bool) or not isinstance(size, int) or size < 0:
            raise ArtifactError(f"invalid artifact-manifest size: {relative_text}")
        normalized_items.append(
            {"path": relative_text, "sha256": checksum, "size_bytes": size}
        )
    actual_manifest_hash = sha256_file(manifest_path)
    completion_manifest_hash = str(completion.get("manifest_sha256", ""))
    if (
        not _SHA256.fullmatch(completion_manifest_hash)
        or completion_manifest_hash != actual_manifest_hash
    ):
        raise ArtifactError("completion marker does not match artifact manifest")
    if completion.get("artifact_count") != len(normalized_items):
        raise ArtifactError("completion marker artifact count does not match manifest")
    expected_paths = {item["path"] for item in normalized_items}
    observed_paths = {record.path for record in inventory_tree(base)}
    unexpected = sorted(observed_paths.difference(expected_paths))
    missing = sorted(expected_paths.difference(observed_paths))
    if unexpected or missing:
        details = []
        if unexpected:
            details.append("unexpected=" + ",".join(unexpected))
        if missing:
            details.append("missing=" + ",".join(missing))
        raise ArtifactError("frozen artifact inventory mismatch: " + "; ".join(details))
    for item in normalized_items:
        path = base / item["path"]
        if path.is_symlink() or not path.is_file():
            raise ArtifactError(f"missing frozen artifact: {item['path']}")
        before = path.stat(follow_symlinks=False)
        if before.st_size != item["size_bytes"]:
            raise ArtifactError(f"size mismatch for frozen artifact: {item['path']}")
        digest = sha256_file(path)
        after = path.stat(follow_symlinks=False)
        if (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
        ) != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
        ):
            raise ArtifactError(f"frozen artifact changed while hashing: {item['path']}")
        if digest != item["sha256"]:
            raise ArtifactError(f"checksum mismatch for frozen artifact: {item['path']}")
    return manifest
