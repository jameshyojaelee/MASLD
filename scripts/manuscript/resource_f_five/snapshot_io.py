#!/usr/bin/env python3
"""Fail-closed filesystem primitives for the Resource candidate snapshot."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import os
import re
import stat
import uuid
from pathlib import Path, PurePosixPath


class SnapshotIOError(RuntimeError):
    """Raised when a snapshot filesystem invariant is violated."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SnapshotIOError(message)


def _identity(value: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _inode_identity(value: os.stat_result) -> tuple[int, int]:
    return (value.st_dev, value.st_ino)


@dataclass(frozen=True)
class _RegularSnapshot:
    identity: tuple[int, int, int, int, int, int]
    mode: int
    size_bytes: int
    sha256: str


@dataclass(frozen=True)
class _ExclusiveCopyResult:
    size_bytes: int
    sha256: str
    inode: tuple[int, int]


@dataclass(frozen=True)
class _TreeEntry:
    entry_type: str
    identity: tuple[int, int, int, int, int, int]
    mode: int
    size_bytes: int
    sha256: str


@dataclass(frozen=True)
class _TreeSnapshot:
    root: _TreeEntry
    entries: tuple[tuple[str, _TreeEntry], ...]


def _directory_flags() -> int:
    _require(hasattr(os, "O_NOFOLLOW"), "O_NOFOLLOW is unavailable")
    _require(hasattr(os, "O_DIRECTORY"), "O_DIRECTORY is unavailable")
    return os.O_RDONLY | os.O_NOFOLLOW | os.O_DIRECTORY


def _open_verified_directory(path: Path) -> tuple[int, os.stat_result]:
    before = os.lstat(path)
    _require(stat.S_ISDIR(before.st_mode), f"not a directory: {path}")
    descriptor = os.open(path, _directory_flags())
    try:
        opened = os.fstat(descriptor)
        _require(
            stat.S_ISDIR(opened.st_mode) and _identity(opened) == _identity(before),
            f"directory changed before open: {path}",
        )
        after = os.lstat(path)
        _require(
            _identity(after) == _identity(before),
            f"directory changed after open: {path}",
        )
        return descriptor, opened
    except Exception:
        os.close(descriptor)
        raise


def _snapshot_regular(path: Path) -> _RegularSnapshot:
    before = os.lstat(path)
    _require(
        stat.S_ISREG(before.st_mode) and before.st_nlink == 1,
        f"not a singly linked regular file: {path}",
    )
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    digest = hashlib.sha256()
    try:
        opened = os.fstat(descriptor)
        _require(
            stat.S_ISREG(opened.st_mode)
            and opened.st_nlink == 1
            and _identity(opened) == _identity(before),
            f"file changed before open or has multiple links: {path}",
        )
        while True:
            block = os.read(descriptor, 4 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
        after_open = os.fstat(descriptor)
        _require(
            after_open.st_nlink == 1
            and _identity(after_open) == _identity(opened),
            f"file changed while hashed or has multiple links: {path}",
        )
    finally:
        os.close(descriptor)
    after = os.lstat(path)
    _require(
        after.st_nlink == 1 and _identity(after) == _identity(before),
        f"file changed after hash or has multiple links: {path}",
    )
    return _RegularSnapshot(
        identity=_identity(opened),
        mode=stat.S_IMODE(opened.st_mode),
        size_bytes=opened.st_size,
        sha256=digest.hexdigest(),
    )


def _snapshot_tree(path: Path) -> _TreeSnapshot:
    root_descriptor, root_opened = _open_verified_directory(path)
    os.close(root_descriptor)
    entries: list[tuple[str, _TreeEntry]] = []

    def visit(directory: Path, prefix: PurePosixPath) -> None:
        before_directory = os.lstat(directory)
        _require(
            stat.S_ISDIR(before_directory.st_mode),
            f"non-directory in publication tree: {directory}",
        )
        for child in sorted(directory.iterdir(), key=lambda item: item.name):
            child_stat = os.lstat(child)
            relative = (prefix / child.name).as_posix()
            if stat.S_ISDIR(child_stat.st_mode):
                child_descriptor, child_opened = _open_verified_directory(child)
                os.close(child_descriptor)
                entries.append(
                    (
                        relative,
                        _TreeEntry(
                            entry_type="directory",
                            identity=_identity(child_opened),
                            mode=stat.S_IMODE(child_opened.st_mode),
                            size_bytes=0,
                            sha256="",
                        ),
                    )
                )
                visit(child, prefix / child.name)
            elif stat.S_ISREG(child_stat.st_mode):
                regular = _snapshot_regular(child)
                entries.append(
                    (
                        relative,
                        _TreeEntry(
                            entry_type="file",
                            identity=regular.identity,
                            mode=regular.mode,
                            size_bytes=regular.size_bytes,
                            sha256=regular.sha256,
                        ),
                    )
                )
            else:
                raise SnapshotIOError(
                    f"symlink or special entry in publication tree: {child}"
                )
            child_after = os.lstat(child)
            _require(
                _identity(child_after) == _identity(child_stat),
                f"tree entry changed during inventory: {child}",
            )
        after_directory = os.lstat(directory)
        _require(
            _identity(after_directory) == _identity(before_directory),
            f"directory changed during inventory: {directory}",
        )

    visit(path, PurePosixPath())
    root_after = os.lstat(path)
    _require(
        _identity(root_after) == _identity(root_opened),
        f"tree root changed during inventory: {path}",
    )
    entries.sort(key=lambda item: item[0])
    return _TreeSnapshot(
        root=_TreeEntry(
            entry_type="directory",
            identity=_identity(root_opened),
            mode=stat.S_IMODE(root_opened.st_mode),
            size_bytes=0,
            sha256="",
        ),
        entries=tuple(entries),
    )


def _tree_content(
    snapshot: _TreeSnapshot,
    *,
    excluded: set[str] | None = None,
) -> tuple[int, tuple[tuple[str, str, int, int, str], ...]]:
    skipped = excluded or set()
    return (
        snapshot.root.mode,
        tuple(
            (
                relative,
                entry.entry_type,
                entry.mode,
                entry.size_bytes,
                entry.sha256,
            )
            for relative, entry in snapshot.entries
            if relative not in skipped
        ),
    )


def _relative_parts(value: str | Path) -> tuple[str, ...]:
    relative = PurePosixPath(str(value))
    _require(
        not relative.is_absolute()
        and bool(relative.parts)
        and all(part not in {"", ".", ".."} for part in relative.parts),
        f"marker must be a safe nonempty relative path: {value}",
    )
    _require(
        all("/" not in part and "\\" not in part and "\x00" not in part for part in relative.parts),
        f"unsafe marker path: {value}",
    )
    return relative.parts


def _reject_existing_destination(path: Path) -> None:
    try:
        existing = os.lstat(path)
    except FileNotFoundError:
        return
    if stat.S_ISLNK(existing.st_mode) or not (
        stat.S_ISREG(existing.st_mode) or stat.S_ISDIR(existing.st_mode)
    ):
        raise SnapshotIOError(
            f"publication destination is symlinked or special: {path}"
        )
    raise FileExistsError(f"publication destination exists: {path}")


def _copy_exclusive_with_identity(
    source: Path,
    destination: Path,
    *,
    mode: int,
) -> _ExclusiveCopyResult:
    before = os.lstat(source)
    _require(stat.S_ISREG(before.st_mode), f"source is not a regular file: {source}")
    source_flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    source_descriptor = os.open(source, source_flags)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination_descriptor: int | None = None
    digest = hashlib.sha256()
    copy_result: _ExclusiveCopyResult | None = None
    try:
        opened = os.fstat(source_descriptor)
        _require(_identity(opened) == _identity(before), f"source changed before open: {source}")
        destination_descriptor = os.open(
            destination,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
            0o600,
        )
        reader = os.fdopen(source_descriptor, "rb", closefd=True)
        source_descriptor = -1
        try:
            writer = os.fdopen(destination_descriptor, "wb", closefd=True)
            destination_descriptor = None
        except Exception:
            reader.close()
            raise
        with reader, writer:
            while True:
                block = reader.read(4 * 1024 * 1024)
                if not block:
                    break
                writer.write(block)
                digest.update(block)
            writer.flush()
            os.fchmod(writer.fileno(), mode)
            os.fsync(writer.fileno())
            copied_open = os.fstat(writer.fileno())
            _require(
                stat.S_ISREG(copied_open.st_mode)
                and copied_open.st_size == before.st_size
                and stat.S_IMODE(copied_open.st_mode) == mode
                and copied_open.st_nlink == 1,
                f"O_EXCL destination descriptor drift: {destination}",
            )
            copy_result = _ExclusiveCopyResult(
                size_bytes=copied_open.st_size,
                sha256=digest.hexdigest(),
                inode=_inode_identity(copied_open),
            )
        after = os.lstat(source)
        _require(_identity(after) == _identity(before), f"source changed while copied: {source}")
        copied = os.lstat(destination)
        _require(
            copy_result is not None
            and stat.S_ISREG(copied.st_mode)
            and _inode_identity(copied) == copy_result.inode
            and copied.st_nlink == 1
            and copied.st_size == copy_result.size_bytes
            and stat.S_IMODE(copied.st_mode) == mode,
            f"named copy differs from O_EXCL-created descriptor: {destination}",
        )
        return copy_result
    finally:
        if source_descriptor >= 0:
            os.close(source_descriptor)
        if destination_descriptor is not None:
            os.close(destination_descriptor)


def copy_exclusive(
    source: Path,
    destination: Path,
    *,
    mode: int = 0o440,
) -> tuple[int, str]:
    """Copy one regular file without following links or replacing a target."""

    copied = _copy_exclusive_with_identity(source, destination, mode=mode)
    return copied.size_bytes, copied.sha256


def _unlink_owned_staging(
    parent_descriptor: int,
    staging_name: str,
    owned_inode: tuple[int, int],
    *,
    expected_mode: int,
    expected_size: int,
    expected_nlink: int,
) -> None:
    current = os.stat(
        staging_name,
        dir_fd=parent_descriptor,
        follow_symlinks=False,
    )
    _require(
        stat.S_ISREG(current.st_mode)
        and _inode_identity(current) == owned_inode
        and stat.S_IMODE(current.st_mode) == expected_mode
        and current.st_size == expected_size
        and current.st_nlink == expected_nlink,
        f"owned staging file identity or link count drift: {staging_name}",
    )
    os.unlink(staging_name, dir_fd=parent_descriptor)
    os.fsync(parent_descriptor)


def _publish_file_by_staged_hardlink(source: Path, destination: Path) -> None:
    source_snapshot = _snapshot_regular(source)
    _reject_existing_destination(destination)
    parent_descriptor, _ = _open_verified_directory(destination.parent)
    staging_name = (
        f".{destination.name}.publish.{os.getpid()}.{uuid.uuid4().hex}.stage"
    )
    staging_path = destination.parent / staging_name
    staging_inode: tuple[int, int] | None = None
    try:
        copied = _copy_exclusive_with_identity(
            source,
            staging_path,
            mode=source_snapshot.mode,
        )
        staging_inode = copied.inode
        staging_snapshot = _snapshot_regular(staging_path)
        staged = os.stat(
            staging_name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        _require(
            stat.S_ISREG(staged.st_mode)
            and staged.st_nlink == 1
            and _inode_identity(staged) == staging_inode
            and staging_snapshot.identity[:2] == staging_inode
            and copied.size_bytes == source_snapshot.size_bytes
            and copied.sha256 == source_snapshot.sha256
            and staging_snapshot.mode == source_snapshot.mode
            and staging_snapshot.size_bytes == source_snapshot.size_bytes
            and staging_snapshot.sha256 == source_snapshot.sha256
            and staged.st_size == source_snapshot.size_bytes
            and stat.S_IMODE(staged.st_mode) == source_snapshot.mode,
            f"publication staging copy or independent snapshot drift: {destination}",
        )
        os.fsync(parent_descriptor)
        os.link(
            staging_name,
            destination.name,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        os.fsync(parent_descriptor)
        staged_after = os.stat(
            staging_name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        published = os.stat(
            destination.name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        _require(
            stat.S_ISREG(staged_after.st_mode)
            and stat.S_ISREG(published.st_mode)
            and _inode_identity(staged_after) == staging_inode
            and _inode_identity(published) == staging_inode
            and staged_after.st_nlink == 2
            and published.st_nlink == 2
            and stat.S_IMODE(staged_after.st_mode) == source_snapshot.mode
            and stat.S_IMODE(published.st_mode) == source_snapshot.mode
            and staged_after.st_size == source_snapshot.size_bytes
            and published.st_size == source_snapshot.size_bytes,
            f"hardlink publication identity, metadata, or link-count drift: {destination}",
        )
        _unlink_owned_staging(
            parent_descriptor,
            staging_name,
            staging_inode,
            expected_mode=source_snapshot.mode,
            expected_size=source_snapshot.size_bytes,
            expected_nlink=2,
        )
        published_after = os.stat(
            destination.name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        _require(
            stat.S_ISREG(published_after.st_mode)
            and _inode_identity(published_after) == staging_inode
            and published_after.st_nlink == 1
            and stat.S_IMODE(published_after.st_mode) == source_snapshot.mode
            and published_after.st_size == source_snapshot.size_bytes,
            f"final marker identity, metadata, or link-count drift: {destination}",
        )
        published_snapshot = _snapshot_regular(destination)
        _require(
            published_snapshot.mode == source_snapshot.mode
            and published_snapshot.size_bytes == source_snapshot.size_bytes
            and published_snapshot.sha256 == source_snapshot.sha256,
            f"published file independent content validation failed: {destination}",
        )
        _require(
            _snapshot_regular(source) == source_snapshot,
            f"publication source mutated: {source}",
        )
    except Exception:
        if staging_inode is not None:
            try:
                staged = os.stat(
                    staging_name,
                    dir_fd=parent_descriptor,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                pass
            else:
                if (
                    stat.S_ISREG(staged.st_mode)
                    and _inode_identity(staged) == staging_inode
                    and stat.S_IMODE(staged.st_mode) == source_snapshot.mode
                    and staged.st_size == source_snapshot.size_bytes
                    and staged.st_nlink in {1, 2}
                ):
                    _unlink_owned_staging(
                        parent_descriptor,
                        staging_name,
                        staging_inode,
                        expected_mode=source_snapshot.mode,
                        expected_size=source_snapshot.size_bytes,
                        expected_nlink=staged.st_nlink,
                    )
        raise
    finally:
        os.close(parent_descriptor)


def publish_file_noreplace(source: Path, destination: Path) -> dict[str, object]:
    """Publish one file via destination-local staging and a no-replace hardlink."""

    source = Path(source)
    destination = Path(destination)
    _publish_file_by_staged_hardlink(source, destination)
    published = _snapshot_regular(destination)
    return {
        "source": str(source),
        "destination": str(destination),
        "size_bytes": published.size_bytes,
        "sha256": published.sha256,
        "source_retained": True,
    }


def _fsync_verified_directory(path: Path) -> None:
    descriptor, _ = _open_verified_directory(path)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _require_reserved_root_named(
    parent_descriptor: int,
    destination_name: str,
    destination_descriptor: int,
    *,
    phase: str,
) -> None:
    held = os.fstat(destination_descriptor)
    named = os.stat(
        destination_name,
        dir_fd=parent_descriptor,
        follow_symlinks=False,
    )
    _require(
        stat.S_ISDIR(held.st_mode)
        and stat.S_ISDIR(named.st_mode)
        and _inode_identity(held) == _inode_identity(named),
        f"reserved destination root changed {phase}: {destination_name}",
    )


def publish_directory_noreplace(
    source: Path,
    destination: Path,
    *,
    marker: str | Path,
) -> dict[str, object]:
    """Copy a tree into a reserved path and hardlink its completion marker last."""

    source = Path(source)
    destination = Path(destination)
    marker_parts = _relative_parts(marker)
    marker_relative = PurePosixPath(*marker_parts).as_posix()
    source_snapshot = _snapshot_tree(source)
    source_entries = dict(source_snapshot.entries)
    _require(
        marker_relative in source_entries
        and source_entries[marker_relative].entry_type == "file",
        f"completion marker is absent or non-regular: {marker_relative}",
    )
    _reject_existing_destination(destination)
    parent_descriptor, _ = _open_verified_directory(destination.parent)
    destination_descriptor: int | None = None
    try:
        try:
            os.mkdir(destination.name, 0o700, dir_fd=parent_descriptor)
        except FileExistsError as error:
            raise FileExistsError(
                f"publication destination won reservation race: {destination}"
            ) from error
        os.fsync(parent_descriptor)
        destination_descriptor = os.open(
            destination.name,
            _directory_flags(),
            dir_fd=parent_descriptor,
        )
        reserved = os.fstat(destination_descriptor)
        named = os.stat(
            destination.name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        _require(
            stat.S_ISDIR(reserved.st_mode)
            and _identity(reserved) == _identity(named),
            f"reserved destination changed: {destination}",
        )

        directories = sorted(
            (
                (relative, entry)
                for relative, entry in source_snapshot.entries
                if entry.entry_type == "directory"
            ),
            key=lambda item: (len(PurePosixPath(item[0]).parts), item[0]),
        )
        for relative, _ in directories:
            target = destination.joinpath(*PurePosixPath(relative).parts)
            try:
                os.mkdir(target, 0o700)
            except FileExistsError as error:
                raise SnapshotIOError(
                    f"destination directory collision: {target}"
                ) from error

        for relative, entry in source_snapshot.entries:
            if entry.entry_type != "file" or relative == marker_relative:
                continue
            source_path = source.joinpath(*PurePosixPath(relative).parts)
            target = destination.joinpath(*PurePosixPath(relative).parts)
            copied_size, copied_hash = copy_exclusive(
                source_path,
                target,
                mode=entry.mode,
            )
            _require(
                copied_size == entry.size_bytes and copied_hash == entry.sha256,
                f"copied artifact drift: {relative}",
            )

        for relative, entry in sorted(
            directories,
            key=lambda item: len(PurePosixPath(item[0]).parts),
            reverse=True,
        ):
            target = destination.joinpath(*PurePosixPath(relative).parts)
            chmod_nofollow(target, entry.mode, directory=True)
            _fsync_verified_directory(target)
        os.fchmod(destination_descriptor, source_snapshot.root.mode)
        os.fsync(destination_descriptor)

        before_marker = _snapshot_tree(destination)
        _require(
            _tree_content(before_marker)
            == _tree_content(source_snapshot, excluded={marker_relative}),
            "destination inventory differs before completion-marker publication",
        )
        _require(
            _snapshot_tree(source) == source_snapshot,
            "source mutated before completion-marker publication",
        )
        _require_reserved_root_named(
            parent_descriptor,
            destination.name,
            destination_descriptor,
            phase="immediately before completion-marker publication",
        )
        marker_source = source.joinpath(*marker_parts)
        marker_destination = destination.joinpath(*marker_parts)
        _publish_file_by_staged_hardlink(marker_source, marker_destination)
        _require(
            _tree_content(_snapshot_tree(destination))
            == _tree_content(source_snapshot),
            "final destination inventory differs from source",
        )
        _require(
            _snapshot_tree(source) == source_snapshot,
            "source mutated during publication",
        )
        _require_reserved_root_named(
            parent_descriptor,
            destination.name,
            destination_descriptor,
            phase="after final inventory verification",
        )
        os.fsync(destination_descriptor)
        os.fsync(parent_descriptor)
    finally:
        if destination_descriptor is not None:
            os.close(destination_descriptor)
        os.close(parent_descriptor)

    return {
        "source": str(source),
        "destination": str(destination),
        "marker": marker_relative,
        "file_count": sum(
            entry.entry_type == "file" for _, entry in source_snapshot.entries
        ),
        "directory_count": 1
        + sum(entry.entry_type == "directory" for _, entry in source_snapshot.entries),
        "source_retained": True,
    }


def chmod_nofollow(path: Path, mode: int, *, directory: bool) -> None:
    """Change mode through a verified descriptor without following links."""

    _require(hasattr(os, "O_NOFOLLOW"), "O_NOFOLLOW is unavailable")
    before = os.lstat(path)
    expected = stat.S_ISDIR if directory else stat.S_ISREG
    _require(expected(before.st_mode), f"unexpected seal target type: {path}")
    flags = os.O_RDONLY | os.O_NOFOLLOW
    if directory:
        _require(hasattr(os, "O_DIRECTORY"), "O_DIRECTORY is unavailable")
        flags |= os.O_DIRECTORY
    descriptor = os.open(path, flags)
    try:
        opened = os.fstat(descriptor)
        _require(
            (opened.st_dev, opened.st_ino) == (before.st_dev, before.st_ino),
            f"seal target changed before open: {path}",
        )
        _require(expected(opened.st_mode), f"unexpected opened seal target: {path}")
        os.fchmod(descriptor, mode)
        sealed = os.fstat(descriptor)
        _require(
            stat.S_IMODE(sealed.st_mode) == mode,
            f"failed to seal target mode: {path}",
        )
    finally:
        os.close(descriptor)


def parse_sha256_manifest(
    path: Path,
    *,
    separator: str,
    expected_rows: int,
) -> list[tuple[str, str]]:
    """Parse a byte-stable lowercase SHA-256 manifest with one fixed separator."""

    raw = path.read_bytes()
    _require(raw.endswith(b"\n") and b"\r" not in raw, f"non-LF manifest: {path}")
    text = raw.decode("utf-8")
    pattern = re.compile(rf"^([0-9a-f]{{64}}){re.escape(separator)}(.+)$")
    rows: list[tuple[str, str]] = []
    seen: set[str] = set()
    for line in text.splitlines():
        match = pattern.fullmatch(line)
        _require(match is not None, f"malformed SHA-256 manifest row: {path}")
        digest, raw_path = match.groups()
        _require(raw_path not in seen, f"duplicate SHA-256 manifest path: {raw_path}")
        seen.add(raw_path)
        rows.append((digest, raw_path))
    _require(len(rows) == expected_rows, f"SHA-256 manifest row-count drift: {path}")
    return rows
