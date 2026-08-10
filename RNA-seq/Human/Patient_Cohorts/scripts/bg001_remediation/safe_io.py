#!/usr/bin/env python3
"""Minimal fail-closed publication helpers for immutable BG-001 records."""

from __future__ import annotations

import os
import secrets
import stat
from pathlib import Path


class SafeIOError(RuntimeError):
    pass


def _contained(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def require_safe_parent(path: Path, run_root: Path) -> None:
    if run_root.is_symlink():
        raise SafeIOError(f"Unsafe symlinked run root: {run_root}")
    run_root = run_root.resolve(strict=True)
    if not run_root.is_dir():
        raise SafeIOError(f"Unsafe run root: {run_root}")
    parent = path.parent
    resolved_parent = parent.resolve(strict=True)
    if not _contained(resolved_parent, run_root):
        raise SafeIOError(f"Publication parent escapes run root: {parent}")
    relative = parent.relative_to(run_root)
    cursor = run_root
    for component in relative.parts:
        cursor = cursor / component
        mode = os.lstat(cursor).st_mode
        if stat.S_ISLNK(mode) or not stat.S_ISDIR(mode):
            raise SafeIOError(f"Publication parent is symlinked or non-directory: {cursor}")


def publish_new_bytes(path: Path, payload: bytes, run_root: Path, mode: int = 0o640) -> None:
    """Publish a new file without following or overwriting a pre-existing target."""

    require_safe_parent(path, run_root)
    try:
        os.lstat(path)
    except FileNotFoundError:
        pass
    else:
        raise SafeIOError(f"Refusing existing publication target: {path}")

    temporary: Path | None = None
    descriptor: int | None = None
    for _ in range(20):
        candidate = path.parent / f".{path.name}.{secrets.token_hex(16)}.tmp"
        try:
            descriptor = os.open(candidate, os.O_WRONLY | os.O_CREAT | os.O_EXCL, mode)
            temporary = candidate
            break
        except FileExistsError:
            continue
    if descriptor is None or temporary is None:
        raise SafeIOError(f"Could not allocate exclusive temporary file for {path}")
    try:
        with os.fdopen(descriptor, "wb") as handle:
            descriptor = None
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        # Atomic no-overwrite publication: link() fails if a target appeared.
        os.link(temporary, path, follow_symlinks=False)
        os.unlink(temporary)
        temporary = None
        directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except BaseException:
        if descriptor is not None:
            os.close(descriptor)
        if temporary is not None:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass
        raise


def require_regular_file(path: Path) -> None:
    try:
        mode = os.lstat(path).st_mode
    except FileNotFoundError as exc:
        raise SafeIOError(f"Missing required regular file: {path}") from exc
    if stat.S_ISLNK(mode) or not stat.S_ISREG(mode):
        raise SafeIOError(f"Required path is symlinked or non-regular: {path}")
