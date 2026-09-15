#!/usr/bin/env python3
"""Bring size companions into line with the files their pin describes.

A pin is rarely one field. Beside ``sha256`` a record often carries
``size_bytes`` describing the same file, and a refresh that moves only the
digest leaves the sibling asserting a size the file no longer has. The pin
contract test does not look at the companion, so the inconsistency is invisible
to it and surfaces later somewhere else.

This only ever writes a value re-derived from the file on disk, and only where
the record's own path field resolves.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

PREFIX = "Analysis/MASLD_Model_Benchmark/"
SIZE_KEYS = ("size_bytes", "bytes", "byte_count")


class CompanionRefreshError(RuntimeError):
    """Raised when a record cannot be reconciled with the file it describes."""


def resolve(root: Path, value: str) -> Path | None:
    relative = value[len(PREFIX):] if value.startswith(PREFIX) else value
    candidate = root / relative
    return candidate if candidate.is_file() else None


def stale_companions(root: Path, path: Path) -> list[dict]:
    """Every size companion in one record file that disagrees with disk."""

    found: list[dict] = []

    def walk(node: object, where: str = "") -> None:
        if isinstance(node, dict):
            named = next(
                (node[k] for k in ("identity", "authority_path", "path")
                 if isinstance(node.get(k), str)),
                None,
            )
            if named:
                target = resolve(root, named)
                if target is not None:
                    actual = target.stat().st_size
                    for key in SIZE_KEYS:
                        if key in node and node[key] != actual:
                            found.append({
                                "json_path": where,
                                "key": key,
                                "target": named,
                                "from": node[key],
                                "to": actual,
                                "digest_current": node.get("sha256") == sha256(
                                    target.read_bytes()).hexdigest(),
                            })
            for key, value in node.items():
                walk(value, f"{where}.{key}")
        elif isinstance(node, list):
            for index, value in enumerate(node):
                walk(value, f"{where}[{index}]")

    walk(json.loads(path.read_text(encoding="utf-8")))
    return found


def apply(root: Path, path: Path) -> list[dict]:
    document = json.loads(path.read_text(encoding="utf-8"))
    changes: list[dict] = []

    def walk(node: object, where: str = "") -> None:
        if isinstance(node, dict):
            named = next(
                (node[k] for k in ("identity", "authority_path", "path")
                 if isinstance(node.get(k), str)),
                None,
            )
            if named:
                target = resolve(root, named)
                if target is not None:
                    actual = target.stat().st_size
                    for key in SIZE_KEYS:
                        if key in node and node[key] != actual:
                            changes.append({"json_path": where, "key": key,
                                            "from": node[key], "to": actual,
                                            "target": named})
                            node[key] = actual
            for key, value in node.items():
                walk(value, f"{where}.{key}")
        elif isinstance(node, list):
            for index, value in enumerate(node):
                walk(value, f"{where}[{index}]")

    walk(document)
    if changes:
        path.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    return changes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--record", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        raise CompanionRefreshError("companion refresh output exists")

    receipt: dict = {
        "schema_version": "masld-bench-pin-companion-refresh-v1",
        "applied": bool(args.apply),
        "records": {},
    }
    for record in args.record:
        path = args.root / record
        if not path.is_file():
            raise CompanionRefreshError(f"record missing: {path}")
        planned = stale_companions(args.root, path)
        receipt["records"][str(record)] = {
            "stale_companions": planned,
            "changes": apply(args.root, path) if args.apply else [],
        }
    receipt["total_stale"] = sum(
        len(v["stale_companions"]) for v in receipt["records"].values()
    )

    # Rewriting a record changes its own digest, so anything pinning it is
    # stranded exactly as it would be by a pin refresh. The first version of
    # this script had no such check and stranded six pins that a branch had
    # just set. Report every holder that must be re-pinned; the caller is
    # responsible for completing the transaction.
    if args.apply:
        changed = [
            args.root / record
            for record, result in receipt["records"].items()
            if result["changes"]
        ]
        strays: list[str] = []
        for path in changed:
            marker = sha256(path.read_bytes()).hexdigest()
            for holder in (args.root / "config").rglob("*.json"):
                if holder == path:
                    continue
                try:
                    text = holder.read_text(encoding="utf-8")
                except OSError:
                    continue
                if path.name in text and marker not in text:
                    strays.append(f"{holder.relative_to(args.root)} may pin {path.name}")
        receipt["holders_to_reverify"] = sorted(set(strays))
        receipt["rewrite_changes_own_digest"] = True
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
