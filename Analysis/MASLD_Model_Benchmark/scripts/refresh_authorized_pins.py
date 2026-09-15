#!/usr/bin/env python3
"""Refresh named stale SHA-256 pins, and refuse to touch anything else.

Refreshing a pin changes the containing bundle's own digest, which breaks every
file pinning that bundle.  So this never scans for work: it acts only on an
explicit authorization file naming each pin, its expected stale value, and the
live value it must become.  Any disagreement aborts before a byte is written.

The same shape drives the later cascade waves, where the ordering matters even
more.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import tomllib
from pathlib import Path
import re


class PinRefreshError(RuntimeError):
    """Raised when the tree disagrees with the authorization."""


DEFERRED = "RECOMPUTED_AFTER_SEEDS"
_INDEX = re.compile(r"^([^\[\]]*)\[(\d+)\]$")


def navigate(document: object, json_path: str) -> object:
    """Resolve a path like ``.audit_binding.registry`` or ``.source_records[0]``."""

    node = document
    for raw in json_path.lstrip(".").split("."):
        if not raw:
            continue
        match = _INDEX.match(raw)
        if match:
            key, index = match.group(1), int(match.group(2))
            if key:
                node = node[key]
            node = node[index]
        else:
            node = node[raw]
    return node


def _toml_table_and_key(json_path: str) -> tuple[str, str]:
    """Split ``.frozen_file_authorities.local_variant_config`` into table+name."""

    parts = [q for q in json_path.lstrip(".").split(".") if q]
    if len(parts) < 2:
        raise PinRefreshError(f"TOML pin path needs a table: {json_path!r}")
    return ".".join(parts), parts[-1]


def read_toml_pin(path: Path, json_path: str, key: str) -> str:
    """Read one key from a TOML table without rewriting anything."""

    table, _ = _toml_table_and_key(json_path)
    document = tomllib.loads(path.read_text(encoding="utf-8"))
    node = document
    for part in table.split("."):
        node = node[part]
    return node[key]


def write_toml_pin(path: Path, json_path: str, key: str, value: str) -> None:
    """Replace a single ``key = "..."`` line inside one TOML table.

    A tomllib/tomli_w round trip would rewrite the whole file and discard
    comments, ordering and spacing from a hand-maintained config. Only the one
    line is touched, and the result is re-parsed to prove it still loads and
    that the value landed where it was meant to.
    """

    table, _ = _toml_table_and_key(json_path)
    header = f"[{table}]"
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    inside = False
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("["):
            inside = stripped == header
            continue
        if inside and stripped.split("=")[0].strip() == key:
            indent = line[: len(line) - len(line.lstrip())]
            lines[index] = f'{indent}{key} = "{value}"\n'
            path.write_text("".join(lines), encoding="utf-8")
            if read_toml_pin(path, json_path, key) != value:
                raise PinRefreshError(f"TOML write did not land: {path} {json_path}")
            return
    raise PinRefreshError(f"{path}: no {key} under {header}")


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def plan_refresh(root: Path, authorization: dict) -> list[dict]:
    """Validate every authorized pin against the live tree before writing."""

    # A bundle that something else pins may only be refreshed as part of a
    # declared, branch-complete transaction. Isolated refresh of a pinned
    # bundle strands its pinners; a branch that names every pinner and updates
    # them in the same run does not. Declaring the closure up front is what
    # makes that property checkable rather than an accident of job structure.
    closure = set(authorization.get("branch_closure") or [])
    touched = {pin["bundle_path"] for pin in authorization["authorized_pins"]}
    for pin in authorization["authorized_pins"]:
        pinners = pin.get("bundle_is_pinned_by") or []
        if not pinners:
            continue
        if not closure:
            raise PinRefreshError(
                f"{pin['bundle_path']} is pinned by {pinners} and no "
                "branch_closure is declared; it belongs to a branch "
                "transaction, not an isolated refresh"
            )
        undeclared = [p for p in pinners if p not in closure]
        if undeclared:
            raise PinRefreshError(
                f"{pin['bundle_path']} is pinned by {undeclared}, which the "
                "declared branch_closure does not cover"
            )
        unrefreshed = [p for p in pinners if p not in touched]
        if unrefreshed:
            raise PinRefreshError(
                f"{pin['bundle_path']} is pinned by {unrefreshed}, which this "
                "run does not re-freeze; a branch must update every pinner in "
                "the same transaction"
            )

    actions: list[dict] = []
    for index, pin in enumerate(authorization["authorized_pins"]):
        bundle = root / pin["bundle_path"]
        if not bundle.is_file():
            raise PinRefreshError(f"bundle missing: {bundle}")
        if bundle.suffix == ".toml":
            recorded = read_toml_pin(bundle, pin["json_path"], pin["key"])
        else:
            document = json.loads(bundle.read_text(encoding="utf-8"))
            recorded = navigate(document, pin["json_path"])[pin["key"]]
        if recorded != pin["expected_stale_sha256"]:
            raise PinRefreshError(
                f"{pin['bundle_path']} {pin['json_path']}.{pin['key']} holds "
                f"{recorded}, not the authorized stale value "
                f"{pin['expected_stale_sha256']}"
            )
        target = root / pin["target_path"]
        if not target.is_file():
            raise PinRefreshError(f"pin target missing: {target}")

        # A pin whose target is itself refreshed earlier in this same run
        # cannot be pre-validated: its live digest does not exist yet. Verify
        # the ordering instead, and compute the value at apply time. Anything
        # whose target is external is still checked against the authorization
        # before a byte is written.
        deferred = pin["expected_live_sha256"] == DEFERRED
        if deferred:
            earlier = [
                q["bundle_path"] for q in authorization["authorized_pins"][:index]
            ]
            if pin["target_path"] not in earlier:
                raise PinRefreshError(
                    f"{pin['bundle_path']} {pin['json_path']} defers its live "
                    f"digest to {pin['target_path']}, which this run does not "
                    "refresh earlier; deferral requires the target to move first"
                )
            live = None
        else:
            live = digest(target)
            if live != pin["expected_live_sha256"]:
                raise PinRefreshError(
                    f"{pin['target_path']} is now {live}, not the authorized "
                    f"{pin['expected_live_sha256']}; the source moved again and "
                    "this authorization is stale"
                )
        actions.append(
            {
                "bundle_path": pin["bundle_path"],
                "json_path": pin["json_path"],
                "key": pin["key"],
                "target_path": pin["target_path"],
                "from": recorded,
                "to": live,
                "deferred": deferred,
                "bundle_sha256_before": digest(bundle),
            }
        )
    return actions


def apply_refresh(root: Path, actions: list[dict]) -> list[dict]:
    for action in actions:
        if action.get("deferred"):
            action["to"] = digest(root / action["target_path"])
        bundle = root / action["bundle_path"]
        if bundle.suffix == ".toml":
            write_toml_pin(bundle, action["json_path"], action["key"], action["to"])
        else:
            document = json.loads(bundle.read_text(encoding="utf-8"))
            navigate(document, action["json_path"])[action["key"]] = action["to"]
            bundle.write_text(
                json.dumps(document, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        action["bundle_sha256_after"] = digest(bundle)
    return actions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--authorization", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        raise PinRefreshError("refresh receipt output exists")

    authorization = json.loads(args.authorization.read_text(encoding="utf-8"))
    if not authorization.get("authorized_by"):
        raise PinRefreshError("authorization names no authorizer")
    actions = plan_refresh(args.root, authorization)

    applied = apply_refresh(args.root, actions) if args.apply else actions
    args.output.mkdir(parents=True, exist_ok=False)
    receipt = {
        "schema_version": "masld-bench-pin-refresh-receipt-v1",
        "authorization_id": authorization["authorization_id"],
        "authorized_by": authorization["authorized_by"],
        "applied": bool(args.apply),
        "pins": applied,
        "pins_refreshed": len(applied) if args.apply else 0,
        "cascade_started": False,
    }
    (args.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
