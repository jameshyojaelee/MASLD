#!/usr/bin/env python3
"""Build reviewable active-panel bindings and exact rendered-token classifications."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

from validate_rendered_figure_labels import (
    INVENTORY_FIELDS,
    TOKEN_FIELDS,
    default_paths,
    expected_token_row,
    inspect_scope,
    load_dynamic_roster,
    load_protected_static,
    load_scope,
    validate_inventory,
)


def atomic_write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--scope", type=Path)
    parser.add_argument("--roster", type=Path)
    parser.add_argument("--inventory", type=Path)
    parser.add_argument("--tokens", type=Path)
    parser.add_argument("--protected", type=Path)
    parser.add_argument(
        "--replace",
        action="store_true",
        help="Explicitly replace generated sidecars after reviewing active PDF changes.",
    )
    args = parser.parse_args()
    project_root = args.project_root.resolve(strict=True)
    defaults = default_paths(project_root)
    scope_path = args.scope or defaults[0]
    roster_path = args.roster or defaults[1]
    inventory_path = args.inventory or defaults[2]
    token_path = args.tokens or defaults[3]
    protected_path = args.protected or defaults[4]
    if not args.replace and (inventory_path.exists() or token_path.exists()):
        raise SystemExit("Refusing existing generated figure-label sidecars without --replace")

    scope_rows = load_scope(project_root, scope_path)
    scope_by_panel = {row["panel_id"]: row for row in scope_rows}
    load_dynamic_roster(project_root, roster_path, scope_by_panel)
    protected_static = load_protected_static(protected_path)
    _, inventory_rows, tokens_by_panel = inspect_scope(project_root, scope_rows)
    token_rows: list[dict[str, str]] = []
    for panel in scope_rows:
        for token in sorted(tokens_by_panel[panel["panel_id"]]):
            token_rows.append(expected_token_row(panel, token, protected_static))
    atomic_write_tsv(inventory_path, INVENTORY_FIELDS, inventory_rows)
    atomic_write_tsv(token_path, TOKEN_FIELDS, token_rows)
    result = validate_inventory(
        project_root,
        scope_path,
        roster_path,
        inventory_path,
        token_path,
        protected_path,
    )
    print(
        f"Wrote and validated {result['panels']} PDFs, {result['classified_rows']} token rows, "
        f"and {result['unique_dynamic_gene_labels']} unique dynamic gene labels"
    )


if __name__ == "__main__":
    main()
