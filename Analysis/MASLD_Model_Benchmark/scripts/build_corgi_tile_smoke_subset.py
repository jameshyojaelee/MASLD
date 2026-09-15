#!/usr/bin/env python3
"""Select a deterministic contig- and cCRE-balanced Corgi smoke tile roster."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence


FOLDS = tuple(range(5))
SEED = "corgi-smoke-20260824-v1"


class CorgiSmokeSubsetError(ValueError):
    """Raised when the frozen tile requirements or balanced selection differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        if not fields:
            raise CorgiSmokeSubsetError(f"{path.name} has no fields")
        return fields, [dict(row) for row in reader]


def write_tsv(path: Path, fields: Sequence[str], rows: list[dict[str, str]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def stable_key(*values: object) -> str:
    return sha256("|".join([SEED, *(str(value) for value in values)]).encode()).hexdigest()


def select_fold(
    tile_rows: list[dict[str, str]],
    window_rows: list[dict[str, str]],
    limit: int,
) -> list[dict[str, str]]:
    if not tile_rows or limit < 1 or len(tile_rows) < limit:
        raise CorgiSmokeSubsetError("fold tile roster cannot satisfy smoke limit")
    tile_ids = {row["tile_id"] for row in tile_rows}
    counts: dict[str, Counter[str]] = {tile_id: Counter() for tile_id in tile_ids}
    for row in window_rows:
        if row["tile_id"] not in counts:
            raise CorgiSmokeSubsetError("window refers outside fold tile roster")
        counts[row["tile_id"]][row["ccre_class"]] += 1
    dominant = {
        tile_id: sorted(values, key=lambda value: (-counts[tile_id][value], value))[0]
        for tile_id, values in counts.items()
    }
    by_stratum: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in tile_rows:
        by_stratum[(row["contig"], dominant[row["tile_id"]])].append(row)
    selected: dict[str, dict[str, str]] = {}
    for stratum in sorted(by_stratum):
        candidates = sorted(
            by_stratum[stratum], key=lambda row: stable_key(stratum[0], stratum[1], row["tile_id"])
        )
        selected[candidates[0]["tile_id"]] = candidates[0]
    if len(selected) > limit:
        raise CorgiSmokeSubsetError("mandatory contig-by-class strata exceed smoke limit")
    by_class: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tile_rows:
        if row["tile_id"] not in selected:
            by_class[dominant[row["tile_id"]]].append(row)
    for value in by_class.values():
        value.sort(key=lambda row: stable_key(dominant[row["tile_id"]], row["tile_id"]))
    selected_by_class = Counter(dominant[tile_id] for tile_id in selected)
    while len(selected) < limit:
        available = [key for key, value in by_class.items() if value]
        if not available:
            raise CorgiSmokeSubsetError("balanced selector exhausted candidates")
        chosen_class = min(available, key=lambda key: (selected_by_class[key], key))
        row = by_class[chosen_class].pop(0)
        selected[row["tile_id"]] = row
        selected_by_class[chosen_class] += 1
    return sorted(
        selected.values(),
        key=lambda row: (row["contig"], int(row["central_start"]), row["tile_id"]),
    )


def build(*, tile_contract: Path, output: Path, tiles_per_fold: int = 128) -> dict[str, object]:
    if output.exists() or tiles_per_fold < 1:
        raise CorgiSmokeSubsetError("invalid output or tile limit")
    tile_contract = tile_contract.resolve(strict=True)
    tile_fields, tiles = read_tsv(tile_contract / "tiles" / "tiles.tsv")
    window_fields, windows = read_tsv(tile_contract / "tiles" / "window_map.tsv")
    output.mkdir(mode=0o750)
    selected_tiles: list[dict[str, str]] = []
    per_fold: dict[str, dict[str, object]] = {}
    for fold in FOLDS:
        fold_tiles = [row for row in tiles if int(row["genomic_fold"]) == fold]
        fold_windows = [row for row in windows if int(row["genomic_fold"]) == fold]
        selected = select_fold(fold_tiles, fold_windows, tiles_per_fold)
        selected_tiles.extend(selected)
        selected_ids = {row["tile_id"] for row in selected}
        selected_windows = [row for row in fold_windows if row["tile_id"] in selected_ids]
        classes = Counter(row["ccre_class"] for row in selected_windows)
        per_fold[str(fold)] = {
            "tiles": len(selected),
            "windows": len(selected_windows),
            "contigs": len({row["contig"] for row in selected}),
            "ccre_class_windows": dict(sorted(classes.items())),
        }
    selected_ids = {row["tile_id"] for row in selected_tiles}
    selected_windows = [row for row in windows if row["tile_id"] in selected_ids]
    if len(selected_tiles) != tiles_per_fold * len(FOLDS) or not selected_windows:
        raise CorgiSmokeSubsetError("smoke tile or window census differs")
    if len({row["tile_id"] for row in selected_tiles}) != len(selected_tiles):
        raise CorgiSmokeSubsetError("smoke tile selection is not unique")
    write_tsv(output / "tiles.tsv", tile_fields, selected_tiles)
    write_tsv(output / "window_map.tsv", window_fields, selected_windows)
    summary: dict[str, object] = {
        "schema_version": "masld-bench-corgi-tile-smoke-subset-v1",
        "status": "pass_outcome_free_smoke_subset",
        "dataset_id": "gse296875",
        "selection_seed": SEED,
        "tiles_per_fold": tiles_per_fold,
        "tiles": len(selected_tiles),
        "windows": len(selected_windows),
        "per_fold": per_fold,
        "outcomes_or_labels_read": False,
        "model_predictions_read": False,
        "tile_contract_artifacts_sha256": digest(tile_contract / "ARTIFACTS.json"),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tile-contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tiles-per-fold", type=int, default=128)
    args = parser.parse_args()
    summary = build(
        tile_contract=args.tile_contract,
        output=args.output,
        tiles_per_fold=args.tiles_per_fold,
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
