#!/usr/bin/env python3
"""Freeze or recheck the immutable Figure 4 v1 byte-level baseline.

The protected census includes the v1 config/results trees, legacy v1 producer
scripts, the complete canonical Figure 4 release tree, and the six accepted
panel producers named by CANONICAL_MAIN_PANELS.tsv.  Candidate-v2 tooling and
Python bytecode are explicitly outside the v1 artifact set.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

from contract_lib import RELEASE_ID, sha256_file, write_tsv


INVENTORY_COLUMNS = (
    "release_id",
    "phase",
    "scope_category",
    "relative_path",
    "bytes",
    "sha256",
    "modified_utc",
)
COMPARISON_COLUMNS = (
    "release_id",
    "relative_path",
    "comparison_state",
    "baseline_bytes",
    "observed_bytes",
    "baseline_sha256",
    "observed_sha256",
)
SCOPE_COLUMNS = ("scope_category", "path_or_rule", "inclusion", "reason")
SUMMARY_COLUMNS = (
    "release_id",
    "phase",
    "status",
    "n_files",
    "n_unchanged",
    "n_changed",
    "n_missing",
    "n_added",
    "aggregate_sha256",
    "created_utc",
)
ANCHOR_COLUMNS = (
    "release_id",
    "n_files",
    "baseline_inventory_sha256",
    "aggregate_sha256",
    "scope_status",
)

# Outcome-blind SP-INT-01 anchors frozen on 2026-08-07.  Updating either value
# is a release-integrity decision, not an ordinary rerun accommodation.
EXPECTED_N_FILES = 591
EXPECTED_BASELINE_INVENTORY_SHA256 = "75d7dd76703f1091b74b709b2331d90e046b53a994ace87f06a1d1641268e956"
EXPECTED_AGGREGATE_SHA256 = "69edc1939fd4818a677fffc7a073d63d6b55aa61e176889bbe6d67910215cf30"

# Parallel Plan 13/30/40 candidate producers are not members of the reviewed
# v1 artifact set. Their exclusion is path-bounded; additions anywhere in the
# canonical config/results/Figure 4 trees remain a hard failure.
CANDIDATE_TOOL_DIRS = ("spatial_context_v2", "genetics_context_v2", "myojin_hlf")


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def utc_mtime(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_panel_producers(project_root: Path) -> list[Path]:
    manifest = project_root / "figures/main/fig4_validation/CANONICAL_MAIN_PANELS.tsv"
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 6:
        raise RuntimeError(f"expected six accepted Figure 4 panels, found {len(rows)}")
    producers = []
    for row in rows:
        if not row["status"].startswith("accepted"):
            raise RuntimeError(f"non-accepted row in canonical panel manifest: {row}")
        path = project_root / row["source"]
        if not path.is_file():
            raise RuntimeError(f"accepted panel producer is missing: {path}")
        producers.append(path)
    return sorted(set(producers))


def protected_paths(project_root: Path) -> tuple[list[tuple[str, Path]], list[dict[str, str]]]:
    roots = (
        ("v1_contract", project_root / "Analysis/Multimodal_Program_Projection/config"),
        ("v1_results", project_root / "Analysis/Multimodal_Program_Projection/results"),
        ("v1_producer", project_root / "Analysis/Multimodal_Program_Projection/scripts"),
        ("figure4_release_tree", project_root / "figures/main/fig4_validation"),
    )
    files: dict[Path, str] = {}
    scope_rows: list[dict[str, str]] = []
    for category, root in roots:
        if not root.is_dir():
            raise RuntimeError(f"v1 protected root is missing: {root}")
        scope_rows.append(
            {
                "scope_category": category,
                "path_or_rule": root.relative_to(project_root).as_posix() + "/**",
                "inclusion": "include",
                "reason": "immutable v1 artifact tree",
            }
        )
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            relative = path.relative_to(project_root).as_posix()
            if category == "v1_producer" and (
                any(f"/{directory}/" in f"/{relative}/" for directory in CANDIDATE_TOOL_DIRS)
                or "/__pycache__/" in f"/{relative}/"
                or path.suffix == ".pyc"
            ):
                continue
            files[path] = category

    scope_rows.extend(
        [
            *[
                {
                    "scope_category": "v2_candidate_tooling",
                    "path_or_rule": f"Analysis/Multimodal_Program_Projection/scripts/{directory}/**",
                    "inclusion": "exclude",
                    "reason": "path-bounded candidate tooling is not a v1 artifact",
                }
                for directory in CANDIDATE_TOOL_DIRS
            ],
            {
                "scope_category": "generated_bytecode",
                "path_or_rule": "**/__pycache__/** and *.pyc",
                "inclusion": "exclude",
                "reason": "non-source cache files are not stable release artifacts",
            },
        ]
    )

    for producer in read_panel_producers(project_root):
        files[producer] = "accepted_panel_producer"
        scope_rows.append(
            {
                "scope_category": "accepted_panel_producer",
                "path_or_rule": producer.relative_to(project_root).as_posix(),
                "inclusion": "include",
                "reason": "producer named by CANONICAL_MAIN_PANELS.tsv",
            }
        )
    return sorted(((category, path) for path, category in files.items()), key=lambda item: item[1].as_posix()), scope_rows


def inventory(project_root: Path, phase: str) -> list[dict[str, object]]:
    protected, _ = protected_paths(project_root)
    rows = []
    for category, path in protected:
        rows.append(
            {
                "release_id": RELEASE_ID,
                "phase": phase,
                "scope_category": category,
                "relative_path": path.relative_to(project_root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "modified_utc": utc_mtime(path),
            }
        )
    return rows


def aggregate_hash(rows: list[dict[str, object]]) -> str:
    payload = "".join(
        f"{row['relative_path']}\0{row['bytes']}\0{row['sha256']}\n"
        for row in sorted(rows, key=lambda row: str(row["relative_path"]))
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def read_inventory(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(INVENTORY_COLUMNS):
            raise RuntimeError(f"unexpected baseline schema in {path}")
        return list(reader)


def freeze(project_root: Path, candidate_root: Path) -> None:
    baseline_path = candidate_root / "v1_preservation_baseline.tsv"
    if baseline_path.exists():
        raise RuntimeError(f"refusing to overwrite immutable baseline: {baseline_path}")
    candidate_root.mkdir(parents=True, exist_ok=True)
    rows = inventory(project_root, "baseline")
    _, scope_rows = protected_paths(project_root)
    write_tsv(candidate_root / "v1_preservation_scope.tsv", SCOPE_COLUMNS, scope_rows)
    write_tsv(baseline_path, INVENTORY_COLUMNS, rows)
    observed_inventory_hash = sha256_file(baseline_path)
    observed_aggregate = aggregate_hash(rows)
    if (
        len(rows) != EXPECTED_N_FILES
        or observed_inventory_hash != EXPECTED_BASELINE_INVENTORY_SHA256
        or observed_aggregate != EXPECTED_AGGREGATE_SHA256
    ):
        raise RuntimeError(
            "new baseline does not match the reviewed SP-INT-01 hard anchors; "
            "do not update anchors merely to accommodate drift"
        )
    write_tsv(
        candidate_root / "v1_preservation_anchor.tsv",
        ANCHOR_COLUMNS,
        [
            {
                "release_id": RELEASE_ID,
                "n_files": len(rows),
                "baseline_inventory_sha256": observed_inventory_hash,
                "aggregate_sha256": observed_aggregate,
                "scope_status": "frozen_outcome_blind_v1_scope",
            }
        ],
    )
    write_tsv(
        candidate_root / "v1_preservation_summary.tsv",
        SUMMARY_COLUMNS,
        [
            {
                "release_id": RELEASE_ID,
                "phase": "baseline",
                "status": "baseline_frozen",
                "n_files": len(rows),
                "n_unchanged": len(rows),
                "n_changed": 0,
                "n_missing": 0,
                "n_added": 0,
                "aggregate_sha256": aggregate_hash(rows),
                "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            }
        ],
    )
    print(f"frozen {len(rows)} v1 files; aggregate={aggregate_hash(rows)}")


def check(project_root: Path, candidate_root: Path) -> None:
    baseline_path = candidate_root / "v1_preservation_baseline.tsv"
    if not baseline_path.is_file():
        raise RuntimeError(f"baseline does not exist: {baseline_path}")
    baseline = read_inventory(baseline_path)
    anchor_path = candidate_root / "v1_preservation_anchor.tsv"
    with anchor_path.open("r", encoding="utf-8", newline="") as handle:
        anchor_reader = csv.DictReader(handle, delimiter="\t")
        if anchor_reader.fieldnames != list(ANCHOR_COLUMNS):
            raise RuntimeError("unexpected v1 preservation anchor schema")
        anchor_rows = list(anchor_reader)
    if len(anchor_rows) != 1:
        raise RuntimeError("v1 preservation anchor must have exactly one row")
    anchor = anchor_rows[0]
    baseline_inventory_hash = sha256_file(baseline_path)
    baseline_aggregate = aggregate_hash(baseline)
    if (
        int(anchor["n_files"]) != EXPECTED_N_FILES
        or anchor["baseline_inventory_sha256"] != EXPECTED_BASELINE_INVENTORY_SHA256
        or anchor["aggregate_sha256"] != EXPECTED_AGGREGATE_SHA256
        or len(baseline) != EXPECTED_N_FILES
        or baseline_inventory_hash != EXPECTED_BASELINE_INVENTORY_SHA256
        or baseline_aggregate != EXPECTED_AGGREGATE_SHA256
    ):
        raise RuntimeError("v1 preservation baseline or anchor drifted from the reviewed hard anchors")
    _, scope_rows = protected_paths(project_root)
    write_tsv(candidate_root / "v1_preservation_scope.tsv", SCOPE_COLUMNS, scope_rows)
    observed = inventory(project_root, "final_check")
    write_tsv(candidate_root / "v1_preservation_final.tsv", INVENTORY_COLUMNS, observed)
    baseline_by_path = {row["relative_path"]: row for row in baseline}
    observed_by_path = {str(row["relative_path"]): row for row in observed}
    comparisons = []
    counts = {"unchanged": 0, "changed": 0, "missing": 0, "added": 0}
    for relative in sorted(set(baseline_by_path) | set(observed_by_path)):
        before = baseline_by_path.get(relative)
        after = observed_by_path.get(relative)
        if before is None:
            state = "added"
        elif after is None:
            state = "missing"
        elif before["sha256"] == after["sha256"] and int(before["bytes"]) == int(after["bytes"]):
            state = "unchanged"
        else:
            state = "changed"
        counts[state] += 1
        comparisons.append(
            {
                "release_id": RELEASE_ID,
                "relative_path": relative,
                "comparison_state": state,
                "baseline_bytes": "" if before is None else before["bytes"],
                "observed_bytes": "" if after is None else after["bytes"],
                "baseline_sha256": "" if before is None else before["sha256"],
                "observed_sha256": "" if after is None else after["sha256"],
            }
        )
    write_tsv(candidate_root / "v1_preservation_comparison.tsv", COMPARISON_COLUMNS, comparisons)
    status = "pass_byte_identical" if not (counts["changed"] or counts["missing"] or counts["added"]) else "failed_v1_drift"
    write_tsv(
        candidate_root / "v1_preservation_summary.tsv",
        SUMMARY_COLUMNS,
        [
            {
                "release_id": RELEASE_ID,
                "phase": "final_check",
                "status": status,
                "n_files": len(observed),
                "n_unchanged": counts["unchanged"],
                "n_changed": counts["changed"],
                "n_missing": counts["missing"],
                "n_added": counts["added"],
                "aggregate_sha256": aggregate_hash(observed),
                "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            }
        ],
    )
    print(f"{status}: {counts}; aggregate={aggregate_hash(observed)}")
    if status != "pass_byte_identical":
        raise RuntimeError("immutable Figure 4 v1 baseline changed")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("baseline", "check"))
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument(
        "--candidate-root",
        type=Path,
        default=None,
        help="defaults to the isolated spatial_context candidate root",
    )
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    candidate_root = args.candidate_root or (
        project_root
        / "Analysis/Multimodal_Program_Projection/candidates"
        / RELEASE_ID
        / "spatial_context"
    )
    try:
        if args.phase == "baseline":
            freeze(project_root, candidate_root)
        else:
            check(project_root, candidate_root)
    except Exception as exc:  # gate scripts should fail closed with a concise message
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
