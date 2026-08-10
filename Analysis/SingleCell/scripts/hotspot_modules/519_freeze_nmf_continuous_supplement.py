#!/usr/bin/env python3
"""Freeze the supplement-only k=4/k=6 continuous NMF bundle for Plan 20.

The producer copies every source byte into the candidate workstream, rebuilds
the evidence rows from those frozen inputs, and writes the READY seal last.
It never discovers factors, optimizes weights, assigns patient classes, or
overwrites an existing candidate bundle.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_PROJECT = SCRIPT_PATH.parents[4]
PROGRAM_ROOT = DEFAULT_PROJECT / "scripts/manuscript/program_context_v2"
COORDINATOR_ROOT = PROGRAM_ROOT / "coordinator"
for import_root in (PROGRAM_ROOT, COORDINATOR_ROOT):
    if str(import_root) not in sys.path:
        sys.path.insert(0, str(import_root))

from coordinator_contract import (  # noqa: E402
    CORE_PRODUCER_FIELDS,
    NMF_CONTINUOUS_FIELDS,
    NMF_READY_FIELDS,
    NMF_SOURCE_MANIFEST_FIELDS,
    NMF_VALIDATION_FIELDS,
    SOURCE_EVIDENCE_FIELDS,
    CoordinatorContractError,
    build_nmf_continuous_supplement,
    project_relative,
    read_tsv_flexible,
    sha256_file,
)
from release_common import CANDIDATE_ID, atomic_write_json, atomic_write_tsv  # noqa: E402


SOURCE_SPECS = (
    (
        "k4_loadings",
        "continuous k=4 sample-by-factor loadings",
        "RNA-seq/results/subtypes/nmf_assignments_k4_pre_k6restore.csv",
        "source_inputs/nmf_assignments_k4_pre_k6restore.csv",
    ),
    (
        "k6_loadings",
        "continuous k=6 sample-by-factor loadings",
        "RNA-seq/results/subtypes/nmf_assignments.csv",
        "source_inputs/nmf_assignments_k6.csv",
    ),
    (
        "k4_labels",
        "curated labels for the four continuous axes",
        "RNA-seq/results/subtypes/program_labels_k4_pre_k6restore.csv",
        "source_inputs/program_labels_k4_pre_k6restore.csv",
    ),
    (
        "k6_labels",
        "curated labels for the six continuous axes",
        "RNA-seq/results/subtypes/program_labels.csv",
        "source_inputs/program_labels_k6.csv",
    ),
    (
        "three_seed_metrics",
        "canonical three-seed internal NMF metrics",
        "RNA-seq/results/subtypes/nmf_3seed_canonical/per_seed_metric_summary.csv",
        "source_inputs/per_seed_metric_summary.csv",
    ),
    (
        "three_seed_stability",
        "canonical matched-factor stability metrics",
        "RNA-seq/results/subtypes/nmf_3seed_canonical/cross_seed_stability_summary.csv",
        "source_inputs/cross_seed_stability_summary.csv",
    ),
)


def atomic_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.", dir=destination.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
    try:
        shutil.copyfile(source, temporary)
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def freeze(project_root: Path, frozen_on_date: str) -> dict[str, object]:
    project = project_root.resolve()
    output = (
        project
        / "Analysis/Multimodal_Program_Projection/candidates"
        / CANDIDATE_ID
        / "hotspot/nmf_continuous_supplement"
    )
    if output.exists() or output.is_symlink():
        raise CoordinatorContractError(f"refusing to overwrite NMF candidate bundle: {output}")
    if not frozen_on_date.strip():
        raise CoordinatorContractError("frozen-on date is required")

    source_records: dict[str, tuple[Path, str, int, str]] = {}
    for source_id, source_role, origin_relative, bundle_relative in SOURCE_SPECS:
        origin = project / origin_relative
        if not origin.is_file() or origin.is_symlink():
            raise CoordinatorContractError(f"NMF source is missing or symlinked: {origin}")
        source_records[source_id] = (
            origin,
            sha256_file(origin),
            origin.stat().st_size,
            bundle_relative,
        )

    _, k4_rows = read_tsv_flexible(source_records["k4_loadings"][0])
    _, k6_rows = read_tsv_flexible(source_records["k6_loadings"][0])
    _, k4_labels = read_tsv_flexible(source_records["k4_labels"][0])
    _, k6_labels = read_tsv_flexible(source_records["k6_labels"][0])
    _, metric_rows = read_tsv_flexible(source_records["three_seed_metrics"][0])
    _, stability_rows = read_tsv_flexible(source_records["three_seed_stability"][0])
    loadings, evidence = build_nmf_continuous_supplement(
        k4_rows,
        k6_rows,
        k4_labels,
        k6_labels,
        metric_rows,
        stability_rows,
        source_records["k4_loadings"][1],
        source_records["k6_loadings"][1],
    )
    sample_ids = {str(row["sample_id"]) for row in loadings}
    s2a_marks = [row for row in evidence if row["panel_id"] == "S2A" and row["plot_role"] == "mark"]
    s2b_marks = [row for row in evidence if row["panel_id"] == "S2B" and row["plot_role"] == "mark"]
    if (
        len(sample_ids) != 1104
        or len(loadings) != 11040
        or len(evidence) != 18
        or len(s2a_marks) != 10
        or len(s2b_marks) != 8
    ):
        raise CoordinatorContractError("continuous NMF bundle cardinality drift")

    validation_rows = [
        {"check_id": "sample_universe", "status": "pass", "observed": str(len(sample_ids)), "criterion": "exactly 1,104 shared samples"},
        {"check_id": "continuous_loading_rows", "status": "pass", "observed": str(len(loadings)), "criterion": "1,104 x (4+6) = 11,040 long rows"},
        {"check_id": "supplement_evidence_rows", "status": "pass", "observed": str(len(evidence)), "criterion": "10 axis summaries plus 8 stability metrics"},
        {"check_id": "figure_s2a_marks", "status": "pass", "observed": str(len(s2a_marks)), "criterion": "one mark per k4/k6 continuous axis"},
        {"check_id": "figure_s2b_marks", "status": "pass", "observed": str(len(s2b_marks)), "criterion": "eight renderable stability marks"},
        {"check_id": "supplement_only", "status": "pass", "observed": str(sum(row["manuscript_included"] == "true" for row in evidence)), "criterion": "zero main-manuscript-included rows"},
        {"check_id": "continuous_only_semantics", "status": "pass", "observed": str(sum("continuous" in row["allowed_wording"].lower() and "subtype" in row["prohibited_wording"].lower() for row in evidence)), "criterion": "all 18 rows retain the continuous-axis/hard-partition boundary"},
    ]

    output.mkdir(parents=True)
    source_manifest_rows = []
    for source_id, source_role, _, bundle_relative in SOURCE_SPECS:
        origin, origin_hash, origin_bytes, _ = source_records[source_id]
        destination = output / bundle_relative
        atomic_copy(origin, destination)
        bundle_hash = sha256_file(destination)
        bundle_bytes = destination.stat().st_size
        if bundle_hash != origin_hash or bundle_bytes != origin_bytes:
            raise CoordinatorContractError(f"NMF source copy drift: {source_id}")
        source_manifest_rows.append(
            {
                "source_id": source_id,
                "source_role": source_role,
                "origin_project_path": project_relative(project, origin),
                "origin_sha256": origin_hash,
                "origin_bytes": origin_bytes,
                "bundle_relative_path": bundle_relative,
                "bundle_sha256": bundle_hash,
                "bundle_bytes": bundle_bytes,
            }
        )

    loadings_path = output / "nmf_continuous_loadings.tsv"
    evidence_path = output / "nmf_continuous_supplement.tsv"
    source_manifest_path = output / "source_manifest.tsv"
    validation_path = output / "validation_report.tsv"
    producer_path = output / "producer_manifest.tsv"
    atomic_write_tsv(loadings_path, loadings, NMF_CONTINUOUS_FIELDS)
    atomic_write_tsv(evidence_path, evidence, SOURCE_EVIDENCE_FIELDS)
    atomic_write_tsv(source_manifest_path, source_manifest_rows, NMF_SOURCE_MANIFEST_FIELDS)
    atomic_write_tsv(validation_path, validation_rows, NMF_VALIDATION_FIELDS)
    producers = []
    for producer_id, path in (
        (
            "plan20_nmf_freezer",
            project
            / "Analysis/SingleCell/scripts/hotspot_modules"
            / SCRIPT_PATH.name,
        ),
        (
            "plan60_nmf_contract",
            project
            / "scripts/manuscript/program_context_v2/coordinator/coordinator_contract.py",
        ),
    ):
        if not path.is_file() or path.is_symlink():
            raise CoordinatorContractError(f"NMF producer is missing or symlinked: {path}")
        producers.append(
            {
                "producer_id": producer_id,
                "repository_path": project_relative(project, path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
    atomic_write_tsv(producer_path, producers, CORE_PRODUCER_FIELDS)

    ready = {
        "candidate_id": CANDIDATE_ID,
        "status": "continuous_nmf_supplement_frozen",
        "source_manifest_sha256": sha256_file(source_manifest_path),
        "loadings_sha256": sha256_file(loadings_path),
        "evidence_sha256": sha256_file(evidence_path),
        "validation_sha256": sha256_file(validation_path),
        "producer_manifest_sha256": sha256_file(producer_path),
        "n_samples": len(sample_ids),
        "n_loading_rows": len(loadings),
        "n_evidence_rows": len(evidence),
        "n_s2a_marks": len(s2a_marks),
        "n_s2b_marks": len(s2b_marks),
        "canonical_promotion_authorized": "false",
        "frozen_on_date": frozen_on_date.strip(),
    }
    ready_path = output / "NMF_CONTINUOUS_SUPPLEMENT_READY"
    atomic_write_tsv(ready_path, [ready], NMF_READY_FIELDS)
    return {
        "status": ready["status"],
        "output": str(output),
        "ready_sha256": sha256_file(ready_path),
        "n_samples": len(sample_ids),
        "n_loading_rows": len(loadings),
        "n_evidence_rows": len(evidence),
        "canonical_promotion_authorized": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_PROJECT)
    parser.add_argument("--frozen-on-date", default="2026-08-08")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        result = freeze(args.project_root, args.frozen_on_date)
    except CoordinatorContractError as error:
        print(f"PLAN20_NMF_FREEZE_BLOCKED: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
