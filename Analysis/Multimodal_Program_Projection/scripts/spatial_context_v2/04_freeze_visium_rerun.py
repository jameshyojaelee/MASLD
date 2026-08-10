#!/usr/bin/env python3
"""Freeze the outcome-blind SP-INT-03 compatible-Visium rerun contract."""

from __future__ import annotations

import argparse
import hashlib
import sys
from decimal import Decimal
from pathlib import Path

from visium_rerun_lib import (
    CORE_OUTPUTS,
    GRAPH_DISTANCE_MULTIPLIER,
    MIN_GRAPH_SPOTS,
    N_NULL,
    N_SENSITIVITY_NULL,
    RELEASE_ID,
    SEED,
    V1_ENGINE_SHA256,
    V1_MEMBERSHIP_SHA256,
    V1_REGISTRY_SHA256,
    V2_MEMBERSHIP_SHA256,
    V2_READY_SHA256,
    V2_REGISTRY_SHA256,
    build_paths,
    read_tsv,
    sha256_file,
    universe_hash,
    utc_now,
    v1_universe,
    v2_universe,
    verify_hotspot_ready,
    verify_v1_anchors,
    write_tsv,
)


FREEZE_FILES = (
    "rerun_specification.tsv",
    "tested_universe_v1.tsv",
    "tested_universe_v2.tsv",
    "expected_v1_outputs.tsv",
    "source_path_preflight.tsv",
)


def bool_text(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def membership_vectors(paths, version: str) -> dict[str, dict[str, Decimal]]:
    out: dict[str, dict[str, Decimal]] = {}
    if version == "v1":
        for row in read_tsv(paths.v1_membership):
            if row["mapped_symbol"].upper() != "TRUE" or not row["gene_symbol"]:
                continue
            vector = out.setdefault(row["program_id"], {})
            symbol = row["gene_symbol"]
            vector[symbol] = vector.get(symbol, Decimal(0)) + Decimal(row["original_l1_weight"])
        return out
    selected = {str(row["program_id"]) for row in v2_universe(paths)}
    for row in read_tsv(paths.hotspot_root / "program_membership_v2.tsv"):
        if row["program_uid"] not in selected:
            continue
        if row["mapped_symbol_status"] != "gencode_v49_unique_symbol_confirmed":
            continue
        vector = out.setdefault(row["program_uid"], {})
        symbol = row["mapped_symbol"]
        vector[symbol] = vector.get(symbol, Decimal(0)) + Decimal(row["original_l1_weight"])
    return out


def vector_hash(vector: dict[str, Decimal]) -> str:
    payload = "".join(f"{gene}\t{weight.normalize()}\n" for gene, weight in sorted(vector.items()))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def vectors_equal(left: dict[str, Decimal], right: dict[str, Decimal]) -> bool:
    if set(left) != set(right):
        return False
    tolerance = Decimal("1e-14")
    return all(abs(left[gene] - right[gene]) <= tolerance for gene in left)


def union_hash(vectors: dict[str, dict[str, Decimal]], ids: set[str]) -> tuple[int, str]:
    genes = sorted({gene for pid in ids for gene in vectors.get(pid, {})})
    return len(genes), hashlib.sha256(("\n".join(genes) + "\n").encode("utf-8")).hexdigest()


def registry_comparison(paths, v1_rows, v2_rows) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    v1_by_key = {(row["cell_type"], str(row["module"])): row for row in v1_rows}
    v2_by_key = {(row["cell_type"], str(row["module"])): row for row in v2_rows}
    v1_vectors = membership_vectors(paths, "v1")
    v2_vectors = membership_vectors(paths, "v2")
    v1_union_n, v1_union_sha = union_hash(v1_vectors, {str(row["program_id"]) for row in v1_rows})
    v2_union_n, v2_union_sha = union_hash(v2_vectors, {str(row["program_id"]) for row in v2_rows})
    family_v1 = universe_hash(v1_rows)
    family_v2 = universe_hash(v2_rows)
    comparison = []
    for cell_type, module in sorted(set(v1_by_key) | set(v2_by_key), key=lambda item: (item[0], int(item[1]))):
        old = v1_by_key.get((cell_type, module))
        new = v2_by_key.get((cell_type, module))
        same_weights = False
        if old and new:
            same_weights = vectors_equal(
                v1_vectors.get(str(old["program_id"]), {}),
                v2_vectors.get(str(new["program_id"]), {}),
            )
        comparison.append(
            {
                "cell_type": cell_type,
                "module": module,
                "v1_selected": bool_text(old is not None),
                "v1_program_id": "" if old is None else old["program_id"],
                "v1_membership_vector_sha256": "" if old is None else vector_hash(v1_vectors[str(old["program_id"])]),
                "v2_selected": bool_text(new is not None),
                "v2_program_uid": "" if new is None else new["program_id"],
                "v2_membership_sha256": "" if new is None else new["membership_sha256"],
                "v2_membership_vector_sha256": "" if new is None else vector_hash(v2_vectors[str(new["program_id"])]),
                "exact_mapped_membership_and_weight": bool_text(same_weights),
                "disposition": (
                    "retained_module_rerun_mapping_pool_and_family_changed"
                    if old and new
                    else "v1_only_not_in_v2_family"
                    if old
                    else "v2_new"
                ),
            }
        )
    reuse = []
    for row in v2_rows:
        legacy = v1_by_key.get((str(row["cell_type"]), str(row["module"])))
        same_weights = bool(
            legacy
            and vectors_equal(
                v1_vectors.get(str(legacy["program_id"]), {}),
                v2_vectors.get(str(row["program_id"]), {}),
            )
        )
        for dataset in paths.datasets:
            reuse.append(
                {
                    "release_id": RELEASE_ID,
                    "dataset": dataset,
                    "program_uid": row["program_id"],
                    "legacy_program_id": "" if legacy is None else legacy["program_id"],
                    "membership_and_weight_equal": bool_text(same_weights),
                    "measured_gene_equality": "deferred_until_candidate_run",
                    "matched_control_exclusion_pool_equal": "FALSE",
                    "tested_bh_family_equal": "FALSE",
                    "null_count_equal": "TRUE",
                    "producer_semantics_equal": "TRUE",
                    "reuse_eligible": "FALSE",
                    "decision": (
                        "rerun_required_mapping_program_union_and_bh_family_changed"
                        if not same_weights
                        else "rerun_required_program_union_and_bh_family_changed"
                    ),
                    "v1_program_union_n_genes": v1_union_n,
                    "v1_program_union_sha256": v1_union_sha,
                    "v2_program_union_n_genes": v2_union_n,
                    "v2_program_union_sha256": v2_union_sha,
                    "v1_tested_family_sha256": family_v1,
                    "v2_tested_family_sha256": family_v2,
                }
            )
    return comparison, reuse


def freeze(paths) -> None:
    verify_v1_anchors(paths)
    verify_hotspot_ready(paths)
    freeze_root = paths.native_root / "freeze"
    owned = [freeze_root / name for name in FREEZE_FILES]
    owned += [
        freeze_root / "freeze_manifest.tsv",
        paths.candidate_root / "registry_comparison.tsv",
        paths.candidate_root / "reuse_eligibility.tsv",
    ]
    existing = [path for path in owned if path.exists()]
    if existing:
        raise RuntimeError("refusing to overwrite SP-INT-03 freeze artifact(s): " + ", ".join(map(str, existing)))
    freeze_root.mkdir(parents=True, exist_ok=True)

    v1_rows = v1_universe(paths)
    v2_rows = v2_universe(paths)
    comparison, reuse = registry_comparison(paths, v1_rows, v2_rows)
    universe_columns = (
        "program_id",
        "legacy_program_id",
        "cell_type",
        "module",
        "program_name",
        "display_order",
        "membership_sha256",
        "primary_direction",
    )
    write_tsv(freeze_root / "tested_universe_v1.tsv", universe_columns, v1_rows)
    write_tsv(freeze_root / "tested_universe_v2.tsv", universe_columns, v2_rows)
    write_tsv(
        paths.candidate_root / "registry_comparison.tsv",
        (
            "cell_type",
            "module",
            "v1_selected",
            "v1_program_id",
            "v1_membership_vector_sha256",
            "v2_selected",
            "v2_program_uid",
            "v2_membership_sha256",
            "v2_membership_vector_sha256",
            "exact_mapped_membership_and_weight",
            "disposition",
        ),
        comparison,
    )
    write_tsv(
        paths.candidate_root / "reuse_eligibility.tsv",
        (
            "release_id",
            "dataset",
            "program_uid",
            "legacy_program_id",
            "membership_and_weight_equal",
            "measured_gene_equality",
            "matched_control_exclusion_pool_equal",
            "tested_bh_family_equal",
            "null_count_equal",
            "producer_semantics_equal",
            "reuse_eligible",
            "decision",
            "v1_program_union_n_genes",
            "v1_program_union_sha256",
            "v2_program_union_n_genes",
            "v2_program_union_sha256",
            "v1_tested_family_sha256",
            "v2_tested_family_sha256",
        ),
        reuse,
    )

    specification = [
        ("release_id", RELEASE_ID),
        ("task_id", "SP-INT-03"),
        ("engine_strategy", "byte_pinned_v1_functions_parameterized_candidate_outputs"),
        ("v1_engine_sha256", V1_ENGINE_SHA256),
        ("v1_registry_sha256", V1_REGISTRY_SHA256),
        ("v1_membership_sha256", V1_MEMBERSHIP_SHA256),
        ("v2_registry_sha256", V2_REGISTRY_SHA256),
        ("v2_membership_sha256", V2_MEMBERSHIP_SHA256),
        ("v2_ready_sha256", V2_READY_SHA256),
        ("v1_program_count", len(v1_rows)),
        ("v2_program_count", len(v2_rows)),
        ("v1_tested_family_sha256", universe_hash(v1_rows)),
        ("v2_tested_family_sha256", universe_hash(v2_rows)),
        ("n_null", N_NULL),
        ("n_sensitivity_null", N_SENSITIVITY_NULL),
        ("random_seed", SEED),
        ("graph_neighbors", 6),
        ("graph_distance_multiplier", GRAPH_DISTANCE_MULTIPLIER),
        ("minimum_graph_spots", MIN_GRAPH_SPOTS),
        ("primary_statistic", "lineage_residualized_donor_collapsed_moran_i"),
        ("gse192741_biological_unit", "four_donors_five_sections_H35_collapsed_first"),
        ("vu_biological_unit", "ten_physical_arrays_source_dependent"),
        ("multiple_testing", "BH_within_dataset_over_complete_prespecified_family"),
        ("external_gene_mapping", "gencode_v49_unique_symbol_confirmed_only_original_L1_denominator_preserved"),
        ("sensitivity_1", "equal_weight"),
        ("sensitivity_2", "leave_highest_weight_gene_out"),
        ("v2_start_gate", "validated_v1_regression_READY"),
        ("canonical_write_allowed", "FALSE"),
        ("yak_outcomes_allowed", "FALSE"),
        ("frozen_at_utc", utc_now()),
    ]
    write_tsv(
        freeze_root / "rerun_specification.tsv",
        ("parameter", "value"),
        [{"parameter": key, "value": value} for key, value in specification],
    )

    expected = []
    for name in CORE_OUTPUTS:
        path = paths.v1_results / name
        with path.open("r", encoding="utf-8") as handle:
            row_count = max(sum(1 for _ in handle) - 1, 0)
        expected.append(
            {
                "relative_path": path.relative_to(paths.project_root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "data_rows": row_count,
            }
        )
    write_tsv(
        freeze_root / "expected_v1_outputs.tsv",
        ("relative_path", "bytes", "sha256", "data_rows"),
        expected,
    )

    script_root = Path(__file__).resolve().parent
    source_paths = [
        (paths.v1_engine, "pinned_scientific_engine", True),
        (script_root / "01_v1_preservation.py", "v1_preservation_gate", True),
        (Path(__file__).resolve(), "candidate_freezer", True),
        (script_root / "visium_rerun_lib.py", "candidate_library", True),
        (script_root / "05_run_visium_rerun.py", "candidate_producer", True),
        (script_root / "06_validate_visium_rerun.py", "candidate_validator", True),
        (script_root / "tests/test_visium_rerun.py", "candidate_semantic_tests", True),
        (script_root / "run_visium_v1_regression.sbatch", "v1_regression_wrapper", True),
        (script_root / "run_visium_v2.sbatch", "v2_candidate_wrapper", True),
        (paths.v1_registry, "v1_tested_registry", True),
        (paths.v1_membership, "v1_program_membership", True),
        (paths.hotspot_root / "program_registry_v2.tsv", "v2_tested_registry", True),
        (paths.hotspot_root / "program_membership_v2.tsv", "v2_program_membership", True),
        (paths.hotspot_root / "READY", "v2_registry_seal", True),
        (paths.gene_metadata, "gene_biotype_metadata", True),
    ] + [(path, f"{dataset}_h5ad", False) for dataset, path in paths.datasets.items()]
    source_rows = []
    for path, role, hash_now in source_paths:
        if not path.is_file():
            raise RuntimeError(f"missing SP-INT-03 source: {path}")
        source_rows.append(
            {
                "relative_path": path.relative_to(paths.project_root).as_posix(),
                "source_role": role,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path) if hash_now else "",
                "sha256_status": "frozen" if hash_now else "required_at_compute_start",
            }
        )
    write_tsv(
        freeze_root / "source_path_preflight.tsv",
        ("relative_path", "source_role", "bytes", "sha256", "sha256_status"),
        source_rows,
    )

    manifest_rows = []
    manifest_paths = [freeze_root / name for name in FREEZE_FILES] + [
        paths.candidate_root / "registry_comparison.tsv",
        paths.candidate_root / "reuse_eligibility.tsv",
    ]
    for path in manifest_paths:
        manifest_rows.append(
            {
                "relative_path": path.relative_to(paths.candidate_root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    write_tsv(
        freeze_root / "freeze_manifest.tsv",
        ("relative_path", "bytes", "sha256"),
        manifest_rows,
    )
    print(
        f"froze SP-INT-03: v1={len(v1_rows)} programs; v2={len(v2_rows)} programs; "
        "reuse=0; v2 requires rerun",
        flush=True,
    )


def check(paths) -> None:
    verify_v1_anchors(paths)
    verify_hotspot_ready(paths)
    manifest_path = paths.native_root / "freeze/freeze_manifest.tsv"
    rows = read_tsv(manifest_path)
    if not rows:
        raise RuntimeError("SP-INT-03 freeze manifest is empty")
    for row in rows:
        path = paths.candidate_root / row["relative_path"]
        if not path.is_file():
            raise RuntimeError(f"frozen SP-INT-03 artifact is missing: {path}")
        if path.stat().st_size != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"frozen SP-INT-03 artifact drifted: {path}")
    for row in read_tsv(paths.native_root / "freeze/source_path_preflight.tsv"):
        path = paths.project_root / row["relative_path"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]):
            raise RuntimeError(f"SP-INT-03 source path or byte size drifted: {path}")
        if row["sha256_status"] == "frozen" and sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"SP-INT-03 frozen source hash drifted: {path}")
    v2_rows = read_tsv(paths.native_root / "freeze/tested_universe_v2.tsv")
    live = v2_universe(paths)
    if universe_hash(v2_rows) != universe_hash(live):
        raise RuntimeError("frozen and live v2 tested universes differ")
    if len(live) != 2:
        raise RuntimeError(f"reviewed SP-INT-03 freeze expected two robust programs, observed {len(live)}")
    reuse = read_tsv(paths.candidate_root / "reuse_eligibility.tsv")
    if len(reuse) != len(paths.datasets) * len(live) or any(row["reuse_eligible"] != "FALSE" for row in reuse):
        raise RuntimeError("reuse decision is incomplete or not fail-closed")
    print("SP-INT-03 freeze check passed", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("write", "check"))
    parser.add_argument("--project-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        paths = build_paths(args.project_root)
        if args.mode == "write":
            freeze(paths)
        else:
            check(paths)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
