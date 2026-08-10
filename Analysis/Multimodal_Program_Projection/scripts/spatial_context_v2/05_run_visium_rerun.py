#!/usr/bin/env python3
"""Run the isolated v1-regression or frozen-v2 compatible Visium projection."""

from __future__ import annotations

import argparse
import gc
import hashlib
import os
import platform
import sys
from pathlib import Path

from visium_rerun_lib import (
    CORE_OUTPUTS,
    N_NULL,
    N_SENSITIVITY_NULL,
    RELEASE_ID,
    V1_ENGINE_SHA256,
    build_paths,
    load_legacy_engine,
    prepare_engine_inputs,
    read_tsv,
    sha256_file,
    universe_hash,
    utc_now,
    verify_hotspot_ready,
    verify_v1_anchors,
    write_tsv,
)


def require_freeze(paths) -> None:
    manifest = paths.native_root / "freeze/freeze_manifest.tsv"
    if not manifest.is_file():
        raise RuntimeError("SP-INT-03 is not frozen; run 04_freeze_visium_rerun.py write first")
    for row in read_tsv(manifest):
        path = paths.candidate_root / row["relative_path"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"SP-INT-03 freeze artifact is missing or drifted: {path}")
    for row in read_tsv(paths.native_root / "freeze/source_path_preflight.tsv"):
        path = paths.project_root / row["relative_path"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]):
            raise RuntimeError(f"SP-INT-03 source path or byte size drifted: {path}")
        if row["sha256_status"] == "frozen" and sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"SP-INT-03 frozen source hash drifted: {path}")


def require_v1_regression_ready(paths) -> str:
    ready_path = paths.native_root / "v1_regression/READY"
    if not ready_path.is_file():
        raise RuntimeError("v2 projection is blocked until the v1 regression validator writes READY")
    rows = read_tsv(ready_path)
    if len(rows) != 1:
        raise RuntimeError("v1 regression READY must contain exactly one row")
    row = rows[0]
    expected = {
        "release_id": RELEASE_ID,
        "registry_version": "v1",
        "status": "pass_v1_regression",
        "v1_engine_sha256": V1_ENGINE_SHA256,
    }
    for key, value in expected.items():
        if row.get(key) != value:
            raise RuntimeError(f"v1 regression READY {key} mismatch")
    report = paths.native_root / "v1_regression/validation_report.tsv"
    if row.get("validation_report_sha256") != sha256_file(report):
        raise RuntimeError("v1 regression validation report drifted after READY")
    return sha256_file(ready_path)


def source_manifest(paths, version: str, producer: Path, library: Path) -> list[dict[str, object]]:
    source_paths = [
        (paths.v1_engine, "pinned_scientific_engine"),
        (producer, "candidate_producer"),
        (library, "candidate_library"),
        (paths.gene_metadata, "gene_biotype_metadata"),
    ]
    if version == "v1":
        source_paths.extend(
            [(paths.v1_registry, "v1_registry"), (paths.v1_membership, "v1_membership")]
        )
    else:
        source_paths.extend(
            [
                (paths.hotspot_root / "program_registry_v2.tsv", "v2_registry"),
                (paths.hotspot_root / "program_membership_v2.tsv", "v2_membership"),
                (paths.hotspot_root / "READY", "v2_registry_ready"),
                (paths.native_root / "v1_regression/READY", "v1_regression_ready_gate"),
            ]
        )
    source_paths.extend((path, f"{dataset}_input_h5ad") for dataset, path in paths.datasets.items())
    rows = []
    for path, role in source_paths:
        if not path.is_file():
            raise RuntimeError(f"missing candidate source {role}: {path}")
        rows.append(
            {
                "relative_path": path.relative_to(paths.project_root).as_posix(),
                "source_role": role,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "observed_utc": utc_now(),
            }
        )
    return rows


def package_versions() -> dict[str, str]:
    import anndata
    import numpy
    import pandas
    import scipy

    return {
        "python": platform.python_version(),
        "anndata": anndata.__version__,
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "scipy": scipy.__version__,
    }


def canonical_columns(paths) -> dict[str, list[str]]:
    out = {}
    for name in CORE_OUTPUTS:
        with (paths.v1_results / name).open("r", encoding="utf-8") as handle:
            out[name] = handle.readline().rstrip("\n").split("\t")
    return out


def build_sample_and_design_audits(paths, version: str):
    import anndata as ad

    sample_rows = []
    design_rows = []
    for dataset, path in paths.datasets.items():
        adata = ad.read_h5ad(path, backed="r")
        try:
            required = {"sample_id", "individual", "condition"}
            missing = required.difference(adata.obs.columns)
            if missing:
                raise RuntimeError(f"{dataset} H5AD lacks sample-design fields: {sorted(missing)}")
            obs = adata.obs[["sample_id", "individual", "condition"]].astype(str).copy()
        finally:
            adata.file.close()
        for sample_id, part in obs.groupby("sample_id", sort=False, observed=False):
            individuals = sorted(set(part["individual"]))
            conditions = sorted(set(part["condition"]))
            if len(individuals) != 1 or len(conditions) != 1:
                raise RuntimeError(f"{dataset} sample {sample_id} has non-unique individual/condition")
            if dataset == "GSE192741":
                biological_unit = "human_donor"
                technical_unit = "Visium_section"
                resolution_note = "five_sections_collapsed_to_four_donors_H35_first"
            else:
                biological_unit = "physical_array"
                technical_unit = "Visium_array"
                resolution_note = "patient_barcode_key_unavailable_no_donor_upgrade"
            sample_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "registry_version": version,
                    "dataset": dataset,
                    "analysis_set_id": "compatible_visium_residual_moran",
                    "biological_id": individuals[0],
                    "technical_id": str(sample_id),
                    "biological_unit": biological_unit,
                    "technical_unit": technical_unit,
                    "n_spots": len(part),
                    "condition": conditions[0],
                    "include_primary": "TRUE",
                    "gate_state": "pass",
                    "resolution_note": resolution_note,
                }
            )
        dataset_rows = [row for row in sample_rows if row["dataset"] == dataset]
        biological_ids = sorted({str(row["biological_id"]) for row in dataset_rows})
        technical_ids = sorted({str(row["technical_id"]) for row in dataset_rows})
        design_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_version": version,
                "dataset": dataset,
                "analysis_set_id": "compatible_visium_residual_moran",
                "biological_unit": dataset_rows[0]["biological_unit"],
                "n_biological": len(biological_ids),
                "technical_unit": dataset_rows[0]["technical_unit"],
                "n_technical": len(technical_ids),
                "biological_ids_sha256": hashlib.sha256(("\n".join(biological_ids) + "\n").encode()).hexdigest(),
                "technical_ids_sha256": hashlib.sha256(("\n".join(technical_ids) + "\n").encode()).hexdigest(),
                "design_status": "pass_source_resolution_preserved",
            }
        )
    return sample_rows, design_rows


def build_program_audits(results, universe, version: str):
    universe_by_id = {str(row["program_id"]): row for row in universe}
    testability_rows = []
    sensitivity_rows = []
    for row in results.to_dict("records"):
        identity = universe_by_id[str(row["program_id"])]
        testable = bool(row["testable"])
        if testable:
            reason = "pass_n_genes_ge_8_and_retained_l1_ge_0.20"
        elif int(row["n_measured"]) < 8:
            reason = "fewer_than_8_measured_genes"
        else:
            reason = "retained_original_l1_weight_below_0.20"
        testability_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_version": version,
                "dataset": row["dataset"],
                "program_id": row["program_id"],
                "legacy_program_id": identity["legacy_program_id"],
                "membership_sha256": identity["membership_sha256"],
                "n_genes_measured": row["n_measured"],
                "retained_l1_weight": row["retained_l1_weight"],
                "testable": "TRUE" if testable else "FALSE",
                "testability_reason": reason,
            }
        )
        if not testable:
            continue
        primary_center = float(row["residual_moran_i"] - row["residual_null_mean"])
        equal_center = float(row["equal_residual_moran_i"] - row["equal_residual_null_mean"])
        leave_center = float(row["leave_top_residual_moran_i"] - row["leave_top_residual_null_mean"])
        for sensitivity_id, centered in (
            ("primary_original_l1_weight", primary_center),
            ("equal_weight", equal_center),
            ("leave_highest_weight_gene_out", leave_center),
        ):
            sensitivity_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "registry_version": version,
                    "dataset": row["dataset"],
                    "program_id": row["program_id"],
                    "membership_sha256": identity["membership_sha256"],
                    "sensitivity_id": sensitivity_id,
                    "centered_residual_moran_i": centered,
                    "direction": "positive" if centered > 0 else "negative" if centered < 0 else "zero",
                    "sign_agree_with_primary": "TRUE" if centered == 0 and primary_center == 0 or centered * primary_center > 0 else "FALSE",
                }
            )
    return testability_rows, sensitivity_rows


def write_output_manifest(staging: Path) -> None:
    rows = []
    excluded = {"output_manifest.tsv", "COMPUTE_COMPLETE", "READY", "validation_report.tsv"}
    for path in sorted(staging.iterdir()):
        if not path.is_file() or path.name in excluded:
            continue
        rows.append(
            {
                "relative_path": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    write_tsv(staging / "output_manifest.tsv", ("relative_path", "bytes", "sha256"), rows)


def run(paths, version: str) -> Path:
    import pandas as pd

    if version not in {"v1", "v2"}:
        raise ValueError("registry version must be v1 or v2")
    verify_v1_anchors(paths)
    require_freeze(paths)
    v1_ready_sha = "not_applicable"
    if version == "v2":
        verify_hotspot_ready(paths)
        v1_ready_sha = require_v1_regression_ready(paths)
    final = paths.native_root / ("v1_regression" if version == "v1" else "v2_candidate")
    if final.exists():
        raise RuntimeError(f"refusing to overwrite candidate output: {final}")
    staging = paths.native_root / f".{final.name}.incomplete.{os.environ.get('SLURM_JOB_ID', os.getpid())}"
    if staging.exists():
        raise RuntimeError(f"staging path already exists: {staging}")
    staging.mkdir(parents=True)

    producer = Path(__file__).resolve()
    library = producer.with_name("visium_rerun_lib.py")
    try:
        sources = source_manifest(paths, version, producer, library)
        write_tsv(
            staging / "source_manifest.tsv",
            ("relative_path", "source_role", "bytes", "sha256", "observed_utc"),
            sources,
        )
        engine = load_legacy_engine(paths)
        registry, membership, universe = prepare_engine_inputs(paths, version)
        if version == "v1" and len(registry) != 22:
            raise RuntimeError(f"v1 regression expected 22 programs, observed {len(registry)}")
        if version == "v2" and len(registry) != 2:
            raise RuntimeError(f"frozen v2 rerun expected two programs, observed {len(registry)}")

        gene_meta = pd.read_csv(
            paths.gene_metadata,
            sep="\t",
            usecols=["gene_name", "gene_biotype"],
            compression="gzip",
        ).drop_duplicates("gene_name")
        biotype = dict(
            zip(gene_meta["gene_name"].astype(str), gene_meta["gene_biotype"].astype(str))
        )
        all_results = []
        all_sections = []
        all_matches = []
        all_graphs = []
        all_nulls = []
        for dataset, path in paths.datasets.items():
            result, section, matching, graph, null = engine.process_dataset(
                dataset, path, registry, membership, biotype
            )
            all_results.append(result)
            all_sections.append(section)
            all_matches.append(matching)
            all_graphs.append(graph)
            all_nulls.append(null)
            gc.collect()

        columns = canonical_columns(paths)
        frames = {
            "spatial_program_results.tsv": pd.concat(all_results, ignore_index=True),
            "spatial_section_results.tsv": pd.concat(all_sections, ignore_index=True),
            "spatial_matching_audit.tsv": pd.concat(all_matches, ignore_index=True),
            "spatial_graph_audit.tsv": pd.concat(all_graphs, ignore_index=True),
            "spatial_null_summary.tsv": pd.concat(all_nulls, ignore_index=True),
        }
        for name, frame in frames.items():
            missing_columns = set(frame.columns).difference(columns[name])
            if missing_columns:
                raise RuntimeError(f"candidate engine emitted unexpected {name} columns: {sorted(missing_columns)}")
            frame.reindex(columns=columns[name]).to_csv(staging / name, sep="\t", index=False)

        sample_rows, design_rows = build_sample_and_design_audits(paths, version)
        write_tsv(
            staging / "native_sample_manifest.tsv",
            (
                "release_id",
                "registry_version",
                "dataset",
                "analysis_set_id",
                "biological_id",
                "technical_id",
                "biological_unit",
                "technical_unit",
                "n_spots",
                "condition",
                "include_primary",
                "gate_state",
                "resolution_note",
            ),
            sample_rows,
        )
        write_tsv(
            staging / "native_design_audit.tsv",
            (
                "release_id",
                "registry_version",
                "dataset",
                "analysis_set_id",
                "biological_unit",
                "n_biological",
                "technical_unit",
                "n_technical",
                "biological_ids_sha256",
                "technical_ids_sha256",
                "design_status",
            ),
            design_rows,
        )
        testability_rows, sensitivity_rows = build_program_audits(frames["spatial_program_results.tsv"], universe, version)
        write_tsv(
            staging / "spatial_testability_audit.tsv",
            (
                "release_id",
                "registry_version",
                "dataset",
                "program_id",
                "legacy_program_id",
                "membership_sha256",
                "n_genes_measured",
                "retained_l1_weight",
                "testable",
                "testability_reason",
            ),
            testability_rows,
        )
        write_tsv(
            staging / "spatial_sensitivity_audit.tsv",
            (
                "release_id",
                "registry_version",
                "dataset",
                "program_id",
                "membership_sha256",
                "sensitivity_id",
                "centered_residual_moran_i",
                "direction",
                "sign_agree_with_primary",
            ),
            sensitivity_rows,
        )

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
        write_tsv(staging / "tested_universe.tsv", universe_columns, universe)
        versions = package_versions()
        execution = [
            ("release_id", RELEASE_ID),
            ("task_id", "SP-INT-03"),
            ("registry_version", version),
            ("engine_strategy", "direct_call_to_byte_pinned_v1_functions"),
            ("v1_engine_sha256", V1_ENGINE_SHA256),
            ("candidate_producer_sha256", sha256_file(producer)),
            ("candidate_library_sha256", sha256_file(library)),
            ("tested_family_sha256", universe_hash(universe)),
            ("program_count", len(universe)),
            ("n_null", N_NULL),
            ("n_sensitivity_null", N_SENSITIVITY_NULL),
            ("dataset_order", ",".join(paths.datasets)),
            ("v1_regression_ready_sha256", v1_ready_sha),
            ("canonical_output_root", str(paths.v1_results)),
            ("candidate_output_root", str(final)),
            ("canonical_write_allowed", "FALSE"),
            ("started_or_completed_utc", utc_now()),
        ] + sorted(versions.items())
        write_tsv(
            staging / "execution_manifest.tsv",
            ("parameter", "value"),
            [{"parameter": key, "value": value} for key, value in execution],
        )
        write_output_manifest(staging)
        write_tsv(
            staging / "COMPUTE_COMPLETE",
            ("release_id", "registry_version", "status", "output_manifest_sha256", "completed_utc"),
            [
                {
                    "release_id": RELEASE_ID,
                    "registry_version": version,
                    "status": "compute_complete_pending_independent_validation",
                    "output_manifest_sha256": sha256_file(staging / "output_manifest.tsv"),
                    "completed_utc": utc_now(),
                }
            ],
        )
        os.replace(staging, final)
    except Exception:
        failure = staging / "FAILED"
        if staging.is_dir() and not failure.exists():
            write_tsv(
                failure,
                ("status", "failed_utc"),
                [{"status": "failed_before_candidate_promotion", "failed_utc": utc_now()}],
            )
        raise
    print(f"SP-INT-03 {version} compute complete: {final}", flush=True)
    return final


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry-version", required=True, choices=("v1", "v2"))
    parser.add_argument("--project-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        run(build_paths(args.project_root), args.registry_version)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
