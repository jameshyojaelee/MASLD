#!/usr/bin/env python3
"""Build candidate-only per-spot map sources for the frozen v2 programs.

One physical section/array per dataset is selected without program outcomes:
the source unit whose eligible spot count is closest to the dataset median,
with a lexical ID tie-break.  Both programs are then displayed in that same
unit.  Scores use the native v2 membership/weights and the pinned spatial
engine's hepatocyte-abundance, library-size, and detected-gene residualization.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from visium_rerun_lib import (
    RELEASE_ID,
    build_paths,
    load_legacy_engine,
    prepare_engine_inputs,
    read_tsv,
    sha256_file,
    utc_now,
    verify_hotspot_ready,
    verify_v1_anchors,
    write_tsv,
)


NATIVE_V2_READY_SHA256 = (
    "fcff53888cd8a29adc817752be4f3798e7fb3bcd7903f02622a90af729089898"
)
MAP_COLUMNS = (
    "release_id",
    "dataset",
    "source_dependence",
    "biological_unit_resolution",
    "reporting_unit_id",
    "source_individual_label",
    "condition",
    "program_uid",
    "legacy_program_id",
    "program_label",
    "membership_sha256",
    "spot_id",
    "x",
    "y",
    "residual_program_score_z",
    "n_genes_measured",
    "retained_l1_weight",
    "graph_eligible",
)


def deterministic_gzip_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(
                    text,
                    fieldnames=list(MAP_COLUMNS),
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="raise",
                )
                writer.writeheader()
                writer.writerows(rows)


def selected_units(sample_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    output = []
    datasets = sorted({row["dataset"] for row in sample_rows})
    for dataset in datasets:
        part = [row for row in sample_rows if row["dataset"] == dataset]
        counts = np.asarray([int(row["n_spots"]) for row in part], dtype=float)
        median = float(np.median(counts))
        chosen = min(
            part,
            key=lambda row: (abs(int(row["n_spots"]) - median), row["technical_id"]),
        )
        output.append(
            {
                "release_id": RELEASE_ID,
                "dataset": dataset,
                "reporting_unit_id": chosen["technical_id"],
                "source_individual_label": chosen["biological_id"],
                "condition": chosen["condition"],
                "n_spots_source": int(chosen["n_spots"]),
                "dataset_median_spot_count": f"{median:.12g}",
                "distance_to_median_spot_count": f"{abs(int(chosen['n_spots']) - median):.12g}",
                "selection_rule": "closest_to_dataset_median_source_spot_count_then_lexical_reporting_unit_id",
                "program_outcomes_used_for_selection": "FALSE",
                "same_unit_for_all_programs": "TRUE",
                "biological_unit_resolution": (
                    "resolved_human_donor" if dataset == "GSE192741" else "unresolved_physical_array"
                ),
            }
        )
    return output


def graph_indices(graphs: dict[str, list[tuple[np.ndarray, np.ndarray, np.ndarray]]]) -> np.ndarray:
    pieces = [component[0] for arrays in graphs.values() for component in arrays]
    if not pieces:
        raise RuntimeError("no graph-eligible tissue spots")
    return np.unique(np.concatenate(pieces).astype(int))


def score_dataset(paths, engine, registry, membership, dataset: str, unit: dict[str, object]):
    import anndata as ad

    adata = ad.read_h5ad(paths.datasets[dataset])
    try:
        adata.var_names = adata.var_names.astype(str)
        obs = adata.obs.copy()
        obs["sample_id"] = obs["sample_id"].astype(str)
        obs["individual"] = obs["individual"].astype(str)
        obs["condition"] = obs["condition"].astype(str)
        counts = engine.dense_counts(adata)
        library = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)
        norm = engine.log_normalize(counts, library)
        detection = np.asarray((counts > 0).mean(axis=0)).ravel()
        gene_to_index = {gene: index for index, gene in enumerate(adata.var_names)}

        factor_names = [
            str(value).replace("means_per_cluster_mu_fg_", "")
            for value in adata.uns["mod"]["factor_names"]
        ]
        if "Hepatocytes" not in factor_names:
            raise RuntimeError(f"{dataset} lacks Hepatocytes cell2location factor")
        abundance = np.asarray(adata.obsm["q05_cell_abundance_w_sf"])[
            :, factor_names.index("Hepatocytes")
        ]
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        graphs, _ = engine.section_graphs(obs, coords, k=6)
        eligible_all = graph_indices(graphs)
        selected_id = str(unit["reporting_unit_id"])
        selected_idx = np.asarray(
            sorted(
                {
                    index
                    for arrays in [graphs.get(selected_id, [])]
                    for component in arrays
                    for index in component[0]
                }
            ),
            dtype=int,
        )
        if len(selected_idx) < 50:
            raise RuntimeError(f"selected unit {dataset}/{selected_id} has <50 graph-eligible spots")

        rows: list[dict[str, object]] = []
        audits = []
        for program in registry.itertuples(index=False):
            program_members = membership.loc[
                membership["program_id"] == program.program_id
            ].copy()
            program_members = program_members.groupby("gene_symbol", as_index=False)[
                "original_l1_weight"
            ].sum()
            program_members["present"] = program_members["gene_symbol"].isin(gene_to_index)
            program_members["detected"] = program_members["gene_symbol"].map(
                lambda gene: (
                    detection[gene_to_index[gene]] >= 0.01
                    if gene in gene_to_index
                    else False
                )
            )
            measured = program_members.loc[
                program_members["present"] & program_members["detected"]
            ].copy()
            retained = float(measured["original_l1_weight"].sum())
            if len(measured) < 8 or retained < 0.20:
                raise RuntimeError(
                    f"map program {dataset}/{program.program_id} fails native testability"
                )
            weights = measured["original_l1_weight"].to_numpy(float) / retained
            gene_indices = np.asarray(
                [gene_to_index[gene] for gene in measured["gene_symbol"]], dtype=int
            )
            score = engine.extract_z(norm, gene_indices) @ weights
            design = engine.design_matrix(obs, abundance)
            residual = np.asarray(engine.residualize(score, design), dtype=float)
            center = float(np.nanmean(residual[eligible_all]))
            spread = float(np.nanstd(residual[eligible_all], ddof=1))
            if not np.isfinite(spread) or spread <= 0:
                raise RuntimeError(f"degenerate map score for {dataset}/{program.program_id}")
            zscore = (residual - center) / spread
            identity = next(
                row
                for row in read_tsv(
                    paths.hotspot_root / "program_registry_v2.tsv"
                )
                if row["program_uid"] == program.program_id
            )
            source_dependence = "independent" if dataset == "GSE192741" else "source_dependent"
            resolution = "resolved_human_donor" if dataset == "GSE192741" else "unresolved_physical_array"
            for index in selected_idx:
                rows.append(
                    {
                        "release_id": RELEASE_ID,
                        "dataset": dataset,
                        "source_dependence": source_dependence,
                        "biological_unit_resolution": resolution,
                        "reporting_unit_id": selected_id,
                        "source_individual_label": str(obs.iloc[index]["individual"]),
                        "condition": str(obs.iloc[index]["condition"]),
                        "program_uid": program.program_id,
                        "legacy_program_id": identity["cell_type"] + "::" + identity["module"],
                        "program_label": identity["module_name"],
                        "membership_sha256": identity["membership_sha256"],
                        "spot_id": str(obs.index[index]),
                        "x": f"{coords[index, 0]:.12g}",
                        "y": f"{coords[index, 1]:.12g}",
                        "residual_program_score_z": f"{zscore[index]:.12g}",
                        "n_genes_measured": len(measured),
                        "retained_l1_weight": f"{retained:.16g}",
                        "graph_eligible": "TRUE",
                    }
                )
            audits.append(
                {
                    "release_id": RELEASE_ID,
                    "dataset": dataset,
                    "reporting_unit_id": selected_id,
                    "program_uid": program.program_id,
                    "membership_sha256": identity["membership_sha256"],
                    "n_graph_eligible_spots": len(selected_idx),
                    "n_genes_measured": len(measured),
                    "retained_l1_weight": f"{retained:.16g}",
                    "residualization": "hepatocyte_abundance_log_library_size_detected_genes",
                    "display_standardization": "z_over_all_graph_eligible_spots_within_dataset",
                    "inferential_use": "FALSE_illustrative_map_only",
                }
            )
        return rows, audits
    finally:
        del adata


def build(project_root: Path | None) -> Path:
    paths = build_paths(project_root)
    verify_v1_anchors(paths)
    verify_hotspot_ready(paths)
    ready = paths.native_root / "v2_candidate/READY"
    if sha256_file(ready) != NATIVE_V2_READY_SHA256:
        raise RuntimeError("native v2 READY drift")
    ready_rows = read_tsv(ready)
    if len(ready_rows) != 1 or ready_rows[0].get("status") != "pass_v2_candidate":
        raise RuntimeError("native v2 candidate is not validated")

    final = paths.candidate_root / "map_source_v2"
    if final.exists():
        raise RuntimeError(f"refusing to overwrite map-source candidate: {final}")
    staging = final.with_name(f".{final.name}.incomplete.{os.environ.get('SLURM_JOB_ID', os.getpid())}")
    if staging.exists():
        raise RuntimeError(f"map-source staging path exists: {staging}")
    staging.mkdir(parents=True)

    engine = load_legacy_engine(paths)
    registry, membership, universe = prepare_engine_inputs(paths, "v2")
    if len(registry) != 2 or len(universe) != 2:
        raise RuntimeError("map source requires the complete two-program v2 family")
    sample_manifest_path = paths.native_root / "v2_candidate/native_sample_manifest.tsv"
    selections = selected_units(read_tsv(sample_manifest_path))
    if {row["dataset"] for row in selections} != set(paths.datasets):
        raise RuntimeError("map selection does not cover both native spatial datasets")

    map_rows = []
    audit_rows = []
    for selection in selections:
        rows, audits = score_dataset(
            paths,
            engine,
            registry,
            membership,
            str(selection["dataset"]),
            selection,
        )
        map_rows.extend(rows)
        audit_rows.extend(audits)
    map_rows.sort(
        key=lambda row: (
            str(row["dataset"]),
            str(row["program_uid"]),
            str(row["spot_id"]),
        )
    )
    deterministic_gzip_tsv(staging / "per_spot_program_map.tsv.gz", map_rows)
    write_tsv(
        staging / "map_selection.tsv",
        (
            "release_id",
            "dataset",
            "reporting_unit_id",
            "source_individual_label",
            "condition",
            "n_spots_source",
            "dataset_median_spot_count",
            "distance_to_median_spot_count",
            "selection_rule",
            "program_outcomes_used_for_selection",
            "same_unit_for_all_programs",
            "biological_unit_resolution",
        ),
        selections,
    )
    write_tsv(
        staging / "map_scoring_audit.tsv",
        (
            "release_id",
            "dataset",
            "reporting_unit_id",
            "program_uid",
            "membership_sha256",
            "n_graph_eligible_spots",
            "n_genes_measured",
            "retained_l1_weight",
            "residualization",
            "display_standardization",
            "inferential_use",
        ),
        audit_rows,
    )
    source_paths = [
        (Path(__file__).resolve(), "map_source_producer"),
        (Path(__file__).resolve().with_name("visium_rerun_lib.py"), "candidate_library"),
        (paths.v1_engine, "pinned_scientific_engine"),
        (ready, "native_v2_READY"),
        (sample_manifest_path, "native_v2_sample_manifest"),
        (paths.hotspot_root / "program_registry_v2.tsv", "v2_registry"),
        (paths.hotspot_root / "program_membership_v2.tsv", "v2_membership"),
    ] + [(path, f"{dataset}_source_h5ad") for dataset, path in paths.datasets.items()]
    write_tsv(
        staging / "source_manifest.tsv",
        ("relative_path", "source_role", "bytes", "sha256"),
        [
            {
                "relative_path": path.relative_to(paths.project_root).as_posix(),
                "source_role": role,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path, role in source_paths
        ],
    )
    write_tsv(
        staging / "execution_manifest.tsv",
        ("parameter", "value"),
        [
            {"parameter": "release_id", "value": RELEASE_ID},
            {"parameter": "task_id", "value": "SP-INT-06-map-source"},
            {"parameter": "selection_rule", "value": selections[0]["selection_rule"]},
            {"parameter": "program_outcomes_used_for_selection", "value": "FALSE"},
            {"parameter": "program_count", "value": len(registry)},
            {"parameter": "map_row_count", "value": len(map_rows)},
            {"parameter": "python", "value": platform.python_version()},
            {"parameter": "numpy", "value": np.__version__},
            {"parameter": "pandas", "value": pd.__version__},
            {"parameter": "built_utc", "value": utc_now()},
        ],
    )
    write_tsv(
        staging / "BUILD_COMPLETE",
        ("release_id", "status", "map_sha256", "row_count"),
        [
            {
                "release_id": RELEASE_ID,
                "status": "built_pending_independent_validation",
                "map_sha256": sha256_file(staging / "per_spot_program_map.tsv.gz"),
                "row_count": len(map_rows),
            }
        ],
    )
    os.replace(staging, final)
    print(f"built v2 spatial map source: rows={len(map_rows)} path={final}")
    return final


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        build(args.project_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
