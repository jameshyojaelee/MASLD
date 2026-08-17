#!/usr/bin/env python3
"""Freeze current program, metadata, atlas, and raw-transport inputs."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import anndata as ad
import pandas as pd

from common import (
    GLOBAL_ATLAS,
    HOTSPOT_ROOT,
    LINEAGES,
    MEMBERSHIP,
    OLD_RAW_EXPORT,
    PROGRAM_ROOT,
    READY,
    REGISTRY,
    SCORE_SCRIPT,
    require,
    refuse_existing,
    sha256,
)


CANDIDATE = Path(os.environ["CAND_ROOT"])
WORK = CANDIDATE / "work"
RESULTS = CANDIDATE / "results"


def load_score_producer():
    specification = importlib.util.spec_from_file_location("hotspot_501", SCORE_SCRIPT)
    require(
        specification is not None and specification.loader is not None,
        "cannot load frozen Hotspot producer",
    )
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def input_row(
    role: str, path: Path, known_hash: str | None = None
) -> dict[str, object]:
    require(path.exists(), f"missing input: {path}")
    return {
        "role": role,
        "path": str(path.resolve()),
        "bytes": int(path.stat().st_size),
        "sha256": known_hash if known_hash is not None else sha256(path),
    }


def main() -> None:
    manifest_path = WORK / "manifest.json"
    source_manifest_path = RESULTS / "input_manifest.tsv"
    refuse_existing(manifest_path)
    refuse_existing(source_manifest_path)

    ready = pd.read_csv(READY, sep="\t", dtype=str)
    require(len(ready) == 1, "READY seal malformed")
    require(ready.loc[0, "status"] == "ready_for_external_testing", "registry not ready")
    require(sha256(REGISTRY) == ready.loc[0, "registry_sha256"], "registry hash drift")
    require(
        sha256(MEMBERSHIP) == ready.loc[0, "membership_table_sha256"],
        "membership hash drift",
    )

    registry = pd.read_csv(REGISTRY, sep="\t")
    membership = pd.read_csv(MEMBERSHIP, sep="\t")
    require(len(registry) == 117, "frozen registry is not 117 programs")
    require(set(registry["cell_type"]) == set(LINEAGES), "registry lineage family drift")
    require(
        not registry[["cell_type", "module"]].duplicated().any(),
        "registry key is not unique",
    )
    require(
        set(membership["program_uid"]) == set(registry["program_uid"]),
        "membership program family drift",
    )

    upstream_work = OLD_RAW_EXPORT / "work"
    upstream_manifest = json.loads((upstream_work / "manifest.json").read_text())
    atlas_checksums = pd.read_csv(
        OLD_RAW_EXPORT / "results/01_input_checksums.tsv", sep="\t", dtype=str
    )
    recorded = atlas_checksums.loc[
        atlas_checksums["path"].str.endswith("results_gpu_v2/integrated_atlas.h5ad")
    ]
    require(len(recorded) == 1, "upstream atlas checksum row missing")
    current_atlas_hash = sha256(GLOBAL_ATLAS)
    require(
        current_atlas_hash == recorded.iloc[0]["sha256"],
        "current atlas differs from checksum-verified raw transport",
    )

    raw_gene_axis = (upstream_work / "genes.txt").read_text().splitlines()
    raw_gene_set = set(raw_gene_axis)
    score_genes = sorted(set(membership["source_gene"].astype(str)))
    missing_genes = sorted(set(score_genes) - raw_gene_set)
    require(not missing_genes, f"program genes absent from raw atlas: {missing_genes[:20]}")
    score_gene_path = WORK / "score_genes.txt"
    refuse_existing(score_gene_path)
    score_gene_path.write_text("\n".join(score_genes) + "\n", encoding="utf-8")

    producer = load_score_producer()
    cell_to_lineage: dict[str, str] = {}
    lineage_census: dict[str, dict[str, object]] = {}
    input_rows = [
        input_row("program_registry", REGISTRY),
        input_row("program_membership", MEMBERSHIP),
        input_row("program_ready_seal", READY),
        input_row("current_global_atlas", GLOBAL_ATLAS, current_atlas_hash),
        input_row("hotspot_score_producer", SCORE_SCRIPT),
    ]
    for lineage in LINEAGES:
        run_metadata_path = HOTSPOT_ROOT / lineage / "run_metadata.json"
        score_path = HOTSPOT_ROOT / lineage / "cell_scores.parquet"
        run_metadata = json.loads(run_metadata_path.read_text())
        atlas_path = producer.load_atlas.__globals__["CELLTYPE_SUBSETS"][lineage]
        atlas = ad.read_h5ad(atlas_path, backed="r")
        keep = pd.Series(True, index=atlas.obs_names)
        excluded = set(run_metadata.get("exclude_datasets") or [])
        if excluded:
            keep = ~atlas.obs["dataset"].astype(str).isin(excluded)
        cell_ids = atlas.obs_names[keep.to_numpy()].astype(str)
        require(
            len(cell_ids) == int(run_metadata["n_cells"]),
            f"frozen cell census drift for {lineage}",
        )
        collision = set(cell_ids).intersection(cell_to_lineage)
        if collision:
            raise RuntimeError(
                f"cell belongs to multiple scoring lineages: {sorted(collision)[0]}"
            )
        cell_to_lineage.update({cell_id: lineage for cell_id in cell_ids})
        lineage_census[lineage] = {
            "n_cells": int(len(cell_ids)),
            "n_programs": int((registry["cell_type"] == lineage).sum()),
            "excluded_datasets": sorted(excluded),
            "atlas_path": str(Path(atlas_path).resolve()),
        }
        input_rows.extend(
            [
                input_row(f"{lineage}_atlas", Path(atlas_path)),
                input_row(f"{lineage}_run_metadata", run_metadata_path),
                input_row(f"{lineage}_stored_scores", score_path),
            ]
        )
        del atlas

    found_cells: set[str] = set()
    dataset_manifest: dict[str, dict[str, object]] = {}
    for dataset in sorted(upstream_manifest["datasets"]):
        upstream_dataset = upstream_work / dataset
        meta = pd.read_csv(upstream_dataset / "meta.csv.gz", dtype=str)
        require("cell_id" in meta and "sample" in meta and "cell_type" in meta, "meta schema drift")
        meta["program_lineage"] = meta["cell_id"].map(cell_to_lineage).fillna("")
        selected = meta["program_lineage"] != ""
        found_cells.update(meta.loc[selected, "cell_id"])
        destination = WORK / dataset / "meta.csv.gz"
        refuse_existing(destination)
        meta.to_csv(destination, index=False)

        dims_path = upstream_dataset / "dims.json"
        dims = json.loads(dims_path.read_text())
        raw_files = [dims_path, upstream_dataset / "meta.csv.gz"]
        for chunk in dims["chunks"]:
            index = int(chunk["chunk"])
            raw_files.extend(
                [
                    upstream_dataset / f"chunk{index}_data.bin",
                    upstream_dataset / f"chunk{index}_indices.bin",
                    upstream_dataset / f"chunk{index}_indptr.bin",
                ]
            )
        for raw_file in raw_files:
            input_rows.append(input_row(f"upstream_raw_transport:{dataset}", raw_file))
        counts = meta.loc[selected, "program_lineage"].value_counts().sort_index()
        dataset_manifest[dataset] = {
            "n_cells": int(len(meta)),
            "n_program_universe_cells": int(selected.sum()),
            "program_lineage_cells": {key: int(value) for key, value in counts.items()},
            "upstream_dims": str(dims_path.resolve()),
            "candidate_meta": str(destination.resolve()),
        }

    require(
        found_cells == set(cell_to_lineage),
        f"raw transport misses {len(set(cell_to_lineage) - found_cells)} scoring cells",
    )

    input_table = pd.DataFrame(input_rows).drop_duplicates(subset=["role", "path"])
    input_table.to_csv(source_manifest_path, sep="\t", index=False)
    manifest = {
        "candidate_id": CANDIDATE.name,
        "seed": 20260815,
        "program_release": str(PROGRAM_ROOT.resolve()),
        "registry_sha256": sha256(REGISTRY),
        "membership_sha256": sha256(MEMBERSHIP),
        "global_atlas_sha256": current_atlas_hash,
        "upstream_raw_export": str(OLD_RAW_EXPORT.resolve()),
        "n_programs": int(len(registry)),
        "n_score_genes": int(len(score_genes)),
        "n_program_universe_cells": int(len(cell_to_lineage)),
        "lineages": lineage_census,
        "datasets": dataset_manifest,
        "analysis_universes": {
            "complete_case_common_universe": "exclude failed decontX datasets from raw and corrected arms",
            "full_universe_passthrough_sensitivity": "retain raw counts in corrected arm for failed datasets",
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "n_programs": manifest["n_programs"],
        "n_score_genes": manifest["n_score_genes"],
        "n_program_universe_cells": manifest["n_program_universe_cells"],
        "datasets": {key: value["n_program_universe_cells"] for key, value in dataset_manifest.items()},
    }, indent=2))


if __name__ == "__main__":
    main()
