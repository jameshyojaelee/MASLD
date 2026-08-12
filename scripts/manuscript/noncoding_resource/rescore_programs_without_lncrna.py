#!/usr/bin/env python3
"""Re-score frozen Hotspot programs after removing their lncRNA members."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path

import hotspot
import numpy as np
import pandas as pd
from hotspot import modules as hotspot_modules


ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
HS_ROOT = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
PROGRAM_ROOT = (
    ROOT
    / "RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/program_annotation/run_19749485"
)
REGISTRY = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv"
)
META = (
    ROOT
    / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
)
SCORE_501 = ROOT / "Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py"
HOTSPOT_IO = ROOT / "Analysis/SingleCell/scripts/hotspot_modules/hotspot_io.py"
PAIRINGS = {
    "GSE244832": ROOT / "data/GSE244832/metadata/donor_pairing.csv",
    "GSE185477": ROOT / "data/GSE185477/metadata/donor_pairing.csv",
    "GSE202379": ROOT / "data/GSE202379/metadata/donor_pairing.csv",
    "GSE136103": ROOT / "data/GSE136103/metadata/donor_pairing.csv",
}
PRIMARY_STAGES = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2}
CELL_TYPES = {"hepatocytes", "macrophages", "fibroblasts", "cholangiocytes", "tcells"}
SEED = 20260811
# Technical reconstruction floor calibrated from the outcome-blind comparison
# between reconstructed and stored frozen program scores. It is not the 0.90
# biological leave-lncRNA-out acceptance threshold.
MINIMUM_BASELINE_REPRODUCTION_R = 0.995


def require(value: bool, message: str) -> None:
    if not value:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_501():
    spec = importlib.util.spec_from_file_location("hotspot_run_501", SCORE_501)
    require(spec is not None and spec.loader is not None, "cannot load 501 producer")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def donor_map() -> dict[str, str]:
    output: dict[str, str] = {}
    for dataset, path in PAIRINGS.items():
        table = pd.read_csv(path, dtype=str)
        for row in table.itertuples(index=False):
            donor = f"{dataset}_{row.donor_id}"
            for sample in str(row.rna_srrs).split(";"):
                if sample.strip():
                    output[sample.strip()] = donor
    return output


def donor_metadata(mapping: dict[str, str]) -> pd.DataFrame:
    table = pd.read_csv(META, sep="\t", dtype=str)
    required = {"sample", "dataset", "disease_stage_coarse", "exclude_stage_analysis"}
    require(required <= set(table), "donor metadata schema drift")
    table = table[list(required)].copy()
    table["donor"] = table["sample"].map(mapping).fillna(table["sample"])
    table["exclude"] = (
        table["exclude_stage_analysis"].str.lower().isin({"true", "t", "1", "yes"})
    )
    rows = []
    for donor, group in table.groupby("donor", sort=False):
        datasets = sorted(set(group["dataset"].dropna()) - {""})
        stages = sorted(set(group["disease_stage_coarse"].dropna()) - {""})
        require(len(datasets) <= 1 and len(stages) <= 1, f"metadata conflict: {donor}")
        rows.append(
            {
                "donor": donor,
                "dataset": datasets[0] if datasets else None,
                "stage": stages[0] if stages else None,
                "exclude": bool(group["exclude"].any()),
            }
        )
    return pd.DataFrame(rows)


def stage_beta(scores: pd.DataFrame, metadata: pd.DataFrame) -> tuple[float, int, int]:
    frame = scores.merge(metadata, on="donor", how="inner", validate="one_to_one")
    frame = frame[(~frame["exclude"]) & frame["stage"].isin(PRIMARY_STAGES)].copy()
    frame["stage_ordinal"] = frame["stage"].map(PRIMARY_STAGES).astype(float)
    datasets = pd.get_dummies(frame["dataset"], drop_first=True, dtype=float)
    design = np.column_stack(
        [
            np.ones(len(frame)),
            frame["stage_ordinal"].to_numpy(float),
            datasets.to_numpy(float),
        ]
    )
    require(
        np.linalg.matrix_rank(design) == design.shape[1], "rank-deficient stage model"
    )
    beta = np.linalg.lstsq(design, frame["score"].to_numpy(float), rcond=None)[0][1]
    return float(beta), len(frame), int(frame["dataset"].nunique())


def score_genes(hs: hotspot.Hotspot, adata, genes: list[str]) -> np.ndarray:
    require(len(genes) >= 2, "fewer than two genes remain after lncRNA removal")
    counts = hs._counts_from_anndata(adata[:, genes], hs.layer_key, dense=True)
    return hotspot_modules.compute_scores(
        counts,
        hs.model,
        hs.umi_counts.values,
        hs.neighbors.values,
        hs.weights.values,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-type", required=True, choices=sorted(CELL_TYPES))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    require(not args.output.exists(), "output already exists")
    np.random.seed(SEED)

    content_path = PROGRAM_ROOT / "program_lncrna_content.tsv"
    annotation_path = PROGRAM_ROOT / "program_membership_lncrna_annotation.tsv"
    content = pd.read_csv(content_path, sep="\t")
    annotation = pd.read_csv(annotation_path, sep="\t", dtype=str)
    registry = pd.read_csv(REGISTRY, sep="\t")
    run_metadata_path = HS_ROOT / args.cell_type / "run_metadata.json"
    run_metadata = json.loads(run_metadata_path.read_text(encoding="utf-8"))
    trigger = content[
        (content["cell_type"] == args.cell_type)
        & (
            content["requires_leave_all_lncrna_out_sensitivity"].astype(str).str.lower()
            == "true"
        )
    ].copy()
    require(len(trigger) > 0, "cell type has no triggered programs")

    run501 = load_501()
    adata = run501.load_atlas(args.cell_type, smoke=False)
    excluded_datasets = run_metadata.get("exclude_datasets") or []
    if excluded_datasets:
        adata = adata[~adata.obs["dataset"].isin(set(excluded_datasets))].copy()
    adata = run501.strip_confounders(adata)
    adata = run501.filter_detected(adata, min_frac=0.01)
    require(adata.n_obs == int(run_metadata["n_cells"]), "frozen run cell census drift")
    require(
        adata.n_vars == int(run_metadata["n_genes_kept"]),
        "frozen run gene census drift",
    )
    run501.ensure_raw_layer(adata)
    latent = run501.resolve_latent(adata)
    hs = hotspot.Hotspot(
        adata,
        layer_key="counts",
        model="danb",
        latent_obsm_key=latent,
        umi_counts_obs_key="n_counts",
    )
    neighbors = 50 if args.cell_type == "macrophages" else 30
    hs.create_knn_graph(weighted_graph=False, n_neighbors=neighbors)

    stored = pd.read_parquet(HS_ROOT / args.cell_type / "cell_scores.parquet")
    require(stored["cell_id"].nunique() == adata.n_obs, "stored cell universe drift")
    mapping = donor_map()
    metadata = donor_metadata(mapping)
    cell_donor = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "donor": adata.obs["sample"]
            .astype(str)
            .map(mapping)
            .fillna(adata.obs["sample"].astype(str)),
        }
    )
    registry_ct = registry[registry["cell_type"] == args.cell_type].set_index("module")

    summary_rows = []
    donor_rows = []
    for row in trigger.itertuples(index=False):
        module = int(row.module)
        members = annotation[
            (annotation["cell_type"] == args.cell_type)
            & (annotation["module"].astype(int) == module)
        ]
        all_genes = members["source_gene"].tolist()
        reduced_genes = members.loc[
            members["gene_type"] != "lncRNA", "source_gene"
        ].tolist()
        require(
            set(all_genes) <= set(adata.var_names.astype(str)),
            f"missing module genes: {module}",
        )
        original = score_genes(hs, adata, all_genes)
        reduced = score_genes(hs, adata, reduced_genes)

        stored_module = stored[stored["module"].astype(int) == module].set_index(
            "cell_id"
        )["score"]
        stored_aligned = stored_module.reindex(adata.obs_names.astype(str))
        require(not stored_aligned.isna().any(), f"stored score join failed: {module}")
        reproduction_r = float(
            np.corrcoef(original, stored_aligned.to_numpy(float))[0, 1]
        )
        max_reproduction_error = float(
            np.max(np.abs(original - stored_aligned.to_numpy(float)))
        )
        print(
            f"REPRODUCTION\t{args.cell_type}\t{module}\t"
            f"r={reproduction_r:.12g}\tmax_abs_error={max_reproduction_error:.12g}",
            flush=True,
        )
        require(
            reproduction_r >= MINIMUM_BASELINE_REPRODUCTION_R,
            f"original Hotspot score not reproduced: {module}",
        )
        # The source serialized exact module scores, but not the neighbor graph.
        # Preserve the exact frozen score and add only the gene-removal delta
        # estimated on the faithfully reconstructed graph.
        reduced_delta_corrected = stored_aligned.to_numpy(float) + (reduced - original)
        reduced_r = float(
            np.corrcoef(stored_aligned.to_numpy(float), reduced_delta_corrected)[0, 1]
        )

        cell_frame = cell_donor.copy()
        cell_frame["original_score"] = stored_aligned.to_numpy(float)
        cell_frame["without_lncrna_score"] = reduced_delta_corrected
        donor = (
            cell_frame.groupby("donor", as_index=False)[
                ["original_score", "without_lncrna_score"]
            ]
            .mean()
            .sort_values("donor")
        )
        original_beta, n_donors, n_datasets = stage_beta(
            donor[["donor", "original_score"]].rename(
                columns={"original_score": "score"}
            ),
            metadata,
        )
        reduced_beta, reduced_n, reduced_datasets = stage_beta(
            donor[["donor", "without_lncrna_score"]].rename(
                columns={"without_lncrna_score": "score"}
            ),
            metadata,
        )
        require(
            (n_donors, n_datasets) == (reduced_n, reduced_datasets),
            "model census drift",
        )
        stored_beta = float(registry_ct.loc[module, "primary_beta"])
        require(
            np.isclose(original_beta, stored_beta, atol=1e-8, rtol=1e-8),
            f"stage beta drift: {module}",
        )
        direction_preserved = np.sign(reduced_beta) == np.sign(original_beta)
        passed = reduced_r >= 0.90 and direction_preserved
        summary_rows.append(
            {
                "program_uid": row.program_uid,
                "cell_type": args.cell_type,
                "module": module,
                "n_original_genes": len(all_genes),
                "n_removed_lncrna": len(all_genes) - len(reduced_genes),
                "lncrna_l1_fraction": row.lncrna_l1_fraction,
                "score_reproduction_r": reproduction_r,
                "score_reproduction_max_abs_error": max_reproduction_error,
                "without_lncrna_score_r": reduced_r,
                "stored_stage_beta": stored_beta,
                "reconstructed_stage_beta": original_beta,
                "without_lncrna_stage_beta": reduced_beta,
                "stage_direction_preserved": direction_preserved,
                "n_stage_donors": n_donors,
                "n_stage_datasets": n_datasets,
                "sensitivity_pass": passed,
            }
        )
        donor.insert(0, "program_uid", row.program_uid)
        donor_rows.append(donor)

    args.output.mkdir(parents=True)
    summary = pd.DataFrame(summary_rows).sort_values("program_uid")
    donor_table = pd.concat(donor_rows, ignore_index=True).sort_values(
        ["program_uid", "donor"]
    )
    summary.to_csv(
        args.output / "program_without_lncrna_sensitivity.tsv", sep="\t", index=False
    )
    donor_table.to_csv(args.output / "donor_scores.tsv", sep="\t", index=False)
    inputs = [
        content_path,
        annotation_path,
        REGISTRY,
        META,
        SCORE_501,
        HOTSPOT_IO,
        HS_ROOT / args.cell_type / "cell_scores.parquet",
        run_metadata_path,
        *PAIRINGS.values(),
    ]
    pd.DataFrame(
        [
            {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in inputs
        ]
    ).to_csv(args.output / "source_manifest.tsv", sep="\t", index=False)
    (args.output / "execution_manifest.json").write_text(
        json.dumps(
            {
                "cell_type": args.cell_type,
                "n_triggered_programs": len(summary),
                "n_cells": int(adata.n_obs),
                "n_genes_after_original_filter": int(adata.n_vars),
                "n_neighbors": neighbors,
                "excluded_datasets": excluded_datasets,
                "score_method": "stored_Hotspot_score_plus_reconstructed_leave_lncrna_out_delta",
                "minimum_reconstructed_baseline_correlation": MINIMUM_BASELINE_REPRODUCTION_R,
                "baseline_gate_status": "technical_threshold_finalized_after_source_reconstruction_audit",
                "module_discovery": False,
                "seed": SEED,
                "slurm_job_id": os.environ.get("SLURM_JOB_ID", "not_slurm"),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    artifacts = sorted(args.output.iterdir())
    pd.DataFrame(
        [
            {
                "relative_path": path.name,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in artifacts
        ]
    ).to_csv(args.output / "output_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
