#!/usr/bin/env python3
"""Re-score all 117 frozen programs on decontX-corrected counts."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import hotspot
import numpy as np
import pandas as pd
import scipy.sparse as sp
from hotspot import modules as hotspot_modules

from common import (
    HOTSPOT_ROOT,
    LINEAGES,
    LINEAGE_TO_LABEL,
    DONOR_METADATA,
    MEMBERSHIP,
    PAIRING_FILES,
    PROGRAM_ROOT,
    READY,
    REGISTRY,
    SCORE_SCRIPT,
    SEED,
    bh_adjust,
    build_sample_to_donor,
    evidence_state,
    load_donor_metadata,
    refuse_existing,
    require,
    sha256,
    stage_fit,
)


CANDIDATE = Path(os.environ["CAND_ROOT"])
WORK = CANDIDATE / "work"
RESULTS = CANDIDATE / "results"
DECONTX = RESULTS / "decontx"
MINIMUM_REPRODUCTION_R = 0.995


def load_score_producer():
    specification = importlib.util.spec_from_file_location("hotspot_501", SCORE_SCRIPT)
    require(specification is not None and specification.loader is not None, "cannot load 501")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def read_matrix(dataset: str, matrix_name: str) -> tuple[pd.DataFrame, sp.csc_matrix]:
    directory = WORK / dataset
    dimensions = json.loads((directory / "dec_dims.json").read_text())
    cells = pd.read_csv(directory / "dec_cells.csv.gz", dtype=str)
    n_cells = int(dimensions["n_cells"])
    n_genes = int(dimensions["n_genes_score"])
    require(len(cells) == n_cells, f"corrected cell index drift: {dataset}")
    prefix = "dec" if matrix_name == "corrected" else "raw"
    nnz_key = "corrected_nnz" if matrix_name == "corrected" else "raw_nnz"
    sum_key = "corrected_sum" if matrix_name == "corrected" else "raw_sum"
    nnz = int(dimensions[nnz_key])
    indices = np.fromfile(directory / f"{prefix}_indices.bin", dtype="<i4")
    indptr = np.fromfile(directory / f"{prefix}_indptr.bin", dtype="<i4")
    values = np.fromfile(directory / f"{prefix}_data.bin", dtype="<f8")
    require(len(indices) == nnz and len(values) == nnz, f"binary nnz drift: {dataset}")
    require(len(indptr) == n_cells + 1, f"binary pointer drift: {dataset}")
    matrix = sp.csc_matrix((values, indices, indptr), shape=(n_genes, n_cells))
    require(
        np.isclose(matrix.data.sum(), float(dimensions[sum_key]), rtol=1e-10),
        f"binary sum drift: {dataset} {matrix_name}",
    )
    return cells, matrix


def assemble_lineage_counts(
    lineage: str, cell_index: pd.Index, score_genes: list[str]
) -> tuple[sp.csr_matrix, sp.csr_matrix, pd.DataFrame]:
    raw_blocks: list[sp.csc_matrix] = []
    corrected_blocks: list[sp.csc_matrix] = []
    cell_frames: list[pd.DataFrame] = []
    manifest = json.loads((WORK / "manifest.json").read_text())
    for dataset in sorted(manifest["datasets"]):
        cells, raw = read_matrix(dataset, "raw")
        _, corrected = read_matrix(dataset, "corrected")
        keep = cells["program_lineage"].eq(lineage).to_numpy()
        if not keep.any():
            continue
        raw_blocks.append(raw[:, keep])
        corrected_blocks.append(corrected[:, keep])
        cell_frames.append(cells.loc[keep].reset_index(drop=True))
    require(raw_blocks, f"no corrected-count blocks for {lineage}")
    all_cells = pd.concat(cell_frames, ignore_index=True)
    require(not all_cells["cell_id"].duplicated().any(), f"duplicate cells for {lineage}")
    position = cell_index.get_indexer(all_cells["cell_id"].astype(str))
    require((position >= 0).all(), f"corrected cells outside frozen universe: {lineage}")
    require(len(position) == len(cell_index), f"incomplete corrected universe: {lineage}")
    inverse = np.empty(len(cell_index), dtype=np.int64)
    inverse[position] = np.arange(len(position))
    raw_all = sp.hstack(raw_blocks, format="csc")[:, inverse].tocsr()
    corrected_all = sp.hstack(corrected_blocks, format="csc")[:, inverse].tocsr()
    all_cells = all_cells.iloc[inverse].reset_index(drop=True)
    require(
        all_cells["cell_id"].astype(str).tolist() == cell_index.astype(str).tolist(),
        f"lineage cell order drift: {lineage}",
    )
    require(raw_all.shape == (len(score_genes), len(cell_index)), "raw assembly shape drift")
    return raw_all, corrected_all, all_cells


def fit_program_universe(
    donor_scores: pd.DataFrame,
    donor_metadata: pd.DataFrame,
    failed_datasets: set[str],
    universe: str,
) -> tuple[dict[str, object], pd.DataFrame]:
    frame = donor_scores.merge(donor_metadata, on="donor", how="inner", validate="one_to_one")
    require(len(frame) == len(donor_scores), "donor score does not join current metadata")
    if universe == "complete_case_common_universe":
        frame = frame.loc[~frame["dataset"].isin(failed_datasets)].copy()
    elif universe != "full_universe_passthrough_sensitivity":
        raise RuntimeError(f"unknown analysis universe: {universe}")
    frame["delta"] = frame["corrected"] - frame["raw"]
    results: dict[str, object] = {}
    for label, column in (("raw", "raw"), ("corrected", "corrected"), ("delta", "delta")):
        try:
            fit = stage_fit(frame, column)
            for key, value in fit.items():
                results[f"{label}_{key}"] = value
            results[f"{label}_failure_reason"] = ""
        except RuntimeError as error:
            for key in (
                "beta", "se", "pvalue", "hc3_se", "hc3_pvalue", "ci_low", "ci_high",
                "hc3_ci_low", "hc3_ci_high", "n", "residual_df", "max_leverage",
                "n_donors", "n_datasets", "n_stage_levels",
            ):
                results[f"{label}_{key}"] = np.nan
            results[f"{label}_failure_reason"] = str(error)
    return results, frame


def build_ambient_burden(
    registry: pd.DataFrame, membership: pd.DataFrame, datasets: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    ambient_tables: dict[str, pd.DataFrame] = {}
    run_status: dict[str, tuple[str, str]] = {}
    for dataset in datasets:
        run = pd.read_csv(DECONTX / f"{dataset}__run_log.tsv", sep="\t", dtype=str)
        require(len(run) == 1, f"run-log cardinality drift: {dataset}")
        run_status[dataset] = (run.loc[0, "correction_status"], run.loc[0, "failure_reason"])
        ambient_tables[dataset] = pd.read_csv(
            DECONTX / f"{dataset}__ambient_by_gene_lineage.tsv.gz", sep="\t"
        )

    rows: list[dict[str, object]] = []
    for program in registry.itertuples(index=False):
        members = membership.loc[
            membership["program_uid"] == program.program_uid,
            ["source_gene", "original_l1_weight"],
        ].copy()
        require(len(members) == int(program.n_source_genes), "program membership size drift")
        members["original_l1_weight"] = members["original_l1_weight"].astype(float)
        require(
            np.isclose(members["original_l1_weight"].sum(), 1.0, atol=1e-10),
            f"program weights do not sum to one: {program.program_uid}",
        )
        native_label = LINEAGE_TO_LABEL[program.cell_type]
        for dataset in datasets:
            status, reason = run_status[dataset]
            data = ambient_tables[dataset]
            native = data.loc[data["cell_type"] == native_label, ["gene", "ambient_fraction"]]
            joined = members.merge(
                native, left_on="source_gene", right_on="gene", how="left", validate="one_to_one"
            )
            finite = joined["ambient_fraction"].notna()
            weight_covered = float(joined.loc[finite, "original_l1_weight"].sum())
            burden = np.nan
            burden_status = status
            if status == "corrected" and weight_covered >= 0.8:
                burden = float(
                    np.average(
                        joined.loc[finite, "ambient_fraction"].astype(float),
                        weights=joined.loc[finite, "original_l1_weight"].astype(float),
                    )
                )
                burden_status = "testable"
            elif status == "corrected" and native.empty:
                burden_status = "no_native_lineage_cells"
            elif status == "corrected":
                burden_status = "insufficient_weight_coverage"
            rows.append(
                {
                    "program_uid": program.program_uid,
                    "cell_type": program.cell_type,
                    "module": int(program.module),
                    "module_name": program.module_name,
                    "dataset": dataset,
                    "native_cell_label": native_label,
                    "correction_status": status,
                    "failure_reason": reason,
                    "burden_status": burden_status,
                    "n_member_genes": int(len(members)),
                    "n_members_with_ambient_estimate": int(finite.sum()),
                    "retained_l1_weight": weight_covered,
                    "l1_weighted_ambient_fraction": burden,
                }
            )
    per_dataset = pd.DataFrame(rows)
    summaries = []
    for program_uid, group in per_dataset.groupby("program_uid", sort=False):
        values = group.loc[
            group["burden_status"] == "testable", "l1_weighted_ambient_fraction"
        ].astype(float)
        summaries.append(
            {
                "program_uid": program_uid,
                "n_testable_datasets": int(len(values)),
                "ambient_burden_median": float(values.median()) if len(values) else np.nan,
                "ambient_burden_min": float(values.min()) if len(values) else np.nan,
                "ambient_burden_max": float(values.max()) if len(values) else np.nan,
            }
        )
    return per_dataset, pd.DataFrame(summaries)


def build_dataset_lineage_qc(
    datasets: list[str], sample_to_donor: dict[str, str]
) -> pd.DataFrame:
    rows = []
    for dataset in datasets:
        run = pd.read_csv(DECONTX / f"{dataset}__run_log.tsv", sep="\t", dtype=str).iloc[0]
        cells = pd.read_csv(WORK / dataset / "dec_cells.csv.gz", dtype=str)
        cells["donor"] = cells["sample"].map(sample_to_donor).fillna(cells["sample"])
        contamination = pd.read_csv(
            DECONTX / f"{dataset}__contamination_per_program_cell.tsv.gz",
            sep="\t",
            dtype={"cell_id": str, "program_lineage": str},
        )
        require(
            len(contamination) == len(cells)
            and contamination["cell_id"].astype(str).tolist()
            == cells["cell_id"].astype(str).tolist(),
            f"contamination/cell index drift: {dataset}",
        )
        contamination["contamination"] = pd.to_numeric(
            contamination["contamination"], errors="coerce"
        )
        for lineage in LINEAGES:
            subset = cells[cells["program_lineage"] == lineage]
            contamination_values = contamination.loc[
                contamination["program_lineage"] == lineage, "contamination"
            ].dropna()
            rows.append(
                {
                    "dataset": dataset,
                    "program_lineage": lineage,
                    "correction_status": run["correction_status"],
                    "failure_reason": run["failure_reason"],
                    "failure_detail": (
                        "decontX aborted because an ambient batch contained only one "
                        "contamination cluster; counts were passed through unchanged"
                        if run["failure_reason"] == "decontX_error" else ""
                    ),
                    "n_cells": int(len(subset)),
                    "n_sequencing_samples": int(subset["sample"].nunique()),
                    "n_biological_donors": int(subset["donor"].nunique()),
                    "n_ambient_batches": int(run["n_ambient_batches"]),
                    "n_pooled_small_samples": int(run["n_pooled_small_samples"]),
                    "n_cells_with_contamination_estimate": int(len(contamination_values)),
                    "contamination_median": (
                        float(contamination_values.median())
                        if len(contamination_values) else np.nan
                    ),
                    "contamination_mean": (
                        float(contamination_values.mean())
                        if len(contamination_values) else np.nan
                    ),
                    "contamination_q10": (
                        float(contamination_values.quantile(0.10))
                        if len(contamination_values) else np.nan
                    ),
                    "contamination_q90": (
                        float(contamination_values.quantile(0.90))
                        if len(contamination_values) else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    outputs = {
        "effects": RESULTS / "ambient_program_effects.tsv",
        "donor_scores": RESULTS / "donor_program_scores.tsv.gz",
        "reproduction": RESULTS / "raw_score_reproduction.tsv",
        "burden": RESULTS / "program_dataset_ambient_burden.tsv",
        "burden_summary": RESULTS / "program_ambient_burden_summary.tsv",
        "dataset_qc": RESULTS / "ambient_dataset_lineage_qc.tsv",
        "donor_roster": RESULTS / "current_source_diagnosis_donor_roster.tsv",
        "provenance": RESULTS / "rescore_provenance.json",
    }
    for path in outputs.values():
        refuse_existing(path)

    ready = pd.read_csv(READY, sep="\t", dtype=str)
    require(len(ready) == 1, "READY seal malformed")
    require(sha256(REGISTRY) == ready.loc[0, "registry_sha256"], "registry hash drift")
    require(
        sha256(MEMBERSHIP) == ready.loc[0, "membership_table_sha256"],
        "membership hash drift",
    )
    registry = pd.read_csv(REGISTRY, sep="\t")
    membership = pd.read_csv(MEMBERSHIP, sep="\t")
    require(len(registry) == 117, "expected 117 frozen programs")
    score_genes = (WORK / "score_genes.txt").read_text().splitlines()
    gene_position = {gene: index for index, gene in enumerate(score_genes)}
    sample_to_donor = build_sample_to_donor()
    donor_metadata = load_donor_metadata()
    manifest = json.loads((WORK / "manifest.json").read_text())
    datasets = sorted(manifest["datasets"])

    run_logs = pd.concat(
        [pd.read_csv(DECONTX / f"{dataset}__run_log.tsv", sep="\t") for dataset in datasets],
        ignore_index=True,
    )
    failed_datasets = set(
        run_logs.loc[run_logs["correction_status"] != "corrected", "dataset"].astype(str)
    )
    require("GSE189600" in set(run_logs["dataset"]), "GSE189600 is absent from QC")

    producer = load_score_producer()
    reproduction_rows: list[dict[str, object]] = []
    donor_score_rows: list[pd.DataFrame] = []
    effect_rows: list[dict[str, object]] = []
    raw_transport_checks: dict[str, dict[str, float]] = {}
    np.random.seed(SEED)

    for lineage in LINEAGES:
        run_metadata = json.loads((HOTSPOT_ROOT / lineage / "run_metadata.json").read_text())
        atlas = producer.load_atlas(lineage, smoke=False)
        excluded = set(run_metadata.get("exclude_datasets") or [])
        if excluded:
            atlas = atlas[~atlas.obs["dataset"].astype(str).isin(excluded)].copy()
        atlas = producer.strip_confounders(atlas)
        atlas = producer.filter_detected(atlas, min_frac=0.01)
        require(atlas.n_obs == int(run_metadata["n_cells"]), f"cell census drift: {lineage}")
        require(atlas.n_vars == int(run_metadata["n_genes_kept"]), f"gene census drift: {lineage}")
        producer.ensure_raw_layer(atlas)
        latent = producer.resolve_latent(atlas)
        neighbors = int(run_metadata["params"]["n_neighbors"])
        hotspot_object = hotspot.Hotspot(
            atlas,
            layer_key="counts",
            model="danb",
            latent_obsm_key=latent,
            umi_counts_obs_key="n_counts",
        )
        hotspot_object.create_knn_graph(weighted_graph=False, n_neighbors=neighbors)
        cell_index = pd.Index(atlas.obs_names.astype(str))
        raw_counts, corrected_counts, cell_table = assemble_lineage_counts(
            lineage, cell_index, score_genes
        )
        lineage_registry = registry[registry["cell_type"] == lineage].set_index("module")
        lineage_membership = membership[membership["cell_type"] == lineage]
        native_layer = atlas.layers["counts"]
        native_data = native_layer.data if sp.issparse(native_layer) else np.asarray(native_layer)
        native_substrate_all_integer = bool(
            np.all(np.asarray(native_data) == np.round(np.asarray(native_data)))
        )
        raw_transport_checks[lineage] = {}
        native_genes = sorted(set(lineage_membership["source_gene"].astype(str)))
        priority_genes = ["IGFBP7", "BICC1", "ALB", "DOCK2", "STAB2"]
        check_genes = [gene for gene in priority_genes if gene in native_genes]
        check_genes += [gene for gene in native_genes if gene not in check_genes][:10]
        require(check_genes, f"no native raw-transport sentinels: {lineage}")
        for gene in check_genes:
            atlas_values = atlas[:, gene].layers["counts"]
            if sp.issparse(atlas_values):
                atlas_values = atlas_values.toarray()
            atlas_values = np.asarray(atlas_values).ravel()
            transported = np.asarray(raw_counts[gene_position[gene], :].todense()).ravel()
            difference = float(np.max(np.abs(atlas_values - transported)))
            raw_transport_checks[lineage][gene] = difference
            if native_substrate_all_integer:
                require(difference < 1e-6, f"raw transport mismatch: {lineage} {gene}")

        stored = pd.read_parquet(HOTSPOT_ROOT / lineage / "cell_scores.parquet")
        stored["module"] = stored["module"].astype(int)
        cell_donors = (
            atlas.obs["sample"].astype(str).map(sample_to_donor).fillna(atlas.obs["sample"].astype(str))
        ).to_numpy()
        cell_datasets = atlas.obs["dataset"].astype(str).to_numpy()

        for module in sorted(int(value) for value in lineage_registry.index):
            program = lineage_registry.loc[module]
            genes = lineage_membership.loc[
                lineage_membership["module"].astype(int) == module, "source_gene"
            ].astype(str).tolist()
            require(len(genes) == int(program["n_source_genes"]), "membership count drift")
            require(set(genes) <= set(atlas.var_names.astype(str)), f"atlas gene missing: {lineage} {module}")
            require(set(genes) <= set(gene_position), f"transport gene missing: {lineage} {module}")
            gene_indices = np.array([gene_position[gene] for gene in genes])
            raw_dense = np.asarray(raw_counts[gene_indices, :].todense(), dtype=np.float64)
            corrected_dense = np.asarray(
                corrected_counts[gene_indices, :].todense(), dtype=np.float64
            )
            if native_substrate_all_integer:
                native_dense = raw_dense
            else:
                native_dense = atlas[:, genes].layers["counts"]
                if sp.issparse(native_dense):
                    native_dense = native_dense.toarray()
                native_dense = np.asarray(native_dense, dtype=np.float64).T
            arguments = (
                hotspot_object.model,
                hotspot_object.umi_counts.values,
                hotspot_object.neighbors.values,
                hotspot_object.weights.values,
            )
            reconstructed_native = hotspot_modules.compute_scores(native_dense, *arguments)
            reconstructed_raw = hotspot_modules.compute_scores(raw_dense, *arguments)
            reconstructed_corrected = hotspot_modules.compute_scores(corrected_dense, *arguments)
            stored_vector = (
                stored[stored["module"] == module]
                .set_index("cell_id")["score"]
                .reindex(cell_index)
            )
            require(not stored_vector.isna().any(), f"stored score join failed: {lineage} {module}")
            stored_vector = stored_vector.to_numpy(float)
            reproduction_r = float(np.corrcoef(reconstructed_native, stored_vector)[0, 1])
            max_error = float(np.max(np.abs(reconstructed_native - stored_vector)))
            require(
                reproduction_r >= MINIMUM_REPRODUCTION_R,
                f"raw score reproduction failed: {lineage} {module} r={reproduction_r}",
            )
            transport_reproduction_r = float(
                np.corrcoef(reconstructed_raw, stored_vector)[0, 1]
            )
            transport_gate_pass = bool(
                transport_reproduction_r >= MINIMUM_REPRODUCTION_R
            )
            corrected_vector = (
                stored_vector + (reconstructed_corrected - reconstructed_raw)
                if transport_gate_pass
                else np.full_like(stored_vector, np.nan)
            )
            reproduction_rows.append(
                {
                    "program_uid": program["program_uid"],
                    "cell_type": lineage,
                    "module": module,
                    "module_name": program["module_name"],
                    "n_cells": int(len(cell_index)),
                    "n_member_genes": int(len(genes)),
                    "n_neighbors": neighbors,
                    "native_substrate_all_integer": native_substrate_all_integer,
                    "reproduction_r": reproduction_r,
                    "max_abs_error": max_error,
                    "transport_reproduction_r": transport_reproduction_r,
                    "transport_gate_pass": transport_gate_pass,
                    "transport_failure_reason": (
                        "complete-atlas integer counts do not reproduce the frozen native score"
                        if not transport_gate_pass else ""
                    ),
                }
            )
            donor = pd.DataFrame(
                {
                    "donor": cell_donors,
                    "dataset_from_cells": cell_datasets,
                    "raw": stored_vector,
                    "corrected": corrected_vector,
                }
            )
            donor = donor.groupby("donor", as_index=False).agg(
                dataset_from_cells=("dataset_from_cells", "first"),
                raw=("raw", "mean"),
                corrected=("corrected", "mean"),
            )
            donor.insert(0, "module", module)
            donor.insert(0, "cell_type", lineage)
            donor.insert(0, "program_uid", program["program_uid"])
            donor_score_rows.append(donor)
            for universe in (
                "complete_case_common_universe",
                "full_universe_passthrough_sensitivity",
            ):
                fits, model_frame = fit_program_universe(
                    donor[["donor", "raw", "corrected"]],
                    donor_metadata,
                    failed_datasets,
                    universe,
                )
                if universe == "full_universe_passthrough_sensitivity":
                    require(
                        np.isclose(float(fits["raw_beta"]), float(program["primary_beta"]), atol=1e-8, rtol=1e-8),
                        f"current registry beta not reproduced: {lineage} {module}",
                    )
                    require(
                        np.isclose(float(fits["raw_se"]), float(program["primary_se"]), atol=1e-8, rtol=1e-8),
                        f"current registry SE not reproduced: {lineage} {module}",
                    )
                    require(
                        np.isclose(float(fits["raw_hc3_se"]), float(program["primary_hc3_se"]), atol=1e-8, rtol=1e-8),
                        f"current registry HC3 SE not reproduced: {lineage} {module}",
                    )
                effect_rows.append(
                    {
                        "program_uid": program["program_uid"],
                        "cell_type": lineage,
                        "module": module,
                        "module_name": program["module_name"],
                        "analysis_universe": universe,
                        "failed_decontx_datasets": ";".join(sorted(failed_datasets)),
                        "primary_selected_frozen": bool(program["primary_selected"]),
                        "robust_display_frozen": bool(program["robust_display"]),
                        "source_stability_frozen": program["source_stability"],
                        "score_reproduction_gate_pass": transport_gate_pass,
                        "score_reproduction_failure_reason": (
                            "complete-atlas integer counts do not reproduce the frozen native score"
                            if not transport_gate_pass else ""
                        ),
                        "n_model_donors": int(
                            model_frame.loc[
                                (~model_frame["exclude"]) & model_frame["stage_ordinal"].notna(),
                                "donor",
                            ].nunique()
                        ),
                        **fits,
                    }
                )
            del raw_dense, corrected_dense, native_dense
        del atlas, hotspot_object, raw_counts, corrected_counts

    effects = pd.DataFrame(effect_rows)
    require(len(effects) == 234, "effect table must contain 117 programs x 2 universes")
    for universe, index in effects.groupby("analysis_universe").groups.items():
        index = list(index)
        require(len(index) == 117, f"incomplete family: {universe}")
        for prefix in ("raw", "corrected", "delta"):
            effects.loc[index, f"{prefix}_qvalue"] = bh_adjust(
                effects.loc[index, f"{prefix}_pvalue"], family_size=117
            )
            effects.loc[index, f"{prefix}_hc3_qvalue"] = bh_adjust(
                effects.loc[index, f"{prefix}_hc3_pvalue"], family_size=117
            )
    effects["evidence_state"] = "sensitivity_only"
    primary = effects["analysis_universe"] == "complete_case_common_universe"
    effects.loc[primary, "evidence_state"] = effects.loc[primary].apply(evidence_state, axis=1)

    control_role = {
        "hotspot_hepatocytes_b1b1287a12c11c43": "positive_control_leukocyte_DOCK2",
        "hotspot_hepatocytes_bece6cd6d8fc9c67": "positive_control_endothelial_STAB2",
        "hotspot_hepatocytes_308bfb88c3a6d001": "intrinsic_control_secretory_ALB",
        "hotspot_hepatocytes_b090a35eafb51328": "intrinsic_control_xenobiotic_CYP",
        "hotspot_hepatocytes_f05c535ae5bbc0b9": "hero_ECM_IGFBP7",
        "hotspot_hepatocytes_48f39dd4d817a10e": "hero_ductular_BICC1",
    }
    effects["display_role"] = effects["program_uid"].map(control_role).fillna("")

    burden, burden_summary = build_ambient_burden(registry, membership, datasets)
    effects = effects.merge(burden_summary, on="program_uid", how="left", validate="many_to_one")
    dataset_qc = build_dataset_lineage_qc(datasets, sample_to_donor)

    effects.to_csv(outputs["effects"], sep="\t", index=False)
    pd.concat(donor_score_rows, ignore_index=True).to_csv(
        outputs["donor_scores"], sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame(reproduction_rows).to_csv(outputs["reproduction"], sep="\t", index=False)
    burden.to_csv(outputs["burden"], sep="\t", index=False)
    burden_summary.to_csv(outputs["burden_summary"], sep="\t", index=False)
    dataset_qc.to_csv(outputs["dataset_qc"], sep="\t", index=False)
    donor_metadata[
        [
            "donor", "dataset", "disease_stage_coarse", "exclude_stage_analysis",
            "exclude", "stage_ordinal",
        ]
    ].sort_values("donor").to_csv(outputs["donor_roster"], sep="\t", index=False)
    reproduction_table = pd.DataFrame(reproduction_rows)
    provenance = {
        "seed": SEED,
        "registry_sha256": sha256(REGISTRY),
        "membership_sha256": sha256(MEMBERSHIP),
        "donor_metadata_path": str(DONOR_METADATA.resolve()),
        "donor_metadata_sha256": sha256(DONOR_METADATA),
        "pairing_file_sha256": {
            dataset: sha256(path) for dataset, path in sorted(PAIRING_FILES.items())
        },
        "n_programs": 117,
        "n_effect_rows": int(len(effects)),
        "failed_decontx_datasets": sorted(failed_datasets),
        "minimum_raw_score_reproduction_r": MINIMUM_REPRODUCTION_R,
        "observed_minimum_raw_score_reproduction_r": float(
            reproduction_table["reproduction_r"].min()
        ),
        "n_transport_testable_programs": int(
            reproduction_table["transport_gate_pass"].astype(bool).sum()
        ),
        "n_transport_untestable_programs": int(
            (~reproduction_table["transport_gate_pass"].astype(bool)).sum()
        ),
        "transport_untestable_program_uids": reproduction_table.loc[
            ~reproduction_table["transport_gate_pass"].astype(bool), "program_uid"
        ].tolist(),
        "raw_transport_checks": raw_transport_checks,
        "model": "program_score ~ stage_ordinal + factor(dataset)",
        "delta_model": "(corrected_score - raw_score) ~ stage_ordinal + factor(dataset)",
        "bh_family_size": 117,
    }
    outputs["provenance"].write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
