"""Common-gene preparation and reference fitting for the external26 sensitivity."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from .config import canonical_json_bytes, write_json_exclusive
from .contracts import ContractError, sha256_path, verify_contract_lock
from .data import _hash_sparse, _hash_strings, verify_prepared_file
from .external_compatibility import read_label_map
from .firewall import validate_program_firewall


class ExternalReferenceError(ContractError):
    """Raised when the common-universe external sensitivity changes identity."""


POLICY_SCHEMA = "masld-cl-external26-common-universe-policy-v21"


def load_external26_policy(config: dict[str, Any], value: str | Path) -> tuple[Path, dict[str, Any]]:
    path = Path(value).resolve()
    expected = Path(config["_config_path"]).resolve().parent / "reference" / "external26_common_universe_policy_v21.json"
    if path != expected:
        raise ExternalReferenceError("external26 requires the source-controlled V21 policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("status") != "sensitivity_only_until_all_original_promotion_gates_pass"
    ):
        raise ExternalReferenceError("external26 policy identity differs")
    return path, policy


def common_gene_mapping(
    canonical_genes: list[str], external_ids: list[str], external_names: list[str],
) -> dict[str, Any]:
    """Resolve one-to-one genes, excluding every collision deterministically."""
    if len(canonical_genes) != len(set(canonical_genes)):
        raise ExternalReferenceError("canonical raw gene IDs are duplicated")
    if len(external_ids) != len(external_names):
        raise ExternalReferenceError("external gene IDs and names differ in length")
    by_id: dict[str, list[int]] = {}
    by_name: dict[str, list[int]] = {}
    for index, (gene_id, gene_name) in enumerate(zip(external_ids, external_names)):
        by_id.setdefault(str(gene_id), []).append(index)
        by_name.setdefault(str(gene_name), []).append(index)
    provisional = []
    methods = Counter()
    for canonical_index, gene in enumerate(canonical_genes):
        if len(by_id.get(gene, [])) == 1:
            provisional.append((canonical_index, by_id[gene][0], "external_var_id"))
        elif len(by_name.get(gene, [])) == 1:
            provisional.append((canonical_index, by_name[gene][0], "external_feature_name"))
    target_counts = Counter(item[1] for item in provisional)
    retained = [item for item in provisional if target_counts[item[1]] == 1]
    for _, _, method in retained:
        methods[method] += 1
    return {
        "canonical_indices": np.asarray([item[0] for item in retained], dtype=np.int64),
        "external_indices": np.asarray([item[1] for item in retained], dtype=np.int64),
        "mapping_methods": dict(sorted(methods.items())),
        "mapped_genes": len(retained),
        "colliding_external_targets_excluded": sum(value > 1 for value in target_counts.values()),
    }


def _metadata_hash(obs) -> str:
    fields = [
        "library_id", "assay_id", "donor_id", "dataset", "preparation_method",
        "native_condition_authoritative", "harmonized_stage_authoritative",
        "technical_batch", "audit_cell_type", "strict_reference",
        "external_reference_candidate", "canonical_cell", "primary_query",
        "query_control", "analysis_eligible",
    ]
    digest = hashlib.sha256()
    digest.update(_hash_strings(obs.index.astype(str)).encode("ascii"))
    for field in fields:
        if field not in obs:
            raise ExternalReferenceError(f"combined metadata field is missing: {field}")
        digest.update(field.encode())
        digest.update(_hash_strings(obs[field].astype(str)).encode("ascii"))
    return digest.hexdigest()


def prepare_external26_common(
    config: dict[str, Any], contract_lock_path: str | Path,
    canonical_prepared_path: str | Path, canonical_prepared_lock_path: str | Path,
    external_h5ad_path: str | Path, external_compatibility_path: str | Path,
    label_map_path: str | Path, policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    """Create one identical 4,000-gene matrix for strict7 and external26."""
    import anndata as ad
    import h5py
    import pandas as pd
    import scanpy as sc
    from anndata.io import sparse_dataset
    from scipy import sparse

    validate_program_firewall(config)
    policy_path, policy = load_external26_policy(config, policy_value)
    sources = policy["sources"]
    checks = {
        "canonical_prepared_sha256": canonical_prepared_path,
        "external_h5ad_sha256": external_h5ad_path,
        "external_compatibility_sha256": external_compatibility_path,
        "external_label_map_sha256": label_map_path,
    }
    for key, path in checks.items():
        if sha256_path(path) != sources[key]:
            raise ExternalReferenceError(f"external26 source hash differs: {key}")
    lock = verify_contract_lock(config, contract_lock_path, full_hash=False)
    prepared_lock = verify_prepared_file(
        canonical_prepared_path, canonical_prepared_lock_path, config, lock
    )
    with Path(external_compatibility_path).open() as handle:
        compatibility = json.load(handle)
    subset_policy = policy["external_subset"]
    if (
        compatibility.get("candidate_donors") != subset_policy["expected_new_donors"]
        or compatibility.get("candidate_cells") != subset_policy["expected_cells"]
        or compatibility.get("eligible_for_identical_4000_gene_comparison") is not False
    ):
        raise ExternalReferenceError("external compatibility census differs")

    canonical_prepared = ad.read_h5ad(canonical_prepared_path, backed="r")
    external = ad.read_h5ad(external_h5ad_path, backed="r")
    atlas_handle = h5py.File(lock["atlas_realpath"], "r")
    try:
        canonical_obs = canonical_prepared.obs.copy()
        if len(canonical_obs) != config["expected_contract"]["descriptive"]["cells"]:
            raise ExternalReferenceError("canonical prepared cell census differs")
        canonical_obs["external_reference_candidate"] = False
        canonical_obs["canonical_cell"] = True
        canonical_genes = list(map(str, atlas_handle["raw/var/_index"].asstr()[:]))
        external_ids = list(map(str, external.raw.var_names))
        external_names = list(map(str, external.raw.var["feature_name"].astype(str)))
        mapping = common_gene_mapping(canonical_genes, external_ids, external_names)
        if mapping["mapped_genes"] < policy["common_feature_rule"]["n_hvg"]:
            raise ExternalReferenceError("common one-to-one gene universe is too small")

        raw = sparse_dataset(atlas_handle["raw/X"])
        strict_positions = np.flatnonzero(canonical_obs["strict_reference"].to_numpy(dtype=bool))
        strict_counts = raw[strict_positions, :].tocsr()[:, mapping["canonical_indices"]]
        donor_values = canonical_obs.iloc[strict_positions]["donor_id"].astype(str).to_numpy()
        detected = np.zeros(strict_counts.shape[1], dtype=np.int16)
        for donor in sorted(set(donor_values)):
            detected += strict_counts[donor_values == donor].getnnz(axis=0) > 0
        eligible = detected >= policy["common_feature_rule"]["minimum_reference_donors_detected"]
        if int(eligible.sum()) < policy["common_feature_rule"]["n_hvg"]:
            raise ExternalReferenceError("fewer than 4,000 common donor-detected genes")
        eligible_positions = np.flatnonzero(eligible)
        reference = ad.AnnData(
            X=strict_counts[:, eligible_positions],
            obs=canonical_obs.iloc[strict_positions].copy(),
            var=pd.DataFrame(
                index=pd.Index(
                    np.asarray(canonical_genes)[mapping["canonical_indices"][eligible_positions]],
                    name="gene",
                )
            ),
        )
        reference.layers["counts"] = reference.X
        sc.pp.highly_variable_genes(
            reference, layer="counts", flavor="seurat_v3",
            n_top_genes=policy["common_feature_rule"]["n_hvg"],
            batch_key="donor_id", subset=False,
        )
        selected_names = reference.var_names[reference.var["highly_variable"]].astype(str).to_numpy()
        if len(selected_names) != 4000:
            raise ExternalReferenceError("common-universe HVG selection did not return 4,000 genes")
        canonical_index = {gene: i for i, gene in enumerate(canonical_genes)}
        external_index = {
            canonical_genes[c_index]: e_index
            for c_index, e_index in zip(mapping["canonical_indices"], mapping["external_indices"])
        }
        selected_canonical = np.asarray([canonical_index[gene] for gene in selected_names], dtype=np.int64)
        selected_external = np.asarray([external_index[gene] for gene in selected_names], dtype=np.int64)
        canonical_counts = raw[:, selected_canonical].tocsr()

        labels = read_label_map(label_map_path)
        included_labels = {key for key, row in labels.items() if row["include_reference"] == "True"}
        external_obs_source = external.obs
        target = (
            external_obs_source["STUDY"].astype(str).isin(subset_policy["included_studies"]).to_numpy()
            & ~external_obs_source["donor_id"].astype(str).isin(subset_policy["overlap_donor_aliases"]).to_numpy()
            & external_obs_source["author_cell_type"].astype(str).isin(included_labels).to_numpy()
        )
        external_positions = np.flatnonzero(target)
        if len(external_positions) != subset_policy["expected_cells"]:
            raise ExternalReferenceError("external retained-cell census differs during preparation")
        retained = external_obs_source.iloc[external_positions].copy()
        external_counts = external.raw.X[external_positions, :].tocsr()[:, selected_external]
        if (
            np.any(external_counts.data < 0)
            or not np.isfinite(external_counts.data).all()
            or not np.all(external_counts.data == np.floor(external_counts.data))
        ):
            raise ExternalReferenceError("external raw matrix is not finite nonnegative integer counts")
        preparation = retained["suspension_type"].astype(str).map(subset_policy["preparation_map"])
        if preparation.isna().any():
            raise ExternalReferenceError("external suspension type is not mapped")
        external_obs = pd.DataFrame(index=pd.Index(
            [f"{subset_policy['cell_prefix']}{value}" for value in retained.index.astype(str)],
            name=canonical_obs.index.name,
        ))
        external_obs["library_id"] = [f"HLiCA|{value}" for value in retained["SAMPLE"].astype(str)]
        external_obs["assay_id"] = external_obs["library_id"]
        external_obs["donor_id"] = [f"{subset_policy['donor_prefix']}{value}" for value in retained["donor_id"].astype(str)]
        external_obs["dataset"] = retained["STUDY"].astype(str).to_numpy()
        external_obs["preparation_method"] = preparation.to_numpy()
        external_obs["native_condition_authoritative"] = "Healthy_external"
        external_obs["harmonized_stage_authoritative"] = "Healthy"
        external_obs["technical_batch"] = [
            f"{study}|{prep}" for study, prep in zip(external_obs["dataset"], preparation)
        ]
        external_obs["audit_cell_type"] = [
            labels[value]["frozen_broad_label"] for value in retained["author_cell_type"].astype(str)
        ]
        external_obs["strict_reference"] = False
        external_obs["external_reference_candidate"] = True
        external_obs["canonical_cell"] = False
        external_obs["primary_query"] = False
        external_obs["query_control"] = False
        external_obs["analysis_eligible"] = False
        for field in canonical_obs.columns:
            if field not in external_obs:
                external_obs[field] = "external_not_applicable"
        external_obs = external_obs[canonical_obs.columns]
        combined_obs = pd.concat([canonical_obs, external_obs], axis=0)
        if len(combined_obs) != 1374354 or not combined_obs.index.is_unique:
            raise ExternalReferenceError("combined cell census or ID uniqueness differs")
        counts = sparse.vstack([canonical_counts, external_counts], format="csr")
        if np.any(counts.data < 0) or not np.all(counts.data == np.floor(counts.data)):
            raise ExternalReferenceError("combined count matrix is not nonnegative integer counts")
        counts.data = counts.data.astype(np.int32, copy=False)
        combined = ad.AnnData(
            X=counts, obs=combined_obs,
            var=pd.DataFrame(index=pd.Index(selected_names, name="gene")),
        )
        combined.layers["counts"] = combined.X
        combined.uns["masld_external26_common"] = {
            "schema_version": "masld-cl-external26-common-prepared-v21",
            "config_sha256": config["_config_sha256"],
            "policy_sha256": sha256_path(policy_path),
            "canonical_raw_counts_sha256": lock["raw_counts_sha256"],
            "canonical_prepared_lock_sha256": prepared_lock["lock_sha256"],
            "common_gene_order_sha256": _hash_strings(selected_names),
            "combined_obs_sha256": _metadata_hash(combined_obs),
            "combined_counts_sha256": _hash_sparse(counts),
            "canonical_cells": len(canonical_obs),
            "external_cells": len(external_obs),
            "strict7_cells": int(combined_obs["strict_reference"].sum()),
            "external_candidate_cells": int(combined_obs["external_reference_candidate"].sum()),
        }
        output = Path(output_value)
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists():
            raise FileExistsError(output)
        combined.write_h5ad(output, compression="gzip")
        result_lock = {
            **combined.uns["masld_external26_common"],
            "schema_version": "masld-cl-external26-common-lock-v21",
            "config_sha256": config["_config_sha256"],
            "policy_realpath": str(policy_path),
            "policy_sha256": sha256_path(policy_path),
            "prepared_realpath": str(output.resolve()),
            "prepared_size_bytes": output.stat().st_size,
            "prepared_file_sha256": sha256_path(output),
            "selected_canonical_indices_sha256": hashlib.sha256(selected_canonical.tobytes()).hexdigest(),
            "selected_external_indices_sha256": hashlib.sha256(selected_external.tobytes()).hexdigest(),
            "common_mapping_summary": {
                key: value for key, value in mapping.items()
                if key not in {"canonical_indices", "external_indices"}
            },
        }
        result_lock["lock_sha256"] = hashlib.sha256(canonical_json_bytes(result_lock)).hexdigest()
        lock_path = output.with_suffix(output.suffix + ".lock.json")
        write_json_exclusive(lock_path, result_lock)
        summary = {
            "schema_version": "masld-cl-external26-common-preparation-summary-v21",
            "output": str(output.resolve()),
            "lock": str(lock_path.resolve()),
            "cells": combined.n_obs,
            "genes": combined.n_vars,
            "canonical_cells": len(canonical_obs),
            "external_cells": len(external_obs),
            "reference_rosters": policy["rosters"],
            "common_mapping_summary": result_lock["common_mapping_summary"],
            "prepared_file_sha256": result_lock["prepared_file_sha256"],
        }
        write_json_exclusive(output.with_suffix(".summary.json"), summary)
        return summary
    finally:
        canonical_prepared.file.close()
        external.file.close()
        atlas_handle.close()


def verify_external26_common(
    config: dict[str, Any], prepared_value: str | Path, lock_value: str | Path,
    policy_value: str | Path,
) -> dict[str, Any]:
    path = Path(prepared_value).resolve()
    policy_path, _ = load_external26_policy(config, policy_value)
    with Path(lock_value).open() as handle:
        lock = json.load(handle)
    payload = {key: value for key, value in lock.items() if key != "lock_sha256"}
    if (
        lock.get("schema_version") != "masld-cl-external26-common-lock-v21"
        or lock.get("config_sha256") != config["_config_sha256"]
        or lock.get("policy_realpath") != str(policy_path)
        or lock.get("policy_sha256") != sha256_path(policy_path)
        or lock.get("prepared_realpath") != str(path)
        or lock.get("prepared_size_bytes") != path.stat().st_size
        or hashlib.sha256(canonical_json_bytes(payload)).hexdigest() != lock.get("lock_sha256")
        or sha256_path(path) != lock.get("prepared_file_sha256")
    ):
        raise ExternalReferenceError("combined external26 prepared lock differs")
    return lock


def external26_roster_mask(obs, roster_name: str) -> np.ndarray:
    strict = obs["strict_reference"].astype(bool).to_numpy()
    external = obs["external_reference_candidate"].astype(bool).to_numpy()
    if roster_name == "common_strict7":
        return strict
    if roster_name == "external_clean26":
        return strict | external
    raise ExternalReferenceError(f"unsupported common-universe reference roster: {roster_name}")


def _enable_inference_without_optimization(model) -> dict[str, Any]:
    """Allow inference from load_query_data without taking an optimization step."""
    before = {
        name: value.detach().cpu().clone()
        for name, value in model.module.state_dict().items()
    }
    model.is_trained = True
    after = model.module.state_dict()
    if set(before) != set(after) or any(
        not np.array_equal(value.numpy(), after[name].detach().cpu().numpy())
        for name, value in before.items()
    ):
        raise ExternalReferenceError(
            "marking the loaded mapper inference-ready changed model state"
        )
    digest = hashlib.sha256()
    for name in sorted(before):
        value = before[name]
        digest.update(name.encode("utf-8"))
        digest.update(str(value.dtype).encode("ascii"))
        digest.update(np.asarray(value.shape, dtype="<i8").tobytes())
        digest.update(value.numpy().tobytes())
    return {
        "load_query_data_weights_optimized": False,
        "is_trained_metadata_flag_set_for_inference": True,
        "state_unchanged_when_flag_set": True,
        "state_sha256": digest.hexdigest(),
    }


def train_external26_reference(
    config: dict[str, Any], prepared_value: str | Path, lock_value: str | Path,
    policy_value: str | Path, output_value: str | Path, roster_name: str, seed: int,
) -> dict[str, Any]:
    """Fit an identical common-universe all-lineage reference model."""
    import anndata as ad
    import torch
    from scvi.model import SCANVI, SCVI

    from .scvi_adapter import assert_runtime_versions
    from .training import (
        _architecture, _capped_indices, _drop_unused_categories,
        _ensure_output, _export_embedding, _external_donor_split,
        _model_bundle_identity, _set_model_labels, _train_kwargs,
        _trained_epochs, _validate_scanvi_label_codes, _write_environment,
        set_all_seeds,
    )

    assert_runtime_versions(config)
    set_all_seeds(seed)
    validate_program_firewall(config)
    common_lock = verify_external26_common(
        config, prepared_value, lock_value, policy_value
    )
    policy_path, policy = load_external26_policy(config, policy_value)
    if seed != policy["model_contract"]["seed"]:
        raise ExternalReferenceError("external26 pilot seed differs from policy")
    if roster_name not in policy["rosters"]:
        raise ExternalReferenceError("external26 roster is absent from policy")
    full = ad.read_h5ad(prepared_value)
    contract = full.uns.get("masld_external26_common", {})
    if (
        full.n_obs != 1374354
        or full.n_vars != 4000
        or contract.get("common_gene_order_sha256") != common_lock["common_gene_order_sha256"]
        or _hash_strings(full.var_names.astype(str)) != common_lock["common_gene_order_sha256"]
        or contract.get("combined_obs_sha256") != common_lock["combined_obs_sha256"]
        or contract.get("combined_counts_sha256") != common_lock["combined_counts_sha256"]
        or "counts" not in full.layers
    ):
        raise ExternalReferenceError("loaded common-universe object differs from its lock")
    counts = full.layers["counts"].data
    if (
        np.any(~np.isfinite(counts)) or np.any(counts < 0)
        or not np.all(counts == np.floor(counts))
    ):
        raise ExternalReferenceError("common-universe counts are invalid at training")
    roster = external26_roster_mask(full.obs, roster_name)
    expected = policy["rosters"][roster_name]
    if (
        int(roster.sum()) != expected["expected_cells"]
        or int(full.obs.loc[roster, "donor_id"].nunique()) != expected["expected_donors"]
    ):
        raise ExternalReferenceError("external26 training roster census differs")
    full.obs["strict_reference"] = roster
    training_indices = _capped_indices(
        full, np.flatnonzero(roster), "all_lineage", config, seed
    )
    reference = full[training_indices].copy()
    batch_key = config["features"]["selection_batch_key"]
    _drop_unused_categories(reference, (batch_key,))
    _set_model_labels(
        reference, query_unknown=False,
        unlabeled=config["features"]["unlabeled_category"],
    )
    _validate_scanvi_label_codes(
        reference, unlabeled=config["features"]["unlabeled_category"]
    )
    external_split = _external_donor_split(reference, seed)
    SCVI.setup_anndata(
        reference, layer="counts", batch_key=batch_key,
        labels_key=config["features"]["labels_key"],
    )
    vae = SCVI(reference, **_architecture(config))
    vae.train(**_train_kwargs(config, "all_lineage", external_split))
    model = SCANVI.from_scvi_model(
        vae, unlabeled_category=config["features"]["unlabeled_category"]
    )
    model.train(**_train_kwargs(config, "all_lineage", external_split))
    output = _ensure_output(output_value)
    model.save(output / "model", save_anndata=True)
    model_bundle = _model_bundle_identity(output / "model")
    np.save(output / "reference_latent.npy", model.get_latent_representation(reference))

    _set_model_labels(
        full, query_unknown=True,
        unlabeled=config["features"]["unlabeled_category"],
    )
    mapping_model = SCANVI.load_query_data(
        full, model, unfrozen=False,
        accelerator="gpu" if torch.cuda.is_available() else "cpu", device=1,
    )
    mapping_inference = _enable_inference_without_optimization(mapping_model)
    mapping_model.save(output / "mapping_model", save_anndata=False)
    mapping_bundle = _model_bundle_identity(output / "mapping_model")
    embedding = _export_embedding(
        mapping_model, full, np.ones(full.n_obs, dtype=bool),
        "all_lineage", config, output,
    )
    result = {
        "schema_version": "masld-cl-external26-common-reference-v21",
        "model_kind": "all_lineage",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "common_prepared_lock": {
            "path": str(Path(lock_value).resolve()), "sha256": sha256_path(lock_value),
            "lock_sha256": common_lock["lock_sha256"],
        },
        "roster_name": roster_name,
        "seed": seed,
        "n_reference_cells_uncapped": int(roster.sum()),
        "n_reference_donors": int(full.obs.loc[roster, "donor_id"].nunique()),
        "n_training_cells": int(reference.n_obs),
        "reference_donors": sorted(full.obs.loc[roster, "donor_id"].astype(str).unique()),
        "genes": list(map(str, reference.var_names)),
        "stopping_epoch": _trained_epochs(model),
        "model_bundle": model_bundle,
        "mapping_model_bundle": mapping_bundle,
        "mapping_inference": mapping_inference,
        "embedding": embedding,
        "query_labels_hidden": True,
        "case_stage_program_hero_gene_umap_cas13_used": False,
        "sensitivity_only": True,
    }
    write_json_exclusive(output / "reference_manifest.json", result)
    _write_environment(output)
    return result


def recover_external26_reference_export(
    config: dict[str, Any], prepared_value: str | Path, lock_value: str | Path,
    policy_value: str | Path, source_value: str | Path, output_value: str | Path,
    roster_name: str, seed: int,
) -> dict[str, Any]:
    """Recover projection from a completed fit that failed before embedding export."""
    import os
    import shutil

    import anndata as ad
    import torch
    from scvi.model import SCANVI

    from .scvi_adapter import assert_runtime_versions
    from .training import (
        _ensure_output, _export_embedding, _model_bundle_identity,
        _set_model_labels, _trained_epochs, _write_environment, set_all_seeds,
    )

    assert_runtime_versions(config)
    set_all_seeds(seed)
    validate_program_firewall(config)
    common_lock = verify_external26_common(
        config, prepared_value, lock_value, policy_value
    )
    policy_path, policy = load_external26_policy(config, policy_value)
    if seed != policy["model_contract"]["seed"] or roster_name not in policy["rosters"]:
        raise ExternalReferenceError("recovery roster or seed differs from policy")
    source = Path(source_value).resolve()
    if (source / "reference_manifest.json").exists():
        raise ExternalReferenceError("recovery source is already a completed reference")
    for relative in ("model/model.pt", "model/adata.h5ad", "mapping_model/model.pt"):
        if not (source / relative).is_file():
            raise ExternalReferenceError(f"recovery source artifact is missing: {relative}")

    full = ad.read_h5ad(prepared_value)
    contract = full.uns.get("masld_external26_common", {})
    if (
        full.n_obs != 1374354 or full.n_vars != 4000
        or contract.get("common_gene_order_sha256") != common_lock["common_gene_order_sha256"]
        or contract.get("combined_obs_sha256") != common_lock["combined_obs_sha256"]
        or contract.get("combined_counts_sha256") != common_lock["combined_counts_sha256"]
    ):
        raise ExternalReferenceError("recovery input differs from common-universe lock")
    roster = external26_roster_mask(full.obs, roster_name)
    expected = policy["rosters"][roster_name]
    if (
        int(roster.sum()) != expected["expected_cells"]
        or int(full.obs.loc[roster, "donor_id"].nunique()) != expected["expected_donors"]
    ):
        raise ExternalReferenceError("recovery reference roster differs")
    full.obs["strict_reference"] = roster
    _set_model_labels(
        full, query_unknown=True,
        unlabeled=config["features"]["unlabeled_category"],
    )
    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    fitted_model = SCANVI.load(source / "model", accelerator=accelerator, device=1)
    mapping_model = SCANVI.load(
        source / "mapping_model", adata=full,
        accelerator=accelerator, device=1,
    )
    mapping_inference = _enable_inference_without_optimization(mapping_model)
    output = _ensure_output(output_value)
    shutil.copytree(source / "model", output / "model", copy_function=os.link)
    mapping_model.save(output / "mapping_model", save_anndata=False)
    embedding = _export_embedding(
        mapping_model, full, np.ones(full.n_obs, dtype=bool),
        "all_lineage", config, output,
    )
    result = {
        "schema_version": "masld-cl-external26-common-reference-v21",
        "model_kind": "all_lineage",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "common_prepared_lock": {
            "path": str(Path(lock_value).resolve()), "sha256": sha256_path(lock_value),
            "lock_sha256": common_lock["lock_sha256"],
        },
        "roster_name": roster_name,
        "seed": seed,
        "n_reference_cells_uncapped": int(roster.sum()),
        "n_reference_donors": int(full.obs.loc[roster, "donor_id"].nunique()),
        "n_training_cells": int(fitted_model.adata.n_obs),
        "reference_donors": sorted(full.obs.loc[roster, "donor_id"].astype(str).unique()),
        "genes": list(map(str, full.var_names)),
        "stopping_epoch": _trained_epochs(fitted_model),
        "model_bundle": _model_bundle_identity(output / "model"),
        "mapping_model_bundle": _model_bundle_identity(output / "mapping_model"),
        "mapping_inference": mapping_inference,
        "embedding": embedding,
        "query_labels_hidden": True,
        "case_stage_program_hero_gene_umap_cas13_used": False,
        "sensitivity_only": True,
        "recovered_after_projection_only_failure": True,
        "recovery_source": {
            "path": str(source),
            "model_bundle": _model_bundle_identity(source / "model"),
            "mapping_model_bundle": _model_bundle_identity(source / "mapping_model"),
        },
    }
    write_json_exclusive(output / "reference_manifest.json", result)
    _write_environment(output)
    return result
