#!/usr/bin/env python3
"""Infer a fixed 50:50 adapted/reference scANVI head consensus."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import torch
from scvi.model import SCANVI

from masld_cl.config import load_config, write_json_exclusive
from masld_cl.contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from masld_cl.embedding import load_embedding
from masld_cl.training import _set_model_labels, set_all_seeds


def _digest(state: dict[str, torch.Tensor], keys: list[str]) -> str:
    result = hashlib.sha256()
    for key in keys:
        value = state[key].detach().cpu().contiguous()
        result.update(key.encode())
        result.update(str(value.dtype).encode())
        result.update(np.asarray(value.shape, dtype=np.int64).tobytes())
        result.update(value.numpy().tobytes())
    return result.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--reference-model", required=True)
    parser.add_argument("--adapted-model", required=True)
    parser.add_argument("--v14-embedding", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    policy_path = Path(args.policy).resolve()
    with policy_path.open() as handle:
        policy = json.load(handle)
    rule = policy.get("consensus_rule", {})
    if (
        policy.get("schema_version") != "masld-cl-dual-head-consensus-policy-v16"
        or policy.get("config_sha256") != config["_config_sha256"]
        or rule.get("adapted_head_probability_weight") != 0.5
        or rule.get("reference_head_probability_weight") != 0.5
        or rule.get("query_labels_used") is not False
        or rule.get("grid_or_tuning") is not False
    ):
        raise ContractError("V16 policy identity or fixed rule differs")

    v14_path = Path(args.v14_embedding).resolve()
    v14_info, _, v14_cells = load_embedding(v14_path)
    if v14_info.get("method") != "geometry_preserving_continual_adapter":
        raise ContractError("V16 V14 source differs")
    set_all_seeds(17)
    reference = SCANVI.load(Path(args.reference_model).resolve())
    adapted = SCANVI.load(Path(args.adapted_model).resolve())
    reference_state = reference.module.state_dict()
    adapted_state = adapted.module.state_dict()
    classifier = sorted(key for key in adapted_state if key.startswith("classifier."))
    reference_classifier = sorted(
        key for key in reference_state if key.startswith("classifier.")
    )
    if (
        not classifier
        or classifier != reference_classifier
        or any(adapted_state[key].shape != reference_state[key].shape for key in classifier)
    ):
        raise ContractError("V16 classifier state rosters differ")
    nonclassifier = sorted(set(adapted_state) - set(classifier))
    nonclassifier_before = _digest(adapted_state, nonclassifier)
    adapted_classifier_digest = _digest(adapted_state, classifier)
    reference_classifier_digest = _digest(reference_state, classifier)

    full = ad.read_h5ad(Path(args.prepared).resolve())
    cell_ids = v14_cells["cell_id"].astype(str).tolist()
    if len(cell_ids) != len(set(cell_ids)) or not set(cell_ids).issubset(set(full.obs_names.astype(str))):
        raise ContractError("V16 embedding and prepared cell rosters differ")
    mapped = full[cell_ids].copy()
    if not np.array_equal(mapped.obs_names.astype(str).to_numpy(), np.asarray(cell_ids)):
        raise ContractError("V16 prepared cell order differs")
    _set_model_labels(
        mapped,
        query_unknown=True,
        unlabeled=config["features"]["unlabeled_category"],
    )
    adapted_probabilities = adapted.predict(
        mapped, soft=True, batch_size=config["architecture"]["batch_size"]
    )
    with torch.no_grad():
        for key in classifier:
            adapted_state[key].copy_(reference_state[key])
    restored_state = adapted.module.state_dict()
    if (
        _digest(restored_state, nonclassifier) != nonclassifier_before
        or _digest(restored_state, classifier) != reference_classifier_digest
    ):
        raise ContractError("V16 restore changed nonclassifier state or failed")
    reference_probabilities = adapted.predict(
        mapped, soft=True, batch_size=config["architecture"]["batch_size"]
    )
    if (
        list(adapted_probabilities.columns) != list(reference_probabilities.columns)
        or adapted_probabilities.shape != reference_probabilities.shape
        or adapted_probabilities.shape[0] != len(v14_cells)
    ):
        raise ContractError("V16 probability outputs differ in class or cell roster")
    probabilities = (
        adapted_probabilities.to_numpy(dtype=np.float64)
        + reference_probabilities.to_numpy(dtype=np.float64)
    ) * 0.5
    if not np.isfinite(probabilities).all():
        raise ContractError("V16 consensus contains non-finite probabilities")
    classes = np.asarray(adapted_probabilities.columns.astype(str))
    predictions = classes[np.argmax(probabilities, axis=1)]
    reference_mask = v14_cells["strict_reference"].to_numpy(dtype=bool)
    predictions[reference_mask] = v14_cells.loc[
        reference_mask, "predicted_cell_type"
    ].astype(str).to_numpy()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    source_latent = v14_path.parent / v14_info["latent_file"]
    latent_path = output / v14_info["latent_file"]
    shutil.copy2(source_latent, latent_path)
    if sha256_path(latent_path) != sha256_path(source_latent):
        raise ContractError("V16 latent copy differs")
    output_cells = v14_cells.copy()
    output_cells["predicted_cell_type"] = predictions
    cells_path = output / v14_info["cells_file"]
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(output_cells.columns)
        writer.writerows(output_cells.itertuples(index=False, name=None))
    embedding = dict(v14_info)
    embedding.update(
        {
            "method": "geometry_preserving_continual_adapter_with_dual_head_consensus",
            "latent_sha256": sha256_path(latent_path),
            "cells_sha256": sha256_path(cells_path),
        }
    )
    write_json_exclusive(output / "embedding_manifest.json", embedding)
    manifest = {
        "schema_version": "masld-cl-dual-head-consensus-v16",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "query_labels_used": False,
        "weights": {"adapted": 0.5, "reference": 0.5},
        "class_order": classes.tolist(),
        "classifier_keys": classifier,
        "adapted_classifier_digest": adapted_classifier_digest,
        "reference_classifier_digest": reference_classifier_digest,
        "nonclassifier_digest_before": nonclassifier_before,
        "nonclassifier_digest_after": _digest(restored_state, nonclassifier),
        "latent_bitwise_identical_to_v14": True,
        "v14_embedding": {"path": str(v14_path), "sha256": sha256_path(v14_path)},
        "prepared": {"path": str(Path(args.prepared).resolve()), "sha256": sha256_path(args.prepared)},
        "embedding": embedding,
    }
    write_json_exclusive(output / "dual_head_manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
