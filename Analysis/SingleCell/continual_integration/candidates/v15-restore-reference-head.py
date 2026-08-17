#!/usr/bin/env python3
"""Restore the strict-reference scANVI head and re-export V14 predictions."""

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


def _tensor_digest(state: dict[str, torch.Tensor], keys: list[str]) -> str:
    digest = hashlib.sha256()
    for key in keys:
        value = state[key].detach().cpu().contiguous()
        digest.update(key.encode())
        digest.update(str(value.dtype).encode())
        digest.update(np.asarray(value.shape, dtype=np.int64).tobytes())
        digest.update(value.numpy().tobytes())
    return digest.hexdigest()


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
    if (
        policy.get("schema_version")
        != "masld-cl-reference-head-preservation-policy-v15"
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("classifier_operation", {}).get("query_labels_used") is not False
    ):
        raise ContractError("V15 policy identity differs")

    v14_path = Path(args.v14_embedding).resolve()
    v14_info, v14_latent, v14_cells = load_embedding(v14_path)
    if (
        v14_info.get("method") != "geometry_preserving_continual_adapter"
        or not v14_cells["strict_reference"].any()
    ):
        raise ContractError("V15 V14 source differs")

    set_all_seeds(17)
    reference = SCANVI.load(Path(args.reference_model).resolve())
    adapted = SCANVI.load(Path(args.adapted_model).resolve())
    reference_state = reference.module.state_dict()
    adapted_state = adapted.module.state_dict()
    reference_classifier = sorted(
        key for key in reference_state if key.startswith("classifier.")
    )
    adapted_classifier = sorted(
        key for key in adapted_state if key.startswith("classifier.")
    )
    if (
        not reference_classifier
        or reference_classifier != adapted_classifier
        or any(
            reference_state[key].shape != adapted_state[key].shape
            for key in reference_classifier
        )
    ):
        raise ContractError("reference and adapted classifier states do not match exactly")
    nonclassifier = sorted(set(adapted_state) - set(adapted_classifier))
    nonclassifier_before = _tensor_digest(adapted_state, nonclassifier)
    classifier_before = _tensor_digest(adapted_state, adapted_classifier)
    classifier_source = _tensor_digest(reference_state, reference_classifier)
    with torch.no_grad():
        for key in adapted_classifier:
            adapted_state[key].copy_(reference_state[key])
    current_state = adapted.module.state_dict()
    nonclassifier_after = _tensor_digest(current_state, nonclassifier)
    classifier_after = _tensor_digest(current_state, adapted_classifier)
    if nonclassifier_before != nonclassifier_after or classifier_after != classifier_source:
        raise ContractError("V15 changed a nonclassifier tensor or failed to restore the head")

    full = ad.read_h5ad(Path(args.prepared).resolve())
    cell_ids = v14_cells["cell_id"].astype(str).tolist()
    if len(cell_ids) != len(set(cell_ids)) or not set(cell_ids).issubset(set(full.obs_names.astype(str))):
        raise ContractError("V15 embedding and prepared cell rosters differ")
    mapped = full[cell_ids].copy()
    if not np.array_equal(mapped.obs_names.astype(str).to_numpy(), np.asarray(cell_ids)):
        raise ContractError("V15 prepared cell order differs")
    _set_model_labels(
        mapped,
        query_unknown=True,
        unlabeled=config["features"]["unlabeled_category"],
    )
    predictions = np.asarray(adapted.predict(mapped), dtype=object).astype(str)
    if len(predictions) != len(v14_cells):
        raise ContractError("V15 prediction roster differs")
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
        raise ContractError("V15 latent copy is not bitwise exact")
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
            "method": "geometry_preserving_continual_adapter_with_reference_head",
            "latent_sha256": sha256_path(latent_path),
            "cells_sha256": sha256_path(cells_path),
        }
    )
    write_json_exclusive(output / "embedding_manifest.json", embedding)
    manifest = {
        "schema_version": "masld-cl-reference-head-preservation-v15",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "query_labels_used": False,
        "classifier_training_or_calibration": False,
        "classifier_keys": adapted_classifier,
        "classifier_digest_before": classifier_before,
        "classifier_source_digest": classifier_source,
        "classifier_digest_after": classifier_after,
        "nonclassifier_digest_before": nonclassifier_before,
        "nonclassifier_digest_after": nonclassifier_after,
        "latent_bitwise_identical_to_v14": True,
        "v14_embedding": {"path": str(v14_path), "sha256": sha256_path(v14_path)},
        "reference_model": str(Path(args.reference_model).resolve()),
        "adapted_model": str(Path(args.adapted_model).resolve()),
        "prepared": {"path": str(Path(args.prepared).resolve()), "sha256": sha256_path(args.prepared)},
        "embedding": embedding,
    }
    write_json_exclusive(output / "reference_head_manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
