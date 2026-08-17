#!/usr/bin/env python3
"""Run the frozen V25 reference-only kNN audit on named sensitivity embeddings."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from masld_cl.config import load_config, write_json_exclusive
from masld_cl.contracts import ContractError, sha256_path
from masld_cl.embedding import load_embedding
from masld_cl.latent_knn_label_audit import (
    _fit_predict_knn,
    _prediction_scores,
    load_latent_knn_label_audit_policy,
)
from masld_cl.training import _capped_indices


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--source", action="append", nargs=2, metavar=("NAME", "MANIFEST"), required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    policy_path, policy = load_latent_knn_label_audit_policy(config, args.policy)
    classifier = policy["classifier"]
    sources = {name: Path(path).resolve() for name, path in args.source}
    if len(sources) != len(args.source):
        raise ContractError("sensitivity method names are duplicated")

    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    results = {}
    query_roster_sha256 = None
    for name, source in sources.items():
        info, latent, cells = load_embedding(source)
        reference = cells["strict_reference"].to_numpy(dtype=bool)
        query = cells["analysis_eligible"].to_numpy(dtype=bool) & ~reference
        if reference.sum() != 216957 or query.sum() != 687559:
            raise ContractError(f"sensitivity cell roster differs: {name}")
        roster_sha256 = hashlib.sha256(
            "\0".join(cells.loc[query, "cell_id"].astype(str)).encode("utf-8")
        ).hexdigest()
        if query_roster_sha256 is None:
            query_roster_sha256 = roster_sha256
        elif roster_sha256 != query_roster_sha256:
            raise ContractError("sensitivity query rosters differ")
        reference_pool = _capped_indices(
            type("ADataView", (), {"obs": cells, "n_obs": len(cells)})(),
            np.flatnonzero(reference), "all_lineage", config,
            classifier["reference_cap_seed"],
        )
        if len(reference_pool) != 61960:
            raise ContractError("sensitivity capped reference roster differs")
        predictions, identity = _fit_predict_knn(
            np.asarray(latent)[reference_pool],
            cells.iloc[reference_pool]["audit_cell_type"].astype(str).to_numpy(),
            np.asarray(latent)[query], classifier,
        )
        scores = _prediction_scores(
            cells.loc[query].reset_index(drop=True), predictions, config["lineages"]
        )
        results[name] = {
            "embedding": str(source),
            "embedding_sha256": sha256_path(source),
            "embedding_method": info.get(
                "method", info.get("training_method", info["model_kind"])
            ),
            "mean_donor_macro_f1": scores["mean_donor_macro_f1"],
            "mean_donor_lineage_f1": scores["mean_donor_lineage_f1"],
            "knn_identity": identity,
        }
        write_json_exclusive(output / f"progress_{name}.json", results[name])

    result = {
        "schema_version": "masld-cl-latent-knn-sensitivity-v25",
        "config_sha256": config["_config_sha256"],
        "v25_policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "reference_only_fit": True,
        "query_labels_available_to_fit": False,
        "query_labels_used_only_for_postfit_scoring": True,
        "sensitivity_only": True,
        "query_roster_sha256": query_roster_sha256,
        "results": results,
    }
    write_json_exclusive(output / "sensitivity.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
