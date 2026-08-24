#!/usr/bin/env python
"""SCimilarity common-lane adapter for the frozen cell-state task.

The encoder is frozen.  Only the common head is trained, and it is trained on
outer-training donors alone, so the held-out donors of the fold never influence
any fitted state.

Preprocessing follows the registered input contract exactly, and it is executed
by SCimilarity's own utilities rather than reimplemented here:

* ``consolidate_duplicate_symbols`` -- duplicate gene symbols are summed
  *before* alignment.
* ``align_dataset(..., gene_overlap_threshold=5000)`` -- the matrix is aligned
  to the archive's exact 28,231-gene order, zero-filling unobserved genes, and
  the call itself enforces the registered ``minimum_gene_overlap`` of 5000.
* ``lognorm_counts`` -- per-cell total 10,000 then natural log1p.

The registry warns that zero-filled genes are *structural missingness before
alignment, not observed biological zeros*.  This adapter therefore records the
observed overlap and the aligned-zero mask fraction in its receipt instead of
letting a low-overlap run look like a well-supported one.

Because the encoder is frozen and label-free, embeddings are computed once in
``prepare`` for every row and reused by ``fit`` and ``predict``.  That is
equivalent to embedding per action -- no fold or label information reaches the
encoder -- and avoids recomputing an identical projection ten times.

The atlas stores Ensembl ids in ``var_names`` and gene symbols in
``source_feature_id``; SCimilarity keys on symbols, so the symbol column is
used and the mapping is asserted rather than assumed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from hvg_pca_logistic import (  # type: ignore[import-not-found]
    ScientificAdapterError,
    _artifact_record,
    _input_by_role,
    _require_sha256,
    _validate_file_artifact,
    _write_json,
    _write_tsv,
    fold_index,
    join_hash,
)

TASK_ID = "cell_state_mapping"
DATASET_ID = "resource_atlas_current"
RUNTIME_ID = "gpu_scimilarity"
RECEIPT_SCHEMA = "masld-bench-adapter-receipt-v1"
ADAPTER_ID = "scimilarity_common_lane_v1"
MODEL_IDS = ("scimilarity_v1_1",)
HEAD_IDS = ("linear", "two_layer_mlp")
EMBEDDING_POLICY = "unit_l2_hypersphere_encoder_latent"
GENE_SYMBOL_COLUMN = "source_feature_id"
MINIMUM_GENE_OVERLAP = 5000
LATENT_DIMENSION = 128
INPUT_DIMENSION = 28231


def _select_device(requested: str) -> str:
    import torch

    if requested == "cpu":
        return "cpu"
    if torch.cuda.is_available():
        return "cuda"
    if requested == "cuda":
        raise ScientificAdapterError("CUDA was required but is unavailable")
    return "cpu"


def _build_head(head_id: str, dimensions: int, classes: int, seed: int):
    import torch
    from torch import nn

    torch.manual_seed(seed)
    if head_id == "linear":
        return nn.Linear(dimensions, classes)
    if head_id == "two_layer_mlp":
        return nn.Sequential(
            nn.Linear(dimensions, 256), nn.ReLU(), nn.Linear(256, classes)
        )
    raise ScientificAdapterError(f"unsupported common head: {head_id}")


def _validate_request(request: Mapping[str, Any], action: str):
    if request.get("schema_version") != "masld-bench-adapter-request-v1":
        raise ScientificAdapterError("unsupported adapter request schema")
    if request.get("action") != action:
        raise ScientificAdapterError("request action differs")
    _require_sha256(request.get("run_id"), "run_id")
    run_spec = request.get("run_spec")
    if not isinstance(run_spec, Mapping):
        raise ScientificAdapterError("request lacks a run_spec")
    for field, expected in (
        ("task_id", TASK_ID),
        ("adaptation_regime", "common_lane"),
        ("runtime_id", RUNTIME_ID),
    ):
        if run_spec.get(field) != expected:
            raise ScientificAdapterError(f"RunSpec {field} differs from {expected}")
    model_id = str(run_spec.get("model_id", ""))
    if model_id not in MODEL_IDS:
        raise ScientificAdapterError(f"unsupported model_id: {model_id}")
    parameters = run_spec.get("hyperparameters")
    if not isinstance(parameters, Mapping):
        raise ScientificAdapterError("run_spec lacks hyperparameters")
    for field in ("join_namespace", "outer_folds", "cell_state_roster"):
        if field not in parameters:
            raise ScientificAdapterError(f"hyperparameters lack {field}")
    if parameters.get("embedding_policy", EMBEDDING_POLICY) != EMBEDDING_POLICY:
        raise ScientificAdapterError("embedding policy differs from the frozen contract")
    head_id = str(parameters.get("common_head_id", ""))
    if head_id not in HEAD_IDS:
        raise ScientificAdapterError(f"unsupported common head: {head_id}")
    roster = tuple(str(item) for item in parameters["cell_state_roster"])
    if len(set(roster)) != len(roster) or not roster:
        raise ScientificAdapterError("cell state roster is invalid")
    return parameters, roster, model_id, head_id


def _model_root(request: Mapping[str, Any], model_id: str) -> Path:
    """Anchor on encoder.ckpt and use its directory as the model path."""

    role = f"model_checkpoint_weights:{model_id}"
    encoder = _validate_file_artifact(_input_by_role(request, role), "SCimilarity encoder")
    if encoder.name != "encoder.ckpt":
        raise ScientificAdapterError("SCimilarity encoder artifact must be encoder.ckpt")
    root = encoder.parent
    for required in ("gene_order.tsv", "layer_sizes.json", "hyperparameters.json"):
        if not (root / required).is_file():
            raise ScientificAdapterError(f"SCimilarity model directory lacks {required}")
    return root


def _atlas_path(request: Mapping[str, Any]) -> Path:
    role = f"dataset_view_data:resource_atlas_geneformer_smoke_1000_v1"
    return _validate_file_artifact(_input_by_role(request, role), "atlas dataset view")


def _load_atlas(path: Path):
    import anndata

    adata = anndata.read_h5ad(path)
    if not {"donor_id", "broad_label"}.issubset(adata.obs.columns):
        raise ScientificAdapterError("Atlas H5AD lacks required observation fields")
    if GENE_SYMBOL_COLUMN not in adata.var.columns:
        raise ScientificAdapterError(
            f"Atlas H5AD lacks the {GENE_SYMBOL_COLUMN!r} gene-symbol column that "
            "SCimilarity's fixed gene order keys on"
        )
    if len(set(map(str, adata.obs_names))) != adata.n_obs:
        raise ScientificAdapterError("Atlas row identifiers are not unique")
    return adata


def _split_rows(adata: Any, parameters: Mapping[str, Any], run_spec: Mapping[str, Any]):
    outer_folds = int(parameters["outer_folds"])
    fold = int(run_spec["fold"])
    if not 0 <= fold < outer_folds:
        raise ScientificAdapterError("outer-fold contract differs")
    namespace = str(parameters["join_namespace"])
    rows = []
    for row_id, donor in zip(
        map(str, adata.obs_names), map(str, adata.obs["donor_id"])
    ):
        assigned = fold_index(donor, seed=20260821, outer_folds=outer_folds)
        rows.append(
            {
                "row_hash": join_hash(namespace, "row", row_id),
                "unit_hash": join_hash(namespace, "unit", donor),
                "fold": assigned,
                "held_out": assigned == fold,
            }
        )
    rows.sort(key=lambda item: item["row_hash"])
    if len({r["row_hash"] for r in rows}) != len(rows):
        raise ScientificAdapterError("row identifiers collide")
    return rows


def _labels_by_hash(request, parameters, roster):
    adata = _load_atlas(_atlas_path(request))
    namespace = str(parameters["join_namespace"])
    index = {label: position for position, label in enumerate(roster)}
    observed = set(map(str, adata.obs["broad_label"]))
    if not observed.issubset(set(roster)):
        raise ScientificAdapterError("Atlas label set differs from the frozen roster")
    return {
        join_hash(namespace, "row", str(name)): index[str(label)]
        for name, label in zip(adata.obs_names, adata.obs["broad_label"])
    }


def _embed_all(adata: Any, model_root: Path, device: str):
    """Align, log-normalise and embed every row with SCimilarity's own code.

    Returns the embedding matrix in ``adata.obs_names`` order plus an
    alignment summary.  No label or fold information is used, so computing this
    once for all rows is equivalent to embedding each action's rows separately.
    """

    import numpy as np
    import scipy.sparse as sp
    from scimilarity import CellEmbedding
    from scimilarity.utils import (
        align_dataset,
        consolidate_duplicate_symbols,
        lognorm_counts,
    )

    working = adata.copy()
    # SCimilarity keys on gene symbols; the atlas stores Ensembl ids in var_names.
    working.var_names = [str(symbol) for symbol in working.var[GENE_SYMBOL_COLUMN]]

    # ``lognorm_counts`` reads layers["counts"], while this atlas carries raw
    # counts in X.  Verify they really are raw UMI counts before promoting them,
    # rather than silently log-normalising something already normalised.
    raw = working.X
    sample = raw[:64].toarray() if sp.issparse(raw) else np.asarray(raw[:64])
    if sample.size:
        if float(sample.min()) < 0.0:
            raise ScientificAdapterError(
                "atlas X holds negative values; the registered input contract "
                "requires raw nonnegative UMI counts"
            )
        if not np.allclose(sample, np.rint(sample)):
            raise ScientificAdapterError(
                "atlas X is not integer-valued, so it is not the raw count matrix "
                "the registered input contract requires"
            )
    working.layers["counts"] = working.X

    embedder = CellEmbedding(model_path=str(model_root))
    gene_order = list(embedder.gene_order)
    if len(gene_order) != INPUT_DIMENSION:
        raise ScientificAdapterError(
            f"SCimilarity gene order is {len(gene_order)}, expected {INPUT_DIMENSION}"
        )

    # Duplicate symbols are summed BEFORE alignment, per the input contract.
    working = consolidate_duplicate_symbols(working)

    observed = {str(name) for name in working.var_names}
    overlap = sum(1 for gene in gene_order if gene in observed)
    if overlap < MINIMUM_GENE_OVERLAP:
        raise ScientificAdapterError(
            f"SCimilarity gene overlap {overlap} is below the registered minimum "
            f"{MINIMUM_GENE_OVERLAP}; the aligned matrix would be mostly structural "
            "zeros rather than observed counts"
        )

    aligned = align_dataset(
        working, gene_order, gene_overlap_threshold=MINIMUM_GENE_OVERLAP
    )
    aligned = lognorm_counts(aligned)
    matrix = aligned.X
    embeddings = embedder.get_embeddings(matrix)
    embeddings = np.asarray(embeddings, dtype=np.float32)
    if embeddings.shape != (adata.n_obs, LATENT_DIMENSION):
        raise ScientificAdapterError(
            f"SCimilarity embeddings are {embeddings.shape}, expected "
            f"{(adata.n_obs, LATENT_DIMENSION)}"
        )

    # Structural-missingness accounting: what fraction of the aligned matrix is
    # a zero that was filled in rather than observed?
    dense_zero_fraction = 1.0 - (overlap / float(len(gene_order)))
    counts = adata.X
    counts = counts.tocsc() if sp.issparse(counts) else counts
    totals = np.asarray(counts.sum(axis=0)).ravel()
    kept = np.array(
        [str(symbol) in set(gene_order) for symbol in adata.var[GENE_SYMBOL_COLUMN]]
    )
    umi_retained = float(totals[kept].sum() / totals.sum()) if totals.sum() else 0.0

    summary = {
        "gene_order_length": len(gene_order),
        "observed_symbols": len(observed),
        "gene_overlap": overlap,
        "gene_overlap_fraction": overlap / float(len(gene_order)),
        "aligned_zero_fill_fraction": dense_zero_fraction,
        "umi_weighted_coverage_retained": umi_retained,
        "minimum_gene_overlap_enforced": MINIMUM_GENE_OVERLAP,
        "structural_missingness_note": (
            "Zero-filled genes are absent from the assay, not observed biological "
            "zeros; aligned_zero_fill_fraction is the share of the model's input "
            "vector that carries no observation."
        ),
    }
    return embeddings, summary


def _receipt(*, action, request, output, artifacts, extra_metadata) -> None:
    run_spec = request["run_spec"]
    payload = {
        "schema_version": RECEIPT_SCHEMA,
        "adapter_id": ADAPTER_ID,
        "action": action,
        "run_id": request["run_id"],
        "model_id": run_spec["model_id"],
        "task_id": run_spec["task_id"],
        "seed": int(run_spec["seed"]),
        "fold": int(run_spec["fold"]),
        "artifacts": [
            _artifact_record(path, relative_to=output) for path in artifacts
        ],
        "metadata": dict(extra_metadata),
    }
    _write_json(output / "adapter_receipt.json", payload)


def _prior_output(request: Mapping[str, Any], action: str) -> Path:
    prior = request.get("prior_action_outputs")
    if not isinstance(prior, Mapping) or action not in prior:
        raise ScientificAdapterError(f"request lacks a prior {action} output")
    path = Path(str(prior[action]))
    if not path.is_dir():
        raise ScientificAdapterError(f"prior {action} output is missing")
    return path


def _read_split(path: Path) -> list[dict[str, Any]]:
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["row_hash", "unit_hash", "fold", "held_out"]:
            raise ScientificAdapterError("split_rows.tsv header differs")
        for line in handle:
            row_hash, unit_hash, fold, held_out = line.rstrip("\n").split("\t")
            rows.append(
                {
                    "row_hash": row_hash,
                    "unit_hash": unit_hash,
                    "fold": int(fold),
                    "held_out": held_out == "true",
                }
            )
    if not rows:
        raise ScientificAdapterError("split_rows.tsv is empty")
    return rows


def _load_embeddings(prepare_dir: Path):
    import numpy as np

    matrix = np.load(prepare_dir / "embeddings.npy")
    order = [
        line.rstrip("\n")
        for line in open(prepare_dir / "embedding_row_order.txt", encoding="utf-8")
    ]
    if matrix.shape[0] != len(order):
        raise ScientificAdapterError("embedding matrix and row order disagree")
    return {row_hash: matrix[index] for index, row_hash in enumerate(order)}


def prepare(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import numpy as np

    parameters, roster, model_id, _ = _validate_request(request, "prepare")
    adata = _load_atlas(_atlas_path(request))
    split = _split_rows(adata, parameters, request["run_spec"])
    namespace = str(parameters["join_namespace"])
    device = _select_device(str(parameters.get("device", "auto")))
    model_root = _model_root(request, model_id)

    embeddings, summary = _embed_all(adata, model_root, device)

    # Reorder to the sorted row_hash order used everywhere downstream.
    by_hash = {
        join_hash(namespace, "row", str(name)): embeddings[index]
        for index, name in enumerate(adata.obs_names)
    }
    order = [row["row_hash"] for row in split]
    matrix = np.stack([by_hash[row_hash] for row_hash in order]).astype(np.float32)

    np.save(output / "embeddings.npy", matrix)
    with open(output / "embedding_row_order.txt", "w", encoding="utf-8") as handle:
        for row_hash in order:
            handle.write(f"{row_hash}\n")
    _write_tsv(
        output / "split_rows.tsv",
        ("row_hash", "unit_hash", "fold", "held_out"),
        [
            {
                "row_hash": row["row_hash"],
                "unit_hash": row["unit_hash"],
                "fold": row["fold"],
                "held_out": "true" if row["held_out"] else "false",
            }
            for row in split
        ],
    )
    _write_json(output / "alignment_summary.json", summary)
    _receipt(
        action="prepare",
        request=request,
        output=output,
        artifacts=[
            output / "embeddings.npy",
            output / "embedding_row_order.txt",
            output / "split_rows.tsv",
            output / "alignment_summary.json",
        ],
        extra_metadata={
            "embedding_policy": EMBEDDING_POLICY,
            "embedding_dimensions": int(matrix.shape[1]),
            "device": device,
            **summary,
        },
    )


def fit(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import torch
    from collections import Counter

    parameters, roster, model_id, head_id = _validate_request(request, "fit")
    prepare_dir = _prior_output(request, "prepare")
    split = _read_split(prepare_dir / "split_rows.tsv")
    embeddings = _load_embeddings(prepare_dir)
    labels = _labels_by_hash(request, parameters, roster)

    training = [row for row in split if not row["held_out"]]
    if not training:
        raise ScientificAdapterError("fold has no training donors")
    held_units = {row["unit_hash"] for row in split if row["held_out"]}
    if held_units & {row["unit_hash"] for row in training}:
        raise ScientificAdapterError("a donor appears in both training and held-out")

    device = _select_device(str(parameters.get("device", "auto")))
    seed = int(request["run_spec"]["seed"])
    torch.manual_seed(seed)

    features = torch.tensor(
        [embeddings[row["row_hash"]] for row in training], dtype=torch.float32
    )
    targets = torch.tensor(
        [labels[row["row_hash"]] for row in training], dtype=torch.long
    )

    donor_counts = Counter(row["unit_hash"] for row in training)
    class_counts = Counter(int(t) for t in targets)
    weights = torch.tensor(
        [
            1.0 / (donor_counts[row["unit_hash"]] * class_counts[int(target)])
            for row, target in zip(training, targets)
        ],
        dtype=torch.float32,
    )
    weights = weights / weights.sum() * len(weights)

    head = _build_head(head_id, features.shape[1], len(roster), seed).to(device)
    optimizer = torch.optim.Adam(
        head.parameters(), lr=float(parameters.get("head_learning_rate", 1e-3))
    )
    epochs = int(parameters.get("head_epochs", 60))
    features_device = features.to(device)
    target_device = targets.to(device)
    weight_device = weights.to(device)
    losses = []
    for _ in range(epochs):
        optimizer.zero_grad()
        logits = head(features_device)
        per_row = torch.nn.functional.cross_entropy(
            logits, target_device, reduction="none"
        )
        loss = (per_row * weight_device).mean()
        loss.backward()
        optimizer.step()
        losses.append(float(loss.detach().cpu()))

    state_path = output / "head_state.pt"
    torch.save(
        {
            "head_id": head_id,
            "state_dict": {k: v.cpu() for k, v in head.state_dict().items()},
            "embedding_dimensions": int(features.shape[1]),
            "roster": list(roster),
            "seed": seed,
        },
        state_path,
    )
    _receipt(
        action="fit",
        request=request,
        output=output,
        artifacts=[state_path],
        extra_metadata={
            "head_id": head_id,
            "training_rows": len(training),
            "training_donors": len(donor_counts),
            "epochs": epochs,
            "final_loss": losses[-1] if losses else None,
            "training_weight_policy": "donor_class_balanced_rescaled_to_n_cells",
            "device": device,
        },
    )


def predict(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import torch

    parameters, roster, model_id, head_id = _validate_request(request, "predict")
    prepare_dir = _prior_output(request, "prepare")
    fit_dir = _prior_output(request, "fit")
    split = _read_split(prepare_dir / "split_rows.tsv")
    embeddings = _load_embeddings(prepare_dir)

    held = [row for row in split if row["held_out"]]
    if not held:
        raise ScientificAdapterError("fold has no held-out rows")

    device = _select_device(str(parameters.get("device", "auto")))
    state = torch.load(fit_dir / "head_state.pt", map_location="cpu", weights_only=True)
    if list(state["roster"]) != list(roster):
        raise ScientificAdapterError("fitted roster differs from the request roster")
    head = _build_head(
        head_id, int(state["embedding_dimensions"]), len(roster), int(state["seed"])
    )
    head.load_state_dict(state["state_dict"], strict=True)
    head.to(device).eval()

    features = torch.tensor(
        [embeddings[row["row_hash"]] for row in held], dtype=torch.float32
    ).to(device)
    with torch.no_grad():
        probabilities = torch.softmax(head(features), dim=1).cpu().numpy()

    prediction_rows = []
    probability_rows = []
    for row, values in zip(held, probabilities):
        # Deterministic tie-break: highest probability, then roster order.
        best = min(range(len(roster)), key=lambda i: (-values[i], roster[i]))
        prediction_rows.append(
            {
                "row_hash": row["row_hash"],
                "unit_hash": row["unit_hash"],
                "predicted_class": roster[best],
            }
        )
        probability_rows.append(
            {
                "row_hash": row["row_hash"],
                "unit_hash": row["unit_hash"],
                **{
                    f"probability::{label}": repr(float(values[index]))
                    for index, label in enumerate(roster)
                },
            }
        )

    _write_tsv(
        output / "predictions.tsv",
        ("row_hash", "unit_hash", "predicted_class"),
        prediction_rows,
    )
    _write_tsv(
        output / "row_ids.tsv",
        ("row_hash", "unit_hash"),
        [
            {"row_hash": row["row_hash"], "unit_hash": row["unit_hash"]}
            for row in held
        ],
    )
    _write_tsv(
        output / "class_probabilities.tsv",
        ("row_hash", "unit_hash", *[f"probability::{label}" for label in roster]),
        probability_rows,
    )
    _receipt(
        action="predict",
        request=request,
        output=output,
        artifacts=[
            output / "predictions.tsv",
            output / "row_ids.tsv",
            output / "class_probabilities.tsv",
        ],
        extra_metadata={
            "head_id": head_id,
            "held_out_rows": len(held),
            "held_out_donors": len({row["unit_hash"] for row in held}),
            "device": device,
        },
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--action", required=True, choices=("prepare", "fit", "predict"))
    parser.add_argument("--request", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args(argv)

    if arguments.output.exists():
        raise ScientificAdapterError(
            f"adapter output already exists: {arguments.output}"
        )
    request = json.loads(arguments.request.read_text(encoding="utf-8"))
    arguments.output.mkdir(parents=True)
    handler = {"prepare": prepare, "fit": fit, "predict": predict}[arguments.action]
    handler(arguments.request, request, arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
