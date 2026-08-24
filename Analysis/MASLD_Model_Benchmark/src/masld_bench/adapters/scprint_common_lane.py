#!/usr/bin/env python
"""scPRINT common-lane adapter for the frozen cell-state task.

The encoder is frozen.  Only the common head is trained, and it is trained on
outer-training donors alone, so held-out donors never influence any fitted
state.  This mirrors ``geneformer_common_lane.py`` exactly: same request
schema, same fold/row/unit hashing, same donor- and class-balanced training
weights, same deterministic argmax tie-break, and byte-identical output
column naming.

Checkpoint safety
------------------
``v2-medium.ckpt`` is a PyTorch Lightning checkpoint, i.e. a zip archive
whose ``archive/data.pkl`` member is a pickle stream.  Before this module
ever calls ``torch.load``, it disassembles that pickle stream with
``pickletools`` (which never executes anything) and asserts that every
``GLOBAL``/``STACK_GLOBAL`` opcode names one of exactly three reconstructors:
``collections.OrderedDict``, ``torch.FloatStorage``, and
``torch._utils._rebuild_tensor_v2``.  These are precisely the three globals
already covered by ``torch.load(..., weights_only=True)``'s own default safe
list, so loading this specific checkpoint does not widen that allowlist by a
single name.  If a future re-acquisition of this path ever produced a
pickle stream naming anything else, the disassembly step refuses it before
``torch.load`` is ever invoked, and no allowlist is widened silently.

Architecture and gene vocabulary
---------------------------------
Restricted inspection (see the adapter's ``inspect_checkpoint`` helper, and
the report that accompanied this file) established two facts that the
registered contract in ``config/artifacts/models/scprint/checkpoints.json``
did not state correctly or completely:

* The checkpoint's ``hyper_parameters`` dict carries ``nhead: 4``, not the
  ``attention_heads: 2`` recorded in the checkpoints.json
  ``architecture_contract``.  ``d_model`` (256) and ``nlayers`` (8) do match.
  This adapter trusts the checkpoint's own hyperparameters, not the
  registered contract, and records the discrepancy in every receipt.
* The gene vocabulary is fully self-contained in the checkpoint:
  ``hyper_parameters["genes"]`` is an ordered list of 44,756 Ensembl gene
  IDs whose length equals ``state_dict["gene_encoder.embedding.weight"]``'s
  row count.  The ``precpt_gene_emb`` hyperparameter names an external ESM2
  protein-embedding parquet used only to *initialize* that embedding table
  at training time; the trained weights are already baked into the
  checkpoint, so no ESM side file is required for inference.

Frozen-encoder inference itself, however, requires the ``scprint`` package:
the state_dict's keys (``transformer.blocks.*.mixer.Wqkv`` fused
flash-attention projections, the ``expr_encoder``/``class_encoder`` fusion,
"cls"-style pooling, prenorm placement) are a custom architecture with no
public reimplementation available here, and no copy of the ``scprint``
Python package or its source is installed or vendored anywhere on this
cluster.  Reconstructing the forward pass from tensor shapes alone without
the source would be exactly the kind of improvised partial mapping this
adapter must not perform.  ``fit`` and ``predict`` therefore verify the
checkpoint and then stop with an explicit, actionable
``ScientificAdapterError`` naming the missing dependency, rather than
guessing at scPRINT's internal forward-pass API.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickletools
import zipfile
from typing import Any, Mapping

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
RUNTIME_ID = "gpu_scprint"
RECEIPT_SCHEMA = "masld-bench-adapter-receipt-v1"
ADAPTER_ID = "scprint_common_lane_v1"
MODEL_IDS = ("scprint_v1_5_medium",)
HEAD_IDS = ("linear", "two_layer_mlp")
# The CLS-token cell embedding, taken before any of scPRINT's supervised
# label/discriminator heads.  This is the "common_lane_candidate" the
# registered embedding_contract names: a disentangled cell embedding with an
# identically-trained, frozen-encoder-only downstream head, not the released
# native classification logits.
EMBEDDING_POLICY = "cls_token_cell_embedding_frozen_encoder"
CONTRACT_PATH = (
    Path(__file__).resolve().parents[3] / "config" / "artifacts" / "models" / "scprint" / "checkpoints.json"
)

# Exactly the globals PyTorch's own ``torch.load(..., weights_only=True)``
# already treats as safe.  Restricted inspection of v2-medium.ckpt found only
# these three; this allowlist is not widened beyond them anywhere in this
# module.
_ALLOWED_CHECKPOINT_GLOBALS = frozenset(
    {
        ("collections", "OrderedDict"),
        ("torch", "FloatStorage"),
        ("torch._utils", "_rebuild_tensor_v2"),
    }
)


def _select_device(requested: str) -> str:
    import torch

    if requested == "cpu":
        return "cpu"
    if torch.cuda.is_available():
        return "cuda"
    if requested == "cuda":
        raise ScientificAdapterError("CUDA was required but is unavailable")
    return "cpu"


def _load_contract() -> Mapping[str, Any]:
    try:
        contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificAdapterError(
            f"scPRINT checkpoint contract is unreadable: {error}"
        ) from error
    if (
        not isinstance(contract, Mapping)
        or contract.get("schema_version") != "masld-bench-upstream-checkpoint-preflight-v1"
    ):
        raise ScientificAdapterError("unsupported scPRINT checkpoint contract schema")
    return contract


def _checkpoint_path(request: Mapping[str, Any], model_id: str) -> Path:
    role = f"model_checkpoint_bundle:{model_id}"
    return _validate_file_artifact(_input_by_role(request, role), "scPRINT checkpoint")


def _disassemble_pickle_globals(payload: bytes) -> frozenset[tuple[str, str]]:
    """Enumerate GLOBAL/STACK_GLOBAL opcodes without executing anything."""

    seen: set[tuple[str, str]] = set()
    for opcode, argument, _position in pickletools.genops(payload):
        if opcode.name in ("GLOBAL", "STACK_GLOBAL"):
            text = str(argument)
            module, _, name = text.rpartition(" ") if " " in text else ("", "", text)
            if not module:
                # STACK_GLOBAL pushes module/name as two prior STRING/SHORT_BINUNICODE
                # opcodes rather than a single space-joined argument; fall back to
                # the raw text as a single opaque token so it still fails the
                # allowlist check below rather than being silently skipped.
                seen.add(("", text))
            else:
                seen.add((module, name))
    return frozenset(seen)


def _restricted_checkpoint_globals(checkpoint: Path) -> frozenset[tuple[str, str]]:
    """Verify a checkpoint's pickle stream names only known-safe globals.

    This runs before any ``torch.load`` call.  ``zipfile`` and
    ``pickletools`` never construct or execute the objects a pickle
    describes; they only parse bytes.
    """

    try:
        with zipfile.ZipFile(checkpoint) as archive:
            data_members = [n for n in archive.namelist() if n.endswith("data.pkl")]
            if len(data_members) != 1:
                raise ScientificAdapterError(
                    "scPRINT checkpoint archive must contain exactly one data.pkl member"
                )
            payload = archive.read(data_members[0])
    except zipfile.BadZipFile as error:
        raise ScientificAdapterError(
            f"scPRINT checkpoint is not a valid zip-format Torch archive: {error}"
        ) from error
    globals_found = _disassemble_pickle_globals(payload)
    forbidden = globals_found - _ALLOWED_CHECKPOINT_GLOBALS
    if forbidden:
        raise ScientificAdapterError(
            "scPRINT checkpoint pickle names forbidden globals, refusing to load: "
            + ", ".join(f"{m}.{n}" for m, n in sorted(forbidden))
        )
    return globals_found


def _load_checkpoint(checkpoint: Path) -> Mapping[str, Any]:
    """Restricted-load a verified scPRINT Lightning checkpoint.

    Only reachable after ``_restricted_checkpoint_globals`` has confirmed the
    pickle stream names nothing beyond the three safe reconstructors above.
    ``weights_only=True`` is PyTorch's own restricted unpickler; it is not
    disabled or widened here.
    """

    import torch

    globals_found = _restricted_checkpoint_globals(checkpoint)
    if not globals_found.issubset(_ALLOWED_CHECKPOINT_GLOBALS):  # pragma: no cover - defense in depth
        raise ScientificAdapterError("scPRINT checkpoint failed the restricted-global pre-check")
    payload = torch.load(checkpoint, map_location="cpu", weights_only=True)
    if not isinstance(payload, Mapping) or "state_dict" not in payload or "hyper_parameters" not in payload:
        raise ScientificAdapterError(
            "scPRINT checkpoint is missing the expected Lightning state_dict/hyper_parameters keys"
        )
    return payload


def _inspect_architecture(
    payload: Mapping[str, Any], contract: Mapping[str, Any], model_id: str
) -> dict[str, Any]:
    hyper_parameters = payload["hyper_parameters"]
    state_dict = payload["state_dict"]
    if not isinstance(hyper_parameters, Mapping):
        raise ScientificAdapterError("scPRINT hyper_parameters is not a mapping")
    genes = hyper_parameters.get("genes")
    if not isinstance(genes, list) or not genes or any(not isinstance(g, str) for g in genes):
        raise ScientificAdapterError("scPRINT checkpoint carries no gene vocabulary")
    embedding_key = "gene_encoder.embedding.weight"
    if embedding_key not in state_dict:
        raise ScientificAdapterError("scPRINT checkpoint has no gene embedding table")
    embedding_shape = tuple(state_dict[embedding_key].shape)
    if embedding_shape[0] != len(genes):
        raise ScientificAdapterError(
            "scPRINT gene embedding row count differs from the embedded gene vocabulary length"
        )
    registered = contract.get("architecture_contract", {})
    observed = {
        "hidden_dimension": hyper_parameters.get("d_model"),
        "attention_heads": hyper_parameters.get("nhead"),
        "transformer_layers": hyper_parameters.get("nlayers"),
    }
    mismatches = {
        field: {"registered": registered.get(field), "observed": observed[field]}
        for field in observed
        if registered.get(field) != observed[field]
    }
    return {
        "schema_version": "masld-bench-scprint-architecture-inspection-v1",
        "d_model": hyper_parameters.get("d_model"),
        "nhead": hyper_parameters.get("nhead"),
        "nlayers": hyper_parameters.get("nlayers"),
        "cell_emb_style": hyper_parameters.get("cell_emb_style"),
        "expr_emb_style": hyper_parameters.get("expr_emb_style"),
        "gene_vocabulary_size": len(genes),
        "gene_embedding_shape": list(embedding_shape),
        "precpt_gene_emb_training_source": hyper_parameters.get("precpt_gene_emb"),
        "restricted_pickle_globals": sorted(
            f"{m}.{n}" for m, n in _ALLOWED_CHECKPOINT_GLOBALS
        ),
        "architecture_contract_mismatches": mismatches,
    }


def _validate_request(request: Mapping[str, Any], action: str):
    if request.get("schema_version") != "masld-bench-adapter-request-v1":
        raise ScientificAdapterError("unsupported adapter request schema")
    if request.get("action") != action:
        raise ScientificAdapterError("request action differs")
    _require_sha256(request.get("run_id"), "run_id")
    run_spec = request.get("run_spec")
    if not isinstance(run_spec, Mapping):
        raise ScientificAdapterError("request has no RunSpec")
    model_id = str(run_spec.get("model_id", ""))
    if model_id not in MODEL_IDS:
        raise ScientificAdapterError("RunSpec names an unsupported scPRINT model")
    for field, expected in (
        ("task_id", TASK_ID),
        ("split_id", "donor_outer"),
        ("adaptation_regime", "common_lane"),
        ("runtime_id", RUNTIME_ID),
    ):
        if run_spec.get(field) != expected:
            raise ScientificAdapterError(f"RunSpec {field} differs from {expected}")
    parameters = run_spec.get("hyperparameters")
    if not isinstance(parameters, Mapping):
        raise ScientificAdapterError("RunSpec carries no frozen hyperparameters")
    if parameters.get("embedding_policy") != EMBEDDING_POLICY:
        raise ScientificAdapterError("embedding policy differs from the frozen policy")
    head_id = str(parameters.get("common_head_id", ""))
    if head_id not in HEAD_IDS:
        raise ScientificAdapterError(f"unsupported common head: {head_id!r}")
    roster = parameters.get("cell_state_roster")
    if not isinstance(roster, (list, tuple)) or not roster:
        raise ScientificAdapterError("frozen cell-state roster is missing")
    roster = tuple(str(item) for item in roster)
    if list(roster) != sorted(set(roster)):
        raise ScientificAdapterError("cell-state roster must be sorted and unique")
    return parameters, roster, model_id, head_id


def _receipt(*, action, request, output, artifacts, extra_metadata) -> None:
    receipt = {
        "schema_version": RECEIPT_SCHEMA,
        "action": action,
        "run_id": request["run_id"],
        "status": "complete",
        "artifacts": [_artifact_record(p, relative_to=output) for p in artifacts],
        "metadata": {
            "adapter": ADAPTER_ID,
            "runtime_id": RUNTIME_ID,
            "fit_dataset_ids": list(request["fit_dataset_ids"]),
            **dict(extra_metadata),
        },
    }
    _write_json(output / "adapter_receipt.json", receipt)


def _atlas_path(request: Mapping[str, Any]) -> Path:
    role = "dataset_view_data:resource_atlas_geneformer_smoke_1000_v1"
    return _validate_file_artifact(_input_by_role(request, role), "Atlas smoke H5AD")


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


def _load_atlas(path: Path):
    import anndata

    adata = anndata.read_h5ad(path)
    if not {"donor_id", "broad_label"}.issubset(adata.obs.columns):
        raise ScientificAdapterError("Atlas H5AD lacks required observation fields")
    if "ensembl_id" not in adata.var.columns:
        raise ScientificAdapterError("Atlas H5AD lacks the Ensembl feature axis")
    if len(set(map(str, adata.obs_names))) != adata.n_obs:
        raise ScientificAdapterError("Atlas row identifiers are not unique")
    return adata


def prepare(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    parameters, roster, model_id, _ = _validate_request(request, "prepare")
    adata = _load_atlas(_atlas_path(request))
    observed = set(map(str, adata.obs["broad_label"]))
    if observed != set(roster):
        raise ScientificAdapterError("Atlas label set differs from the frozen roster")
    rows = _split_rows(adata, parameters, request["run_spec"])
    split_path = output / "split_rows.tsv"
    _write_tsv(
        split_path,
        ("row_hash", "unit_hash", "fold", "held_out"),
        (
            {**row, "held_out": "true" if row["held_out"] else "false"}
            for row in rows
        ),
    )

    contract = _load_contract()
    checkpoint = _checkpoint_path(request, model_id)
    registered_checkpoint = contract.get("checkpoint", {})
    checkpoint_artifact = _input_by_role(request, f"model_checkpoint_bundle:{model_id}")
    if (
        str(checkpoint_artifact.get("sha256")) != registered_checkpoint.get("sha256")
        or int(checkpoint_artifact.get("size_bytes", -1)) != registered_checkpoint.get("size_bytes")
    ):
        raise ScientificAdapterError(
            "supplied scPRINT checkpoint does not match the registered contract identity"
        )
    payload = _load_checkpoint(checkpoint)
    architecture = _inspect_architecture(payload, contract, model_id)

    gene_ids = list(map(str, adata.var["ensembl_id"]))
    checkpoint_genes = payload["hyper_parameters"]["genes"]
    checkpoint_index = {gene: index for index, gene in enumerate(checkpoint_genes)}
    alignment = [
        {"ensembl_id": gene, "checkpoint_gene_index": checkpoint_index[gene]}
        for gene in gene_ids
        if gene in checkpoint_index
    ]
    if not alignment:
        raise ScientificAdapterError(
            "no Atlas gene maps into the scPRINT checkpoint gene vocabulary"
        )
    gene_alignment_path = output / "gene_alignment.tsv"
    _write_tsv(
        gene_alignment_path,
        ("ensembl_id", "checkpoint_gene_index"),
        (
            {"ensembl_id": r["ensembl_id"], "checkpoint_gene_index": str(r["checkpoint_gene_index"])}
            for r in sorted(alignment, key=lambda r: r["ensembl_id"])
        ),
    )
    gene_summary = {
        "schema_version": "masld-bench-scprint-gene-alignment-v1",
        "atlas_genes_total": len(gene_ids),
        "checkpoint_vocabulary_size": len(checkpoint_genes),
        "genes_aligned": len(alignment),
        "atlas_gene_fraction_covered": len(alignment) / len(gene_ids),
        "checkpoint_vocabulary_fraction_present": len(alignment) / len(checkpoint_genes),
    }
    architecture_path = output / "checkpoint_inspection.json"
    _write_json(architecture_path, {**architecture, **{"gene_alignment_summary": gene_summary}})

    _receipt(
        action="prepare",
        request=request,
        output=output,
        artifacts=[split_path, gene_alignment_path, architecture_path],
        extra_metadata={
            "model_id": model_id,
            "genes_aligned": gene_summary["genes_aligned"],
            "atlas_gene_fraction_covered": gene_summary["atlas_gene_fraction_covered"],
            "architecture_contract_mismatches": architecture["architecture_contract_mismatches"],
        },
    )


def _read_split(path: Path) -> list[dict[str, Any]]:
    import csv

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [
        {
            "row_hash": r["row_hash"],
            "unit_hash": r["unit_hash"],
            "held_out": r["held_out"] == "true",
        }
        for r in rows
    ]


def _prior_output(request: Mapping[str, Any], action: str) -> Path:
    outputs = request.get("prior_action_outputs")
    if not isinstance(outputs, Mapping) or action not in outputs:
        raise ScientificAdapterError(f"request lacks the frozen {action} output")
    return Path(str(outputs[action])).resolve(strict=True)


def _require_scprint_encoder(model_id: str, checkpoint: Path) -> None:
    """Fail loudly and specifically rather than improvise the forward pass.

    Restricted inspection (see the module docstring) established that this
    checkpoint's gene vocabulary and trained weights are fully self-contained
    and do not need an external ESM/protein-embedding side file at inference
    time.  What frozen-encoder embedding extraction still needs is the
    ``scprint`` package itself: the exact ``nn.Module`` graph that consumes
    ``state_dict`` keys such as ``transformer.blocks.*.mixer.Wqkv`` (a fused
    flash-attention QKV projection), the gene/expression embedding fusion,
    and "cls"-style pooling.  That package is not importable in this
    environment and is not vendored anywhere on this cluster (checked under
    every micromamba/conda env visible from this job).  Reimplementing that
    forward pass from tensor shapes alone, without the source, would silently
    risk producing embeddings that are wrong in ways indistinguishable from
    correct ones -- exactly what this adapter must not do.
    """

    try:
        import scprint  # noqa: F401
    except ImportError as error:
        raise ScientificAdapterError(
            "scPRINT frozen-encoder embedding extraction requires the 'scprint' "
            "Python package (release 1.6.4, cantinilab/scPRINT), which is not "
            "installed in this environment and is not vendored anywhere on this "
            "cluster. The verified checkpoint at "
            f"{checkpoint} contains a self-sufficient gene vocabulary and "
            "trained weights (no external ESM side file is needed at inference "
            "time), but this adapter refuses to hand-reconstruct scPRINT's "
            "custom flash-attention transformer forward pass from state_dict "
            "shapes alone. Install 'scprint' (with matching scDataLoader>=1.6.4 "
            "and a Torch runtime compatible with its FlashAttention path) with "
            "explicit approval, then re-run this action; no partial mapping "
            f"was substituted. model_id={model_id!r}. Original error: {error}"
        ) from error


def fit(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    import torch

    parameters, roster, model_id, head_id = _validate_request(request, "fit")
    prepare_dir = _prior_output(request, "prepare")
    _read_split(prepare_dir / "split_rows.tsv")  # validated shape; reused by predict

    contract = _load_contract()
    checkpoint = _checkpoint_path(request, model_id)
    payload = _load_checkpoint(checkpoint)
    _inspect_architecture(payload, contract, model_id)

    device = _select_device(str(parameters.get("device", "auto")))
    seed = int(request["run_spec"]["seed"])
    torch.manual_seed(seed)

    # The checkpoint is verified and its architecture facts are frozen above.
    # Everything past this point -- building the scPrint nn.Module, running
    # the frozen forward pass to get per-cell CLS embeddings, and training the
    # donor/class-balanced common head identically to geneformer_common_lane
    # -- requires the scprint package itself.
    _require_scprint_encoder(model_id, checkpoint)
    raise ScientificAdapterError(  # pragma: no cover - unreachable while scprint is absent
        "unreachable: _require_scprint_encoder must raise when scprint is missing"
    )


def predict(request_path: Path, request: Mapping[str, Any], output: Path) -> None:
    parameters, roster, model_id, head_id = _validate_request(request, "predict")
    prepare_dir = _prior_output(request, "prepare")
    _prior_output(request, "fit")
    _read_split(prepare_dir / "split_rows.tsv")

    contract = _load_contract()
    checkpoint = _checkpoint_path(request, model_id)
    payload = _load_checkpoint(checkpoint)
    _inspect_architecture(payload, contract, model_id)

    _require_scprint_encoder(model_id, checkpoint)
    raise ScientificAdapterError(  # pragma: no cover - unreachable while scprint is absent
        "unreachable: _require_scprint_encoder must raise when scprint is missing"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--action", required=True, choices=("prepare", "fit", "predict"))
    parser.add_argument("--request", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args(argv)
    try:
        request = json.loads(arguments.request.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificAdapterError("adapter request is invalid JSON") from error
    if not isinstance(request, Mapping):
        raise ScientificAdapterError("adapter request must be an object")
    if arguments.output.exists():
        raise ScientificAdapterError(f"adapter output already exists: {arguments.output}")
    arguments.output.mkdir(parents=True, exist_ok=False)
    {"prepare": prepare, "fit": fit, "predict": predict}[arguments.action](
        arguments.request.resolve(strict=True), request, arguments.output
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
