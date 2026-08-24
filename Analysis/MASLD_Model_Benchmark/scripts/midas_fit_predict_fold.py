#!/usr/bin/env python3
"""Fit one donor-held MIDAS fold and translate held-donor RNA to ATAC."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import random
import sys
from typing import Any, Mapping, Sequence


MODEL_ID = "midas_inductive"


class MIDASFitError(ValueError):
    """Raised when the donor-held MIDAS execution contract differs."""


def load_helper(path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("midas_multivi_export_helper", path)
    if spec is None or spec.loader is None:
        raise MIDASFitError("prediction helper spec differs")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.MODEL_ID = MODEL_ID
    return module


def chromosome_chunks(peaks: Sequence[Mapping[str, str]]) -> list[int]:
    groups: list[tuple[str, int]] = []
    for row in peaks:
        chromosome = row["chromosome"]
        if not groups or groups[-1][0] != chromosome:
            groups.append((chromosome, 1))
        else:
            groups[-1] = (chromosome, groups[-1][1] + 1)
    if len(groups) < 2 or len({item[0] for item in groups}) != len(groups):
        raise MIDASFitError("peak chromosome blocks are not contiguous")
    counts = [item[1] for item in groups]
    if sum(counts) != len(peaks):
        raise MIDASFitError("peak chromosome census differs")
    return counts


def project_config() -> dict[str, Any]:
    from scmidas.config import load_config

    config = dict(load_config())
    config.update(
        {
            "num_workers": 0,
            "pin_memory": False,
            "persistent_workers": False,
        }
    )
    return config


def state_hash(state: Mapping[str, Any]) -> str:
    digest = sha256()
    for name in sorted(state):
        tensor = state[name].detach().cpu().contiguous()
        digest.update(name.encode())
        digest.update(str(tensor.dtype).encode())
        digest.update(json.dumps(list(tensor.shape), separators=(",", ":")).encode())
        digest.update(tensor.numpy().tobytes(order="C"))
    return digest.hexdigest()


def build_training_mudata(
    training_rna: Any,
    training_atac: Any,
    training_rows: Sequence[Mapping[str, str]],
    genes: Sequence[Mapping[str, str]],
    peaks: Sequence[Mapping[str, str]],
) -> Any:
    import anndata as ad
    import mudata as mu
    import pandas as pd

    cells = [row["cell_id"] for row in training_rows]
    obs = pd.DataFrame({"batch": ["paired_training"] * len(cells)}, index=cells)
    rna = ad.AnnData(
        X=training_rna,
        obs=obs.copy(),
        var=pd.DataFrame(index=[row["ensembl_id"] for row in genes]),
    )
    atac = ad.AnnData(
        X=training_atac,
        obs=obs.copy(),
        var=pd.DataFrame(index=[row["peak_id"] for row in peaks]),
    )
    return mu.MuData({"rna": rna, "atac": atac})


def query_dataset(
    query_rna: Any,
    query_rows: Sequence[Mapping[str, str]],
    genes: Sequence[Mapping[str, str]],
) -> Any:
    import anndata as ad
    import pandas as pd
    from scmidas.data import MultiModalDataset

    cells = [row["cell_id"] for row in query_rows]
    adata = ad.AnnData(
        X=query_rna,
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=[row["ensembl_id"] for row in genes]),
    )
    return MultiModalDataset(
        {"rna": adata},
        {"rna": 0, "joint": 0},
        {"rna": "anndata"},
        None,
        None,
    )


def predict_query(model: Any, dataset: Any, seed: int) -> Any:
    import numpy as np
    import torch

    model.datalist = [dataset]
    model.batch_names = ["rna_query"]
    model.combs = [["rna"]]
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    predicted = model.predict(
        return_in_memory=True,
        joint_latent=False,
        mod_latent=False,
        impute=False,
        batch_correct=False,
        translate=True,
        input=False,
        verbose=False,
    )
    try:
        values = np.asarray(
            predicted["rna_query"]["x_trans"]["rna_to_atac"], dtype=np.float64
        )
    except (KeyError, TypeError) as error:
        raise MIDASFitError("MIDAS RNA-to-ATAC output is absent") from error
    if (
        values.ndim != 2
        or len(values) != len(dataset)
        or np.any(~np.isfinite(values))
        or np.any(values < 0)
    ):
        raise MIDASFitError("MIDAS RNA-to-ATAC output differs")
    return values


def rewrite_bundle(path: Path) -> None:
    bundle = json.loads(path.read_text(encoding="utf-8"))
    metadata = bundle["metadata"]
    metadata.update(
        {
            "query_atac_tensor_created": False,
            "query_atac_tensor_nonzero_values": 0,
            "query_rna_used_for_training": False,
            "query_rna_used_after_fit_only": True,
            "training_datalist_excluded_query": True,
            "inductive_held_donor_inference": True,
            "exact_distribution": "scmidas 0.3.0",
            "rna_conditioned_atac_eligible": True,
            "sealed_inference_eligible": False,
        }
    )
    path.write_text(
        json.dumps(bundle, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )


def run(
    prepared: Path,
    helper_path: Path,
    output: Path,
    fold: int,
    seed: int,
    epochs: int,
) -> dict[str, Any]:
    import lightning as lightning
    import numpy as np
    from scipy import sparse
    import scmidas
    from scmidas import MIDAS
    import torch

    if output.exists() or fold not in range(5) or seed < 0 or epochs < 1:
        raise MIDASFitError("output, fold, seed, or epoch contract differs")
    if scmidas.__version__ != "0.3.0" or torch.__version__ != "2.6.0+cu124":
        raise MIDASFitError("MIDAS runtime version differs")
    if torch.cuda.device_count() != 1 or "L40S" not in torch.cuda.get_device_name(0):
        raise MIDASFitError("exactly one L40S is required")
    prepare_receipt = json.loads((prepared / "prepare_receipt.json").read_text())
    if (
        prepare_receipt.get("status") != "pass"
        or prepare_receipt.get("pairing_topology") != "same_nucleus"
        or prepare_receipt.get("held_atac_exported_from_prepare") is not False
        or prepare_receipt.get("fit_predict_source_hdf5_required") is not False
        or prepare_receipt.get("outcomes_read") is not False
    ):
        raise MIDASFitError("prepared same-nucleus contract differs")
    helper = load_helper(helper_path)
    root = prepared / f"fold_{fold}"
    training_rna = sparse.load_npz(root / "training_rna.npz").tocsr()
    training_atac = sparse.load_npz(root / "training_atac.npz").tocsr()
    query_rna = sparse.load_npz(root / "query_rna.npz").tocsr()
    train_fields, train_rows = helper.read_tsv(root / "training_rows.tsv")
    query_fields, query_rows = helper.read_tsv(root / "query_rows.tsv")
    gene_fields, genes = helper.read_tsv(root / "selected_genes.tsv")
    peak_fields, peaks = helper.read_tsv(root / "selected_peaks.tsv")
    if (
        train_fields != ("cell_id", "donor_id", "lineage", "rna_state", "atac_state")
        or query_fields != train_fields
        or gene_fields != ("selected_index", "source_index", "ensembl_id", "gene_name")
        or peak_fields
        != (
            "selected_index",
            "source_index",
            "peak_id",
            "chromosome",
            "bed_start",
            "bed_end",
        )
        or training_rna.shape != (len(train_rows), 2000)
        or training_atac.shape != (len(train_rows), 10_000)
        or query_rna.shape != (len(query_rows), 2000)
        or len(genes) != 2000
        or len(peaks) != 10_000
        or any(row["rna_state"] != "observed" for row in train_rows + query_rows)
        or any(row["atac_state"] != "observed" for row in train_rows)
        or any(row["atac_state"] != "structurally_missing" for row in query_rows)
    ):
        raise MIDASFitError("prepared fold axes or missing states differ")
    train_donors = {row["donor_id"] for row in train_rows}
    query_donors = {row["donor_id"] for row in query_rows}
    if train_donors & query_donors:
        raise MIDASFitError("donor crosses fit and query")
    if any(value <= 0 for value in np.asarray(training_rna.sum(axis=1)).ravel()):
        raise MIDASFitError("training RNA contains an empty nucleus")
    if any(value <= 0 for value in np.asarray(training_atac.sum(axis=1)).ravel()):
        raise MIDASFitError("training ATAC contains an empty nucleus")
    if any(value <= 0 for value in np.asarray(query_rna.sum(axis=1)).ravel()):
        raise MIDASFitError("query RNA contains an empty nucleus")
    chunks = chromosome_chunks(peaks)
    mdata = build_training_mudata(training_rna, training_atac, train_rows, genes, peaks)
    MIDAS.setup_mudata(
        mdata,
        batch_key="batch",
        dims_x={"rna": [len(genes)], "atac": chunks},
    )
    output.mkdir(parents=True, mode=0o750)
    lightning.seed_everything(seed, workers=True)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True)
    config = project_config()
    model = MIDAS(
        mdata,
        configs=config,
        batch_size=128,
        n_save=100_000,
        save_model_path=str(output / "unused_checkpoints"),
    )
    if len(model.datalist) != 1 or model.combs != [["rna", "atac"]]:
        raise MIDASFitError("training datalist is not paired-only")
    model.train(
        max_epochs=epochs,
        accelerator="gpu",
        devices=1,
        logger=False,
        enable_checkpointing=False,
        enable_progress_bar=False,
        num_sanity_val_steps=0,
        deterministic=True,
    )
    state = {name: value.detach().cpu() for name, value in model.state_dict().items()}
    state_digest = state_hash(state)
    state_path = output / "state_dict.pt"
    torch.save(state, state_path)
    restored = MIDAS(
        mdata,
        configs=config,
        batch_size=128,
        n_save=100_000,
        save_model_path=str(output / "unused_checkpoints_resumed"),
    ).cuda()
    restored.load_state_dict(
        torch.load(state_path, map_location="cuda", weights_only=True), strict=True
    )
    if state_hash(restored.state_dict()) != state_digest:
        raise MIDASFitError("strict MIDAS tensor-state resume differs")
    held_query = query_dataset(query_rna, query_rows, genes)
    cell_profiles = predict_query(restored, held_query, seed + 1)
    repeated = predict_query(restored, held_query, seed + 1)
    if not np.array_equal(cell_profiles, repeated):
        raise MIDASFitError("repeated MIDAS query translation differs")
    cell_profiles = helper.normalize_profiles(cell_profiles)
    fold_output = output / f"fold_{fold}"
    fold_output.mkdir(mode=0o750)
    run_id = helper.canonical_hash(
        {
            "model_id": MODEL_ID,
            "fold": fold,
            "seed": seed,
            "epochs": epochs,
            "prepared_receipt_sha256": helper.sha256_file(
                prepared / "prepare_receipt.json"
            ),
            "state_dict_sha256": state_digest,
        }
    )
    bundle_path = helper.export_bundle(
        fold=fold,
        run_id=run_id,
        query_rows=query_rows,
        peaks=peaks,
        cell_profiles=cell_profiles,
        output=fold_output,
        seed=seed,
        training_nuclei=len(train_rows),
        state_hash=state_digest,
    )
    rewrite_bundle(bundle_path)
    fold_receipt = {
        "run_id": run_id,
        "training_nuclei": len(train_rows),
        "query_nuclei": len(query_rows),
        "training_donors": len(train_donors),
        "query_donors": len(query_donors),
        "epochs": epochs,
        "state_dict_sha256": state_digest,
        "state_dict_file_sha256": helper.sha256_file(state_path),
        "prediction_bundle": bundle_path.relative_to(output).as_posix(),
        "prediction_bundle_sha256": helper.sha256_file(bundle_path),
        "strict_state_resume": True,
        "repeated_query_prediction_bit_identical": True,
        "query_rna_used_for_training": False,
        "held_atac_read": False,
    }
    receipt = {
        "schema_version": "masld-bench-midas-inductive-fit-predict-v1",
        "status": "pass",
        "model_id": MODEL_ID,
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "exact_distribution": f"scmidas {scmidas.__version__}",
        "pairing_topology": "same_nucleus_training",
        "outer_unit": "donor",
        "folds": {str(fold): fold_receipt},
        "base_seed": seed,
        "epochs": epochs,
        "atac_chromosome_chunks": chunks,
        "official_lightning_trainer_used": True,
        "training_datalist_paired_only": True,
        "query_rna_used_for_training": False,
        "query_rna_used_after_fit_only": True,
        "query_atac_state": "structurally_missing",
        "held_atac_input_exposed": False,
        "source_hdf5_available_to_fit_predict": False,
        "outcomes_read": False,
        "metrics_calculated_by_adapter": False,
        "test_outcomes_read": False,
        "sealed_inference_eligible": False,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--helper-path", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fold", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--epochs", type=int, default=30)
    arguments = parser.parse_args()
    run(
        arguments.prepared,
        arguments.helper_path,
        arguments.output,
        arguments.fold,
        arguments.seed,
        arguments.epochs,
    )


if __name__ == "__main__":
    main()
