#!/usr/bin/env python3
"""Synthetic L40S probe for exact scGLUE 0.4.1 latent integration."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


class SCGLUEProbeError(ValueError):
    """Raised when the synthetic scGLUE requirement differs."""


def _array_sha256(value: Any) -> str:
    import numpy as np

    array = np.ascontiguousarray(value)
    digest = sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(json.dumps(list(array.shape)).encode("ascii"))
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _state_sha256(state: dict[str, Any]) -> str:
    digest = sha256()
    for name in sorted(state):
        value = state[name].detach().contiguous().cpu()
        digest.update(name.encode("utf-8"))
        digest.update(_array_sha256(value.numpy()).encode("ascii"))
    return digest.hexdigest()


def _device_contract(torch: Any) -> str:
    if torch.cuda.device_count() != 1:
        raise SCGLUEProbeError("exactly one GPU is required")
    name = torch.cuda.get_device_name(0)
    if "L40S" not in name:
        raise SCGLUEProbeError("exactly one L40S is required")
    return name


def _make_inputs() -> tuple[Any, Any, Any]:
    import anndata
    import networkx as nx
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA, TruncatedSVD

    rng = np.random.default_rng(20260824)
    n_cells, n_genes, n_peaks = 48, 16, 20
    rna_counts = rng.poisson(1.6, size=(n_cells, n_genes)).astype(np.float32)
    atac_counts = rng.poisson(0.35, size=(n_cells, n_peaks)).astype(np.float32)
    rna_counts[:, 0] += 1
    atac_counts[:, 0] += 1
    cell_ids = pd.Index([f"SAME_NUCLEUS_{index:03d}" for index in range(n_cells)])
    batches = pd.Categorical(["batch_a"] * 24 + ["batch_b"] * 24)
    obs = pd.DataFrame({"batch": batches}, index=cell_ids)
    genes = pd.Index([f"gene_{index:03d}" for index in range(n_genes)])
    peaks = pd.Index([f"peak_{index:03d}" for index in range(n_peaks)])
    rna = anndata.AnnData(
        X=rna_counts,
        obs=obs.copy(),
        var=pd.DataFrame({"highly_variable": True}, index=genes),
    )
    atac = anndata.AnnData(
        X=atac_counts,
        obs=obs.copy(),
        var=pd.DataFrame({"highly_variable": True}, index=peaks),
    )
    rna.obsm["X_pca"] = PCA(n_components=8, random_state=20260824).fit_transform(
        np.log1p(rna_counts)
    )
    atac.obsm["X_lsi"] = TruncatedSVD(
        n_components=8, random_state=20260824
    ).fit_transform(np.log1p(atac_counts))

    graph = nx.MultiDiGraph()
    for feature in genes.append(peaks):
        graph.add_edge(feature, feature, weight=1.0, sign=1.0)
    for peak_index, peak in enumerate(peaks):
        gene = genes[peak_index % n_genes]
        graph.add_edge(gene, peak, weight=0.5, sign=1.0)
        graph.add_edge(peak, gene, weight=0.5, sign=1.0)
    return rna, atac, graph


def _configure(scglue: Any, rna: Any, atac: Any) -> None:
    scglue.models.configure_dataset(
        rna,
        "NB",
        use_rep="X_pca",
        use_highly_variable=True,
        use_batch="batch",
        use_obs_names=True,
    )
    scglue.models.configure_dataset(
        atac,
        "NB",
        use_rep="X_lsi",
        use_highly_variable=True,
        use_batch="batch",
        use_obs_names=True,
    )


def _run_model(
    scglue: Any,
    torch: Any,
    model_class: Any,
    model_name: str,
    rna: Any,
    atac: Any,
    graph: Any,
    workdir: Path,
) -> dict[str, Any]:
    vertices = sorted(graph.nodes)
    model = model_class(
        {"rna": rna, "atac": atac},
        vertices,
        latent_dim=4,
        h_depth=1,
        h_dim=16,
        dropout=0.0,
        random_seed=20260824,
    )
    if model_name == "paired_scglue":
        model.compile(lam_joint_cross=0.02, lam_real_cross=0.02, lam_cos=0.02)
    else:
        model.compile()
    model.fit(
        {"rna": rna, "atac": atac},
        graph,
        data_batch_size=16,
        graph_batch_size=64,
        align_burnin=1,
        safe_burnin=False,
        max_epochs=3,
        patience=None,
        reduce_lr_patience=None,
        directory=workdir,
    )
    rna_embedding = model.encode_data("rna", rna, batch_size=16)
    atac_embedding = model.encode_data("atac", atac, batch_size=16)
    graph_embedding = model.encode_graph(graph)
    if (
        rna_embedding.shape != (48, 4)
        or atac_embedding.shape != (48, 4)
        or graph_embedding.shape != (36, 4)
    ):
        raise SCGLUEProbeError(f"{model_name} embedding axes differ")

    query = rna[:4].copy()
    expanded = rna[:7].copy()
    query_embedding = model.encode_data("rna", query, batch_size=16)
    expanded_embedding = model.encode_data("rna", expanded, batch_size=16)[:4]
    composition_max_abs_difference = float(
        abs(query_embedding - expanded_embedding).max()
    )
    if composition_max_abs_difference > 1.0e-7:
        raise SCGLUEProbeError(f"{model_name} query composition changes embeddings")

    state = {
        name: value.detach().cpu().clone()
        for name, value in model.net.state_dict().items()
    }
    resumed = model_class(
        {"rna": rna, "atac": atac},
        vertices,
        latent_dim=4,
        h_depth=1,
        h_dim=16,
        dropout=0.0,
        random_seed=20260824,
    )
    resumed.net.load_state_dict(state, strict=True)
    resumed_embedding = resumed.encode_data("rna", query, batch_size=16)
    resume_max_abs_difference = float(abs(query_embedding - resumed_embedding).max())
    if resume_max_abs_difference != 0.0:
        raise SCGLUEProbeError(f"{model_name} state-dict resume differs")
    if not all(
        torch.isfinite(value).all().item() for value in model.net.state_dict().values()
    ):
        raise SCGLUEProbeError(f"{model_name} contains non-finite state")
    return {
        "model_class": model_class.__name__,
        "training_rows_per_modality": 48,
        "latent_dimension": 4,
        "rna_embedding_sha256": _array_sha256(rna_embedding),
        "atac_embedding_sha256": _array_sha256(atac_embedding),
        "graph_embedding_sha256": _array_sha256(graph_embedding),
        "state_dict_sha256": _state_sha256(state),
        "state_key_count": len(state),
        "query_composition_max_abs_difference": composition_max_abs_difference,
        "query_composition_tolerance": 1.0e-7,
        "state_dict_resume_max_abs_difference": resume_max_abs_difference,
        "fit_executed": True,
        "rna_encode_executed": True,
        "atac_encode_executed": True,
        "graph_encode_executed": True,
        "state_dict_resume_executed": True,
        "external_dill_deserialization_executed": False,
        "decode_data_executed": False,
    }


def run(output: Path) -> dict[str, Any]:
    import scglue
    import torch

    if output.exists():
        raise SCGLUEProbeError("output already exists")
    if scglue.__version__ != "0.4.1":
        raise SCGLUEProbeError("scGLUE version differs")
    torch.manual_seed(20260824)
    torch.cuda.manual_seed_all(20260824)
    gpu = _device_contract(torch)
    scglue.config.CPU_ONLY = False
    scglue.config.DATALOADER_NUM_WORKERS = 0
    scglue.config.DATALOADER_PIN_MEMORY = False
    rna, atac, graph = _make_inputs()
    if not rna.obs_names.equals(atac.obs_names):
        raise SCGLUEProbeError("same-nucleus row identity differs")
    _configure(scglue, rna, atac)
    output.mkdir(parents=True, mode=0o750)
    models = {
        "scglue": _run_model(
            scglue,
            torch,
            scglue.models.SCGLUEModel,
            "scglue",
            rna,
            atac,
            graph,
            output / "training" / "scglue",
        ),
        "paired_scglue": _run_model(
            scglue,
            torch,
            scglue.models.PairedSCGLUEModel,
            "paired_scglue",
            rna,
            atac,
            graph,
            output / "training" / "paired_scglue",
        ),
    }
    receipt = {
        "schema_version": "masld-bench-scglue-no-outcome-probe-v1",
        "status": "pass",
        "scglue": scglue.__version__,
        "torch": torch.__version__,
        "torch_cuda_build": torch.version.cuda,
        "gpu": gpu,
        "models": models,
        "synthetic_inputs": True,
        "pairing_topology": "same_nucleus",
        "ordered_same_nucleus_ids_verified": True,
        "gse244832_false_pairing": False,
        "gse281367_false_pairing": False,
        "project_data_read": False,
        "outcomes_read": False,
        "labels_read": False,
        "latent_integration_only": True,
        "rna_to_atac_prediction_executed": False,
        "sealed_rna_conditioned_atac_eligible": False,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    run(arguments.output)


if __name__ == "__main__":
    main()
