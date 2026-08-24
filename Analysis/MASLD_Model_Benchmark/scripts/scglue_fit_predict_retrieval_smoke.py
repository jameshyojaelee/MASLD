#!/usr/bin/env python3
"""Fit scGLUE/PairedSCGLUE on training donors and encode blinded held modalities."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence


MODELS = ("scglue", "paired_scglue", "linear_cca")


class SCGLUERetrievalError(ValueError):
    """Raised when project scGLUE retrieval execution differs."""


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SCGLUERetrievalError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def array_hash(value: Any) -> str:
    import numpy as np

    array = np.ascontiguousarray(value)
    digest = sha256()
    digest.update(str(array.dtype).encode())
    digest.update(json.dumps(list(array.shape)).encode())
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def state_hash(state: Mapping[str, Any]) -> str:
    digest = sha256()
    for name in sorted(state):
        tensor = state[name].detach().contiguous().cpu()
        digest.update(name.encode())
        digest.update(array_hash(tensor.numpy()).encode())
    return digest.hexdigest()


def tfidf_representations(train: Any, query: Any, seed: int) -> tuple[Any, Any]:
    import numpy as np
    from sklearn.decomposition import TruncatedSVD

    train_depth = np.asarray(train.sum(axis=1)).ravel()
    query_depth = np.asarray(query.sum(axis=1)).ravel()
    if np.any(train_depth <= 0) or np.any(query_depth <= 0):
        raise SCGLUERetrievalError("ATAC representation contains empty rows")
    detection = np.asarray((train > 0).sum(axis=0)).ravel()
    inverse = np.log1p(train.shape[0] / np.maximum(detection, 1))
    train_tfidf = train.multiply((10_000.0 / train_depth)[:, None]).multiply(inverse)
    query_tfidf = query.multiply((10_000.0 / query_depth)[:, None]).multiply(inverse)
    svd = TruncatedSVD(n_components=50, random_state=seed, n_iter=7)
    return svd.fit_transform(train_tfidf), svd.transform(query_tfidf)


def rna_representations(train: Any, query: Any, seed: int) -> tuple[Any, Any]:
    import numpy as np
    from sklearn.decomposition import TruncatedSVD

    train_depth = np.asarray(train.sum(axis=1)).ravel()
    query_depth = np.asarray(query.sum(axis=1)).ravel()
    if np.any(train_depth <= 0) or np.any(query_depth <= 0):
        raise SCGLUERetrievalError("RNA representation contains empty rows")
    train_norm = train.multiply((10_000.0 / train_depth)[:, None]).tocsr()
    query_norm = query.multiply((10_000.0 / query_depth)[:, None]).tocsr()
    train_norm.data = np.log1p(train_norm.data)
    query_norm.data = np.log1p(query_norm.data)
    svd = TruncatedSVD(n_components=50, random_state=seed, n_iter=7)
    return svd.fit_transform(train_norm), svd.transform(query_norm)


def make_anndata(
    matrix: Any,
    identifiers: Sequence[str],
    features: Sequence[str],
    representation: Any,
    representation_key: str,
) -> Any:
    import anndata
    import pandas as pd

    if matrix.shape != (len(identifiers), len(features)):
        raise SCGLUERetrievalError("AnnData axes differ")
    obs = pd.DataFrame(
        {"batch": pd.Categorical(["batch0"] * len(identifiers))},
        index=pd.Index(identifiers),
    )
    var = pd.DataFrame(
        {"highly_variable": True},
        index=pd.Index(features),
    )
    adata = anndata.AnnData(X=matrix, obs=obs, var=var)
    adata.obsm[representation_key] = representation
    return adata


def configure(scglue: Any, rna: Any, atac: Any) -> None:
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


def read_graph(path: Path, vertices: set[str]) -> Any:
    import networkx as nx

    fields, rows = read_tsv(path)
    if fields != ("source", "target", "weight", "sign"):
        raise SCGLUERetrievalError("guidance edge schema differs")
    graph = nx.MultiDiGraph()
    for row in rows:
        source, target = row["source"], row["target"]
        weight, sign = float(row["weight"]), float(row["sign"])
        if (
            source not in vertices
            or target not in vertices
            or not math.isfinite(weight)
            or weight <= 0
            or sign != 1.0
        ):
            raise SCGLUERetrievalError("guidance edge differs")
        graph.add_edge(source, target, weight=weight, sign=sign)
    if set(graph.nodes) != vertices or graph.number_of_edges() <= len(vertices):
        raise SCGLUERetrievalError("guidance graph vertex or edge census differs")
    return graph


def save_embeddings(
    path: Path,
    rna_ids: Sequence[str],
    rna_embedding: Any,
    atac_ids: Sequence[str],
    atac_embedding: Any,
) -> dict[str, Any]:
    import numpy as np

    rna = np.asarray(rna_embedding, dtype=np.float32)
    atac = np.asarray(atac_embedding, dtype=np.float32)
    if (
        rna.shape != (len(rna_ids), 16)
        or atac.shape != (len(atac_ids), 16)
        or np.any(~np.isfinite(rna))
        or np.any(~np.isfinite(atac))
    ):
        raise SCGLUERetrievalError("query embedding axes or values differ")
    with path.open("xb") as handle:
        np.savez_compressed(
            handle,
            rna_ids=np.asarray(rna_ids),
            rna_embedding=rna,
            atac_ids=np.asarray(atac_ids),
            atac_embedding=atac,
        )
    return {
        "path": path.name,
        "rna_embedding_sha256": array_hash(rna),
        "atac_embedding_sha256": array_hash(atac),
        "rna_rows": len(rna_ids),
        "atac_rows": len(atac_ids),
        "latent_dimension": 16,
    }


def fit_scglue_model(
    *,
    scglue: Any,
    torch: Any,
    model_name: str,
    train_rna: Any,
    train_atac: Any,
    query_rna: Any,
    query_atac: Any,
    graph: Any,
    output: Path,
    seed: int,
) -> dict[str, Any]:
    model_class = (
        scglue.models.SCGLUEModel
        if model_name == "scglue"
        else scglue.models.PairedSCGLUEModel
    )
    vertices = sorted(graph.nodes)
    model = model_class(
        {"rna": train_rna, "atac": train_atac},
        vertices,
        latent_dim=16,
        h_depth=1,
        h_dim=128,
        dropout=0.1,
        random_seed=seed,
    )
    if model_name == "paired_scglue":
        model.compile(lam_joint_cross=0.02, lam_real_cross=0.02, lam_cos=0.02)
    else:
        model.compile()
    training = output / "training"
    model.fit(
        {"rna": train_rna, "atac": train_atac},
        graph,
        data_batch_size=128,
        graph_batch_size=4096,
        align_burnin=2,
        safe_burnin=False,
        max_epochs=30,
        patience=None,
        reduce_lr_patience=None,
        directory=training,
    )
    state = {name: value.detach().cpu().clone() for name, value in model.net.state_dict().items()}
    state_sha = state_hash(state)
    state_path = output / "state_dict.pt"
    torch.save(state, state_path)
    resumed = model_class(
        {"rna": train_rna, "atac": train_atac},
        vertices,
        latent_dim=16,
        h_depth=1,
        h_dim=128,
        dropout=0.1,
        random_seed=seed,
    )
    resumed.net.load_state_dict(
        torch.load(state_path, map_location="cpu", weights_only=True), strict=True
    )
    first_rna = resumed.encode_data("rna", query_rna, batch_size=128)
    first_atac = resumed.encode_data("atac", query_atac, batch_size=128)
    second_rna = resumed.encode_data("rna", query_rna, batch_size=128)
    second_atac = resumed.encode_data("atac", query_atac, batch_size=128)
    if not __import__("numpy").array_equal(first_rna, second_rna) or not __import__(
        "numpy"
    ).array_equal(first_atac, second_atac):
        raise SCGLUERetrievalError(f"{model_name} repeated query encoding differs")
    embedding = save_embeddings(
        output / "query_embeddings.npz",
        query_rna.obs_names.tolist(),
        first_rna,
        query_atac.obs_names.tolist(),
        first_atac,
    )
    return {
        **embedding,
        "model_class": model_class.__name__,
        "state_dict_sha256": state_sha,
        "state_key_count": len(state),
        "fit_epochs": 30,
        "strict_state_resume": True,
        "repeat_query_encoding_bit_identical": True,
        "hidden_pair_map_read": False,
        "retrieval_metrics_calculated": False,
    }


def run(inputs: Path, output: Path, seed: int) -> dict[str, Any]:
    import numpy as np
    from scipy import sparse
    from sklearn.cross_decomposition import CCA
    import scglue
    import torch

    if output.exists() or seed < 0 or scglue.__version__ != "0.4.1":
        raise SCGLUERetrievalError("output, seed, or scGLUE version differs")
    if torch.__version__ != "2.6.0+cu124":
        raise SCGLUERetrievalError("torch version differs")
    if torch.cuda.device_count() != 1 or "L40S" not in torch.cuda.get_device_name(0):
        raise SCGLUERetrievalError("exactly one L40S is required")
    scglue.config.CPU_ONLY = False
    scglue.config.DATALOADER_NUM_WORKERS = 0
    scglue.config.DATALOADER_PIN_MEMORY = False
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    output.mkdir(parents=True, mode=0o750)
    fold_receipts: dict[str, Any] = {}
    for fold in range(5):
        fold_seed = seed + fold
        root = inputs / f"fold_{fold}"
        train_rna_matrix = sparse.load_npz(root / "training_rna.npz").tocsr()
        train_atac_matrix = sparse.load_npz(root / "training_atac.npz").tocsr()
        query_rna_matrix = sparse.load_npz(root / "query_rna.npz").tocsr()
        query_atac_matrix = sparse.load_npz(root / "query_atac.npz").tocsr()
        train_fields, train_rows = read_tsv(root / "training_rows.tsv")
        rna_fields, rna_rows = read_tsv(root / "query_rna_rows.tsv")
        atac_fields, atac_rows = read_tsv(root / "query_atac_rows.tsv")
        _gene_fields, genes = read_tsv(root / "genes.tsv")
        _peak_fields, peaks = read_tsv(root / "peaks.tsv")
        if (
            train_fields != ("training_id", "rna_state", "atac_state")
            or rna_fields != ("rna_query_id", "rna_state")
            or atac_fields != ("atac_query_id", "atac_state")
            or len(train_rows) != train_rna_matrix.shape[0]
            or train_atac_matrix.shape != (len(train_rows), 10_000)
            or query_rna_matrix.shape != (len(rna_rows), 2000)
            or query_atac_matrix.shape != (len(atac_rows), 10_000)
            or len({row["rna_query_id"] for row in rna_rows}) != len(rna_rows)
            or len({row["atac_query_id"] for row in atac_rows}) != len(atac_rows)
        ):
            raise SCGLUERetrievalError("fold matrix, row, or identity axes differ")
        train_rna_rep, query_rna_rep = rna_representations(
            train_rna_matrix, query_rna_matrix, fold_seed
        )
        train_atac_rep, query_atac_rep = tfidf_representations(
            train_atac_matrix, query_atac_matrix, fold_seed
        )
        train_ids = [row["training_id"] for row in train_rows]
        rna_ids = [row["rna_query_id"] for row in rna_rows]
        atac_ids = [row["atac_query_id"] for row in atac_rows]
        gene_ids = [row["gene_id"] for row in genes]
        peak_ids = [row["peak_id"] for row in peaks]
        train_rna = make_anndata(
            train_rna_matrix, train_ids, gene_ids, train_rna_rep, "X_pca"
        )
        train_atac = make_anndata(
            train_atac_matrix, train_ids, peak_ids, train_atac_rep, "X_lsi"
        )
        query_rna = make_anndata(
            query_rna_matrix, rna_ids, gene_ids, query_rna_rep, "X_pca"
        )
        query_atac = make_anndata(
            query_atac_matrix, atac_ids, peak_ids, query_atac_rep, "X_lsi"
        )
        for pair in ((train_rna, train_atac), (query_rna, query_atac)):
            configure(scglue, *pair)
        graph = read_graph(root / "guidance_edges.tsv", set(gene_ids + peak_ids))
        fold_output = output / f"fold_{fold}"
        fold_output.mkdir(mode=0o750)
        models: dict[str, Any] = {}
        for model_name in ("scglue", "paired_scglue"):
            model_output = fold_output / model_name
            model_output.mkdir(mode=0o750)
            models[model_name] = fit_scglue_model(
                scglue=scglue,
                torch=torch,
                model_name=model_name,
                train_rna=train_rna,
                train_atac=train_atac,
                query_rna=query_rna,
                query_atac=query_atac,
                graph=graph,
                output=model_output,
                seed=fold_seed,
            )
        cca = CCA(n_components=16, max_iter=1000, tol=1.0e-6, scale=True)
        cca.fit(train_rna_rep, train_atac_rep)
        query_rna_cca, query_atac_cca = cca.transform(query_rna_rep, query_atac_rep)
        baseline_output = fold_output / "linear_cca"
        baseline_output.mkdir(mode=0o750)
        models["linear_cca"] = {
            **save_embeddings(
                baseline_output / "query_embeddings.npz",
                rna_ids,
                query_rna_cca,
                atac_ids,
                query_atac_cca,
            ),
            "fit_on_training_pairs_only": True,
            "hidden_pair_map_read": False,
            "retrieval_metrics_calculated": False,
        }
        fold_receipts[str(fold)] = {
            "training_nuclei": len(train_rows),
            "query_rna_nuclei": len(rna_rows),
            "query_atac_nuclei": len(atac_rows),
            "guidance_vertices": graph.number_of_nodes(),
            "guidance_edges": graph.number_of_edges(),
            "models": models,
        }
    receipt = {
        "schema_version": "masld-bench-scglue-retrieval-fit-predict-v1",
        "status": "pass",
        "model_roster": list(MODELS),
        "dataset_id": "gse296875",
        "pairing_topology": "same_nucleus_training_blinded_query",
        "outer_unit": "donor",
        "folds": fold_receipts,
        "scglue": scglue.__version__,
        "torch": torch.__version__,
        "gpu": torch.cuda.get_device_name(0),
        "hidden_pair_map_read": False,
        "query_modality_ids_disjoint": True,
        "retrieval_metrics_calculated": False,
        "rna_to_atac_prediction_executed": False,
        "sealed_rna_conditioned_atac_eligible": False,
        "smoke_only": True,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=1103)
    arguments = parser.parse_args()
    run(arguments.inputs, arguments.output, arguments.seed)


if __name__ == "__main__":
    main()
