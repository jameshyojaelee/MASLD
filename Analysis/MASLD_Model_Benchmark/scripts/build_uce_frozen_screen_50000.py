#!/usr/bin/env python3
"""Build deterministic, label-blind UCE sentences for the 50k Atlas screen."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
from typing import Any, Mapping

from scripts.uce_frozen_screen_50000_common import (
    ARCH,
    MODEL_IDS,
    PE_EMBEDDING_ROWS,
    UCEActivationError,
    inspect_state_dict,
    load_activation_contract,
    load_restricted_plain_mapping,
    sha256_file,
    tokenize_cell,
    verify_record,
)


ALLOWED_COUNT_OBS = {"row_id", "donor_id", "dataset_id", "outer_fold", "assay"}
FORBIDDEN_OUTCOMES = {
    "broad_label",
    "source_cell_type",
    "cell_type",
    "label",
    "fibrosis",
    "nas",
    "mash",
    "masld",
    "sex",
    "age",
    "stage",
}


def _strings(values: Any) -> list[str]:
    result: list[str] = []
    for value in values:
        result.append(value.decode("utf-8") if isinstance(value, bytes) else str(value))
    return result


def _read_h5_column(group: Any, key: str) -> list[str]:
    import h5py
    import numpy as np

    if key not in group:
        raise UCEActivationError(f"H5AD source lacks required field {key}")
    node = group[key]
    if isinstance(node, h5py.Dataset):
        return _strings(node[:])
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        categories = _strings(node["categories"][:])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        if np.any(codes < 0) or np.any(codes >= len(categories)):
            raise UCEActivationError(f"H5AD categorical field {key} has invalid codes")
        return [categories[int(code)] for code in codes]
    raise UCEActivationError(f"H5AD field {key} has an unsupported encoding")


def read_symbol_axis(path: Path) -> tuple[list[str], list[str]]:
    """Read only feature identifiers from the label-bearing source H5AD."""

    import h5py

    with h5py.File(path, "r") as handle:
        if "var" not in handle:
            raise UCEActivationError("symbol-axis H5AD lacks var")
        var = handle["var"]
        ensembl = _read_h5_column(var, "ensembl_id")
        symbols = _read_h5_column(var, "source_feature_id")
    if len(ensembl) != len(symbols) or not ensembl:
        raise UCEActivationError("symbol-axis feature arrays differ")
    return ensembl, symbols


def read_split(path: Path) -> dict[str, tuple[str, str, int]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["row_id", "donor_id", "dataset", "outer_fold"]:
            raise UCEActivationError("study-held-out split schema differs")
        result: dict[str, tuple[str, str, int]] = {}
        for row in reader:
            row_id = row["row_id"]
            if row_id in result:
                raise UCEActivationError("study-held-out row identifiers collide")
            fold = int(row["outer_fold"])
            if fold not in range(5):
                raise UCEActivationError("study-held-out fold is outside 0..4")
            result[row_id] = (row["donor_id"], row["dataset"], fold)
    return result


def validate_firewall(
    row_ids: list[str],
    donors: list[str],
    studies: list[str],
    split: Mapping[str, tuple[str, str, int]],
    *,
    expected_rows: int,
    expected_donors: int,
    expected_studies: int,
) -> list[int]:
    if (
        len(row_ids) != expected_rows
        or len(set(row_ids)) != expected_rows
        or len(row_ids) != len(donors)
        or len(row_ids) != len(studies)
        or len(set(donors)) != expected_donors
        or len(set(studies)) != expected_studies
        or set(split) != set(row_ids)
    ):
        raise UCEActivationError("row, donor, study, or split cardinality differs")
    donor_folds: dict[str, int] = {}
    study_folds: dict[str, int] = {}
    folds: list[int] = []
    for row_id, donor, study in zip(row_ids, donors, studies, strict=True):
        split_donor, split_study, fold = split[row_id]
        if (donor, study) != (split_donor, split_study):
            raise UCEActivationError("count and study-split metadata join differs")
        if donor_folds.setdefault(donor, fold) != fold:
            raise UCEActivationError("one donor crosses outer folds")
        if study_folds.setdefault(study, fold) != fold:
            raise UCEActivationError("one study crosses outer folds")
        folds.append(fold)
    if set(folds) != set(range(5)):
        raise UCEActivationError("study-held-out split lacks an outer fold")
    for held in range(5):
        train = {study for study, fold in study_folds.items() if fold != held}
        test = {study for study, fold in study_folds.items() if fold == held}
        if not test or train & test:
            raise UCEActivationError("outer train/test studies are not disjoint")
    return folds


def build_gene_universe(
    gene_symbols: list[str],
    gene_passes_min_cells: Any,
    protein_embeddings_path: Path,
    species_offsets_path: Path,
    species_chrom_path: Path,
) -> tuple[Any, Any, Any, Any, dict[str, Any]]:
    import numpy as np
    import pandas as pd
    import torch

    offsets = load_restricted_plain_mapping(species_offsets_path)
    human_offset = offsets.get("human")
    if not isinstance(human_offset, int):
        raise UCEActivationError("UCE species offsets lack an integer human offset")
    table = pd.read_csv(species_chrom_path)
    if not {"gene_symbol", "chromosome", "start", "species"}.issubset(table.columns):
        raise UCEActivationError("UCE species chromosome table schema differs")
    categories = pd.Categorical(
        table["species"].astype(str) + "_" + table["chromosome"].astype(str)
    )
    human_mask = table["species"] == "human"
    human_symbols = table.loc[human_mask, "gene_symbol"].astype(str).str.upper()
    if human_symbols.duplicated().any():
        raise UCEActivationError("UCE chromosome table has duplicate human gene symbols")
    chrom_by_symbol = dict(
        zip(human_symbols, categories.codes[human_mask.to_numpy()], strict=True)
    )
    start_by_symbol = dict(
        zip(human_symbols, table.loc[human_mask, "start"].to_numpy(), strict=True)
    )
    protein = torch.load(
        protein_embeddings_path, map_location="cpu", weights_only=True, mmap=True
    )
    if not isinstance(protein, Mapping) or not protein:
        raise UCEActivationError("UCE human protein embedding mapping differs")
    ordered_symbols = [str(symbol).upper() for symbol in protein]
    if len(set(ordered_symbols)) != len(ordered_symbols):
        raise UCEActivationError("UCE protein symbols collide under uppercasing")
    protein_index = {symbol: index for index, symbol in enumerate(ordered_symbols)}

    symbols = [str(symbol).upper() for symbol in gene_symbols]
    token_rows = np.full(len(symbols), -1, dtype=np.int64)
    chromosome_codes = np.full(len(symbols), -1, dtype=np.int64)
    starts = np.zeros(len(symbols), dtype=np.int64)
    usable = np.asarray(gene_passes_min_cells, dtype=bool).copy()
    for index, symbol in enumerate(symbols):
        if not usable[index] or symbol not in protein_index or symbol not in chrom_by_symbol:
            usable[index] = False
            continue
        token_rows[index] = protein_index[symbol] + human_offset
        chromosome_codes[index] = int(chrom_by_symbol[symbol])
        starts[index] = int(start_by_symbol[symbol])
    if not usable.any():
        raise UCEActivationError("no filtered Atlas gene maps to the UCE human universe")
    if (
        np.any(token_rows[usable] < 0)
        or np.any(token_rows[usable] >= ARCH["chromosome_token_offset"])
        or np.any(chromosome_codes[usable] < 0)
        or np.any(
            chromosome_codes[usable] + ARCH["chromosome_token_offset"]
            >= PE_EMBEDDING_ROWS
        )
    ):
        raise UCEActivationError("UCE token or chromosome indices exceed the embedding table")
    summary = {
        "source_genes": len(symbols),
        "genes_passing_official_min_cells": int(np.asarray(gene_passes_min_cells).sum()),
        "genes_usable_after_protein_and_chromosome_mapping": int(usable.sum()),
        "genes_dropped_after_all_filters": int(len(symbols) - usable.sum()),
        "human_species_offset": human_offset,
        "chromosome_categories_total": int(len(categories.categories)),
        "dataset_symbol_duplicates_retained_as_official_feature_rows": int(
            len(symbols) - len(set(symbols))
        ),
    }
    return token_rows, chromosome_codes, starts, usable, summary


def inspect_checkpoints(contract: Mapping[str, Any]) -> dict[str, Any]:
    import torch

    all_tokens_path = verify_record(
        contract, contract["model_files"]["all_tokens"], "UCE all-token table"
    )
    all_tokens = torch.load(
        all_tokens_path, map_location="cpu", weights_only=True, mmap=True
    )
    if not isinstance(all_tokens, torch.Tensor) or tuple(all_tokens.shape) != (
        PE_EMBEDDING_ROWS,
        ARCH["token_dim"],
    ):
        raise UCEActivationError("UCE all-token table shape differs")
    result: dict[str, Any] = {}
    for model_id in MODEL_IDS:
        checkpoint = verify_record(
            contract, contract["checkpoints"][model_id], f"{model_id} checkpoint"
        )
        state = torch.load(checkpoint, map_location="cpu", weights_only=True, mmap=True)
        inspection = inspect_state_dict(state, model_id)
        checkpoint_tokens = state["pe_embedding.weight"]
        for start in range(0, PE_EMBEDDING_ROWS, 4096):
            stop = min(start + 4096, PE_EMBEDDING_ROWS)
            if not torch.equal(all_tokens[start:stop], checkpoint_tokens[start:stop]):
                raise UCEActivationError(
                    f"{model_id} token table differs from shared all_tokens"
                )
        inspection["checkpoint_sha256"] = contract["checkpoints"][model_id]["sha256"]
        inspection["shared_token_table_equal"] = True
        result[model_id] = inspection
        del state
    return result


def build(
    count_source: Path,
    symbol_source: Path,
    split_path: Path,
    protein_embeddings_path: Path,
    species_offsets_path: Path,
    species_chrom_path: Path,
    output: Path,
    *,
    expected_rows: int = 50_000,
    expected_donors: int = 102,
    expected_studies: int = 7,
    gene_min_cells: int = 10,
    cell_min_genes: int = 25,
    run_seed: int = 20260824,
) -> dict[str, Any]:
    import anndata
    import numpy as np
    from scipy import sparse

    if output.exists() or min(expected_rows, expected_donors, expected_studies) < 1:
        raise UCEActivationError("UCE fixture output or cardinality contract differs")
    if min(gene_min_cells, cell_min_genes) < 1:
        raise UCEActivationError("UCE official filter thresholds must be positive")
    adata = anndata.read_h5ad(count_source)
    observed_obs = set(map(str, adata.obs.columns))
    if (
        not observed_obs.issubset(ALLOWED_COUNT_OBS)
        or observed_obs.intersection(FORBIDDEN_OUTCOMES)
        or not {"donor_id", "dataset_id"}.issubset(observed_obs)
    ):
        raise UCEActivationError("outcome-blind count observation firewall differs")
    matrix = adata.X
    if not sparse.issparse(matrix):
        raise UCEActivationError("UCE count input is not sparse")
    matrix = sparse.csr_matrix(matrix, copy=True)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix.sort_indices()
    if (
        matrix.shape[0] != expected_rows
        or matrix.nnz < 1
        or not np.isfinite(matrix.data).all()
        or np.any(matrix.data < 0.0)
        or not np.allclose(matrix.data, np.rint(matrix.data), rtol=0.0, atol=1e-6)
    ):
        raise UCEActivationError("UCE input is not the frozen raw nonnegative UMI matrix")
    row_ids = list(map(str, adata.obs_names))
    donors = list(map(str, adata.obs["donor_id"]))
    studies = list(map(str, adata.obs["dataset_id"]))
    count_ensembl = list(map(str, adata.var_names))
    symbol_ensembl, symbols = read_symbol_axis(symbol_source)
    if count_ensembl != symbol_ensembl:
        raise UCEActivationError("outcome-blind count and symbol feature axes differ")
    split = read_split(split_path)
    folds = validate_firewall(
        row_ids,
        donors,
        studies,
        split,
        expected_rows=expected_rows,
        expected_donors=expected_donors,
        expected_studies=expected_studies,
    )

    gene_pass = np.asarray((matrix > 0).getnnz(axis=0) >= gene_min_cells).ravel()
    official_filtered = matrix[:, gene_pass]
    raw_genes_per_cell = np.asarray((official_filtered > 0).getnnz(axis=1)).ravel()
    if np.any(raw_genes_per_cell < cell_min_genes):
        raise UCEActivationError(
            "frozen row roster contains a cell removed by official UCE min-genes filtering"
        )
    token_rows, chromosome_codes, starts, usable, universe = build_gene_universe(
        symbols,
        gene_pass,
        protein_embeddings_path,
        species_offsets_path,
        species_chrom_path,
    )
    usable_indices = np.flatnonzero(usable)
    mapped = matrix[:, usable_indices].tocsr()
    mapped_genes_per_cell = np.asarray((mapped > 0).getnnz(axis=1)).ravel()
    if np.any(mapped_genes_per_cell < 1):
        raise UCEActivationError("a frozen row has no count after UCE protein mapping")
    token_rows = token_rows[usable_indices]
    chromosome_codes = chromosome_codes[usable_indices]
    starts = starts[usable_indices]

    stage = output.with_name(f".{output.name}.{os.getpid()}.staging")
    if stage.exists():
        raise UCEActivationError("UCE fixture staging path exists")
    stage.mkdir(parents=True)
    token_path = stage / "tokens.npy"
    tokens = np.lib.format.open_memmap(
        token_path,
        mode="w+",
        dtype="int32",
        shape=(expected_rows, ARCH["pad_length"]),
    )
    content_lengths = np.empty(expected_rows, dtype=np.int16)
    sampled_unique = np.empty(expected_rows, dtype=np.int16)
    for row_index, row_id in enumerate(row_ids):
        row = mapped.getrow(row_index)
        sentence, length, unique = tokenize_cell(
            row.data,
            token_rows[row.indices],
            chromosome_codes[row.indices],
            starts[row.indices],
            row_id=row_id,
            run_seed=run_seed,
        )
        tokens[row_index] = sentence
        content_lengths[row_index] = length
        sampled_unique[row_index] = unique
    tokens.flush()
    del tokens
    np.save(stage / "content_lengths.npy", content_lengths, allow_pickle=False)
    np.save(stage / "outer_folds.npy", np.asarray(folds, dtype=np.int8), allow_pickle=False)
    (stage / "row_ids.txt").write_text("\n".join(row_ids) + "\n", encoding="utf-8")
    row_order_sha = sha256(("\n".join(row_ids) + "\n").encode("utf-8")).hexdigest()
    receipt = {
        "schema_version": "masld-bench-uce-frozen-screen-fixture-v1",
        "status": "pass_outcome_blind_uce_fixture",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "split_id": "resource_atlas_study_outer_5fold_v1",
        "rows": expected_rows,
        "donors": len(set(donors)),
        "studies": len(set(studies)),
        "outer_folds": 5,
        "source_obs_fields_read": sorted(observed_obs),
        "symbol_source_fields_read": ["var.ensembl_id", "var.source_feature_id"],
        "symbol_source_observation_fields_read": [],
        "evaluation_label_columns_read": [],
        "histology_columns_read": [],
        "sealed_outcomes_read": False,
        "official_additional_filter": True,
        "gene_min_cells": gene_min_cells,
        "cell_min_genes_before_protein_mapping": cell_min_genes,
        "cells_removed": 0,
        "row_roster_preserved": True,
        "sampling_seed": run_seed,
        "sampling_size": ARCH["sample_size"],
        "sampling_with_replacement": True,
        "sampling_weight": "log1p_raw_count",
        "sampling_rng": "per_row_numpy_default_rng_from_sha256_identity",
        "cell_order_invariant": True,
        "pad_length": ARCH["pad_length"],
        "truncation_policy": "fail_if_sentence_exceeds_1536;never_truncate",
        "cells_truncated": 0,
        "cells_with_more_than_1024_expressed_usable_genes": int(
            np.sum(mapped_genes_per_cell > ARCH["sample_size"])
        ),
        "raw_filtered_genes_per_cell": {
            "min": int(raw_genes_per_cell.min()),
            "median": float(np.median(raw_genes_per_cell)),
            "max": int(raw_genes_per_cell.max()),
        },
        "mapped_genes_per_cell": {
            "min": int(mapped_genes_per_cell.min()),
            "median": float(np.median(mapped_genes_per_cell)),
            "max": int(mapped_genes_per_cell.max()),
        },
        "sampled_unique_genes": {
            "min": int(sampled_unique.min()),
            "median": float(np.median(sampled_unique)),
            "max": int(sampled_unique.max()),
        },
        "content_lengths": {
            "min": int(content_lengths.min()),
            "median": float(np.median(content_lengths)),
            "max": int(content_lengths.max()),
            "at_pad_length": int(np.sum(content_lengths == ARCH["pad_length"])),
        },
        "gene_universe": universe,
        "row_order_sha256": row_order_sha,
        "tokens_sha256": sha256_file(token_path),
        "downstream_head_fit": False,
        "head_fitting_delegated_to_hardened_common_lane": True,
    }
    (stage / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    os.rename(stage, output)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inspect-checkpoints", action="store_true")
    arguments = parser.parse_args()
    contract = load_activation_contract(arguments.contract, arguments.project_root)
    runtime = contract["runtime"]
    runtime_root = (arguments.project_root / runtime["root"]).resolve(strict=True)
    if sha256_file(runtime_root / "ARTIFACTS.json") != runtime["artifacts_sha256"]:
        raise UCEActivationError("frozen runtime ARTIFACTS SHA-256 differs")
    inputs = contract["inputs"]
    files = {
        key: verify_record(contract, record, f"UCE {key}")
        for key, record in contract["model_files"].items()
        if key != "all_tokens"
    }
    checkpoint_inspection = inspect_checkpoints(contract) if arguments.inspect_checkpoints else {}
    receipt = build(
        verify_record(contract, inputs["outcome_blind_counts"], "outcome-blind counts"),
        verify_record(contract, inputs["symbol_axis_source"], "symbol-axis source"),
        verify_record(contract, inputs["study_split"], "study split"),
        files["human_protein_embeddings"],
        files["species_offsets"],
        files["species_chrom"],
        arguments.output,
        gene_min_cells=contract["tokenization"]["gene_min_cells"],
        cell_min_genes=contract["tokenization"]["cell_min_genes_before_protein_mapping"],
        run_seed=contract["tokenization"]["common_lane_seed"],
    )
    if checkpoint_inspection:
        receipt_path = arguments.output / "receipt.json"
        observed = json.loads(receipt_path.read_text(encoding="utf-8"))
        observed["checkpoint_inspection"] = checkpoint_inspection
        receipt_path.write_text(
            json.dumps(observed, sort_keys=True, indent=2) + "\n", encoding="utf-8"
        )
        receipt = observed
    print(json.dumps({key: receipt[key] for key in ("status", "rows", "content_lengths")}))


if __name__ == "__main__":
    main()
