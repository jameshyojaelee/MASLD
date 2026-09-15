#!/usr/bin/env python3
"""Build a label-free, SCimilarity-aligned fixture for the frozen 50k screen."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
from typing import Any


ROWS = 50_000
DONORS = 102
STUDIES = 7
INPUT_DIMENSION = 28_231
MINIMUM_GENE_OVERLAP = 5_000
DATASET_VIEW_ID = "resource_atlas_frozen_screen_50000_v1"
SPLIT_ID = "resource_atlas_study_outer_5fold_v1"
ALLOWED_OBS = {"row_id", "donor_id", "dataset_id", "outer_fold", "assay"}
FORBIDDEN_NAMES = {
    "age",
    "broad_label",
    "cell_type",
    "fibrosis",
    "label",
    "mash",
    "masld",
    "nas",
    "sex",
    "source_cell_type",
    "stage",
}
EXPOSURE_BY_STUDY = {
    "GSE136103": "reference_only",
    "GSE174748": "clean_declared",
    "GSE185477": "encoder_seen",
    "GSE189600": "clean_declared",
    "GSE202379": "clean_declared",
    "GSE244832": "clean_declared",
    "Liver_Atlas": "clean_declared",
}
AGGREGATE_EXPOSURE_STATUS = "encoder_seen"


class SCimilarityFixtureError(ValueError):
    """Raised when a label-free SCimilarity fixture requirement differs."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _strings(values: Any) -> list[str]:
    return [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values]


def _read_h5_column(group: Any, key: str) -> list[str]:
    import h5py
    import numpy as np

    if key not in group:
        raise SCimilarityFixtureError(f"required H5AD column {key!r} is absent")
    node = group[key]
    if isinstance(node, h5py.Dataset):
        return _strings(node[:])
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        categories = _strings(node["categories"][:])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        if np.any(codes < 0) or np.any(codes >= len(categories)):
            raise SCimilarityFixtureError(f"categorical H5AD column {key!r} is invalid")
        return [categories[int(code)] for code in codes]
    raise SCimilarityFixtureError(f"unsupported H5AD encoding for {key!r}")


def _read_symbol_crosswalk(path: Path) -> tuple[list[str], list[str]]:
    """Read feature metadata only; observation/outcome groups are never opened."""

    import h5py

    with h5py.File(path, "r") as handle:
        if "var" not in handle:
            raise SCimilarityFixtureError("symbol source lacks H5AD var metadata")
        var = handle["var"]
        ensembl_ids = _read_h5_column(var, "ensembl_id")
        symbols = _read_h5_column(var, "source_feature_id")
    if len(ensembl_ids) != len(symbols) or len(set(ensembl_ids)) != len(ensembl_ids):
        raise SCimilarityFixtureError("symbol crosswalk identity differs")
    return ensembl_ids, symbols


def _read_split(path: Path) -> dict[str, tuple[str, str, int]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["row_id", "donor_id", "dataset", "outer_fold"]:
            raise SCimilarityFixtureError("study-held-out split schema differs")
        result: dict[str, tuple[str, str, int]] = {}
        for record in reader:
            row_id = record["row_id"]
            if row_id in result:
                raise SCimilarityFixtureError("study-held-out row identifiers collide")
            fold = int(record["outer_fold"])
            if fold not in range(5):
                raise SCimilarityFixtureError("study-held-out outer fold differs")
            result[row_id] = (record["donor_id"], record["dataset"], fold)
    return result


def _load_gene_order(path: Path, expected_dimension: int) -> list[str]:
    genes = path.read_text(encoding="utf-8").splitlines()
    if len(genes) != expected_dimension or len(set(genes)) != len(genes) or any(not gene for gene in genes):
        raise SCimilarityFixtureError("SCimilarity gene order differs")
    return genes


def _study_observability(
    library_manifest: Path,
    library_root: Path,
    gene_order: list[str],
    datasets: set[str],
) -> tuple[list[str], Any, list[dict[str, Any]]]:
    """Derive conservative study panels from every analysis-eligible library."""

    import h5py
    import numpy as np

    panels: dict[str, list[set[str]]] = {dataset: [] for dataset in sorted(datasets)}
    panel_hashes: dict[str, set[str]] = {dataset: set() for dataset in sorted(datasets)}
    panel_sizes: dict[str, set[int]] = {dataset: set() for dataset in sorted(datasets)}
    with library_manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"dataset", "library_id", "analysis_eligible"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise SCimilarityFixtureError("library manifest schema differs")
        for record in reader:
            dataset = record["dataset"]
            if dataset not in datasets or record["analysis_eligible"] != "True":
                continue
            path = library_root / f"{record['library_id']}.h5ad"
            if not path.is_file():
                raise SCimilarityFixtureError(f"analysis-eligible library is absent: {path.name}")
            with h5py.File(path, "r") as source:
                if "var" not in source:
                    raise SCimilarityFixtureError(f"library lacks var metadata: {path.name}")
                var = source["var"]
                index_key = var.attrs.get("_index", "_index")
                if isinstance(index_key, bytes):
                    index_key = index_key.decode("utf-8")
                names = _read_h5_column(var, str(index_key))
            if len(names) != len(set(names)):
                raise SCimilarityFixtureError(f"library feature symbols are not unique: {path.name}")
            panels[dataset].append(set(names))
            panel_hashes[dataset].add(sha256(("\n".join(names) + "\n").encode("utf-8")).hexdigest())
            panel_sizes[dataset].add(len(names))
    studies = sorted(datasets)
    masks = np.zeros((len(studies), len(gene_order)), dtype=np.bool_)
    records: list[dict[str, Any]] = []
    for index, dataset in enumerate(studies):
        values = panels[dataset]
        if not values:
            raise SCimilarityFixtureError(f"study has no analysis-eligible library panels: {dataset}")
        intersection = set.intersection(*values)
        masks[index] = np.asarray([gene in intersection for gene in gene_order], dtype=np.bool_)
        records.append(
            {
                "study_index": index,
                "dataset_id": dataset,
                "n_libraries": len(values),
                "distinct_library_panel_hashes": len(panel_hashes[dataset]),
                "library_panel_gene_counts": ",".join(map(str, sorted(panel_sizes[dataset]))),
                "model_genes_in_every_library_panel": int(masks[index].sum()),
            }
        )
    return studies, masks, records


def _align_and_lognormalize(matrix: Any, symbols: list[str], gene_order: list[str], minimum_overlap: int):
    import numpy as np
    from scipy import sparse

    if not sparse.issparse(matrix):
        raise SCimilarityFixtureError("outcome-blind counts are not sparse")
    matrix = sparse.csr_matrix(matrix, dtype=np.float32, copy=True)
    matrix.sort_indices()
    if (
        matrix.nnz < 1
        or float(matrix.data.min()) < 0.0
        or not np.allclose(matrix.data, np.rint(matrix.data), rtol=0.0, atol=1e-6)
    ):
        raise SCimilarityFixtureError("outcome-blind matrix is not raw nonnegative UMI counts")
    target = {gene: index for index, gene in enumerate(gene_order)}
    source_indices: list[int] = []
    target_indices: list[int] = []
    for index, symbol in enumerate(symbols):
        target_index = target.get(symbol)
        if target_index is not None:
            source_indices.append(index)
            target_indices.append(target_index)
    observed_targets = sorted(set(target_indices))
    if len(observed_targets) < minimum_overlap:
        raise SCimilarityFixtureError(
            f"SCimilarity gene overlap {len(observed_targets)} is below {minimum_overlap}"
        )
    mapping = sparse.csr_matrix(
        (
            np.ones(len(source_indices), dtype=np.float32),
            (np.asarray(source_indices), np.asarray(target_indices)),
        ),
        shape=(matrix.shape[1], len(gene_order)),
    )
    aligned = sparse.csr_matrix(matrix @ mapping, dtype=np.float32)
    aligned.sum_duplicates()
    aligned.sort_indices()
    row_totals = np.asarray(aligned.sum(axis=1)).ravel()
    if np.any(~np.isfinite(row_totals)) or np.any(row_totals <= 0.0):
        raise SCimilarityFixtureError("a row has no counts in the SCimilarity gene order")
    scalers = (10_000.0 / row_totals).astype(np.float32, copy=False)
    aligned.data *= np.repeat(scalers, np.diff(aligned.indptr))
    np.log1p(aligned.data, out=aligned.data)
    observed_mask = np.zeros(len(gene_order), dtype=np.bool_)
    observed_mask[observed_targets] = True
    return aligned, observed_mask, len(source_indices)


def build(
    source: Path,
    symbol_source: Path,
    split_path: Path,
    gene_order_path: Path,
    library_manifest: Path,
    library_root: Path,
    output: Path,
    *,
    expected_source_sha256: str,
    expected_symbol_source_sha256: str,
    expected_split_sha256: str,
    expected_gene_order_sha256: str,
    expected_library_manifest_sha256: str,
    expected_rows: int = ROWS,
    expected_donors: int = DONORS,
    expected_studies: int = STUDIES,
    expected_dimension: int = INPUT_DIMENSION,
    minimum_overlap: int = MINIMUM_GENE_OVERLAP,
) -> dict[str, Any]:
    import anndata
    import numpy as np
    from scipy import sparse

    if output.exists() or any(value < 1 for value in (expected_rows, expected_donors, expected_studies)):
        raise SCimilarityFixtureError("fixture output/cardinality contract differs")
    expected_files = (
        (source, expected_source_sha256, "outcome-blind count source"),
        (symbol_source, expected_symbol_source_sha256, "feature-symbol source"),
        (split_path, expected_split_sha256, "study-held-out split"),
        (gene_order_path, expected_gene_order_sha256, "SCimilarity gene order"),
        (library_manifest, expected_library_manifest_sha256, "library manifest"),
    )
    for path, expected, name in expected_files:
        if sha256_file(path) != expected:
            raise SCimilarityFixtureError(f"{name} hash differs")

    adata = anndata.read_h5ad(source)
    observed_obs = set(map(str, adata.obs.columns))
    if not observed_obs.issubset(ALLOWED_OBS) or observed_obs.intersection(FORBIDDEN_NAMES):
        raise SCimilarityFixtureError("outcome-blind observation firewall differs")
    if not {"donor_id", "dataset_id"}.issubset(observed_obs):
        raise SCimilarityFixtureError("outcome-blind donor/study join is absent")
    if adata.n_obs != expected_rows or len(set(map(str, adata.obs_names))) != expected_rows:
        raise SCimilarityFixtureError("outcome-blind row identity differs")
    row_ids = list(map(str, adata.obs_names))
    donors = list(map(str, adata.obs["donor_id"]))
    datasets = list(map(str, adata.obs["dataset_id"]))
    if len(set(donors)) != expected_donors or len(set(datasets)) != expected_studies:
        raise SCimilarityFixtureError("outcome-blind donor/study cardinality differs")
    if not set(datasets).issubset(EXPOSURE_BY_STUDY) or (
        expected_studies == STUDIES and set(datasets) != set(EXPOSURE_BY_STUDY)
    ):
        raise SCimilarityFixtureError("SCimilarity development exposure roster differs")
    ensembl_ids = list(map(str, adata.var_names))
    mapping_ids, symbols = _read_symbol_crosswalk(symbol_source)
    if ensembl_ids != mapping_ids:
        raise SCimilarityFixtureError("label-free counts and feature-symbol order differ")
    gene_order = _load_gene_order(gene_order_path, expected_dimension)

    split = _read_split(split_path)
    if set(split) != set(row_ids):
        raise SCimilarityFixtureError("count rows and study split rows differ")
    folds: list[int] = []
    donor_folds: dict[str, int] = {}
    study_folds: dict[str, int] = {}
    for row_id, donor, dataset in zip(row_ids, donors, datasets, strict=True):
        split_donor, split_dataset, fold = split[row_id]
        if (donor, dataset) != (split_donor, split_dataset):
            raise SCimilarityFixtureError("count and study-split metadata join differs")
        if donor_folds.setdefault(donor, fold) != fold:
            raise SCimilarityFixtureError("one donor crosses outer folds")
        if study_folds.setdefault(dataset, fold) != fold:
            raise SCimilarityFixtureError("one study crosses outer folds")
        folds.append(fold)
    if set(folds) != set(range(5)):
        raise SCimilarityFixtureError("study-held-out split lacks an outer fold")

    normalized, observed_mask, matched_source_features = _align_and_lognormalize(
        adata.X, symbols, gene_order, minimum_overlap
    )
    study_order, library_panel_masks, observability_records = _study_observability(
        library_manifest, library_root, gene_order, set(datasets)
    )
    study_observability = library_panel_masks & observed_mask[None, :]
    for record, mask in zip(observability_records, study_observability, strict=True):
        record["model_genes_available_in_integrated_input"] = int(mask.sum())
    global_panel = bool(np.all(study_observability == study_observability[0]))
    if normalized.shape != (expected_rows, expected_dimension) or not sparse.isspmatrix_csr(normalized):
        raise SCimilarityFixtureError("aligned SCimilarity matrix shape differs")
    stage = output.with_name(f".{output.name}.{os.getpid()}.staging")
    if stage.exists():
        raise SCimilarityFixtureError("fixture staging path already exists")
    stage.mkdir(parents=True)
    sparse.save_npz(stage / "normalized_counts.npz", normalized, compressed=True)
    np.save(stage / "observed_gene_mask.npy", observed_mask, allow_pickle=False)
    np.save(stage / "study_gene_observability.npy", study_observability, allow_pickle=False)
    np.save(stage / "outer_folds.npy", np.asarray(folds, dtype=np.int8), allow_pickle=False)
    (stage / "embedding_row_order.txt").write_text("\n".join(row_ids) + "\n", encoding="utf-8")
    with (stage / "row_contract.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("row_position", "row_id", "donor_id", "dataset_id", "outer_fold", "exposure_state"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for position, (row_id, donor, dataset, fold) in enumerate(
            zip(row_ids, donors, datasets, folds, strict=True)
        ):
            writer.writerow(
                {
                    "row_position": position,
                    "row_id": row_id,
                    "donor_id": donor,
                    "dataset_id": dataset,
                    "outer_fold": fold,
                    "exposure_state": EXPOSURE_BY_STUDY[dataset],
                }
            )
    with (stage / "study_gene_observability.tsv").open("x", encoding="utf-8", newline="") as handle:
        fieldnames = list(observability_records[0])
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(observability_records)
    receipt = {
        "schema_version": "masld-bench-scimilarity-v1-1-frozen-screen-fixture-v1",
        "status": "pass_outcome_blind_fixture",
        "model_id": "scimilarity_v1_1",
        "dataset_view_id": DATASET_VIEW_ID,
        "split_id": SPLIT_ID,
        "rows": expected_rows,
        "donors": len(set(donors)),
        "studies": len(set(datasets)),
        "outer_folds": 5,
        "input_dimension": expected_dimension,
        "gene_overlap": int(observed_mask.sum()),
        "gene_overlap_fraction": float(observed_mask.mean()),
        "matched_source_features_including_duplicates": matched_source_features,
        "structurally_missing_genes": int((~observed_mask).sum()),
        "structural_missingness_encoded_as_zero_fill": True,
        "normalization": "per_cell_total_10000_over_retained_model_genes_then_natural_log1p",
        "duplicate_gene_symbols": "summed_before_alignment",
        "minimum_gene_overlap_enforced": minimum_overlap,
        "source_obs_fields_read": sorted(observed_obs),
        "feature_symbol_source_groups_read": ["var/ensembl_id", "var/source_feature_id"],
        "feature_symbol_source_observations_read": False,
        "evaluation_label_columns_read": [],
        "histology_columns_read": [],
        "sealed_outcomes_read": False,
        "released_reference_index_used": False,
        "aggregate_exposure_status": AGGREGATE_EXPOSURE_STATUS,
        "exposure_by_study": EXPOSURE_BY_STUDY,
        "study_gene_observability_order": study_order,
        "study_gene_observability_shape": list(study_observability.shape),
        "study_gene_observability_is_global": global_panel,
        "study_gene_observability_mask_semantics": (
            "true iff the exact SCimilarity gene symbol is present in every analysis-eligible "
            "input library panel for that study and retained in the integrated raw-count axis"
        ),
        "native_model_accepts_observability_mask": False,
        "native_fixed_vocabulary_exception": (
            "SCimilarity requires a dense fixed 28,231-gene vector and zero-fills the 7,622 "
            "globally unavailable genes. The preserved mask records input compatibility; "
            "those zeros are not interpreted as observed biological zeros."
        ),
        "source_sha256": expected_source_sha256,
        "symbol_source_sha256": expected_symbol_source_sha256,
        "split_sha256": expected_split_sha256,
        "gene_order_sha256": expected_gene_order_sha256,
        "library_manifest_sha256": expected_library_manifest_sha256,
        "row_order_sha256": sha256(("\n".join(row_ids) + "\n").encode("utf-8")).hexdigest(),
    }
    (stage / "fixture_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    os.rename(stage, output)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--symbol-source", type=Path, required=True)
    parser.add_argument("--split", type=Path, required=True)
    parser.add_argument("--gene-order", type=Path, required=True)
    parser.add_argument("--library-manifest", type=Path, required=True)
    parser.add_argument("--library-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-symbol-source-sha256", required=True)
    parser.add_argument("--expected-split-sha256", required=True)
    parser.add_argument("--expected-gene-order-sha256", required=True)
    parser.add_argument("--expected-library-manifest-sha256", required=True)
    arguments = parser.parse_args()
    result = build(
        arguments.source,
        arguments.symbol_source,
        arguments.split,
        arguments.gene_order,
        arguments.library_manifest,
        arguments.library_root,
        arguments.output,
        expected_source_sha256=arguments.expected_source_sha256,
        expected_symbol_source_sha256=arguments.expected_symbol_source_sha256,
        expected_split_sha256=arguments.expected_split_sha256,
        expected_gene_order_sha256=arguments.expected_gene_order_sha256,
        expected_library_manifest_sha256=arguments.expected_library_manifest_sha256,
    )
    print(json.dumps({key: result[key] for key in ("status", "rows", "gene_overlap")}))


if __name__ == "__main__":
    main()
