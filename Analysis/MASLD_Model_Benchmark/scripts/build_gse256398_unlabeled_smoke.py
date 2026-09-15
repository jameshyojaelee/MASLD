#!/usr/bin/env python3
"""Build a donor-balanced, outcome-blind GSE256398 snRNA compatibility view."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import re

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


VIEW_ID = "gse256398_rna_unlabeled_smoke_1040_v1"
SELECTION_SEED = 20260824
ROWS_PER_DONOR = 40
EXPECTED_DONORS = 26
EXPECTED_ROWS = ROWS_PER_DONOR * EXPECTED_DONORS
EXPECTED_SOURCE_ARTIFACTS = (
    "f54ac799b22cb02912b82432929255e6400d871cc0fd99dbf9b4756f8f977a20"
)
EXPECTED_CROSSWALK_ARTIFACTS = (
    "54ca7af31f6ae223defdee4a7c1d7931f7cba4362f621fe9646fd08530c4c6da"
)
EXPECTED_QC_ARTIFACTS = (
    "3499eed58c6da15e644f5114da33b2874e3252a86fad6d6937960424e781ecbc"
)
EXPECTED_SENSITIVITY_ARTIFACTS = (
    "77b7b62b30c49fdf306039a3c31930040a5a140d23a64dd24875370cc11e91ad"
)
H5_PATTERN = re.compile(
    r"^(GSM809\d{4})_(S\d+)_CB_raw_feature_bc_matrix_filtered\.h5$"
)


class GSE256398SmokeError(RuntimeError):
    """Raised when source identities or smoke-view invariants differ."""


def decode_many(values: object) -> list[str]:
    return [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values]


def read_10x_counts(path: Path) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    import h5py
    import numpy as np
    from scipy import sparse

    with h5py.File(path, "r") as handle:
        source = handle["matrix"]
        shape = tuple(int(value) for value in source["shape"][...])
        counts = sparse.csc_matrix(
            (
                np.asarray(source["data"][...]),
                np.asarray(source["indices"][...]),
                np.asarray(source["indptr"][...]),
            ),
            shape=shape,
        ).T.tocsr()
        barcodes = decode_many(source["barcodes"][...])
        names = decode_many(source["features/name"][...])
    if counts.shape != (len(barcodes), len(names)) or counts.shape[1] != 36_601:
        raise GSE256398SmokeError("source 10x matrix axis differs")
    return counts, barcodes, names


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def priority(row_id: str, seed: int = SELECTION_SEED) -> str:
    return sha256(f"{seed}:{row_id}".encode("utf-8")).hexdigest()


def select_balanced(
    source_rows: list[dict[str, str]], sensitivity_row_ids: set[str]
) -> list[dict[str, str]]:
    by_donor: dict[str, list[dict[str, str]]] = {}
    for row in source_rows:
        if row["row_id"] in sensitivity_row_ids:
            by_donor.setdefault(row["source_sample_id"], []).append(row)
    if len(by_donor) != EXPECTED_DONORS:
        raise GSE256398SmokeError("eligible donor census differs")
    selected: list[dict[str, str]] = []
    for donor in sorted(by_donor, key=lambda value: int(value.removeprefix("S"))):
        rows = sorted(by_donor[donor], key=lambda row: (priority(row["row_id"]), row["row_id"]))
        if len(rows) < ROWS_PER_DONOR:
            raise GSE256398SmokeError(f"donor {donor} lacks 40 eligible nuclei")
        selected.extend(rows[:ROWS_PER_DONOR])
    if len(selected) != EXPECTED_ROWS or len({row["row_id"] for row in selected}) != EXPECTED_ROWS:
        raise GSE256398SmokeError("balanced selection census differs")
    return selected


def load_membership(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"row_id", "gsm", "source_sample_id", "barcode"}
    if not rows or set(rows[0]) != required or len({row["row_id"] for row in rows}) != len(rows):
        raise GSE256398SmokeError(f"membership schema or row identity differs: {path}")
    return rows


def load_allowed_genes(path: Path) -> tuple[np.ndarray, pd.DataFrame]:
    import numpy as np
    import pandas as pd

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    allowed = [row for row in rows if row["allowed_project_input"] == "true"]
    indices = np.asarray([int(row["source_feature_index"]) for row in allowed], dtype=np.int64)
    stable = [row["source_stable_gene_id"] for row in allowed]
    if (
        len(rows) != 36_601
        or len(allowed) != 35_455
        or len(set(stable)) != len(stable)
        or np.any(indices[1:] <= indices[:-1])
    ):
        raise GSE256398SmokeError("allowed GENCODE v49 gene axis differs")
    var = pd.DataFrame(
        {
            "ensembl_id": stable,
            "gencode49_gene_id": [row["gencode49_gene_id"] for row in allowed],
            "gencode49_gene_name": [row["gencode49_gene_name"] for row in allowed],
            "source_feature_id": [row["source_feature_id"] for row in allowed],
            "source_feature_name": [row["source_feature_name"] for row in allowed],
            "mapping_state": [row["mapping_state"] for row in allowed],
        },
        index=pd.Index(stable, name="stable_ensembl_gene_id"),
    )
    return indices, var


def materialize(
    *, source: Path, crosswalk: Path, selected: list[dict[str, str]], output: Path
) -> dict[str, object]:
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    allowed_indices, var = load_allowed_genes(crosswalk / "gene_crosswalk.tsv")
    by_donor: dict[str, list[dict[str, str]]] = {}
    for row in selected:
        by_donor.setdefault(row["source_sample_id"], []).append(row)
    h5_by_donor: dict[str, Path] = {}
    for path in sorted((source / "h5").glob("*.h5")):
        match = H5_PATTERN.fullmatch(path.name)
        if match is None:
            raise GSE256398SmokeError("source H5 filename differs")
        h5_by_donor[match.group(2)] = path
    if set(h5_by_donor) != set(by_donor):
        raise GSE256398SmokeError("selected donors do not match source H5 roster")

    chunks: list[sparse.csr_matrix] = []
    ordered: list[dict[str, str]] = []
    for donor in sorted(by_donor, key=lambda value: int(value.removeprefix("S"))):
        donor_rows = by_donor[donor]
        counts, barcodes, _ = read_10x_counts(h5_by_donor[donor])
        barcode_to_index = {barcode: index for index, barcode in enumerate(barcodes)}
        if len(barcode_to_index) != len(barcodes):
            raise GSE256398SmokeError("source barcode axis is duplicated")
        try:
            row_indices = [barcode_to_index[row["barcode"]] for row in donor_rows]
        except KeyError as error:
            raise GSE256398SmokeError(f"selected barcode is absent for {donor}") from error
        chunks.append(counts[row_indices, :][:, allowed_indices].tocsr())
        ordered.extend(donor_rows)
    matrix = sparse.vstack(chunks, format="csr")
    if (
        matrix.shape != (EXPECTED_ROWS, 35_455)
        or matrix.nnz < EXPECTED_ROWS
        or np.any(matrix.data <= 0)
        or np.any(matrix.data != np.floor(matrix.data))
    ):
        raise GSE256398SmokeError("materialized raw-count matrix differs")
    n_counts = np.asarray(matrix.sum(axis=1)).ravel().astype(np.int64)
    detected = np.diff(matrix.indptr).astype(np.int64)
    obs = pd.DataFrame(
        {
            "row_id": [row["row_id"] for row in ordered],
            "gsm": [row["gsm"] for row in ordered],
            "source_sample_id": [row["source_sample_id"] for row in ordered],
            "barcode": [row["barcode"] for row in ordered],
            "n_counts": n_counts,
            "detected_allowed_genes": detected,
            "cell_state_label_state": "structurally_missing",
        },
        index=pd.Index([row["row_id"] for row in ordered], name="nucleus_row_id"),
    )
    value = ad.AnnData(X=matrix, obs=obs, var=var)
    value.uns["masld_bench_dataset_view"] = {
        "view_id": VIEW_ID,
        "selection_seed": SELECTION_SEED,
        "rows_per_donor": ROWS_PER_DONOR,
        "selection_outcomes_used": False,
        "sealed_outcomes_used": False,
        "normalization_fitted": False,
        "variable_feature_selection_fitted": False,
    }
    value.write_h5ad(output, compression="gzip")
    reloaded = ad.read_h5ad(output)
    if reloaded.shape != value.shape or list(reloaded.obs_names) != list(value.obs_names):
        raise GSE256398SmokeError("written H5AD identity differs")
    if np.any(np.asarray(reloaded.X.sum(axis=1)).ravel() != reloaded.obs["n_counts"]):
        raise GSE256398SmokeError("written H5AD counts do not rederive")
    return {
        "shape": list(matrix.shape),
        "nnz": int(matrix.nnz),
        "matrix_data_sha256": sha256(matrix.data.tobytes()).hexdigest(),
        "matrix_indices_sha256": sha256(matrix.indices.tobytes()).hexdigest(),
        "matrix_indptr_sha256": sha256(matrix.indptr.tobytes()).hexdigest(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--qc", type=Path, required=True)
    parser.add_argument("--sensitivity", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise GSE256398SmokeError("smoke output exists")
    expected = (
        (args.source, EXPECTED_SOURCE_ARTIFACTS),
        (args.crosswalk, EXPECTED_CROSSWALK_ARTIFACTS),
        (args.qc, EXPECTED_QC_ARTIFACTS),
        (args.sensitivity, EXPECTED_SENSITIVITY_ARTIFACTS),
    )
    for root, digest in expected:
        verify_frozen_tree(root)
        if sha256_file(root / "ARTIFACTS.json") != digest:
            raise GSE256398SmokeError(f"frozen input differs: {root}")
    source_rows = load_membership(args.qc / "retained_membership.tsv")
    sensitivity_rows = load_membership(
        args.sensitivity / "expected_rate_capped_membership.tsv"
    )
    selected = select_balanced(source_rows, {row["row_id"] for row in sensitivity_rows})
    args.output.mkdir(parents=True, exist_ok=False)
    selection_path = args.output / "selection.tsv"
    with selection_path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "cell_id",
                "donor_id",
                "broad_label",
                "row_id",
                "gsm",
                "source_sample_id",
                "barcode",
                "selection_priority",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in selected:
            writer.writerow(
                {
                    "cell_id": row["row_id"],
                    "donor_id": f"gse256398:{row['source_sample_id']}",
                    "broad_label": "unlabeled",
                    **row,
                    "selection_priority": priority(row["row_id"]),
                }
            )
    ontology = {
        "schema_version": "masld-bench-cell-ontology-v1",
        "ontology_id": "unlabeled_snRNA_v1",
        "source_field": "barcode_level_cell_label_state",
        "classes": {"unlabeled": ["structurally_missing"]},
        "barcode_level_cell_labels": "structurally_missing",
        "supervised_cell_state_scoring_allowed": False,
        "selection_outcomes_used": False,
        "sealed_outcomes_used": False,
    }
    ontology_path = args.output / "unlabeled_snRNA_v1.json"
    write_json_exclusive(ontology_path, ontology)
    data_path = args.output / "gse256398_rna_unlabeled_smoke_1040.h5ad"
    matrix = materialize(
        source=args.source,
        crosswalk=args.crosswalk,
        selected=selected,
        output=data_path,
    )
    derivation = {
        "schema_version": "masld-bench-unlabeled-snrna-smoke-subset-v1",
        "status": "pass_compatibility_smoke_only",
        "dataset_id": "gse256398",
        "subset_id": VIEW_ID,
        "view_id": VIEW_ID,
        "row_count": EXPECTED_ROWS,
        "biological_unit_count": EXPECTED_DONORS,
        "rows_per_donor": ROWS_PER_DONOR,
        "modalities": ["single_nucleus_rna"],
        "pairing_levels": ["same_study_unpaired"],
        "selection": {
            "cell_budget": EXPECTED_ROWS,
            "counts_by_class": {"unlabeled": EXPECTED_ROWS},
            "donors": EXPECTED_DONORS,
            "rows_per_donor": ROWS_PER_DONOR,
            "selection_policy": "intersection_of_source_exact_and_expected_rate_capped_QC_then_40_per_donor_SHA256_row_priority",
            "selection_seed": SELECTION_SEED,
            "selection_outcomes_used": False,
            "sealed_outcomes_used": False,
        },
        "donor_disease_age_sex_histology_in_model_input": False,
        "barcode_level_cell_labels": "structurally_missing",
        "features_native": 36_601,
        "features_allowed_stable_id_exact_gencode49": 35_455,
        "features_masked": 1_146,
        "normalization_or_hvg_fitted": False,
        "biological_evaluation_or_model_selection_allowed": False,
        "source_artifacts_sha256": EXPECTED_SOURCE_ARTIFACTS,
        "crosswalk_artifacts_sha256": EXPECTED_CROSSWALK_ARTIFACTS,
        "qc_artifacts_sha256": EXPECTED_QC_ARTIFACTS,
        "qc_sensitivity_artifacts_sha256": EXPECTED_SENSITIVITY_ARTIFACTS,
        "matrix": matrix,
        "artifacts": {
            "h5ad": {
                "path": data_path.resolve().as_posix(),
                "sha256": sha256_file(data_path),
                "size_bytes": data_path.stat().st_size,
            },
            "selection": {
                "path": selection_path.resolve().as_posix(),
                "sha256": sha256_file(selection_path),
                "size_bytes": selection_path.stat().st_size,
            },
            "ontology": {
                "path": ontology_path.resolve().as_posix(),
                "sha256": sha256_file(ontology_path),
                "size_bytes": ontology_path.stat().st_size,
            },
        },
        "unlabeled_input_contract": {
            "matrix": "raw_integer_CellBender_counts_on_stable_ID_exact_GENCODE_v49_genes",
            "barcode_cell_labels": "structurally_missing",
            "donor_covariates_in_model_input": False,
            "ready": True,
        },
    }
    write_json_exclusive(args.output / "smoke_subset_manifest.json", derivation)
    freeze_tree(
        args.output,
        {
            "artifact_class": "unlabeled_snrna_smoke_subset",
            "dataset_id": "gse256398",
            "subset_id": VIEW_ID,
            "rows": EXPECTED_ROWS,
            "donors": EXPECTED_DONORS,
            "outcome_independent": True,
            "compatibility_smoke_only": True,
        },
    )
    verify_frozen_tree(args.output)
    print(json.dumps(derivation, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
