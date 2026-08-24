#!/usr/bin/env python3
"""Build a deterministic donor-balanced Geneformer smoke subset."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import heapq
import json
from pathlib import Path
from typing import Any, Iterable


ONTOLOGY: dict[str, tuple[str, ...]] = {
    "cholangiocyte": ("Cholangiocytes",),
    "endothelial": ("Endothelial cells",),
    "hepatocyte": ("Hepatocytes",),
    "immune": (
        "B cells",
        "Basophils",
        "Circulating NK/NKT",
        "Macrophages",
        "Mono+mono derived cells",
        "Neutrophils",
        "Plasma cells",
        "Resident NK",
        "T cells",
        "cDC1s",
        "cDC2s",
        "pDCs",
    ),
    "mesenchymal_stromal": ("Fibroblasts",),
}
CELL_BUDGET = 1_000
CELLS_PER_CLASS = CELL_BUDGET // len(ONTOLOGY)
SELECTION_SEED = 20260821
PER_DONOR_CLASS_BUFFER = 16


class SmokeSubsetError(ValueError):
    """Raised when the source atlas cannot satisfy the frozen smoke contract."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def write_json_exclusive(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(
            value,
            handle,
            sort_keys=True,
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        handle.write("\n")


def freeze_output_tree(output: Path, *, subset_id: str, manifest_id: str) -> None:
    records = []
    for path in sorted(output.rglob("*")):
        if path.is_symlink():
            raise SmokeSubsetError(f"derived smoke tree contains a symlink: {path}")
        if path.is_file() and path.name not in {"ARTIFACTS.json", "COMPLETE"}:
            records.append(
                {
                    "path": path.relative_to(output).as_posix(),
                    "sha256": sha256_file(path),
                    "size_bytes": path.stat().st_size,
                }
            )
    artifacts = {
        "schema_version": "masld-bench-artifacts-v1",
        "metadata": {
            "artifact_class": "geneformer_smoke_subset",
            "artifact_id": manifest_id,
            "subset_id": subset_id,
        },
        "artifacts": records,
    }
    artifacts_path = output / "ARTIFACTS.json"
    write_json_exclusive(artifacts_path, artifacts)
    complete = {
        "schema_version": "masld-bench-complete-v1",
        "manifest_sha256": sha256_file(artifacts_path),
        "artifact_count": len(records),
    }
    complete_path = output / "COMPLETE"
    write_json_exclusive(complete_path, complete)
    for record in records:
        path = output / record["path"]
        if sha256_file(path) != record["sha256"] or path.stat().st_size != record["size_bytes"]:
            raise SmokeSubsetError(f"derived artifact changed while freezing: {path}")
        path.chmod(0o440)
    artifacts_path.chmod(0o440)
    complete_path.chmod(0o440)
    output.chmod(0o550)


def priority(*parts: object) -> int:
    text = "\x1f".join(map(str, (SELECTION_SEED, *parts)))
    return int.from_bytes(hashlib.sha256(text.encode("utf-8")).digest(), "big")


def _raw_to_broad() -> dict[str, str]:
    result: dict[str, str] = {}
    for broad, raw_labels in ONTOLOGY.items():
        for raw in raw_labels:
            if raw in result:
                raise SmokeSubsetError(f"raw label occurs twice in ontology: {raw}")
            result[raw] = broad
    return result


def select_cells(cell_manifest: Path) -> list[dict[str, str | int]]:
    raw_to_broad = _raw_to_broad()
    pools: dict[str, dict[str, list[tuple[int, int, dict[str, str | int]]]]] = {
        broad: {} for broad in ONTOLOGY
    }
    opener = gzip.open if cell_manifest.suffix == ".gz" else open
    with opener(cell_manifest, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "cell_id",
            "source_cell_id",
            "library_id",
            "donor_id",
            "dataset",
            "cell_type",
            "analysis_eligible",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise SmokeSubsetError("cell manifest lacks required source fields")
        for row_index, row in enumerate(reader):
            if row["analysis_eligible"] != "True":
                continue
            broad = raw_to_broad.get(row["cell_type"])
            if broad is None:
                raise SmokeSubsetError(
                    f"analysis-eligible cell has unmapped type {row['cell_type']!r}"
                )
            record: dict[str, str | int] = {
                "source_row_index": row_index,
                "cell_id": row["cell_id"],
                "source_cell_id": row["source_cell_id"],
                "library_id": row["library_id"],
                "donor_id": row["donor_id"],
                "dataset": row["dataset"],
                "source_cell_type": row["cell_type"],
                "broad_label": broad,
            }
            score = priority("cell", broad, row["donor_id"], row["cell_id"])
            donor_pool = pools[broad].setdefault(row["donor_id"], [])
            item = (-score, -row_index, record)
            if len(donor_pool) < PER_DONOR_CLASS_BUFFER:
                heapq.heappush(donor_pool, item)
            elif item > donor_pool[0]:
                heapq.heapreplace(donor_pool, item)

    selected: list[dict[str, str | int]] = []
    for broad in sorted(ONTOLOGY):
        donor_rows = {
            donor: [
                item[2]
                for item in sorted(
                    heap,
                    key=lambda item: (-item[0], -item[1]),
                )
            ]
            for donor, heap in pools[broad].items()
        }
        donor_order = sorted(
            donor_rows,
            key=lambda donor: (priority("donor", broad, donor), donor),
        )
        if not donor_order:
            raise SmokeSubsetError(f"no eligible donors for {broad}")
        offset = 0
        class_rows: list[dict[str, str | int]] = []
        while len(class_rows) < CELLS_PER_CLASS:
            progressed = False
            for donor in donor_order:
                values = donor_rows[donor]
                if offset < len(values):
                    class_rows.append(values[offset])
                    progressed = True
                    if len(class_rows) == CELLS_PER_CLASS:
                        break
            if not progressed:
                raise SmokeSubsetError(
                    f"donor-balanced buffer cannot supply {CELLS_PER_CLASS} {broad} cells"
                )
            offset += 1
        selected.extend(class_rows)
    if len(selected) != CELL_BUDGET:
        raise SmokeSubsetError("smoke selection did not produce exactly 1,000 cells")
    if len({str(row["cell_id"]) for row in selected}) != len(selected):
        raise SmokeSubsetError("smoke selection repeats a cell")
    selected.sort(key=lambda row: int(row["source_row_index"]))
    return selected


def _strings(dataset: Any) -> list[str]:
    try:
        values: Iterable[Any] = dataset.asstr()[:]
    except (AttributeError, TypeError, ValueError):
        values = dataset[:]
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def materialize_subset(
    *, atlas: Path, selected: list[dict[str, str | int]], output_h5ad: Path
) -> dict[str, Any]:
    import anndata as ad
    import h5py
    import numpy as np
    import pandas as pd
    from scipy import sparse

    row_indices = np.asarray(
        [int(row["source_row_index"]) for row in selected], dtype=np.int64
    )
    with h5py.File(atlas, "r") as handle:
        raw_x = handle["raw/X"]
        indptr = np.asarray(raw_x["indptr"][:], dtype=np.int64)
        source_cell_ids = _strings(handle["obs/_index"])
        raw_var = handle["raw/var"]
        source_features = _strings(raw_var["_index"])
        source_gene_ids = _strings(raw_var["gene_ids"])
        if len(source_features) != len(source_gene_ids):
            raise SmokeSubsetError("raw variable axes differ")
        ensembl_ids = [value.split(".", 1)[0] for value in source_gene_ids]
        if (
            any(not value.startswith("ENSG") for value in ensembl_ids)
            or len(set(ensembl_ids)) != len(ensembl_ids)
        ):
            raise SmokeSubsetError("source gene_ids cannot form a unique Ensembl axis")

        data_parts: list[Any] = []
        index_parts: list[Any] = []
        subset_indptr = [0]
        for row_index, row in zip(row_indices, selected, strict=True):
            if source_cell_ids[int(row_index)] != str(row["source_cell_id"]):
                raise SmokeSubsetError(
                    f"cell-manifest order differs at source row {row_index}"
                )
            start = int(indptr[int(row_index)])
            end = int(indptr[int(row_index) + 1])
            data = raw_x["data"][start:end]
            indices = raw_x["indices"][start:end]
            data_parts.append(data)
            index_parts.append(indices)
            subset_indptr.append(subset_indptr[-1] + len(data))

    data = np.concatenate(data_parts)
    indices = np.concatenate(index_parts)
    if np.any(data < 0) or np.any(data != np.floor(data)):
        raise SmokeSubsetError("selected raw counts are not nonnegative integers")
    matrix = sparse.csr_matrix(
        (data, indices, np.asarray(subset_indptr, dtype=np.int64)),
        shape=(len(selected), len(ensembl_ids)),
    )
    n_counts = np.asarray(matrix.sum(axis=1)).ravel().astype(np.int64)
    if np.any(n_counts <= 0):
        raise SmokeSubsetError("selected cell has no raw counts")

    obs = pd.DataFrame(selected).set_index("cell_id", drop=True)
    obs["n_counts"] = n_counts
    var = pd.DataFrame(
        {
            "ensembl_id": ensembl_ids,
            "source_feature_id": source_features,
            "source_gene_id": source_gene_ids,
        },
        index=pd.Index(ensembl_ids, name="ensembl_id_index"),
    )
    value = ad.AnnData(X=matrix, obs=obs, var=var)
    value.uns["masld_bench_smoke_contract"] = {
        "cell_budget": CELL_BUDGET,
        "cells_per_class": CELLS_PER_CLASS,
        "selection_seed": SELECTION_SEED,
        "selection_policy": "donor_round_robin_then_sha256_cell_priority",
    }
    value.write_h5ad(output_h5ad, compression="gzip")

    reloaded = ad.read_h5ad(output_h5ad)
    if reloaded.shape != value.shape:
        raise SmokeSubsetError("written smoke H5AD shape changed")
    if "n_counts" not in reloaded.obs or "ensembl_id" not in reloaded.var:
        raise SmokeSubsetError("written smoke H5AD lacks Geneformer fields")
    if np.any(np.asarray(reloaded.X.sum(axis=1)).ravel() != reloaded.obs["n_counts"]):
        raise SmokeSubsetError("written smoke H5AD n_counts do not rederive")

    content_identity = {
        "data_sha256": hashlib.sha256(matrix.data.tobytes()).hexdigest(),
        "indices_sha256": hashlib.sha256(matrix.indices.tobytes()).hexdigest(),
        "indptr_sha256": hashlib.sha256(matrix.indptr.tobytes()).hexdigest(),
        "cell_ids_sha256": canonical_sha256(list(map(str, value.obs_names))),
        "ensembl_ids_sha256": canonical_sha256(ensembl_ids),
        "n_counts_sha256": hashlib.sha256(n_counts.tobytes()).hexdigest(),
    }
    return {
        "shape": list(value.shape),
        "nnz": int(value.X.nnz),
        "content_identity": content_identity,
    }


def build(
    *, atlas: Path, contract_lock: Path, cell_manifest: Path, output: Path
) -> None:
    output.mkdir(parents=True, exist_ok=False)
    contract = json.loads(contract_lock.read_text(encoding="utf-8"))
    expected_manifest_hash = contract["manifest_sha256"][cell_manifest.name]
    if sha256_file(cell_manifest) != expected_manifest_hash:
        raise SmokeSubsetError("cell manifest differs from the production contract")
    if atlas.resolve(strict=True).as_posix() != Path(
        contract["atlas_realpath"]
    ).resolve(strict=True).as_posix():
        raise SmokeSubsetError("atlas path differs from the production contract")
    if atlas.stat().st_size != int(contract["atlas_size_bytes"]):
        raise SmokeSubsetError("atlas size differs from the production contract")

    selected = select_cells(cell_manifest)
    ontology_payload = {
        "schema_version": "masld-bench-cell-ontology-v1",
        "ontology_id": "broad_liver_five_v1",
        "source_field": "cell_type",
        "classes": {key: list(value) for key, value in sorted(ONTOLOGY.items())},
        "selection_outcomes_used": False,
        "sealed_outcomes_used": False,
    }
    ontology_payload["ontology_sha256"] = canonical_sha256(ontology_payload)
    ontology_path = output / "broad_liver_five_v1.json"
    write_json_exclusive(ontology_path, ontology_payload)

    selection_path = output / "selection.tsv"
    selection_fields = (
        "source_row_index",
        "cell_id",
        "source_cell_id",
        "library_id",
        "donor_id",
        "dataset",
        "source_cell_type",
        "broad_label",
    )
    with selection_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=selection_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    output_h5ad = output / "resource_atlas_geneformer_smoke_1000.h5ad"
    matrix_summary = materialize_subset(
        atlas=atlas, selected=selected, output_h5ad=output_h5ad
    )
    by_class: dict[str, int] = {}
    by_dataset: dict[str, int] = {}
    donors_by_class: dict[str, set[str]] = {}
    for row in selected:
        broad = str(row["broad_label"])
        dataset = str(row["dataset"])
        donor = str(row["donor_id"])
        by_class[broad] = by_class.get(broad, 0) + 1
        by_dataset[dataset] = by_dataset.get(dataset, 0) + 1
        donors_by_class.setdefault(broad, set()).add(donor)
    if set(by_class) != set(ONTOLOGY) or set(by_class.values()) != {CELLS_PER_CLASS}:
        raise SmokeSubsetError("written selection is not exactly class-balanced")

    manifest = {
        "schema_version": "masld-bench-geneformer-smoke-subset-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "dataset_id": "resource_atlas_current",
        "subset_id": "resource_atlas_geneformer_smoke_1000_v1",
        "source": {
            "atlas_path": atlas.resolve(strict=True).as_posix(),
            "atlas_size_bytes": atlas.stat().st_size,
            "contract_lock_path": contract_lock.resolve(strict=True).as_posix(),
            "contract_lock_file_sha256": sha256_file(contract_lock),
            "contract_lock_id": contract["lock_sha256"],
            "raw_counts_sha256": contract["raw_counts_sha256"],
            "cell_manifest_path": cell_manifest.resolve(strict=True).as_posix(),
            "cell_manifest_sha256": expected_manifest_hash,
        },
        "selection": {
            "cell_budget": CELL_BUDGET,
            "cells_per_class": CELLS_PER_CLASS,
            "selection_seed": SELECTION_SEED,
            "selection_policy": "donor_round_robin_then_sha256_cell_priority",
            "sealed_outcomes_used": False,
            "selection_outcomes_used": False,
            "counts_by_class": dict(sorted(by_class.items())),
            "donors_by_class": {
                key: len(value) for key, value in sorted(donors_by_class.items())
            },
            "counts_by_dataset": dict(sorted(by_dataset.items())),
        },
        "artifacts": {
            "h5ad": {
                "path": output_h5ad.resolve().as_posix(),
                "sha256": sha256_file(output_h5ad),
                "size_bytes": output_h5ad.stat().st_size,
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
        "matrix": matrix_summary,
        "geneformer_input_contract": {
            "matrix": "X_raw_unselected_counts",
            "obs_n_counts": True,
            "var_ensembl_id": True,
            "ready": True,
        },
    }
    manifest["manifest_sha256"] = canonical_sha256(manifest)
    write_json_exclusive(output / "smoke_subset_manifest.json", manifest)
    freeze_output_tree(
        output,
        subset_id=manifest["subset_id"],
        manifest_id=manifest["manifest_sha256"],
    )


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--atlas", required=True, type=Path)
    value.add_argument("--contract-lock", required=True, type=Path)
    value.add_argument("--cell-manifest", required=True, type=Path)
    value.add_argument("--output", required=True, type=Path)
    return value


def main() -> int:
    args = parser().parse_args()
    build(
        atlas=args.atlas,
        contract_lock=args.contract_lock,
        cell_manifest=args.cell_manifest,
        output=args.output,
    )
    print((args.output / "smoke_subset_manifest.json").resolve().as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
