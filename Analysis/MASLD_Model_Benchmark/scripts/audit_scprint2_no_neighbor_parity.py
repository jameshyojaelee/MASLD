#!/usr/bin/env python3
"""Audit scPRINT-2 no-neighbor inputs against released collation semantics."""

from __future__ import annotations

import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy import sparse

from masld_bench.artifacts import freeze_tree, write_json_exclusive


SOURCE_SHA256 = "d6fa563b11570f37d25d9a62dd4626e8f63f0624adea873bb102a550fedc8854"
FIXTURE_ARTIFACTS_SHA256 = "9c2ecbecb610ea35687085ad10188e6c3f6418469cee229c9d83f6ebc7e43d9f"
ALIGNED_SHA256 = "462d7cdb79d630022c2ad1462ba7f0666300f70b4f31d9f9f4c6513a2a3b35e1"
ROWS = 1000
FORBIDDEN = ("label", "fibrosis", "nas", "mash", "masld", "outcome", "stage")


class Scprint2ParityAuditError(RuntimeError):
    """Raised when a bound input or released-source requirement differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def row_totals(matrix: Any) -> np.ndarray:
    values = sparse.csr_matrix(matrix).sum(axis=1)
    return np.asarray(values).ravel().astype(np.float64)


def embedder_defaults(path: Path) -> dict[str, Any]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == "Embedder":
            initializer = next(
                (
                    item
                    for item in node.body
                    if isinstance(item, ast.FunctionDef) and item.name == "__init__"
                ),
                None,
            )
            if initializer is None:
                break
            names = [argument.arg for argument in initializer.args.args]
            defaults = [ast.literal_eval(item) for item in initializer.args.defaults]
            mapping = dict(zip(names[-len(defaults) :], defaults, strict=True))
            return {
                "how": mapping["how"],
                "max_len": mapping["max_len"],
                "pred_embedding": mapping["pred_embedding"],
                "use_knn": mapping["use_knn"],
            }
    raise Scprint2ParityAuditError("released Embedder defaults were not found")


def build(
    *, source: Path, fixture_root: Path, cell_embedder_source: Path, output: Path
) -> dict[str, Any]:
    if output.exists():
        raise Scprint2ParityAuditError("audit output already exists")
    if digest(source) != SOURCE_SHA256:
        raise Scprint2ParityAuditError("registered source H5AD differs")
    if digest(fixture_root / "ARTIFACTS.json") != FIXTURE_ARTIFACTS_SHA256:
        raise Scprint2ParityAuditError("registered scPRINT-2 fixture differs")
    aligned_path = fixture_root / "fixture/aligned_counts.h5ad"
    observed_aligned_sha = digest(aligned_path)
    if observed_aligned_sha != ALIGNED_SHA256:
        raise Scprint2ParityAuditError("aligned count matrix differs")

    source_adata = ad.read_h5ad(source)
    aligned_adata = ad.read_h5ad(aligned_path)
    allowed_source_obs = {"row_id", "donor_id", "dataset_id", "outer_fold", "assay"}
    if set(source_adata.obs.columns) != allowed_source_obs:
        raise Scprint2ParityAuditError("source observation firewall differs")
    if any(
        any(token in column.lower() for token in FORBIDDEN)
        for column in source_adata.obs.columns
    ):
        raise Scprint2ParityAuditError("source contains a forbidden outcome field")
    if (
        source_adata.n_obs != ROWS
        or aligned_adata.n_obs != ROWS
        or source_adata.obs_names.astype(str).tolist()
        != aligned_adata.obs_names.astype(str).tolist()
    ):
        raise Scprint2ParityAuditError("source and aligned row identities differ")
    if not sparse.issparse(source_adata.X) or not sparse.issparse(aligned_adata.X):
        raise Scprint2ParityAuditError("count inputs are not sparse")

    source_totals = row_totals(source_adata.X)
    aligned_totals = row_totals(aligned_adata.X)
    lost = source_totals - aligned_totals
    if np.any(lost < 0) or np.any(source_totals <= 0):
        raise Scprint2ParityAuditError("aligned counts exceed source or source is empty")
    loss_fraction = lost / source_totals
    defaults = embedder_defaults(cell_embedder_source)
    if defaults != {
        "how": "random expr",
        "max_len": 2000,
        "pred_embedding": ["all"],
        "use_knn": True,
    }:
        raise Scprint2ParityAuditError("released Embedder defaults differ")

    mismatch_rows = int(np.count_nonzero(lost))
    result = {
        "schema_version": "masld-bench-scprint2-no-neighbor-parity-audit-v1",
        "status": "fail_depth_parity" if mismatch_rows else "pass_depth_parity",
        "rows": ROWS,
        "source_features": source_adata.n_vars,
        "aligned_features": aligned_adata.n_vars,
        "rows_with_depth_mismatch": mismatch_rows,
        "rows_with_depth_mismatch_fraction": mismatch_rows / ROWS,
        "source_total_counts": int(source_totals.sum()),
        "aligned_total_counts": int(aligned_totals.sum()),
        "counts_removed_before_collator_depth": int(lost.sum()),
        "counts_removed_fraction": float(lost.sum() / source_totals.sum()),
        "per_row_removed_fraction_quantiles": {
            str(quantile): float(np.quantile(loss_fraction, quantile))
            for quantile in (0.0, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0)
        },
        "released_embedder_defaults": defaults,
        "active_lane_max_len": 3200,
        "active_lane_neighbor_policy": "no_neighbor",
        "interpretation": (
            "Released Collator records depth before valid_genes filtering. The active "
            "fixture filters to checkpoint genes before collation, so every nonzero "
            "removed-count row receives a non-native req_depth. The 3200-token path is "
            "checkpoint-max random-expressed, not the released Embedder default of 2000."
        ),
        "evaluation_allowed": mismatch_rows == 0,
        "evaluation_label_columns_read": [],
        "sealed_outcomes_read": False,
        "source_h5ad_sha256": SOURCE_SHA256,
        "fixture_artifacts_sha256": FIXTURE_ARTIFACTS_SHA256,
        "aligned_h5ad_sha256": observed_aligned_sha,
        "released_cell_embedder_source_sha256": digest(cell_embedder_source),
    }
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "parity_audit.json", result)
    freeze_tree(
        output,
        {
            "artifact_class": "scprint2_no_neighbor_exact_source_parity_audit",
            "evaluation_labels_read": False,
            "rows": ROWS,
            "sealed_outcomes_read": False,
            "status": result["status"],
        },
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--fixture-root", required=True, type=Path)
    parser.add_argument("--cell-embedder-source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = build(
        source=args.source,
        fixture_root=args.fixture_root,
        cell_embedder_source=args.cell_embedder_source,
        output=args.output,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
