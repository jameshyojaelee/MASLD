#!/usr/bin/env python3
"""Create source-replicate-by-lineage raw-count pseudobulks from the deposited OS-HLO object."""

from __future__ import annotations

import gzip
import json
from collections import Counter

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from public_functional_common import CANDIDATE_ROOT, read_tsv, require_sealed, sha256_file, write_tsv


def main() -> None:
    seal = require_sealed()
    gates = read_tsv(CANDIDATE_ROOT / "source_gate_status.tsv")
    gate = [row for row in gates if row["dataset_id"] == "GSE207889"]
    if len(gate) != 1 or gate[0]["inference_authorized"] != "true":
        raise RuntimeError("GSE207889 source gate does not authorize pseudobulk construction")
    sample_rows = [
        row for row in read_tsv(CANDIDATE_ROOT / "sample_manifest.tsv")
        if row["dataset_id"] == "GSE207889" and row["include_in_inference"] == "true"
    ]
    sample_map = {row["sample_id"]: row for row in sample_rows}
    if len(sample_map) != 12:
        raise RuntimeError(f"Expected 12 prespecified OS-HLO source samples, found {len(sample_map)}")

    source = CANDIDATE_ROOT / "sources/GSE207889/GSE207889_os.h5ad"
    data = ad.read_h5ad(source)
    required_obs = {"sample", "detailed_condition", "celltype_category"}
    if not required_obs.issubset(data.obs.columns) or data.raw is None:
        raise RuntimeError("Deposited OS-HLO source fields/raw matrix drifted")
    obs = data.obs.copy()
    selected = obs["sample"].astype(str).isin(sample_map)
    obs = obs.loc[selected].copy()
    raw = data.raw.X[selected.values]
    genes = np.asarray(data.raw.var_names.astype(str))
    if len(set(genes)) != len(genes):
        raise RuntimeError("Raw OS-HLO gene symbols are not unique")
    if not sparse.issparse(raw):
        raw = sparse.csr_matrix(raw)
    raw = raw.tocsr()
    sampled_values = raw.data[: min(1_000_000, raw.data.size)]
    if sampled_values.size and not np.allclose(sampled_values, np.rint(sampled_values)):
        raise RuntimeError("OS-HLO raw.X is not integer-like count data")

    group_keys = list(zip(obs["sample"].astype(str), obs["celltype_category"].astype(str)))
    unique_groups = sorted(set(group_keys))
    pseudobulk = np.zeros((len(genes), len(unique_groups)), dtype=np.int64)
    manifest = []
    total_by_sample = Counter(obs["sample"].astype(str))
    for column, (sample_id, lineage) in enumerate(unique_groups):
        mask = np.asarray([(left == sample_id and right == lineage) for left, right in group_keys])
        counts = np.asarray(raw[mask].sum(axis=0)).ravel()
        if not np.allclose(counts, np.rint(counts)):
            raise RuntimeError(f"Noninteger pseudobulk counts for {sample_id}/{lineage}")
        pseudobulk[:, column] = np.rint(counts).astype(np.int64)
        source_row = sample_map[sample_id]
        n_cells = int(mask.sum())
        manifest.append(
            {
                "pseudobulk_id": f"{sample_id}__{lineage}",
                "sample_id": sample_id,
                "biological_unit_id": source_row["biological_unit_id"],
                "condition": source_row["condition"],
                "treatment": source_row["treatment"],
                "replicate": source_row["replicate"],
                "lineage": lineage,
                "n_cells": n_cells,
                "sample_total_cells": int(total_by_sample[sample_id]),
                "cell_fraction": n_cells / total_by_sample[sample_id],
                "include_in_inference": "true",
                "biological_n_semantics": "source_replicate; cells are technical observations",
                "specification_sha256": seal["specification_sha256"],
            }
        )
    if len(unique_groups) != 48 or {row["lineage"] for row in manifest} != {"Hepatocytes", "Cholangiocytes", "Fibroblasts", "Hepatic stellate cells"}:
        raise RuntimeError("Expected 12 samples x 4 deposited lineages")

    output_root = CANDIDATE_ROOT / "preprocessed/GSE207889"
    output_root.mkdir(parents=True, exist_ok=True)
    count_path = output_root / "pseudobulk_counts.tsv.gz"
    with gzip.open(count_path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("gene_symbol\t" + "\t".join(row["pseudobulk_id"] for row in manifest) + "\n")
        for index, gene in enumerate(genes):
            handle.write(gene + "\t" + "\t".join(map(str, pseudobulk[index])) + "\n")
    write_tsv(output_root / "pseudobulk_manifest.tsv", manifest, list(manifest[0]))
    gene_audit = [
        {"source": "raw.var_names", "n_genes": len(genes), "n_unique": len(set(genes)), "count_matrix_sha256": sha256_file(count_path), "source_h5ad_sha256": sha256_file(source)}
    ]
    write_tsv(output_root / "gene_mapping_audit.tsv", gene_audit, list(gene_audit[0]))
    print(json.dumps({"n_cells": int(selected.sum()), "n_pseudobulks": len(manifest), "n_genes": len(genes), "lineages": sorted({row["lineage"] for row in manifest})}, indent=2))


if __name__ == "__main__":
    main()
