#!/usr/bin/env python3
"""Freeze outcome-independent gene length, biotype, and cell-specificity covariates."""

from __future__ import annotations

import csv
import gzip
import re

from bridge_common import CANDIDATE_ROOT, PROJECT_ROOT, require_validated_seal, write_tsv


GTF = "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"


def main() -> None:
    require_validated_seal()
    atlas: dict[str, dict[str, str]] = {}
    with (PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            atlas[row["human_symbol"]] = {
                "atlas_gene_biotype": row.get("gene_biotype", ""),
                "sc_is_celltype_specific": row.get("sc_is_celltype_specific", ""),
                "bulk_tstat": row.get("bulk_tstat", ""),
            }
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    with gzip.open(GTF, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            match = re.search(r'gene_name "([^"]+)"', fields[8])
            if not match or match.group(1) in seen:
                continue
            symbol = match.group(1); seen.add(symbol)
            type_match = re.search(r'gene_type "([^"]+)"', fields[8])
            extra = atlas.get(symbol, {})
            rows.append({
                "gene_symbol": symbol,
                "gene_length_bp": int(fields[4]) - int(fields[3]) + 1,
                "log_gene_length": __import__("math").log(int(fields[4]) - int(fields[3]) + 1),
                "gene_biotype": type_match.group(1) if type_match else extra.get("atlas_gene_biotype", ""),
                "sc_is_celltype_specific": extra.get("sc_is_celltype_specific", ""),
                "canonical_bulk_tstat": extra.get("bulk_tstat", ""),
                "source": "GENCODE_v49_gene_span_and_frozen_multi_evidence_atlas",
            })
    if len(rows) < 70000:
        raise RuntimeError(f"Unexpectedly small GENCODE gene universe: {len(rows)}")
    write_tsv(CANDIDATE_ROOT / "protein/gene_covariates.tsv", rows, list(rows[0]))
    print("GENE_COVARIATES_COMPLETE")


if __name__ == "__main__":
    main()
