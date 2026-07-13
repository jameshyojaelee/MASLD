#!/usr/bin/env python3
"""Build leakage-free, chromosome-fold-specific Currin SNV score inputs."""

from __future__ import annotations

import csv
import gzip
import json
import os
from pathlib import Path

import pandas as pd
from pyfaidx import Fasta

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
BASE = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2"
TRUTH = ROOT / "GWAS/finemapping/results/seqfunc/chrombpnet_caqtl/v2_truth"
FASTA = Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")


def main() -> None:
    out = BASE / "currin_calibration"
    inputs = out / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    pos = pd.read_csv(TRUTH / "currin_positive_1kb.tsv.gz", sep="\t", dtype={"chr": str})
    neg = pd.read_csv(TRUTH / "currin_source_background.tsv.gz", sep="\t", dtype={"chr": str})
    source_counts = {"positive_rows": len(pos), "background_rows": len(neg)}
    pos = pos.loc[pos["variant_class"].eq("SNV")].copy()
    neg = neg.loc[neg["variant_class"].eq("SNV")].copy()
    pos["label"] = 1
    neg["label"] = 0
    keep = ["variant_id_hg38", "chr", "pos_hg38", "ref", "alt", "label", "block_1mb"]
    positive_extra = ["beta_alt", "q_value", "p_nominal", "peak_id"]
    pos = pos[keep + positive_extra]
    for col in positive_extra:
        neg[col] = pd.NA
    neg = neg[keep + positive_extra]
    combined = pd.concat([pos, neg], ignore_index=True)
    combined["chr"] = combined["chr"].astype(str)
    combined = combined.loc[combined["chr"].isin([f"chr{i}" for i in range(1, 23)])]
    # Positive precedence prevents a source-significant lead from appearing in
    # the negative class if source tables ever overlap.
    combined = combined.sort_values(["variant_id_hg38", "label"], ascending=[True, False])
    n_pre_dedup = len(combined)
    combined = combined.drop_duplicates("variant_id_hg38", keep="first").copy()

    genome = Fasta(str(FASTA), as_raw=True, sequence_always_upper=True)
    ref_match = []
    for row in combined.itertuples(index=False):
        ref_match.append(genome[row.chr][int(row.pos_hg38) - 1:int(row.pos_hg38)] == row.ref)
    combined["reference_match"] = ref_match
    n_ref_mismatch = int((~combined["reference_match"]).sum())
    combined = combined.loc[combined["reference_match"]].copy()

    fold_map = {}
    with (BASE / "fold_manifest.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for chrom in row["test_chromosomes"].split(","):
                if chrom in fold_map:
                    raise SystemExit(f"chromosome assigned twice: {chrom}")
                fold_map[chrom] = int(row["fold"])
    if set(fold_map) != {f"chr{i}" for i in range(1, 23)}:
        raise SystemExit("fold map does not cover chr1-chr22 exactly once")
    combined["heldout_fold"] = combined["chr"].map(fold_map)
    combined.to_csv(out / "currin_scoring_metadata.tsv.gz", sep="\t", index=False, compression="gzip")
    fold_counts = {}
    for fold in range(5):
        subset = combined.loc[combined["heldout_fold"].eq(fold), ["chr", "pos_hg38", "ref", "alt", "variant_id_hg38"]]
        subset.to_csv(inputs / f"fold{fold}.variants.tsv", sep="\t", index=False, header=False)
        fold_counts[str(fold)] = {"total": len(subset), "positive": int(combined.loc[combined["heldout_fold"].eq(fold), "label"].sum())}
        if len(subset) == 0:
            raise SystemExit(f"fold{fold} Currin input is empty")
    report = {
        "source_counts": source_counts,
        "primary_truth": "official Currin GRCh38 1-kb lead caQTLs",
        "negative_truth": "official source-tested non-caPeak beta-adjusted-p>0.5 background",
        "variant_scope": "biallelic SNVs only; indels excluded because allele normalization differs from variant-scorer schema",
        "rows_before_variant_dedup": n_pre_dedup,
        "duplicate_variant_rows_removed": n_pre_dedup - len(combined) - n_ref_mismatch,
        "reference_mismatches_removed": n_ref_mismatch,
        "scored_unique_snvs": len(combined),
        "positive_unique_snvs": int(combined["label"].sum()),
        "background_unique_snvs": int((combined["label"] == 0).sum()),
        "fold_counts": fold_counts,
        "fold_rule": "score each variant exactly once with the model whose test group contains its chromosome",
    }
    with (out / "input_contract.json").open("w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

