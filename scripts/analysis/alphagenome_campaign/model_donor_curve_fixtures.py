#!/usr/bin/env python3
"""Cut the donor learning-curve panel into one fixture per training size.

The earlier donor-context result rested on 128 training regions, and a null there
bounds the training support as much as the model. The curve panel selected 8,192
training regions and a fixed 1,024-region held panel outcome-independently, with
nested training sets, so the question becomes whether context prediction improves
as training regions grow while the evaluated regions stay identical.

Rather than teach the fitting code a new notion of training subset, this emits
three fixtures in the shape it already reads. Each holds the same 1,024 held
regions and the training regions for its size, so every curve point evaluates the
same output regions in the same held donors, and the only thing that moves is how
many training regions the model saw.

The gene axis is carried over from the original 128-region fixture after checking
the gene identities match exactly, and the reverse-complement-invariant
dinucleotide comparator is recomputed with the original producer's own function.
No target value, variance, significance or previous prediction quality is read.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

import review_donor as D

PROJ = Path(__file__).resolve().parents[3]
BASE = PROJ / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model"
ORIGINAL_FIXTURE = BASE / "donor_fixture_21773589"
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
SIZES = (128, 1024, 8192)


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Compute-node execution required")
    args.out.mkdir(parents=True, exist_ok=False)
    reg = pd.read_csv(args.panel / "regions.tsv", sep="\t")
    with np.load(args.panel / "paired.npz", allow_pickle=True) as f:
        ids = f["participant_ids"].astype(str)
        fold = f["donor_fold"]
        rna = f["rna_log2cpm_full_library"]
        target = f["target_log2cpm_full_library"]
        gene_ids = f["gene_ids"].astype(str)
        region_key = f["region_key"].astype(str)
    if not np.array_equal(region_key, reg.region_key.to_numpy(dtype=str)):
        raise ValueError("Panel region order disagrees with its own paired arrays")
    if target.shape != (len(ids), len(reg)):
        raise ValueError("Panel target axis disagrees with its region count")

    with np.load(ORIGINAL_FIXTURE / "paired.npz", allow_pickle=True) as f:
        if not np.array_equal(f["gene_ids"].astype(str), gene_ids):
            raise ValueError("Gene axis differs from the original fixture; cannot carry coordinates")
        gene_chrom = f["gene_chrom"].astype(str)
        gene_start = f["gene_start"]
        gene_end = f["gene_end"]

    features = None
    if args.sequence_features:
        with np.load(args.sequence_features, allow_pickle=False) as f:
            if not np.array_equal(f["region_key"].astype(str), region_key):
                raise ValueError("Sequence feature region identities differ from the panel")
            features = {"region_key": f["region_key"], "strand_features": f["strand_features"],
                        "mean_strands": f["mean_strands"]}

    genome = pysam.FastaFile(FASTA)
    try:
        dinucleotides = np.asarray(
            [D.sequence_features(genome.fetch(r.chrom, int(r.start0), int(r.end0)).upper())
             for r in reg.itertuples()], dtype=float)
    finally:
        genome.close()
    if dinucleotides.shape[0] != len(reg) or not np.isfinite(dinucleotides).all():
        raise ValueError("Dinucleotide comparator is incomplete for the panel")

    held = reg.region_role.eq("held").to_numpy()
    if int(held.sum()) != args.held_regions:
        raise ValueError(f"Expected {args.held_regions} held regions, found {int(held.sum())}")
    written = []
    for size in SIZES:
        column = f"training_{size}"
        train = reg[column].to_numpy(dtype=bool)
        if int(train.sum()) != size:
            raise ValueError(f"{column} marks {int(train.sum())} regions, not {size}")
        if (train & held).any():
            raise ValueError("A held region is marked for training")
        keep = np.flatnonzero(train | held)
        folder = args.out / f"size_{size}"
        folder.mkdir()
        arrays = {"participant_ids": ids, "donor_fold": fold,
                  "rna_log2cpm_full_library": rna,
                  "target_log2cpm_full_library": target[:, keep],
                  "dinucleotide_features": dinucleotides[keep],
                  "gene_ids": gene_ids, "gene_chrom": gene_chrom,
                  "gene_start": gene_start, "gene_end": gene_end}
        np.savez_compressed(folder / "paired.npz", **arrays)
        subset = reg.iloc[keep].reset_index(drop=True)
        subset.to_csv(folder / "regions.tsv", sep="\t", index=False)
        if features is not None:
            np.savez_compressed(folder / "features.npz",
                                region_key=features["region_key"][keep],
                                strand_features=features["strand_features"][keep],
                                mean_strands=features["mean_strands"][keep])
        # The held panel must be byte-identical across curve points, or the curve
        # compares evaluation sets rather than training sizes.
        held_keys = subset.loc[subset.region_role.eq("held"), "region_key"].to_numpy(dtype=str)
        written.append({"training_regions": size, "regions": int(len(keep)),
                        "held_regions": int(len(held_keys)),
                        "inner_train": int(subset.inner_region_role.eq("inner_train").sum()),
                        "inner_valid": int(subset.inner_region_role.eq("inner_valid").sum()),
                        "input_components": int(subset.input_component.nunique()),
                        "held_key_sha256": hashlib.sha256(
                            "\n".join(sorted(held_keys)).encode()).hexdigest(),
                        "fixture": str(folder.relative_to(PROJ))})
    if len({w["held_key_sha256"] for w in written}) != 1:
        raise ValueError("Held panel differs between curve points")

    receipt = {"status": "curve_fixtures_written",
               "panel": str(args.panel.relative_to(PROJ)),
               "participants": int(len(ids)), "RNA_features": int(rna.shape[1]),
               "training_sizes": list(SIZES), "held_regions": int(held.sum()),
               "nested_training_sets": all(
                   set(reg.index[reg[f"training_{a}"]]) <= set(reg.index[reg[f"training_{b}"]])
                   for a, b in zip(SIZES, SIZES[1:])),
               "identical_held_panel_across_points": True,
               "sequence_features": str(args.sequence_features) if args.sequence_features else None,
               "dinucleotide_comparator": "review_donor.sequence_features, the original producer's own function",
               "gene_axis": "carried from the original 128-region fixture after exact gene-identity check",
               "outcomes_read": False, "target_values_used_for_selection": False,
               "fixtures": written,
               "panel_regions_sha256": hashlib.sha256(
                   (args.panel / "regions.tsv").read_bytes()).hexdigest()}
    (args.out / "curve_fixtures_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--sequence-features", type=Path)
    parser.add_argument("--held-regions", type=int, default=1024)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
