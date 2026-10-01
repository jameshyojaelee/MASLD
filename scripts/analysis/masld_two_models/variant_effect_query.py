#!/usr/bin/env python3
"""Query a measured Currin variant-peak pair and its matched development predictions.

The adapter was trained on association-selected variant-peak leads and does not
use an arbitrary target interval as an input. A request outside a source pair
therefore fails instead of inventing target-specific prediction.
"""
import argparse
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
LABELS = ROOT/"GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz"
VARIANT = re.compile(r"^GRCh38:chr(?:[1-9]|1[0-9]|2[0-2]|X|Y):[1-9][0-9]*:[ACGT]:[ACGT]$")
INTERVAL = re.compile(r"^(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)):([0-9]+)-([0-9]+)$")


def query(variant, interval, assay, predictions):
    if assay != "caQTL":
        raise ValueError("Only the measured Currin caQTL endpoint is supported")
    if not VARIANT.fullmatch(variant):
        raise ValueError("Expected GRCh38:chrN:1-based-position:REF:ALT for a single SNV")
    target = INTERVAL.fullmatch(interval)
    if target is None or int(target.group(2)) >= int(target.group(3)):
        raise ValueError("Expected chrN:start-end with positive source interval bounds")
    # Identify the row using metadata alone before parsing any effect fields.
    labels = pd.read_csv(LABELS, sep="\t", usecols=["lead_variant_id", "heldout_fold"])
    source_key = variant.removeprefix("GRCh38:")
    matched = labels.loc[labels.lead_variant_id.eq(source_key)]
    if len(matched) != 1:
        raise ValueError("Variant is not one unique measured Currin lead; arbitrary variant-target effects are unsupported")
    if int(matched.iloc[0].heldout_fold) not in range(1, 5):
        raise ValueError("Spent fold 0 is unavailable for further queries")
    source_index = int(matched.index[0])
    source = pd.read_csv(LABELS, sep="\t",
                         skiprows=lambda line: line != 0 and line != source_index+1,
                         usecols=["lead_variant_id", "chr", "pos_hg38", "ref", "alt",
                                  "peak_id", "peak_start_hg38", "peak_stop_hg38",
                                  "beta_alt", "heldout_fold"]).iloc[0]
    source_interval = f"{source.chr}:{int(source.peak_start_hg38)}-{int(source.peak_stop_hg38)}"
    if interval != source_interval:
        raise ValueError(f"Target interval is not this variant's measured source peak ({source_interval})")
    frame = pd.read_csv(predictions, sep="\t")
    if frame.key.duplicated().any() or frame.fold.eq(0).any():
        raise ValueError("Prediction table has duplicate identity or spent fold 0")
    row = frame.loc[frame.key.eq(source_key.removeprefix("chr"))]
    if len(row) != 1:
        raise ValueError("Pair lacks all matched development predictions; unsupported")
    row = row.iloc[0]
    if int(row.fold) != int(source.heldout_fold) or row.ref != source.ref or row.alt != source.alt:
        raise ValueError("Source fold or allele identity differs")
    if abs(float(row.beta_alt)-float(source.beta_alt)) > 1e-6:
        raise ValueError("Measured signed beta differs from prediction table")
    return {
        "variant": variant, "target_source_peak": interval, "source_peak_id": source.peak_id,
        "variant_coordinate_system": "GRCh38, 1-based SNV position",
        "target_coordinate_system": "GRCh38 source BED peak bounds, 0-based start and exclusive end; unchanged from the Currin source BED",
        "assay": "Currin liver caQTL", "allele_contrast": "ALT minus REF",
        "effect_unit": "source FastQTL additive ALT dosage beta",
        "effect_interpretation": "Association coefficient per additional ALT allele, not an isolated causal allele-editing effect",
        "measured_beta_alt": float(source.beta_alt),
        "predicted_beta_alt": {
            "native_ATAC_2048bp": float(row.native2k),
            "fixed_head_2048bp": float(row.fixed_head),
            "adapter_2048bp_five_seed": float(row.adapter),
            "native_ATAC_1048576bp": float(row.native1m),
            "native_plus_adapter_combination": float(row.combined)},
        "prediction_state": "held-chromosome-fold development for this exact measured pair; fold 0 excluded",
        "target_specificity": "The native models read out 501 bp centered on the variant; the adapter pools a 384-bp central embedding interval within 2-kb sequence. Neither uses source peak coordinates as a target input. The peak identifies the measured association, not an independently modeled target interval.",
        "target_gene": None, "target_gene_evidence": "unsupported",
        "single_variant_predictive_interval": None,
        "uncertainty_reason": "Only aggregate paired intervals are available; no calibrated individual interval",
        "recommended_experiment": "Test the allele substitution at the variant and measure accessibility of the source peak; assess target-gene links separately",
        "clinical_or_MASLD_effect": "unsupported",
    }


def main(args):
    if args.out.exists():
        raise FileExistsError(args.out)
    result = query(args.variant, args.target_interval, args.assay, args.predictions)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", required=True)
    parser.add_argument("--target-interval", required=True)
    parser.add_argument("--assay", default="caQTL")
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
