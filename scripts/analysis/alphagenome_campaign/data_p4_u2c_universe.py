#!/usr/bin/env python3
"""Build the U2c genetic-overlap region universe, before any of its scores is read.

Step 41 left U2c empty: it ran at 15:22 on 2026-09-09 and logged that Track 0's
variant-peak overlap table was not present; that table was written at 23:27 the
same day. The original universe is not reopened. This registers the overlap
regions as a separate universe with its own prespecification
(data_p4_u2c_prespec.json), retrieved separately and tested as its own family.

Admission uses lineage and peak coordinate only. The overlap table also carries
hosted Atlas variant-effect columns; reading those to admit a region and then
testing a model score over that region would be circular, so this refuses to run
if those columns are loaded.

Output: <out>/tables/saturation_universe.tsv.gz in the schema step 42 consumes,
plus universe_receipt.json. Point AGA_OUT_ROOT at <out> to retrieve it.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
CTX = ROOT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
OVERLAPS = ROOT / "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/tables/variant_peak_overlaps.tsv.gz"
ORIGINAL = ROOT / "GWAS/finemapping/results/alphagenome_atlas/p4-saturation-20260909T191917Z/tables/saturation_universe.tsv.gz"
PRESPEC = Path(__file__).with_name("data_p4_u2c_prespec.json")
UNIVERSE_NAME = "U2c_genetic_overlap"
ADMISSION_COLUMNS = ["lineage", "peak", "queried", "signal_uid"]
FORBIDDEN = ["atlas_atac_liver_median_raw", "atlas_atac_liver_median_quantile",
             "atlas_dnase_liver_median_raw"]
COLUMNS = ["universe", "region_key", "chrom", "start0", "end", "lineage", "u2a_da_state",
           "u2b_sig_either", "u2c_genetic_overlap", "u2d_program_linked", "width",
           "n_overlapping_signals", "any_queried_variant"]


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Compute-node execution required")
    args.out.mkdir(parents=True, exist_ok=False)
    prespec = json.loads(PRESPEC.read_text())
    if prespec["family"] != "F6_genetic_overlap":
        raise ValueError("Prespecification is not the registered genetic-overlap family")

    overlaps = pd.read_csv(OVERLAPS, sep="\t", usecols=ADMISSION_COLUMNS, dtype=str)
    if any(c in overlaps.columns for c in FORBIDDEN):
        raise ValueError("Hosted variant-effect columns must not be loaded for admission")

    manifest = pd.read_csv(CTX / "consensus_peak_manifest.tsv", sep="\t", dtype=str)
    manifest = manifest.loc[~manifest.blacklist_overlap.str.upper().eq("TRUE")]
    eligible = set(manifest.lineage + "|" + manifest.chrom + ":" + manifest.start0 + "-" + manifest.end)
    overlaps = overlaps.loc[(overlaps.lineage + "|" + overlaps.peak).isin(eligible)]

    original = pd.read_csv(ORIGINAL, sep="\t", usecols=["universe", "region_key"])
    registered = set(original.loc[original.universe.eq("U2_snatac"), "region_key"])
    if len(registered) != 6930:
        raise ValueError("Original U2 coordinate universe changed")

    grouped = overlaps.groupby("peak")
    stats = pd.DataFrame({
        "n_overlapping_signals": grouped.signal_uid.nunique(),
        "any_queried_variant": grouped.queried.apply(lambda s: bool((s.str.lower() == "true").any())),
        "lineage": grouped.lineage.apply(lambda s: ";".join(sorted(set(s)))),
    })
    new = sorted(set(stats.index) - registered)

    rows = []
    for key in new:
        chrom, span = key.split(":")
        start0, end = (int(x) for x in span.split("-"))
        rows.append({"universe": UNIVERSE_NAME, "region_key": key, "chrom": chrom, "start0": start0,
                     "end": end, "lineage": stats.at[key, "lineage"], "u2a_da_state": "",
                     "u2b_sig_either": False, "u2c_genetic_overlap": True,
                     "u2d_program_linked": False, "width": end - start0,
                     "n_overlapping_signals": int(stats.at[key, "n_overlapping_signals"]),
                     "any_queried_variant": bool(stats.at[key, "any_queried_variant"])})
    frame = pd.DataFrame(rows, columns=COLUMNS)
    if frame.region_key.duplicated().any() or len(frame) != len(new):
        raise ValueError("Duplicate or dropped overlap coordinate")
    tables = args.out / "tables"
    tables.mkdir(parents=True)
    path = tables / "saturation_universe.tsv.gz"
    frame.to_csv(path, sep="\t", index=False)

    receipt = {
        "status": "u2c_universe_defined_before_any_of_its_scores_was_read",
        "universe": UNIVERSE_NAME,
        "overlap_coordinates_total": int(len(stats)),
        "already_registered_in_U2": int(len(set(stats.index) & registered)),
        "regions_to_retrieve": len(frame),
        "widths": sorted(frame.width.unique().tolist()),
        "total_bp": int(frame.width.sum()),
        "regions_with_a_queried_variant": int(frame.any_queried_variant.sum()),
        "admission_columns": ADMISSION_COLUMNS,
        "hosted_variant_effect_columns_not_loaded": FORBIDDEN,
        "original_universe_unchanged": True,
        "no_saturation_values_read": True,
        "prespecification": str(PRESPEC.relative_to(ROOT)),
        "prespecification_sha256": hashlib.sha256(PRESPEC.read_bytes()).hexdigest(),
        "universe_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "overlaps_sha256": hashlib.sha256(OVERLAPS.read_bytes()).hexdigest(),
    }
    (args.out / "universe_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
