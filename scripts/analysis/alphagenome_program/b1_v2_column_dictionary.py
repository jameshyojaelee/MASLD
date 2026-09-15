#!/usr/bin/env python3
"""Correction 1, completed: one deposited row per column of the response layer, naming its track panel.

Renaming the exported quantile columns to `*_primary_liver_*` fixes the columns this package computes. Three
columns are carried verbatim from Track 0 and cannot be renamed without breaking the correspondence with
`direction_variant_level.tsv.gz` and the P6d dossier: `pred_allele1_liver_rna`, `predicted_direction` and
`direction_concordant`. Those three rest on the Track 0 liver panel, which is `primary_liver` (5 RNA tracks:
3 adult, 1 child, 1 embryonic) POOLED WITH `hepatocyte` = CL:0000182 (2 tracks, embryonic) -- a different
panel from every other liver column in the same row. This dictionary says so, per column, in the deposit.

Life stages are read from the Atlas build's own raw/scorer_metadata, never inferred from a track name.

Usage: b1_v2_column_dictionary.py <v2 deposit directory>
"""

from __future__ import annotations

import csv
import gzip
import json
import pathlib
import sys

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
TRACK0 = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z"
META = TRACK0 / "raw/scorer_metadata"

TRACK0_DIRECTION_TRACKS = [
    ("RNA_SEQ", "UBERON:0001114 total RNA-seq"), ("RNA_SEQ", "UBERON:0001115 total RNA-seq"),
    ("RNA_SEQ", "UBERON:0002107 polyA plus RNA-seq"), ("RNA_SEQ", "UBERON:0002107 total RNA-seq"),
    ("RNA_SEQ", "UBERON:0001114 gtex Liver polyA plus RNA-seq"),
    ("RNA_SEQ", "CL:0000182 total RNA-seq"), ("RNA_SEQ", "CL:0000182 polyA plus RNA-seq")]

# column -> (what it is, unit, sign convention or the column that states it)
NOTES = {
    "block_1mb": ("resampling unit for every variant-level test in this program", "chrom~floor(pos/1e6)", ""),
    "block_1mb_recovered": ("block_1mb, or the block of the recovered hg38 identity when the variant was "
                            "deferred and has no variant_uid", "chrom~floor(pos/1e6)", ""),
    "analysis_block": ("Track 0 locus block; a DIFFERENT unit from block_1mb, empty where the signal has "
                       "none", "Track 0 block id", ""),
    "weight": ("posterior weight within the signal, ranked over the full posterior and never renormalised; "
               "exported rows are floored at 0.01 plus the rank-1 variant, so exported weights do not sum "
               "to total_mass", "posterior probability", ""),
    "ase_gse281367_mean_log2_alt_over_ref": ("measured allelic imbalance over heterozygous donors",
                                             "log2(ALT/REF)", "ase_sign_refers_to"),
    "ase_gse244832_mean_log2_alt_over_ref": ("measured allelic imbalance over heterozygous donors",
                                             "log2(ALT/REF)", "ase_sign_refers_to"),
    "ase_gse281367_n_het_donors": ("heterozygous donors behind the mean; donors are the unit, the SITE is "
                                   "the unit of any median over this column", "donors", ""),
    "ase_gse244832_n_het_donors": ("heterozygous donors behind the mean; donors are the unit, the SITE is "
                                   "the unit of any median over this column", "donors", ""),
    "measured_atac_logFC_gse244832": ("donor-level differential accessibility of the overlapping consensus "
                                      "peak", "log2 fold change",
                                      "measured_atac_logFC_gse244832_sign_refers_to"),
    "measured_atac_logFC_gse281367": ("donor-level differential accessibility of the overlapping consensus "
                                      "peak", "log2 fold change",
                                      "measured_atac_logFC_gse281367_sign_refers_to"),
    "gwas_beta_allele1": ("GWAS effect carried from Track 0", "trait-specific",
                          "direction_sign_refers_to"),
    "eqtl_beta_allele1": ("liver eQTL effect carried from Track 0", "eQTL beta",
                          "direction_sign_refers_to"),
    "model_api_rna_log2": ("hosted model-API rescue effect; own scale, never pooled with an Atlas quantile",
                           "log2(ALT/REF)", "model_api_sign_refers_to"),
    "model_api_splice_log2": ("hosted model-API rescue effect", "log2(ALT/REF)", "model_api_sign_refers_to"),
    "model_api_atac_log2": ("hosted model-API rescue effect", "log2(ALT/REF)", "model_api_sign_refers_to"),
    "model_api_dnase_log2": ("hosted model-API rescue effect", "log2(ALT/REF)", "model_api_sign_refers_to"),
    "model_api_h3k27ac_log2": ("hosted model-API rescue effect", "log2(ALT/REF)", "model_api_sign_refers_to"),
    "avi_quantile": ("Atlas AVI feature score quantile; a prediction, not a liver panel", "quantile", ""),
    "avi_raw": ("Atlas AVI feature score", "model units", ""),
}


def life_stage_map():
    out = {}
    for p in sorted(META.glob("*.track_metadata.tsv")):
        scorer = p.name.split(".track_metadata.tsv")[0]
        with open(p) as h:
            for r in csv.DictReader(h, delimiter="\t"):
                out[(scorer, r["name"])] = r.get("biosample_life_stage", "")
    return out


def main() -> None:
    dep = pathlib.Path(sys.argv[1]).resolve()
    with gzip.open(dep / "response_layer_variants.tsv.gz", "rt") as h:
        header = next(csv.reader(h, delimiter="\t"))
    panel = {}
    with open(dep / "track_panel_composition.tsv") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            panel[r["exported_column"]] = r
    stages = life_stage_map()
    d_stage = [stages[k] for k in TRACK0_DIRECTION_TRACKS]
    track0_panel = ("Track 0 liver RNA median over 7 tracks: primary_liver 5 ("
                    f"{d_stage[:5].count('adult')} adult, {d_stage[:5].count('child')} child, "
                    f"{d_stage[:5].count('embryonic')} embryonic) POOLED WITH hepatocyte CL:0000182 2 ("
                    f"{d_stage[5:].count('embryonic')} embryonic). NOT the panel behind the "
                    "*_primary_liver_* columns in the same row.")

    rows = []
    for c in header:
        p = panel.get(c)
        note, unit, sign = NOTES.get(c, ("", "", ""))
        if c in ("pred_allele1_liver_rna", "predicted_direction", "direction_concordant",
                 "measured_direction"):
            tp = track0_panel
            note = note or ("carried verbatim from Track 0 direction_variant_level.tsv.gz; the column name "
                            "is Track 0's own and is kept so the join stays checkable")
            sign = sign or "direction_sign_refers_to"
        elif p is not None:
            tp = (f"q|primary_liver| group, {p['n_tracks']} tracks ({p['life_stages']}): {p['tracks']}")
        else:
            tp = ""
        rows.append({"column": c, "what_it_is": note, "unit": unit, "sign_convention_column": sign,
                     "track_panel": tp})
    dest = dep / "column_dictionary.tsv"
    if dest.exists():
        raise SystemExit(f"refusing to overwrite {dest}")
    with dest.open("w") as h:
        w = csv.DictWriter(h, fieldnames=["column", "what_it_is", "unit", "sign_convention_column",
                                          "track_panel"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(json.dumps({"n_columns": len(rows),
                      "n_columns_with_a_track_panel": sum(1 for r in rows if r["track_panel"]),
                      "track0_direction_panel": track0_panel}, indent=1))


if __name__ == "__main__":
    main()
