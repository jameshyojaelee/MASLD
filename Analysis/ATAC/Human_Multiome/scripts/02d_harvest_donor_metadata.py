#!/usr/bin/env python
"""
02d_harvest_donor_metadata.py

Build a curated donor metadata table for the GSE244832 (Wang/Guo 2024 Genome
Med) multiome ATAC pipeline. Bridges the ATAC donor IDs (MM_*) to the RNA
donor IDs (JB_*) to SRR accessions to the canonical F-stage axis
(`F_stage_augmented`) sourced from donor_metadata_extended.tsv.

Inputs:
  - data/GSE244832/metadata/donor_pairing.csv
      MM_*/JB_* sample pairing, ATAC SRR, RNA SRRs (;-joined), condition
      (NORMAL/MASL/MASH) — already curated from GEO Sample_title parsing.
  - data/GSE244832/metadata/ena_runinfo.tsv
      ENA SRR <-> JB_* mapping (audit only).
  - Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv
      Canonical scRNA donor metadata with disease_stage_coarse,
      F_stage_documented, F_stage_inferred, F_stage_source, F_stage_augmented,
      age, sex, nas_score at SRR-row granularity.

Note on GEO characteristics:
  GSE244832 GEO Sample_characteristics_ch1 only carries `tissue: liver`
  (no F-stage, NAS, age, or sex). This is documented in
  docs/dataset_labeling_and_harmonization.md. F-stage harvest therefore
  flows through donor_metadata_extended.tsv (which sources Andrews
  documented F-stage anchors + scVI projection across all 117 GSE244832
  SRRs via Script 343m).

For each MM_* donor we:
  1. Look up the paired JB_* RNA donor's SRR list.
  2. Pull all matching SRR rows from donor_metadata_extended.tsv.
  3. Aggregate F_stage_augmented / F_stage_documented / F_stage_inferred /
     disease_stage_coarse / age / sex / nas_score across those SRRs
     (modal aggregation — same donor across multiplexed runs should
     return identical values; any inconsistency is flagged).
  4. Emit one row per ATAC donor.

Output:
  Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv
  Schema: donor_id | donor_id_atac | donor_id_rna | atac_srr | rna_srrs |
          condition | disease_stage_coarse | fibrosis_stage_documented |
          F_stage_inferred | F_stage_source | F_stage_augmented |
          nas_score | age | sex | n_rna_srrs_matched | aggregation_notes

Run: ~5 sec; no GEO fetch (GEO already cached under data/GSE244832/metadata/).
"""

from __future__ import annotations

import os
from collections import Counter
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)

DONOR_PAIRING = PROJECT_ROOT / "data/GSE244832/metadata/donor_pairing.csv"
ENA_RUNINFO = PROJECT_ROOT / "data/GSE244832/metadata/ena_runinfo.tsv"
DONOR_META_EXT = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
)
OUT_DIR = PROJECT_ROOT / "Analysis/ATAC/Human_Multiome/metadata"
OUT_TSV = OUT_DIR / "donor_metadata_curated.tsv"


def _mode_or_none(values):
    """Return modal non-null value (and a note if non-unanimous), else None."""
    clean = [v for v in values if pd.notna(v) and str(v) not in {"", "nan", "NA"}]
    if not clean:
        return None, "all_na"
    counts = Counter(clean)
    top_val, top_n = counts.most_common(1)[0]
    if len(counts) == 1:
        return top_val, ""
    return top_val, f"non_unanimous({dict(counts)})"


def _coerce_int(x):
    if x is None:
        return None
    try:
        return int(float(str(x).lstrip("F").lstrip("f")))
    except Exception:
        return None


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load pairing
    pairing = pd.read_csv(DONOR_PAIRING)
    print(f"[input] donor_pairing: {len(pairing)} ATAC donors")

    # 2. Load donor_metadata_extended; restrict to GSE244832 SRR rows
    ext = pd.read_csv(DONOR_META_EXT, sep="\t", low_memory=False)
    ext_w = ext[ext["dataset"] == "GSE244832"].copy()
    print(f"[input] donor_metadata_extended GSE244832 rows: {len(ext_w)}")
    print(f"[input] unique SRRs: {ext_w['sample'].nunique()}")

    # 3. Sanity audit against ena_runinfo (which SRR <-> JB_* it says)
    ena = pd.read_csv(ENA_RUNINFO, sep="\t")
    ena["jb_id"] = ena["sample_title"].str.extract(r"(JB_\d+)")
    ena_jb_to_srr = (
        ena.dropna(subset=["jb_id"])
        .groupby("jb_id")["run_accession"]
        .apply(list)
        .to_dict()
    )

    # 4. Per ATAC donor, aggregate axes across JB_* RNA SRRs.
    out_rows = []
    for _, r in pairing.iterrows():
        mm = r["atac_sample"]
        jb = r["rna_sample"]
        atac_srr = r["atac_srr"]
        rna_srrs_decl = [s.strip() for s in str(r["rna_srrs"]).split(";") if s.strip()]
        cond = r["condition"]

        # Audit: ENA-resolved SRRs for this JB
        ena_resolved = set(ena_jb_to_srr.get(jb, []))
        decl_set = set(rna_srrs_decl)
        # Use the union of declared + ENA-resolved for maximum coverage
        # (donor_pairing.csv is the curated truth, but ENA gives a fallback).
        srr_pool = list(decl_set | ena_resolved)

        matched = ext_w[ext_w["sample"].isin(srr_pool)]
        n_matched = len(matched)

        # Aggregate
        notes = []
        if n_matched == 0:
            notes.append(f"no_srrs_matched_pool={srr_pool}")
            stage_coarse, _ = (None, "no_match")
            f_doc, _ = (None, "no_match")
            f_inf, _ = (None, "no_match")
            f_src, _ = (None, "no_match")
            f_aug, _ = (None, "no_match")
            nas, _ = (None, "no_match")
            age, _ = (None, "no_match")
            sex, _ = (None, "no_match")
        else:
            stage_coarse, n_sc = _mode_or_none(matched["disease_stage_coarse"].tolist())
            f_doc, n_fd = _mode_or_none(matched["F_stage_documented"].tolist())
            f_inf, n_fi = _mode_or_none(matched["F_stage_inferred"].tolist())
            f_src, n_fs = _mode_or_none(matched["F_stage_source"].tolist())
            f_aug, n_fa = _mode_or_none(matched["F_stage_augmented"].tolist())
            nas, n_n = _mode_or_none(matched["nas_score"].tolist())
            age, n_a = _mode_or_none(matched["age"].tolist())
            sex, n_s = _mode_or_none(matched["sex"].tolist())
            for tag, txt in (
                ("stage", n_sc),
                ("fdoc", n_fd),
                ("finf", n_fi),
                ("fsrc", n_fs),
                ("faug", n_fa),
                ("nas", n_n),
                ("age", n_a),
                ("sex", n_s),
            ):
                if txt:
                    notes.append(f"{tag}:{txt}")

        out_rows.append(
            {
                "donor_id": r["donor_id"],
                "donor_id_atac": mm,
                "donor_id_rna": jb,
                "atac_srr": atac_srr,
                "rna_srrs": ";".join(sorted(rna_srrs_decl)),
                "condition": cond,
                "disease_stage_coarse": stage_coarse,
                "fibrosis_stage_documented": _coerce_int(f_doc),
                "F_stage_inferred": _coerce_int(f_inf),
                "F_stage_source": f_src,
                "F_stage_augmented": _coerce_int(f_aug),
                "nas_score": nas,
                "age": age,
                "sex": sex,
                "n_rna_srrs_matched": n_matched,
                "aggregation_notes": ";".join(notes) if notes else "",
            }
        )

    out = pd.DataFrame(out_rows)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\n[output] {OUT_TSV}")
    print(f"[output] {len(out)} ATAC donors written")

    # 5. Report
    n_placeholder = (
        (out["condition"].isna())
        | (out["condition"].astype(str).str.upper() == "PLACEHOLDER")
        | (out["condition"].astype(str) == "")
    ).sum()
    print(f"\n[verify] PLACEHOLDER/NA condition rows: {n_placeholder} (expected 0)")

    print("\n[F_stage_augmented distribution among 18 ATAC donors]")
    print(out["F_stage_augmented"].value_counts(dropna=False).sort_index().to_string())

    print("\n[disease_stage_coarse distribution]")
    print(out["disease_stage_coarse"].value_counts(dropna=False).sort_index().to_string())

    print("\n[condition distribution]")
    print(out["condition"].value_counts(dropna=False).sort_index().to_string())

    print("\n[F_stage_source distribution]")
    print(out["F_stage_source"].value_counts(dropna=False).sort_index().to_string())

    # 6. Flag any aggregation issues
    issues = out[out["aggregation_notes"].astype(str) != ""]
    if len(issues) > 0:
        print(
            f"\n[note] {len(issues)} donors have aggregation notes (see column "
            "'aggregation_notes' in TSV)"
        )

    # 7. Halt if anything still PLACEHOLDER
    if n_placeholder > 0:
        print(
            "\n[HALT] Some donors still PLACEHOLDER. Manual transcription from "
            "Guo et al. 2024 Genome Med supplementary Table 1 required for:"
        )
        missing = out[
            (out["condition"].isna())
            | (out["condition"].astype(str).str.upper() == "PLACEHOLDER")
            | (out["condition"].astype(str) == "")
        ]
        print(missing[["donor_id", "donor_id_atac", "donor_id_rna"]].to_string(index=False))


if __name__ == "__main__":
    main()
