#!/usr/bin/env python3
"""C1 acquisition check for Model C (PRESPEC_C section 3). No image is downloaded
and no RNA value is read: counts GTEx liver slides in the NCI Imaging Data Commons
that match a GTEx v8 liver RNA sample, their total size, and the donor seal.

Liver RNA samples come from the open-access v8 sample attributes (SMTSD == Liver).
A slide matches when its GTEx tissue-sample ID (GTEX-XXXXX-NNNN) equals the RNA
sample's first three ID fields. The 20% seal ranks donors by
sha256("MASLD-C-2026-09-23" + SUBJID).
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
from idc_index import IDCClient

TISSUE_ID = re.compile(r"(GTEX-[A-Z0-9]+-\d{4})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-attributes", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    attr = pd.read_csv(a.sample_attributes, sep="\t", low_memory=False)
    liver = attr[(attr["SMTSD"] == "Liver") & attr["SAMPID"].str.startswith("GTEX-")
                 & (attr["SMAFRZE"] == "RNASEQ")].copy()
    liver["tissue_id"] = liver["SAMPID"].str.split("-").str[:3].str.join("-")
    liver["SUBJID"] = liver["SAMPID"].str.split("-").str[:2].str.join("-")

    client = IDCClient()
    client.fetch_index("sm_index")  # slide-microscopy index: ContainerIdentifier = GTEx tissue ID
    idx = client.index[["SeriesInstanceUID", "collection_id", "PatientID", "series_size_MB"]]
    gtex = client.sm_index.merge(idx, on="SeriesInstanceUID")
    gtex = gtex[gtex["collection_id"].str.lower() == "gtex"].copy()
    gtex["tissue_id"] = gtex["ContainerIdentifier"].astype(str).str.extract(TISSUE_ID)[0]
    keep_cols = ["PatientID", "SeriesInstanceUID", "ContainerIdentifier", "tissue_id",
                 "primaryAnatomicStructure_CodeMeaning", "min_PixelSpacing_2sf", "ObjectiveLensPower",
                 "max_TotalPixelMatrixColumns", "max_TotalPixelMatrixRows", "series_size_MB"]
    gtex[keep_cols].to_csv(out / "idc_gtex_series.tsv", sep="\t", index=False)

    # Liver slides and liver RNA aliquots carry different tissue codes (0911/0914 matched 1 of 251
    # by tissue ID), so match at the donor: the unit of inference is the donor anyway.
    liver_site = gtex[gtex["primaryAnatomicStructure_CodeMeaning"].astype(str).str.contains("iver", na=False)]
    matched = liver_site[liver_site["PatientID"].isin(set(liver["SUBJID"]))].copy()
    matched = matched.sort_values(["PatientID", "series_size_MB"], ascending=[True, False]).drop_duplicates("PatientID")
    donors = sorted(matched["PatientID"])
    rank = sorted(donors, key=lambda s: hashlib.sha256(("MASLD-C-2026-09-23" + s).encode()).hexdigest())
    n_seal = round(0.2 * len(rank))
    sealed = pd.DataFrame({"SUBJID": rank, "sealed": [i < n_seal for i in range(len(rank))]})
    sealed.to_csv(out / "donor_seal.tsv", sep="\t", index=False)
    matched.to_csv(out / "liver_slides_matched.tsv", sep="\t", index=False)

    summary = {
        "match_rule": "donor-level: one liver-site slide per donor (largest series) among donors with v8 liver RNA-seq",
        "idc_liver_site_series": int(len(liver_site)),
        "matched_pixel_spacing_mm": matched["min_PixelSpacing_2sf"].astype(str).value_counts().head(5).to_dict(),
        "matched_objective_power": matched["ObjectiveLensPower"].astype(str).value_counts().head(5).to_dict(),
        "idc_gtex_series": int(len(gtex)),
        "idc_gtex_series_with_tissue_id": int(gtex["tissue_id"].notna().sum()),
        "liver_rna_samples_v8_rnaseq": int(len(liver)),
        "liver_slides_matched_to_rna": int(len(matched)),
        "liver_donors_matched": len(donors),
        "sealed_donors": n_seal,
        "matched_size_GB": float(matched["series_size_MB"].sum() / 1024),
        "idc_version": str(getattr(client, "get_idc_version", lambda: "unknown")()),
    }
    (out / "c1_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
