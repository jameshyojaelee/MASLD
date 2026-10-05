#!/usr/bin/env python3
"""C1b acquisition check for the Model C external test (PRESPEC_C Amendment 6).
No image is downloaded and no RNA value is read.

Counts TCGA-LIHC Solid Tissue Normal (sample type 11) H&E slides in the NCI
Imaging Data Commons whose sample barcode (TCGA-XX-XXXX-11A) also has a column in
the UCSC Xena GDC-hub STAR TPM matrix. Only the header line of that matrix is read.
Age, sex and tissue source site come from the Xena GDC clinical table. One slide per
donor (the largest series). Writes a seal file with no sealed donor so
c3_tiles_embed.py can be reused unchanged.
"""
import argparse
import gzip
import json
from pathlib import Path

import pandas as pd
import requests
from idc_index import IDCClient

HUB = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"


def tpm_header():
    with requests.get(f"{HUB}/TCGA-LIHC.star_tpm.tsv.gz", stream=True, timeout=120) as r:
        r.raise_for_status()
        with gzip.GzipFile(fileobj=r.raw) as g:
            return g.readline().decode().rstrip("\n").split("\t")[1:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    client = IDCClient()
    client.fetch_index("sm_index")
    idx = client.index[["SeriesInstanceUID", "collection_id", "PatientID", "series_size_MB"]]
    sm = client.sm_index.merge(idx, on="SeriesInstanceUID")
    sm = sm[sm["collection_id"] == "tcga_lihc"].copy()
    cid = sm["ContainerIdentifier"].astype(str)
    sm["sample"] = cid.str[:16]
    sm["sample_type"] = cid.str[13:15]
    sm["slide_code"] = cid.str[17:]
    normal = sm[sm["sample_type"] == "11"].copy()
    normal.to_csv(out / "idc_tcga_lihc_normal_series.tsv", sep="\t", index=False)

    rna = pd.Series(tpm_header())
    rna_normal = set(rna[rna.str[13:15] == "11"])
    matched = normal[normal["sample"].isin(rna_normal)].copy()
    matched = matched.sort_values(["PatientID", "series_size_MB"], ascending=[True, False]).drop_duplicates("PatientID")
    matched["tss"] = matched["PatientID"].str.split("-").str[1]

    clin = pd.read_csv(f"{HUB}/TCGA-LIHC.clinical.tsv.gz", sep="\t", low_memory=False)
    clin = clin[["sample", "gender.demographic", "age_at_index.demographic", "preservation_method.samples",
                 "is_ffpe.samples"]].rename(columns={"gender.demographic": "gender", "age_at_index.demographic": "age",
                                                     "preservation_method.samples": "preservation",
                                                     "is_ffpe.samples": "is_ffpe"})
    matched = matched.merge(clin, on="sample", how="left")
    matched.to_csv(out / "liver_slides_matched.tsv", sep="\t", index=False)
    pd.DataFrame({"SUBJID": matched["PatientID"], "sealed": False}).to_csv(out / "donor_seal_none.tsv", sep="\t", index=False)

    summary = {
        "match_rule": "same TCGA sample barcode (TCGA-XX-XXXX-11A) for the slide and the RNA-seq column; one slide per donor (largest series)",
        "idc_tcga_lihc_series": int(len(sm)),
        "idc_normal_series": int(len(normal)),
        "idc_normal_donors": int(normal["PatientID"].nunique()),
        "rna_normal_samples": len(rna_normal),
        "donors_matched": int(matched["PatientID"].nunique()),
        "donors_per_tissue_source_site": matched["tss"].value_counts().to_dict(),
        "slide_code": matched["slide_code"].value_counts().to_dict(),
        "embedding_medium": matched["embeddingMedium_CodeMeaning"].astype(str).value_counts().to_dict(),
        "pixel_spacing_mm": matched["min_PixelSpacing_2sf"].astype(str).value_counts().to_dict(),
        "objective_power": matched["ObjectiveLensPower"].astype(str).value_counts().to_dict(),
        "age_missing": int(matched["age"].isna().sum()),
        "gender": matched["gender"].value_counts(dropna=False).to_dict(),
        "matched_size_GB": float(matched["series_size_MB"].sum() / 1024),
        "idc_version": str(client.get_idc_version()),
        "rna_source": f"{HUB}/TCGA-LIHC.star_tpm.tsv.gz (header line only)",
        "non_tumor_histology": "not in the Xena clinical table; GDC API returned 503 on 2026-10-03",
    }
    (out / "c1b_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
