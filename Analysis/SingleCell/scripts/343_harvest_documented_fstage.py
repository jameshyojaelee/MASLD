#!/usr/bin/env python
"""
343_harvest_documented_fstage.py

Harvest documented donor-level fibrosis stage (F0-F4) from GEO series_matrix
files and published supplementary tables for the 7 scRNA datasets in the
MASLD atlas. NO classifier projection — only labels that were explicitly
published.

Inputs:
  - donor_metadata.tsv (sample, dataset, ...) — drives the donor roster
  - unified_metadata.csv (bulk-RNA-seq metadata, where some scRNA donors
    overlap and carry fibrosis_stage)
  - GEO series_matrix.txt.gz (downloaded on demand)

Outputs:
  - Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv
    Columns: sample, dataset, F_stage_documented (0-4 or NA), source

Note: this script has two parts.
  (1) Auto-harvest from GEO series_matrix characteristics_ch1 fields.
  (2) Fallback overlay from bulk unified_metadata.csv (sample-ID-matched).

Coverage estimate after run:
  - GSE244832 (Wang) -> if SRA / supp has fibrosis_stage, ~80-100 of 117 donors
  - GSE202379 (Andrews) -> partial (cirrhosis-heavy, F-stage rare)
  - Liver_Atlas (Guilliams) -> some donors via supp Table S1
  - Others -> mostly NA
"""

from __future__ import annotations
import gzip
import io
import os
import re
import sys
import urllib.request
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

DONOR_META = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
BULK_META  = PROJECT_ROOT / "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"
OUT_DIR    = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
GEO_CACHE  = OUT_DIR / "geo_cache"
OUT_TSV    = OUT_DIR / "donor_fstage_documented.tsv"

GEO_CACHE.mkdir(parents=True, exist_ok=True)

# GEO Series accessions -> series_matrix URLs
GEO_URLS = {
    "GSE136103": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/matrix/GSE136103_series_matrix.txt.gz",
    "GSE174748": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174748/matrix/GSE174748_series_matrix.txt.gz",
    "GSE185477": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185477/matrix/GSE185477_series_matrix.txt.gz",
    "GSE189600": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189600/matrix/GSE189600_series_matrix.txt.gz",
    "GSE202379": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202379/matrix/GSE202379_series_matrix.txt.gz",
    "GSE244832": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244832/matrix/GSE244832_series_matrix.txt.gz",
    # Liver_Atlas = Guilliams GSE192742 (mouse + human); use GSE192740 (human only)
    "Liver_Atlas": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE192nnn/GSE192740/matrix/GSE192740_series_matrix.txt.gz",
}

# Regex patterns to extract fibrosis stage from GEO characteristics strings.
# Matches: "fibrosis_stage: F2", "Kleiner fibrosis stage: F3", "saf score: S2A3F3"
FSTAGE_PATTERNS = [
    re.compile(r"saf\s*score[:\s]*S\d+A\d+F(\d)", re.I),
    re.compile(r"fibrosis[ _]*stage[\":\s]*F?(\d)", re.I),
    re.compile(r"fibrosis[\":\s]+F?(\d)", re.I),
    re.compile(r"Kleiner[ _]*fibrosis[\":\s]*F?(\d)", re.I),
    re.compile(r"\bf[ _]*stage[\":\s]*F?(\d)", re.I),
    re.compile(r"\bF(\d)\b"),  # bare "F3" anywhere (last-resort, low precedence)
]

# "saf score: healthy control" -> F=0; "saf score: end stage" -> F=4
SAF_LABEL_MAP = {
    "healthy control": 0,
    "end stage": 4,
}

def fetch_geo_index(acc: str, url_template: str) -> list[str]:
    """List actual *_series_matrix.txt.gz files in the GEO matrix dir.

    GEO sometimes returns 404 for the bare GSE<acc>_series_matrix.txt.gz
    pattern when the series spans multiple platforms; the actual files
    are named GSE<acc>-GPL<plat>_series_matrix.txt.gz.
    """
    # Drop the filename, fetch the parent dir listing.
    dir_url = url_template.rsplit("/", 1)[0] + "/"
    try:
        html = urllib.request.urlopen(dir_url, timeout=30).read().decode(errors="ignore")
    except Exception as e:
        print(f"[GEO] WARN: dir listing failed for {acc}: {e}")
        return []
    files = re.findall(r'href="([^"]+_series_matrix\.txt\.gz)"', html)
    return [dir_url + f for f in files]

def download_geo(acc: str, url: str) -> list[Path]:
    """Download series_matrix.txt.gz(s) to cache. Returns list of local paths."""
    locals_out = []
    # First try the canonical URL
    local = GEO_CACHE / f"{acc}_series_matrix.txt.gz"
    if local.exists() and local.stat().st_size > 1000:
        return [local]
    try:
        print(f"[GEO] downloading {acc} from {url}")
        urllib.request.urlretrieve(url, local)
        if local.stat().st_size > 1000:
            return [local]
    except Exception as e:
        print(f"[GEO] WARN: canonical URL failed for {acc}: {e}; trying dir listing")

    # Fall back to dir-listing for multi-platform series
    urls = fetch_geo_index(acc, url)
    for u in urls:
        fname = u.rsplit("/", 1)[-1]
        local_i = GEO_CACHE / fname
        if local_i.exists() and local_i.stat().st_size > 1000:
            locals_out.append(local_i)
            continue
        try:
            print(f"[GEO] downloading {fname}")
            urllib.request.urlretrieve(u, local_i)
            locals_out.append(local_i)
        except Exception as e:
            print(f"[GEO] WARN: failed to download {fname}: {e}")
    return locals_out

def parse_geo_characteristics(path: Path) -> pd.DataFrame:
    """Parse GSM accession -> dict of characteristic field -> value."""
    rows = []
    with gzip.open(path, "rt", errors="ignore") as fh:
        lines = fh.readlines()

    # !Sample_geo_accession provides the GSM order
    gsms = []
    for L in lines:
        if L.startswith("!Sample_geo_accession"):
            parts = L.strip().split("\t")[1:]
            gsms = [p.strip('"') for p in parts]
            break
    if not gsms:
        return pd.DataFrame()

    char_lines = [L for L in lines if L.startswith("!Sample_characteristics_ch1")]
    char_rows = []
    for L in char_lines:
        parts = L.strip().split("\t")[1:]
        char_rows.append([p.strip('"') for p in parts])

    if not char_rows:
        return pd.DataFrame({"sample": gsms})

    # Each row is a separate characteristic across all samples
    n_samples = len(gsms)
    df = pd.DataFrame({"sample": gsms})
    for i, row in enumerate(char_rows):
        if len(row) != n_samples:
            continue
        df[f"char_{i}"] = row

    return df

def extract_fstage(row: pd.Series) -> int | None:
    """Run regex + label-map over all characteristic strings; return 0-4 or None."""
    # Label-map pass (e.g. "saf score: healthy control" -> 0)
    for col in row.index:
        if not col.startswith("char_"):
            continue
        val = str(row[col]).lower()
        for label, fval in SAF_LABEL_MAP.items():
            if label in val:
                return fval
    # Regex pass
    for col in row.index:
        if not col.startswith("char_"):
            continue
        val = str(row[col])
        for pat in FSTAGE_PATTERNS:
            m = pat.search(val)
            if m:
                try:
                    f = int(m.group(1))
                    if 0 <= f <= 4:
                        return f
                except ValueError:
                    pass
    return None

def parse_geo_bioproject(path: Path) -> str | None:
    """Extract the BioProject (PRJNA*) accession from a GEO series_matrix."""
    with gzip.open(path, "rt", errors="ignore") as fh:
        for L in fh:
            if L.startswith("!Series_relation") and "BioProject" in L:
                m = re.search(r"(PRJ[A-Z]{2}\d+)", L)
                if m:
                    return m.group(1)
            if L.startswith("!Series_relation") and "study" in L:
                m = re.search(r"(SRP\d+)", L)
                if m:
                    return m.group(1)
    return None

def resolve_gsm_to_srr(prj: str) -> dict[str, str]:
    """Use ENA portal to map GSM -> SRR for a BioProject. Returns {gsm: srr}."""
    url = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={prj}&result=read_run"
        f"&fields=run_accession,sample_alias&format=tsv"
    )
    try:
        body = urllib.request.urlopen(url, timeout=60).read().decode()
    except Exception as e:
        print(f"[ENA] WARN: {prj} resolution failed: {e}")
        return {}
    lines = body.strip().split("\n")
    if len(lines) < 2:
        return {}
    out = {}
    for L in lines[1:]:
        parts = L.split("\t")
        if len(parts) >= 2:
            srr, gsm = parts[0], parts[1]
            if gsm.startswith("GSM") and srr.startswith("SRR"):
                out[gsm] = srr
    return out

def harvest_geo() -> pd.DataFrame:
    """Return DataFrame [sample, dataset, F_stage_documented, source].

    Two-pass:
      (1) Parse GSM-level F-stages from series_matrix characteristics.
      (2) Resolve each dataset's BioProject -> ENA SRR mapping; emit BOTH
          GSM- and SRR-keyed rows so the merge with donor_metadata.tsv
          (which uses SRR) works regardless of which ID convention the
          donor table uses.
    """
    rows = []
    for acc, url in GEO_URLS.items():
        locals_ = download_geo(acc, url)
        if not locals_:
            continue
        # Pick one platform file for the BioProject (they share PRJ)
        prj = parse_geo_bioproject(locals_[0])
        srr_map = resolve_gsm_to_srr(prj) if prj else {}
        print(f"[ENA] {acc} -> PRJ={prj}, resolved {len(srr_map)} GSM->SRR")

        for local in locals_:
            chars = parse_geo_characteristics(local)
            if chars.empty:
                continue
            chars["F_stage_documented"] = chars.apply(extract_fstage, axis=1)
            for _, r in chars.iterrows():
                gsm = r["sample"]
                fstage = r["F_stage_documented"]
                src = f"GEO_series_matrix:{local.name}"
                # Emit GSM-keyed row
                rows.append({"sample": gsm, "dataset": acc,
                             "F_stage_documented": fstage, "source": src})
                # Emit SRR-keyed row when ENA resolved this GSM
                srr = srr_map.get(gsm)
                if srr:
                    rows.append({"sample": srr, "dataset": acc,
                                 "F_stage_documented": fstage,
                                 "source": f"{src}+ENA:{gsm}"})

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("F_stage_documented", na_position="last")
        df = df.drop_duplicates(subset=["sample", "dataset"], keep="first")
    return df

def overlay_bulk_fstage(harvest: pd.DataFrame, donor: pd.DataFrame) -> pd.DataFrame:
    """For donors with NA F-stage, try bulk unified_metadata.csv sample match."""
    if not BULK_META.exists():
        return harvest
    bulk = pd.read_csv(BULK_META)
    # Try several plausible F-stage column names
    fcols = [c for c in bulk.columns if c.lower() in
             ("fibrosis_stage", "fibrosis", "kleiner_fibrosis", "f_stage")]
    if not fcols:
        return harvest
    fcol = fcols[0]
    bulk_lookup = bulk[["sample_id" if "sample_id" in bulk.columns else "sample", fcol]].copy()
    bulk_lookup.columns = ["sample", "F_stage_bulk"]

    # Coerce numeric
    def coerce(x):
        try:
            s = str(x).upper().lstrip("F").strip()
            return int(float(s)) if s.replace(".", "").isdigit() else None
        except Exception:
            return None
    bulk_lookup["F_stage_bulk"] = bulk_lookup["F_stage_bulk"].apply(coerce)

    # Merge into the donor list, fill NA from bulk
    donor_keys = donor[["sample", "dataset"]].drop_duplicates()
    merged = donor_keys.merge(harvest, on=["sample", "dataset"], how="left")
    merged = merged.merge(bulk_lookup, on="sample", how="left")
    fill_mask = merged["F_stage_documented"].isna() & merged["F_stage_bulk"].notna()
    merged.loc[fill_mask, "F_stage_documented"] = merged.loc[fill_mask, "F_stage_bulk"]
    merged.loc[fill_mask, "source"] = "bulk_unified_metadata"
    return merged[["sample", "dataset", "F_stage_documented", "source"]]

def main():
    donor = pd.read_csv(DONOR_META, sep="\t")
    print(f"[input] {len(donor)} scRNA donors across {donor['dataset'].nunique()} datasets")

    harvest = harvest_geo()
    print(f"[GEO]   {len(harvest)} GSM entries with characteristics parsed")
    print(f"[GEO]   {harvest['F_stage_documented'].notna().sum()} carry a documented F-stage")

    final = overlay_bulk_fstage(harvest, donor)
    # Coerce to nullable Int (F0-F4) for tidy output
    final["F_stage_documented"] = pd.to_numeric(final["F_stage_documented"], errors="coerce").astype("Int64")

    print()
    print("Coverage per dataset:")
    cov = final.groupby("dataset").agg(
        n_donors=("sample", "count"),
        n_with_fstage=("F_stage_documented", lambda s: s.notna().sum()),
    )
    print(cov.to_string())
    print()
    print("F-stage distribution (all datasets):")
    print(final["F_stage_documented"].value_counts(dropna=False).sort_index().to_string())

    final.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\n[output] {OUT_TSV}")

if __name__ == "__main__":
    main()
