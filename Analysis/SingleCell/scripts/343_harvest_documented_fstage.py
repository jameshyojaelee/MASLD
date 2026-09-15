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
    Columns: sample, dataset, F_stage_documented (0-4 or NA), source,
             saf_S, saf_A, saf_F (deposited SAF grades, or NA),
             f_stage_source (which pass produced F_stage_documented)

    saf_S / saf_A were previously DISCARDED: the SAF regex captured only the
    F group. They are now kept.

    2026-08-30 -- LABEL-MAP PRECEDENCE CORRECTED, F_stage_documented CHANGES.
    extract_fstage() used to run SAF_LABEL_MAP before the regex pass and
    return on first hit, so a free-text "disease status: Healthy control"
    overwrote a deposited "saf score: S1A1F1" on the same sample. The order is
    now measurement-first: SAF_FULL_PATTERN, then the remaining FSTAGE
    patterns, then SAF_LABEL_MAP as a last resort. f_stage_source records
    which pass actually produced the value.

    Blast radius of the correction, recomputed from all 18 cached series
    matrices (288 GSM records): exactly ONE donor moves --
      GSE202379 P98 (GSM6112244, run SRR19129274): F0 -> F1, cause
      "disease status: Healthy control" short-circuiting "saf score: S1A1F1".
    P30 (GSM6112243, "saf score: S0A1F0") hits the identical short-circuit and
    does not move, because its measured grade is F0 either way.
    The 7 donors whose stage comes from a label with no SAF string at all
    (PHL1, PHL2 -> F0; PCL16, PCL17, PCL18, PCL103, PCL104 -> F4) are
    unchanged: for them the label map is now the last resort rather than the
    first, and it still supplies the only available value.

    The pre-correction release
      donor_fstage_documented.tsv sha256 9548296f...485735c6
    is left in place and readable. The corrected harvest is written beside it
    as a NEW timestamped artifact; nothing named is overwritten.

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
import argparse
import datetime as _dt
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
# The named pre-correction release. NEVER written to by this script.
FROZEN_TSV = OUT_DIR / "donor_fstage_documented.tsv"
FROZEN_SHA = "9548296ff43888a583e6d39552e87f77b37520ca96bfca798daa8f07485735c6"

# Every run writes a NEW timestamped artifact beside the frozen release.
# Override with --out to pin an exact path (the frozen path is refused).
RUN_STAMP  = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
OUT_TSV    = OUT_DIR / f"donor_fstage_documented_{RUN_STAMP}.tsv"

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

# FULL SAF pattern. The previous version of this file matched
#     r"saf\s*score[:\s]*S\d+A\d+F(\d)"
# which captured ONLY the F group and DISCARDED the S (steatosis) and A
# (activity) grades that are deposited alongside it. Those two histology
# aspects are the only NAS-like signal anywhere in the single-cell atlas, so
# throwing them away cost two axes for nothing. This pattern captures all
# three; extract_saf() below returns them.
#
# Grade ranges as deposited (GSE202379): S is 0-3 -- FOUR levels, an S0 is
# present. Any downstream code assuming a 1-3 steatosis scale is wrong.
SAF_FULL_PATTERN = re.compile(r"saf\s*score[:\s]*S(\d)A(\d)F(\d)", re.I)

# Regex patterns to extract fibrosis stage from GEO characteristics strings.
# Matches: "fibrosis_stage: F2", "Kleiner fibrosis stage: F3", "saf score: S2A3F3"
FSTAGE_PATTERNS = [
    SAF_FULL_PATTERN,
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

def extract_fstage_with_source(row: pd.Series) -> tuple[int | None, str]:
    """Return (F stage 0-4 or None, name of the pass that produced it).

    PRECEDENCE, corrected 2026-08-30. A deposited SAF grade is a MEASUREMENT;
    "healthy control" / "end stage" are free-text CATEGORY LABELS. When both
    are present on the same record the measurement wins:

      1. SAF_FULL_PATTERN     explicit S/A/F  -> "saf_regex"
      2. remaining FSTAGE_*   explicit F      -> "other_fstage_regex"
      3. SAF_LABEL_MAP        free text       -> "label_map:<label>"

    The previous order ran the label map FIRST and returned on first hit, so
    "disease status: Healthy control" short-circuited a "saf score: S1A1F1"
    deposited on the same sample. That overwrote GSE202379 donor P98's real
    F1 with F0 (donor P30 hit the identical short-circuit and agreed only
    because its SAF is F0 anyway). See the header note for the frozen-value
    consequence.

    SAF_FULL_PATTERN captures three groups (S, A, F) and the fibrosis grade is
    its LAST group; every other pattern captures F alone, so m.groups()[-1] is
    the F value under either shape.
    """
    # ---- pass 1: the deposited SAF measurement
    for col in row.index:
        if not col.startswith("char_"):
            continue
        m = SAF_FULL_PATTERN.search(str(row[col]))
        if m:
            try:
                f = int(m.groups()[-1])
            except ValueError:
                continue
            if 0 <= f <= 4:
                return f, "saf_regex"
    # ---- pass 2: any other explicitly deposited F
    for col in row.index:
        if not col.startswith("char_"):
            continue
        val = str(row[col])
        for pat in FSTAGE_PATTERNS:
            if pat is SAF_FULL_PATTERN:
                continue          # already tried, at higher precedence
            m = pat.search(val)
            if m:
                try:
                    f = int(m.groups()[-1])
                except ValueError:
                    continue
                if 0 <= f <= 4:
                    return f, "other_fstage_regex"
    # ---- pass 3: free-text label, last resort only
    for col in row.index:
        if not col.startswith("char_"):
            continue
        val = str(row[col]).lower()
        for label, fval in SAF_LABEL_MAP.items():
            if label in val:
                return fval, f"label_map:{label}"
    return None, "none"


def extract_fstage(row: pd.Series) -> int | None:
    """Value-only wrapper over extract_fstage_with_source()."""
    return extract_fstage_with_source(row)[0]


def extract_saf(row: pd.Series) -> tuple[int | None, int | None, int | None]:
    """Return the deposited (S, A, F) grades, or (None, None, None).

    This is the pass that stops S and A being discarded. It runs ONLY the SAF
    regex, so it reports what was measured.

    Since the 2026-08-30 precedence correction, extract_fstage() also puts the
    SAF measurement first, so saf_F and F_stage_documented now agree on every
    record that carries a SAF string. They can still differ in the other
    direction: a record with no SAF string has saf_F NA while
    F_stage_documented may come from SAF_LABEL_MAP or another F regex.
    f_stage_source names which pass supplied the value in every case.
    """
    for col in row.index:
        if not col.startswith("char_"):
            continue
        m = SAF_FULL_PATTERN.search(str(row[col]))
        if m:
            try:
                s_, a_, f_ = (int(g) for g in m.groups())
            except ValueError:
                continue
            if 0 <= s_ <= 4 and 0 <= a_ <= 4 and 0 <= f_ <= 4:
                return s_, a_, f_
    return None, None, None


def label_map_hit(row: pd.Series) -> str | None:
    """Which SAF_LABEL_MAP label (if any) short-circuits this row."""
    for col in row.index:
        if not col.startswith("char_"):
            continue
        val = str(row[col]).lower()
        for label in SAF_LABEL_MAP:
            if label in val:
                return label
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

def ena_map_from_harvest(path: Path) -> dict[str, str]:
    """Recover a GSM -> SRR map from a previous harvest's `source` column.

    resolve_gsm_to_srr() builds a dict, so a GSM with several ENA runs keeps
    whichever run the portal happens to list last; a refetch can therefore
    return a different run for the same GSM (observed 2026-08-30 for
    GSM6112251/2/3 against the frozen release). Pinning the map to a previous
    harvest keeps a re-run a 1-to-1 diff of that harvest instead of mixing a
    genuine correction with run-accession churn.
    """
    prev = pd.read_csv(path, sep="\t")
    out: dict[str, str] = {}
    for _, r in prev.iterrows():
        m = re.search(r"\+ENA:(GSM\d+)$", str(r["source"]))
        if m:
            out[m.group(1)] = r["sample"]
    if not out:
        raise SystemExit(f"recovered ZERO GSM->SRR pairs from {path}")
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
            fs = chars.apply(extract_fstage_with_source, axis=1)
            chars["F_stage_documented"] = [t[0] for t in fs]
            chars["f_stage_pass"] = [t[1] for t in fs]
            saf = chars.apply(extract_saf, axis=1)
            chars["saf_S"] = [t[0] for t in saf]
            chars["saf_A"] = [t[1] for t in saf]
            chars["saf_F"] = [t[2] for t in saf]
            chars["label_map_label"] = chars.apply(label_map_hit, axis=1)
            for _, r in chars.iterrows():
                gsm = r["sample"]
                fstage = r["F_stage_documented"]
                src = f"GEO_series_matrix:{local.name}"
                # pd.notna(), not "is not None": these came out of a DataFrame
                # column and an absent grade is NaN, not None.
                has_saf = pd.notna(r["saf_F"])
                has_lab = pd.notna(r["label_map_label"])
                # f_stage_source names the pass that ACTUALLY produced the
                # value (r["f_stage_pass"]), then annotates the cases where
                # precedence was decisive so they stay auditable:
                #   *_over_label_map:<label>  a label was present but demoted
                #   label_map:<label>_no_saf_string  no measurement existed
                fsrc = str(r["f_stage_pass"])
                if fsrc.startswith("label_map:") and not has_saf:
                    fsrc = f"{fsrc}_no_saf_string"
                elif fsrc in ("saf_regex", "other_fstage_regex") and has_lab:
                    fsrc = f"{fsrc}_over_label_map:{r['label_map_label']}"
                grades = {"saf_S": r["saf_S"], "saf_A": r["saf_A"], "saf_F": r["saf_F"],
                          "f_stage_source": fsrc}
                # Emit GSM-keyed row
                rows.append({"sample": gsm, "dataset": acc,
                             "F_stage_documented": fstage, "source": src, **grades})
                # Emit SRR-keyed row when ENA resolved this GSM
                srr = srr_map.get(gsm)
                if srr:
                    rows.append({"sample": srr, "dataset": acc,
                                 "F_stage_documented": fstage,
                                 "source": f"{src}+ENA:{gsm}", **grades})

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
    keep = ["sample", "dataset", "F_stage_documented", "source"]
    keep += [c for c in ("saf_S", "saf_A", "saf_F", "f_stage_source") if c in merged.columns]
    return merged[keep]

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(OUT_TSV),
                    help="output TSV path (default: a new timestamped file)")
    ap.add_argument("--ena-map-from", default=None, metavar="TSV",
                    help="pin the GSM->SRR map to a previous harvest's "
                         "`source` column instead of refetching from ENA "
                         "(see ena_map_from_harvest)")
    ap.add_argument("--geo-cache-only", action="store_true",
                    help="use only the on-disk geo_cache; never download")
    args = ap.parse_args()
    if args.ena_map_from:
        pinned = ena_map_from_harvest(Path(args.ena_map_from))
        print(f"[ENA] pinned {len(pinned)} GSM->SRR from {args.ena_map_from}")
        globals()["resolve_gsm_to_srr"] = lambda prj: pinned
        globals()["parse_geo_bioproject"] = lambda p: "PINNED"
    if args.geo_cache_only:
        def _cached_only(acc, url, _c=GEO_CACHE):
            pref = "GSE192740" if acc == "Liver_Atlas" else acc
            hits = sorted(_c.glob(f"{pref}*_series_matrix.txt.gz"))
            if not hits:
                raise SystemExit(f"--geo-cache-only: nothing cached for {acc}")
            return hits
        globals()["download_geo"] = _cached_only
    out_tsv = Path(args.out)
    if out_tsv.resolve() == FROZEN_TSV.resolve():
        raise SystemExit(
            f"REFUSING to write the named release {FROZEN_TSV}. "
            "Pass --out with a new timestamped path.")

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

    final.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n[output] {out_tsv}")
    print(f"[frozen] pre-correction release left untouched: {FROZEN_TSV}")

if __name__ == "__main__":
    main()
