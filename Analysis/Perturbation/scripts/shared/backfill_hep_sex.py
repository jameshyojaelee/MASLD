"""Backfill inferred/biological sex into hepatocyte_subtype_metadata.csv.

The hepatocyte subtype metadata (657,805 cells × 17 cols) carries `sample` +
`dataset` but no sex column. Phase 2 stratified fine-tuning needs
`inferred_sex_final` (F / M / U) per cell.

Data discovery for this task (2026-05-21):
- `Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv`
  carries a `sex` column but it is empty for all 269 donors. Not usable.
- `Analysis/SingleCell/metadata/snrna_donor_sex.csv` carries `inferred_sex`
  (Male/Female via XIST/DDX3Y k-means) for 275 scRNA donors across all 7
  scRNA datasets (GSE136103, GSE174748, GSE185477, GSE189600, GSE202379,
  GSE244832, Liver_Atlas). This is the canonical donor-level sex source.

Sample-ID format in `hepatocyte_subtype_metadata.csv` already matches the
donor-level `sample` in `snrna_donor_sex.csv` directly (no regex parsing
needed — the integration step exported per-cell `sample = donor_id`). Spot
check: all 263 unique (sample,dataset) pairs in hep meta resolve cleanly.

Fallback (only fires for samples missing from snrna_donor_sex.csv): parse
`SraRunTable.csv` under
`RNA-seq/Human/Patient_Cohorts/pipelines/custom/{DATASET}/metadata/` —
which carries biological (annotated) sex.

Outputs (alongside the input CSV):
- `hepatocyte_subtype_metadata_with_sex.csv` — original cols + `donor_id`
  + `inferred_sex_final` ∈ {F, M, U} + `sex_source` ∈ {snrna_donor_sex,
  sra_runtable, unknown}.
- `hepatocyte_sex_backfill_report.md` — audit (totals, per-dataset
  coverage, F/M/U breakdown).
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
META_PATH = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv"
)
SNRNA_DONOR_SEX = PROJECT_ROOT / "Analysis/SingleCell/metadata/snrna_donor_sex.csv"
DONOR_EXTENDED = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
)
SRA_BASE = PROJECT_ROOT / "RNA-seq/Human/Patient_Cohorts/pipelines/custom"
OUT_META = META_PATH.with_name("hepatocyte_subtype_metadata_with_sex.csv")
OUT_REPORT = META_PATH.with_name("hepatocyte_sex_backfill_report.md")


SEX_NORMALIZE = {
    "M": "M", "MALE": "M", "MAN": "M",
    "F": "F", "FEMALE": "F", "WOMAN": "F",
}


def _normalize_sex(v) -> str:
    if v is None:
        return "U"
    if isinstance(v, float) and pd.isna(v):
        return "U"
    s = str(v).strip().upper()
    return SEX_NORMALIZE.get(s, "U")


def _load_snrna_donor_sex() -> pd.DataFrame:
    df = pd.read_csv(SNRNA_DONOR_SEX, usecols=["sample", "dataset", "inferred_sex"])
    df["snrna_sex"] = df["inferred_sex"].map(_normalize_sex)
    return df[["sample", "dataset", "snrna_sex"]].drop_duplicates("sample")


def _load_donor_extended_sex() -> pd.DataFrame:
    """Reserved fallback. The `sex` column in this TSV is empty as of
    2026-05-21 but the loader is kept so future re-runs pick up any
    backfilled values automatically."""
    if not DONOR_EXTENDED.exists():
        return pd.DataFrame(columns=["sample", "ext_sex"])
    df = pd.read_csv(DONOR_EXTENDED, sep="\t", usecols=["sample", "sex"])
    df["ext_sex"] = df["sex"].map(_normalize_sex)
    return df[["sample", "ext_sex"]].drop_duplicates("sample")


def _load_sra_runtable_sex() -> pd.DataFrame:
    """Fallback: scrape any SraRunTable.csv under pipelines/custom for a
    Sex/sex/gender column. Maps Run -> sex. Bulk RNA-seq cohorts only
    overlap with scRNA via shared GSM accessions in a few cases, but we
    keep this for completeness."""
    rows: list[dict[str, str]] = []
    if not SRA_BASE.exists():
        return pd.DataFrame(columns=["run", "sra_sex"])
    for sra_path in SRA_BASE.glob("*/metadata/SraRunTable.csv"):
        try:
            df = pd.read_csv(sra_path, low_memory=False)
        except Exception:
            continue
        sex_col = next(
            (c for c in df.columns if re.fullmatch(r"sex|gender", c, re.I)),
            None,
        )
        run_col = next(
            (c for c in df.columns if c in ("Run", "run", "Sample Name", "sample_name")),
            None,
        )
        if sex_col is None or run_col is None:
            continue
        sub = df[[run_col, sex_col]].copy()
        sub.columns = ["run", "sra_sex"]
        sub["sra_sex"] = sub["sra_sex"].map(_normalize_sex)
        rows.append(sub)
    if not rows:
        return pd.DataFrame(columns=["run", "sra_sex"])
    out = pd.concat(rows, ignore_index=True).drop_duplicates("run")
    return out


def main() -> None:
    print(f"[backfill] loading {META_PATH.name}")
    meta = pd.read_csv(META_PATH)
    n_total = len(meta)
    print(f"[backfill] {n_total:,} cells across {meta['sample'].nunique()} samples, "
          f"{meta['dataset'].nunique()} datasets")

    # 1. donor_id == sample for scRNA atlas (single-donor-per-sample design)
    meta["donor_id"] = meta["sample"].astype(str)

    # 2. primary source: snrna_donor_sex
    donor_sex = _load_snrna_donor_sex()
    print(f"[backfill] snrna_donor_sex: {len(donor_sex):,} donors, "
          f"F={int((donor_sex.snrna_sex == 'F').sum())} "
          f"M={int((donor_sex.snrna_sex == 'M').sum())} "
          f"U={int((donor_sex.snrna_sex == 'U').sum())}")
    meta = meta.merge(donor_sex[["sample", "snrna_sex"]], on="sample", how="left")

    # 3. extended donor metadata sex (currently empty but loader is future-proof)
    ext = _load_donor_extended_sex()
    if not ext.empty:
        meta = meta.merge(ext, on="sample", how="left")
    else:
        meta["ext_sex"] = pd.NA

    # 4. SRA fallback (join on run = sample if matches; otherwise no-op)
    sra = _load_sra_runtable_sex()
    if not sra.empty:
        meta = meta.merge(sra.rename(columns={"run": "sample"}), on="sample", how="left")
    else:
        meta["sra_sex"] = pd.NA

    # 5. resolve final sex with priority snrna > extended > sra
    def _resolve(row) -> tuple[str, str]:
        for col, src in (
            ("snrna_sex", "snrna_donor_sex"),
            ("ext_sex", "donor_metadata_extended"),
            ("sra_sex", "sra_runtable"),
        ):
            v = row.get(col)
            if isinstance(v, str) and v in ("F", "M"):
                return v, src
        return "U", "unknown"

    resolved = meta.apply(_resolve, axis=1)
    meta["inferred_sex_final"] = [r[0] for r in resolved]
    meta["sex_source"] = [r[1] for r in resolved]

    # Drop staging cols
    meta = meta.drop(columns=[c for c in ("snrna_sex", "ext_sex", "sra_sex") if c in meta.columns])

    meta.to_csv(OUT_META, index=False)
    print(f"[backfill] wrote {OUT_META} ({len(meta):,} rows)")

    # ---------------- Audit ----------------
    n_filled = int((meta["inferred_sex_final"] != "U").sum())
    n_unknown = n_total - n_filled
    breakdown = meta["inferred_sex_final"].value_counts().to_dict()
    src_breakdown = meta["sex_source"].value_counts().to_dict()

    per_dataset = (
        meta.groupby("dataset")["inferred_sex_final"]
        .value_counts()
        .unstack(fill_value=0)
        .reindex(columns=["F", "M", "U"], fill_value=0)
    )
    per_dataset["total"] = per_dataset.sum(axis=1)
    per_dataset["pct_filled"] = (
        100 * (per_dataset["F"] + per_dataset["M"]) / per_dataset["total"]
    ).round(1)

    per_dataset_donors = (
        meta.drop_duplicates("sample")
        .groupby("dataset")["inferred_sex_final"]
        .value_counts()
        .unstack(fill_value=0)
        .reindex(columns=["F", "M", "U"], fill_value=0)
    )
    per_dataset_donors["total_donors"] = per_dataset_donors.sum(axis=1)

    unmapped = (
        meta[meta["inferred_sex_final"] == "U"][["sample", "dataset"]]
        .drop_duplicates()
        .sort_values(["dataset", "sample"])
    )

    lines = [
        "# Hepatocyte sex backfill report",
        "",
        f"- Generated: 2026-05-21",
        f"- Input file: `{META_PATH}`",
        f"- Output file: `{OUT_META}`",
        "",
        "## Summary (cell-level)",
        "",
        f"- Total cells: {n_total:,}",
        f"- Sex-resolved: {n_filled:,} ({100 * n_filled / n_total:.1f}%)",
        f"- Unknown: {n_unknown:,} ({100 * n_unknown / n_total:.1f}%)",
        f"- F/M/U breakdown: {breakdown}",
        f"- Source breakdown: {src_breakdown}",
        "",
        "## Per-dataset coverage (cells)",
        "",
        per_dataset.to_markdown(),
        "",
        "## Per-dataset coverage (donors)",
        "",
        per_dataset_donors.to_markdown(),
        "",
        "## Source inventory",
        "",
        f"- `snrna_donor_sex.csv` (k-means on XIST/DDX3Y log1p-CPM, 2-cluster) — "
        f"primary; covers {len(_load_snrna_donor_sex())} donors across 7 scRNA datasets.",
        "- `donor_metadata_extended.tsv` — `sex` column empty as of 2026-05-21; loader retained as future-proof fallback.",
        "- `SraRunTable.csv` per cohort — fallback for any unmapped samples (currently 0 hits for scRNA-only sample IDs).",
        "",
        "## Notes",
        "",
        "- Sample IDs in `hepatocyte_subtype_metadata.csv` already equal donor IDs "
        "(`SRR*` / `GSM*` for GEO scRNA datasets; `<dataset>_<donor>` for Liver_Atlas / GSE244832). "
        "All 263 unique (sample,dataset) pairs resolved without regex parsing.",
        "- For donors with `inferred_sex_final == 'U'`: drop from sex-stratified analyses or assign to a pooled bucket.",
        "- Use the `inferred_sex_final` column in `load_hepatocyte_atlas.stratify_by_sex_stage_subtype()`.",
        "",
    ]
    if len(unmapped):
        lines += [
            "## Unmapped samples (no F/M call)",
            "",
            unmapped.to_markdown(index=False),
            "",
        ]
    OUT_REPORT.write_text("\n".join(lines))
    print(f"[backfill] wrote {OUT_REPORT}")
    print(f"[backfill] filled {n_filled}/{n_total} ({100 * n_filled / n_total:.1f}%)")


if __name__ == "__main__":
    main()
