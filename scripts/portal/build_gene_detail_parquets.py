#!/usr/bin/env python3
"""
build_gene_detail_parquets.py
=============================
Long-format, symbol-keyed gene-detail parquets for the MASLD Atlas v2 web app.

Replaces the per-gene-JSON approach of ``generate_gene_profiles.py``: instead of
writing 33,943 ``genes/{SYMBOL}.json`` files, this appends rows to Arrow tables
and writes 5 long-format parquets keyed + sorted by ``symbol`` (enables DuckDB
predicate pushdown / row-group skipping on ``WHERE symbol = ?``).

Builds EXACTLY to the data contract §3 (docs/superpowers/specs/atlas-web-data-contract.md):
  1. gene_per_cohort_de.parquet   (§3.1) per-study legacy + per-study-contrast DE
  2. gene_trajectories.parquet    (§3.2) fibrosis + nas staging trajectories
  3. gene_pseudobulk_de.parquet   (§3.4) per-cell-type pseudobulk DE (padj < 0.1)
  4. gene_lincs.parquet           (§3.5) LINCS/CGP reversal compounds per target
  5. gene_drugs.parquet           (§3.6) clinical drug validation + target_class

Does NOT build §3.3 gene_coloc_by_gwas.parquet — the gen-genetics agent
materializes the shared coloc_by_gwas.parquet once and the gene page queries it.

The source-path resolution + per-file parsing / shim logic is REUSED VERBATIM by
importing the loader functions from ``masld-atlas-v2/scripts/generate_gene_profiles.py``
(read-only). This script only reshapes their (ensembl/symbol -> records) maps into
long-format tables and resolves Ensembl-keyed rows to gene symbols via the atlas.

Usage:
  micromamba run -n spatial python scripts/portal/build_gene_detail_parquets.py \
      --output-dir masld-atlas-v2/public/data
"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Import the loader functions verbatim from generate_gene_profiles.py
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
_GEN_SCRIPTS_DIR = PROJECT_ROOT / "masld-atlas-v2" / "scripts"
sys.path.insert(0, str(_GEN_SCRIPTS_DIR))

import generate_gene_profiles as ggp  # noqa: E402  (path inserted above)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def build_ensembl_to_symbol(atlas: pd.DataFrame) -> dict:
    """Map version-stripped Ensembl id -> human_symbol from the atlas.

    Used to resolve Ensembl-keyed per-study / trajectory / pseudobulk records to
    the canonical gene symbol (contract §0.5). First occurrence wins.
    """
    a = atlas[["ensembl_id", "human_symbol"]].copy()
    a = a[a["human_symbol"].notna()]
    a["ens_base"] = a["ensembl_id"].astype(str).map(ggp.strip_version)
    a["sym"] = a["human_symbol"].astype(str).str.strip()
    a = a[(a["ens_base"].ne("nan")) & (a["ens_base"].ne("")) & (a["sym"].ne("nan")) & (a["sym"].ne(""))]
    a = a.drop_duplicates("ens_base", keep="first")
    return dict(zip(a["ens_base"], a["sym"]))


def write_parquet(df: pd.DataFrame, path: Path, columns, float_cols=(), int_cols=()):
    """Sort by symbol, clean ±Inf -> NaN on numeric cols, write pyarrow parquet.

    Contract §0.4: engine pyarrow, index=False, ±Inf -> NaN before write,
    symbol-keyed long-format parquets sorted by symbol.
    """
    if df.empty:
        df = pd.DataFrame(columns=list(columns))
    df = df[list(columns)].copy()
    for c in float_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").replace([np.inf, -np.inf], np.nan)
    for c in int_cols:
        vals = pd.to_numeric(df[c], errors="coerce").replace([np.inf, -np.inf], np.nan)
        df[c] = vals.astype("Int64")  # nullable integer
    # mergesort = stable, so within-symbol order (already stage/celltype-sorted) is preserved
    df = df.sort_values("symbol", kind="mergesort").reset_index(drop=True)
    df.to_parquet(path, engine="pyarrow", index=False)
    n_sym = df["symbol"].nunique() if not df.empty else 0
    print(
        f"  wrote {path.name}: {len(df):,} rows x {len(df.columns)} cols "
        f"({n_sym:,} distinct symbols) -> {path}"
    )
    print(f"    columns: {list(df.columns)}")
    return df


# ---------------------------------------------------------------------------
# Builders (one per contract sub-section)
# ---------------------------------------------------------------------------


def build_per_cohort_de(ens2sym: dict) -> pd.DataFrame:
    """§3.1 gene_per_cohort_de.parquet — per-cohort contrast DE only.

    Contrasts = {disease_vs_control, MASH_vs_control, MASL_vs_control,
    MASH_vs_MASL} via ggp.load_per_cohort_contrasts() (drops PRJNA512027).

    The retired `legacy_dream` single-contrast arm (per_study/*_de_results.csv,
    the dream-based method) is NOT shipped — project convention forbids retired-
    method labels/data in web artifacts (team-lead 2026-07-08).
    """
    rows = []
    # per-cohort per-contrast DE (per_study_contrasts/{contrast}/*_de.csv)
    per_cohort = ggp.load_per_cohort_contrasts()
    for contrast, gene_map in per_cohort.items():
        for ens, recs in gene_map.items():
            sym = ens2sym.get(ens)
            if not sym:
                continue
            for rec in recs:
                rows.append(
                    {
                        "symbol": sym,
                        "contrast": contrast,
                        "dataset": rec["dataset"],
                        "logfc": rec["logfc"],
                        "pval": rec["pval"],
                        "padj": rec.get("padj"),
                    }
                )
    return pd.DataFrame(rows)


def build_trajectories(ens2sym: dict) -> pd.DataFrame:
    """§3.2 gene_trajectories.parquet — fibrosis + nas staging trajectories.

    One row per (symbol, axis, reference, stage). n_cohorts null where the
    source is a one-vs-rest dream table rather than a per-cohort meta shim.
    """
    # (axis, reference, loader-result-map, per-record stage key)
    specs = [
        ("fibrosis", "vs_rest", ggp.load_fibrosis_stage(), "stage"),
        ("fibrosis", "vs_F0", ggp.load_stage_vs_control(), "stage"),
        ("fibrosis", "vs_healthy", ggp.load_stage_vs_healthy(), "stage"),
        ("nas", "vs_rest", ggp.load_nas_stage(), "nas"),
        ("nas", "vs_nas0", ggp.load_nas_vs_nas0(), "nas"),
        ("nas", "vs_healthy", ggp.load_nas_vs_healthy(), "nas"),
    ]
    rows = []
    for axis, reference, mapping, stage_key in specs:
        for ens, recs in mapping.items():
            sym = ens2sym.get(ens)
            if not sym:
                continue
            for rec in recs:
                rows.append(
                    {
                        "symbol": sym,
                        "axis": axis,
                        "reference": reference,
                        "stage": str(rec[stage_key]),
                        "logfc": rec.get("logfc"),
                        "padj": rec.get("padj"),
                        "n_cohorts": rec.get("n_cohorts"),
                    }
                )
    return pd.DataFrame(rows)


# Plain per-cell-type disease files carry contrast "MASLD_vs_Healthy"; normalize
# to "disease_vs_control" per the team-lead §3.4 extension (2026-07-08). Stage-
# transition files keep their transition name (Cirrhosis_vs_Steatohepatitis, …).
_PSEUDOBULK_CONTRAST_REMAP = {"MASLD_vs_Healthy": "disease_vs_control"}
_PSEUDOBULK_REQUIRED = {"gene", "cell_type", "logFC", "padj"}


def build_pseudobulk_de(ens2sym: dict) -> pd.DataFrame:
    """§3.4 gene_pseudobulk_de.parquet — per-cell-type pseudobulk DE (padj < 0.1),
    WITH a `contrast` column (team-lead §3.4 extension 2026-07-08 to disambiguate
    the plain disease files from the 3 stage-transition contrasts).

    Local extension of ggp.load_pseudobulk_de(): the verbatim loader drops the
    source `contrast` column via usecols, so we re-read it here. Inclusion criterion
    is identical to verbatim (skip any file missing gene/cell_type/logFC/padj — i.e.
    allcell / F-transition bayesprism files); padj<0.1 filter, ENSG/symbol split,
    strip_version + rounding all reuse ggp.
    """
    pbdir = ggp.PSEUDOBULK_DIR
    print(f"Loading pseudobulk DE (with contrast) from {pbdir}")
    frames = []
    for fp in sorted(pbdir.glob("*_de.csv")):
        if fp.name.startswith("archive"):
            continue
        try:
            chunk = pd.read_csv(fp)
        except Exception as e:
            print(f"  WARNING: skipping {fp.name}: {e}")
            continue
        if not _PSEUDOBULK_REQUIRED.issubset(chunk.columns):
            # mirrors the verbatim usecols skip (missing cell_type / logFC)
            continue
        if "contrast" in chunk.columns:
            contrast = chunk["contrast"].astype(str)
        else:
            # no explicit contrast column -> treat as the plain disease contrast
            contrast = pd.Series("disease_vs_control", index=chunk.index)
        sub = chunk[["gene", "cell_type", "logFC", "padj"]].copy()
        sub["contrast"] = contrast.replace(_PSEUDOBULK_CONTRAST_REMAP)
        frames.append(sub)

    cols = ["symbol", "cell_type", "contrast", "logfc", "padj"]
    if not frames:
        print("  WARNING: no pseudobulk DE files found")
        return pd.DataFrame(columns=cols)

    df = pd.concat(frames, ignore_index=True)
    df = df[df["padj"].notna() & (df["padj"] < 0.1)].copy()
    print(f"  {len(df):,} significant (padj<0.1) pseudobulk DE records")

    # ENSG* rows -> resolve via atlas; everything else is already a symbol (verbatim split)
    is_ens = df["gene"].astype(str).str.startswith("ENSG")
    df["symbol"] = np.where(
        is_ens,
        df["gene"].astype(str).map(ggp.strip_version).map(ens2sym),
        df["gene"].astype(str).str.strip(),
    )
    df = df[df["symbol"].notna() & (df["symbol"].astype(str) != "")]
    df["logfc"] = df["logFC"].map(lambda v: ggp.safe_round(v, 3))
    df["padj"] = df["padj"].map(lambda v: ggp.safe_round(v, 4))
    out = df[cols].copy()
    out["cell_type"] = out["cell_type"].astype(str)
    out["contrast"] = out["contrast"].astype(str)
    # stable pre-sort so write_parquet's symbol sort leaves tidy within-symbol order
    return out.sort_values(["symbol", "cell_type", "contrast"], kind="mergesort").reset_index(drop=True)


def build_lincs() -> pd.DataFrame:
    """§3.5 gene_lincs.parquet — LINCS/CGP compounds per target gene.

    ggp.load_lincs_ranked() returns target_symbol -> [{name, score, moa?}], with
    targets already parsed from target.x / dgidb_targets.
    """
    lincs = ggp.load_lincs_ranked()
    rows = []
    for sym, recs in lincs.items():
        sym = str(sym).strip()
        if not sym or sym.lower() == "nan":
            continue
        if not any(ch.isalnum() for ch in sym):
            # drop quoted-empty / junk target tokens (e.g. '""""') from target.x
            continue
        for rec in recs:
            rows.append(
                {
                    "symbol": sym,
                    "compound": rec.get("name"),
                    "score": rec.get("score"),
                    "moa": rec.get("moa"),
                }
            )
    return pd.DataFrame(rows)


def build_drugs() -> pd.DataFrame:
    """§3.6 gene_drugs.parquet — clinical drug validation + progression target_class.

    ggp.load_clinical_drugs() gives target_symbol -> [{drug, stage, moa, support}];
    ggp.load_drug_progression() gives target_symbol -> {..., target_class} (one row
    per target). target_class joined by symbol (nullable).
    """
    clinical = ggp.load_clinical_drugs()
    prog = ggp.load_drug_progression()
    rows = []
    for sym, recs in clinical.items():
        sym = str(sym).strip()
        if not sym or sym.lower() == "nan":
            continue
        target_class = prog.get(sym, {}).get("target_class")
        for rec in recs:
            rows.append(
                {
                    "symbol": sym,
                    "drug": rec.get("drug"),
                    "stage": rec.get("stage"),
                    "moa": rec.get("moa"),
                    "support": rec.get("support"),
                    "target_class": target_class,
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Build long-format, symbol-keyed gene-detail parquets (contract §3)."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(PROJECT_ROOT / "masld-atlas-v2" / "public" / "data"),
        help="Output directory for the parquet files.",
    )
    parser.add_argument(
        "--only",
        type=str,
        default=None,
        help="Comma-separated subset to (re)build: "
        "per_cohort_de,trajectories,pseudobulk_de,lincs,drugs. Default: all 5.",
    )
    args = parser.parse_args()

    all_keys = ["per_cohort_de", "trajectories", "pseudobulk_de", "lincs", "drugs"]
    if args.only:
        only = {s.strip() for s in args.only.split(",") if s.strip()}
        bad = only - set(all_keys)
        if bad:
            parser.error(f"--only has unknown keys: {sorted(bad)}; valid: {all_keys}")
    else:
        only = set(all_keys)

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {out_dir}")
    print(f"Building: {[k for k in all_keys if k in only]}\n")

    t0 = time.time()

    # Atlas -> ensembl_base -> symbol map (only needed for symbol-resolving builds).
    ens2sym = {}
    if only & {"per_cohort_de", "trajectories", "pseudobulk_de"}:
        atlas = ggp.load_atlas()  # asserts bulk_* present
        ens2sym = build_ensembl_to_symbol(atlas)
        print(f"  ensembl->symbol map: {len(ens2sym):,} entries\n")

    if "per_cohort_de" in only:
        print("[per_cohort_de] gene_per_cohort_de.parquet (§3.1)")
        df1 = build_per_cohort_de(ens2sym)
        write_parquet(
            df1,
            out_dir / "gene_per_cohort_de.parquet",
            columns=["symbol", "contrast", "dataset", "logfc", "pval", "padj"],
            float_cols=["logfc", "pval", "padj"],
        )
        print(f"    contrasts: {sorted(df1['contrast'].unique()) if not df1.empty else []}\n")

    if "trajectories" in only:
        print("[trajectories] gene_trajectories.parquet (§3.2)")
        df2 = build_trajectories(ens2sym)
        write_parquet(
            df2,
            out_dir / "gene_trajectories.parquet",
            columns=["symbol", "axis", "reference", "stage", "logfc", "padj", "n_cohorts"],
            float_cols=["logfc", "padj"],
            int_cols=["n_cohorts"],
        )
        if not df2.empty:
            combos = df2[["axis", "reference"]].drop_duplicates().sort_values(["axis", "reference"])
            print("    axis x reference:")
            for _, r in combos.iterrows():
                print(f"      {r['axis']:9s} {r['reference']}")
        print()

    if "pseudobulk_de" in only:
        print("[pseudobulk_de] gene_pseudobulk_de.parquet (§3.4)")
        df3 = build_pseudobulk_de(ens2sym)
        write_parquet(
            df3,
            out_dir / "gene_pseudobulk_de.parquet",
            columns=["symbol", "cell_type", "contrast", "logfc", "padj"],
            float_cols=["logfc", "padj"],
        )
        if not df3.empty:
            print(f"    cell types: {sorted(df3['cell_type'].unique())}")
            print(f"    contrasts:  {sorted(df3['contrast'].unique())}")
        print()

    if "lincs" in only:
        print("[lincs] gene_lincs.parquet (§3.5)")
        df4 = build_lincs()
        write_parquet(
            df4,
            out_dir / "gene_lincs.parquet",
            columns=["symbol", "compound", "score", "moa"],
            float_cols=["score"],
        )
        print()

    if "drugs" in only:
        print("[drugs] gene_drugs.parquet (§3.6)")
        df5 = build_drugs()
        write_parquet(
            df5,
            out_dir / "gene_drugs.parquet",
            columns=["symbol", "drug", "stage", "moa", "support", "target_class"],
        )
        print()

    print(f"{'=' * 64}")
    n = len(only)
    print(f"Done in {time.time() - t0:.1f}s. {n} gene-detail parquet(s) in {out_dir}")


if __name__ == "__main__":
    main()
