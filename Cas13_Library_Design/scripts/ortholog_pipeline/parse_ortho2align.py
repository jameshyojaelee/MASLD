#!/usr/bin/env python
"""
parse_ortho2align.py — Parse ortho2align output into L9_ortho2align.tsv

Reads the ortho2align run_pipeline outputs:
  - bestSignificant.annotation.tsv      (query->mouse gene annotation with Jaccard/OC)
  - bestSignificant.subject_orthologs.tsv (mouse coordinates + q-values)
  - bestSignificant.query_orthologs.tsv   (human coordinates + alignment score)
  - stats.txt                            (pipeline statistics)

Joins with GENCODE gene lookups for symbol/biotype, then compares to Phase J
(L1.5_phasej_synteny.tsv) for recovery/disambiguation statistics.

Output: data/external/orthologs/layers/L9_ortho2align.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np

# ============================================================================
# Paths
# ============================================================================
PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
WORK = PROJECT / "Cas13_Library_Design/scripts/ortholog_pipeline/ortho2align_work"
LAYERS = PROJECT / "data/external/orthologs/layers"

DEFAULT_RUNDIR = WORK / "full_run"
MOUSE_LOOKUP = WORK / "mouse_gene_lookup.tsv"
HUMAN_LOOKUP = WORK / "human_lncrna_lookup.tsv"
PHASEJ_LAYER = LAYERS / "L1.5_phasej_synteny.tsv"
OUT_LAYER = LAYERS / "L9_ortho2align.tsv"


def load_gene_lookup(path: Path) -> pd.DataFrame:
    """Load ensembl_id -> gene_name, gene_biotype lookup."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["ensembl_id"] = df["ensembl_id"].str.split(".").str[0]
    return df.drop_duplicates(subset=["ensembl_id"]).set_index("ensembl_id")


def parse_annotation(rundir: Path) -> pd.DataFrame:
    """Parse bestSignificant.annotation.tsv.

    Columns: Query, Orthologs, Query_length, Orthologs_lengths, JI, OC
    Orthologs can be comma-separated (multiple mouse genes overlapping the
    ortholog region) or 'NotAnnotated'.
    """
    ann_path = rundir / "annotation_files" / "bestSignificant.annotation.tsv"
    if not ann_path.exists():
        raise FileNotFoundError(f"Annotation file not found: {ann_path}")

    df = pd.read_csv(ann_path, sep="\t", dtype=str)
    print(f"  Loaded {len(df)} best-significant orthologs from annotation file")
    return df


def parse_subject_bed(rundir: Path) -> pd.DataFrame:
    """Parse bestSignificant.subject_orthologs.tsv (BED12 + q-values).

    Col 4 = query gene name, Col 5 = alignment score (sum of HSP raw scores),
    Col 13 = comma-separated q-values per block.
    """
    tsv_path = rundir / "bestSignificant.subject_orthologs.tsv"
    if not tsv_path.exists():
        raise FileNotFoundError(f"Subject orthologs TSV not found: {tsv_path}")

    cols = ["chrom", "start", "end", "query_gene", "score", "strand",
            "thick_start", "thick_end", "item_rgb", "block_count",
            "block_sizes", "block_starts", "qvalues"]
    df = pd.read_csv(tsv_path, sep="\t", header=None, names=cols, dtype=str)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")

    # Parse q-values: take the minimum q-value across blocks
    def min_qval(qv_str):
        try:
            vals = [float(x) for x in qv_str.split(",") if x.strip()]
            return min(vals) if vals else np.nan
        except (ValueError, AttributeError):
            return np.nan

    df["min_qvalue"] = df["qvalues"].apply(min_qval)

    # Compute total aligned length from block_sizes
    def total_aligned(bs_str):
        try:
            return sum(int(x) for x in bs_str.split(",") if x.strip())
        except (ValueError, AttributeError):
            return 0

    df["aligned_length"] = df["block_sizes"].apply(total_aligned)

    # Subject coordinates
    df["subject_chrom"] = df["chrom"]
    df["subject_start"] = pd.to_numeric(df["start"], errors="coerce")
    df["subject_end"] = pd.to_numeric(df["end"], errors="coerce")

    return df[["query_gene", "score", "min_qvalue", "aligned_length",
               "subject_chrom", "subject_start", "subject_end"]].copy()


def parse_query_bed(rundir: Path) -> pd.DataFrame:
    """Parse bestSignificant.query_orthologs.tsv — for alignment coverage info."""
    tsv_path = rundir / "bestSignificant.query_orthologs.tsv"
    if not tsv_path.exists():
        return pd.DataFrame()

    cols = ["chrom", "start", "end", "query_gene", "score", "strand",
            "thick_start", "thick_end", "item_rgb", "block_count",
            "block_sizes", "block_starts", "qvalues"]
    df = pd.read_csv(tsv_path, sep="\t", header=None, names=cols, dtype=str)

    def total_aligned(bs_str):
        try:
            return sum(int(x) for x in bs_str.split(",") if x.strip())
        except (ValueError, AttributeError):
            return 0

    df["query_aligned_length"] = df["block_sizes"].apply(total_aligned)
    df["query_total_length"] = pd.to_numeric(df["end"], errors="coerce") - \
                                pd.to_numeric(df["start"], errors="coerce")
    df["alignment_coverage"] = df["query_aligned_length"] / df["query_total_length"]
    df["alignment_coverage"] = df["alignment_coverage"].clip(0, 1)

    return df[["query_gene", "query_aligned_length",
               "query_total_length", "alignment_coverage"]].copy()


def build_layer(rundir: Path) -> pd.DataFrame:
    """Build L9 ortho2align layer table."""
    # Load gene lookups
    mouse_lk = load_gene_lookup(MOUSE_LOOKUP)
    human_lk = load_gene_lookup(HUMAN_LOOKUP)

    # Parse ortho2align outputs
    ann_df = parse_annotation(rundir)
    subj_df = parse_subject_bed(rundir)
    query_df = parse_query_bed(rundir)

    # Merge annotation with subject BED scores
    merged = ann_df.merge(subj_df, left_on="Query", right_on="query_gene",
                          how="left")

    # Merge with query BED coverage
    if not query_df.empty:
        merged = merged.merge(query_df, left_on="Query", right_on="query_gene",
                              how="left", suffixes=("", "_q"))

    # Expand multi-gene annotations: pick the BEST mouse gene per query
    # Best = highest Jaccard index (or highest overlap coefficient if tied)
    rows = []
    for _, row in merged.iterrows():
        human_id = str(row["Query"]).strip()
        orthologs_str = str(row.get("Orthologs", "NotAnnotated")).strip()

        # Scores from subject BED
        score = float(row.get("score", 0)) if pd.notna(row.get("score")) else 0.0
        min_qval = float(row.get("min_qvalue", np.nan)) if pd.notna(row.get("min_qvalue")) else np.nan
        aligned_len = int(row.get("aligned_length", 0)) if pd.notna(row.get("aligned_length")) else 0
        coverage = float(row.get("alignment_coverage", np.nan)) if pd.notna(row.get("alignment_coverage")) else np.nan

        if orthologs_str == "NotAnnotated" or not orthologs_str:
            # Ortholog found but doesn't overlap any known mouse gene
            rows.append({
                "human_ensembl": human_id,
                "mouse_ensembl": pd.NA,
                "ortho2align_score": score,
                "ortho2align_min_qvalue": min_qval,
                "ortho2align_aligned_bp": aligned_len,
                "ortho2align_coverage": coverage,
                "ortho2align_jaccard": 0.0,
                "ortho2align_overlap_coeff": 0.0,
                "ortho2align_annotated": False,
            })
            continue

        # Multiple overlapping genes: split and pick best by Jaccard
        gene_ids = [g.strip() for g in orthologs_str.split(",")]
        ji_str = str(row.get("JI", "0"))
        oc_str = str(row.get("OC", "0"))
        jis = [float(x) for x in ji_str.split(",") if x.strip()]
        ocs = [float(x) for x in oc_str.split(",") if x.strip()]

        # Pad if needed
        while len(jis) < len(gene_ids):
            jis.append(0.0)
        while len(ocs) < len(gene_ids):
            ocs.append(0.0)

        # Pick the best by Jaccard, break ties by OC
        best_idx = 0
        best_ji = jis[0]
        best_oc = ocs[0]
        for i in range(1, len(gene_ids)):
            if jis[i] > best_ji or (jis[i] == best_ji and ocs[i] > best_oc):
                best_idx = i
                best_ji = jis[i]
                best_oc = ocs[i]

        rows.append({
            "human_ensembl": human_id,
            "mouse_ensembl": gene_ids[best_idx],
            "ortho2align_score": score,
            "ortho2align_min_qvalue": min_qval,
            "ortho2align_aligned_bp": aligned_len,
            "ortho2align_coverage": coverage,
            "ortho2align_jaccard": best_ji,
            "ortho2align_overlap_coeff": best_oc,
            "ortho2align_annotated": True,
            "n_overlapping_genes": len(gene_ids),
        })

    result = pd.DataFrame(rows)
    print(f"  Built {len(result)} ortholog rows")
    print(f"    Annotated (overlaps known mouse gene): "
          f"{result['ortho2align_annotated'].sum()}")
    print(f"    Not annotated (novel syntenic region): "
          f"{(~result['ortho2align_annotated']).sum()}")

    # Strip Ensembl versions
    result["human_ensembl"] = result["human_ensembl"].astype(str).str.split(".").str[0]
    if "mouse_ensembl" in result.columns:
        result["mouse_ensembl"] = result["mouse_ensembl"].astype(str).str.split(".").str[0]
        result["mouse_ensembl"] = result["mouse_ensembl"].replace({"nan": pd.NA})

    # Join human symbol/biotype
    result = result.merge(
        human_lk[["gene_name", "gene_biotype"]].rename(
            columns={"gene_name": "human_symbol", "gene_biotype": "human_biotype"}),
        left_on="human_ensembl", right_index=True, how="left")

    # Join mouse symbol/biotype
    result = result.merge(
        mouse_lk[["gene_name", "gene_biotype"]].rename(
            columns={"gene_name": "mouse_symbol", "gene_biotype": "mouse_biotype"}),
        left_on="mouse_ensembl", right_index=True, how="left")

    # Add tier/provenance columns
    result["tier_M_ortho2align"] = 1
    result["confidence_tier"] = "M"
    result["provenance_sources"] = "ortho2align_v1.0.5"

    # Reorder columns
    col_order = [
        "mouse_ensembl", "mouse_symbol", "mouse_biotype",
        "human_ensembl", "human_symbol", "human_biotype",
        "tier_M_ortho2align",
        "ortho2align_score", "ortho2align_min_qvalue",
        "ortho2align_aligned_bp", "ortho2align_coverage",
        "ortho2align_jaccard", "ortho2align_overlap_coeff",
        "ortho2align_annotated", "n_overlapping_genes",
        "confidence_tier", "provenance_sources",
    ]
    for c in col_order:
        if c not in result.columns:
            result[c] = pd.NA
    result = result[col_order]

    return result


def compare_to_phasej(l9: pd.DataFrame) -> None:
    """Compare ortho2align results to Phase J synteny layer."""
    if not PHASEJ_LAYER.exists():
        print("\n  Phase J layer not found; skipping comparison")
        return

    pj = pd.read_csv(PHASEJ_LAYER, sep="\t", dtype=str)
    pj_human = set(pj["human_ensembl"].dropna().str.split(".").str[0].unique())
    l9_human = set(l9["human_ensembl"].dropna().unique())

    overlap = pj_human & l9_human
    pj_only = pj_human - l9_human
    l9_only = l9_human - pj_human

    print(f"\n  === Comparison to Phase J (L1.5_phasej_synteny) ===")
    print(f"  Phase J human anchors:       {len(pj_human):,}")
    print(f"  ortho2align human hits:      {len(l9_human):,}")
    print(f"  Recovery (PJ anchors in L9): {len(overlap):,} / {len(pj_human):,} "
          f"({100*len(overlap)/len(pj_human):.1f}%)")
    print(f"  Phase J-only (no L9 hit):    {len(pj_only):,}")
    print(f"  L9-only (new discoveries):   {len(l9_only):,}")

    # For annotated L9 hits that overlap Phase J: does ortho2align pick the same
    # mouse gene as Phase J listed first?
    l9_annotated = l9.dropna(subset=["mouse_ensembl"]).copy()
    l9_annotated["mouse_ensembl_clean"] = l9_annotated["mouse_ensembl"].astype(str).str.split(".").str[0]

    pj_pairs = pj[["human_ensembl", "mouse_ensembl"]].copy()
    pj_pairs["human_ensembl"] = pj_pairs["human_ensembl"].astype(str).str.split(".").str[0]
    pj_pairs["mouse_ensembl"] = pj_pairs["mouse_ensembl"].astype(str).str.split(".").str[0]

    # Phase J first mouse gene per human anchor
    pj_first = pj_pairs.drop_duplicates(subset=["human_ensembl"], keep="first")
    pj_first = pj_first.set_index("human_ensembl")["mouse_ensembl"]

    # All Phase J mouse genes per human anchor
    pj_all = pj_pairs.groupby("human_ensembl")["mouse_ensembl"].apply(set).to_dict()

    # Count partners per anchor in Phase J
    pj_counts = pj_pairs.groupby("human_ensembl")["mouse_ensembl"].nunique()

    match_first = 0
    match_any = 0
    different = 0
    total = 0

    for _, row in l9_annotated.iterrows():
        h = row["human_ensembl"]
        m = row["mouse_ensembl_clean"]
        if h not in pj_all:
            continue
        total += 1
        if h in pj_first.index and pj_first[h] == m:
            match_first += 1
        if m in pj_all.get(h, set()):
            match_any += 1
        else:
            different += 1

    print(f"\n  --- Annotated ortho2align hits vs Phase J partners ---")
    print(f"  Total compared:              {total}")
    print(f"  Same as PJ first-listed:     {match_first} ({100*match_first/max(total,1):.1f}%)")
    print(f"  Among ANY PJ partner:        {match_any} ({100*match_any/max(total,1):.1f}%)")
    print(f"  New (not in PJ set):         {different} ({100*different/max(total,1):.1f}%)")

    # Disambiguation: Phase J anchors with 10+ partners that got a single best
    big_anchors = pj_counts[pj_counts >= 10].index
    big_in_l9 = l9_annotated[l9_annotated["human_ensembl"].isin(big_anchors)]
    print(f"\n  --- Disambiguation of multi-partner Phase J anchors ---")
    print(f"  PJ anchors with >=10 partners: {len(big_anchors):,}")
    print(f"  Of those with L9 annotated hit: {big_in_l9['human_ensembl'].nunique():,}")


def spot_check_canonical(l9: pd.DataFrame) -> None:
    """Spot-check canonical lncRNA pairs."""
    canonical = {
        "MALAT1": "Malat1",
        "NEAT1": "Neat1",
        "MEG3": "Meg3",
        "H19": "H19",
        "XIST": "Xist",
        "HOTAIR": "Hotair",
    }
    print(f"\n  === Canonical lncRNA spot-check ===")
    for human_name, mouse_name in canonical.items():
        hit = l9[l9["human_symbol"] == human_name]
        if len(hit) == 0:
            print(f"  {human_name}: NOT FOUND in ortho2align output")
        else:
            row = hit.iloc[0]
            mouse_match = row.get("mouse_symbol", "")
            score = row.get("ortho2align_score", 0)
            qval = row.get("ortho2align_min_qvalue", "NA")
            ji = row.get("ortho2align_jaccard", 0)
            ann = row.get("ortho2align_annotated", False)
            expected = "MATCH" if str(mouse_match) == mouse_name else f"MISMATCH (got {mouse_match})"
            print(f"  {human_name} -> {mouse_match} [{expected}] | "
                  f"score={score:.1f}, q={qval}, JI={ji:.4f}, annotated={ann}")


def main():
    parser = argparse.ArgumentParser(
        description="Parse ortho2align output into L9_ortho2align.tsv")
    parser.add_argument("--rundir", type=str, default=str(DEFAULT_RUNDIR),
                        help="ortho2align output directory")
    parser.add_argument("--output", type=str, default=str(OUT_LAYER),
                        help="Output L9 layer TSV path")
    args = parser.parse_args()

    rundir = Path(args.rundir)
    print(f"Parsing ortho2align output from: {rundir}")

    l9 = build_layer(rundir)

    # Save
    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    l9.to_csv(outpath, sep="\t", index=False)
    print(f"\n  Wrote L9 layer: {outpath} ({len(l9):,} rows)")

    # Comparison and spot-checks
    compare_to_phasej(l9)
    spot_check_canonical(l9)

    # Summary statistics
    print(f"\n  === Summary ===")
    print(f"  Total orthologs:          {len(l9):,}")
    ann_count = l9["ortho2align_annotated"].sum()
    print(f"  Annotated (known gene):   {ann_count:,}")
    print(f"  Unannotated (novel):      {len(l9) - ann_count:,}")
    if "ortho2align_score" in l9.columns:
        scores = l9["ortho2align_score"].dropna()
        print(f"  Score: median={scores.median():.1f}, "
              f"mean={scores.mean():.1f}, max={scores.max():.1f}")
    if "ortho2align_coverage" in l9.columns:
        cov = l9["ortho2align_coverage"].dropna()
        print(f"  Coverage: median={cov.median():.3f}, "
              f"mean={cov.mean():.3f}")

    # Mouse biotype breakdown
    if "mouse_biotype" in l9.columns:
        bt = l9.dropna(subset=["mouse_biotype"])["mouse_biotype"].value_counts()
        print(f"\n  --- Mouse gene biotype breakdown ---")
        for b, n in bt.head(10).items():
            print(f"    {b}: {n:,}")


if __name__ == "__main__":
    main()
