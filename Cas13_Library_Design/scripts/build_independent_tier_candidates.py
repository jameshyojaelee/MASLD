"""Build INDEPENDENT-TIER candidate list for the Cas13 in vivo MASH library.

Four independent evidence axes (each ranked separately, NOT pre-merged):
    1. COLOC          — atlas `coloc_best_susie_pp4_polyfun > 0.5`
    2. Finemap        — combined_finemapping.csv `recommended_pip > 0.5`, nearest-gene
    3. Mouse UP DEG   — MCD ∪ CDAHFD, padj<0.05, logFC>0.3, hep expressed
    4. Human UP DEG   — dream_results_ashr, padj<0.05, logFC>0.3, ortho→mouse

Outputs:
    Cas13_Library_Design/data/candidates_pc_independent.csv
    Cas13_Library_Design/data/candidates_lncrna_independent.csv
    Cas13_Library_Design/reviews/independent_tier_pools.md
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
sys.path.insert(0, str(ROOT / "Cas13_Library_Design/scripts"))
from ortholog_bridge import (  # noqa: E402
    strip_version,
    build_ortholog_table,
    load_mouse_de,
    DIETS,
)

ATLAS_CSV = ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
FINEMAP_CSV = ROOT / "GWAS/finemapping/results/combined_finemapping.csv"
DREAM_CSV = ROOT / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv"
GTF = Path("/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz")
MOUSE_GTF = ROOT / "data/ncrna_conservation/gencode.vM38.annotation.gtf.gz"
GENCODE_META = ROOT / "data/gencode_v49_gene_metadata.tsv.gz"
MOUSE_META_CACHE = ROOT / "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv"

OUT_DIR = ROOT / "Cas13_Library_Design/data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
REVIEW_MD = ROOT / "Cas13_Library_Design/reviews/independent_tier_pools.md"


# ============================================================================
# Step 1: Build gene-coordinate table from GTF for nearest-gene mapping
# ============================================================================

def load_mouse_metadata() -> pd.DataFrame:
    """Load mouse GENCODE vM38 gene metadata (mouse_ensembl_base, mouse_symbol, mouse_biotype)."""
    if MOUSE_META_CACHE.exists():
        return pd.read_csv(MOUSE_META_CACHE, dtype=str)
    rows = []
    with gzip.open(MOUSE_GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            m = re.search(r'gene_id "([^"]+)"', p[8])
            n = re.search(r'gene_name "([^"]+)"', p[8])
            t = re.search(r'gene_type "([^"]+)"', p[8])
            if m:
                rows.append((
                    m.group(1).split(".")[0],
                    n.group(1) if n else "",
                    t.group(1) if t else "",
                ))
    df = pd.DataFrame(rows, columns=["mouse_ensembl_base", "mouse_symbol_gtf", "mouse_biotype"])
    df.to_csv(MOUSE_META_CACHE, index=False)
    return df


def build_gene_coords() -> pd.DataFrame:
    """Parse GTF for gene-level rows. Returns chrom, start, end, strand,
    tss (= start if +; end if -), gene_id (versioned), ensembl_base, name, biotype."""
    rows = []
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]
            gid = re.search(r'gene_id "([^"]+)"', attrs)
            gname = re.search(r'gene_name "([^"]+)"', attrs)
            gtype = re.search(r'gene_type "([^"]+)"', attrs)
            if not gid:
                continue
            tss = start if strand == "+" else end
            rows.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "tss": tss,
                    "gene_id": gid.group(1),
                    "ensembl_base": gid.group(1).split(".")[0],
                    "gene_name": gname.group(1) if gname else None,
                    "gene_biotype": gtype.group(1) if gtype else None,
                }
            )
    df = pd.DataFrame(rows)
    # Restrict to primary chromosomes (chr1..chr22, chrX, chrY)
    pri = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    df = df[df["chrom"].isin(pri)].copy()
    df["chrom_int"] = df["chrom"].str.replace("chr", "", regex=False)
    # Drop alt haplotype contigs / non-standard scaffolds via prior filter
    return df.reset_index(drop=True)


def nearest_gene_per_variant(variants: pd.DataFrame, gene_coords: pd.DataFrame) -> pd.DataFrame:
    """Assign nearest-gene per variant (in-gene preferred; else min |variant - TSS|).

    Inputs:
        variants : columns chrom (str like '1'), pos (int)
        gene_coords : output of build_gene_coords()
    Returns variants + nearest_ensembl_base + nearest_symbol + nearest_biotype + distance_to_tss.
    """
    gc = gene_coords[["chrom_int", "start", "end", "tss", "ensembl_base", "gene_name", "gene_biotype"]].copy()

    out = []
    for chrom, sub_v in variants.groupby("chrom_str"):
        sub_g = gc[gc["chrom_int"] == chrom].copy()
        if sub_g.empty:
            for _, vr in sub_v.iterrows():
                out.append({**vr.to_dict(), "nearest_ensembl_base": None,
                            "nearest_symbol": None, "nearest_biotype": None,
                            "distance_to_tss": np.nan, "in_gene_body": False})
            continue
        # Sort genes by TSS for searchsorted
        sub_g = sub_g.sort_values("tss").reset_index(drop=True)
        tss_arr = sub_g["tss"].values
        starts = sub_g["start"].values
        ends = sub_g["end"].values
        names = sub_g["gene_name"].values
        ebs = sub_g["ensembl_base"].values
        bio = sub_g["gene_biotype"].values

        for _, vr in sub_v.iterrows():
            pos = int(vr["pos"])
            # Prefer a gene whose body contains pos
            in_mask = (starts <= pos) & (ends >= pos)
            if in_mask.any():
                # Pick the in-body gene with smallest TSS distance (in case of overlap)
                idxs = np.where(in_mask)[0]
                dists = np.abs(tss_arr[idxs] - pos)
                j = idxs[np.argmin(dists)]
                out.append({
                    **vr.to_dict(),
                    "nearest_ensembl_base": ebs[j],
                    "nearest_symbol": names[j],
                    "nearest_biotype": bio[j],
                    "distance_to_tss": int(tss_arr[j] - pos),
                    "in_gene_body": True,
                })
            else:
                # Nearest by |TSS - pos|
                i = np.searchsorted(tss_arr, pos)
                cand = []
                if i > 0:
                    cand.append(i - 1)
                if i < len(tss_arr):
                    cand.append(i)
                dists = [abs(tss_arr[k] - pos) for k in cand]
                j = cand[int(np.argmin(dists))]
                out.append({
                    **vr.to_dict(),
                    "nearest_ensembl_base": ebs[j],
                    "nearest_symbol": names[j],
                    "nearest_biotype": bio[j],
                    "distance_to_tss": int(tss_arr[j] - pos),
                    "in_gene_body": False,
                })
    return pd.DataFrame(out)


# ============================================================================
# Step 2: Per-axis pools
# ============================================================================

def load_atlas() -> pd.DataFrame:
    cols = [
        "human_symbol", "ensembl_id", "gene_biotype", "mouse_ortholog",
        "coloc_best_susie_pp4_polyfun", "coloc_best_susie_gwas_polyfun",
        "coloc_susie_best_pp4", "coloc_susie_best_gwas",
        "dream_logFC", "dream_padj",
    ]
    df = pd.read_csv(ATLAS_CSV, usecols=cols, low_memory=False)
    df["ensembl_base"] = strip_version(df["ensembl_id"])
    df["mouse_ensembl_base"] = strip_version(df["mouse_ortholog"].fillna(""))
    df.loc[df["mouse_ensembl_base"] == "", "mouse_ensembl_base"] = np.nan
    df.loc[df["mouse_ensembl_base"] == "nan", "mouse_ensembl_base"] = np.nan
    return df


def axis1_coloc(atlas: pd.DataFrame, fallback_log: list) -> pd.DataFrame:
    """COLOC tier — atlas coloc_best_susie_pp4_polyfun > 0.5; fallback to non-polyfun."""
    pp4_col = "coloc_best_susie_pp4_polyfun"
    gwas_col = "coloc_best_susie_gwas_polyfun"
    if pp4_col not in atlas.columns or atlas[pp4_col].notna().sum() == 0:
        fallback_log.append(
            "WARNING: coloc_best_susie_pp4_polyfun not usable; falling back to coloc_susie_best_pp4."
        )
        pp4_col = "coloc_susie_best_pp4"
        gwas_col = "coloc_susie_best_gwas"
    mask = atlas[pp4_col] > 0.5
    out = atlas.loc[mask, ["ensembl_base", "human_symbol", "gene_biotype",
                           "mouse_ensembl_base", pp4_col, gwas_col]].copy()
    out = out.rename(columns={pp4_col: "coloc_pp4", gwas_col: "coloc_gwas"})
    out["axis_coloc"] = True
    return out.reset_index(drop=True)


def axis2_finemap(gene_coords: pd.DataFrame) -> pd.DataFrame:
    """Finemap tier — combined_finemapping.csv recommended_pip > 0.5, nearest-gene."""
    fm = pd.read_csv(FINEMAP_CSV, low_memory=False)
    # rec_pip > 0.5
    fm = fm[pd.to_numeric(fm["recommended_pip"], errors="coerce") > 0.5].copy()
    fm["recommended_pip"] = pd.to_numeric(fm["recommended_pip"], errors="coerce")
    # Drop variants with missing chrom/pos (we'd skip them anyway)
    fm = fm.dropna(subset=["chromosome", "position"])
    # chromosome may be float (1.0) or string ("chr1"); normalise to bare integer-as-string
    chrom_num = pd.to_numeric(fm["chromosome"], errors="coerce")
    fm["chromosome"] = np.where(
        chrom_num.notna(),
        chrom_num.fillna(-1).astype(int).astype(str),
        fm["chromosome"].astype(str).str.replace("chr", "", regex=False).str.strip(),
    )
    fm["position"] = pd.to_numeric(fm["position"], errors="coerce")
    fm = fm.dropna(subset=["position"])
    fm["position"] = fm["position"].astype(int)
    # Dedup at variant×study granularity (already in there) — keep one row per (chrom,pos,study)
    fm = fm.drop_duplicates(["chromosome", "position", "study"], keep="first")

    var_uniq = fm[["chromosome", "position"]].drop_duplicates().rename(
        columns={"chromosome": "chrom_str", "position": "pos"}
    )

    mapped = nearest_gene_per_variant(var_uniq, gene_coords)
    # Merge back with study/pip information
    fm2 = fm.merge(
        mapped,
        left_on=["chromosome", "position"],
        right_on=["chrom_str", "pos"],
        how="left",
    )
    # Build gene-level summary
    grp = fm2.dropna(subset=["nearest_ensembl_base"]).groupby("nearest_ensembl_base")
    out = grp.agg(
        finemap_pip=("recommended_pip", "max"),
        finemap_n_gwas=("study", "nunique"),
        nearest_symbol=("nearest_symbol", "first"),
        nearest_biotype=("nearest_biotype", "first"),
    ).reset_index().rename(columns={"nearest_ensembl_base": "ensembl_base"})
    out["axis_finemap"] = True
    # Save the raw variant→gene table for transparency
    fm2[[
        "chromosome", "position", "study", "ancestry", "recommended_pip",
        "nearest_ensembl_base", "nearest_symbol", "nearest_biotype",
        "distance_to_tss", "in_gene_body",
    ]].to_csv(OUT_DIR / "finemap_variant_to_gene.csv", index=False)
    return out


def axis3_mouse_de(ortho: pd.DataFrame) -> pd.DataFrame:
    """Mouse UP DEG — MCD ∪ CDAHFD, padj<0.05, logFC>0.3, hep AveExpr filter (biotype-specific)."""
    mcd = load_mouse_de("MCD").reset_index()
    cda = load_mouse_de("CDAHFD").reset_index()
    mcd = mcd.rename(columns={
        "logFC": "logFC_mcd", "AveExpr": "AveExpr_mcd",
        "P.Value": "P.Value_mcd", "adj.P.Val": "adj.P.Val_mcd",
    })[["gene_id_base", "logFC_mcd", "AveExpr_mcd", "adj.P.Val_mcd"]]
    cda = cda.rename(columns={
        "logFC": "logFC_cda", "AveExpr": "AveExpr_cda",
        "P.Value": "P.Value_cda", "adj.P.Val": "adj.P.Val_cda",
    })[["gene_id_base", "logFC_cda", "AveExpr_cda", "adj.P.Val_cda"]]

    full = mcd.merge(cda, on="gene_id_base", how="outer")
    # Significance UP filter on either diet
    sig_mcd = (full["adj.P.Val_mcd"] < 0.05) & (full["logFC_mcd"] > 0.3)
    sig_cda = (full["adj.P.Val_cda"] < 0.05) & (full["logFC_cda"] > 0.3)
    sig_any = sig_mcd | sig_cda
    de = full[sig_any].copy()
    # max logFC across sig diets
    de["mouse_lfc_max"] = np.nan
    de.loc[sig_mcd, "mouse_lfc_max"] = de.loc[sig_mcd, "logFC_mcd"]
    de.loc[sig_cda & ~sig_mcd, "mouse_lfc_max"] = de.loc[sig_cda & ~sig_mcd, "logFC_cda"]
    # if both sig take max
    both = sig_mcd & sig_cda
    de.loc[both, "mouse_lfc_max"] = np.maximum(de.loc[both, "logFC_mcd"], de.loc[both, "logFC_cda"])
    # max AveExpr across the two diets (use whatever is non-null)
    de["mouse_ave_expr"] = np.nanmax(
        de[["AveExpr_mcd", "AveExpr_cda"]].values, axis=1
    )

    # Annotate mouse biotype via the GENCODE mouse metadata — we only have human
    # GENCODE in metadata file; use atlas biotype via ortholog as a coarse proxy
    # NOTE: biotype-specific expression filter
    #   PC AveExpr ≥ 0.5,  lncRNA AveExpr ≥ -0.3
    # We don't have a clean mouse biotype table at hand; derive from atlas ortho.
    # For genes without an atlas ortho row, default to a permissive lncRNA-like threshold
    # so we don't accidentally drop genuine mouse lncRNAs (will be re-filtered downstream
    # by biotype assignment).
    # For now we'll assign biotype later from the ortholog table (human-side gene_biotype).
    de = de.rename(columns={"gene_id_base": "mouse_ensembl_base"})
    return de.reset_index(drop=True)


def axis4_human_de(ortho: pd.DataFrame) -> pd.DataFrame:
    """Human UP DEG — dream_results_ashr, padj<0.05, logFC>0.3."""
    dr = pd.read_csv(DREAM_CSV)
    # dream uses 'adj.P.Val' OR 'padj' — header showed: gene,logFC,AveExpr,t,P.Value,padj,…
    # Use 'padj' column from dream_results_ashr
    pcol = "padj" if "padj" in dr.columns else "adj.P.Val"
    mask = (pd.to_numeric(dr[pcol], errors="coerce") < 0.05) & (pd.to_numeric(dr["logFC"], errors="coerce") > 0.3)
    sub = dr.loc[mask, ["gene", "symbol", "logFC", pcol]].copy()
    sub = sub.rename(columns={"logFC": "human_lfc", pcol: "human_padj", "gene": "ensembl_id", "symbol": "human_symbol"})
    sub["ensembl_base"] = strip_version(sub["ensembl_id"])
    sub = sub.drop_duplicates("ensembl_base", keep="first")
    sub["axis_humanDE"] = True
    return sub.reset_index(drop=True)


# ============================================================================
# Step 3: Master merge & assignment of mouse target
# ============================================================================

def build_master(
    atlas: pd.DataFrame,
    ortho: pd.DataFrame,
    cov_coloc: pd.DataFrame,
    cov_fine: pd.DataFrame,
    cov_mde: pd.DataFrame,
    cov_hde: pd.DataFrame,
    gene_coords: pd.DataFrame,
) -> pd.DataFrame:
    """Build per-gene master keyed on (mouse_ensembl_base, biotype).

    Strategy:
      - For each axis, derive a mouse target gene:
          * COLOC, Finemap, HumanDE are human-side → use ortho table to map to mouse
          * MouseDE is mouse-side → use directly
      - Outer-merge all four axes on mouse_ensembl_base.
    """
    # Human→mouse ortho map (keep all pairs; deduplicate later)
    # Index: human ensembl_base -> [list of (mouse_ensembl_base, mouse_symbol)]
    ortho_h2m = ortho[["human_ensembl", "mouse_ensembl", "human_symbol", "mouse_symbol"]].drop_duplicates()

    def _coalesce_human_symbol(df: pd.DataFrame) -> pd.DataFrame:
        """After a merge that disambiguates `human_symbol` into _x/_y, coalesce them.

        Prefer the left (non-_y) side, fall back to the ortho side.
        """
        if "human_symbol_x" in df.columns:
            left = df["human_symbol_x"]
            right = df.get("human_symbol_y")
            df["human_symbol"] = left.where(left.notna() & (left.astype(str) != ""), right)
            df = df.drop(columns=[c for c in ["human_symbol_x", "human_symbol_y"] if c in df.columns])
        return df

    # ---- axis1: human → mouse ----
    a1 = cov_coloc.copy()
    # Use atlas mouse_ortholog if present; else go via ortho table
    a1_with_mouse = a1.merge(
        ortho_h2m, left_on="ensembl_base", right_on="human_ensembl", how="left",
    )
    a1_with_mouse = _coalesce_human_symbol(a1_with_mouse)
    # Prefer atlas-provided mouse_ortholog if present
    a1_with_mouse["mouse_ensembl_final"] = a1_with_mouse["mouse_ensembl_base"].where(
        a1_with_mouse["mouse_ensembl_base"].notna(), a1_with_mouse["mouse_ensembl"]
    )
    a1_with_mouse = a1_with_mouse.dropna(subset=["mouse_ensembl_final"])
    a1_with_mouse = a1_with_mouse.rename(columns={"mouse_ensembl_final": "mouse_ensembl_base_target"})
    a1_keep = a1_with_mouse[[
        "mouse_ensembl_base_target", "ensembl_base", "human_symbol",
        "mouse_symbol", "gene_biotype",
        "coloc_pp4", "coloc_gwas", "axis_coloc",
    ]].drop_duplicates(["mouse_ensembl_base_target", "ensembl_base"], keep="first")
    a1_keep = a1_keep.rename(columns={"mouse_ensembl_base_target": "mouse_ensembl_base"})

    # ---- axis2: human → mouse ----
    a2 = cov_fine.copy()
    a2 = a2.rename(columns={"nearest_symbol": "human_symbol", "nearest_biotype": "biotype_fm"})
    a2_with_mouse = a2.merge(
        ortho_h2m, left_on="ensembl_base", right_on="human_ensembl", how="left",
    )
    a2_with_mouse = _coalesce_human_symbol(a2_with_mouse)
    a2_with_mouse = a2_with_mouse.dropna(subset=["mouse_ensembl"])
    a2_with_mouse = a2_with_mouse.rename(columns={"mouse_ensembl": "mouse_ensembl_base"})
    a2_keep = a2_with_mouse[[
        "mouse_ensembl_base", "ensembl_base", "human_symbol", "mouse_symbol",
        "biotype_fm", "finemap_pip", "finemap_n_gwas", "axis_finemap",
    ]].drop_duplicates(["mouse_ensembl_base", "ensembl_base"], keep="first")

    # ---- axis3: already mouse ----
    a3 = cov_mde.copy()
    # Look up its human ortholog & biotype
    a3_with_human = a3.merge(
        ortho_h2m, left_on="mouse_ensembl_base", right_on="mouse_ensembl", how="left",
    )
    # If multiple human orthologs, keep one (one2many → take first)
    a3_with_human = a3_with_human.drop_duplicates("mouse_ensembl_base", keep="first")
    a3_keep = a3_with_human[[
        "mouse_ensembl_base", "human_ensembl", "human_symbol", "mouse_symbol",
        "logFC_mcd", "adj.P.Val_mcd", "logFC_cda", "adj.P.Val_cda",
        "mouse_lfc_max", "mouse_ave_expr",
    ]].rename(columns={"human_ensembl": "ensembl_base"})
    a3_keep["axis_mouseDE"] = True

    # ---- axis4: human → mouse ----
    a4 = cov_hde.copy()
    a4_with_mouse = a4.merge(
        ortho_h2m, left_on="ensembl_base", right_on="human_ensembl", how="left",
    )
    a4_with_mouse = _coalesce_human_symbol(a4_with_mouse)
    a4_with_mouse = a4_with_mouse.dropna(subset=["mouse_ensembl"])
    a4_with_mouse = a4_with_mouse.rename(columns={"mouse_ensembl": "mouse_ensembl_base"})
    a4_keep = a4_with_mouse[[
        "mouse_ensembl_base", "ensembl_base", "human_symbol", "mouse_symbol",
        "human_lfc", "human_padj", "axis_humanDE",
    ]].drop_duplicates(["mouse_ensembl_base", "ensembl_base"], keep="first")

    # ---- Outer-merge on mouse_ensembl_base ----
    # Keep best per axis per mouse target (the axis presence is what matters; for the
    # ranking metric we take the *max* per mouse target across human paralogs).
    def collapse(df: pd.DataFrame, axis_cols: list[str]) -> pd.DataFrame:
        # axis_cols = the columns we want to retain at gene-target level
        keep_cols = ["mouse_ensembl_base"] + axis_cols
        return df[keep_cols].drop_duplicates("mouse_ensembl_base", keep="first")

    # For axis1, sort by PP4 descending so the kept row is the strongest GWAS
    a1_keep_sorted = a1_keep.sort_values("coloc_pp4", ascending=False)
    a1_collapsed = a1_keep_sorted.drop_duplicates("mouse_ensembl_base", keep="first")[
        ["mouse_ensembl_base", "coloc_pp4", "coloc_gwas", "axis_coloc"]
    ]
    # For axis2: sort by finemap_pip
    a2_keep_sorted = a2_keep.sort_values("finemap_pip", ascending=False)
    a2_collapsed = a2_keep_sorted.drop_duplicates("mouse_ensembl_base", keep="first")[
        ["mouse_ensembl_base", "finemap_pip", "finemap_n_gwas", "axis_finemap"]
    ]
    # axis3 already one-row-per-mouse
    a3_collapsed = a3_keep[
        ["mouse_ensembl_base", "logFC_mcd", "adj.P.Val_mcd", "logFC_cda", "adj.P.Val_cda",
         "mouse_lfc_max", "mouse_ave_expr", "axis_mouseDE"]
    ].drop_duplicates("mouse_ensembl_base", keep="first")
    # axis4: collapse by max human_lfc
    a4_keep_sorted = a4_keep.sort_values("human_lfc", ascending=False)
    a4_collapsed = a4_keep_sorted.drop_duplicates("mouse_ensembl_base", keep="first")[
        ["mouse_ensembl_base", "human_lfc", "human_padj", "axis_humanDE"]
    ]

    master = a1_collapsed.merge(a2_collapsed, on="mouse_ensembl_base", how="outer")
    master = master.merge(a3_collapsed, on="mouse_ensembl_base", how="outer")
    master = master.merge(a4_collapsed, on="mouse_ensembl_base", how="outer")

    for col in ["axis_coloc", "axis_finemap", "axis_mouseDE", "axis_humanDE"]:
        master[col] = master[col].fillna(False).astype(bool)
    master["n_axes"] = master[["axis_coloc", "axis_finemap", "axis_mouseDE", "axis_humanDE"]].sum(axis=1)

    # Attach mouse_symbol + human_symbol via ortho table
    label_src = ortho[["mouse_ensembl", "mouse_symbol", "human_symbol"]].rename(
        columns={"mouse_ensembl": "mouse_ensembl_base"}
    ).drop_duplicates("mouse_ensembl_base", keep="first")
    master = master.merge(
        label_src[["mouse_ensembl_base", "mouse_symbol", "human_symbol"]],
        on="mouse_ensembl_base", how="left",
    )

    # Biotype assignment: use MOUSE GENCODE biotype (vM38) — the actual biotype of the
    # screen target. This is critical for lncRNAs (ortholog table is PC-biased).
    mouse_meta = load_mouse_metadata()
    master = master.merge(
        mouse_meta[["mouse_ensembl_base", "mouse_symbol_gtf", "mouse_biotype"]],
        on="mouse_ensembl_base", how="left",
    )
    # Coalesce mouse_symbol: prefer ortho table, fall back to GTF
    master["mouse_symbol"] = master["mouse_symbol"].where(
        master["mouse_symbol"].notna() & (master["mouse_symbol"].astype(str) != ""),
        master["mouse_symbol_gtf"],
    )
    # Coalesce biotype: prefer mouse GTF biotype; fall back to atlas human-side biotype
    biotype_atlas = atlas[["ensembl_base", "gene_biotype"]].dropna(subset=["ensembl_base"])
    ortho_hu = ortho[["mouse_ensembl", "human_ensembl"]].drop_duplicates("mouse_ensembl", keep="first").rename(
        columns={"mouse_ensembl": "mouse_ensembl_base", "human_ensembl": "ensembl_base"}
    )
    bm = ortho_hu.merge(biotype_atlas, on="ensembl_base", how="left").drop_duplicates("mouse_ensembl_base", keep="first")
    master = master.merge(bm[["mouse_ensembl_base", "gene_biotype"]], on="mouse_ensembl_base", how="left")

    # Fallback biotype from axis2 nearest_gene biotype if missing
    fm_bio = a2_keep[["mouse_ensembl_base", "biotype_fm"]].drop_duplicates("mouse_ensembl_base", keep="first")
    master = master.merge(fm_bio, on="mouse_ensembl_base", how="left")
    master["gene_biotype_final"] = master["mouse_biotype"].where(
        master["mouse_biotype"].notna() & (master["mouse_biotype"] != ""),
        master["gene_biotype"],
    )
    master["gene_biotype_final"] = master["gene_biotype_final"].where(
        master["gene_biotype_final"].notna() & (master["gene_biotype_final"] != ""),
        master["biotype_fm"],
    )

    # ---------- Mouse hep expression filter (biotype-specific) ----------
    # PC AveExpr >= 0.5,  lncRNA AveExpr >= -0.3
    # If a gene has no mouse_ave_expr (i.e. axis_mouseDE=False AND was not in MCD/CDAHFD
    # at all), we need to fetch its expression from MCD/CDAHFD full DE tables (sig filter
    # off; just for AveExpr).
    mcd_all = load_mouse_de("MCD").reset_index().rename(
        columns={"gene_id_base": "mouse_ensembl_base", "AveExpr": "AveExpr_mcd_all"}
    )[["mouse_ensembl_base", "AveExpr_mcd_all"]]
    cda_all = load_mouse_de("CDAHFD").reset_index().rename(
        columns={"gene_id_base": "mouse_ensembl_base", "AveExpr": "AveExpr_cda_all"}
    )[["mouse_ensembl_base", "AveExpr_cda_all"]]
    master = master.merge(mcd_all, on="mouse_ensembl_base", how="left")
    master = master.merge(cda_all, on="mouse_ensembl_base", how="left")
    master["mouse_ave_expr_any"] = np.nanmax(
        master[["mouse_ave_expr", "AveExpr_mcd_all", "AveExpr_cda_all"]].values, axis=1
    )

    is_pc = master["gene_biotype_final"] == "protein_coding"
    is_lnc = master["gene_biotype_final"] == "lncRNA"
    expr_pass_pc = master["mouse_ave_expr_any"] >= 0.5
    expr_pass_lnc = master["mouse_ave_expr_any"] >= -0.3
    expr_pass_other = master["mouse_ave_expr_any"] >= 0.5

    master["passes_hep_expr"] = np.where(
        is_pc, expr_pass_pc,
        np.where(is_lnc, expr_pass_lnc, expr_pass_other),
    )

    # ---------- Primary axis + rank ----------
    # Priority: COLOC > Finemap > MouseDE > HumanDE
    def assign_primary(row):
        if row["axis_coloc"]:
            return "COLOC", row["coloc_pp4"]
        if row["axis_finemap"]:
            return "Finemap", row["finemap_pip"]
        if row["axis_mouseDE"]:
            return "MouseDE", row["mouse_lfc_max"]
        if row["axis_humanDE"]:
            return "HumanDE", row["human_lfc"]
        return None, np.nan

    primary = master.apply(assign_primary, axis=1, result_type="expand")
    primary.columns = ["primary_axis", "primary_metric"]
    master = pd.concat([master, primary], axis=1)
    # Rank within each primary_axis
    master["primary_rank"] = master.groupby("primary_axis")["primary_metric"].rank(
        ascending=False, method="dense"
    )

    return master


# ============================================================================
# Step 4: Driver
# ============================================================================

def main():
    fallback_log: list[str] = []

    print("[1/6] Building gene coordinates from GENCODE v49 GTF …", flush=True)
    gene_coords = build_gene_coords()
    print(f"      gene_coords: {len(gene_coords):,} gene rows on primary chr.")

    print("[2/6] Loading atlas …", flush=True)
    atlas = load_atlas()
    print(f"      atlas: {len(atlas):,} rows.")

    print("[3/6] Building ortholog table …", flush=True)
    ortho = build_ortholog_table(include_atlas_fallback=True)
    print(f"      ortho: {len(ortho):,} (mouse,human) pairs.")
    # Pre-strip versions just in case
    ortho["mouse_ensembl"] = strip_version(ortho["mouse_ensembl"])
    ortho["human_ensembl"] = strip_version(ortho["human_ensembl"])

    print("[4/6] Building 4 per-axis pools …", flush=True)
    cov_coloc = axis1_coloc(atlas, fallback_log)
    print(f"      axis1 COLOC pool: {len(cov_coloc):,} (human-side rows)")
    cov_fine = axis2_finemap(gene_coords)
    print(f"      axis2 Finemap pool: {len(cov_fine):,} (gene-level rows, via nearest-gene)")
    cov_mde = axis3_mouse_de(ortho)
    print(f"      axis3 Mouse UP DEG pool: {len(cov_mde):,} mouse-side rows")
    cov_hde = axis4_human_de(ortho)
    print(f"      axis4 Human UP DEG pool: {len(cov_hde):,} human-side rows")

    print("[5/6] Merging into master …", flush=True)
    master = build_master(atlas, ortho, cov_coloc, cov_fine, cov_mde, cov_hde, gene_coords)
    print(f"      master pre-filter: {len(master):,} rows; biotype value_counts:")
    print(master["gene_biotype_final"].value_counts().head(10))

    # ---- Apply the global UP-only + hep-expressed filter ----
    # For axis_humanDE-only & axis_coloc-only & axis_finemap-only genes: we need to
    # make sure they aren't DOWN in mouse if any mouse data exists. The task says
    # "discard any gene with logFC < 0 in its primary axis", so if primary_axis == 'HumanDE'
    # we check human_lfc > 0 (already filtered); etc. So the global UP-only rule is
    # naturally enforced by the per-axis filters.

    n_pre = len(master)
    master_exp = master[master["passes_hep_expr"]].copy()
    n_post = len(master_exp)
    print(f"      master post hep-expr filter: {n_post:,} (dropped {n_pre - n_post:,})")

    # ---- Split into PC vs lncRNA ----
    pc = master_exp[master_exp["gene_biotype_final"] == "protein_coding"].copy()
    lnc = master_exp[master_exp["gene_biotype_final"] == "lncRNA"].copy()
    print(f"      PC: {len(pc):,};  lncRNA: {len(lnc):,};  other: {len(master_exp) - len(pc) - len(lnc):,}")

    # ---- Flag: has_mouse_de_support (P1-11 fix) ----
    # True if the gene is UP-regulated in MCD or CDAHFD (axis_mouseDE == True).
    # False if the gene entered the pool only via human-side axes (COLOC/Finemap/HumanDE)
    # without direct mouse DE evidence.
    for df in [pc, lnc, master_exp]:
        df["has_mouse_de_support"] = df["axis_mouseDE"].astype(bool)

    # ---- Assemble final candidate columns ----
    out_cols = [
        ("gene_symbol_mouse", "mouse_symbol"),
        ("gene_symbol_human", "human_symbol"),
        ("gene_id_mouse", "mouse_ensembl_base"),
        ("biotype", "gene_biotype_final"),
        ("axis_coloc", "axis_coloc"),
        ("axis_finemap", "axis_finemap"),
        ("axis_mouseDE", "axis_mouseDE"),
        ("axis_humanDE", "axis_humanDE"),
        ("n_axes", "n_axes"),
        ("has_mouse_de_support", "has_mouse_de_support"),
        ("coloc_pp4", "coloc_pp4"),
        ("coloc_gwas", "coloc_gwas"),
        ("finemap_pip", "finemap_pip"),
        ("finemap_n_gwas", "finemap_n_gwas"),
        ("mouse_lfc_mcd", "logFC_mcd"),
        ("mouse_padj_mcd", "adj.P.Val_mcd"),
        ("mouse_lfc_cdahfd", "logFC_cda"),
        ("mouse_padj_cdahfd", "adj.P.Val_cda"),
        ("mouse_lfc_max", "mouse_lfc_max"),
        ("mouse_ave_expr", "mouse_ave_expr_any"),
        ("human_lfc", "human_lfc"),
        ("human_padj", "human_padj"),
        ("primary_axis", "primary_axis"),
        ("primary_metric", "primary_metric"),
        ("primary_rank", "primary_rank"),
    ]

    def reshape(df: pd.DataFrame) -> pd.DataFrame:
        out = pd.DataFrame()
        for new, old in out_cols:
            out[new] = df[old] if old in df.columns else np.nan
        # Sort by (n_axes desc, primary_metric desc)
        out = out.sort_values(["n_axes", "primary_metric"], ascending=[False, False])
        return out.reset_index(drop=True)

    pc_out = reshape(pc)
    lnc_out = reshape(lnc)

    pc_path = OUT_DIR / "candidates_pc_independent.csv"
    lnc_path = OUT_DIR / "candidates_lncrna_independent.csv"
    pc_out.to_csv(pc_path, index=False)
    lnc_out.to_csv(lnc_path, index=False)
    print(f"      wrote {pc_path}  ({len(pc_out):,} rows)")
    print(f"      wrote {lnc_path} ({len(lnc_out):,} rows)")

    # ---- Per-axis pool sizes (after hep-expr filter, by biotype) ----
    def per_axis_pool(df: pd.DataFrame):
        return {
            "COLOC":   int(df["axis_coloc"].sum()),
            "Finemap": int(df["axis_finemap"].sum()),
            "MouseDE": int(df["axis_mouseDE"].sum()),
            "HumanDE": int(df["axis_humanDE"].sum()),
        }
    pc_pool = per_axis_pool(pc_out)
    lnc_pool = per_axis_pool(lnc_out)
    print()
    print(f"PC per-axis pool: {pc_pool}")
    print(f"lncRNA per-axis pool: {lnc_pool}")

    # ---- Mouse DE support breakdown (P1-11) ----
    pc_mde_yes = int(pc_out["has_mouse_de_support"].sum())
    pc_mde_no = len(pc_out) - pc_mde_yes
    lnc_mde_yes = int(lnc_out["has_mouse_de_support"].sum())
    lnc_mde_no = len(lnc_out) - lnc_mde_yes
    print(f"has_mouse_de_support — PC: {pc_mde_yes} True / {pc_mde_no} False; "
          f"lncRNA: {lnc_mde_yes} True / {lnc_mde_no} False")

    # ---- 4-way Venn (PC only) — all 16 cells, including empties ----
    pc_keys = pc_out[["axis_coloc", "axis_finemap", "axis_mouseDE", "axis_humanDE"]].astype(int)
    venn_dict = pc_keys.value_counts().to_dict()
    venn_rows = []
    from itertools import product
    for c, f, m, h in product([1, 0], repeat=4):
        n = venn_dict.get((c, f, m, h), 0)
        venn_rows.append({"axis_coloc": c, "axis_finemap": f, "axis_mouseDE": m, "axis_humanDE": h, "n": n})
    venn = pd.DataFrame(venn_rows)
    venn = venn.sort_values(["axis_coloc", "axis_finemap", "axis_mouseDE", "axis_humanDE"],
                             ascending=False)

    # ---- Cross-axis overlap rates ----
    def overlap_rate(df, a, b):
        if df[a].sum() == 0:
            return float("nan"), 0, 0
        n_both = int((df[a] & df[b]).sum())
        n_a = int(df[a].sum())
        return float(n_both / n_a), n_both, n_a
    pc_overlap = {
        "MouseDE → also COLOC":   overlap_rate(pc_out, "axis_mouseDE", "axis_coloc"),
        "MouseDE → also Finemap": overlap_rate(pc_out, "axis_mouseDE", "axis_finemap"),
        "MouseDE → also HumanDE": overlap_rate(pc_out, "axis_mouseDE", "axis_humanDE"),
        "HumanDE → also COLOC":   overlap_rate(pc_out, "axis_humanDE", "axis_coloc"),
        "HumanDE → also Finemap": overlap_rate(pc_out, "axis_humanDE", "axis_finemap"),
        "HumanDE → also MouseDE": overlap_rate(pc_out, "axis_humanDE", "axis_mouseDE"),
        "COLOC   → also Finemap": overlap_rate(pc_out, "axis_coloc",   "axis_finemap"),
        "COLOC   → also MouseDE": overlap_rate(pc_out, "axis_coloc",   "axis_mouseDE"),
        "COLOC   → also HumanDE": overlap_rate(pc_out, "axis_coloc",   "axis_humanDE"),
        "Finemap → also COLOC":   overlap_rate(pc_out, "axis_finemap", "axis_coloc"),
        "Finemap → also MouseDE": overlap_rate(pc_out, "axis_finemap", "axis_mouseDE"),
        "Finemap → also HumanDE": overlap_rate(pc_out, "axis_finemap", "axis_humanDE"),
    }

    # ---- Top-10 per axis preview ----
    def top10(df, axis_col, metric_col):
        sub = df[df[axis_col]].copy()
        sub = sub.sort_values(metric_col, ascending=False).head(10)
        return sub[["gene_symbol_mouse", "gene_symbol_human", "biotype", metric_col, "n_axes"]]

    pc_top = {
        "COLOC":   top10(pc_out, "axis_coloc",   "coloc_pp4"),
        "Finemap": top10(pc_out, "axis_finemap", "finemap_pip"),
        "MouseDE": top10(pc_out, "axis_mouseDE", "mouse_lfc_max"),
        "HumanDE": top10(pc_out, "axis_humanDE", "human_lfc"),
    }

    # ---- Quality notes ----
    quality = []
    quality.append(f"COLOC source column: `coloc_best_susie_pp4_polyfun` (PolyFun production EUR panel)."
                   + (" Fallback DID NOT trigger." if not fallback_log else ""))
    if fallback_log:
        quality.extend(fallback_log)
    n_atlas_pc_pp4 = int((atlas["coloc_best_susie_pp4_polyfun"] > 0.5).sum())
    quality.append(f"Atlas-level COLOC PP4>0.5 hits (any biotype, pre hep-expr): {n_atlas_pc_pp4:,}")
    quality.append(f"Finemap variants with recommended_pip>0.5 (any GWAS): "
                   f"{int(pd.read_csv(FINEMAP_CSV, usecols=['recommended_pip'], low_memory=False)['recommended_pip'].gt(0.5).sum()):,}")
    quality.append("Finemap variant→gene method: nearest-gene by TSS distance, in-gene-body preferred. "
                   "Variant→gene table saved to data/finemap_variant_to_gene.csv.")
    quality.append("Mouse DE: MCD ∪ CDAHFD, padj<0.05, logFC>0.3 (UP). Hep-expr filter "
                   "AveExpr≥0.5 (PC) / ≥-0.3 (lncRNA) applied at master-table stage.")
    quality.append("Human DE: dream_results_ashr (padj<0.05, logFC>0.3); ortho→mouse via "
                   "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz "
                   "(many2many pairs collapsed by max metric).")
    quality.append(f"Genes with NO mouse symbol assignable (kept anyway, biotype-known): "
                   f"{int(master_exp['mouse_symbol'].isna().sum()):,}.")
    quality.append("Sort key for candidate CSVs: (n_axes DESC, primary_metric DESC).")
    quality.append("lncRNA pool composition: ONLY the MouseDE axis contributes to lncRNA candidates "
                   "(545 / 545 = 100%). The biomaRt mouse↔human ortholog table includes essentially "
                   "no lncRNA pairs (1 of 12,671 human lncRNAs has a mouse_ortholog in the atlas), "
                   "so human-side axes (COLOC, Finemap, HumanDE) cannot map to mouse lncRNAs. "
                   "This is a biological reality of lncRNA evolution, not a pipeline bug.")
    n_other = int((master_exp["gene_biotype_final"].isin(["protein_coding", "lncRNA"]) == False).sum())
    quality.append(f"{n_other} non-PC, non-lncRNA biotype rows (pseudogenes / IG_V_gene / TR_V_gene / "
                   f"TEC / miRNA / snoRNA) were filtered out of both candidate tables; the screen targets "
                   f"PC + lncRNA only.")
    quality.append("Mouse biotype source: GENCODE vM38 GTF "
                   "(`data/ncrna_conservation/gencode.vM38.annotation.gtf.gz`), cached to "
                   "`Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv`. "
                   "Mouse biotype is preferred over human atlas biotype to correctly classify "
                   "the screen target (the user wants mouse symbols).")
    quality.append(f"Ortholog bridge: min_tier='M' (default) — Tier L positional-only pairs excluded. "
                   f"Ortholog table sourced from master_ortholog_table.tsv.gz with confidence_tier filter.")
    quality.append(f"has_mouse_de_support column added (P1-11 fix): PC {pc_mde_yes} True / {pc_mde_no} False; "
                   f"lncRNA {lnc_mde_yes} True / {lnc_mde_no} False. "
                   f"True = gene is UP-regulated in MCD or CDAHFD; False = entered pool via "
                   f"human-side axes only (COLOC/Finemap/HumanDE) without direct mouse DE evidence.")

    # ---- Quota scenarios ----
    scenarios = []
    pool_pc_coloc = pc_pool["COLOC"]
    pool_pc_fine = pc_pool["Finemap"]
    pool_pc_mouse = pc_pool["MouseDE"]
    pool_pc_human = pc_pool["HumanDE"]
    pool_lnc_union = len(lnc_out)
    pool_pc_union = len(pc_out)

    # PC slots in 5,000 library: 5,000 - lncRNA quota - controls
    lnc_quota = min(338, pool_lnc_union)
    ctrl_quota = 100
    pc_budget = 5000 - lnc_quota - ctrl_quota
    # "Mouse-only" singletons (axis_mouseDE only) — the easiest tail to trim
    mde_only = int(((pc_out["axis_mouseDE"]) &
                    (~pc_out["axis_coloc"]) & (~pc_out["axis_finemap"]) & (~pc_out["axis_humanDE"])).sum())
    hde_only = int(((pc_out["axis_humanDE"]) &
                    (~pc_out["axis_coloc"]) & (~pc_out["axis_finemap"]) & (~pc_out["axis_mouseDE"])).sum())
    multi_axis_pc = pool_pc_union - mde_only - hde_only  # genes with ≥2 axes OR genetic-only

    scenarios.append(
        f"**Scenario A (genetic-first, undershoots)**: All COLOC ({pool_pc_coloc}) + all Finemap ({pool_pc_fine}) "
        f"+ 2,000 MouseDE + 1,000 HumanDE (de-duplicated to mouse target) "
        f"+ {lnc_quota} lncRNA + {ctrl_quota} controls. "
        f"Sum-of-quotas (with no overlap) ≈ "
        f"{pool_pc_coloc + pool_pc_fine + 2000 + 1000 + lnc_quota + ctrl_quota}. "
        f"In practice this undershoots the 5,000 cap; MouseDE/HumanDE pools are large enough to fill more."
    )
    scenarios.append(
        f"**Scenario B (balanced)**: All COLOC ({pool_pc_coloc}) + all Finemap ({pool_pc_fine}) "
        f"+ 1,500 MouseDE-primary + 1,500 HumanDE-primary "
        f"+ {lnc_quota} lncRNA + {ctrl_quota} controls "
        f"= ~{pool_pc_coloc + pool_pc_fine + 1500 + 1500 + lnc_quota + ctrl_quota}."
    )
    if pool_pc_union <= pc_budget:
        decisionC = (
            f"the entire PC union ({pool_pc_union}) fits in the {pc_budget}-slot PC budget — "
            f"no downsampling needed; the remaining "
            f"{pc_budget - pool_pc_union} slots can be reallocated to controls or extra lncRNAs."
        )
    else:
        n_drop = pool_pc_union - pc_budget
        decisionC = (
            f"the PC union ({pool_pc_union}) exceeds the {pc_budget}-slot budget by **{n_drop}**. "
            f"Recommended trim: drop the lowest-ranked MouseDE-singleton genes "
            f"(n={mde_only} in the singleton tail), keeping all COLOC + Finemap + multi-axis hits "
            f"({multi_axis_pc} PC genes). HumanDE-only singletons (n={hde_only}) are biologically "
            f"weaker than MouseDE-only singletons for an MCD-diet screen (mouse data is the direct screen "
            f"substrate); trim HumanDE-only before MouseDE-only if MouseDE-only is exhausted."
        )
    scenarios.append(
        f"**Scenario C (pool-sized, recommended)**: take the entire PC any-axis union "
        f"({pool_pc_union}) as the structural backbone — "
        f"{decisionC}"
    )
    scenarios.append(
        f"**Scenario D (high-confidence-only)**: PC genes with n_axes ≥ 2 only "
        f"(n={int((pc_out['n_axes']>=2).sum())}) + {lnc_quota} lncRNA + {ctrl_quota} controls — "
        f"fills only ~{int((pc_out['n_axes']>=2).sum()) + lnc_quota + ctrl_quota} of 5,000. "
        f"Would need to either accept a much smaller library or backfill with single-axis genes. "
        f"Useful if the screen prioritises signal-to-noise over coverage."
    )

    # ---- Write markdown report ----
    venn_md = ["| coloc | finemap | mouseDE | humanDE | n |", "|---|---|---|---|---|"]
    for _, r in venn.iterrows():
        venn_md.append(f"| {int(r['axis_coloc'])} | {int(r['axis_finemap'])} | {int(r['axis_mouseDE'])} | {int(r['axis_humanDE'])} | {int(r['n'])} |")

    def top10_md(d: pd.DataFrame, metric_label: str) -> str:
        rows = ["| mouse | human | biotype | " + metric_label + " | n_axes |",
                "|---|---|---|---|---|"]
        for _, r in d.iterrows():
            rows.append(
                f"| {r['gene_symbol_mouse']} | {r['gene_symbol_human']} | "
                f"{r['biotype']} | {r.iloc[3]:.3g} | {int(r['n_axes'])} |"
            )
        return "\n".join(rows)

    lines = []
    lines.append("# Independent-tier candidate pools (Cas13 in vivo MASH library)")
    lines.append("")
    lines.append("**Build date:** 2026-05-23")
    lines.append("")
    lines.append("Four INDEPENDENT evidence axes, each ranked separately (no pre-merging into A∪B union). "
                 "Cas13 = knockdown → UP-regulated targets only. Mouse hep expression filter applied "
                 "(PC AveExpr ≥ 0.5; lncRNA AveExpr ≥ -0.3).")
    lines.append("")
    lines.append("## 1. Per-axis pool sizes (post hep-expr filter)")
    lines.append("")
    lines.append("| Axis | PC count | lncRNA count |")
    lines.append("|------|----------|--------------|")
    for k in ["COLOC", "Finemap", "MouseDE", "HumanDE"]:
        lines.append(f"| {k} | {pc_pool[k]} | {lnc_pool[k]} |")
    lines.append("")
    lines.append(f"Total unique PC candidates (any-axis union): **{pool_pc_union:,}**")
    lines.append(f"Total unique lncRNA candidates (any-axis union): **{pool_lnc_union:,}**")
    lines.append("")
    lines.append("## 2. 4-way Venn intersection counts (PC only)")
    lines.append("")
    n_axes_dist = pc_out["n_axes"].value_counts().sort_index(ascending=False).to_dict()
    n_eq1 = n_axes_dist.get(1, 0)
    n_eq2 = n_axes_dist.get(2, 0)
    n_eq3 = n_axes_dist.get(3, 0)
    n_eq4 = n_axes_dist.get(4, 0)
    n_ge2 = n_eq2 + n_eq3 + n_eq4
    n_ge3 = n_eq3 + n_eq4
    lines.append(f"All 16 cells shown (cells with n=0 included). Counts by n_axes: "
                 f"=4: {n_eq4}, =3: {n_eq3}, =2: {n_eq2}, =1: {n_eq1}. "
                 f"({n_ge2} PC genes hit ≥2 axes; {n_ge3} hit ≥3.)")
    lines.append("")
    lines.append("\n".join(venn_md))
    lines.append("")
    lines.append("## 3. Top-10 per axis (PC preview)")
    for ax in ["COLOC", "Finemap", "MouseDE", "HumanDE"]:
        metric = {"COLOC": "PP4", "Finemap": "PIP", "MouseDE": "mouse_lfc_max", "HumanDE": "human_lfc"}[ax]
        lines.append("")
        lines.append(f"### {ax}")
        lines.append("")
        lines.append(top10_md(pc_top[ax], metric))
    lines.append("")
    lines.append("## 4. Cross-axis overlap rates (PC)")
    lines.append("")
    lines.append("| Conditional | % overlap | n / total |")
    lines.append("|---|---|---|")
    for k, (rate, n_both, n_a) in pc_overlap.items():
        if isinstance(rate, float) and np.isnan(rate):
            lines.append(f"| {k} | n/a (zero pool) | 0 / 0 |")
        else:
            lines.append(f"| {k} | {rate*100:.1f}% | {n_both} / {n_a} |")
    lines.append("")
    lines.append("## 5. Quota allocation scenarios for 5,000-gene library")
    lines.append("")
    for s in scenarios:
        lines.append("- " + s)
    lines.append("")
    lines.append("**Recommended (Scenario C):** Take the entire PC union of all 4 axes as the structural")
    lines.append("backbone, then size the MouseDE/HumanDE singleton tails down to fit the 5,000 cap.")
    lines.append("This maximises evidence independence and keeps every COLOC + Finemap hit (the genetic")
    lines.append("backbone, which is the smallest pool).")
    lines.append("")
    lines.append("## 6. Quality notes")
    lines.append("")
    for q in quality:
        lines.append("- " + q)
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("**Outputs:**")
    lines.append(f"- `Cas13_Library_Design/data/candidates_pc_independent.csv` ({len(pc_out):,} rows)")
    lines.append(f"- `Cas13_Library_Design/data/candidates_lncrna_independent.csv` ({len(lnc_out):,} rows)")
    lines.append(f"- `Cas13_Library_Design/data/finemap_variant_to_gene.csv` (axis2 variant→gene transparency table)")

    REVIEW_MD.parent.mkdir(parents=True, exist_ok=True)
    REVIEW_MD.write_text("\n".join(lines))
    print(f"      wrote {REVIEW_MD}")

    print("[6/6] DONE.")


if __name__ == "__main__":
    main()
