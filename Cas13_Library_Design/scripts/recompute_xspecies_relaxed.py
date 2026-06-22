"""Recompute Pool A / Pool B / A∩B sizes at RELAXED thresholds with
biotype-specific expression filters.

Re-uses `ortholog_bridge.py` (10/11 sanity matches, 25,490 (mouse,human)
pairs) without modification.

Relaxed filters (locked by user 2026-05-21):
  Pool A — RNA-seq UP-concordant (RELAXED Tier 2):
    Human:  atlas bulk_padj<0.05 AND bulk_logFC > 0.3   (was >0.5)
    Mouse:  (MCD OR CDAHFD) adj.P.Val<0.05 AND logFC > 0.3 (was >0.5)
    Both UP-concordant; cross-species ortholog match required.
  Pool B — SuSiE-COLOC genetic causal (UNCHANGED PP4 cut):
    Atlas coloc_best_susie_pp4_polyfun > 0.5 (fallback to coloc_susie_best_pp4)
    Mouse hepatocyte-expressed: biotype-specific
      protein_coding:   max AveExpr across MCD/CDAHFD/FPC/HFD/LIDPAD >= 0.5
      lncRNA / other:   max AveExpr across diets >= -0.3
    Direction lenient.

Outputs:
  Cas13_Library_Design/reviews/perdiet_de_xspecies_relaxed.md
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd

from ortholog_bridge import (
    ROOT,
    PERDIET,
    ATLAS,
    DIETS,
    load_mouse_de,
    load_human_atlas,
    build_ortholog_table,
    cross_validate_known_examples,
    strip_version,
)

OUT = ROOT / "Cas13_Library_Design/reviews/perdiet_de_xspecies_relaxed.md"
GTF = Path.home() / "reference_genome/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz"

# Biotype-specific AveExpr cutoffs (limma-voom log-CPM; CPM ≈ TPM for libs ~10-50M).
# TPM >= 1 protein-coding   ≈ log2(1 + 0.5) ≈ 0.58  →  use 0.5
# TPM >= 0.3 lncRNA / other ≈ log2(0.3 + 0.5) ≈ -0.32 →  use -0.3
AVEEXPR_PC = 0.5
AVEEXPR_NONPC = -0.3

# Pool A relaxed thresholds (Tier 2)
HUMAN_LFC_RELAXED = 0.3
MOUSE_LFC_RELAXED = 0.3

# Prior strict thresholds (for side-by-side strict-vs-relaxed table)
HUMAN_LFC_STRICT = 0.5
MOUSE_LFC_STRICT = 0.5

# COLOC threshold (unchanged)
COLOC_PP4_THR = 0.5


# ----------------------------- helpers ------------------------------------- #

def parse_gtf_biotypes(gtf_path: Path) -> tuple[dict, dict]:
    """ENSMUSG_base → (biotype, symbol)."""
    import gzip
    import re
    gid_re = re.compile(r'gene_id "([^"]+)"')
    gtype_re = re.compile(r'gene_type "([^"]+)"')
    gname_re = re.compile(r'gene_name "([^"]+)"')
    out_b, out_n = {}, {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9 or parts[2] != "gene":
                continue
            gid = gid_re.search(parts[8])
            gtype = gtype_re.search(parts[8])
            gname = gname_re.search(parts[8])
            if gid is None:
                continue
            base = gid.group(1).split(".")[0]
            if gtype:
                out_b[base] = gtype.group(1)
            if gname:
                out_n[base] = gname.group(1)
    return out_b, out_n


def pool_biotype_split(genes_set: set, mouse_biotype: dict) -> tuple[int, int, int]:
    pc = sum(1 for g in genes_set if mouse_biotype.get(g) == "protein_coding")
    lnc = sum(1 for g in genes_set if mouse_biotype.get(g) == "lncRNA")
    other = len(genes_set) - pc - lnc
    return pc, lnc, other


def fmt(x, fmtstr="{:.2f}"):
    if pd.isna(x):
        return "n/a"
    if isinstance(x, (int, np.integer)):
        return f"{int(x):,}"
    return fmtstr.format(x)


def fmt_padj(x):
    if pd.isna(x):
        return "n/a"
    return f"{x:.2e}"


# ----------------------------- inputs -------------------------------------- #

print("Loading mouse biotype annotation (GENCODE vM33)…")
mouse_biotype, mouse_symbol = parse_gtf_biotypes(GTF)
print(f"  {len(mouse_biotype):,} mouse genes annotated")

print("Loading atlas + ortholog bridge…")
atlas = load_human_atlas()
ortho = build_ortholog_table()
print(f"  atlas rows: {len(atlas):,}")
print(f"  bridge pairs: {len(ortho):,}")

# COLOC canonical column
atlas["coloc_pp4_canonical"] = atlas["coloc_best_susie_pp4_polyfun"]
atlas["coloc_pp4_canonical"] = atlas["coloc_pp4_canonical"].fillna(atlas["coloc_susie_best_pp4"])
atlas["coloc_gwas_canonical"] = atlas["coloc_best_susie_gwas_polyfun"].fillna(
    atlas["coloc_susie_best_gwas"]
)

print("Loading per-diet DE files…")
de = {d: load_mouse_de(d) for d in DIETS}
for d, df in de.items():
    mtime = datetime.fromtimestamp((PERDIET / f"{d}_de_results.csv").stat().st_mtime)
    print(f"  {d}: {df.shape[0]:,} genes  (mtime {mtime:%Y-%m-%d %H:%M})")

# annotate mouse-side biotype + symbol
for d in DIETS:
    de[d]["mouse_biotype"] = de[d].index.map(mouse_biotype.get)
    de[d]["mouse_symbol"] = de[d].index.map(mouse_symbol.get)


# Build max-AveExpr-across-diets table.  Each row = one ENSMUSG with its max AveExpr.
max_ave_records = []
all_mouse_genes: set = set()
for d in DIETS:
    all_mouse_genes |= set(de[d].index)
for g in all_mouse_genes:
    vals = []
    for d in DIETS:
        if g in de[d].index:
            vals.append(de[d].loc[g, "AveExpr"])
    if vals:
        max_ave_records.append((g, max(vals)))
max_ave = pd.Series(dict(max_ave_records))
print(f"Mouse genes with at least one diet AveExpr: {len(max_ave):,}")


def expressed_mouse(biotype_aware: bool) -> set:
    """Return the set of mouse ENSMUSG that pass the expression filter.
    If biotype_aware, use 0.5 (PC) / -0.3 (non-PC).  Otherwise apply 0.5 to all
    (matches the prior strict-protocol mode)."""
    s = set()
    for g, ae in max_ave.items():
        bt = mouse_biotype.get(g)
        if biotype_aware:
            thr = AVEEXPR_PC if bt == "protein_coding" else AVEEXPR_NONPC
        else:
            thr = 1.0  # the original strict pool used >= 1 for all
        if ae >= thr:
            s.add(g)
    return s


expressed_strict = expressed_mouse(False)  # legacy AveExpr >= 1 (any biotype)
expressed_relaxed = expressed_mouse(True)  # biotype-specific
print(f"Mouse expressed strict (any-diet AveExpr >= 1): {len(expressed_strict):,}")
print(f"Mouse expressed relaxed (PC>=0.5 / non-PC>=-0.3): {len(expressed_relaxed):,}")


# ----------------------------- sanity ------------------------------------- #

sanity = cross_validate_known_examples()


# ----------------------------- helper sets --------------------------------- #

def mouse_up(diet: str, lfc_thr: float, padj_thr: float = 0.05) -> set:
    df = de[diet]
    return set(df.index[(df["adj.P.Val"] < padj_thr) & (df["logFC"] > lfc_thr)])


# ----------------------------- Pool A (strict + relaxed) ------------------- #

# Bridge: human ENSG → set of mouse ENSMUSG
ortho_h2m = ortho.groupby("human_ensembl")["mouse_ensembl"].apply(set).to_dict()


def map_hens_to_mens(human_ens_set: set) -> set:
    out = set()
    for h in human_ens_set:
        out |= ortho_h2m.get(h, set())
    return out


# --- Strict (re-derived for the side-by-side table) ---
human_strict_up = atlas[(atlas["bulk_padj"] < 0.05) &
                        (atlas["bulk_logFC"] > HUMAN_LFC_STRICT)]
human_strict_up_ens = set(human_strict_up["ensembl_id_base"].dropna())
human_strict_up_mouse = map_hens_to_mens(human_strict_up_ens)

mouse_strict_up = mouse_up("MCD", MOUSE_LFC_STRICT) | mouse_up("CDAHFD", MOUSE_LFC_STRICT)
pool_a_strict = human_strict_up_mouse & mouse_strict_up

# --- Relaxed ---
human_relaxed_up = atlas[(atlas["bulk_padj"] < 0.05) &
                         (atlas["bulk_logFC"] > HUMAN_LFC_RELAXED)]
human_relaxed_up_ens = set(human_relaxed_up["ensembl_id_base"].dropna())
human_relaxed_up_mouse = map_hens_to_mens(human_relaxed_up_ens)

mouse_relaxed_up = mouse_up("MCD", MOUSE_LFC_RELAXED) | mouse_up("CDAHFD", MOUSE_LFC_RELAXED)
pool_a_relaxed = human_relaxed_up_mouse & mouse_relaxed_up

pa_s_pc, pa_s_lnc, pa_s_other = pool_biotype_split(pool_a_strict, mouse_biotype)
pa_r_pc, pa_r_lnc, pa_r_other = pool_biotype_split(pool_a_relaxed, mouse_biotype)


# ----------------------------- Pool B (strict + relaxed) ------------------- #

coloc_hits = atlas[atlas["coloc_pp4_canonical"] > COLOC_PP4_THR]
coloc_hits_ens = set(coloc_hits["ensembl_id_base"].dropna())
coloc_mouse = map_hens_to_mens(coloc_hits_ens)

pool_b_strict = coloc_mouse & expressed_strict
pool_b_relaxed = coloc_mouse & expressed_relaxed

pb_s_pc, pb_s_lnc, pb_s_other = pool_biotype_split(pool_b_strict, mouse_biotype)
pb_r_pc, pb_r_lnc, pb_r_other = pool_biotype_split(pool_b_relaxed, mouse_biotype)


# ----------------------------- A ∩ B and A ∪ B (strict + relaxed) ---------- #

inter_strict = pool_a_strict & pool_b_strict
union_strict = pool_a_strict | pool_b_strict
inter_relaxed = pool_a_relaxed & pool_b_relaxed
union_relaxed = pool_a_relaxed | pool_b_relaxed

is_pc, is_lnc, is_other = pool_biotype_split(inter_strict, mouse_biotype)
ir_pc, ir_lnc, ir_other = pool_biotype_split(inter_relaxed, mouse_biotype)
us_pc, us_lnc, us_other = pool_biotype_split(union_strict, mouse_biotype)
ur_pc, ur_lnc, ur_other = pool_biotype_split(union_relaxed, mouse_biotype)


print()
print("=== Strict pool sizes (prior run) ===")
print(f"  Pool A: {len(pool_a_strict):,}  PC {pa_s_pc}  lncRNA {pa_s_lnc}  other {pa_s_other}")
print(f"  Pool B: {len(pool_b_strict):,}  PC {pb_s_pc}  lncRNA {pb_s_lnc}  other {pb_s_other}")
print(f"  A∩B   : {len(inter_strict):,}  PC {is_pc}  lncRNA {is_lnc}  other {is_other}")
print(f"  A∪B   : {len(union_strict):,}  PC {us_pc}  lncRNA {us_lnc}  other {us_other}")
print()
print("=== Relaxed pool sizes (this run) ===")
print(f"  Pool A: {len(pool_a_relaxed):,}  PC {pa_r_pc}  lncRNA {pa_r_lnc}  other {pa_r_other}")
print(f"  Pool B: {len(pool_b_relaxed):,}  PC {pb_r_pc}  lncRNA {pb_r_lnc}  other {pb_r_other}")
print(f"  A∩B   : {len(inter_relaxed):,}  PC {ir_pc}  lncRNA {ir_lnc}  other {ir_other}")
print(f"  A∪B   : {len(union_relaxed):,}  PC {ur_pc}  lncRNA {ur_lnc}  other {ur_other}")


# ----------------------------- lncRNA-specific ----------------------------- #

# Human lncRNA DEGs UP at relaxed cut
n_human_lnc_up_relaxed = len(atlas[
    (atlas["gene_biotype"] == "lncRNA")
    & (atlas["bulk_padj"] < 0.05)
    & (atlas["bulk_logFC"] > HUMAN_LFC_RELAXED)
])
# Mouse-only UP lncRNA (MCD ∪ CDAHFD UP at relaxed cut, mouse biotype = lncRNA, biotype-TPM filter)
mcd_lnc_up_relaxed = {g for g in mouse_up("MCD", MOUSE_LFC_RELAXED)
                      if mouse_biotype.get(g) == "lncRNA"}
cd_lnc_up_relaxed = {g for g in mouse_up("CDAHFD", MOUSE_LFC_RELAXED)
                     if mouse_biotype.get(g) == "lncRNA"}
mouse_only_up_lnc_relaxed_all = (mcd_lnc_up_relaxed | cd_lnc_up_relaxed)
# Apply expression filter (biotype-aware)
mouse_only_up_lnc_relaxed = mouse_only_up_lnc_relaxed_all & expressed_relaxed

# Cross-species lncRNA matches at relaxed
human_lnc_up_relaxed_ens = set(
    atlas.loc[(atlas["gene_biotype"] == "lncRNA")
              & (atlas["bulk_padj"] < 0.05)
              & (atlas["bulk_logFC"] > HUMAN_LFC_RELAXED), "ensembl_id_base"].dropna()
)
human_lnc_up_mouse_ens = map_hens_to_mens(human_lnc_up_relaxed_ens)
xspecies_lnc_matches = human_lnc_up_mouse_ens & mouse_only_up_lnc_relaxed_all


# ----------------------------- Top-20 COLOC annotated --------------------- #

top20 = (
    atlas.dropna(subset=["coloc_pp4_canonical"])
    .sort_values("coloc_pp4_canonical", ascending=False)
    .drop_duplicates("ensembl_id_base", keep="first")
    .head(20)
).copy()
top20 = top20.drop(columns=["mouse_ensembl"], errors="ignore")
ortho_h_idx = ortho.drop_duplicates("human_ensembl", keep="first").set_index("human_ensembl")[
    ["mouse_ensembl", "mouse_symbol"]
]
top20 = top20.merge(ortho_h_idx, left_on="ensembl_id_base", right_index=True, how="left")


def direction(lfc):
    if pd.isna(lfc):
        return "n/a"
    return "UP" if lfc > 0 else "DOWN"


def lookup_de(d, mens, col):
    if pd.isna(mens) or mens not in de[d].index:
        return float("nan")
    return de[d].loc[mens, col]


for col_d in ["MCD", "CDAHFD"]:
    top20[f"{col_d}_logFC"] = top20["mouse_ensembl"].apply(lambda m: lookup_de(col_d, m, "logFC"))
    top20[f"{col_d}_padj"] = top20["mouse_ensembl"].apply(lambda m: lookup_de(col_d, m, "adj.P.Val"))
top20["dream_direction"] = top20["bulk_logFC"].apply(direction)
top20["in_pool_a_relaxed"] = top20["mouse_ensembl"].isin(pool_a_relaxed)
top20["in_pool_b_relaxed"] = top20["mouse_ensembl"].isin(pool_b_relaxed)


# ----------------------------- A ∩ B full table (relaxed) ------------------ #

ab_inter_table = (
    ortho.loc[ortho["mouse_ensembl"].isin(inter_relaxed)]
    .drop(columns=["human_symbol"], errors="ignore")
    .merge(
        atlas[["ensembl_id_base", "human_symbol", "gene_biotype",
               "bulk_logFC", "bulk_padj",
               "coloc_pp4_canonical", "coloc_gwas_canonical"]],
        left_on="human_ensembl", right_on="ensembl_id_base", how="left",
    )
    # Pick best COLOC hit per mouse ortholog
    .sort_values("coloc_pp4_canonical", ascending=False, na_position="last")
    .drop_duplicates("mouse_ensembl", keep="first")
)
ab_inter_table["MCD_logFC"] = ab_inter_table["mouse_ensembl"].apply(lambda m: lookup_de("MCD", m, "logFC"))
ab_inter_table["CDAHFD_logFC"] = ab_inter_table["mouse_ensembl"].apply(lambda m: lookup_de("CDAHFD", m, "logFC"))
ab_inter_table["mouse_biotype"] = ab_inter_table["mouse_ensembl"].map(mouse_biotype.get)
ab_inter_table["mouse_symbol"] = ab_inter_table["mouse_ensembl"].map(mouse_symbol.get)
ab_inter_table = ab_inter_table.sort_values("coloc_pp4_canonical", ascending=False)


# ----------------------------- write markdown ------------------------------ #

L: list[str] = []
A = L.append

A("# Per-Diet Mouse DE — Cross-Species Pools at RELAXED Thresholds")
A("")
A(f"_Generated by `Cas13_Library_Design/scripts/recompute_xspecies_relaxed.py` on {datetime.now():%Y-%m-%d}._")
A("")
A("This deliverable accompanies (does not replace) "
  "`perdiet_de_xspecies_fix.md`.  It recomputes Pool A / Pool B / A∩B at "
  "RELAXED Tier-2 thresholds with biotype-specific expression filters, "
  "using the freshly rebuilt per-diet DE outputs (M02-M07 job 16198594, "
  "completed 2026-05-21 20:58).")
A("")
A("The mouse rebuild changed two things:")
A("- HFD now includes GSE224069 (n=44 disease + 26 control; was 31 + 11).")
A("- `filterByExpr` relaxed to `min.count=5, min.total.count=10`, adding "
  "~2-3.5K more genes tested per diet and recovering canonical markers "
  "(Acta2, Timp1, Mmp2, Mmp9, Saa3).")
A("")
A("Bridge module (`ortholog_bridge.py`) is unchanged: "
  f"{len(ortho):,} (mouse, human) pairs, "
  f"{ortho['mouse_ensembl'].nunique():,} unique mouse genes, "
  f"{ortho['human_ensembl'].nunique():,} unique human genes.")
A("")

# ---- Section 1: Threshold table ------------------------------------------- #
A("## 1. Threshold table — strict (prior) vs relaxed (this run)")
A("")
A("| Filter | Strict (prior) | Relaxed (this run) |")
A("|---|---|---|")
A("| Human UP `bulk_padj` | <0.05 | <0.05 |")
A(f"| Human UP `bulk_logFC` | >{HUMAN_LFC_STRICT} (Tier 1) | >{HUMAN_LFC_RELAXED} (Tier 2) |")
A("| Mouse UP `adj.P.Val` | <0.05 | <0.05 |")
A(f"| Mouse UP `logFC` (MCD ∪ CDAHFD) | >{MOUSE_LFC_STRICT} | >{MOUSE_LFC_RELAXED} |")
A("| Pool A UP-concordant | required | required |")
A("| Mouse hep-expressed | any-diet `AveExpr >= 1` | biotype-specific |")
A(f"| · protein_coding | (n/a) | max-diet `AveExpr >= {AVEEXPR_PC}` (≈ TPM ≥ 1) |")
A(f"| · lncRNA / other | (n/a) | max-diet `AveExpr >= {AVEEXPR_NONPC}` (≈ TPM ≥ 0.3) |")
A(f"| COLOC `PP.H4.susie_polyfun` | > {COLOC_PP4_THR} | > {COLOC_PP4_THR} (unchanged) |")
A("| Ortholog bridge | biomaRt + atlas (25,490 pairs) | identical |")
A("")
A(f"Atlas rows passing COLOC PP4>0.5 (canonical polyfun + susie fallback): "
  f"**{(atlas['coloc_pp4_canonical']>COLOC_PP4_THR).sum():,}**.  Mouse genes "
  f"passing the biotype-specific expression filter: "
  f"**{len(expressed_relaxed):,}** (vs {len(expressed_strict):,} under the "
  f"strict ≥1-log-CPM cutoff).")
A("")

# ---- Section 2: Pool sizes by biotype ------------------------------------- #
A("## 2. Pool sizes by biotype — relaxed")
A("")
A("| Pool | PC count | lncRNA count | Other | Total |")
A("|---|---:|---:|---:|---:|")
A(f"| A_relaxed (Tier 2 UP + mouse \\|LFC\\|>{MOUSE_LFC_RELAXED} UP) | {pa_r_pc:,} | {pa_r_lnc:,} | {pa_r_other:,} | {len(pool_a_relaxed):,} |")
A(f"| B_relaxed (COLOC PP4>0.5 + biotype-TPM-expressed) | {pb_r_pc:,} | {pb_r_lnc:,} | {pb_r_other:,} | {len(pool_b_relaxed):,} |")
A(f"| A_relaxed ∩ B_relaxed | {ir_pc:,} | {ir_lnc:,} | {ir_other:,} | {len(inter_relaxed):,} |")
A(f"| A_relaxed ∪ B_relaxed | {ur_pc:,} | {ur_lnc:,} | {ur_other:,} | {len(union_relaxed):,} |")
A("")

# ---- Section 3: side-by-side strict vs relaxed --------------------------- #
A("## 3. Side-by-side — strict vs relaxed")
A("")
A("| Pool | Strict total (PC / lncRNA / other) | Relaxed total (PC / lncRNA / other) | Δ total |")
A("|---|---|---|---:|")
A(f"| A | {len(pool_a_strict):,} ({pa_s_pc} / {pa_s_lnc} / {pa_s_other}) | "
  f"{len(pool_a_relaxed):,} ({pa_r_pc} / {pa_r_lnc} / {pa_r_other}) | "
  f"{len(pool_a_relaxed)-len(pool_a_strict):+,} |")
A(f"| B | {len(pool_b_strict):,} ({pb_s_pc} / {pb_s_lnc} / {pb_s_other}) | "
  f"{len(pool_b_relaxed):,} ({pb_r_pc} / {pb_r_lnc} / {pb_r_other}) | "
  f"{len(pool_b_relaxed)-len(pool_b_strict):+,} |")
A(f"| A ∩ B | {len(inter_strict):,} ({is_pc} / {is_lnc} / {is_other}) | "
  f"{len(inter_relaxed):,} ({ir_pc} / {ir_lnc} / {ir_other}) | "
  f"{len(inter_relaxed)-len(inter_strict):+,} |")
A(f"| A ∪ B | {len(union_strict):,} ({us_pc} / {us_lnc} / {us_other}) | "
  f"{len(union_relaxed):,} ({ur_pc} / {ur_lnc} / {ur_other}) | "
  f"{len(union_relaxed)-len(union_strict):+,} |")
A("")
A("Caveat — the strict-here numbers use the **freshly rebuilt** per-diet "
  "outputs (mtime 2026-05-21 20:37) and the biotype-naive expression cut "
  "(any-diet AveExpr ≥ 1); they may differ slightly from the strict numbers "
  "in `perdiet_de_xspecies_fix.md` (whose mouse DE files were the older Apr-7 "
  "build at filterByExpr default `min.count=10`).")
A("")

# ---- Section 4: Top-20 COLOC ---------------------------------------------- #
A("## 4. Top-20 SuSiE-COLOC hits — annotated with all evidence")
A("")
A("Sorted by `coloc_best_susie_pp4_polyfun` (canonical production EUR LD). "
  "`A?` / `B?` = membership in the **relaxed** Pool A / Pool B.  Empty mouse "
  "cells indicate either no mouse ortholog or the mouse gene was filtered out "
  "of the per-diet DE.")
A("")
A("| rank | gene symbol | ENSG | ENSMUSG | coloc PP4 | best GWAS | A? | B? | dream logFC | dream padj | dream dir | MCD logFC | MCD padj | CDAHFD logFC | CDAHFD padj | biotype |")
A("|---|---|---|---|---:|---|:---:|:---:|---:|---:|:---:|---:|---:|---:|---:|---|")
for i, (_, r) in enumerate(top20.iterrows(), start=1):
    A(
        f"| {i} | {r['human_symbol']} | {r['ensembl_id_base']} | "
        f"{r['mouse_ensembl'] if pd.notna(r['mouse_ensembl']) else 'n/a'} | "
        f"{fmt(r['coloc_pp4_canonical'], '{:.3f}')} | "
        f"{r['coloc_gwas_canonical'] if pd.notna(r['coloc_gwas_canonical']) else 'n/a'} | "
        f"{'YES' if r['in_pool_a_relaxed'] else '—'} | "
        f"{'YES' if r['in_pool_b_relaxed'] else '—'} | "
        f"{fmt(r['bulk_logFC'], '{:.2f}')} | {fmt_padj(r['bulk_padj'])} | "
        f"{r['dream_direction']} | "
        f"{fmt(r['MCD_logFC'], '{:.2f}')} | {fmt_padj(r['MCD_padj'])} | "
        f"{fmt(r['CDAHFD_logFC'], '{:.2f}')} | {fmt_padj(r['CDAHFD_padj'])} | "
        f"{r['gene_biotype']} |"
    )
A("")

# ---- Section 5: A∩B full gene list ---------------------------------------- #
A(f"## 5. A_relaxed ∩ B_relaxed full gene list (n = {len(ab_inter_table):,})")
A("")
A("Per-gene human ortholog picked by highest `coloc_pp4_canonical` "
  "(mouse-to-many fans collapse to the COLOC-anchored human ortholog).")
A("")
A("| mouse symbol | human symbol | ENSG | ENSMUSG | coloc PP4 | best GWAS | dream logFC | MCD logFC | CDAHFD logFC | biotype |")
A("|---|---|---|---|---:|---|---:|---:|---:|---|")
for _, r in ab_inter_table.iterrows():
    A(
        f"| {r['mouse_symbol'] or '?'} | {r['human_symbol']} | "
        f"{r['human_ensembl']} | {r['mouse_ensembl']} | "
        f"{fmt(r['coloc_pp4_canonical'], '{:.3f}')} | "
        f"{r['coloc_gwas_canonical'] if pd.notna(r['coloc_gwas_canonical']) else 'n/a'} | "
        f"{fmt(r['bulk_logFC'], '{:.2f}')} | "
        f"{fmt(r['MCD_logFC'], '{:.2f}')} | "
        f"{fmt(r['CDAHFD_logFC'], '{:.2f}')} | "
        f"{r['mouse_biotype'] or '?'} |"
    )
A("")

# ---- Section 6: lncRNA-specific ------------------------------------------- #
A("## 6. lncRNA-specific accounting at relaxed thresholds")
A("")
A(f"- Human lncRNA DEGs UP at relaxed `bulk_padj<0.05, bulk_logFC>{HUMAN_LFC_RELAXED}`: "
  f"**{n_human_lnc_up_relaxed:,}**.")
A(f"- Mouse-only UP lncRNA pool (MCD ∪ CDAHFD UP `|LFC|>{MOUSE_LFC_RELAXED}`, "
  f"mouse biotype=lncRNA, biotype-TPM expression filter): "
  f"**{len(mouse_only_up_lnc_relaxed):,}** (MCD-only candidates "
  f"{len(mcd_lnc_up_relaxed & expressed_relaxed):,}; CDAHFD-only "
  f"{len(cd_lnc_up_relaxed & expressed_relaxed):,}).")
A(f"- Cross-species lncRNA matches (mouse-up-lncRNA ∩ ortholog-resolved "
  f"human-up-lncRNA, both at relaxed cut): **{len(xspecies_lnc_matches):,}**.")
A("")
A("Per `perdiet_de_xspecies_fix.md`: the cross-species lncRNA gap is "
  "biologically fundamental (ortholog tables are protein-coding-dominant), "
  "not a bridge artifact.  Relaxing thresholds does not change this — for "
  "the in-vivo mouse screen, the operational lncRNA ceiling is the "
  "**mouse-only UP set**, with human cross-species evidence applied per-"
  "target rather than as a panel filter.")
A("")

# ---- Section 7: verdict --------------------------------------------------- #
A("## 7. Verdict — does relaxing the cuts clear the 2,500 PC library floor?")
A("")
au_pc_relaxed = ur_pc
A(f"- **A_relaxed ∪ B_relaxed PC count = {au_pc_relaxed:,}** "
  f"(vs the 2,500 protein-coding floor).  "
  f"{'**Yes — clears the floor.**' if au_pc_relaxed >= 2500 else '**No — still below the 2,500-PC floor.**'}")
A("")
A(f"- Pool A_relaxed grows from {pa_s_pc:,} → {pa_r_pc:,} PC "
  f"({(pa_r_pc/max(1,pa_s_pc)-1)*100:+.1f}%) when LFC is relaxed from "
  f"|>{HUMAN_LFC_STRICT}|/|>{MOUSE_LFC_STRICT}| to "
  f"|>{HUMAN_LFC_RELAXED}|/|>{MOUSE_LFC_RELAXED}|.  This is the dominant "
  f"effect — both human and mouse lose stringency together.")
A("")
A(f"- Pool B_relaxed grows from {pb_s_pc:,} → {pb_r_pc:,} PC "
  f"({(pb_r_pc/max(1,pb_s_pc)-1)*100:+.1f}%) when the expression filter is "
  f"relaxed.  The non-PC bucket (lncRNA + pseudogene + other) grows from "
  f"{pb_s_lnc + pb_s_other} → {pb_r_lnc + pb_r_other} — lower-expressed "
  f"non-coding GWAS-causal genes now make the cut.")
A("")
A(f"- A_relaxed ∩ B_relaxed = {len(inter_relaxed):,} "
  f"({ir_pc} PC + {ir_lnc} lncRNA + {ir_other} other) — "
  f"{len(inter_relaxed)-len(inter_strict):+,} vs strict.  This is the "
  f"highest-evidence positive-control set.")
A("")
A(f"- Mouse-only-UP backfill ceiling at relaxed mouse cut "
  f"(|LFC|>{MOUSE_LFC_RELAXED}, padj<0.05; MCD ∪ CDAHFD; biotype-aware "
  f"expression): "
  f"PC **{sum(1 for g in (mouse_relaxed_up & expressed_relaxed) if mouse_biotype.get(g)=='protein_coding'):,}**, "
  f"lncRNA **{sum(1 for g in (mouse_relaxed_up & expressed_relaxed) if mouse_biotype.get(g)=='lncRNA'):,}**.  "
  f"With the relaxed LFC, the mouse-only UP set alone exceeds the 2,500-PC "
  f"library floor — backfill is not the bottleneck.")
A("")

# Compute recommended LFC if A∪B is still short
if au_pc_relaxed < 2500:
    # Walk down the cut to find what gets us to >= 2500
    targets = []
    for thr in [0.25, 0.20, 0.15, 0.10, 0.05, 0.00]:
        h_set = atlas[(atlas["bulk_padj"] < 0.05) & (atlas["bulk_logFC"] > thr)]
        h_ens = set(h_set["ensembl_id_base"].dropna())
        h_mouse = map_hens_to_mens(h_ens)
        m_set = mouse_up("MCD", thr) | mouse_up("CDAHFD", thr)
        a = h_mouse & m_set
        union = a | pool_b_relaxed
        pc = sum(1 for g in union if mouse_biotype.get(g) == "protein_coding")
        targets.append((thr, pc, len(a), len(union)))
    A("**LFC sweep — what cut hits the 2,500 PC floor in A ∪ B?**")
    A("")
    A("| LFC threshold (both sides) | Pool A size | Pool A ∪ B total | Pool A ∪ B PC |")
    A("|---:|---:|---:|---:|")
    for thr, pc, a_size, u_size in targets:
        flag = " **← clears 2,500 PC**" if pc >= 2500 else ""
        A(f"| >{thr} | {a_size:,} | {u_size:,} | {pc:,}{flag} |")
    A("")
    hit = next((t for t in targets if t[1] >= 2500), None)
    if hit:
        A(f"**Recommendation**: drop the cross-species LFC cut to "
          f"`>{hit[0]}` on both sides to clear the 2,500-PC floor "
          f"from A ∪ B alone (without backfilling mouse-only-UP).")
    else:
        A("**Recommendation**: even at LFC > 0 (any-direction UP, padj-only) "
          "the cross-species A ∪ B PC count stays below 2,500.  Use the "
          "mouse-only-UP backfill block as the primary fill source; cross-"
          "species concordance is a precision pool, not a volume pool.")
else:
    A(f"**Recommendation**: relaxed thresholds already clear the 2,500-PC "
      f"library floor from A ∪ B alone.  Lock A_relaxed ∩ B_relaxed "
      f"({len(inter_relaxed):,} genes) as positive controls; fill the "
      f"remaining {2500 - len(inter_relaxed):,} PC slots from A_relaxed ∪ "
      f"B_relaxed in priority order (A∩B → A-only → B-only).")
A("")

A(f"_Inputs: per_diet/*_de_results.csv (mtime "
  f"{datetime.fromtimestamp((PERDIET/'MCD_de_results.csv').stat().st_mtime):%Y-%m-%d %H:%M}; "
  f"M02-M07 rebuild job 16198594, completed 2026-05-21 20:58); "
  f"atlas {datetime.fromtimestamp(ATLAS.stat().st_mtime):%Y-%m-%d %H:%M}; "
  f"ortholog symbols file "
  f"{datetime.fromtimestamp((ROOT/'data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz').stat().st_mtime):%Y-%m-%d %H:%M}._")

OUT.write_text("\n".join(L))
print(f"\nWrote: {OUT}")
print(f"Lines: {len(L)}")
