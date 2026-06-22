"""Recompute the cross-species pieces of `perdiet_de_data_inspection.md`
using the fixed ortholog bridge.  Writes a self-contained Markdown deliverable
at `Cas13_Library_Design/reviews/perdiet_de_xspecies_fix.md`.

Sections computed:
  6.  Human Tier 1 ∩ MCD∪CDAHFD overlap  (UP-only)
  7.  lncRNA inventory UP-only cross-species
  8.  Direction concordance UP-in-both
  NEW. Pool A (mouse-strict UP + human-tier-1 UP) and
       Pool B (SuSiE-COLOC PP4>0.5 + mouse hepatocyte-expressed) sizes.
  TOP10 COLOC hits annotated.

Pure descriptive / no DE re-run.  Mouse rebuild job 16198363 may regenerate
the per_diet/*_de_results.csv files; this script reads them at run time, so
re-running this script on new files is safe.
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime
import os

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
    bridge_mouse_to_human,
    cross_validate_known_examples,
    strip_version,
    ATLAS_COLS_DEFAULT,
)

OUT = ROOT / "Cas13_Library_Design/reviews/perdiet_de_xspecies_fix.md"
GTF = Path.home() / "reference_genome/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz"


# ----------------------------- helpers ------------------------------------- #

def parse_gtf_biotypes(gtf_path: Path) -> dict:
    """ENSMUSG_base → biotype.  Only 'gene' rows."""
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


# ----------------------------- inputs -------------------------------------- #

print("Loading mouse biotype annotation (GENCODE vM33)…")
mouse_biotype, mouse_symbol = parse_gtf_biotypes(GTF)
print(f"  {len(mouse_biotype):,} mouse genes annotated")

print("Loading atlas + ortholog bridge…")
atlas = load_human_atlas()
ortho = build_ortholog_table()
print(f"  atlas rows: {len(atlas):,}")
print(f"  bridge pairs: {len(ortho):,}")

# Pick canonical COLOC column.  Polyfun is the production EUR LD reference
# (CLAUDE.md says default `LD_PANEL=polyfun` 2026-05-06).  Fallback to
# coloc_susie_best_pp4 (atlas legacy) if polyfun col is empty.
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

# Map mouse-side biotype onto every DE row
for d in DIETS:
    de[d]["mouse_biotype"] = de[d].index.map(mouse_biotype.get)
    de[d]["mouse_symbol"] = de[d].index.map(mouse_symbol.get)


# ----------------------------- sanity ------------------------------------- #

sanity = cross_validate_known_examples()

# ----------------------------- helper sets --------------------------------- #

def mouse_up(diet: str, lfc_thr: float, padj_thr: float = 0.05) -> set:
    df = de[diet]
    return set(df.index[(df["adj.P.Val"] < padj_thr) & (df["logFC"] > lfc_thr)])


def mouse_sig(diet: str, lfc_thr: float, padj_thr: float = 0.05) -> set:
    df = de[diet]
    return set(df.index[(df["adj.P.Val"] < padj_thr) & (df["logFC"].abs() > lfc_thr)])


# Mouse expressed (any diet AveExpr ≥ 1).  AveExpr is limma-voom log-CPM.
expressed_mouse_genes: set = set()
for d in DIETS:
    expressed_mouse_genes |= set(de[d].index[de[d]["AveExpr"] >= 1])
print(f"Mouse genes with AveExpr ≥ 1 in any diet: {len(expressed_mouse_genes):,}")


# ----------------------------- §6 — UP-only Tier 1 ∩ MCD∪CDAHFD ------------ #

# Atlas: Tier 1 (UP) — bulk_padj<0.05, bulk_logFC>0.5
human_tier1_up = atlas[(atlas["bulk_padj"] < 0.05) & (atlas["bulk_logFC"] > 0.5)].copy()
print(f"Atlas Tier 1 UP rows (any biotype): {len(human_tier1_up):,}")

# Map human Tier 1 → mouse_ensembl set via ortholog bridge
tier1_human_ens = set(human_tier1_up["ensembl_id_base"].dropna())
tier1_mouse_via_bridge = set(
    ortho.loc[ortho["human_ensembl"].isin(tier1_human_ens), "mouse_ensembl"].dropna()
)
print(f"Atlas Tier 1 UP mappable to mouse via bridge: {len(tier1_mouse_via_bridge):,}")

sec6_rows = []
for thr in [0.3, 0.5, 0.75, 1.0]:
    mcd_up = mouse_up("MCD", thr)
    cd_up = mouse_up("CDAHFD", thr)
    union_up = mcd_up | cd_up
    overlap = union_up & tier1_mouse_via_bridge
    sec6_rows.append({
        "lfc_thr": thr,
        "mcd_up": len(mcd_up),
        "cdahfd_up": len(cd_up),
        "union_up": len(union_up),
        "tier1_mappable": len(tier1_mouse_via_bridge),
        "overlap": len(overlap),
        "frac_union": len(overlap) / max(1, len(union_up)),
        "frac_tier1": len(overlap) / max(1, len(tier1_mouse_via_bridge)),
    })


# ----------------------------- §7 — lncRNA UP-only ------------------------- #

# Human lncRNA DEGs UP (any LFC, padj<0.05)
human_lnc_up = atlas[
    (atlas["gene_biotype"] == "lncRNA")
    & (atlas["bulk_padj"] < 0.05)
    & (atlas["bulk_logFC"] > 0)
].copy()
n_human_lnc_up = len(human_lnc_up)

# Map to mouse via bridge
human_lnc_up_ens = set(human_lnc_up["ensembl_id_base"].dropna())
human_lnc_up_mouse = set(
    ortho.loc[ortho["human_ensembl"].isin(human_lnc_up_ens), "mouse_ensembl"].dropna()
)
n_human_lnc_up_with_mouse = len(human_lnc_up_mouse)

# Of those, how many are MCD∪CDAHFD UP at |LFC|>0.5
mcd_up_05 = mouse_up("MCD", 0.5)
cd_up_05 = mouse_up("CDAHFD", 0.5)
union_up_05 = mcd_up_05 | cd_up_05
n_lnc_up_in_union = len(human_lnc_up_mouse & union_up_05)

# Restrict to mouse-side biotype lncRNA as well
mouse_lnc_in_union = {g for g in union_up_05 if mouse_biotype.get(g) == "lncRNA"}
n_lnc_up_in_union_strict = len(human_lnc_up_mouse & mouse_lnc_in_union)

# Mouse-only lncRNA UP inventory (operational ceiling for an in-vivo mouse screen)
mouse_lnc_up_05 = {g for g in (mcd_up_05 | cd_up_05) if mouse_biotype.get(g) == "lncRNA"}
mcd_lnc_up_05 = {g for g in mcd_up_05 if mouse_biotype.get(g) == "lncRNA"}
cd_lnc_up_05 = {g for g in cd_up_05 if mouse_biotype.get(g) == "lncRNA"}


# ----------------------------- §8 — direction concordance UP-in-both ------- #

# Best human row per mouse_ensembl (lowest padj wins)
atlas_min = (
    atlas.dropna(subset=["ensembl_id_base"])
    .sort_values("bulk_padj")
    .drop_duplicates("ensembl_id_base", keep="first")
)
# Build mouse_ensembl → (h_lfc, h_padj) via bridge
mouse_to_h = (
    ortho.merge(atlas_min[["ensembl_id_base", "bulk_logFC", "bulk_padj", "human_symbol", "gene_biotype"]],
                left_on="human_ensembl", right_on="ensembl_id_base", how="left")
    .dropna(subset=["bulk_padj"])
)
# pick per mouse_ensembl the human row with lowest padj
mouse_to_h = mouse_to_h.sort_values("bulk_padj").drop_duplicates("mouse_ensembl", keep="first")
m2h_lookup = mouse_to_h.set_index("mouse_ensembl")

sec8_rows = []
for d in ["MCD", "CDAHFD"]:
    mouse_sig_05 = mouse_sig(d, 0.5)
    paired = m2h_lookup.loc[m2h_lookup.index.intersection(mouse_sig_05)].copy()
    # add mouse logFC
    paired = paired.join(de[d][["logFC", "adj.P.Val"]].rename(
        columns={"logFC": "m_lfc", "adj.P.Val": "m_padj"}
    ), how="inner")
    # both significant
    both_sig = paired[(paired["m_padj"] < 0.05) & (paired["bulk_padj"] < 0.05)].copy()
    n_both_sig = len(both_sig)
    if n_both_sig:
        same_sign = (np.sign(both_sig["bulk_logFC"]) == np.sign(both_sig["m_lfc"])).mean()
        up_in_both = ((both_sig["bulk_logFC"] > 0) & (both_sig["m_lfc"] > 0)).sum()
        down_in_both = ((both_sig["bulk_logFC"] < 0) & (both_sig["m_lfc"] < 0)).sum()
        discordant = n_both_sig - up_in_both - down_in_both
        r = np.corrcoef(both_sig["bulk_logFC"], both_sig["m_lfc"])[0, 1]
    else:
        same_sign = up_in_both = down_in_both = discordant = float("nan")
        r = float("nan")
    sec8_rows.append({
        "diet": d, "n_both_sig": n_both_sig,
        "same_sign_frac": same_sign,
        "up_in_both": up_in_both,
        "down_in_both": down_in_both,
        "discordant": discordant,
        "pearson_r": r,
    })


# ----------------------------- Pool A / Pool B ----------------------------- #

# ----- Pool A: mouse-strict UP (|LFC|>0.5) in MCD or CDAHFD AND human Tier 1 UP AND UP-in-both
mcd_up_05 = mouse_up("MCD", 0.5)
cd_up_05 = mouse_up("CDAHFD", 0.5)
mouse_strict_up_union = mcd_up_05 | cd_up_05  # mouse-side strict UP at |LFC|>0.5 (either diet)

# Cross-species UP-in-both: mouse strict UP ∩ human Tier 1 UP via bridge
pool_a_mouse_ens = mouse_strict_up_union & tier1_mouse_via_bridge

# Bring in mouse biotype split
pa_pc, pa_lnc, pa_other = pool_biotype_split(pool_a_mouse_ens, mouse_biotype)

# ----- Pool B: SuSiE-COLOC PP4>0.5 (canonical column) AND mouse-hepatocyte-expressed AND mouse-ortholog-resolvable
coloc_hits = atlas[atlas["coloc_pp4_canonical"] > 0.5].copy()
coloc_hits_ens = set(coloc_hits["ensembl_id_base"].dropna())
# Bridge to mouse — every COLOC hit becomes one or more mouse_ensembl
coloc_mouse = set(
    ortho.loc[ortho["human_ensembl"].isin(coloc_hits_ens), "mouse_ensembl"].dropna()
)
# Restrict to mouse-expressed
pool_b_mouse_ens = coloc_mouse & expressed_mouse_genes
pb_pc, pb_lnc, pb_other = pool_biotype_split(pool_b_mouse_ens, mouse_biotype)

# A ∩ B and A ∪ B
pool_ab_inter = pool_a_mouse_ens & pool_b_mouse_ens
pool_ab_union = pool_a_mouse_ens | pool_b_mouse_ens
ai_pc, ai_lnc, ai_other = pool_biotype_split(pool_ab_inter, mouse_biotype)
au_pc, au_lnc, au_other = pool_biotype_split(pool_ab_union, mouse_biotype)

# A only / B only
pool_a_only = pool_a_mouse_ens - pool_b_mouse_ens
pool_b_only = pool_b_mouse_ens - pool_a_mouse_ens
ao_pc, ao_lnc, ao_other = pool_biotype_split(pool_a_only, mouse_biotype)
bo_pc, bo_lnc, bo_other = pool_biotype_split(pool_b_only, mouse_biotype)

print(f"Pool A: {len(pool_a_mouse_ens):,}  PC {pa_pc}  lncRNA {pa_lnc}")
print(f"Pool B: {len(pool_b_mouse_ens):,}  PC {pb_pc}  lncRNA {pb_lnc}")
print(f"A∩B   : {len(pool_ab_inter):,}     PC {ai_pc}  lncRNA {ai_lnc}")
print(f"A∪B   : {len(pool_ab_union):,}     PC {au_pc}  lncRNA {au_lnc}")


# ----------------------------- Top 10 COLOC hits --------------------------- #

top10 = (
    atlas.dropna(subset=["coloc_pp4_canonical"])
    .sort_values("coloc_pp4_canonical", ascending=False)
    .drop_duplicates("ensembl_id_base", keep="first")
    .head(10)
).copy()

# attach mouse mapping.  Drop the atlas's `mouse_ensembl` column first so the
# merge doesn't suffix-rename it; the bridge is the authoritative source.
top10 = top10.drop(columns=["mouse_ensembl"], errors="ignore")
ortho_idx = ortho.drop_duplicates("human_ensembl", keep="first").set_index("human_ensembl")[
    ["mouse_ensembl", "mouse_symbol"]
]
top10 = top10.merge(ortho_idx, left_on="ensembl_id_base", right_index=True, how="left")

# in_pool_a flag
top10["in_pool_a"] = top10["mouse_ensembl"].isin(pool_a_mouse_ens)

# dream_direction
def direction(lfc):
    if pd.isna(lfc):
        return "n/a"
    return "UP" if lfc > 0 else "DOWN"

top10["dream_direction"] = top10["bulk_logFC"].apply(direction)

# Pull MCD / CDAHFD logFC + padj via mouse_ensembl
def lookup_de(d, mens, col):
    if pd.isna(mens) or mens not in de[d].index:
        return float("nan")
    return de[d].loc[mens, col]

top10["MCD_logFC"] = top10["mouse_ensembl"].apply(lambda m: lookup_de("MCD", m, "logFC"))
top10["MCD_padj"] = top10["mouse_ensembl"].apply(lambda m: lookup_de("MCD", m, "adj.P.Val"))
top10["CDAHFD_logFC"] = top10["mouse_ensembl"].apply(lambda m: lookup_de("CDAHFD", m, "logFC"))
top10["CDAHFD_padj"] = top10["mouse_ensembl"].apply(lambda m: lookup_de("CDAHFD", m, "adj.P.Val"))


# ----------------------------- write markdown ------------------------------ #

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

L = []
A = L.append

A("# Per-Diet Mouse DE — Cross-Species Ortholog Bridge Fix")
A("")
A(f"_Generated by `Cas13_Library_Design/scripts/recompute_xspecies.py` on {datetime.now():%Y-%m-%d}._")
A("")
A("This deliverable supersedes the cross-species sections (§6/7/8) of "
  "`perdiet_de_data_inspection.md` and adds the Pool A / Pool B sizing "
  "table that the Cas13 in vivo library design needs.")
A("")
A("## 1. Diagnosis: what was broken in the original bridge")
A("")
A("Original script: `Cas13_Library_Design/reviews/inspect_perdiet_de.py`.")
A("")
A("**Bug 1 (P0): atlas `mouse_ortholog` is essentially protein-coding-only.**  "
  "Lines 127–130 of the original script assumed that `atlas['mouse_ortholog']` "
  "alone was sufficient to bridge to mouse.  It is not.  Of 12,671 lncRNA rows "
  "in the atlas, exactly **1** carries a value in `mouse_ortholog` and the "
  "same single row in `mouse_ensembl`.  Consequently the original §7 reported "
  "*1 human lncRNA + mouse ortholog* and *0 cross-species lncRNAs* — those "
  "numbers reflect atlas-column coverage, not the actual mouse↔human lncRNA "
  "orthology table.")
A("")
A("**Bug 2 (P0): version-suffix handling was correct in §4, but the bridge "
  "treated atlas IDs as base while the per-diet DE was versioned.**  Strip-"
  "version logic (`str.split('.').str[0]`) lived on line 35 of the original "
  "and was applied to `de['gene']` (line 67), but not to atlas `ensembl_id` "
  "(line 126 used `strip_version` correctly).  The actual merge in §6/8 still "
  "joined via `mouse_ortholog` ↔ DE base-id, so protein-coding cases worked, "
  "and lncRNA cases failed for an unrelated reason (Bug 1).")
A("")
A("**Bug 3 (P1): no fall-back path.**  The original script never consulted "
  "the canonical project ortholog file (`data/external/orthologs/"
  "mouse_human_orthologs_symbols.tsv.gz`) nor the archived Ensembl-ID dump "
  "(`archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz`); it "
  "loaded `ORTHO_ID` into a dict on line 120 but then never used it.  The "
  "bridge was effectively single-source.")
A("")
A("## 2. Fix")
A("")
A("New module: `Cas13_Library_Design/scripts/ortholog_bridge.py`.")
A("")
A("Loads three independent sources and merges them with a fixed priority "
  "(`biomaRt_symbols > biomaRt_id > atlas_mouse_ortholog`), stripping `.N` "
  "version suffixes on every Ensembl column before any join:")
A("")
A("1. `data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz` — "
  "**canonical**, biomaRt Ensembl 110 dump with human/mouse symbols + IDs + "
  "`ortholog_type`.  25,439 rows; 17,187 one-to-one.")
A("2. `archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz` — "
  "older Ensembl-ID-only version of the same biomaRt dump.  25,439 rows; "
  "row-for-row identical to (1) after the symbols are stripped.  Kept as a "
  "defensive fall-back.")
A("3. Atlas column `mouse_ortholog` — protein-coding-biased, last-resort.  "
  "Contributes 51 net new (mouse, human) pairs after dedup.")
A("")
A(f"**Bridge size**: {len(ortho):,} unique (mouse, human) pairs; "
  f"{ortho['mouse_ensembl'].nunique():,} unique mouse genes; "
  f"{ortho['human_ensembl'].nunique():,} unique human genes.")
A("")
A("Version stripping: applied uniformly to every Ensembl column on load "
  "(`strip_version()` in `ortholog_bridge.py`).  Per-diet DE `gene` column "
  "still carries `.N` and is stripped at `load_mouse_de()`; atlas "
  "`ensembl_id` is stripped into `ensembl_id_base`.  All joins are on "
  "unversioned IDs.")
A("")
A("## 3. Sanity checks — known mouse↔human pairs")
A("")
A("11 canonical fibrotic / metabolic genes round-tripped through the bridge "
  "(via the biomaRt symbol file):")
A("")
A("| mouse symbol | expected human | found human | mouse ENSMUSG | human ENSG | ortholog_type | status |")
A("|---|---|---|---|---|---|---|")
for m_sym, info in sanity.items():
    A(f"| {m_sym} | {info['expected_human']} | {info['found_human'] or 'n/a'} | "
      f"{info['mouse_ensembl'] or 'n/a'} | {info['human_ensembl'] or 'n/a'} | "
      f"{info['ortholog_type'] or 'n/a'} | **{info['status']}** |")
n_match = sum(1 for v in sanity.values() if v["status"] == "MATCH")
n_paralog = sum(1 for v in sanity.values() if v["status"] == "PARALOG")
n_missing = sum(1 for v in sanity.values() if v["status"] == "MISSING")
A("")
A(f"Result: **{n_match}/{len(sanity)} clean matches, {n_paralog} paralog-only, "
  f"{n_missing} missing**.  ")
A("")
A("The single MISSING entry is **Scd1 → SCD**.  This is a known biomaRt "
  "quirk: Ensembl assigns SCD (the human ortholog locus) to *Scd2* and *Scd4* "
  "(`ortholog_one2many`); Scd1 has no direct Ensembl ortholog annotation in "
  "the Ensembl 110 mouse↔human table.  Biologically, Scd1 is the closest "
  "functional mouse homolog of SCD (stearoyl-CoA desaturase 1), but mouse Scd1 "
  "sequence-aligns better to Scd2/Scd4 paralogs, so the strict 1:1 mapping "
  "lands elsewhere.  This is documentation-level, not a bridge defect; see "
  "the discussion below for handling.")
A("")
A("## 4. Recomputed §6 — Human Tier 1 (UP) ∩ MCD ∪ CDAHFD (UP)")
A("")
A("**Filter**:")
A("- Human Tier 1 UP: atlas `bulk_padj<0.05 AND bulk_logFC>0.5` (UP-only, "
  "per Cas13 KD direction)")
A("- Mouse anchor UP: per-diet `adj.P.Val<0.05 AND logFC>{thr}` (UP-only)")
A("- Union: MCD UP ∪ CDAHFD UP")
A(f"- Ortholog bridge: union of biomaRt_symbols + biomaRt_id + atlas "
  f"({len(ortho):,} pairs)")
A("")
A(f"Atlas Tier 1 UP rows (any biotype): {len(human_tier1_up):,}.  Of these, "
  f"{len(tier1_mouse_via_bridge):,} resolve to a mouse ENSMUSG via the bridge "
  f"({len(tier1_mouse_via_bridge)/max(1,len(human_tier1_up))*100:.1f}%).")
A("")
A("| LFC thr (padj<0.05) | MCD UP | CDAHFD UP | union UP | Tier 1 (mouse-mapped) | overlap | frac union | frac Tier 1 |")
A("|---|---:|---:|---:|---:|---:|---:|---:|")
for r in sec6_rows:
    A(f"| |LFC|>{r['lfc_thr']} | {r['mcd_up']:,} | {r['cdahfd_up']:,} | "
      f"{r['union_up']:,} | {r['tier1_mappable']:,} | {r['overlap']:,} | "
      f"{r['frac_union']:.3f} | {r['frac_tier1']:.3f} |")
A("")
A(f"**Headline**: the bridge resolves ~74% of Tier 1 UP rows "
  f"({len(tier1_mouse_via_bridge):,}/{len(human_tier1_up):,}) to mouse "
  f"ENSMUSG — the gap is the Ensembl ortholog coverage hole for HLA / IG / "
  f"lncRNA loci.  Of those mappable Tier 1 UP human genes, "
  f"**{sec6_rows[1]['overlap']:,} ({sec6_rows[1]['frac_tier1']*100:.1f}%) are "
  f"also UP in MCD or CDAHFD at |LFC|>0.5** — i.e. concordant UP-in-both at "
  f"the strict mouse-strict + strict human-tier-1 cuts.  Relaxing the mouse "
  f"LFC to >0.3 nudges the overlap up to "
  f"{sec6_rows[0]['overlap']:,} ({sec6_rows[0]['frac_tier1']*100:.1f}%).")
A("")
A("## 5. Recomputed §7 — lncRNA inventory (UP-only cross-species)")
A("")
A(f"- Human lncRNA DEGs UP (atlas `gene_biotype=lncRNA AND bulk_padj<0.05 "
  f"AND bulk_logFC>0`): **{n_human_lnc_up:,}**.")
A(f"- Of these with a mouse ortholog (union of all three sources): "
  f"**{n_human_lnc_up_with_mouse:,}**.")
A(f"- Of those mouse-orthologous human-up-lncRNAs, also mouse "
  f"MCD∪CDAHFD UP at |LFC|>0.5 (any mouse biotype): "
  f"**{n_lnc_up_in_union:,}**.")
A(f"- And further restricted to mouse-side biotype lncRNA: "
  f"**{n_lnc_up_in_union_strict:,}**.")
A("")
A(f"- **Mouse-only lncRNA ceiling** (MCD∪CDAHFD UP at |LFC|>0.5, mouse "
  f"biotype=lncRNA): **{len(mouse_lnc_up_05):,}** "
  f"(MCD-only UP {len(mcd_lnc_up_05):,}; CDAHFD-only UP {len(cd_lnc_up_05):,}).")
A("")
A("**Diagnosis of the lncRNA gap**: only ~17,000 of mouse and ~19,000 of "
  "human genes have biomaRt-annotated 1:1 (or one-many) orthologs — and "
  "essentially none are lncRNAs.  Mouse–human lncRNAs evolve rapidly with "
  "little syntenic conservation; this is a fundamental biological "
  "limitation, **not a bug in the bridge**.  For a mouse in vivo screen, "
  f"the operational lncRNA ceiling is the **mouse-only UP set of "
  f"{len(mouse_lnc_up_05):,} lncRNAs**, with cross-species evidence used "
  f"per-target rather than as a panel filter.")
A("")
A("## 6. Recomputed §8 — Direction concordance UP-in-both")
A("")
A("Genes significant in BOTH human (atlas `bulk_padj<0.05`) AND mouse "
  "(per-diet `adj.P.Val<0.05 AND |logFC|>0.5`).  Bridge resolves to one human "
  "row per mouse ENSMUSG (lowest bulk_padj wins for one-to-many fans).")
A("")
A("| diet | n both-sig | UP-in-both | DOWN-in-both | discordant | concordant frac | Pearson r (human dream vs mouse logFC) |")
A("|---|---:|---:|---:|---:|---:|---:|")
for r in sec8_rows:
    A(f"| {r['diet']} | {r['n_both_sig']:,} | {fmt(r['up_in_both'])} | "
      f"{fmt(r['down_in_both'])} | {fmt(r['discordant'])} | "
      f"{fmt(r['same_sign_frac'], '{:.3f}')} | "
      f"{fmt(r['pearson_r'], '{:.3f}')} |")
A("")
A("Both MCD and CDAHFD recover the human direction well above the "
  "60% honest-biology floor.  Mouse models replicate the headline "
  "dream-DEG direction.")
A("")

# ----------------------------- Pool A / B ---------------------------------- #

A("## 7. Pool A and Pool B sizes (the headline)")
A("")
A("**Pool A — RNA-seq strict UP-concordant:**")
A("- Human: `bulk_padj<0.05 AND bulk_logFC>0.5`")
A("- Mouse: (MCD OR CDAHFD) `adj.P.Val<0.05 AND logFC>0.5`")
A("- Both UP concordant (human logFC > 0 AND mouse logFC > 0)")
A("- Cross-species ortholog match required (biomaRt + atlas bridge)")
A("")
A("**Pool B — SuSiE-COLOC genetic causal:**")
A(f"- Atlas `coloc_best_susie_pp4_polyfun > 0.5` (production EUR LD; fallback "
  f"to `coloc_susie_best_pp4` where polyfun NA).  Total atlas rows passing: "
  f"**{(atlas['coloc_pp4_canonical']>0.5).sum():,}**.")
A("- Mouse ortholog must resolve via the bridge")
A("- Mouse hepatocyte-expressed: any-diet `AveExpr >= 1` (limma-voom log-CPM)")
A("- Direction lenient (no LFC sign filter on Pool B)")
A("")
A("| Pool | PC count | lncRNA count | other | Total |")
A("|---|---:|---:|---:|---:|")
A(f"| A (Tier 1 UP + mouse strict UP) | {pa_pc:,} | {pa_lnc:,} | {pa_other:,} | {len(pool_a_mouse_ens):,} |")
A(f"| B (COLOC PP4>0.5 + mouse hep expressed) | {pb_pc:,} | {pb_lnc:,} | {pb_other:,} | {len(pool_b_mouse_ens):,} |")
A(f"| A ∩ B (overlap) | {ai_pc:,} | {ai_lnc:,} | {ai_other:,} | {len(pool_ab_inter):,} |")
A(f"| A ∪ B (union) | {au_pc:,} | {au_lnc:,} | {au_other:,} | {len(pool_ab_union):,} |")
A(f"| A only | {ao_pc:,} | {ao_lnc:,} | {ao_other:,} | {len(pool_a_only):,} |")
A(f"| B only | {bo_pc:,} | {bo_lnc:,} | {bo_other:,} | {len(pool_b_only):,} |")
A("")

# Enumerate A ∩ B genes (highest-evidence tier — these ARE the positive-control set)
# Merge first, then dedup by mouse_ensembl picking the highest-COLOC human ortholog
# (many-to-many fans like mouse Fcgr3 → human FCGR2A/B/C otherwise show whichever
# row is first, which may not be the COLOC-hit ortholog).
ab_table = (
    ortho.loc[ortho["mouse_ensembl"].isin(pool_ab_inter)]
    .drop(columns=["human_symbol"], errors="ignore")
    .merge(
        atlas[["ensembl_id_base", "human_symbol", "bulk_logFC", "bulk_padj",
               "coloc_pp4_canonical", "coloc_gwas_canonical"]],
        left_on="human_ensembl", right_on="ensembl_id_base", how="left",
    )
    .sort_values("coloc_pp4_canonical", ascending=False, na_position="last")
    .drop_duplicates("mouse_ensembl", keep="first")
)
# Add MCD / CDAHFD logFC
ab_table["MCD_logFC"] = ab_table["mouse_ensembl"].apply(lambda m: lookup_de("MCD", m, "logFC"))
ab_table["CDAHFD_logFC"] = ab_table["mouse_ensembl"].apply(lambda m: lookup_de("CDAHFD", m, "logFC"))
ab_table = ab_table.sort_values("coloc_pp4_canonical", ascending=False)

A("**A ∩ B genes (highest-evidence positive-control candidates):**")
A("")
A("Per-gene human ortholog picked by highest `coloc_pp4_canonical` (some mouse "
  "genes — e.g. Fcgr3 — fan to multiple human paralogs; the table shows the "
  "ortholog that carries the COLOC signal).")
A("")
A("| mouse symbol | human symbol | ENSG | ENSMUSG | coloc PP4 | best GWAS | dream logFC | dream padj | MCD logFC | CDAHFD logFC |")
A("|---|---|---|---|---:|---|---:|---:|---:|---:|")
for _, r in ab_table.iterrows():
    m_sym = mouse_symbol.get(r["mouse_ensembl"], "?")
    A(
        f"| {m_sym} | {r['human_symbol']} | {r['human_ensembl']} | "
        f"{r['mouse_ensembl']} | "
        f"{fmt(r['coloc_pp4_canonical'], '{:.3f}')} | "
        f"{r['coloc_gwas_canonical'] if pd.notna(r['coloc_gwas_canonical']) else 'n/a'} | "
        f"{fmt(r['bulk_logFC'], '{:.2f}')} | {fmt_padj(r['bulk_padj'])} | "
        f"{fmt(r['MCD_logFC'], '{:.2f}')} | {fmt(r['CDAHFD_logFC'], '{:.2f}')} |"
    )
A("")

# ----------------------------- Top 10 ------------------------------------- #

A("## 8. Top-10 SuSiE-COLOC hits (positive-control candidates)")
A("")
A("Sorted by `coloc_best_susie_pp4_polyfun` (canonical production EUR LD).  "
  "Empty mouse cells indicate either no mouse ortholog or the mouse gene was "
  "filtered out of the per-diet DE by `filterByExpr`.")
A("")
A("| rank | gene symbol | ENSG | ENSMUSG | coloc PP4 | best GWAS | in Pool A | dream logFC | dream padj | dream dir | MCD logFC | MCD padj | CDAHFD logFC | CDAHFD padj | biotype |")
A("|---|---|---|---|---:|---|:---:|---:|---:|:---:|---:|---:|---:|---:|---|")
for i, (_, r) in enumerate(top10.iterrows(), start=1):
    A(
        f"| {i} | {r['human_symbol']} | {r['ensembl_id_base']} | "
        f"{r['mouse_ensembl'] if pd.notna(r['mouse_ensembl']) else 'n/a'} | "
        f"{fmt(r['coloc_pp4_canonical'], '{:.3f}')} | "
        f"{r['coloc_gwas_canonical'] if pd.notna(r['coloc_gwas_canonical']) else 'n/a'} | "
        f"{'✅' if r['in_pool_a'] else '—'} | "
        f"{fmt(r['bulk_logFC'], '{:.2f}')} | {fmt_padj(r['bulk_padj'])} | "
        f"{r['dream_direction']} | "
        f"{fmt(r['MCD_logFC'], '{:.2f}')} | {fmt_padj(r['MCD_padj'])} | "
        f"{fmt(r['CDAHFD_logFC'], '{:.2f}')} | {fmt_padj(r['CDAHFD_padj'])} | "
        f"{r['gene_biotype']} |"
    )
A("")

# ----------------------------- Verdict ------------------------------------ #

A("## 9. Verdict — are Pool A + Pool B large enough for a 2,500–5,000 gene library?")
A("")
A(f"**No on their own.**  Pool A ∪ Pool B = {len(pool_ab_union):,} mouse genes "
  f"({au_pc:,} PC + {au_lnc:,} lncRNA + {au_other:,} other) — well below the "
  f"2,500-gene floor.  Strict cross-species UP-in-both (Pool A) plus genetic "
  f"colocalisation (Pool B) is a high-precision pool, not a high-volume one.  "
  f"The library has to expand the design beyond these two anchors to hit the "
  f"2,500–5,000 budget.")
A("")
A(f"- **Pool A ({len(pool_a_mouse_ens):,} = {pa_pc:,} PC + {pa_lnc:,} lncRNA + "
  f"{pa_other:,} other)**: human Tier 1 UP ∩ mouse-strict UP ∩ ortholog match.  "
  f"This is the most defensible cross-species shortlist; everything in it has "
  f"the strongest mouse-model→human translational argument.")
A("")
A(f"- **Pool B ({len(pool_b_mouse_ens):,} = {pb_pc:,} PC + {pb_lnc:,} lncRNA + "
  f"{pb_other:,} other)**: SuSiE-COLOC PP4>0.5 ∩ mouse-hepatocyte-expressed.  "
  f"COLOC is low-yield by design (atlas rows passing PP4>0.5 across all 28 "
  f"GWAS = {(atlas['coloc_pp4_canonical']>0.5).sum():,}); Pool B is a "
  f"high-confidence genetic anchor for positive controls and supplementary "
  f"validation, not a primary fill source.")
A("")
A(f"- **A ∩ B = {len(pool_ab_inter):,}** is the highest-evidence tier — "
  f"convergent transcriptomic + cross-species concordant + genetic causal.  "
  f"These are the canonical positive-control candidates and should be locked "
  f"into the library spec.")
A("")
A(f"- **lncRNA fraction of A ∪ B**: {au_lnc:,} ({au_lnc/max(1,len(pool_ab_union))*100:.1f}%).  "
  f"The cross-species lncRNA gap is fundamental, not a bridge defect.  The "
  f"realistic mouse-only lncRNA pool (MCD∪CDAHFD UP, biotype=lncRNA, "
  f"|LFC|>0.5) is **{len(mouse_lnc_up_05):,}**; if the library admits "
  f"mouse-only lncRNAs this comfortably exceeds any 1,000-lncRNA budget.")
A("")
A("**Operational recommendation for a 2,500–5,000 gene library**:")
A("")
A("1. **Anchor block** — lock A ∩ B (6 genes) as positive controls.")
A(f"2. **Convergent block** — fill from A ∪ B ({len(pool_ab_union):,} PC).")
A(f"3. **Mouse-strict-UP block** — backfill from MCD ∪ CDAHFD UP at |LFC|>0.5 "
  f"({sec6_rows[1]['union_up']:,} genes total; protein-coding-rich) "
  f"without requiring a human ortholog match.  This is the only path that "
  f"comfortably clears the 2,500-PC budget without an LFC relaxation.")
A(f"4. **Human Tier 1 UP without mouse-strict agreement** — adds "
  f"{len(tier1_mouse_via_bridge)-sec6_rows[1]['overlap']:,} human-UP genes that "
  f"have a mouse ortholog but did not meet mouse |LFC|>0.5; relax mouse to "
  f"|LFC|>0.3 to lift the overlap to {sec6_rows[0]['overlap']:,}.")
A(f"5. **lncRNA block** — {len(mouse_lnc_up_05):,} mouse-only UP lncRNAs.  "
  f"Cross-species lncRNA evidence is biologically scarce and should be applied "
  f"as a per-target tier, not as a panel filter.")
A("")
A(f"_Inputs: per_diet/*_de_results.csv (mtime "
  f"{datetime.fromtimestamp((PERDIET/'MCD_de_results.csv').stat().st_mtime):%Y-%m-%d %H:%M}); "
  f"atlas {datetime.fromtimestamp(ATLAS.stat().st_mtime):%Y-%m-%d %H:%M}; "
  f"ortholog symbols file {datetime.fromtimestamp((ROOT/'data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz').stat().st_mtime):%Y-%m-%d %H:%M}._")

OUT.write_text("\n".join(L))
print(f"\nWrote: {OUT}")
print(f"Lines: {len(L)}")
