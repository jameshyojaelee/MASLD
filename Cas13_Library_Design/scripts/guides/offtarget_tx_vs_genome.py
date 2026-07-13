#!/usr/bin/env python
"""Transcriptome-level vs genome-level off-target profile for the Cas13 library guides.
Joins the library guide roster to sfriedman's .3.transcriptome (hard-filter pass list)
and .5.genome (all-guides + c0:c1:c2 annotation) files."""
import polars as pl, re

LIB = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design/data/guides/cas13_library_guides_vM38.csv"
SF  = "/gpfs/commons/groups/sanjana_lab/sfriedman/track-cas13/outputs/mouse_GENCODEvM38_v1/sgrna"
TX  = f"{SF}/gencode.vM38.3.transcriptome.csv"
GN  = f"{SF}/gencode.vM38.5.genome.csv"

strip = lambda c: pl.col(c).str.replace(r"\.\d+$", "")

lib = pl.read_csv(LIB)
lib = lib.with_columns((strip("gene_id_mouse") + "|" + pl.col("guide_seq")).alias("key"))
keys = set(lib["key"].to_list())
N = lib.height
print(f"library guides: {N}  (unique keys: {len(keys)})")

def pull(path, cols):
    lf = pl.scan_csv(path, separator=",", infer_schema_length=2000, ignore_errors=True)
    lf = lf.with_columns((strip("gene_id") + "|" + pl.col("guide_seq")).alias("key"))
    df = lf.filter(pl.col("key").is_in(keys)).select(["key"] + cols).collect(engine="streaming")
    return df.unique(subset=["key"])

# ---- transcriptome level (.3): the hard-filter PASS list ----
tx = pull(TX, ["offtarget_pseudogene_only", "gene_family_hit_genes", "gene_family_hit_transcripts"])
tx = tx.rename({c: f"tx_{c}" for c in tx.columns if c != "key"})
print(f"\n[.3 transcriptome] library guides found in tx PASS list: {tx.height}/{N} "
      f"({100*tx.height/N:.1f}%)  -> guides absent here FAILED the transcriptome filter")

# ---- genome level (.5): all guides + mismatch triple ----
gn = pull(GN, ["mismatch", "any_indel", "offtarget_pseudogene_only", "gene_family_hit_transcripts"])
gn = gn.rename({"offtarget_pseudogene_only": "gn_pseudo_only",
                "gene_family_hit_transcripts": "gn_gf_tx"})
print(f"[.5 genome]        library guides found: {gn.height}/{N}")

m = lib.join(gn, on="key", how="left").join(tx, on="key", how="left")

def c012(colexpr):
    parts = colexpr.str.split(":")
    return (parts.list.get(0).cast(pl.Int32, strict=False),
            parts.list.get(1).cast(pl.Int32, strict=False),
            parts.list.get(2).cast(pl.Int32, strict=False))

c0, c1, c2 = c012(pl.col("mismatch_right") if "mismatch_right" in m.columns else pl.col("mismatch"))
mcol = "mismatch_right" if "mismatch_right" in m.columns else "mismatch"
m = m.with_columns([c012(pl.col(mcol))[i].alias(n) for i, n in enumerate(["g0","g1","g2"])])

# cross-check: does the genome mismatch from .5 match what our library table already carries?
same = (m["mismatch"] == m[mcol]).sum() if mcol != "mismatch" else N
print(f"\n=== VERIFY genome-level 43.9% ===")
tot = m.height
g0pos = (m["g0"] > 0).sum(); g1pos = (m["g1"] > 0).sum(); g2pos = (m["g2"] > 0).sum()
zero = ((m["g0"]==0) & (m["g1"]==0) & (m["g2"]==0)).sum()
print(f"library-table mismatch == .5.genome.csv mismatch (provenance): {same}/{tot}")
print(f"genome 0:0:0 (no OT)          : {zero}  ({100*zero/tot:.1f}%)")
print(f"genome >=1 perfect (0mm) OT   : {g0pos}  ({100*g0pos/tot:.1f}%)   <-- the 43.9% claim")
print(f"genome >=1 @1mm               : {g1pos}  ({100*g1pos/tot:.1f}%)")
print(f"genome >=1 @2mm               : {g2pos}  ({100*g2pos/tot:.1f}%)")

# bool columns (offtarget_pseudogene_only) are inferred as Boolean by polars
as_bool = lambda s: m[s].cast(pl.Boolean, strict=False).fill_null(False)
nonempty = lambda s: m[s].cast(pl.Utf8, strict=False).fill_null("") != ""

# how many of the genome perfect-match OTs are pseudogene-only or gene-family (tolerable)?
gf_or_pseudo = ((m["g0"]>0) & (as_bool("gn_pseudo_only") | nonempty("gn_gf_tx"))).sum()
print(f"  of the {g0pos} genome-0mm-OT guides, pseudogene-only or gene-family flagged: {gf_or_pseudo}")

print(f"\n=== TRANSCRIPTOME-level off-target burden (.3 pass list) ===")
present = m["tx_offtarget_pseudogene_only"].is_not_null().sum()
print(f"guides present in tx PASS list           : {present}/{tot}  (rest failed tx filter -> not in library-eligible pool)")
tx_pseudo = as_bool("tx_offtarget_pseudogene_only").sum()
tx_gf     = nonempty("tx_gene_family_hit_transcripts").sum()
tx_clean  = ((~as_bool("tx_offtarget_pseudogene_only")) &
             (~nonempty("tx_gene_family_hit_transcripts")) &
             (~nonempty("tx_gene_family_hit_genes"))).sum()
print(f"transcriptome pseudogene-only OT (tolerated) : {tx_pseudo}  ({100*tx_pseudo/tot:.1f}%)")
print(f"transcriptome gene-family/paralog tx hit     : {tx_gf}  ({100*tx_gf/tot:.1f}%)")
print(f"transcriptome FULLY CLEAN (no tx OT at all)  : {tx_clean}  ({100*tx_clean/tot:.1f}%)")
print(f"transcriptome REAL disqualifying OT          : 0 by construction (hard-filtered out upstream)")

# side-by-side burden
print(f"\n=== GENOME vs TRANSCRIPTOME (same {tot} library guides) ===")
print(f"  >=1 perfect-match OT   genome: {g0pos} ({100*g0pos/tot:.1f}%)   transcriptome(tolerated only): {tx_pseudo+tx_gf} ({100*(tx_pseudo+tx_gf)/tot:.1f}%)")

# ---- persist merged per-guide table for the figure ----
out = m.select([
    "guide_id", "gene_symbol_mouse", "gene_id_mouse", "biotype", "tier", "library_arm", "region",
    "n_target", "single_isoform", "tiger_score", "cas13_score", "combined_score",
    "mismatch", "g0", "g1", "g2", "any_indel", "gn_pseudo_only",
    "tx_offtarget_pseudogene_only", "tx_gene_family_hit_transcripts", "tx_gene_family_hit_genes",
]).rename({"mismatch": "genome_mismatch", "gn_pseudo_only": "genome_pseudo_only",
           "any_indel": "genome_any_indel"})
OUTP = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design/data/guides/cache/offtarget_merged.parquet"
out.write_parquet(OUTP)
print(f"\nmerged per-guide table -> {OUTP}  ({out.height} rows, {out.width} cols)")
