#!/usr/bin/env python
"""
76_decima_celltype.py  --  Phase 5, Track 2

Cell-type-resolved variant effect with Decima (Genentech, Nat Methods 2026;
s41592-026-03102-0). Decima is a single-cell-resolution sequence->expression
model with published weights (HF Genentech/decima-model, Zenodo 15092691).
It predicts a signed variant effect on gene expression *conditioned on cell
type* (8,856 pseudobulk tasks; 210 of them liver). We use it to say WHICH
liver cell type each of our eQTL-absent nomination + coding-effector variants
acts in, tying the genetic arm to the single-cell atlas.

MODALITY NOTE: Decima is a cell-type EXPRESSION variant-effect predictor
(modality-appropriate for scoring DNA variants). It is NOT a perturbation
foundation model; do not conflate.

ADDITIVE ONLY -- new script, reads existing seqfunc outputs, writes
SF/decima_celltype.tsv. Does not edit any existing script.

Anchor positive control: SORT1 rs12740374 (chr1 hg38 109274968 G>T) should act
most strongly in HEPATOCYTE (hepatic regulatory eQTL).
"""
import os, sys, glob
# isolate from ~/.local user-site (a decima dep leaked there during env build)
os.environ.setdefault("PYTHONNOUSERSITE", "1")
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

import numpy as np
import pandas as pd

REPO = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SF   = f"{REPO}/GWAS/finemapping/results/seqfunc"
ENVD = f"{SF}/decima_env"
META = f"{ENVD}/metadata.h5ad"                       # pre-downloaded HF metadata
ATAC = f"{REPO}/Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2"
OUT  = f"{SF}/decima_celltype.tsv"
# Confidence floor for the cell-type argmax: |effect| above the noise floor AND a clear
# top1/top2 separation. The all-variant argmax split is argmax-over-near-noise for most
# variants (median |effect| ~0.0017); report the split among confidently-assigned variants.
CONF_ABS    = 0.005
CONF_MARGIN = 1.2
# hg38 fasta (chr-prefixed, .fai-indexed) — pass a PATH so decima/grelu reads it
# directly instead of trying to genomepy-download "hg38". Same assembly the
# borzoi seqfunc arm used.
GENOME = "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

# HF cache local to seqfunc so we never touch $HOME
os.environ.setdefault("HF_HOME", f"{SF}/hf_cache")

# ---- coarse liver compartment mapping (Decima fine cell_type -> MASLD compartment)
COARSE = {
    "hepatocyte": "hepatocyte",
    "endothelial cell of hepatic sinusoid": "LSEC_endothelial",
    "blood vessel endothelial cell": "LSEC_endothelial",
    "capillary endothelial cell": "LSEC_endothelial",
    "endothelial cell of vascular tree": "LSEC_endothelial",
    "endothelial cell of lymphatic vessel": "LSEC_endothelial",
    "vein endothelial cell": "LSEC_endothelial",
    "fibroblast": "stellate_fibroblast",
    "myofibroblast cell": "stellate_fibroblast",
    "smooth muscle cell": "stellate_fibroblast",
    "vascular associated smooth muscle cell": "stellate_fibroblast",
    "macrophage": "kupffer_macrophage",
    "classical monocyte": "kupffer_macrophage",
    "non-classical monocyte": "kupffer_macrophage",
    "conventional dendritic cell": "kupffer_macrophage",
    "plasmacytoid dendritic cell": "kupffer_macrophage",
    "mast cell": "kupffer_macrophage",
    "cholangiocyte": "cholangiocyte",
    "epithelial cell of bile duct": "cholangiocyte",
    # everything else lymphoid (T/NK/B/plasma/ILC)
}
def to_coarse(ct):
    ct = str(ct)
    if ct in COARSE:
        return COARSE[ct]
    lym = ("T cell", "NK", "NKT", "B cell", "plasma", "lymphoid", "lymphocyte",
           "CD4", "CD8", "regulatory T", "innate lymphoid", "mucosal invariant",
           "plasmablast", "IgG", "memory B", "naive B", "germinal")
    if any(k.lower() in ct.lower() for k in lym):
        return "lymphoid_T_NK_B"
    return "other"

# ATAC compartment -> peak bed basename
ATAC_BED = {
    "hepatocyte": "Hepatocytes_peaks.bed",
    "LSEC_endothelial": "Endothelial_cells_peaks.bed",
    "stellate_fibroblast": "Fibroblasts_peaks.bed",
    "kupffer_macrophage": "Macrophages_peaks.bed",
    "cholangiocyte": "Cholangiocytes_peaks.bed",
    "lymphoid_T_NK_B": ["T_cells_peaks.bed", "Circulating_NK_NKT_peaks.bed",
                         "Resident_NK_peaks.bed", "Plasma_cells_peaks.bed"],
}


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- input assembly
def load_variants(decima_genes):
    """Assemble the union of eQTL-absent nominations + magnitude leads +
    coding effectors + SORT1 anchor, in hg38 chrom/pos/ref/alt/gene."""
    sub = pd.read_csv(f"{SF}/variant_substrate_hg38.tsv", sep="\t", dtype=str)
    sub_by_id = sub.drop_duplicates("variant_id_hg19").set_index("variant_id_hg19")

    rows = []

    # 1) eQTL-absent nominations (hg19 credible_variant -> hg38 via substrate)
    nom = pd.read_csv(f"{SF}/eqtl_absent_nominations.tsv", sep="\t", dtype=str)
    for _, r in nom.iterrows():
        vid = r["credible_variant"]
        if vid in sub_by_id.index:
            s = sub_by_id.loc[vid]
            rows.append(dict(chrom=str(s["chr"]), pos=int(float(s["pos_hg38"])),
                             ref=s["ref_hg38"], alt=s["alt_hg38"],
                             gene=r["nominated_gene"], var_set="nomination",
                             variant_id_hg19=vid))
    # 2) borzoi magnitude leads (already hg38)
    lead = pd.read_csv(f"{SF}/borzoi_magnitude_leads.tsv", sep="\t", dtype=str)
    for _, r in lead.iterrows():
        g = r.get("gene_substrate", "")
        if not isinstance(g, str) or g in ("", "-", "nan"):
            g = np.nan
        rows.append(dict(chrom=str(r["chr"]), pos=int(float(r["pos_hg38"])),
                         ref=r["ref_hg38"], alt=r["alt_hg38"], gene=g,
                         var_set="magnitude_lead", variant_id_hg19=r["variant_id_hg19"]))
    # 3) coding effectors (variant_hg38 = chr:pos:ref:alt already hg38)
    cod = pd.read_csv(f"{SF}/coding_hardening.tsv", sep="\t", dtype=str)
    for _, r in cod.iterrows():
        c, p, ref, alt = r["variant_hg38"].split(":")
        rows.append(dict(chrom=("chr"+c) if not str(c).startswith("chr") else str(c),
                         pos=int(p), ref=ref, alt=alt, gene=r["gene"],
                         var_set="coding_effector", variant_id_hg19=r["variant"]))
    # 4) SORT1 anchor positive control
    rows.append(dict(chrom="chr1", pos=109274968, ref="G", alt="T", gene="SORT1",
                     var_set="anchor", variant_id_hg19="rs12740374"))

    df = pd.DataFrame(rows)
    df["chrom"] = df["chrom"].where(df["chrom"].str.startswith("chr"), "chr"+df["chrom"])
    # merge var_set memberships across duplicated coordinates
    key = ["chrom", "pos", "ref", "alt"]
    df["gene"] = df["gene"].replace({"": np.nan, "-": np.nan, "nan": np.nan})
    agg = (df.sort_values("var_set")
             .groupby(key, dropna=False)
             .agg(gene=("gene", lambda s: next((x for x in s if isinstance(x, str)), np.nan)),
                  var_set=("var_set", lambda s: ",".join(sorted(set(s)))),
                  variant_id_hg19=("variant_id_hg19", "first"))
             .reset_index())
    agg["gene_in_decima"] = agg["gene"].isin(decima_genes)
    log(f"[input] {len(agg)} unique variants "
        f"({agg['gene_in_decima'].sum()} with gene in decima, "
        f"{(~agg['gene_in_decima']).sum()} auto-overlap fallback)")
    return agg


# ---------------------------------------------------------------- decima run
def run_decima(agg, liver_tasks, device="cuda"):
    from decima import predict_variant_effect

    def _call(df_in, gene_col):
        d = df_in[["chrom", "pos", "ref", "alt"]].copy()
        if gene_col:
            d["gene"] = df_in["gene"].values
        d = d.reset_index(drop=True)
        res = predict_variant_effect(
            df_variant=d,
            output_pq=None,
            tasks=None,                    # all tasks; we subset liver ourselves
            model="ensemble",
            metadata_anndata=META,
            genome=GENOME,
            gene_col=("gene" if gene_col else None),
            device=device,
            num_workers=4,
            batch_size=1,
        )
        return res

    parts = []
    with_gene = agg[agg["gene_in_decima"]]
    no_gene   = agg[~agg["gene_in_decima"]]

    if len(with_gene):
        log(f"[decima] scoring {len(with_gene)} variants WITH nominated gene ...")
        r = _call(with_gene, gene_col=True)
        r["_mode"] = "nominated_gene"
        parts.append(r)
    if len(no_gene):
        log(f"[decima] scoring {len(no_gene)} variants via auto TSS-overlap ...")
        r = _call(no_gene, gene_col=False)
        r["_mode"] = "auto_overlap"
        parts.append(r)

    res = pd.concat(parts, ignore_index=True)
    log(f"[decima] raw prediction rows: {len(res)} (cols incl tasks: {res.shape[1]})")
    return res


# ---------------------------------------------------------------- post-process
def summarize(res, liver_meta):
    """res: per (variant, gene) rows with one column per task. Collapse liver
    tasks to per-coarse-compartment signed + abs effect; pick argmax gene per
    variant (largest liver abs effect); return per-variant summary."""
    liver_tasks = list(liver_meta.index)
    present = [t for t in liver_tasks if t in res.columns]
    log(f"[post] liver tasks present in output: {len(present)}/{len(liver_tasks)}")

    # coarse compartment for each present task
    coarse = liver_meta.loc[present, "coarse"]
    comps = sorted(set(coarse) - {"other"})

    key = ["chrom", "pos", "ref", "alt"]
    # per-row (variant,gene): coarse signed mean + abs mean, overall liver abs
    for comp in comps:
        cols = [t for t in present if coarse.loc[t] == comp]
        block = res[cols].astype(float)
        res[f"eff_{comp}"] = block.mean(axis=1)
        res[f"abs_{comp}"] = block.abs().mean(axis=1)
    abs_cols = [f"abs_{c}" for c in comps]
    res["liver_abs_max"] = res[abs_cols].max(axis=1)

    # choose the gene-row with the largest liver abs effect per variant
    res = res.sort_values("liver_abs_max", ascending=False)
    top = res.groupby(key, as_index=False).first()

    # argmax compartment
    top["argmax_celltype"] = top[abs_cols].idxmax(axis=1).str.replace("abs_", "", regex=False)
    top["argmax_abs_effect"] = top[abs_cols].max(axis=1)
    # signed effect in the argmax compartment
    top["argmax_signed_effect"] = top.apply(
        lambda r: r[f"eff_{r['argmax_celltype']}"], axis=1)
    # top1/top2 margin + confidence floor. Most variants sit near the noise floor
    # (median |effect| ~0.0017); the all-120 argmax split is argmax-over-near-noise.
    # A magnitude+margin floor (|effect|>=CONF_ABS AND top1/top2>=CONF_MARGIN) marks the
    # variants whose compartment assignment is above noise — report the split among THESE.
    def _second_max(r):
        v = sorted((abs(r[c]) for c in abs_cols), reverse=True)
        return v[1] if len(v) > 1 else 0.0
    top["argmax_margin"] = top.apply(
        lambda r: (r["argmax_abs_effect"] / _second_max(r)) if _second_max(r) > 0 else np.inf,
        axis=1)
    top["confident_call"] = (top["argmax_abs_effect"] >= CONF_ABS) & (top["argmax_margin"] >= CONF_MARGIN)
    return top, comps


# ---------------------------------------------------------------- snATAC corroboration
def load_atac_peaks():
    peaks = {}   # compartment -> dict(chrom -> np.array of [start,end])
    for comp, beds in ATAC_BED.items():
        beds = beds if isinstance(beds, list) else [beds]
        by_chrom = {}
        for b in beds:
            fp = os.path.join(ATAC, b)
            if not os.path.exists(fp):
                continue
            bed = pd.read_csv(fp, sep="\t", header=None, usecols=[0, 1, 2],
                              names=["chrom", "start", "end"])
            for ch, g in bed.groupby("chrom"):
                arr = g[["start", "end"]].values
                by_chrom.setdefault(ch, []).append(arr)
        peaks[comp] = {ch: np.vstack(v) for ch, v in by_chrom.items()}
    return peaks


def atac_overlap(chrom, pos, peaks_comp):
    """Return set of compartments whose accessible peaks contain pos (hg38, 1-based)."""
    hits = set()
    for comp, by_chrom in peaks_comp.items():
        arr = by_chrom.get(chrom)
        if arr is None:
            continue
        # bed is 0-based half-open; variant pos 1-based -> pos-1
        p = pos - 1
        if np.any((arr[:, 0] <= p) & (p < arr[:, 1])):
            hits.add(comp)
    return hits


# ---------------------------------------------------------------- main
def main():
    import anndata as ad
    log("[meta] loading liver task metadata ...")
    a = ad.read_h5ad(META, backed="r")
    obs = a.obs
    liver = obs[(obs["organ"].astype(str).str.contains("liver", case=False, na=False)) |
                (obs["tissue"].astype(str).str.contains("liver", case=False, na=False))].copy()
    liver["coarse"] = liver["cell_type"].map(to_coarse)
    liver_meta = liver[["cell_type", "coarse"]]
    decima_genes = set(map(str, a.var_names))
    log(f"[meta] {len(liver)} liver tasks; compartments: "
        f"{liver['coarse'].value_counts().to_dict()}")

    agg = load_variants(decima_genes)

    res = run_decima(agg, liver_meta)
    # attach input metadata back onto res rows
    key = ["chrom", "pos", "ref", "alt"]
    res = res.merge(agg[key + ["gene", "var_set", "variant_id_hg19", "gene_in_decima"]],
                    on=key, how="left", suffixes=("", "_input"))

    top, comps = summarize(res, liver_meta)
    # var_set / variant_id_hg19 / gene_in_decima already carried through res->summarize;
    # merge them fresh only if summarize dropped them (avoid _x/_y suffix collision).
    for c in ["var_set", "variant_id_hg19", "gene_in_decima"]:
        if c not in top.columns:
            top = top.merge(agg[key + [c]], on=key, how="left")

    # snATAC corroboration
    log("[atac] loading per-cell-type peak sets ...")
    peaks_comp = load_atac_peaks()
    atac_hits, argmax_match = [], []
    for _, r in top.iterrows():
        hits = atac_overlap(r["chrom"], int(r["pos"]), peaks_comp)
        atac_hits.append(",".join(sorted(hits)) if hits else "")
        argmax_match.append(r["argmax_celltype"] in hits)
    top["atac_accessible_celltypes"] = atac_hits
    top["argmax_matches_atac"] = argmax_match
    top["atac_any_peak"] = [h != "" for h in atac_hits]

    # tidy output (explicit, ordered)
    if "gene" not in top.columns and "gene_input" in top.columns:
        top["gene"] = top["gene_input"]
    out_cols = (["chrom", "pos", "ref", "alt", "gene", "variant_id_hg19",
                 "var_set", "gene_in_decima"]
                + [f"eff_{c}" for c in comps] + [f"abs_{c}" for c in comps]
                + ["argmax_celltype", "argmax_abs_effect", "argmax_signed_effect",
                   "argmax_margin", "confident_call",
                   "atac_accessible_celltypes", "argmax_matches_atac", "atac_any_peak"])
    top = top[[c for c in out_cols if c in top.columns]].sort_values(
        "argmax_abs_effect", ascending=False)
    top.to_csv(OUT, sep="\t", index=False, float_format="%.5f")
    log(f"[done] wrote {OUT}  ({len(top)} variants)")

    # ---- report anchors + distribution
    log("\n===== SORT1 ANCHOR =====")
    s = top[top["gene"] == "SORT1"]
    if len(s):
        r = s.iloc[0]
        log(f"SORT1 rs12740374: argmax={r['argmax_celltype']} "
            f"(abs={r['argmax_abs_effect']:.4f}); "
            f"atac_accessible={r['atac_accessible_celltypes']}")
        for c in comps:
            log(f"    {c:22s} eff={r[f'eff_{c}']:+.4f} abs={r[f'abs_{c}']:.4f}")
    else:
        log("SORT1 not scored (check gene/overlap).")

    log("\n===== argmax cell-type distribution (ALL variants) =====")
    log(top["argmax_celltype"].value_counts().to_string())
    hep_all = int((top["argmax_celltype"] == "hepatocyte").sum())
    log(f"ALL: {hep_all}/{len(top)} hepatocyte = {100*hep_all/len(top):.1f}% "
        f"hep / {100-100*hep_all/len(top):.1f}% non-parenchymal")

    log(f"\n===== ROBUSTNESS: confident subset (|effect|>={CONF_ABS} & "
        f"margin>={CONF_MARGIN}) =====")
    conf = top[top["confident_call"]]
    if len(conf):
        hep_c = int((conf["argmax_celltype"] == "hepatocyte").sum())
        log(top.loc[top["confident_call"], "argmax_celltype"].value_counts().to_string())
        log(f"CONFIDENT: {hep_c}/{len(conf)} hepatocyte = {100*hep_c/len(conf):.1f}% "
            f"hep / {100-100*hep_c/len(conf):.1f}% non-parenchymal")
        log("NOTE: the all-variant 69% non-parenchymal split is NOT robust to a confidence "
            "floor — among confidently-assigned variants the split is ~balanced. Report the "
            "cell-type localization as nomination/mechanism-framed, NOT a quantitative "
            "hep-vs-nonparenchymal proportion; SORT1 anchor + snATAC check are the load-bearing "
            "corroboration.")
    else:
        log("no variants pass the confidence floor.")
    log(f"\nsnATAC corroboration: argmax matches an accessible cell type in "
        f"{top['argmax_matches_atac'].sum()}/{top['atac_any_peak'].sum()} "
        f"variants that overlap ANY peak "
        f"({top['argmax_matches_atac'].sum()}/{len(top)} of all).")


if __name__ == "__main__":
    main()
