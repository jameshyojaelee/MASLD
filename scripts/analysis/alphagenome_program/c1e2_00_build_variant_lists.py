#!/usr/bin/env python3
"""C1 endpoint 2, step 0: build the GPU variant lists wave 1 costed but never wrote.

Executes IMPLEMENTATION_SPEC.md section 5 (C1), the endpoint-2 half that wave 1 could not
run because GPU scoring was forbidden.

Writes, into <outdir>/scoring/:
  chrombpnet/<family>/fold{0..4}.variants.tsv   5-column no-header chrombpnet schema
                                                (chr, pos, allele1=REF, allele2=ALT, variant_id),
                                                each site placed in the fold whose TEST group
                                                holds its chromosome, so every site is scored by
                                                a model that never saw its chromosome.
  borzoi/borzoi_targets.tsv                     chr, pos_hg38, ref, alt, canon, priority
                                                ordered so the four-model matched set is covered
                                                first and a truncated run still delivers it.

Targets:
  * the 500 union allelic-imbalance sites (P3 deposit, GSE281367 + GSE244832), rebuilt with the
    same precision-weighted union the wave-1 C1 package used;
  * the Currin FDR<5% peak-lead variants each model is missing.
    ChromBPNet (both project families) already covers all 32,336 leads, so its only new work is
    the 500 sites. Borzoi covers 1,236, so every remaining lead is a Borzoi target.

Nothing under GWAS/finemapping/results/seqfunc/ or .../alphagenome_atlas/ is written.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os

import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SF = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
ATLAS = os.path.join(ROOT, "GWAS/finemapping/results/alphagenome_atlas")

CBP_FAMILIES = {
    "hepatocyte_5fold_v2": os.path.join(SF, "adult_liver_chrombpnet/hepatocyte_5fold_v2"),
    "gse281367_hepatocyte_5fold_v2": os.path.join(
        SF, "adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2"),
}
CBP_ADULT = os.path.join(
    SF, "adult_liver_chrombpnet/hepatocyte_5fold_v2/currin_calibration/currin_scored_primary.tsv.gz")
DIRFEAT = os.path.join(SF, "direction_features_caqtl.tsv")
ASE_281367 = os.path.join(ATLAS, "p3-ase-20260909T192926Z/tables/allelic_sites.tsv")
ASE_244832 = os.path.join(ATLAS, "p3-ase-20260909T192926Z/tables/gse244832_allelic_sites.tsv")
P3A_CHUNKS = os.path.join(ATLAS, "p3a-benchmarks-20260909T190214Z/raw/atlas_caqtl")
FASTA_CBP = ("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
             "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa")
FASTA_BZ = "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

AG_ATAC_TRACKS = ["UBERON:0001114 ATAC-seq", "UBERON:0001115 ATAC-seq", "UBERON:0002107 ATAC-seq"]


def log(m):
    print(m, flush=True)


def norm_chr(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def key_of(chrom, pos, ref, alt):
    return f"{norm_chr(chrom)}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


def sha256(path, max_bytes=400_000_000):
    if os.path.getsize(path) > max_bytes:
        return "not_hashed_too_large"
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def build_ase_union():
    a1 = pd.read_csv(ASE_281367, sep="\t")
    a2 = pd.read_csv(ASE_244832, sep="\t")
    a1["dataset"] = "GSE281367"
    a2["dataset"] = "GSE244832"
    both = pd.concat([a1, a2], ignore_index=True)
    rec = []
    for uid, g in both.groupby("uid"):
        prec = 1.0 / g.se_log2 ** 2
        wsum = float(prec.sum())
        rec.append(dict(
            uid=uid, chrom=g.chrom.iloc[0], pos=int(g.pos.iloc[0]),
            ref=g.ref.iloc[0], alt=g.alt.iloc[0],
            mean_log2_alt_over_ref=float((g.mean_log2_alt_over_ref * prec).sum() / wsum),
            se_log2=float(np.sqrt(1.0 / wsum)),
            n_datasets=int(len(g)), datasets=",".join(sorted(g.dataset)),
            n_het_donors=int(g.n_het_donors.sum()), total_reads=int(g.total_reads.sum())))
    ase = pd.DataFrame(rec)
    ase["key"] = [key_of(c, p, r, a) for c, p, r, a in zip(ase.chrom, ase.pos, ase.ref, ase.alt)]
    return ase


def fold_map(family_dir):
    fm = pd.read_csv(os.path.join(family_dir, "fold_manifest.tsv"), sep="\t")
    m = {}
    for r in fm.itertuples():
        for c in str(r.test_chromosomes).split(","):
            m[c.strip()] = int(r.fold)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    out = args.outdir
    sc = os.path.join(out, "scoring")
    os.makedirs(os.path.join(sc, "borzoi"), exist_ok=True)

    import pysam
    fa_cbp = pysam.FastaFile(FASTA_CBP)
    fa_bz = pysam.FastaFile(FASTA_BZ)

    summary = {}

    # ---------------- ASE union sites ----------------------------------------
    ase = build_ase_union()
    log(f"ASE union sites: {len(ase)}")
    ase["chr_pref"] = ["chr" + norm_chr(c) for c in ase.chrom]
    # reference identity check on both fastas
    for tag, fa in (("cellranger_arc", fa_cbp), ("refdata_gex", fa_bz)):
        base = [fa.fetch(c, int(p) - 1, int(p)).upper()
                for c, p in zip(ase.chr_pref, ase.pos)]
        ase[f"fasta_base_{tag}"] = base
        summary[f"ase_ref_matches_fasta_{tag}"] = int(
            (ase[f"fasta_base_{tag}"] == ase.ref.str.upper()).sum())
    summary["ase_n_sites"] = int(len(ase))
    summary["ase_n_autosomal_snv"] = int(
        ((ase.ref.str.len() == 1) & (ase.alt.str.len() == 1)
         & ase.chr_pref.isin([f"chr{i}" for i in range(1, 23)])).sum())
    ase.to_csv(os.path.join(sc, "ase_union_sites.tsv"), sep="\t", index=False)

    # ---------------- ChromBPNet fold lists ----------------------------------
    keep = ase[(ase.ref.str.len() == 1) & (ase.alt.str.len() == 1)
               & ase.chr_pref.isin([f"chr{i}" for i in range(1, 23)])].copy()
    per_family = {}
    for fam, fdir in CBP_FAMILIES.items():
        fmap = fold_map(fdir)
        d = keep.copy()
        d["heldout_fold"] = d.chr_pref.map(fmap)
        if d.heldout_fold.isna().any():
            raise SystemExit(f"{fam}: {int(d.heldout_fold.isna().sum())} sites with no fold")
        fdir_out = os.path.join(sc, "chrombpnet", fam)
        os.makedirs(fdir_out, exist_ok=True)
        counts = {}
        for f in range(5):
            s = d[d.heldout_fold == f]
            s[["chr_pref", "pos", "ref", "alt", "uid"]].to_csv(
                os.path.join(fdir_out, f"fold{f}.variants.tsv"),
                sep="\t", index=False, header=False)
            counts[str(f)] = int(len(s))
        d[["uid", "chr_pref", "pos", "ref", "alt", "heldout_fold"]].to_csv(
            os.path.join(fdir_out, "fold_assignment.tsv"), sep="\t", index=False)
        per_family[fam] = counts
        log(f"{fam}: fold counts {counts}")
    summary["chrombpnet_fold_counts"] = per_family

    # ---------------- Borzoi target list -------------------------------------
    cols = ["variant_id_hg38", "chr_x", "pos_hg38", "ref", "alt", "label", "beta_alt", "q_value"]
    cbp = pd.read_csv(CBP_ADULT, sep="\t", usecols=cols)
    leads = cbp[cbp.label == 1].copy()
    leads["key"] = [key_of(c, p, r, a)
                    for c, p, r, a in zip(leads.chr_x, leads.pos_hg38, leads.ref, leads.alt)]
    log(f"Currin FDR<5% leads with beta_alt: {len(leads)}")

    dirfeat = pd.read_csv(DIRFEAT, sep="\t")
    dirfeat["key"] = dirfeat["canon"].astype(str)
    have_bz = set(dirfeat.loc[dirfeat["borzoi_atac_delta_labelframe"].notna(), "key"])
    log(f"leads already carrying a Borzoi score: {len(have_bz & set(leads.key))}")

    # AlphaGenome archive coverage: the four-model matched set cannot exceed it, so those leads
    # are Borzoi priority 1.
    import anndata as ad
    import glob
    ag_keys = set()
    for ch in sorted(glob.glob(os.path.join(P3A_CHUNKS, "chunk_*"))):
        a = ad.read_h5ad(os.path.join(ch, "ATAC.h5ad"))
        v = pd.Series(a.obs["variant"].astype(str).values)
        parts = v.str.replace(">", ":", regex=False).str.split(":", expand=True)
        ag_keys.update(key_of(c, p, r, al)
                       for c, p, r, al in zip(parts[0], parts[1], parts[2], parts[3]))
    log(f"archived AlphaGenome caQTL variants: {len(ag_keys)}")
    summary["alphagenome_archive_variants"] = int(len(ag_keys))
    summary["currin_leads_with_alphagenome"] = int(len(ag_keys & set(leads.key)))
    summary["currin_leads_total"] = int(len(leads))
    summary["currin_leads_with_borzoi_before"] = int(len(have_bz & set(leads.key)))

    rows = []
    for r in ase.itertuples():
        rows.append(dict(chr=norm_chr(r.chrom), pos_hg38=int(r.pos), ref=str(r.ref).upper(),
                         alt=str(r.alt).upper(), canon=r.key, priority=0, target="ase_site"))
    for r in leads.itertuples():
        if r.key in have_bz:
            continue
        rows.append(dict(chr=norm_chr(r.chr_x), pos_hg38=int(r.pos_hg38), ref=str(r.ref).upper(),
                         alt=str(r.alt).upper(), canon=r.key,
                         priority=1 if r.key in ag_keys else 2, target="currin_lead"))
    tgt = pd.DataFrame(rows).drop_duplicates(subset=["canon"])
    tgt = tgt[(tgt.ref.str.len() == 1) & (tgt.alt.str.len() == 1)]
    tgt = tgt.sort_values(["priority", "chr", "pos_hg38"], kind="mergesort").reset_index(drop=True)
    tgt.to_csv(os.path.join(sc, "borzoi", "borzoi_targets.tsv"), sep="\t", index=False)
    summary["borzoi_targets_total"] = int(len(tgt))
    summary["borzoi_targets_by_priority"] = {
        str(k): int(v) for k, v in tgt.priority.value_counts().sort_index().items()}
    log(f"Borzoi targets: {len(tgt)}  by priority {summary['borzoi_targets_by_priority']}")

    # ---------------- manifest ------------------------------------------------
    man = [dict(path=p, bytes=os.path.getsize(p), sha256=sha256(p))
           for p in (ASE_281367, ASE_244832, CBP_ADULT, DIRFEAT,
                     os.path.join(CBP_FAMILIES["hepatocyte_5fold_v2"], "fold_manifest.tsv"),
                     os.path.join(CBP_FAMILIES["gse281367_hepatocyte_5fold_v2"],
                                  "fold_manifest.tsv"))]
    pd.DataFrame(man).to_csv(os.path.join(sc, "MANIFEST_step0.tsv"), sep="\t", index=False)
    with open(os.path.join(sc, "build_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
