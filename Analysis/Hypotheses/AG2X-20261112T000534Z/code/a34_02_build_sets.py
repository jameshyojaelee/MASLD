#!/usr/bin/env python3
"""Build the A3b (sign) and A4 (fine-mapped variant recovery) variant sets from the extracted eQTL PIPs.

Rules (fixed in PRESPEC_AG2X.md):
- Gene start = 5' end of the gene in GENCODE (cellranger-arc GRCh38-2024-A genes.gtf.gz), strand-aware.
- hg38 REF is the GRCh38 FASTA base at the lifted position; the Broadaway NEA/EA pair must contain it.
  REF == NEA -> ALT = EA (effect allele is ALT, sign +1); REF == EA -> ALT = NEA (sign -1); else excluded.
- Eligible variant: SNV, one marginal row at that gene x position, unique liftover to the same chromosome
  on the + strand, REF check passed, gene start within 450,000 bp (inside the 1,048,576-bp scored window).
- Top variant: the gene's maximum PIP, PIP >= 0.9 (primary) or >= 0.5 (sensitivity), unique maximum.
- A3b additionally drops palindromic (A/T, C/G) tops, whose strand cannot be checked.
- A4 controls: same gene fit, PIP < 0.01, eligible, not the top; MAF within 0.05 of the top and
  |log10(d + 1000) - log10(d_top + 1000)| <= 0.2 (d = distance to gene start). Up to K controls, nearest by
  (dMAF/0.05)^2 + (dlogd/0.2)^2, ties broken by a seeded random key (seed 20261111). Genes with 0 controls
  leave A4. Unit: 1-Mb block of the top variant.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

GTF = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/genes/genes.gtf.gz"
IN = L.PKG / "inputs"
MAX_TSS = 450_000
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def gene_starts():
    rows = []
    rx = re.compile(r'gene_id "([^"]+)"')
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] != "gene":
                continue
            gid = rx.search(f[8]).group(1).split(".")[0]
            start = int(f[3]) if f[6] == "+" else int(f[4])
            rows.append((gid, f[0], start, f[6]))
    g = pd.DataFrame(rows, columns=["ensg", "gene_chr", "gene_start", "gene_strand"])
    return g.drop_duplicates("ensg", keep=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True, help="maximum controls per gene")
    ap.add_argument("--dry", action="store_true", help="counts only; write nothing but the count json")
    a = ap.parse_args()
    counts = {"k_max_controls": a.k}
    d = pd.read_csv(IN / "eqtl_pip_candidates_hg38.tsv.gz", sep="\t", low_memory=False)
    counts["rows_in"] = int(len(d))
    gs = gene_starts()
    d = d.merge(gs, on="ensg", how="left")
    d["is_snv"] = d.NEA.astype(str).str.len().eq(1) & d.EA.astype(str).str.len().eq(1)
    ok = (d.is_snv & d.n_rows_at_position.eq(1) & d.liftover_hits.eq(1) & d.strand_hg38.eq("+")
          & d.chr_hg38.eq("chr" + d.chr.astype(str)) & d.gene_chr.eq(d.chr_hg38))
    d["pre_ok"] = ok
    import pyfaidx
    fa = pyfaidx.Fasta(str(L.FASTA), as_raw=True, sequence_always_upper=True)
    ref = []
    for c, p, o in zip(d.chr_hg38, d.pos_hg38, d.pre_ok):
        ref.append(fa[c][int(p) - 1] if o else None)
    d["ref_hg38"] = ref
    nea, ea = d.NEA.astype(str).str.upper(), d.EA.astype(str).str.upper()
    d["ea_is_alt"] = np.where(d.ref_hg38 == nea, 1, np.where(d.ref_hg38 == ea, -1, 0))
    d["alt_hg38"] = np.where(d.ea_is_alt == 1, ea, np.where(d.ea_is_alt == -1, nea, None))
    d["palindromic"] = [COMP.get(x) == y for x, y in zip(nea, ea)]
    d["dist_tss"] = (d.pos_hg38 - d.gene_start).abs()
    d["maf"] = np.minimum(d.EAF, 1 - d.EAF)
    d["eligible"] = d.pre_ok & d.ea_is_alt.ne(0) & d.dist_tss.le(MAX_TSS)
    d["uid"] = np.where(d.eligible, d.chr_hg38.astype(str) + ":" + d.pos_hg38.astype("Int64").astype(str) + ":" + d.ref_hg38.astype(str) + ":" + d.alt_hg38.astype(str), None)
    d["block_1mb"] = np.where(d.pre_ok, d.chr_hg38.astype(str) + "~" + (d.pos_hg38 // 1_000_000).astype("Int64").astype(str), None)
    tops = d[d.is_top.astype(str).str.upper().eq("TRUE")].copy()
    counts["genes_top_ge_0.5"] = int(tops.ensg.nunique())
    tied = tops.ensg.value_counts()
    counts["genes_with_tied_maximum_pip_resolved_by_smallest_p"] = int((tied > 1).sum())
    trng = np.random.default_rng(L.SEED + 1)
    tops = tops.sort_values(["ensg", "pos_hg19"]).copy()
    tops["tiebreak"] = trng.random(len(tops))
    tops = tops.sort_values(["ensg", "PVAL", "tiebreak"]).drop_duplicates("ensg", keep="first").drop(columns="tiebreak")
    counts["tops_not_snv"] = int((~tops.is_snv).sum())
    counts["tops_multirow_position"] = int((tops.is_snv & tops.n_rows_at_position.ne(1)).sum())
    counts["tops_liftover_fail"] = int((tops.is_snv & tops.n_rows_at_position.eq(1) & ~tops.pre_ok).sum())
    counts["tops_ref_mismatch"] = int((tops.pre_ok & tops.ea_is_alt.eq(0)).sum())
    counts["tops_gene_start_beyond_450kb"] = int((tops.pre_ok & tops.ea_is_alt.ne(0) & tops.dist_tss.gt(MAX_TSS)).sum())
    counts["tops_no_gene_start_in_gtf"] = int(tops.gene_start.isna().sum())
    counts["ref_equals_nea_share_among_prechecked_tops"] = float((tops[tops.pre_ok].ea_is_alt == 1).mean())
    el = tops[tops.eligible].copy()
    for thr in (0.9, 0.5):
        s = el[el.pip >= thr]
        counts[f"A3b_genes_pip_ge_{thr}"] = int((~s.palindromic).sum())
        counts[f"A3b_blocks_pip_ge_{thr}"] = int(s[~s.palindromic].block_1mb.nunique())
        counts[f"eligible_tops_pip_ge_{thr}_incl_palindromic"] = int(len(s))

    # ---- A4 controls (PIP >= 0.9 tops only) ----
    rng = np.random.default_rng(L.SEED)
    cand = d[d.eligible & (d.pip < 0.01)].sort_values(["ensg", "pos_hg38"]).copy()
    cand["tiebreak"] = rng.random(len(cand))
    cg = {g: x for g, x in cand.groupby("ensg")}
    pairs = []
    n_ctrl = {}
    for _, t in el[el.pip >= 0.9].iterrows():
        c = cg.get(t.ensg)
        if c is None:
            n_ctrl[t.ensg] = 0
            continue
        c = c[c.uid != t.uid]
        dm = (c.maf - t.maf).abs()
        dl = (np.log10(c.dist_tss + 1000) - np.log10(t.dist_tss + 1000)).abs()
        c = c[(dm <= 0.05) & (dl <= 0.2)].copy()
        c["score"] = (dm[c.index] / 0.05) ** 2 + (dl[c.index] / 0.2) ** 2
        c = c.sort_values(["score", "tiebreak"]).head(a.k)
        n_ctrl[t.ensg] = int(len(c))
        for _, r in c.iterrows():
            pairs.append({"ensg": t.ensg, "gene": t.gene, "top_uid": t.uid, "ctrl_uid": r.uid, "top_pip": t.pip, "ctrl_pip": r.pip,
                          "top_maf": t.maf, "ctrl_maf": r.maf, "top_dist_tss": int(t.dist_tss), "ctrl_dist_tss": int(r.dist_tss),
                          "block_1mb": t.block_1mb})
    nc = pd.Series(n_ctrl)
    counts["A4_genes_pip_ge_0.9"] = int(len(nc))
    counts["A4_genes_with_ge1_control"] = int((nc > 0).sum())
    counts["A4_genes_with_k_controls"] = int((nc == a.k).sum())
    counts["A4_controls_total"] = int(nc.sum())
    pairs = pd.DataFrame(pairs)
    counts["A4_blocks"] = int(pairs.block_1mb.nunique()) if len(pairs) else 0
    counts["A4_distinct_control_uids"] = int(pairs.ctrl_uid.nunique()) if len(pairs) else 0
    if not a.dry:
        cols = ["ensg", "gene", "chr_hg38", "pos_hg38", "ref_hg38", "alt_hg38", "uid", "NEA", "EA", "EAF", "maf", "Beta", "SE",
                "PVAL", "ea_is_alt", "palindromic", "pip", "pos_hg19", "gene_start", "gene_strand", "dist_tss", "block_1mb"]
        el[cols].to_csv(IN / "a3b_a4_tops.tsv", sep="\t", index=False)
        pairs.to_csv(IN / "a4_pairs.tsv", sep="\t", index=False)
        json.dump(counts, open(IN / "a34_set_counts.json", "w"), indent=1)
    print(json.dumps(counts, indent=1))


if __name__ == "__main__":
    main()
