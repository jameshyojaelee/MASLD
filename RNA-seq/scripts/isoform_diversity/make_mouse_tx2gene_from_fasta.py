#!/usr/bin/env python3
"""
Build the mouse tx2gene from the SAME transcript FASTA that the vM38 kallisto index
was built from (guarantees every index target is covered -> tximport dtuScaledTPM
won't fail), merging gene_biotype / canonical tags from the GTF-derived table where
available. Overwrites tx2gene_vM38_primary.tsv.gz.

FASTA header (GENCODE): >ENSMUST..|ENSMUSG..|OTT..|OTT..|txname|gene_symbol|len|tx_biotype|
"""
import gzip, csv, sys

IDX_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/isoform_diversity/_index"
FASTA   = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/reference/kallisto/gencode.vM38.transcripts.fa.gz"
GTF_T2G = f"{IDX_DIR}/tx2gene_vM38_primary.tsv.gz"          # GTF-derived (has tags)
OUT     = f"{IDX_DIR}/tx2gene_vM38_primary.tsv.gz"          # overwrite

# load GTF-derived tags keyed by txname
gtf = {}
with gzip.open(GTF_T2G, "rt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        gtf[r["txname"]] = r
sys.stderr.write(f"[mouse_tx2gene] GTF-derived rows: {len(gtf)}\n")

n = 0
with gzip.open(FASTA, "rt") as fh, gzip.open(OUT + ".tmp", "wt") as out:
    out.write("txname\tgeneid\tgene_symbol\ttx_biotype\tgene_biotype\t"
              "seqname\tis_mane_select\tis_ensembl_canonical\tappris\n")
    for line in fh:
        if not line.startswith(">"):
            continue
        f = line[1:].rstrip("\n").split("|")
        tx, gid, sym, txbt = f[0], f[1], (f[5] if len(f) > 5 else ""), (f[7] if len(f) > 7 else "")
        g = gtf.get(tx)
        if g:
            gnbt = g["gene_biotype"]; seqn = g["seqname"]
            mane = g["is_mane_select"]; canon = g["is_ensembl_canonical"]; appr = g["appris"]
        else:
            gnbt = txbt; seqn = "NA"; mane = "0"; canon = "0"; appr = ""
        out.write(f"{tx}\t{gid}\t{sym}\t{txbt}\t{gnbt}\t{seqn}\t{mane}\t{canon}\t{appr}\n")
        n += 1

import os
os.replace(OUT + ".tmp", OUT)
sys.stderr.write(f"[mouse_tx2gene] wrote {n} transcripts (from FASTA) -> {OUT}\n")
