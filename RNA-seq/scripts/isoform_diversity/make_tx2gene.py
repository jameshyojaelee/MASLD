#!/usr/bin/env python3
"""
Parse a GENCODE GTF into a transcript-level annotation table (tx2gene) restricted
to the PRIMARY assembly (chr1-22, X, Y, M), carrying canonical-isoform tags used
downstream for DTU / diversity annotation.

Outputs (to --outdir):
  tx2gene_<tag>.tsv.gz   columns: txname, geneid, gene_symbol, tx_biotype,
                                  gene_biotype, seqname, is_mane_select,
                                  is_ensembl_canonical, appris
  primary_enst_<tag>.txt one versioned ENST per line (the index whitelist)

Versioned IDs are retained (kallisto FASTA headers are versioned).
"""
import argparse, gzip, re, sys

PRIMARY = re.compile(r'^chr([0-9]+|X|Y|M)$')

def attr(field, key):
    m = re.search(key + r' "([^"]+)"', field)
    return m.group(1) if m else ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", required=True, help="label, e.g. v49_primary")
    args = ap.parse_args()

    opn = gzip.open if args.gtf.endswith(".gz") else open
    t2g_path = f"{args.outdir}/tx2gene_{args.tag}.tsv.gz"
    wl_path  = f"{args.outdir}/primary_enst_{args.tag}.txt"

    n_in = n_keep = 0
    with opn(args.gtf, "rt") as fh, \
         gzip.open(t2g_path, "wt") as out, \
         open(wl_path, "w") as wl:
        out.write("txname\tgeneid\tgene_symbol\ttx_biotype\tgene_biotype\t"
                  "seqname\tis_mane_select\tis_ensembl_canonical\tappris\n")
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            n_in += 1
            seqname = f[0]
            if not PRIMARY.match(seqname):
                continue
            a = f[8]
            txname = attr(a, "transcript_id")
            geneid = attr(a, "gene_id")
            sym    = attr(a, "gene_name")
            txbt   = attr(a, "transcript_type")
            gnbt   = attr(a, "gene_type")
            tags   = re.findall(r'tag "([^"]+)"', a)
            is_mane = "1" if "MANE_Select" in tags else "0"
            is_can  = "1" if "Ensembl_canonical" in tags else "0"
            appris  = next((t for t in tags if t.startswith("appris_")), "")
            out.write(f"{txname}\t{geneid}\t{sym}\t{txbt}\t{gnbt}\t{seqname}\t"
                      f"{is_mane}\t{is_can}\t{appris}\n")
            wl.write(txname + "\n")
            n_keep += 1

    sys.stderr.write(f"[make_tx2gene] transcripts seen={n_in} kept(primary)={n_keep}\n")
    sys.stderr.write(f"[make_tx2gene] wrote {t2g_path} and {wl_path}\n")

if __name__ == "__main__":
    main()
