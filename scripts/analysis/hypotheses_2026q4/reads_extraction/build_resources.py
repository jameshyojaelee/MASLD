#!/usr/bin/env python3
"""Build the resource tables read by extract_reads_one_library.py (GRCh38, GENCODE v49 chr names).

Writes into RESOURCE_DIR:
  variants_hg38.tsv            the seven named variants, checked against the FASTA and the local dbSNP VCF
  hotspots_foxo1.tsv           FOXO1 codon 21/22 positions derived from the MANE Select CDS
  foxo1_window.tsv             201 bp window around c.65 used for background error
  ig_constant_exons.tsv        exon union per IGH constant gene from the GENCODE GTF
  redi_v3_slim.tsv.gz          REDIportal v3 hg38 columns used here (all rows)
  redi_v3_stats.tsv            site counts by type x nSamples threshold x dbSNP flag
  editing_panel_hg38.tsv.gz    the A-to-I panel (human readable)
  editing_panel_hg38.npz       the same panel as sorted arrays per chromosome

Usage: build_resources.py RESOURCE_DIR REDI_V3_GZ FASTA GTF BCFTOOLS DBSNP_VCF MIN_NSAMPLES [EXISTING_SLIM_TSV_GZ]
FASTA and GTF are the ones the STAR index was built from (alt contigs included).
The REDIportal file is untrusted input: it is only parsed as text here.
"""
import gzip
import os
import subprocess
import sys
from collections import Counter

import numpy as np
import pysam

res_dir, redi_gz, fasta_path, gtf_path, bcftools, dbsnp_vcf, min_ns = sys.argv[1:8]
min_ns = int(min_ns)
existing_slim = sys.argv[8] if len(sys.argv) > 8 else None
os.makedirs(res_dir, exist_ok=True)
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
NC = {"chr1": "NC_000001.11", "chr2": "NC_000002.12", "chr4": "NC_000004.12", "chr13": "NC_000013.11",
      "chr14": "NC_000014.9", "chr19": "NC_000019.10", "chr22": "NC_000022.11"}
fa = pysam.FastaFile(fasta_path)


def out(name):
    p = os.path.join(res_dir, name)
    if os.path.exists(p):
        sys.exit(f"refusing to overwrite {p}")
    return p


# ---------------------------------------------------------------- named variants
# Positions are GRCh38, 1-based, plus strand. ref/alt follow the local dbSNP VCF
# (GCF_000001405.40). rs72613567 is stored at the left-anchored T>TA form.
VARIANTS = [
    # rsid, gene, chrom, pos, ref, alt, kind, note, hypothesis, alt-contig copy (contig, pos) or None
    ("rs72613567", "HSD17B13", "chr4", 87310240, "T", "TA", "insertion",
     "A insertion next to the exon 6 splice donor; gene on minus strand; Ensembl: splice_donor_variant; intronic, few RNA reads", "H3", None),
    ("rs738409", "PNPLA3", "chr22", 43928847, "C", "G", "snv", "p.Ile148Met (alt G = 148M)", "H3", None),
    ("rs58542926", "TM6SF2", "chr19", 19268740, "C", "T", "snv", "p.Glu167Lys (alt T = 167K); gene on minus strand", "H3", None),
    ("rs2642438", "MTARC1", "chr1", 220796686, "A", "G", "snv",
     "ref A encodes Thr165 (Ensembl minor allele A, MAF 0.19); alt G encodes Ala165", "H3", None),
    ("rs1260326", "GCKR", "chr2", 27508073, "T", "C", "snv",
     "ref T encodes Leu446; alt C encodes Pro446 (Ensembl L/P at 446)", "H3", None),
    ("rs17580", "SERPINA1", "chr14", 94380925, "T", "A", "snv",
     "Pi*S; Ensembl canonical p.Glu288Val (E264V in mature numbering); gene on minus strand; copy on alt contig KI270847.1", "H31",
     ("KI270847.1", 1497497)),
    ("rs28929474", "SERPINA1", "chr14", 94378610, "C", "T", "snv",
     "Pi*Z; Ensembl canonical p.Glu366Lys (E342K in mature numbering); gene on minus strand; copy on alt contig KI270847.1", "H31",
     ("KI270847.1", 1495182)),
]
rows = []
for rsid, gene, chrom, pos, ref, alt, kind, note, hyp, altc in VARIANTS:
    seq = fa.fetch(chrom, pos - 1, pos + 1).upper()
    if kind == "snv":
        fasta_ok = seq[0] == ref
    else:
        fasta_ok = seq[0] == ref[0] and seq[1] == "A"
    if altc:
        # the alt-contig copy must carry the same reference base and the same 21 bp context
        a_ctx = fa.fetch(altc[0], altc[1] - 11, altc[1] + 10).upper()
        m_ctx = fa.fetch(chrom, pos - 11, pos + 10).upper()
        if a_ctx != m_ctx:
            sys.exit(f"{rsid}: alt-contig context differs {a_ctx} vs {m_ctx}")
    q = subprocess.run([bcftools, "query", "-r", f"{NC[chrom]}:{pos}-{pos}", "-f", "%ID\t%REF\t%ALT\n", dbsnp_vcf],
                       check=True, capture_output=True, text=True).stdout.strip().splitlines()
    hit = [x.split("\t") for x in q if x.split("\t")[0] == rsid]
    if not hit:
        sys.exit(f"{rsid} not found at {chrom}:{pos} in the local dbSNP VCF")
    d_ref, d_alt = hit[0][1], hit[0][2].split(",")
    vcf_ok = d_ref == ref and alt in d_alt
    if not (fasta_ok and vcf_ok):
        sys.exit(f"{rsid}: FASTA/VCF mismatch fasta={seq} vcf={hit[0]}")
    rows.append([rsid, gene, chrom, pos, ref, alt, kind, hyp, altc[0] if altc else "NA", altc[1] if altc else "NA",
                 "FASTA_ok", f"dbSNP_VCF {d_ref}>{','.join(d_alt)}", note])
with open(out("variants_hg38.tsv"), "w") as fh:
    fh.write("rsid\tgene\tchrom\tpos\tref\talt\tkind\thypothesis\talt_contig\talt_pos\tfasta_check\tdbsnp_local\tnote\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")

# ---------------------------------------------------------------- GTF: FOXO1 CDS and IG exons
IG_GENES = ["IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2"]
foxo_cds, foxo_tid = [], None
ig_exons = {g: set() for g in IG_GENES}


def attr(a, key):
    i = a.find(key + ' "')
    if i < 0:
        return None
    j = a.find('"', i + len(key) + 2)
    return a[i + len(key) + 2:j]


opener = gzip.open if gtf_path.endswith(".gz") else open
with opener(gtf_path, "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        a = f[8]
        gn = attr(a, "gene_name")
        if gn == "FOXO1" and f[2] == "transcript" and "MANE_Select" in a:
            foxo_tid = attr(a, "transcript_id")
        if gn in ig_exons and f[2] == "exon":
            ig_exons[gn].add((f[0], int(f[3]), int(f[4])))
fo_chrom = None
with opener(gtf_path, "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] == "CDS" and f'transcript_id "{foxo_tid}"' in f[8]:
            foxo_cds.append((int(f[3]), int(f[4])))
            fo_chrom, fo_strand = f[0], f[6]
assert fo_strand == "-" and foxo_cds
foxo_cds.sort(reverse=True)  # minus strand: first CDS base is the highest coordinate


def cds_pos_to_genome(c):
    """c is 1-based CDS index; returns the 1-based genomic coordinate (minus strand)."""
    left = c
    for s, e in foxo_cds:
        n = e - s + 1
        if left <= n:
            return e - (left - 1)
        left -= n
    raise ValueError(c)


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


codon21 = revcomp(fa.fetch(fo_chrom, cds_pos_to_genome(63) - 1, cds_pos_to_genome(61)).upper())
codon22 = revcomp(fa.fetch(fo_chrom, cds_pos_to_genome(66) - 1, cds_pos_to_genome(64)).upper())
assert codon21 in {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}, codon21  # Arg
assert codon22 == "TCG", codon22  # Ser; TCG>TGG is the only single-base route to Trp
p62, p65 = cds_pos_to_genome(62), cds_pos_to_genome(65)
assert (p62, p65) == (40666151, 40666148), (p62, p65)
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
hot = [
    ("R21L", fo_chrom, p62, comp["G"], comp["T"], "c.62G>T", "p.Arg21Leu", foxo_tid),  # plus strand C>A
    ("S22W", fo_chrom, p65, comp["C"], comp["G"], "c.65C>G", "p.Ser22Trp", foxo_tid),  # plus strand G>C
    ("S22X", fo_chrom, p65, comp["C"], comp["A"], "c.65C>A", "p.Ser22Ter", foxo_tid),  # plus strand G>T
]
with open(out("hotspots_foxo1.tsv"), "w") as fh:
    fh.write("name\tchrom\tpos\tref_plus\talt_plus\tcdna\tprotein\ttranscript\tcodon21\tcodon22\tsource\n")
    for h in hot:
        assert fa.fetch(h[1], h[2] - 1, h[2]).upper() == h[3], h
        fh.write("\t".join(map(str, h)) + f"\t{codon21}\t{codon22}\tNg 2021 (PMID 34646017) amino acid changes; coordinates from MANE CDS, checked with Ensembl VEP\n")
wstart, wend = p65 - 100, p65 + 100
with open(out("foxo1_window.tsv"), "w") as fh:
    fh.write("chrom\tpos\tref_plus\n")
    seq = fa.fetch(fo_chrom, wstart - 1, wend).upper()
    for i, b in enumerate(seq):
        fh.write(f"{fo_chrom}\t{wstart + i}\t{b}\n")

with open(out("ig_constant_exons.tsv"), "w") as fh:
    fh.write("gene\tchrom\tstart\tend\n")
    for g in IG_GENES:
        merged = []
        for chrom, s_, e in sorted(ig_exons[g]):
            if merged and merged[-1][0] == chrom and s_ <= merged[-1][2] + 1:
                merged[-1][2] = max(merged[-1][2], e)
            else:
                merged.append([chrom, s_, e])
        assert {m[0] for m in merged} == {"chr14", "KI270846.1"}, (g, {m[0] for m in merged})
        for chrom, s_, e in merged:
            fh.write(f"{g}\t{chrom}\t{s_}\t{e}\n")

# ---------------------------------------------------------------- REDIportal v3
want = ["Region", "Position", "Ref", "Ed", "Strand", "db", "type", "dbsnp", "repeat", "Gene.refGene", "nSamples", "nTissues"]
if existing_slim:
    slim_path = existing_slim
else:
    slim_path = out("redi_v3_slim.tsv.gz")
    n_rows = 0
    with gzip.open(redi_gz, "rt") as fh, gzip.open(slim_path, "wt") as outfh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = [head.index(w) for w in want]
        top = max(idx) + 1
        outfh.write("\t".join(want) + "\n")
        for line in fh:
            f = line.split("\t", top)
            outfh.write("\t".join(f[i] for i in idx) + "\n")
            n_rows += 1
    print("REDIportal v3 rows:", n_rows)

stats = Counter()
panel = []
with gzip.open(slim_path, "rt") as fh:
    fh.readline()
    for line in fh:
        chrom, pos, ref, ed, strand, db, typ, dbsnp, rep, gene, ns, nt = line.rstrip("\n").split("\t")
        ok_sub = (ref, ed, strand) in {("A", "G", "+"), ("T", "C", "-")}
        if not ok_sub or chrom not in CHROMS:
            stats[("not_A2I_or_nonstandard_chrom", "", "", "")] += 1
            continue
        try:
            nsi = int(ns)
        except ValueError:
            nsi = 0
        snp = "dbsnp" if dbsnp not in ("-", "") else "no_dbsnp"
        for thr in (1, 10, 50, 100, 500, 1000):
            if nsi >= thr:
                stats[(typ, f"nSamples>={thr}", snp, "")] += 1
        if nsi >= min_ns and snp == "no_dbsnp":
            panel.append((chrom, int(pos), ref, ed, strand, typ, rep, gene, nsi, int(nt) if nt.isdigit() else 0))
with open(out("redi_v3_stats.tsv"), "w") as fh:
    fh.write("type\tthreshold\tdbsnp\tn_sites\n")
    for (typ, thr, snp, _), n in sorted(stats.items()):
        fh.write(f"{typ}\t{thr}\t{snp}\t{n}\n")

order = {c: i for i, c in enumerate(CHROMS)}
panel.sort(key=lambda r: (order[r[0]], r[1]))
with gzip.open(out("editing_panel_hg38.tsv.gz"), "wt") as fh:
    fh.write("site\tchrom\tpos\tref_plus\talt_plus\tsite_strand\ttype\trepeat\tgene\tnSamples\tnTissues\n")
    for i, r in enumerate(panel):
        fh.write(f"{i}\t" + "\t".join(map(str, r[:2])) + f"\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t{r[7]}\t{r[8]}\t{r[9]}\n")
arrays = {}
sites = np.arange(len(panel), dtype=np.int32)
chrom_arr = np.array([r[0] for r in panel])
pos_arr = np.array([r[1] for r in panel], dtype=np.int32)
ref_code = np.array([{"A": 0, "T": 3}[r[2]] for r in panel], dtype=np.int8)
alt_code = np.array([{"G": 2, "C": 1}[r[3]] for r in panel], dtype=np.int8)
strand_arr = np.array([1 if r[4] == "+" else -1 for r in panel], dtype=np.int8)
type_arr = np.array([0 if r[5] == "ALU" else 1 for r in panel], dtype=np.int8)
for c in CHROMS:
    m = chrom_arr == c
    arrays[f"{c}_pos"] = pos_arr[m]
    arrays[f"{c}_site"] = sites[m]
arrays["ref_code"], arrays["alt_code"], arrays["strand"], arrays["is_nonalu"] = ref_code, alt_code, strand_arr, type_arr
np.savez(out("editing_panel_hg38.npz"), **arrays)
print(f"panel: {len(panel)} sites (nSamples>={min_ns}, not in dbSNP); ALU {int((type_arr == 0).sum())}, other {int((type_arr == 1).sum())}")
