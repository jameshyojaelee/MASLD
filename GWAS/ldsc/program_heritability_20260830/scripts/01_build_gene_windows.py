#!/usr/bin/env python3
"""
Turn the 117 liver Hotspot gene programs into genomic windows on GRCh37/hg19.

Each program becomes one genomic annotation: the union of 100 kb windows drawn
around the transcribed region (annotated gene start to gene end) of every gene
that belongs to the program. 100 kb is the window size used for gene-set
annotations in stratified LD score regression (Finucane et al., Nature Genetics
2018), so the numbers here are comparable to published gene-set results.

Two kinds of background annotation are also built, because a program must be
compared against the pool of genes it was actually drawn from, not against the
whole genome:

  background::<lineage>     every gene that entered the Hotspot autocorrelation
                            test in that cell type (the genes that could in
                            principle have joined a program there)
  background::all_lineages  the union of those five gene universes

The program gene lists are read only. Nothing in this script writes to
Analysis/SingleCell/.

Coordinates come from GENCODE v46 mapped back to GRCh37 (gencode.v46lift37),
because the 1000 Genomes Phase 3 reference panel and the baselineLD v2.2 model
are both on GRCh37. Only autosomes are kept: the LD score regression reference
panel has no sex chromosomes.

Outputs (all under gene_sets/):
  annotation_windows.bed        chrom, start, end, annotation  (0-based, merged)
  annotation_manifest.tsv       one row per annotation, with mapping counts
  unmapped_genes.tsv            genes with no usable autosomal GRCh37 locus
"""

import gzip
import os
import sys
from collections import defaultdict

REPO = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
HOTSPOT = os.path.join(REPO, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
OUTDIR = os.path.join(REPO, "GWAS/ldsc/program_heritability_20260830/gene_sets")
GENES = os.path.join(OUTDIR, "gencode_v46lift37_genes_hg19.tsv")

LINEAGES = ["hepatocytes", "fibroblasts", "macrophages", "cholangiocytes", "tcells"]
WINDOW = 100_000
AUTOSOMES = ["chr%d" % c for c in range(1, 23)]
# GRCh37 autosome lengths, used only to report the genome fraction each
# annotation covers.
HG19_LEN = {
    "chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276,
    "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022,
    "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753,
    "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
    "chr21": 48129895, "chr22": 51304566,
}
AUTOSOME_BP = sum(HG19_LEN.values())


def load_gene_coordinates():
    """symbol -> [(chrom, start, end)], ensembl_id -> [(chrom, start, end)]."""
    by_symbol = defaultdict(list)
    by_ensembl = defaultdict(list)
    with open(GENES) as fh:
        for line in fh:
            chrom, start, end, _strand, ensg, symbol, _biotype = line.rstrip("\n").split("\t")
            locus = (chrom, int(start) - 1, int(end))  # GTF is 1-based inclusive
            if symbol:
                by_symbol[symbol].append(locus)
            if ensg:
                by_ensembl[ensg].append(locus)
    return by_symbol, by_ensembl


def read_program_genes():
    """(lineage, module) -> [gene identifiers], preserving the file exactly."""
    programs = defaultdict(list)
    for lineage in LINEAGES:
        path = os.path.join(HOTSPOT, lineage, "module_genes.tsv")
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            assert header[:3] == ["gene", "module", "weight"], (path, header)
            for line in fh:
                gene, module, _weight = line.rstrip("\n").split("\t")
                programs[(lineage, module)].append(gene)
    return programs


def read_gene_universe():
    """lineage -> set of genes that entered the Hotspot autocorrelation test."""
    universe = {}
    for lineage in LINEAGES:
        path = os.path.join(HOTSPOT, lineage, "autocorr.tsv")
        genes = set()
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            assert header[0] == "gene", (path, header)
            for line in fh:
                genes.add(line.split("\t", 1)[0])
        universe[lineage] = genes
    return universe


def windows_for(genes, by_symbol, by_ensembl):
    """Return merged windows and the identifiers that could not be placed."""
    raw = defaultdict(list)
    unmapped = []
    mapped = 0
    for gene in genes:
        loci = by_symbol.get(gene)
        if not loci and gene.startswith("ENSG"):
            loci = by_ensembl.get(gene.split(".")[0])
        if not loci:
            unmapped.append((gene, "no_grch37_locus"))
            continue
        autosomal = [l for l in loci if l[0] in HG19_LEN]
        if not autosomal:
            unmapped.append((gene, "sex_chromosome_only"))
            continue
        mapped += 1
        for chrom, start, end in autosomal:
            lo = max(0, start - WINDOW)
            hi = min(HG19_LEN[chrom], end + WINDOW)
            raw[chrom].append((lo, hi))

    merged = {}
    for chrom, spans in raw.items():
        spans.sort()
        out = []
        for lo, hi in spans:
            if out and lo <= out[-1][1]:
                out[-1][1] = max(out[-1][1], hi)
            else:
                out.append([lo, hi])
        merged[chrom] = out
    return merged, mapped, unmapped


def main():
    by_symbol, by_ensembl = load_gene_coordinates()
    programs = read_program_genes()
    universe = read_gene_universe()

    if len(programs) != 117:
        sys.exit("expected 117 programs, found %d" % len(programs))

    annotations = []  # (name, kind, lineage, module, genes)
    for lineage in LINEAGES:
        modules = sorted(
            (m for (l, m) in programs if l == lineage),
            key=lambda m: int(m) if m.lstrip("-").isdigit() else 10**9,
        )
        for module in modules:
            annotations.append((
                "%s::%s" % (lineage, module), "program", lineage, module,
                programs[(lineage, module)],
            ))
    for lineage in LINEAGES:
        annotations.append((
            "background::%s" % lineage, "background", lineage, "",
            sorted(universe[lineage]),
        ))
    union = sorted(set().union(*universe.values()))
    annotations.append(("background::all_lineages", "background", "", "", union))

    bed_path = os.path.join(OUTDIR, "annotation_windows.bed")
    manifest_path = os.path.join(OUTDIR, "annotation_manifest.tsv")
    unmapped_path = os.path.join(OUTDIR, "unmapped_genes.tsv")

    with open(bed_path, "w") as bed, \
         open(manifest_path, "w") as man, \
         open(unmapped_path, "w") as unm:
        man.write("\t".join([
            "annotation", "file_slug", "kind", "lineage", "module",
            "n_genes_listed", "n_genes_placed", "n_intervals",
            "bp_covered", "frac_autosome",
        ]) + "\n")
        unm.write("annotation\tgene\treason\n")

        for name, kind, lineage, module, genes in annotations:
            merged, mapped, unmapped = windows_for(genes, by_symbol, by_ensembl)
            n_int = 0
            bp = 0
            for chrom in AUTOSOMES:
                for lo, hi in merged.get(chrom, []):
                    bed.write("%s\t%d\t%d\t%s\n" % (chrom, lo, hi, name))
                    n_int += 1
                    bp += hi - lo
            slug = (("%s-M%s" % (lineage, module)) if kind == "program"
                    else ("background-%s" % (lineage or "all_lineages")))
            man.write("\t".join([
                name, slug, kind, lineage, module,
                str(len(genes)), str(mapped), str(n_int),
                str(bp), "%.5f" % (bp / AUTOSOME_BP),
            ]) + "\n")
            for gene, reason in unmapped:
                unm.write("%s\t%s\t%s\n" % (name, gene, reason))

    print("wrote %d annotations (%d programs, %d backgrounds)"
          % (len(annotations), sum(1 for a in annotations if a[1] == "program"),
             sum(1 for a in annotations if a[1] == "background")))
    print("  " + bed_path)
    print("  " + manifest_path)
    print("  " + unmapped_path)


if __name__ == "__main__":
    main()
