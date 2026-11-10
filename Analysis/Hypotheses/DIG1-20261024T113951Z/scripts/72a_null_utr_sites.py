# DIG1 job 72a (post hoc control for job 71b P1; written before any control-site count was read).
# Null sites for reference-allele mapping bias: common SNVs (1kGP high-coverage EUR AF 0.05-0.50, biallelic,
# non-palindromic) in the last exon of the Ensembl-canonical transcript of protein-coding genes whose GTEx v8 liver
# median TPM is 3-30 (PANX1 6.3) and which are tested but not liver eGenes (GTEx liver egenes qval > 0.5). One SNV per
# gene, 400 genes sampled with seed 20261025. Output: BED-like target file for mpileup.
# Deviation (first run, no counts read): GTEx liver qval never exceeds 0.46, so 'not an eGene' is pval_beta > 0.2.
# Usage: python -I 72a_null_utr_sites.py <bed_of_candidates_from_bcftools> <out_targets>  (two-pass; see run_72.sbatch)
import sys, gzip, random, re
mode = sys.argv[1]
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
if mode == "regions":
    tpm = {}
    with gzip.open(f"{ROOT}/data/external/gtex_v8_gene_median_tpm/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", "rt") as f:
        next(f); next(f); h = next(f).rstrip("\n").split("\t"); j = h.index("Liver")
        for l in f:
            r = l.rstrip("\n").split("\t"); tpm[r[0].split(".")[0]] = float(r[j])
    ok = set()
    with gzip.open(f"{ROOT}/data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz", "rt") as f:
        h = next(f).rstrip("\n").split("\t"); q = h.index("pval_beta")
        for l in f:
            r = l.rstrip("\n").split("\t"); g = r[0].split(".")[0]
            if float(r[q]) > 0.2 and 3 <= tpm.get(g, 0) <= 30: ok.add(g)
    last = {}
    for l in open("/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/gencode.v49.primary_chr.gtf"):
        if l.startswith("#"): continue
        r = l.split("\t")
        if r[2] != "exon" or "Ensembl_canonical" not in r[8] or 'gene_type "protein_coding"' not in r[8]: continue
        g = re.search(r'gene_id "([^".]+)', r[8]).group(1)
        if g not in ok: continue
        s, e = int(r[3]), int(r[4])
        if g not in last or (r[6] == "+" and s > last[g][1]) or (r[6] == "-" and s < last[g][1]): last[g] = (r[0], s, e)
    with open(sys.argv[2], "w") as o:
        for g, (c, s, e) in sorted(last.items(), key=lambda x: (x[1][0], x[1][1])): o.write(f"{c}\t{s}\t{e}\t{g}\n")
    print("genes", len(ok), "regions", len(last))
else:  # pick: argv[2] = bcftools query output (CHROM POS REF ALT AF_EUR), argv[3] = regions, argv[4] = out targets
    import bisect
    reg = {}
    for l in open(sys.argv[3]):
        c, s, e, g = l.split(); reg.setdefault(c, []).append((int(s), int(e), g))
    pal = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
    by = {}
    for l in open(sys.argv[2]):
        c, p, ref, alt, af = l.split(); p = int(p)
        if len(ref) != 1 or len(alt) != 1 or "," in alt or (ref, alt) in pal: continue
        try: af = float(af)
        except ValueError: continue
        if not 0.05 <= af <= 0.5: continue
        for s, e, g in reg.get(c, []):
            if s <= p <= e: by.setdefault(g, []).append((c, p, ref, alt, af))
    rnd = random.Random(20261025); gs = sorted(by); rnd.shuffle(gs)
    with open(sys.argv[4], "w") as o:
        for g in gs[:400]:
            c, p, ref, alt, af = rnd.choice(by[g]); o.write(f"{c}\t{p}\t{p}\t{g}\t{ref}\t{alt}\t{af}\n")
    print("genes with SNV", len(by), "written", min(400, len(gs)))
