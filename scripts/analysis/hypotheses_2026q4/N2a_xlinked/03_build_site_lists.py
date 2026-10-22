"""Build the N2a site list: common biallelic dbSNP SNVs in exons, chrX (non-PAR, escape/variable genes removed)
and chr19-22 (autosomal technical-noise null). Reads no BAM.
Usage: python -I 03_build_site_lists.py RUN_ROOT
Writes RUN_ROOT/inputs/sites.tsv and RUN_ROOT/inputs/site_build_summary.json (refuses to overwrite)."""
import sys, os, re, gzip, json, bisect, random, subprocess, collections, tempfile

root = sys.argv[1]
SEED = 20260923
AUTO_CAP = 60000
MIN_MAF = 0.05   # TOPMed, else gnomAD genomes, else 1000 Genomes allele frequency from dbSNP FREQ
GTF = "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
VCF = "/gpfs/commons/home/jameslee/reference_genome/dbsnp/GCF_000001405.40.gz"
BCFTOOLS = "/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools"
NC = {"chrX": "NC_000023.11", "chr19": "NC_000019.10", "chr20": "NC_000020.11", "chr21": "NC_000021.9", "chr22": "NC_000022.11"}
NC_INV = {v: k for k, v in NC.items()}
PAR1_END, PAR2_START = 2781479, 155701383   # GRC GRCh38: PAR1 X:10,001-2,781,479; PAR2 X:155,701,383-156,030,895
out_sites = os.path.join(root, "inputs", "sites.tsv")
out_json = os.path.join(root, "inputs", "site_build_summary.json")
for p in (out_sites, out_json):
    if os.path.exists(p):
        sys.exit("refusing to overwrite " + p)

# Tukiainen status by Ensembl stem
status = {}
with open(os.path.join(root, "inputs", "tukiainen_xci_status.tsv")) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        s = r["combined_status"] or ("unknown" if r["reported_t13"] == "Unknown" else "")
        if s:
            status[r["ensg"]] = s

# exons
exons = collections.defaultdict(list)   # chrom -> [(start, end, gene_id, gene_name)]
rx_gid = re.compile(r'gene_id "([^"]+)"'); rx_gn = re.compile(r'gene_name "([^"]+)"')
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line[0] == "#":
            continue
        f = line.split("\t")
        if f[2] != "exon" or f[0] not in NC:
            continue
        gid = rx_gid.search(f[8]).group(1).split(".")[0]
        gn = rx_gn.search(f[8]); gn = gn.group(1) if gn else gid
        exons[f[0]].append((int(f[3]), int(f[4]), gid, gn))
for c in exons:
    exons[c].sort()
starts = {c: [e[0] for e in exons[c]] for c in exons}
maxlen = {c: max(e[1] - e[0] + 1 for e in exons[c]) for c in exons}

# merged region file for bcftools (1-based inclusive, tab-delimited)
tmpdir = tempfile.mkdtemp(prefix="n2a_regions_", dir=os.path.join(root, "logs"))
reg = os.path.join(tmpdir, "exon_regions.txt")
with open(reg, "w") as fh:
    for c in NC:
        cur = None
        for s, e, _, _ in exons[c]:
            if cur and s <= cur[1] + 1:
                cur[1] = max(cur[1], e)
            else:
                if cur:
                    fh.write("%s\t%d\t%d\n" % (NC[c], cur[0], cur[1]))
                cur = [s, e]
        fh.write("%s\t%d\t%d\n" % (NC[c], cur[0], cur[1]))

proc = subprocess.run([BCFTOOLS, "query", "-R", reg, "-i", "INFO/COMMON=1", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/FREQ\n", VCF],
                      capture_output=True, text=True, check=True)
bases = set("ACGT")
rx_freq = {k: re.compile(r"(?:^|\|)" + re.escape(k) + r":([0-9.eE-]+),([0-9.eE-]+)(?:\||$)") for k in ("TOPMED", "GnomAD_genomes", "1000Genomes")}
counts = collections.Counter()
x_rows, a_rows = [], []
for line in proc.stdout.splitlines():
    ch, pos, ref, alt, rsid, freq = line.split("\t")
    counts["vcf_records_in_exon_regions"] += 1
    if len(ref) != 1 or len(alt) != 1 or ref not in bases or alt not in bases:
        counts["not_biallelic_snv"] += 1
        continue
    maf = None
    for study in ("TOPMED", "GnomAD_genomes", "1000Genomes"):
        m = rx_freq[study].search(freq)
        if m:
            maf = min(float(m.group(1)), float(m.group(2)))
            break
    if maf is None:
        counts["no_frequency"] += 1
        continue
    if maf < MIN_MAF:
        counts["maf_below_threshold"] += 1
        continue
    c = NC_INV[ch]; p = int(pos)
    hits = []
    i = bisect.bisect_right(starts[c], p) - 1
    while i >= 0 and starts[c][i] >= p - maxlen[c]:
        s, e, gid, gn = exons[c][i]
        if s <= p <= e:
            hits.append((gid, gn))
        i -= 1
    if not hits:
        counts["no_exon_overlap"] += 1
        continue
    genes = sorted(set(hits))
    if c == "chrX":
        if p <= PAR1_END or p >= PAR2_START:
            counts["x_par"] += 1
            continue
        st = [status.get(g[0], "unclassified") for g in genes]
        if any(g[1] == "XIST" for g in genes):
            counts["x_xist"] += 1
            continue
        if any(s in ("escape", "variable") for s in st):
            counts["x_escape_or_variable"] += 1
            continue
        cls = "X_inactive" if all(s == "inactive" for s in st) else "X_other"
        x_rows.append((c, p, ref, alt, rsid, cls, ",".join(g[1] for g in genes), ",".join(st)))
    else:
        a_rows.append((c, p, ref, alt, rsid, "AUTO", ",".join(g[1] for g in genes), ""))
counts["x_kept"] = len(x_rows); counts["x_kept_inactive_only"] = sum(r[5] == "X_inactive" for r in x_rows)
counts["auto_before_cap"] = len(a_rows)
rng = random.Random(SEED)
if len(a_rows) > AUTO_CAP:
    a_rows = sorted(rng.sample(a_rows, AUTO_CAP), key=lambda r: (r[0], r[1]))
counts["auto_kept"] = len(a_rows)
counts["auto_by_chrom"] = dict(collections.Counter(r[0] for r in a_rows))
x_genes = collections.Counter(r[6].split(",")[0] for r in x_rows)
with open(out_sites, "x") as fh:
    fh.write("chrom\tpos\tref\talt\trsid\tcls\tgenes\tstatus\n")
    for r in x_rows + a_rows:
        fh.write("\t".join(map(str, r)) + "\n")
with open(out_json, "x") as fh:
    json.dump(dict(counts=counts, n_x_genes_with_sites=len(x_genes), seed=SEED, auto_cap=AUTO_CAP, min_maf=MIN_MAF, par1_end=PAR1_END, par2_start=PAR2_START,
                   n_tukiainen_status=len(status)), fh, indent=1)
print(json.dumps(dict(counts), indent=1)); print("X genes with sites:", len(x_genes))
