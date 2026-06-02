#!/bin/bash
#SBATCH --job-name=gwas_var2gene
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=90:00:00
#SBATCH --output=logs/gwas_var2gene_%j.out
#SBATCH --error=logs/gwas_var2gene_%j.err

# Phase 1b: variant -> nearest-TSS gene mapping for traits 12a (100kb) + 12b (500kb).
#
# 1. Aggregate combined_finemapping.csv variants by variant_id (max recommended_pip)
# 2. Build variants.bed (chr, pos-1, pos, variant_id, max_pip)
# 3. Build gene_tss.bed from GENCODE v49 GTF (TSS = start for + strand, end for -)
# 4. bedtools closest -d -t first  -> variant_id, gene, distance_bp
# 5. For 100kb: filter dist <= 100000; aggregate gene-level max(pip); filter > 0.1
# 6. For 500kb: same with 500000
# 7. Append two .gs traits (variant-level 100kb / 500kb) to _phase1b_traits.gs
# 8. Update GWAS_anchored_scdrs_summary.csv with the two new rows.
#
# Final merge happens via:
#   Rscript Analysis/SingleCell/scripts/400_mcp/build_gwas_gs.R --merge-only
# which concatenates _phase1_traits.gs + _phase1b_traits.gs -> GWAS_anchored_scdrs.gs

set -eo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$PROJECT_ROOT"
mkdir -p logs

module load BEDTools/2.31.0-GCC-12.3.0

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

GTF="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
FINEMAP="GWAS/finemapping/results/combined_finemapping.csv"
OUT_DIR="Analysis/SingleCell/results_gpu_v2/disease_signatures"
WORK_DIR="${OUT_DIR}/_phase1b_workdir"
mkdir -p "$WORK_DIR"

PHASE1B_GS="${OUT_DIR}/_phase1b_traits.gs"
SUMMARY_CSV="${OUT_DIR}/GWAS_anchored_scdrs_summary.csv"

VARIANTS_BED="${WORK_DIR}/variants.bed"
VARIANTS_BED_SORTED="${WORK_DIR}/variants.sorted.bed"
GENE_TSS_BED="${WORK_DIR}/gene_tss.bed"
GENE_TSS_BED_SORTED="${WORK_DIR}/gene_tss.sorted.bed"
CLOSEST_TSV="${WORK_DIR}/variant_to_gene.tsv"
GENE_PIP_TSV="${WORK_DIR}/gene_pip.tsv"

echo "[$(date)] Step 1: aggregate variant_id -> max(recommended_pip) and emit BED"

# Use awk to:
#   - skip header
#   - skip rows where variant_id is empty or "NA:NA:NA:NA"
#   - skip rows where recommended_pip is empty/NA
#   - aggregate max pip per variant_id
#   - emit BED rows
# combined_finemapping.csv col 20 = variant_id, col 21 = recommended_pip
awk -F',' '
  BEGIN {OFS="\t"}
  NR==1 {
    for (i=1; i<=NF; i++) {
      if ($i=="variant_id") vid=i
      if ($i=="recommended_pip") pip=i
    }
    if (!vid || !pip) {print "ERROR: header missing variant_id or recommended_pip" > "/dev/stderr"; exit 1}
    next
  }
  {
    v=$vid; p=$pip
    if (v=="" || v=="NA:NA:NA:NA") next
    if (p=="" || p=="NA") next
    if (!(v in best) || p+0 > best[v]) best[v] = p+0
  }
  END {
    n_kept=0; n_dropped=0
    for (v in best) {
      n=split(v, a, ":")
      if (n != 4) {n_dropped++; continue}
      chr=a[1]; pos=a[2]
      if (chr !~ /^[0-9XY]+$/) {n_dropped++; continue}
      if (pos !~ /^[0-9]+$/)   {n_dropped++; continue}
      # Add chr prefix to match GENCODE GTF chromosome naming
      print "chr"chr, pos-1, pos, v, best[v]
      n_kept++
    }
    printf("aggregated variants: kept=%d dropped=%d\n", n_kept, n_dropped) > "/dev/stderr"
  }
' "$FINEMAP" > "$VARIANTS_BED"

echo "[$(date)] variants.bed lines: $(wc -l < "$VARIANTS_BED")"

echo "[$(date)] Step 2: build gene TSS BED from GTF"
# Extract TSS per gene: + strand -> start (col 4); - strand -> end (col 5).
# Output: chr, tss-1, tss, gene_name (HGNC symbol; falls back to gene_id when name missing).
zcat "$GTF" \
  | awk -F'\t' '
      BEGIN {OFS="\t"}
      $0 ~ /^#/ { next }
      $3 == "gene" {
        chr=$1; start=$4; end=$5; strand=$7; attr=$9
        # Skip alt/patch/scaffold contigs (only keep chr1-22, chrX, chrY)
        if (chr !~ /^chr([0-9]+|X|Y)$/) next
        if (strand == "+") tss=start
        else if (strand == "-") tss=end
        else next
        # Parse gene_name from attribute string
        gname=""
        n=split(attr, parts, ";")
        for (i=1; i<=n; i++) {
          if (match(parts[i], /gene_name "[^"]+"/)) {
            gname=substr(parts[i], RSTART+11, RLENGTH-12)
            break
          }
        }
        if (gname == "") next
        print chr, tss-1, tss, gname
      }
    ' > "$GENE_TSS_BED"

echo "[$(date)] gene_tss.bed lines: $(wc -l < "$GENE_TSS_BED")"

echo "[$(date)] Step 3: sort BED files"
sort -k1,1 -k2,2n "$VARIANTS_BED" > "$VARIANTS_BED_SORTED"
sort -k1,1 -k2,2n "$GENE_TSS_BED" > "$GENE_TSS_BED_SORTED"

echo "[$(date)] Step 4: bedtools closest"
# -t first: pick first tie among equally-close TSS (deterministic)
# -d: report distance in bp
# Output cols: <variants.bed cols (5)> <gene_tss.bed cols (4)> <distance>
bedtools closest \
  -a "$VARIANTS_BED_SORTED" \
  -b "$GENE_TSS_BED_SORTED" \
  -d -t first \
  > "$CLOSEST_TSV"

echo "[$(date)] closest.tsv lines: $(wc -l < "$CLOSEST_TSV")"

echo "[$(date)] Step 5: build per-gene PIPs at 100kb / 500kb proximity"
# Input cols: 1=chr 2=v_start 3=v_end 4=variant_id 5=pip 6=tss_chr 7=tss_start
#             8=tss_end 9=gene_name 10=distance
awk -F'\t' '
  BEGIN {OFS="\t"; print "gene","pip_100kb","pip_500kb"}
  {
    pip=$5+0; gene=$9; d=$10+0
    if (gene == "" || gene == ".") next
    # 500kb gene-level
    if (d <= 500000) {
      if (!(gene in g500) || pip > g500[gene]) g500[gene] = pip
    }
    if (d <= 100000) {
      if (!(gene in g100) || pip > g100[gene]) g100[gene] = pip
    }
  }
  END {
    for (g in g500) {
      a = (g in g100) ? g100[g] : "NA"
      print g, a, g500[g]
    }
  }
' "$CLOSEST_TSV" > "$GENE_PIP_TSV"

echo "[$(date)] gene_pip.tsv lines: $(wc -l < "$GENE_PIP_TSV")"

echo "[$(date)] Step 6: write Phase 1b .gs traits via R"

Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({ library(data.table) })

PROJECT_ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT_DIR <- file.path(PROJECT_ROOT,
                     "Analysis/SingleCell/results_gpu_v2/disease_signatures")
WORK_DIR    <- file.path(OUT_DIR, "_phase1b_workdir")
GENE_PIP_FP <- file.path(WORK_DIR, "gene_pip.tsv")
PHASE1B_GS  <- file.path(OUT_DIR, "_phase1b_traits.gs")
SUMMARY_CSV <- file.path(OUT_DIR, "GWAS_anchored_scdrs_summary.csv")

GENE_PIP_THRESH <- 0.1
SOURCE_DESC_100KB <- "GWAS/finemapping/results/combined_finemapping.csv + GENCODE v49 TSS [<=100kb]"
SOURCE_DESC_500KB <- "GWAS/finemapping/results/combined_finemapping.csv + GENCODE v49 TSS [<=500kb]"

dt <- fread(GENE_PIP_FP, na.strings = c("NA", ""))
cat("Loaded gene_pip.tsv: ", nrow(dt), " gene rows\n", sep = "")

format_geneset_row <- function(trait, genes, weights) {
  ok <- !is.na(genes) & nzchar(genes) & !is.na(weights) & is.finite(weights) &
        weights > 0
  genes <- genes[ok]; weights <- weights[ok]
  if (!length(genes)) {
    warning(sprintf("[%s] no valid genes", trait))
    return(NULL)
  }
  d <- data.table(g = genes, w = weights)
  d <- d[, .(w = max(w)), by = g]
  setorder(d, -w)
  paste0(trait, "\t",
         paste0(d$g, ":", signif(d$w, 6), collapse = ","))
}

# 100kb trait
t100 <- dt[!is.na(pip_100kb) & pip_100kb > GENE_PIP_THRESH,
           .(gene, pip = pip_100kb)]
t100 <- t100[, .(pip = max(pip)), by = gene]
cat(sprintf("[gwas_finemapping_variants_100kb] n_genes=%d weight=[%.3f, %.3f]\n",
            nrow(t100), min(t100$pip), max(t100$pip)))

# 500kb trait
t500 <- dt[!is.na(pip_500kb) & pip_500kb > GENE_PIP_THRESH,
           .(gene, pip = pip_500kb)]
t500 <- t500[, .(pip = max(pip)), by = gene]
cat(sprintf("[gwas_finemapping_variants_500kb] n_genes=%d weight=[%.3f, %.3f]\n",
            nrow(t500), min(t500$pip), max(t500$pip)))

rows <- list(
  format_geneset_row("gwas_finemapping_variants_100kb", t100$gene, t100$pip),
  format_geneset_row("gwas_finemapping_variants_500kb", t500$gene, t500$pip)
)
rows <- Filter(Negate(is.null), rows)

con <- file(PHASE1B_GS, "w")
writeLines("TRAIT\tGENESET", con)
for (r in rows) writeLines(r, con)
close(con)
cat(sprintf("Wrote %d Phase 1b traits to %s\n", length(rows), PHASE1B_GS))

# Append to summary CSV
new_summary <- data.table(
  trait        = c("gwas_finemapping_variants_100kb",
                   "gwas_finemapping_variants_500kb"),
  n_genes      = c(nrow(t100), nrow(t500)),
  weight_min   = c(min(t100$pip), min(t500$pip)),
  weight_max   = c(max(t100$pip), max(t500$pip)),
  weight_median= c(median(t100$pip), median(t500$pip)),
  source_files = c(SOURCE_DESC_100KB, SOURCE_DESC_500KB)
)
if (file.exists(SUMMARY_CSV)) {
  prev <- fread(SUMMARY_CSV)
  # drop any existing rows for these two trait names (idempotent)
  prev <- prev[!trait %in% new_summary$trait]
  combined <- rbind(prev, new_summary, use.names = TRUE)
} else {
  combined <- new_summary
}
fwrite(combined, SUMMARY_CSV)
cat(sprintf("Updated summary CSV at %s (now %d rows)\n",
            SUMMARY_CSV, nrow(combined)))
RSCRIPT

echo "[$(date)] Step 7: merge Phase 1 + Phase 1b into final .gs"
Rscript "${PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp/build_gwas_gs.R" --merge-only

echo "[$(date)] DONE."
