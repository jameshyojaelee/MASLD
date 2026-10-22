# Re-annotate every gated credible-set member site with the same three deposited
# contexts the noncoding-context companion uses, keyed on position (allele-independent):
#   promoter: GENCODE v49 TSS within 2 kb of the hg38 position (unique 1:1 liftOver)
#   accessible: hg38 position inside any cell_type_peak_sets_v2 peak (BED -> 1-based)
#   ABC: hg19 position inside a liver/hepatocyte/HepG2 ABC enhancer with isSelfPromoter FALSE
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
})
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out_dir <- file.path(root, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R02-02-coverage")
mem <- fread(file.path(root, "GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv"),
             select = c("chromosome", "position", "pip_sum_gate_passed"), colClasses = "character")
sites <- unique(mem[, .(chromosome, position = as.integer(position))])
cat("unique member sites (all rows):", nrow(sites), "\n")

chain <- import.chain(file.path(root, "data/broadaway_eqtl/hg19ToHg38.over.chain"))
gr19 <- GRanges(paste0("chr", sites$chromosome), IRanges(sites$position, width = 1))
lifted <- liftOver(gr19, chain)
nmap <- lengths(lifted)
cat("liftOver mapping counts (0 / 1 / >1):", sum(nmap == 0), sum(nmap == 1), sum(nmap > 1), "\n")
sites[, `:=`(chr_hg38 = NA_character_, pos_hg38 = NA_integer_)]
ok <- which(nmap == 1)
g38 <- unlist(lifted[ok])
sites[ok, `:=`(chr_hg38 = as.character(seqnames(g38)), pos_hg38 = start(g38))]

# Promoter: GENCODE v49 canonical-chromosome TSS within +/- 2,000 bp
gi <- fread(file.path(root, "RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/annotation/GENCODE_v49/run_19749483/gene_identity.tsv"),
            select = c("chromosome", "start_1based", "end_1based", "strand", "is_canonical_chromosome"))
gi <- gi[tolower(as.character(is_canonical_chromosome)) == "true"]
gi[, tss := ifelse(strand == "+", start_1based, end_1based)]
tss <- GRanges(gi$chromosome, IRanges(gi$tss, width = 1))
s38 <- sites[!is.na(pos_hg38)]
q38 <- GRanges(s38$chr_hg38, IRanges(s38$pos_hg38, width = 1))
prom_hit <- overlapsAny(q38, tss, maxgap = 2000L)
sites[!is.na(pos_hg38), promoter_new := prom_hit]

# Accessibility: union of all cell-type peak sets
pdir <- file.path(root, "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2")
pf <- list.files(pdir, pattern = "_peaks\\.bed$", full.names = TRUE)
peaks <- rbindlist(lapply(pf, function(f) { b <- fread(f, header = FALSE); if (ncol(b) == 4) b <- b[, 2:4]; setnames(b, c("chr", "start", "end")); b }))
pk <- GRanges(peaks$chr, IRanges(peaks$start + 1L, peaks$end))
sites[!is.na(pos_hg38), accessible_new := overlapsAny(q38, pk)]

# ABC (hg19): non-self-promoter liver / hepatocyte / HepG2 enhancers
abc <- fread(file.path(root, "data/external/abc_liver/abc_liver_hepg2_enhancers.tsv.gz"))
abc <- abc[toupper(as.character(isSelfPromoter)) == "FALSE"]
ag <- GRanges(abc$chr, IRanges(abc$start + 1L, abc$end))
sites[, abc_new := overlapsAny(gr19, ag)]

fwrite(sites, file.path(out_dir, "site_context_reannotated.tsv"), sep = "\t")
cat("sites with hg38:", sum(!is.na(sites$pos_hg38)), " promoter:", sum(sites$promoter_new, na.rm = TRUE),
    " accessible:", sum(sites$accessible_new, na.rm = TRUE), " abc:", sum(sites$abc_new), "\n")
