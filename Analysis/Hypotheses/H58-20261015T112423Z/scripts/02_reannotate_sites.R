# H58 step 02: re-annotate every unique gated credible-set member site by position.
# Positions only; computes no fraction, difference or estimate. PRESPEC_H58.md
# "Context repair" step 1.
#   promoter: |GRCh38 pos - TSS| <= 2,000 bp inclusive (GENCODE v49 canonical chromosomes;
#             TSS = start_1based on +, end_1based on -). Implemented as a 4,001-bp query
#             window [pos-2000, pos+2000] overlapping a 1-bp TSS (no maxgap).
#   accessible: GRCh38 pos inside any of the 9 cell-type peak BEDs (start+1 .. end).
#   ABC: GRCh37 pos inside a liver/hepatocyte/HepG2 ABC enhancer, isSelfPromoter FALSE.
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
})
args <- commandArgs(trailingOnly = TRUE)
out <- args[which(args == "--out") + 1]
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
base <- file.path(root, "GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321")

mem <- fread(file.path(base, "credible_set_members_annotated.tsv"),
             select = c("credible_set_uid", "trait_scope", "chromosome", "position", "pip_sum_gate_passed"),
             colClasses = "character")
arch <- fread(file.path(base, "credible_set_pip_architecture.tsv"),
              select = c("credible_set_uid", "trait_scope", "pip_sum_gate_passed"), colClasses = "character")
scopes <- c("tier1_direct_masld_pdff", "tier2_liver_enzyme")
mem <- mem[trait_scope %in% scopes]
arch <- arch[trait_scope %in% scopes]
# Gate: same rule as the producer's apply_gate (a set failing in either table is dropped).
bad <- union(arch[tolower(pip_sum_gate_passed) != "true", credible_set_uid],
             mem[tolower(pip_sum_gate_passed) != "true", credible_set_uid])
mem <- mem[!(credible_set_uid %in% bad)]
sites <- unique(mem[, .(chromosome, position = as.integer(position))])
cat("gated member rows:", nrow(mem), " unique sites:", nrow(sites), "\n")

chain <- import.chain(file.path(root, "data/broadaway_eqtl/hg19ToHg38.over.chain"))
gr19 <- GRanges(paste0("chr", sites$chromosome), IRanges(sites$position, width = 1))
lifted <- liftOver(gr19, chain)
nmap <- lengths(lifted)
cat("liftOver map counts: 0 =", sum(nmap == 0), " 1 =", sum(nmap == 1), " >1 =", sum(nmap > 1), "\n")
sites[, n_map := as.integer(nmap)]
sites[, `:=`(chr_hg38 = NA_character_, pos_hg38 = NA_integer_)]
ok <- which(nmap == 1)
g38 <- unlist(lifted[ok])
sites[ok, `:=`(chr_hg38 = as.character(seqnames(g38)), pos_hg38 = start(g38))]

gi <- fread(file.path(root, "RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/annotation/GENCODE_v49/run_19749483/gene_identity.tsv"),
            select = c("chromosome", "start_1based", "end_1based", "strand", "is_canonical_chromosome"))
gi <- gi[tolower(as.character(is_canonical_chromosome)) == "true"]
stopifnot(all(gi$strand %in% c("+", "-")))
gi[, tss := ifelse(strand == "+", start_1based, end_1based)]
tss <- GRanges(gi$chromosome, IRanges(gi$tss, width = 1))
s38 <- sites[!is.na(pos_hg38)]
qwin <- GRanges(s38$chr_hg38, IRanges(s38$pos_hg38 - 2000L, s38$pos_hg38 + 2000L))
q38 <- GRanges(s38$chr_hg38, IRanges(s38$pos_hg38, width = 1))
sites[, promoter_new := NA]
sites[!is.na(pos_hg38), promoter_new := overlapsAny(qwin, tss)]
# boundary self-check on synthetic points: 2,000 bp in, 2,001 bp out
tchk <- GRanges("chrZ", IRanges(c(10000L), width = 1))
qchk <- GRanges("chrZ", IRanges(c(12000L, 12001L) - 2000L, c(12000L, 12001L) + 2000L))
stopifnot(identical(overlapsAny(qchk, tchk), c(TRUE, FALSE)))

pdir <- file.path(root, "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2")
pf <- sort(list.files(pdir, pattern = "_peaks\\.bed$", full.names = TRUE))
stopifnot(length(pf) == 9)
peaks <- rbindlist(lapply(pf, function(f) {
  b <- fread(f, header = FALSE); stopifnot(ncol(b) == 3); setnames(b, c("chr", "start", "end")); b }))
pk <- GRanges(peaks$chr, IRanges(peaks$start + 1L, peaks$end))
sites[, accessible_new := NA]
sites[!is.na(pos_hg38), accessible_new := overlapsAny(q38, pk)]

abc <- fread(file.path(root, "data/external/abc_liver/abc_liver_hepg2_enhancers.tsv.gz"))
stopifnot(all(abc$CellType %in% c("hepatocyte-ENCODE", "HepG2-Roadmap", "liver-ENCODE")))
abc <- abc[toupper(as.character(isSelfPromoter)) == "FALSE"]
ag <- GRanges(abc$chr, IRanges(abc$start + 1L, abc$end))
sites[, abc_new := overlapsAny(gr19, ag)]

fwrite(sites, file.path(out, "site_reannotated.tsv"), sep = "\t", na = "NA")
cat("sites with 1:1 hg38:", sum(sites$n_map == 1), " promoter TRUE:", sum(sites$promoter_new, na.rm = TRUE),
    " accessible TRUE:", sum(sites$accessible_new, na.rm = TRUE), " ABC TRUE:", sum(sites$abc_new), "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_02.txt"))
