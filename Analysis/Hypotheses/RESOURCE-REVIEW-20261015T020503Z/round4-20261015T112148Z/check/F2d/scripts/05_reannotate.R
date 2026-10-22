# Independent position re-annotation of gated credible-set member sites.
# GRCh37 -> GRCh38 liftOver (1:1 kept); promoter = |pos38 - TSS| <= 2000 for any
# GENCODE v49 gene on a canonical chromosome; accessible = pos38 in any of 9
# cell-type peak sets (BED start+1..end); ABC = pos37 in any enhancer with
# isSelfPromoter FALSE (start+1..end).
suppressPackageStartupMessages({
  library(rtracklayer); library(GenomicRanges); library(data.table)
})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT <- file.path(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2d/out")
sites <- fread(file.path(OUT, "gated_member_sites.tsv"))  # chrom, pos (GRCh37)
cat("unique sites", nrow(sites), "\n")
gr37 <- GRanges(paste0("chr", sites$chrom), IRanges(sites$pos, width = 1))
gr37$site <- seq_len(nrow(sites))
chain <- import.chain(file.path(ROOT, "data/broadaway_eqtl/hg19ToHg38.over.chain"))
lo <- unlist(liftOver(gr37, chain))
nmap <- tabulate(lo$site, nbins = nrow(sites))
cat("map counts: 0 =", sum(nmap == 0), " 1 =", sum(nmap == 1), " >1 =", sum(nmap > 1), "\n")
one <- lo[lo$site %in% which(nmap == 1)]
sites$pos38 <- NA_integer_
sites$chr38 <- NA_character_
sites$pos38[one$site] <- start(one)
sites$chr38[one$site] <- as.character(seqnames(one))
cat("1:1 maps to a different chromosome:", sum(!is.na(sites$chr38) & sites$chr38 != paste0("chr", sites$chrom)), "\n")

# promoter
g <- fread(file.path(ROOT, "RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/annotation/GENCODE_v49/run_19749483/gene_identity.tsv"))
g <- g[tolower(as.character(is_canonical_chromosome)) == "true"]
g[, tss := ifelse(strand == "+", start_1based, end_1based)]
tss <- GRanges(g$chromosome, IRanges(g$tss, width = 1))
ok <- !is.na(sites$pos38)
s38 <- GRanges(sites$chr38[ok], IRanges(sites$pos38[ok], width = 1))
win <- resize(s38, width = 4001, fix = "center")  # pos-2000 .. pos+2000 inclusive
sites$promoter <- NA
sites$promoter[ok] <- overlapsAny(win, tss)
# accessible
bfiles <- list.files(file.path(ROOT, "Analysis/ATAC/Human_Multiome/results/label_transfer/cell_type_peak_sets_v2"), pattern = "_peaks.bed$", full.names = TRUE)
cat("peak files", length(bfiles), "\n")
pk <- rbindlist(lapply(bfiles, function(f) fread(f, header = FALSE, select = 1:3)))
pkr <- GRanges(pk$V1, IRanges(pk$V2 + 1L, pk$V3))
sites$accessible <- NA
sites$accessible[ok] <- overlapsAny(s38, pkr)
# ABC (GRCh37)
abc <- fread(file.path(ROOT, "data/external/abc_liver/abc_liver_hepg2_enhancers.tsv.gz"))
abc <- abc[toupper(as.character(isSelfPromoter)) == "FALSE"]
abr <- GRanges(abc$chr, IRanges(abc$start + 1L, abc$end))
sites$abc <- overlapsAny(gr37, abr)
fwrite(sites, file.path(OUT, "sites_reannotated_check.tsv"), sep = "\t")
cat("promoter TRUE", sum(sites$promoter, na.rm = TRUE), " accessible TRUE", sum(sites$accessible, na.rm = TRUE), " abc TRUE", sum(sites$abc), "\n")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_05.txt"))
