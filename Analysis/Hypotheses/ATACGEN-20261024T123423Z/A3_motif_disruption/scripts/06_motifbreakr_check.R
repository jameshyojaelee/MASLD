#!/usr/bin/env Rscript
# A3 step 6: independent re-scoring of the Q4 top examples with motifbreakR
# (MotifDb jaspar2024 matrices for the same motif IDs, method "log", p < 5e-5).

suppressPackageStartupMessages({
  library(motifbreakR)
  library(MotifDb)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(data.table)
  library(BiocParallel)
})
set.seed(42)
A3 <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A3_motif_disruption"
RES <- file.path(A3, "results")
WORK <- file.path(A3, "work")

top <- fread(file.path(RES, "q4_top_examples.tsv"))
cl <- fread(file.path(RES, "motif_clusters.tsv"))
bg <- fread(file.path(WORK, "bg.tsv"))
bkg <- setNames(bg$freq, bg$base)

top[, snp_name := paste0("chr", chr, ":", pos38, ":", ref, ":", alt)]
bed <- top[, .(chrom = paste0("chr", chr), start = pos38 - 1L, end = pos38, name = snp_name, score = 0, strand = "+")]
bed_file <- file.path(WORK, "q4_top_examples.bed")
fwrite(unique(bed), bed_file, sep = "\t", col.names = FALSE)
snps <- snps.from.file(bed_file, search.genome = BSgenome.Hsapiens.UCSC.hg38, format = "bed")

ids <- unique(cl[cluster_id %in% top$cluster_id, motif_id])
pwm <- subset(MotifDb, dataSource == "jaspar2024" & providerId %in% ids)
pwm <- pwm[!duplicated(mcols(pwm)$providerId)]
message("motifs used: ", length(pwm), " of ", length(ids))

res <- motifbreakR(snpList = snps, pwmList = pwm, threshold = 5e-5, filterp = TRUE,
                   method = "log", bkg = bkg, show.neutral = TRUE, BPPARAM = SerialParam())
out <- as.data.table(as.data.frame(res, row.names = NULL))
out <- out[, .(SNP_id, REF, ALT, providerId, geneSymbol, scoreRef, scoreAlt, alleleDiff, effect)]
fwrite(out, file.path(RES, "q4_motifbreakr_raw.tsv"), sep = "\t")

# agreement on the example's best motif: sign of alleleDiff vs our delta
cmp <- merge(top[, .(snp_name, best_motif, best_motif_name, delta_alt_minus_ref, cluster_id)],
             out, by.x = c("snp_name", "best_motif"), by.y = c("SNP_id", "providerId"), all.x = TRUE)
cmp[, same_sign := sign(alleleDiff) == sign(delta_alt_minus_ref)]
fwrite(cmp, file.path(RES, "q4_motifbreakr_check.tsv"), sep = "\t")
message("best-motif rows found by motifbreakR: ", sum(!is.na(cmp$alleleDiff)), " of ", nrow(cmp),
        "; same sign: ", sum(cmp$same_sign, na.rm = TRUE))
writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_motifbreakr_env.txt"))
