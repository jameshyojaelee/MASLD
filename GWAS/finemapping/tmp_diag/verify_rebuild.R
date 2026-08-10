#!/usr/bin/env Rscript
# Is the rebuilt LD silently wrong?
#
# The consumer (02_run_susie_locus.R) addresses the matrix by BIM ROW NUMBER
# (ld_index). If the rebuilt matrix's row/column order does not match the block
# .bim, every variant is mis-mapped while the matrix still looks perfectly PSD
# and well-conditioned -- the worst possible failure because nothing downstream
# would notice.
#
# Test: recompute specific (i,j) correlations DIRECTLY from genotypes and compare
# to rebuilt[i,j]. Independent of the rebuild code path.
suppressPackageStartupMessages({ library(data.table) })
P <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
PLINK2 <- "/nfs/sw/plink/plink-2.0a5.13/plink"
POP <- "afr"; CHR <- 8L; BLK <- "9154694.9640787"

bim_p <- file.path(P, "data/ld_ref", paste0("1kg_", POP), paste0("chr", CHR), BLK, paste0(BLK, ".bim"))
gram_p <- file.path(P, "data/ld_ref", paste0("1kg_", POP, "_gram"), paste0("chr", CHR), BLK, paste0(BLK, ".ld"))
ship_p <- file.path(P, "data/ld_ref", paste0("1kg_", POP), paste0("chr", CHR), BLK, paste0(BLK, ".ld"))
chrpre <- file.path(P, "data/ld_ref", paste0("1kg_", POP), sprintf("chr%d_%s", CHR, POP))

bim <- fread(bim_p, header = FALSE)
cat(sprintf("block %s: %d variants\n", BLK, nrow(bim)))

# Independent genotype extraction (does NOT reuse the rebuild's code path).
tmp <- tempfile(); writeLines(as.character(bim$V2), paste0(tmp, ".snps"))
st <- system2(PLINK2, c("--bfile", chrpre, "--extract", paste0(tmp, ".snps"),
                        "--export", "A", "--out", tmp), stdout = FALSE, stderr = FALSE)
stopifnot(st == 0L)
raw <- fread(paste0(tmp, ".raw"))
G <- as.matrix(raw[, -(1:6)])
colnames(G) <- sub("_[^_]*$", "", colnames(G))
G <- G[, match(as.character(bim$V2), colnames(G)), drop = FALSE]   # BIM order
for (j in seq_len(ncol(G))) { m <- is.na(G[, j]); if (any(m)) G[m, j] <- mean(G[, j], na.rm = TRUE) }

gram <- as.matrix(fread(gram_p, header = FALSE))
ship <- as.matrix(fread(ship_p, header = FALSE))
cat(sprintf("gram %dx%d | shipped %dx%d | bim %d\n", nrow(gram), ncol(gram), nrow(ship), ncol(ship), nrow(bim)))
stopifnot(nrow(gram) == nrow(bim))

set.seed(1)
pairs <- cbind(sample(nrow(bim), 40), sample(nrow(bim), 40))
pairs <- pairs[pairs[,1] != pairs[,2], , drop = FALSE]
direct <- apply(pairs, 1, function(ij) suppressWarnings(cor(G[, ij[1]], G[, ij[2]])))
fromgram <- gram[pairs]
fromship <- ship[pairs]
ok <- is.finite(direct) & is.finite(fromgram)

cat("\n=== ORDER CHECK: direct genotype cor(i,j) vs rebuilt[i,j] ===\n")
cat(sprintf("  pairs tested        : %d\n", sum(ok)))
cat(sprintf("  max |direct - gram| : %.3e   <- must be ~1e-6 (6-dp rounding)\n", max(abs(direct[ok] - fromgram[ok]))))
cat(sprintf("  correlation         : %.10f\n", cor(direct[ok], fromgram[ok])))
cat(sprintf("  VERDICT: %s\n", if (max(abs(direct[ok] - fromgram[ok])) < 1e-5) "ORDER CORRECT" else "*** ORDER WRONG / MIS-MAPPED ***"))

# A permuted matrix would still be PSD, so prove the test can detect one.
perm <- sample(nrow(gram))
gp <- gram[perm, perm]
cat(sprintf("\n  control (deliberately permuted): max|direct - permuted| = %.3f (test detects mis-mapping: %s)\n",
            max(abs(direct[ok] - gp[pairs][ok])), if (max(abs(direct[ok] - gp[pairs][ok])) > 1e-5) "YES" else "NO"))

cat("\n=== how different are the correlations themselves? ===\n")
okb <- is.finite(fromship) & is.finite(fromgram)
d <- abs(fromship[okb] - fromgram[okb])
cat(sprintf("  |shipped - rebuilt| off-diagonal: median %.4f  mean %.4f  max %.4f  (n=%d)\n",
            median(d), mean(d), max(d), sum(okb)))
cat(sprintf("  correlation between shipped and rebuilt entries: %.6f\n", cor(fromship[okb], fromgram[okb])))
unlink(paste0(tmp, c(".snps", ".raw", ".log")))
