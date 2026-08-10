#!/usr/bin/env Rscript
# Decisive test: is the non-EUR ridge requirement caused by SAMPLE SIZE or by
# how the shipped .ld was ESTIMATED?
#
# Recompute LD for one AFR block from the SAME genotypes, mean-imputing missing
# calls so the result is a genuine Gram matrix (PSD by construction), and compare
# its minimum eigenvalue to the shipped matrix.
#   shipped min eig ~ -0.36  ->  needs lambda 0.3
#   if rebuilt min eig ~ 0   ->  estimation method, and the loci are RECOVERABLE
#   if rebuilt min eig < 0   ->  something else; sample size argument survives
suppressPackageStartupMessages({ library(data.table) })
P <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PREFIX <- file.path(P, "GWAS/finemapping/data/ld_ref/1kg_afr/chr8/9154694.9640787/9154694.9640787")
PLINK2 <- "/nfs/sw/plink/plink-2.0a5.13/plink"
CHRPREF <- file.path(P, "GWAS/finemapping/data/ld_ref/1kg_afr/chr8_afr")

bim <- fread(paste0(PREFIX, ".bim"), header = FALSE)
cat(sprintf("block: %d variants\n", nrow(bim)))

# --- shipped matrix -------------------------------------------------------
Rs <- as.matrix(fread(paste0(PREFIX, ".ld"), header = FALSE))
Rs <- (Rs + t(Rs)) / 2
fin <- rowSums(!is.finite(Rs)) == 0 & colSums(!is.finite(Rs)) == 0
cat(sprintf("shipped: %d x %d, %d non-finite rows dropped\n", nrow(Rs), ncol(Rs), sum(!fin)))
Rs <- Rs[fin, fin, drop = FALSE]; diag(Rs) <- 1
es <- eigen(Rs, symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("SHIPPED  min eig = %+.6f   (rank-def would give ~0)\n", min(es)))

# --- rebuild from genotypes ----------------------------------------------
tmp <- tempfile()
ids <- bim$V2[fin]
writeLines(ids, paste0(tmp, ".snps"))
st <- system2(PLINK2, c("--bfile", CHRPREF, "--extract", paste0(tmp, ".snps"),
                        "--export", "A", "--out", tmp), stdout = FALSE, stderr = FALSE)
stopifnot(st == 0L)
raw <- fread(paste0(tmp, ".raw"))
G <- as.matrix(raw[, -(1:6)])
cat(sprintf("rebuilt from genotypes: %d samples x %d variants; missing calls %.4f%%\n",
            nrow(G), ncol(G), 100 * mean(is.na(G))))
# Mean-impute: makes the covariance a genuine Gram matrix of a complete dataset.
for (j in seq_len(ncol(G))) { m <- is.na(G[, j]); if (any(m)) G[m, j] <- mean(G[, j], na.rm = TRUE) }
keep <- apply(G, 2, sd) > 0
G <- G[, keep, drop = FALSE]
Rr <- cor(G)
er <- eigen(Rr, symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("REBUILT  min eig = %+.6f   (%d monomorphic variants dropped)\n", min(er), sum(!keep)))

lad <- c(1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1)
need <- function(M) { for (l in lad) if (tryCatch({ invisible(chol(M + l * diag(nrow(M)))); TRUE },
                                                  error = function(e) FALSE)) return(l); NA }
cat(sprintf("\nridge needed  SHIPPED = %s   REBUILT = %s\n",
            format(need(Rs)), format(need(Rr))))
unlink(paste0(tmp, c(".snps", ".raw", ".log")))
