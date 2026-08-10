#!/usr/bin/env Rscript
# Is "shipped values are accurate, only definiteness is broken" true at HIGH LD?
# The 40-random-pair check is dominated by low-LD pairs. If the plink --r
# estimator diverges where |r| is large, then previously-PASSING genes have
# wrong LD and the coverage-loss-only framing is wrong.
suppressPackageStartupMessages({ library(data.table) })
P <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
POP <- "afr"; CHR <- 8L; BLK <- "9154694.9640787"
g <- file.path(P,"data/ld_ref",paste0("1kg_",POP,"_gram"),paste0("chr",CHR),BLK,paste0(BLK,".ld"))
s <- file.path(P,"data/ld_ref",paste0("1kg_",POP),paste0("chr",CHR),BLK,paste0(BLK,".ld"))
G <- as.matrix(fread(g, header=FALSE)); S <- as.matrix(fread(s, header=FALSE))
cat(sprintf("matrices %dx%d\n", nrow(G), ncol(G)))
set.seed(7); idx <- sort(sample(nrow(G), 700))   # 700x700 submatrix = 244k pairs
Gs <- G[idx,idx]; Ss <- S[idx,idx]
ut <- upper.tri(Gs)
gv <- Gs[ut]; sv <- Ss[ut]; d <- abs(gv - sv)
ok <- is.finite(gv) & is.finite(sv)
gv<-gv[ok]; sv<-sv[ok]; d<-d[ok]
cat(sprintf("off-diagonal pairs compared: %d\n\n", length(d)))
brk <- c(-1e-9,0.05,0.2,0.4,0.6,0.8,0.95,1.0001)
bin <- cut(abs(gv), brk, include.lowest=TRUE)
cat(sprintf("%-14s %9s %10s %10s %10s\n","|r| bin (gram)","n","median|d|","mean|d|","max|d|"))
for (b in levels(bin)) { m <- bin==b; if(!any(m)) next
  cat(sprintf("%-14s %9d %10.5f %10.5f %10.5f\n", b, sum(m), median(d[m]), mean(d[m]), max(d[m]))) }
cat(sprintf("\nOVERALL  cor(shipped,rebuilt) = %.8f   max|d| = %.5f\n", cor(gv,sv), max(d)))
cat(sprintf("pairs with |d| > 0.01 : %d (%.4f%%)\n", sum(d>0.01), 100*mean(d>0.01)))
cat(sprintf("pairs with |d| > 0.05 : %d (%.4f%%)\n", sum(d>0.05), 100*mean(d>0.05)))
ev <- function(M){ e <- eigen((M+t(M))/2, symmetric=TRUE, only.values=TRUE)$values; range(e) }
cat(sprintf("\n700x700 submatrix eigenvalue range  shipped [%.6f, %.4f]   rebuilt [%.6f, %.4f]\n",
            ev(Ss)[1], ev(Ss)[2], ev(Gs)[1], ev(Gs)[2]))
