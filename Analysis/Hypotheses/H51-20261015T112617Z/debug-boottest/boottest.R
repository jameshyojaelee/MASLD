# Debug only: is the GSE213621 bin x sex bootstrap slope centred on the observed slope? Uses SCRAMBLED labels (as stage 2 DEBUG).
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H51-20261015T112617Z/h51_common.R")
D <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H51-20261015T112617Z/debug-20261015T113244Z"
X <- load_inputs(); g <- X$g; S <- readRDS(file.path(D, "stage1/stage1_objects.rds"))
D0 <- fread(file.path(D, "stage1/D0_discovery_grouped_per_gene.tsv.gz"))
set.seed(99L); lab <- as.character(g$bin)
for (s in sort(unique(g$stratum))) { i <- which(g$stratum == s); lab[i] <- lab[i][sample.int(length(i))] }
g[, bin := factor(lab, levels = BINS)]
dgG <- subdge(X$dge_all, g$sample_id, genes = S$KEEP_G)
gs <- S$HCk$F3F4; x <- D0[contrast == "binF3F4"][match(gs, gene), logFC]
set.seed(5L)
idxs <- lapply(1:47, function(k) { cl <- sort(unique(g$cell)); unlist(lapply(cl, function(cc) { i <- which(g$cell == cc); i[sample.int(length(i), length(i), replace = TRUE)] })) })
idxs <- c(list(seq_len(nrow(g))), idxs)
res <- mclapply(seq_along(idxs), function(k) { idx <- idxs[[k]]; d <- dgG[, idx]; colnames(d) <- make.unique(colnames(d))
  f <- fit_qw(d, des_G0(g[idx])); v <- voom(d, des_G0(g[idx])); f2 <- eBayes(lmFit(v, des_G0(g[idx])))
  c(qw = ols_slope(f$coefficients[match(gs, rownames(f$coefficients)), "binF3F4"], x), plain = ols_slope(f2$coefficients[match(gs, rownames(f2$coefficients)), "binF3F4"], x)) }, mc.cores = 16)
r <- do.call(rbind, res); print(r[1, ]); print(colMeans(r[-1, ])); print(apply(r[-1, ], 2, sd))
