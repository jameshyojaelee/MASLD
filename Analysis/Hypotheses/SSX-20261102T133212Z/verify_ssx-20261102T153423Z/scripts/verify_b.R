# Independent recompute of SSX-B EISA (limma on donor log2 exon/intron ratio) and replicated-gene classes.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
wd <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z"
out <- file.path(wd, "verify_ssx-20261102T153423Z")
ss <- file.path(wd, "../SS-20261025T135237Z")
set.seed(1)
samples <- readLines(file.path(wd, "scripts/ssx_b_samples.txt"))
cat("samples:", length(samples), "\n")

# intronic matrix
il <- lapply(samples, function(s) fread(file.path(wd, "out_b_counts", paste0(s, ".intron.counts.gz")), header = FALSE,
                                         col.names = c("gene", "count")))
genes_i <- il[[1]]$gene
stopifnot(all(vapply(il, function(z) identical(z$gene, genes_i), TRUE)))
I <- sapply(il, `[[`, "count"); rownames(I) <- genes_i; colnames(I) <- samples
cat("intronic features:", nrow(I), "\n")

# exonic matrix
dge <- readRDS(file.path(ss, "out_s2_qc/dge_included.rds"))
E515 <- dge$counts
E <- E515[, samples]

# stage
st <- fread(cmd = "tr -d '\\r' < /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv")
stg <- setNames(as.integer(st$nash_crn_kleiner_fibrosis_stage), st$sample_alias)[samples]
stopifnot(!anyNA(stg), all(stg %in% 1:4))
bin <- factor(c("F0F1", "F2", "F3F4", "F3F4")[stg], levels = c("F0F1", "F2", "F3F4"))
print(table(bin))
smp <- dge$samples[samples, ]
batch <- factor(smp$batch); sex <- factor(smp$inferred_sex)

# TMM within count type on the full matrix of that type (493 units)
libE <- colSums(E); nfE <- calcNormFactors(DGEList(E))$samples$norm.factors
libI <- colSums(I); nfI <- calcNormFactors(DGEList(I))$samples$norm.factors
# alternative: exon TMM on all 515 included libraries
nfE515 <- setNames(calcNormFactors(DGEList(E515))$samples$norm.factors, colnames(E515))[samples]
fwrite(data.table(sample_alias = samples, lib_exon = libE, nf_exon = nfE, nf_exon_515 = nfE515,
                  lib_intron = libI, nf_intron = nfI), file.path(out, "b_library_totals.tsv"), sep = "\t")

# gene filter
common <- intersect(rownames(I), rownames(E))
cat("intronic features also in exon matrix:", length(common), "\n")
Ic <- I[common, ]; Ec <- E[common, ]
half <- 0.5 * length(samples)
keep_joint <- rowSums(Ic >= 10 & Ec >= 10) >= half
keep_sep <- rowSums(Ic >= 10) >= half & rowSums(Ec >= 10) >= half
cat("filter joint (both >=10 in same unit, >=50% units):", sum(keep_joint),
    "; separate (each >=10 in >=50% units):", sum(keep_sep), "\n")

run_eisa <- function(keep, nfE_use, label) {
  g <- common[keep]
  le <- log2(t(t(Ec[g, ] + 0.5) / (libE * nfE_use)) * 1e6)
  li <- log2(t(t(Ic[g, ] + 0.5) / (libI * nfI)) * 1e6)
  r <- le - li
  fit <- eBayes(lmFit(r, model.matrix(~ bin)), trend = TRUE)
  # set Amean to mean exon log2 CPM before eBayes trend
  fit0 <- lmFit(r, model.matrix(~ bin)); fit0$Amean <- rowMeans(le)
  fit <- eBayes(fit0, trend = TRUE)
  tt <- topTable(fit, coef = "binF3F4", number = Inf, sort.by = "none")
  tt$q <- p.adjust(tt$P.Value, "BH")
  cat(label, ": genes", length(g), "; hits BH<0.05", sum(tt$q < 0.05),
      "(pos", sum(tt$q < 0.05 & tt$logFC > 0), "; neg", sum(tt$q < 0.05 & tt$logFC < 0), ")\n")
  data.table(gene = g, dEmI = tt$logFC, P = tt$P.Value, q = tt$q)
}
e_joint <- run_eisa(keep_joint, nfE, "EISA joint filter, exon TMM 493")
e_sep <- if (!identical(keep_joint, keep_sep)) run_eisa(keep_sep, nfE, "EISA separate filter, exon TMM 493") else NULL
e_515 <- run_eisa(keep_joint, nfE515, "EISA joint filter, exon TMM 515")
fwrite(e_joint, file.path(out, "b_eisa_joint.tsv.gz"), sep = "\t")

# dE, dI: voomWithQualityWeights ~ batch + sex + bin, TMM factors within count type (full-matrix factors)
g <- e_joint$gene
des <- model.matrix(~ batch + sex + bin)
fit_vq <- function(M, lib, nf) {
  d <- DGEList(M[g, ], lib.size = lib, norm.factors = nf)
  v <- voomWithQualityWeights(d, des)
  tt <- topTable(eBayes(lmFit(v, des)), coef = "binF3F4", number = Inf, sort.by = "none")
  data.table(gene = g, lfc = tt$logFC, q = p.adjust(tt$P.Value, "BH"))
}
dE <- fit_vq(Ec, libE, nfE); dI <- fit_vq(Ic, libI, nfI)
cat("Spearman(dI, dE):", cor(dE$lfc, dI$lfc, method = "spearman"), "\n")

# classes for replicated genes
sc <- fread(file.path(ss, "out_s3b-20261028T035515Z/steatosite_catalog_sidecar.tsv.gz"))
rep <- sc[steatosite_stage_state == "replicated", gene_id_versioned]
cat("replicated genes:", length(rep), "\n")
cl <- data.table(gene = rep)
cl[, class := ifelse(!gene %in% rownames(I), "no intronic feature", ifelse(!gene %in% g, "intronic signal too low", NA_character_))]
m <- merge(merge(e_joint[, .(gene, eq = q)], dE[, .(gene, dE = lfc, dEq = q)], by = "gene"), dI[, .(gene, dI = lfc, dIq = q)], by = "gene")
cl <- merge(cl, m, by = "gene", all.x = TRUE)
cl[is.na(class), class := ifelse(eq < 0.05, "post-transcriptional",
                          ifelse(sign(dI) == sign(dE) & dIq < 0.05, "transcriptional", "undetermined"))]
print(cl[, .N, by = class][order(-N)])
fwrite(cl, file.path(out, "b_replicated_classes_recomputed.tsv.gz"), sep = "\t")

# compare with reported per-gene tables
rep_all <- fread(file.path(wd, "out_b-20261102T144813Z/eisa_all_filtered_genes.tsv.gz"))
cat("reported eisa table columns:", paste(names(rep_all), collapse = ","), "\n")
rc <- fread(file.path(wd, "out_b-20261102T144813Z/replicated_gene_classes.tsv.gz"))
x <- merge(cl[, .(gene, class_v = class, dE_v = dE, dI_v = dI, dEmI_v = NA_real_)], rc, by = "gene", all = TRUE)
x <- merge(x, e_joint[, .(gene, dEmI_v2 = dEmI, eq_v = q)], by = "gene", all.x = TRUE)
norm_class <- function(z) sub(" \\(overlap rule\\)", "", z)
cat("class agreement with reported (replicated genes):", sum(norm_class(x$class) == x$class_v, na.rm = TRUE), "of", nrow(x), "\n")
print(table(reported = norm_class(x$class), recomputed = x$class_v, useNA = "ifany"))
ok <- !is.na(x$dEmI) & !is.na(x$dEmI_v2)
cat("dEmI max rel diff (replicated, filtered):", max(abs(x$dEmI_v2[ok] - x$dEmI[ok]) / abs(x$dEmI[ok])),
    "; eisa_FDR max rel diff:", max(abs(x$eq_v[ok] - x$eisa_FDR[ok]) / x$eisa_FDR[ok]),
    "; dE max abs diff:", max(abs(x$dE_v[ok] - x$dE[ok])), "; dI max abs diff:", max(abs(x$dI_v[ok] - x$dI[ok])),
    "; dE_FDR max rel diff:", max(abs(merge(dE, rc, by = "gene")[, abs(q - dE_FDR) / dE_FDR])),
    "; dI_FDR max rel diff:", max(abs(merge(dI, rc, by = "gene")[, abs(q - dI_FDR) / dI_FDR])), "\n")
# reported-column self-consistency: apply class rules to reported columns
y <- rc[!is.na(eisa_FDR)]
y[, class_rule := ifelse(eisa_FDR < 0.05, "post-transcriptional",
                   ifelse(sign(dI) == sign(dE) & dI_FDR < 0.05, "transcriptional", "undetermined"))]
cat("rules applied to reported columns reproduce reported class:", sum(y$class_rule == y$class), "of", nrow(y), "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_b.txt"))
