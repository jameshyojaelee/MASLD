# SSX-B (PRESPEC_SSX_B.md; exploratory). Exon-intron split (EISA, Gaidatzis 2015) of the F3-F4 vs F0-F1 change in
# SteatoSITE. dE and dI: separate voomWithQualityWeights fits ~ batch + inferred_sex + bin (S3a helpers), TMM within
# count type on the full matrix of that type. EISA: edgeR QL (legacy = FALSE) on stacked exon and intron columns with
# design donor + exon + exon:F2 + exon:F3F4. This is the identifiable form of the prespecified ~ 0 + donor + type:bin:
# the exon:F3F4 coefficient equals (exon:F3F4 - intron:F3F4), i.e. the stage change in log exon/intron. BH over genes
# passing the filter (>= 10 intronic and >= 10 exonic reads in >= 50% of units).
# Controls: Spearman(dI, dE) > 0 (else stop); 200 permutations of bin within batch x sex for the EISA hit count
# (seed 20261103); intron/exon proxy covariate in dE, dI; dE vs S3a F3-F4 logFC.
# D-B-1: EISA statistic computed with limma on the donor log2 exon/intron ratio; edgeR stacked fit kept as a
# 300-gene agreement check. Usage: Rscript ssx_b_eisa.R --out <new dir>
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_common.R")
SSX <- file.path(ROOT, "Analysis/Hypotheses/SSX-20261102T133212Z")
a <- parse_args(commandArgs(trailingOnly = TRUE), "out", defaults = list(`time-test` = "0"))
out <- new_outdir(a$out); n_test <- as.integer(a$`time-test`)
stopifnot(system(paste("cd", SSX, "&& sha256sum -c --quiet PRESPEC_SSX.sha256")) == 0)
NAMED <- c("GNMT", "MAT1A", "CYP2C19", "COL1A1", "LUM", "IGFBP7", "EPCAM", "TNFSF10")

U <- fread(file.path(SS, "out_s2_qc/units.tsv"))[toupper(as.character(include)) == "TRUE"]
U[, sample_alias := as.character(sample_alias)]
stg <- read_stage(file.path(ROOT, "data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv"))
U <- merge(U, stg, by = "sample_alias")[!is.na(stage_code)]
setorder(U, sample_alias); stopifnot(nrow(U) == 493L)
U <- add_labels(U, U$stage_code)
info <- U[, .(sample_alias, batch, inferred_sex, bin, intron_exon_ratio)]

d <- readRDS(file.path(SS, "out_s2_qc/dge_included.rds"))
E <- (if (inherits(d, "DGEList")) d$counts else as.matrix(d))[, U$sample_alias]; rm(d)
cdir <- file.path(SSX, "out_b_counts")
I <- do.call(cbind, lapply(U$sample_alias, function(s) {
  z <- fread(file.path(cdir, paste0(s, ".intron.counts.gz")), header = FALSE, col.names = c("gene", "n"))
  setNames(z$n, z$gene) }))
colnames(I) <- U$sample_alias
stopifnot(all(rownames(I) %in% rownames(E)))
nf_E <- calcNormFactors(DGEList(E))$samples$norm.factors; lib_E <- colSums(E)
nf_I <- calcNormFactors(DGEList(I))$samples$norm.factors; lib_I <- colSums(I)
G <- rownames(I)
keep <- rowMeans(E[G, ] >= 10) >= 0.5 & rowMeans(I[G, ] >= 10) >= 0.5
K <- G[keep]
ann <- load_ann(); sym <- setNames(ann$gene_name, ann$gene_id_versioned)
log_msg("genes with intronic features:", length(G), "; passing filter:", sum(keep), "; used:", length(K))
fwrite(data.table(sample_alias = U$sample_alias, lib_exon = lib_E, nf_exon = nf_E, lib_intron = lib_I, nf_intron = nf_I,
                  intron_share = lib_I / (lib_I + lib_E)), file.path(out, "library_totals.tsv"), sep = "\t")

# ---- dE, dI (S3a model), and with the intron/exon proxy covariate
de_fit <- function(M, lib, nf, extra = character()) {
  y <- DGEList(M[K, ], lib.size = lib, norm.factors = nf)
  f <- fit_qw(y, design_of(info, "bin", extra))
  gene_tab(f, "binF3F4") }
dE <- de_fit(E, lib_E, nf_E); dI <- de_fit(I, lib_I, nf_I)
dEp <- de_fit(E, lib_E, nf_E, "intron_exon_ratio"); dIp <- de_fit(I, lib_I, nf_I, "intron_exon_ratio")
rho_EI <- cor(dE$logFC, dI$logFC, method = "spearman")
log_msg("Spearman(dI, dE):", rho_EI)

# ---- EISA (D-B-1): donor-level log2 exon/intron ratio, limma-trend; edgeR stacked fit kept as a check on 300 genes
n <- nrow(U)
lratio <- function(g) log2(sweep(E[g, , drop = FALSE] + .5, 2, lib_E * nf_E, "/") * 1e6) -
                      log2(sweep(I[g, , drop = FALSE] + .5, 2, lib_I * nf_I, "/") * 1e6)
RR <- lratio(K)
AE <- rowMeans(log2(sweep(E[K, ] + .5, 2, lib_E * nf_E, "/") * 1e6))   # exon abundance for the eBayes trend
eisa_fit <- function(bin, extra = FALSE) {
  inf <- copy(info); inf$bin <- factor(bin, levels = BINS)
  X <- if (extra) design_of(inf, "bin") else model.matrix(~ bin, as.data.frame(inf))
  f <- lmFit(RR, X); f$Amean <- AE; f <- eBayes(f, trend = TRUE)
  t <- topTable(f, coef = "binF3F4", number = Inf, sort.by = "none")
  data.table(gene = rownames(RR), dEmI = t$logFC, t = t$t, P = t$P.Value, FDR = p.adjust(t$P.Value, "BH")) }
EI <- eisa_fit(U$bin); EIb <- eisa_fit(U$bin, extra = TRUE)
log_msg("EISA (limma log ratio) hits BH < 0.05:", sum(EI$FDR < 0.05))
# edgeR check on the 300-gene subset of the timing test (same seed and sampling call)
set.seed(1); K300 <- sample(G[keep], 300)
Y <- cbind(E[K300, ], I[K300, ])
yS <- DGEList(Y, lib.size = c(lib_E, lib_I), norm.factors = c(nf_E, nf_I))
ex <- rep(c(1, 0), each = n); b2 <- rep(as.character(U$bin), 2)
XS <- cbind(model.matrix(~ 0 + factor(rep(U$sample_alias, 2))), exon = ex, exF2 = ex * (b2 == "F2"), exF3F4 = ex * (b2 == "F3F4"))
t0 <- Sys.time(); fS <- glmQLFit(yS, XS, legacy = FALSE); tS <- glmQLFTest(fS, coef = "exF3F4")$table
t_fit <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
ER <- data.table(gene = rownames(tS), edgeR_logFC = tS$logFC, edgeR_FDR300 = p.adjust(tS$PValue, "BH"))
ER <- merge(ER, EI[, .(gene, dEmI, P)], by = "gene"); ER[, limma_FDR300 := p.adjust(P, "BH")]
fwrite(ER, file.path(out, "edgeR_vs_limma_300.tsv"), sep = "\t")
if (n_test > 0) quit(save = "no")

# ---- permutation null for the EISA hit count (bin permuted within batch x sex)
set.seed(20261103L)
strat <- paste(info$batch, info$inferred_sex)
perms <- lapply(1:200, function(i) { b <- as.character(U$bin)
  for (s in unique(strat)) { j <- which(strat == s); b[j] <- b[j][sample.int(length(j))] }; b })
null_hits <- unlist(mclapply(perms, function(b) {
  sum(eisa_fit(b)$FDR < 0.05) }, mc.cores = NCORES, mc.preschedule = FALSE))

# ---- classes for replicated S3b genes
SC <- fread(file.path(SS, "out_s3b-20261028T035515Z/steatosite_catalog_sidecar.tsv.gz"))
R <- data.table(gene = K, symbol = sym[K])
R <- merge(R, dE[, .(gene, dE = logFC, dE_FDR = FDR)], by = "gene")
R <- merge(R, dI[, .(gene, dI = logFC, dI_FDR = FDR)], by = "gene")
R <- merge(R, EI[, .(gene, dEmI, eisa_FDR = FDR)], by = "gene")
R <- merge(R, EIb[, .(gene, dEmI_batch_sex = dEmI, eisa_FDR_batch_sex = FDR)], by = "gene")
R <- merge(R, dEp[, .(gene, dE_proxy = logFC)], by = "gene"); R <- merge(R, dIp[, .(gene, dI_proxy = logFC)], by = "gene")
R[, class := fifelse(eisa_FDR < 0.05, "post-transcriptional",
             fifelse(sign(dI) == sign(dE) & dI_FDR < 0.05, "transcriptional", "undetermined"))]
SC[, gene := gene_id_versioned]
REP <- SC[steatosite_stage_state == "replicated", .(gene, symbol, primary_evidence_class, S3a_F3F4 = steatosite_F3F4_log2FC)]
REP <- merge(REP, R[, !"symbol"], by = "gene", all.x = TRUE)
REP[is.na(class), class := fifelse(gene %in% G, "intronic signal too low", "no intronic feature (overlap rule)")]
CL <- REP[, .N, by = class][order(-N)]
CL_dir <- REP[class == "post-transcriptional", .(n = .N), by = .(dE_dir = sign(dE), dEmI_dir = sign(dEmI))]
chk <- merge(R[, .(gene, dE)], SC[, .(gene, S3a = steatosite_F3F4_log2FC)], by = "gene")
NM <- rbindlist(lapply(NAMED, function(s) { g <- ann$gene_id_versioned[ann$gene_name == s][1]
  r <- REP[gene == g]; if (!nrow(r)) r <- R[gene == g]
  data.table(symbol = s, gene = g, in_intron_gtf = g %in% G, passes_filter = g %in% K,
             class = if (nrow(r) && "class" %in% names(r)) r$class[1] else NA_character_,
             dE = if (g %in% K) R[gene == g, dE] else NA_real_, dI = if (g %in% K) R[gene == g, dI] else NA_real_,
             dEmI = if (g %in% K) R[gene == g, dEmI] else NA_real_, eisa_FDR = if (g %in% K) R[gene == g, eisa_FDR] else NA_real_) }))

fwrite(R, file.path(out, "eisa_all_filtered_genes.tsv.gz"), sep = "\t")
fwrite(REP, file.path(out, "replicated_gene_classes.tsv.gz"), sep = "\t")
fwrite(CL, file.path(out, "class_counts.tsv"), sep = "\t")
fwrite(NM, file.path(out, "named_genes.tsv"), sep = "\t")
fwrite(data.table(perm = 1:200, hits = null_hits), file.path(out, "perm_null_hits.tsv"), sep = "\t")
summ <- c(
  paste("donors:", n, "; genes with intronic features:", length(G), "; passing filter:", length(K)),
  paste("method check Spearman(dI, dE):", sprintf("%.3f", rho_EI), if (rho_EI > 0) "PASS" else "FAIL: not interpretable"),
  paste("proxy covariate: Spearman(dE, dE_proxy)", sprintf("%.3f", cor(R$dE, R$dE_proxy, method = "spearman")),
        "; Spearman(dI, dI_proxy)", sprintf("%.3f", cor(R$dI, R$dI_proxy, method = "spearman"))),
  paste("dE vs S3a F3F4 logFC: Pearson", sprintf("%.4f", cor(chk$dE, chk$S3a)), "; median |diff|", signif(median(abs(chk$dE - chk$S3a)), 3)),
  paste("EISA hits BH < 0.05:", sum(R$eisa_FDR < 0.05), "(exon/intron up", sum(R$eisa_FDR < 0.05 & R$dEmI > 0), "; down",
        sum(R$eisa_FDR < 0.05 & R$dEmI < 0), ") ; permutation null median", median(null_hits), "95th pct",
        quantile(null_hits, 0.95), "max", max(null_hits), "; p", signif((1 + sum(null_hits >= sum(R$eisa_FDR < 0.05))) / 201, 3)),
  paste("beside, EISA with exon:batch + exon:sex: hits", sum(R$eisa_FDR_batch_sex < 0.05),
        "; Spearman of dEmI with primary", sprintf("%.3f", cor(R$dEmI, R$dEmI_batch_sex, method = "spearman"))),
  paste("edgeR check (300 genes): Pearson logFC", sprintf("%.3f", cor(ER$edgeR_logFC, ER$dEmI)), "; calls BH<0.05 within 300: edgeR",
        sum(ER$edgeR_FDR300 < 0.05), "limma", sum(ER$limma_FDR300 < 0.05), "both", sum(ER$edgeR_FDR300 < 0.05 & ER$limma_FDR300 < 0.05),
        "; edgeR fit minutes", round(t_fit, 1)))
writeLines(summ, file.path(out, "summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 200); cat(summ, sep = "\n"); print(CL); print(CL_dir); print(NM, digits = 3)
