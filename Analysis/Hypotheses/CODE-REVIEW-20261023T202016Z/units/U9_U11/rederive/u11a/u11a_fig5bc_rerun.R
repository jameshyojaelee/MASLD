#!/usr/bin/env Rscript
# U11a independent re-derivation of Figure 5B (protein vs mRNA) and 5C (within-MASLD
# protein-histology partial Spearman, BH over 125). Written from METHODS.md
# (Proteomic analyses) and the 5B/5C legends. Deterministic (no resampling here;
# the 5C bootstrap intervals do not enter any call and are not re-derived).
suppressPackageStartupMessages({ library(data.table); library(limma) })
set.seed(20260923L)
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- Sys.getenv("U11A_OUT"); stopifnot(nzchar(out)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
inp <- c(
  raw = file.path(R, "data/PXD051911/liver_protein_quant.txt"),
  meta = file.path(R, "data/PXD051911/meta_data.txt"),
  deg = file.path(R, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"),
  universe = file.path(R, "Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt"),
  fixed25 = file.path(R, "Analysis/Multimodal_Program_Projection/results/proteomics/fixed_panel4c_rows.tsv"),
  prod_de = file.path(R, "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv"),
  prod_5B = file.path(R, "figures/candidates/panel-deep-fixes-20261017T035633Z/fig5/source_tables/5B.tsv"),
  prod_5C = file.path(R, "figures/candidates/panel-deep-fixes-20261017T035633Z/fig5/source_tables/5C.tsv")
)
fwrite(data.table(input = names(inp), path = inp,
                  sha256 = sapply(inp, function(f) sub(" .*", "", system2("sha256sum", f, stdout = TRUE)))),
       file.path(out, "input_sha256.tsv"), sep = "\t")

# ---- protein matrix
raw <- fread(inp[["raw"]]); md <- fread(inp[["meta"]])
cols <- setdiff(names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% cols], by = "liver_proteomics_filename")
stopifnot(nrow(md) == 58L)
md[, `:=`(sid = liver_proteomics_filename,
          masld = as.integer(saf_diagnosis != "No_MASLD"),
          batch2020 = as.integer(!grepl("^2019", liver_proteomics_filename)),
          age = as.numeric(alder), bmi = as.numeric(bmi), male = as.integer(gender == "Male"),
          Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)), Steatosis = as.numeric(steatosis_score),
          Ballooning = as.numeric(hepatocellular_ballooning_score), Inflammation = as.numeric(lobular_inflammation_score),
          NAS = as.numeric(nafld_activity_score))]
cat("participants", nrow(md), "MASLD", sum(md$masld), "control", sum(1 - md$masld), "\n")
cat("gender levels:", paste(unique(md$gender), collapse = ","), "; NA age/bmi:", sum(is.na(md$age)), sum(is.na(md$bmi)), "\n")
print(md[, .N, by = .(masld, batch2020)])
cat("raw protein groups", nrow(raw), "\n")
r1 <- raw[!is.na(Genes) & nzchar(Genes) & !grepl(";", Genes, fixed = TRUE)]
cat("single-gene groups", nrow(r1), "\n")
E <- as.matrix(r1[, md$sid, with = FALSE]); storage.mode(E) <- "double"
E[!is.finite(E) | E <= 0] <- NA; E <- log2(E); rownames(E) <- r1$Genes
E <- E[rowMeans(is.na(E)) <= 0.5, ]
cat("groups with <=50% missing", nrow(E), "\n")
E <- normalizeBetweenArrays(E, method = "quantile")
G <- do.call(rbind, lapply(split(seq_len(nrow(E)), rownames(E)), function(i)
  if (length(i) == 1L) E[i, ] else apply(E[i, , drop = FALSE], 2, median, na.rm = TRUE)))
G[!is.finite(G)] <- NA
cat("proteins (symbols) in DE family", nrow(G), "\n")

# ---- differential abundance: ~ group + batch + age + BMI + sex, BH over all proteins
X <- model.matrix(~ masld + batch2020 + age + bmi + male, data = md)
de <- list()
for (rob in c(TRUE, FALSE)) {
  f <- eBayes(lmFit(G, X), robust = rob)
  tt <- topTable(f, coef = "masld", number = Inf, sort.by = "none")
  de[[as.character(rob)]] <- data.table(gene = rownames(tt), logFC = tt$logFC, p = tt$P.Value,
                                        q = p.adjust(tt$P.Value, "BH"), robust = rob)
}
cat("protein BH<0.05: robust=TRUE", de[["TRUE"]][q < 0.05, .N], "; robust=FALSE", de[["FALSE"]][q < 0.05, .N], "\n")
pde <- fread(inp[["prod_de"]])
cmpde <- merge(de[["TRUE"]], pde[, .(gene, logFC_p = logFC, padj_p = padj)], by = "gene")
cat("producer DE rows", nrow(pde), "matched", nrow(cmpde), "max|dlogFC|", max(abs(cmpde$logFC - cmpde$logFC_p)),
    "max|dq|", max(abs(cmpde$q - cmpde$padj_p), na.rm = TRUE), "\n")
P <- de[["TRUE"]]

# ---- Figure 5B: join to 774-unit bulk by gene symbol (symbol unique in the bulk fit)
deg <- fread(inp[["deg"]], select = c("gene", "symbol", "logFC", "padj"))
cat("bulk rows", nrow(deg), "\n")
deg <- deg[!is.na(symbol) & nzchar(symbol)]
dupsym <- deg[, .N, by = symbol][N > 1, symbol]
cat("bulk symbols duplicated:", length(dupsym), "\n")
deg1 <- deg[!symbol %in% dupsym]
uni <- unique(trimws(readLines(inp[["universe"]]))); uni <- uni[nzchar(uni)]
cat("universe genes", length(uni), "\n")
J <- merge(P, deg1[, .(gene = symbol, rna_logFC = logFC, rna_padj = padj)], by = "gene")
J <- J[is.finite(logFC) & is.finite(rna_logFC)]
Ju <- J[gene %in% uni]
cat("universe proteins measured", length(intersect(P$gene, uni)), "; with bulk row", nrow(Ju), "\n")
rho <- cor(Ju$logFC, Ju$rna_logFC, method = "spearman")
sig <- Ju[q < 0.05]
res5b <- data.table(n_genes = nrow(Ju), spearman_rho = rho, sign_agree = mean(sign(Ju$logFC) == sign(Ju$rna_logFC)),
                    n_protein_sig = nrow(sig), n_same = sig[sign(logFC) == sign(rna_logFC), .N],
                    n_opp = sig[sign(logFC) != sign(rna_logFC), .N],
                    rnasig_same = sig[rna_padj < 0.05 & sign(logFC) == sign(rna_logFC), .N], rnasig_n = sig[rna_padj < 0.05, .N],
                    rnans_same = sig[rna_padj >= 0.05 & sign(logFC) == sign(rna_logFC), .N], rnans_n = sig[rna_padj >= 0.05, .N],
                    n_all_measured = nrow(J), rho_all = cor(J$logFC, J$rna_logFC, method = "spearman"))
print(res5b); fwrite(res5b, file.path(out, "fig5b_rederived.tsv"), sep = "\t")
p5b <- fread(inp[["prod_5B"]])
p5u <- p5b[mark == "scatter_point" | TRUE]
cat("producer 5B rows", nrow(p5b), "\n"); print(p5b[, .N, by = .(class)])
only_me <- setdiff(Ju$gene, p5b[class != "rest", gene]); only_prod <- setdiff(p5b[class != "rest", gene], Ju$gene)
cat("genes only in rederived universe set:", length(only_me), paste(head(only_me, 20), collapse = ","), "\n")
cat("genes only in producer universe set:", length(only_prod), paste(head(only_prod, 20), collapse = ","), "\n")
fwrite(data.table(side = c(rep("rederived_only", length(only_me)), rep("producer_only", length(only_prod))),
                  gene = c(only_me, only_prod)), file.path(out, "fig5b_gene_set_diff.tsv"), sep = "\t")
fwrite(Ju, file.path(out, "fig5b_rederived_genes.tsv"), sep = "\t")

# ---- Figure 5C: 25 fixed proteins x 5 histology scores within MASLD; BH over 125
fx <- fread(inp[["fixed25"]]); stopifnot(nrow(fx) == 25L, all(fx$gene %in% rownames(G)))
feats <- c("Steatosis", "Inflammation", "Ballooning", "Fibrosis", "NAS")
pc <- function(x, y, d) {
  Z <- cbind(1, d$age, d$bmi)
  if (length(unique(d$batch2020)) > 1) Z <- cbind(Z, d$batch2020)
  if (length(unique(d$male)) > 1) Z <- cbind(Z, d$male)
  rx <- qr.resid(qr(Z), rank(x)); ry <- qr.resid(qr(Z), rank(y))
  r <- cor(rx, ry); dfree <- length(x) - qr(Z)$rank - 1
  tt <- r * sqrt(dfree / (1 - r^2)); c(n = length(x), rho = r, p = 2 * pt(-abs(tt), dfree))
}
run5c <- function(keep) rbindlist(lapply(fx$gene, function(g) rbindlist(lapply(feats, function(f) {
  x <- G[g, md$sid]; y <- md[[f]]; ok <- keep & is.finite(x) & is.finite(y)
  s <- pc(x[ok], y[ok], md[ok]); data.table(gene = g, feature = f, n = s[["n"]], rho = s[["rho"]], p = s[["p"]])
}))))
m5 <- run5c(md$masld == 1L); m5[, q := p.adjust(p, "BH")]
a5 <- run5c(rep(TRUE, 58L)); a5[, q := p.adjust(p, "BH")]
cat("5C within MASLD: tests", nrow(m5), "q<0.05", m5[q < 0.05, .N], "; n range", min(m5$n), max(m5$n), "\n")
cat("5C all 58: q<0.05", a5[q < 0.05, .N], "\n")
p5c <- fread(inp[["prod_5C"]])
cc <- merge(m5, p5c[, .(gene, feature, rho_p = rho, q_p = padj_bh125, n_p = n)], by = c("gene", "feature"))
cat("5C producer match rows", nrow(cc), "max|drho|", max(abs(cc$rho - cc$rho_p)), "max|dq|", max(abs(cc$q - cc$q_p)),
    "n equal", all(cc$n == cc$n_p), "call flips", cc[(q < 0.05) != (q_p < 0.05), .N], "\n")
fwrite(m5, file.path(out, "fig5c_rederived_masld.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("U11A_DONE\n")
