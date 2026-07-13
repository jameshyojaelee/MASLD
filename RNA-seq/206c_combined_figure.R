#!/usr/bin/env Rscript
# 206c_combined_figure.R
# Generate one combined supplementary figure with all published DEG comparison
# panels: Hoang apple-to-apple + investigation + Govaere pairwise.
# Reads pre-computed CLM results from 206b; re-runs limma contrasts (fast).

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
INTB <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PUB  <- file.path(BASE, "data/published_degs")
OUTDIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 206c: Combined Published DEG Comparison Figure ===\n")

# ── Helpers ─────────────────────────────────────────────────────────────────
annot <- fread(file.path(INT, "gene_annotation/human_ensg_to_symbol.tsv"))
annot_pc <- annot[gene_type == "protein_coding"]
sym2ens <- setNames(annot_pc$gene_base, annot_pc$symbol)
ens2sym <- setNames(annot_pc$symbol, annot_pc$gene_base)
extra <- annot[!gene_base %in% annot_pc$gene_base & !symbol %in% names(sym2ens)]
extra <- extra[!duplicated(symbol)]
sym2ens <- c(sym2ens, setNames(extra$gene_base, extra$symbol))
sym2ens_upper <- setNames(sym2ens, toupper(names(sym2ens)))

map_symbols <- function(symbols) {
  symbols <- as.character(symbols)
  mapped <- sym2ens[symbols]
  miss <- is.na(mapped)
  if (any(miss)) mapped[miss] <- sym2ens_upper[toupper(symbols[miss])]
  mapped
}

lfc_cor <- function(dt1, dt2, l1, l2, id = "gene_base") {
  m <- merge(dt1[, c(id, l1), with = FALSE], dt2[, c(id, l2), with = FALSE], by = id)
  setnames(m, c(l1, l2), c("lfc1", "lfc2"))
  m <- m[is.finite(lfc1) & is.finite(lfc2)]
  n <- nrow(m)
  if (n < 10) return(list(n=n, r=NA, rho=NA, dir=NA, data=m))
  list(n=n, r=cor(m$lfc1,m$lfc2), rho=cor(m$lfc1,m$lfc2,method="spearman"),
       dir=mean(sign(m$lfc1)==sign(m$lfc2))*100, data=m)
}

scatter <- function(data, xlab, ylab, title) {
  anno <- sprintf("r = %.3f\nrho = %.3f\ndir = %.1f%%\nn = %s",
                  cor(data$lfc1,data$lfc2), cor(data$lfc1,data$lfc2,method="spearman"),
                  mean(sign(data$lfc1)==sign(data$lfc2))*100,
                  formatC(nrow(data), big.mark=","))
  ggplot(data, aes(x=lfc1, y=lfc2)) +
    geom_point(alpha=0.12, size=0.2, color=masld_colors$ns) +
    geom_hline(yintercept=0, linewidth=0.25, linetype="dashed", color="grey50") +
    geom_vline(xintercept=0, linewidth=0.25, linetype="dashed", color="grey50") +
    geom_smooth(method="lm", se=FALSE, linewidth=0.4, color=masld_colors$up) +
    geom_abline(slope=1, intercept=0, linewidth=0.25, linetype="dotted", color="grey40") +
    annotate("text", x=-Inf, y=Inf, label=anno, hjust=-0.05, vjust=1.2,
             size=1.6, color="grey30") +
    labs(x=xlab, y=ylab, title=title) +
    theme_masld(base_size=6)
}

recovery_rate <- function(ref, ql) sapply(ql, function(q) mean(ref %in% q)*100)

# ── Load data ───────────────────────────────────────────────────────────────
cat("Loading data...\n")
counts <- readRDS(file.path(INTB, "results/integration/merged_counts_raw.rds"))
sample_meta <- readRDS(file.path(INTB, "results/integration/meta_matched.rds"))
qc <- fread(file.path(INTB, "qc/sample_qc_report.csv"))
sample_meta <- sample_meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

run_limma <- function(sample_ids, design_mat, coef_name) {
  idx <- colnames(counts) %in% sample_ids
  dge <- DGEList(counts = counts[, idx])
  dge <- calcNormFactors(dge, method = "TMM")
  keep <- filterByExpr(dge, design = design_mat)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  v <- voom(dge, design_mat, plot = FALSE)
  fit <- lmFit(v, design_mat)
  if (is.character(coef_name) && !coef_name %in% colnames(fit$coefficients)) {
    fit <- contrasts.fit(fit, coef_name); fit <- eBayes(fit)
    tt <- topTable(fit, number = Inf, sort.by = "none")
  } else {
    fit <- eBayes(fit)
    tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  }
  tt$gene <- rownames(tt)
  dt <- as.data.table(tt); dt[, gene_base := sub("\\.\\d+$", "", gene)]; dt
}

# ── Hoang published ─────────────────────────────────────────────────────────
cat("Loading Hoang published DEGs...\n")
std_h <- function(dt) {
  setnames(dt, names(dt), tolower(gsub("\\s+","_",names(dt))))
  if ("gene_symbol" %in% names(dt)) setnames(dt, "gene_symbol", "symbol", skip_absent=TRUE)
  if ("adj_p" %in% names(dt)) setnames(dt, "adj_p", "adj_p_val", skip_absent=TRUE)
  dt[, gene_base := map_symbols(symbol)]; dt[!is.na(gene_base)]
}
hoang_nas <- std_h(as.data.table(read_excel(file.path(PUB,"GSE130970/MOESM2.xlsx"),
                                             sheet="NAS ordinal regression")))
hoang_fib <- std_h(as.data.table(read_excel(file.path(PUB,"GSE130970/MOESM2.xlsx"),
                                             sheet="fibrosis ordinal regression")))

# ── Hoang: our results ──────────────────────────────────────────────────────
cat("Computing Hoang apple-to-apple contrasts...\n")
m130 <- sample_meta[dataset=="GSE130970" & !is.na(nas_score)]
design_nas <- model.matrix(~ nas_score + age + sex, data=m130)
our_nas_ord <- run_limma(m130$sample_id, design_nas, "nas_score")

fib_ord_all <- fread(file.path(INT,"disease_signatures/fibrosis_ordinal_per_study.csv"))
fib_ord_all[, gene_base := sub("\\.\\d+$","",gene)]
our_fib_ord <- fib_ord_all[dataset=="GSE130970"]

dream <- fread(file.path(INT,"integration/dream_results.csv"))
dream[, gene_base := sub("\\.\\d+$","",gene)]
dream_fib_ord <- fread(file.path(INT,"../results/progression/c17_fibrosis_ordinal_dream.csv"))
dream_fib_ord[, gene_base := sub("\\.\\d+$","",gene)]

# ── Hoang: CLM results (from 206b) ─────────────────────────────────────────
cat("Loading pre-computed CLM results...\n")
clm_nas <- fread(file.path(OUTDIR, "clm_nas_results.csv"))
clm_fib <- fread(file.path(OUTDIR, "clm_fib_results.csv"))

# ── Govaere published ───────────────────────────────────────────────────────
cat("Loading Govaere supplementary tables...\n")
govaere_supp <- list()
govaere_xlsx <- file.path(PUB, "GSE135251/aba4448_supplementary_tables.xlsx")
supp_info <- data.table(
  sheet = paste("Table", c("S3","S4","S5","S6","S7","S8")),
  contrast = c("NASH_F2_vs_NAFL","NASH_F3_vs_NAFL","NASH_F4_vs_NAFL",
               "NASH_F3_vs_NASH_F01","NASH_F4_vs_NASH_F01","NAS_ge4"),
  expected_n = c(50, 907, 1369, 434, 1194, 369)
)
avail <- excel_sheets(govaere_xlsx)
for (i in seq_len(nrow(supp_info))) {
  sn <- supp_info$sheet[i]
  if (!sn %in% avail) next
  raw <- as.data.table(read_excel(govaere_xlsx, sheet=sn, skip=1))
  setnames(raw, names(raw), tolower(gsub("\\s+","_",names(raw))))
  nms <- names(raw)
  ens_col <- grep("ensembl|gene_id", nms, value=TRUE)[1]
  lfc_col <- grep("log2?fc|logfc", nms, value=TRUE)[1]
  if (!is.na(ens_col)) raw[, gene_base := sub("\\.\\d+$","",as.character(raw[[ens_col]]))]
  if (!is.na(lfc_col)) raw[, pub_logFC := as.numeric(raw[[lfc_col]])]
  govaere_supp[[supp_info$contrast[i]]] <- list(
    data=raw[!is.na(gene_base)], contrast=supp_info$contrast[i],
    has_lfc=!is.na(lfc_col))
}

# ── Govaere: our reproduced contrasts ───────────────────────────────────────
cat("Computing Govaere reproduced contrasts...\n")
m135 <- sample_meta[dataset=="GSE135251" & !is.na(fibrosis_stage)]
m135[, govaere_group := fcase(
  condition=="Control","Control", condition=="NAFL","NAFL",
  condition %in% "NASH" & fibrosis_stage %in% 0:1,"NASH_F01",
  condition %in% "NASH" & fibrosis_stage==2,"NASH_F2",
  condition %in% "NASH_Fibrosis" & fibrosis_stage==3,"NASH_F3",
  condition %in% "NASH_Fibrosis" & fibrosis_stage==4,"NASH_F4",
  default=NA_character_)]

govaere_contrasts <- list(
  NASH_F2_vs_NAFL=list(g1="NASH_F2",g2="NAFL"),
  NASH_F3_vs_NAFL=list(g1="NASH_F3",g2="NAFL"),
  NASH_F4_vs_NAFL=list(g1="NASH_F4",g2="NAFL"),
  NASH_F3_vs_NASH_F01=list(g1="NASH_F3",g2="NASH_F01"),
  NASH_F4_vs_NASH_F01=list(g1="NASH_F4",g2="NASH_F01"),
  NAS_ge4="NAS_binary")

our_govaere <- list()
for (cn in names(govaere_contrasts)) {
  info <- govaere_contrasts[[cn]]
  if (identical(info,"NAS_binary")) {
    ms <- m135[!is.na(nas_score)]
    ms[, nas_ge4 := factor(ifelse(nas_score>=4,"high","low"), levels=c("low","high"))]
    d <- model.matrix(~ nas_ge4 + inferred_sex, data=ms)
    our_govaere[[cn]] <- run_limma(ms$sample_id, d, "nas_ge4high")
  } else {
    ms <- m135[govaere_group %in% c(info$g1, info$g2)]
    ms[, grp := factor(govaere_group, levels=c(info$g2, info$g1))]
    d <- model.matrix(~ grp + inferred_sex, data=ms)
    our_govaere[[cn]] <- run_limma(ms$sample_id, d, paste0("grp",info$g1))
  }
}

# ── Govaere 25-gene signature ──────────────────────────────────────────────
ps_135251 <- fread(file.path(INT,"per_study/GSE135251_de_results.csv"))
ps_135251[, gene_base := sub("\\.\\d+$","",gene)]
govaere_25 <- data.table(
  symbol = c("AKR1B10","ANKRD29","CCL20","CFAP221","CLIC6","COL1A1","COL1A2",
             "DTNA","DUSP8","EPB41L4A","FERMT1","GDF15","HECW1","IL32","ITGBL1",
             "LTBP2","PDGFA","PPAPDC1A","RGS4","SCTR","STMN2","THY1",
             "TNFRSF12A","TYMS","HSD17B14"),
  pub_dir = c(rep("up",24),"down"))
govaere_25[, gene_base := map_symbols(symbol)]
g25 <- merge(govaere_25[!is.na(gene_base)],
             ps_135251[, .(gene_base, logFC, adj.P.Val)], by="gene_base", all.x=TRUE)
g25 <- merge(g25, dream[, .(gene_base, dream_lfc=logFC, dream_padj=padj)],  # C2-OK-sensitivity: dream_results.csv retired sensitivity arm, labeled "Dream" comparator
             by="gene_base", all.x=TRUE)

# ══════════════════════════════════════════════════════════════════════════
# BUILD ALL PANELS
# ══════════════════════════════════════════════════════════════════════════
cat("\nBuilding panels...\n")
PADJ <- 0.05; DPADJ <- 0.1
hoang_nas_degs <- hoang_nas[adj_p_val < 0.01, gene_base]
hoang_fib_degs <- hoang_fib[adj_p_val < 0.01, gene_base]

# --- Row 1: Hoang per-study CLM (same method as Hoang) + Dream ---
# CLM range_log2FC: same metric, same method, different quantification
clm_nas_rng <- merge(hoang_nas[!is.na(gene_base), .(gene_base, lfc1=range_log2fc)],
                     clm_nas[!is.na(gene_base), .(gene_base, lfc2=range_log2FC)], by="gene_base")
clm_nas_rng <- clm_nas_rng[is.finite(lfc1) & is.finite(lfc2)]
clm_fib_rng <- merge(hoang_fib[!is.na(gene_base), .(gene_base, lfc1=range_log2fc)],
                     clm_fib[!is.na(gene_base), .(gene_base, lfc2=range_log2FC)], by="gene_base")
clm_fib_rng <- clm_fib_rng[is.finite(lfc1) & is.finite(lfc2)]

p_h1 <- scatter(clm_nas_rng, "Hoang NAS ordinal","Our NAS ordinal (CLM)","NAS ordinal per-study")
p_h2 <- scatter(clm_fib_rng, "Hoang Fib ordinal","Our Fib ordinal (CLM)","Fib ordinal per-study")

h3 <- lfc_cor(hoang_nas, dream, "range_log2fc","logFC")
h4 <- lfc_cor(hoang_fib, dream_fib_ord, "range_log2fc","logFC")
p_h3 <- scatter(h3$data, "Hoang NAS ordinal","Dream disease vs ctrl","NAS ordinal vs Dream")
p_h4 <- scatter(h4$data, "Hoang Fib ordinal","Dream fib ordinal","Fib ordinal vs Dream")

# Recovery using CLM DEGs for per-study
clm_nas_degs <- clm_nas[!is.na(adj_P) & adj_P < 0.01, gene_base]
clm_fib_degs <- clm_fib[!is.na(adj_P) & adj_P < 0.01, gene_base]
rec_nas <- recovery_rate(hoang_nas_degs, list("Per-study"=clm_nas_degs,
                                               "Dream"=dream[padj<DPADJ,gene_base]))
rec_fib <- recovery_rate(hoang_fib_degs, list("Per-study"=clm_fib_degs,
                                               "Dream"=dream_fib_ord[padj<DPADJ,gene_base]))
rec_dt <- data.table(published=rep(c("NAS\n(2,970)","Fib\n(1,656)"),each=2),
                     method=rep(c("Per-study","Dream"),2),
                     recovery=c(rec_nas, rec_fib))
rec_dt[, method := factor(method, levels=c("Per-study","Dream"))]
p_rec <- ggplot(rec_dt, aes(x=published, y=recovery, fill=method)) +
  geom_col(position=position_dodge(0.7), width=0.6) +
  geom_text(aes(label=sprintf("%.0f%%",recovery)), position=position_dodge(0.7),
            vjust=-0.3, size=1.6) +
  scale_fill_manual(values=c("Per-study"=masld_colors$down, "Dream"=masld_colors$up), name=NULL) +
  labs(x=NULL, y="Recovery (%)", title="Hoang DEG recovery") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  theme_masld(base_size=6) + theme(legend.key.size=unit(0.3,"cm"))

# --- Row 2: Discordance analysis + variance decomposition ---
clm_fib_rng[, concordant := sign(lfc1)==sign(lfc2)]
idx <- colnames(counts) %in% m130$sample_id
dge_tmp <- DGEList(counts=counts[,idx])
dge_tmp <- calcNormFactors(dge_tmp, method="TMM")
keep_tmp <- filterByExpr(dge_tmp)
dge_tmp <- dge_tmp[keep_tmp,,keep.lib.sizes=FALSE]
l2cpm <- cpm(dge_tmp, log=TRUE, prior.count=1)
expr_dt <- data.table(gene_base=sub("\\.\\d+$","",rownames(l2cpm)), mean_expr=rowMeans(l2cpm))
disc_dt <- merge(clm_fib_rng, expr_dt, by="gene_base")
disc_dt[, expr_q := cut(mean_expr, quantile(mean_expr,0:4/4), include.lowest=TRUE,
                        labels=c("Q1\n(low)","Q2","Q3","Q4\n(high)"))]
disc_by_q <- disc_dt[, .(n=.N, n_disc=sum(!concordant), pct=100*mean(!concordant)), by=expr_q]

p_disc <- ggplot(disc_by_q, aes(x=expr_q, y=pct)) +
  geom_col(fill=masld_colors$up, width=0.6) +
  geom_text(aes(label=sprintf("%.1f%%",pct)), vjust=-0.3, size=1.6) +
  labs(x="Expression quartile", y="Direction discordance (%)", title="Discordance by expression") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  theme_masld(base_size=6)

# Variance decomposition
r_clm_fib <- cor(clm_fib_rng$lfc1, clm_fib_rng$lfc2)
r_clm_nas <- cor(clm_nas_rng$lfc1, clm_nas_rng$lfc2)
decomp <- data.table(
  comparison=c("Govaere\nlimma vs limma\n(HT-Seq vs featureCounts)",
               "Hoang CLM\nsame method\n(Salmon vs featureCounts)"),
  r=c(0.955, mean(c(r_clm_fib, r_clm_nas))),
  quant=c("HT-Seq vs\nfeatureCounts","Salmon vs\nfeatureCounts"))
decomp[, comparison := factor(comparison, levels=rev(comparison))]
p_decomp <- ggplot(decomp, aes(x=r, y=comparison)) +
  geom_col(width=0.5, fill=masld_colors$down) +
  geom_text(aes(label=sprintf("%.3f",r)), hjust=-0.1, size=1.8) +
  labs(x="Pearson r (same method, different quantification)", y=NULL,
       title="Quantification gap") +
  xlim(0,1.1) +
  theme_masld(base_size=6)

# --- Row 3-4: Govaere pairwise (j-o) ---
contrast_labels <- c(NASH_F2_vs_NAFL="NASH F2 vs NAFL", NASH_F3_vs_NAFL="NASH F3 vs NAFL",
                     NASH_F4_vs_NAFL="NASH F4 vs NAFL", NASH_F3_vs_NASH_F01="NASH F3 vs F0/1",
                     NASH_F4_vs_NASH_F01="NASH F4 vs F0/1", NAS_ge4="NAS >= 4")
gov_panels <- list()
gov_summary <- list()
for (cn in names(contrast_labels)) {
  if (!cn %in% names(govaere_supp) || !cn %in% names(our_govaere)) next
  pub <- govaere_supp[[cn]]; our <- our_govaere[[cn]]
  res <- lfc_cor(pub$data, our, "pub_logFC","logFC")
  gov_panels[[cn]] <- scatter(res$data,
    paste("Govaere",contrast_labels[cn],"(logFC)"), "Our reproduced (logFC)",
    contrast_labels[cn])
  gov_summary[[cn]] <- data.table(contrast=contrast_labels[cn], r=res$r, pub_n=nrow(pub$data))
}

# --- Row 5: 25-gene + Govaere correlation summary ---
g25_plot <- melt(g25[, .(symbol, pub_dir, logFC, adj.P.Val, dream_lfc, dream_padj)],  # C2-OK-sensitivity: dream arm of published comparison figure
                 id.vars=c("symbol","pub_dir"),
                 measure.vars=list(lfc=c("logFC","dream_lfc"), padj=c("adj.P.Val","dream_padj")))  # C2-OK-sensitivity
g25_plot[, source := fifelse(variable==1, "Per-study","Dream")]
g25_plot[, sig := padj < 0.05]
g25_plot[, symbol := factor(symbol, levels=g25[order(logFC)]$symbol)]

p_25gene <- ggplot(g25_plot, aes(x=source, y=symbol)) +
  geom_point(aes(color=ifelse(lfc>0,"up","down"), size=pmin(-log10(padj),10), shape=sig)) +
  scale_color_manual(values=c(up=masld_colors$up, down=masld_colors$down), guide="none") +
  scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1), labels=c("NS","Sig"), name=NULL) +
  scale_size_continuous(range=c(0.4,2.5), name=expression(-log[10](padj))) +
  labs(x=NULL, y=NULL, title="Govaere 25-gene signature") +
  theme_masld(base_size=6) + theme(axis.text.y=element_text(size=4),
                                    legend.key.size=unit(0.3,"cm"))

sum_dt <- rbindlist(gov_summary)
sum_dt[, contrast := factor(contrast, levels=rev(contrast))]
p_corsum <- ggplot(sum_dt, aes(x=r, y=contrast)) +
  geom_point(aes(size=pub_n), color=masld_colors$up) +
  geom_text(aes(label=sprintf("%.3f",r)), hjust=-0.3, size=1.8) +
  scale_size_continuous(range=c(1,3), name="Published\nDEGs") +
  labs(x="Pearson r (published vs reproduced)", y=NULL, title="Govaere correlation summary") +
  xlim(0.85,1.05) +
  theme_masld(base_size=6) + theme(legend.key.size=unit(0.3,"cm"))

# ══════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ══════════════════════════════════════════════════════════════════════════
cat("Assembling combined figure...\n")

# Row 1: Hoang per-study CLM + Dream (5 panels)
row1 <- p_h1 | p_h2 | p_h3 | p_h4 | p_rec

# Row 2: Discordance + variance decomp (3 panels)
row2 <- p_disc | p_decomp | plot_spacer()

# Row 3: Govaere pairwise top (3 panels)
row3 <- gov_panels[["NASH_F2_vs_NAFL"]] | gov_panels[["NASH_F3_vs_NAFL"]] | gov_panels[["NASH_F4_vs_NAFL"]]

# Row 4: Govaere pairwise bottom (3 panels)
row4 <- gov_panels[["NASH_F3_vs_NASH_F01"]] | gov_panels[["NASH_F4_vs_NASH_F01"]] | gov_panels[["NAS_ge4"]]

# Row 5: 25-gene + correlation summary
row5 <- p_25gene | p_corsum | plot_spacer()

fig <- (row1 / row2 / row3 / row4 / row5) +
  plot_layout(heights = c(1, 0.8, 1, 1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 7, face = "bold"))

save_fig(fig, file.path(OUTDIR, "published_deg_combined.pdf"),
         width = fig_full_width, height = 12)
cat(sprintf("  Saved %s\n", file.path(OUTDIR, "published_deg_combined.pdf")))

# Also save metrics table
metrics <- list()
for (cn in names(gov_summary)) {
  pub <- govaere_supp[[cn]]; our <- our_govaere[[cn]]
  res <- lfc_cor(pub$data, our, "pub_logFC","logFC")
  metrics[[cn]] <- data.table(study="Govaere", comparison=cn, r=res$r, rho=res$rho,
                               dir_conc=res$dir, n=res$n)
}
metrics[["H_NAS_CLM"]] <- data.table(study="Hoang", comparison="NAS_ordinal_CLM_per_study",
                                      r=r_clm_nas, rho=cor(clm_nas_rng$lfc1,clm_nas_rng$lfc2,method="spearman"),
                                      dir_conc=mean(sign(clm_nas_rng$lfc1)==sign(clm_nas_rng$lfc2))*100, n=nrow(clm_nas_rng))
metrics[["H_Fib_CLM"]] <- data.table(study="Hoang", comparison="Fib_ordinal_CLM_per_study",
                                      r=r_clm_fib, rho=cor(clm_fib_rng$lfc1,clm_fib_rng$lfc2,method="spearman"),
                                      dir_conc=mean(clm_fib_rng$concordant)*100, n=nrow(clm_fib_rng))
metrics[["H_NAS_dream"]] <- data.table(study="Hoang", comparison="NAS_ordinal_vs_dream",
                                        r=h3$r, rho=h3$rho, dir_conc=h3$dir, n=h3$n)
metrics[["H_Fib_dream"]] <- data.table(study="Hoang", comparison="Fib_ordinal_vs_dream",
                                        r=h4$r, rho=h4$rho, dir_conc=h4$dir, n=h4$n)
fwrite(rbindlist(metrics), file.path(OUTDIR, "published_vs_perstudy_metrics.csv"))
fwrite(g25, file.path(OUTDIR, "govaere_25gene_detail.csv"))
cat("  Saved metrics and 25-gene detail\n")

cat("\n=== Done ===\n")
