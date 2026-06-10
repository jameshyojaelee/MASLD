#!/usr/bin/env Rscript
# figS_degx_mega_validation.R
# Full-data between-method concordance for figS_methods_validation/mega_validation/,
# rebuilt from the degx R-exact fits (real R/Bioconductor methods on the canonical
# 5-cohort disease_vs_control mega data). Replaces the legacy 6-engine panels.
#
# mega_validation's angle = EFFECT-SIZE concordance (Spearman rho on logFC) +
# consensus core, complementary to multimethod_validation's SET-overlap panels
# (Jaccard / CAT). Same FE/HK exclusion as the comparison panels (FE invalid under
# heterogeneity; HK degenerate) -> 15 methods.
#
# Publication theme; control/neutral = #9E9E9E; cairo_pdf; ASCII. CPU.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

REX_DIR <- file.path(Sys.getenv("HOME"), "degx", "runs", "Rexact")
CONTRAST <- "disease_vs_control"
META_CSV <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/meta_analysis_results.csv")
OUT <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DEGX_VERSION <- "0.0.1"
MF_DROP <- c("metafor_deseq2_hk","metafor_voom_hk","metafor_deseq2_fe","metafor_voom_fe")

rex_labels <- c(deseq2_wald="DESeq2 (Wald)", deseq2_lrt="DESeq2 (LRT)", deseq2_apeglm="DESeq2 (apeglm)",
  deseq2_ashr="DESeq2 (ashr)", edger_qlf="edgeR (QLF)", edger_qlf_robust="edgeR (QLF-robust)",
  edger_exact="edgeR (exact)", limma_voom="limma-voom", limma_trend="limma-trend",
  limma_voom_qw="limma-voom (QW)", dream="dream", metafor_deseq2_re="metafor DESeq2 (RE)",
  metafor_voom_re="metafor voom (RE)", combatseq_deseq2="ComBat-seq+DESeq2", sva_limma="sva+limma")
pmr <- function(x){x<-as.character(x); ifelse(x %in% names(rex_labels), rex_labels[x], x)}
green_sc <- colorRampPalette(c("#dae7c7","#7fb069","#193c1e"))(12)
prov <- function(...) NULL   # on-panel Source captions removed (provenance lives in the README)
written <- character(0); mark <- function(p){written<<-c(written,p); message("  wrote ",basename(p))}

# ---- load the 15 R-exact method tables (gene, logFC, pval, padj) ----------
files <- list.files(REX_DIR, pattern=paste0("^Rexact_",CONTRAST,"_.*\\.csv$"), full.names=TRUE)
files <- files[!grepl("manifest", files)]
meth <- sub(paste0("^Rexact_",CONTRAST,"_"), "", sub("\\.csv$","", basename(files)))
keep <- !meth %in% MF_DROP
files <- files[keep]; meth <- meth[keep]
tabs <- setNames(lapply(files, function(f) read.csv(f, stringsAsFactors=FALSE)), meth)
message(sprintf("[mega] %d methods: %s", length(meth), paste(meth, collapse=", ")))

# wide logFC matrix (gene x method)
lfc <- Reduce(function(a,b) full_join(a,b,by="gene"),
              lapply(meth, function(m) tabs[[m]][,c("gene","logFC")] %>% rename(!!m := logFC)))
LM <- as.matrix(lfc[,-1]); rownames(LM) <- lfc$gene

# ===========================================================================
# A. Effect-size concordance: 15x15 Spearman rho on logFC, clustered
# ===========================================================================
rho <- cor(LM, method="spearman", use="pairwise.complete.obs")
hc <- hclust(as.dist(1-rho), method="average"); ord <- rownames(rho)[hc$order]
long <- as.data.frame(as.table(rho)); names(long) <- c("mi","mj","rho")
long <- long %>% mutate(mi=factor(pmr(mi),levels=pmr(ord)), mj=factor(pmr(mj),levels=pmr(ord)))
pA <- ggplot(long, aes(mi,mj,fill=rho)) +
  geom_tile(color="white", linewidth=0.3) +
  scale_fill_gradientn(colours=green_sc, limits=c(min(rho),1), oob=scales::squish, name="Spearman\nrho") +
  coord_fixed() +
  labs(title="Effect-size concordance across methods (Spearman rho on log2FC)",
       subtitle=sprintf("Full-data mega-analysis; all pairs >= %.2f -> methods agree on effect-size rank", min(rho[rho<1])),
       x=NULL, y=NULL, caption=prov("15 methods, clustered")) +
  theme_masld()+theme_pub()+
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="right", legend.key.width=unit(0.18,"cm"))
save_fig(pA, file.path(OUT,"panelI3_logfc_concordance.pdf"), width=6.6, height=6.2)
mark(file.path(OUT,"panelI3_logfc_concordance.pdf"))

# ===========================================================================
# B. Consensus: # genes called DEG by k of 15 methods, stacked by up/down direction
#    (direction = sign of dream's log2FC). Two thresholds: padj<0.05; +|log2FC|>0.2.
# ===========================================================================
nM <- length(meth)
dd <- tabs[["dream"]]; dream_sign <- setNames(sign(dd$logFC), dd$gene)
consensus_tbl <- function(dsets){
  allg <- unique(unlist(dsets))
  kc  <- rowSums(sapply(dsets, function(s) allg %in% s))
  dir <- dream_sign[allg]
  dir <- ifelse(is.na(dir) | dir == 0, NA_character_, ifelse(dir > 0, "up", "down"))
  ok  <- kc >= 1 & !is.na(dir)
  tb  <- as.data.frame(table(k = factor(kc[ok], levels = 1:nM),
                             dir = factor(dir[ok], levels = c("up","down"))))
  names(tb) <- c("k","dir","n"); tb$k <- as.integer(as.character(tb$k)); tb
}
plot_consensus <- function(tb, ttl){
  ggplot(tb, aes(k, n, fill = dir)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = c(up = masld_colors$up, down = masld_colors$down),
                      name = NULL, labels = c(up = "up", down = "down")) +
    scale_x_continuous(breaks = 1:nM) +
    labs(title = ttl, x = "# methods calling the gene a DEG", y = "# genes") +
    theme_masld() + theme_pub() + theme(legend.position = "top")
}
ds05     <- lapply(tabs, function(d) d$gene[!is.na(d$padj) & d$padj < 0.05])
ds05_lfc <- lapply(tabs, function(d) d$gene[!is.na(d$padj) & d$padj < 0.05 & !is.na(d$logFC) & abs(d$logFC) > 0.2])
save_fig(plot_consensus(consensus_tbl(ds05), "Consensus DEGs (padj < 0.05)"),
         file.path(OUT,"panelI4_consensus_core.pdf"), width=fig_col_width, height=3.6)
mark(file.path(OUT,"panelI4_consensus_core.pdf"))
save_fig(plot_consensus(consensus_tbl(ds05_lfc), "Consensus DEGs (padj < 0.05, |log2FC| > 0.2)"),
         file.path(OUT,"panelI4c_consensus_lfc02.pdf"), width=fig_col_width, height=3.6)
mark(file.path(OUT,"panelI4c_consensus_lfc02.pdf"))

# I4b -- UpSet of the 15 DEG sets (padj<0.05), top intersections by size
if (requireNamespace("UpSetR", quietly=TRUE)) {
  suppressPackageStartupMessages(library(UpSetR))
  ds_pretty <- setNames(ds05, pmr(names(ds05)))
  tmpf <- tempfile(fileext=".pdf"); fin <- file.path(OUT,"panelI4b_upset.pdf")
  grDevices::pdf(tmpf, width=9.5, height=5.2)
  print(UpSetR::upset(fromList(ds_pretty), nsets=length(ds_pretty), nintersects=25,
        order.by="freq", mb.ratio=c(0.62,0.38), text.scale=0.85, point.size=1.6, line.size=0.4,
        mainbar.y.label="genes in intersection", sets.x.label="DEGs per method",
        main.bar.color="#518dc9", sets.bar.color="#9E9E9E", matrix.color="#1565C0"))
  dev.off()
  # UpSetR's print() emits a leading blank page; the plot is page 2 -> extract it
  gsbin <- Sys.which("gs"); ok <- FALSE
  if (nzchar(gsbin)) {
    rc <- tryCatch(system2(gsbin, c("-q","-dBATCH","-dNOPAUSE","-dFirstPage=2","-dLastPage=2",
            "-sDEVICE=pdfwrite", paste0("-sOutputFile=",fin), tmpf)), error=function(e) 1L)
    ok <- (identical(rc, 0L) && file.exists(fin))
  }
  if (!ok) file.copy(tmpf, fin, overwrite=TRUE)
  written <- c(written, fin); message("  wrote panelI4b_upset.pdf")
}

# ===========================================================================
# C. dream vs metafor-RE effect sizes (the conservative family), I2-coloured
# ===========================================================================
if (all(c("dream","metafor_voom_re") %in% meth)) {
  mg <- merge(tabs[["dream"]][,c("gene","logFC")], tabs[["metafor_voom_re"]][,c("gene","logFC")],
              by="gene", suffixes=c("_dream","_meta"))
  i2 <- NULL
  if (file.exists(META_CSV)) {
    md <- tryCatch(read.csv(META_CSV, stringsAsFactors=FALSE), error=function(e) NULL)
    gc <- intersect(c("gene","gene_name","gene_id"), names(md))[1]
    i2c <- intersect(c("meta_I2","I2","i2"), names(md))[1]
    if (!is.na(gc) && !is.na(i2c)) { md$gene <- md[[gc]]; mg <- left_join(mg, md[,c("gene",i2c)] %>% rename(I2=all_of(i2c)), by="gene"); i2 <- TRUE }
  }
  rr <- suppressWarnings(cor(mg$logFC_dream, mg$logFC_meta, method="spearman", use="complete.obs"))
  lim <- as.numeric(quantile(abs(c(mg$logFC_dream, mg$logFC_meta)), 0.999, na.rm=TRUE))
  base <- ggplot(mg, aes(logFC_dream, logFC_meta))
  pt <- if (!is.null(i2)) geom_point(aes(color=I2), size=0.3, alpha=0.4) else geom_point(size=0.3, alpha=0.3, color="#1565C0")
  pC <- base + rasterize_layer(pt) +
    geom_abline(slope=1,intercept=0, linetype="dashed", linewidth=0.3, color="grey45") +
    coord_fixed(xlim=c(-lim,lim), ylim=c(-lim,lim)) +
    annotate("text", x=-lim*0.95, y=lim*0.9, hjust=0, size=PUB_GEOM_TEXT+0.3,
             label=sprintf("Spearman rho = %.3f\nn = %s genes", rr, format(nrow(mg),big.mark=","))) +
    labs(title="dream vs metafor random-effects: effect-size agreement",
         subtitle="The conservative metafor-RE calls fewer DEGs but agrees with dream on log2FC",
         x="dream log2FC", y="metafor voom (RE) log2FC", caption=prov("R-exact")) +
    theme_masld()+theme_pub()
  if (!is.null(i2)) pC <- pC + scale_color_gradientn(colours=colorRampPalette(c("#dae7c7","#193c1e"))(12), name="I2 (%)", na.value="grey85")
  save_fig(pC, file.path(OUT,"panelI5_dream_vs_metafor.pdf"), width=fig_half_width, height=fig_half_width)
  mark(file.path(OUT,"panelI5_dream_vs_metafor.pdf"))
}

message("\n[degx mega_validation] wrote ", length(written), " panel(s) to:\n  ", OUT)
invisible(lapply(written, function(p) cat("   -", basename(p), "\n")))
