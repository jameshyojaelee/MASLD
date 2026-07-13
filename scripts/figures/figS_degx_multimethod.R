#!/usr/bin/env Rscript
# figS_degx_multimethod.R
# Head-to-head DE-method comparison panels for figS_methods_validation/multimethod_validation/,
# regenerated from the degx GPU robustness battery (replaces the legacy 3-method
# dream/DESeq2/metafor panels, archived under archive_3method_pre_degx_2026-06-04/).
#
# DESIGN LANGUAGE for the 11-19 method matrix (the old 3-method designs do not scale):
#   * N x N clustered heatmaps replace pairwise hexbins (concordance, agreement).
#   * Sorted horizontal bars / lollipops, coloured by METHOD FAMILY, replace
#     dodged 3-colour bars. No per-bar value labels (colour + axis carry it).
#   * Calibration panels use a diverging colourbar centred on the nominal alpha
#     (blue = conservative -> neutral gray = calibrated -> magenta = inflated).
#   * Family (voom / NB-GLM / dream / metafor; R-exact also DESeq2 / edgeR /
#     limma / batch-corr) is the categorical axis -- never 15+ raw colours.
#
# Colours from ~/publication_color_themes.R via publication_theme.R. Control/
# neutral = #9E9E9E. PDF only (cairo_pdf, useDingbats=FALSE). ASCII text. CPU job.
# Data-driven: panels whose degx inputs are missing are skipped with a message.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
# Shared engine/correction/k_sv label lookups. The local method_labels / rex_labels
# below are kept (they already carry the NB-sim + R-exact variant names) but the
# canonical engine names match scripts/figures/method_correction_labels.R.
source(file.path(BASE, "scripts/figures/method_correction_labels.R"))

DEGX_RUNS <- Sys.getenv("DEGX_RUNS", file.path(Sys.getenv("HOME"), "degx", "runs", "robustness"))
REX_DIR   <- file.path(dirname(DEGX_RUNS), "Rexact")
DEGX_VERSION <- "0.0.1"
CONTRAST  <- "disease_vs_control"; CONTRAST_LAB <- "disease vs control"
NOMINAL   <- 0.05; CAL_LO <- 0.04; CAL_HI <- 0.06

OUT <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- pretty method labels (sim / resample set) -----------------------------
method_labels <- c(
  voom="voom", voom_tmm="voom (TMM)", voom_robust="voom (robust)", voom_ashr="voom (ashr)",
  nb_glm="NB-GLM (Wald)", nb_glm_lrt="NB-GLM (LRT)", nb_glm_tmm="NB-GLM (TMM)", nb_glm_ashr="NB-GLM (ashr)",
  dream="dream", dream_tmm="dream (TMM)", dream_ashr="dream (ashr)",
  metafor="metafor (RE)", metafor_deseq2="metafor (DESeq2)", metafor_fe="metafor (FE)", metafor_hk="metafor (HK)")
pm <- function(x){x<-as.character(x); ifelse(x %in% names(method_labels), method_labels[x], x)}

# --- family maps + palettes (guideline colours) ----------------------------
fam_resample <- function(m){m<-as.character(m)
  ifelse(grepl("^dream",m),"dream", ifelse(grepl("^nb_glm",m),"NB-GLM",
  ifelse(grepl("^voom",m),"voom", ifelse(grepl("^metafor",m),"metafor","other"))))}
fam_cols <- c(voom="#4aa2c2", "NB-GLM"="#C9265E", dream="#40b499", metafor="#9b75d6")

fam_rex <- function(m){m<-as.character(m)
  ifelse(grepl("^deseq2",m),"DESeq2", ifelse(grepl("^edger",m),"edgeR",
  ifelse(grepl("^limma",m),"limma", ifelse(grepl("^dream",m),"dream",
  ifelse(grepl("^metafor",m),"metafor", ifelse(grepl("combatseq|^sva",m),"batch-corr","other"))))))}
fam_cols_rex <- c(DESeq2="#4aa2c2", edgeR="#518dc9", limma="#9b75d6",
                  dream="#40b499", metafor="#e37faf", "batch-corr"="#f7bf87")
rex_labels <- c(deseq2_wald="DESeq2 (Wald)", deseq2_lrt="DESeq2 (LRT)", deseq2_apeglm="DESeq2 (apeglm)",
  deseq2_ashr="DESeq2 (ashr)", edger_qlf="edgeR (QLF)", edger_qlf_robust="edgeR (QLF-robust)",
  edger_exact="edgeR (exact)", limma_voom="limma-voom", limma_trend="limma-trend",
  limma_voom_qw="limma-voom (QW)", dream="dream", metafor_deseq2_re="metafor DESeq2 (RE)",
  metafor_deseq2_fe="metafor DESeq2 (FE)", metafor_deseq2_hk="metafor DESeq2 (HK)",
  metafor_voom_re="metafor voom (RE)", metafor_voom_fe="metafor voom (FE)",
  metafor_voom_hk="metafor voom (HK)", combatseq_deseq2="ComBat-seq+DESeq2", sva_limma="sva+limma")
pmr <- function(x){x<-as.character(x); ifelse(x %in% names(rex_labels), rex_labels[x], x)}

# --- shared scales ---------------------------------------------------------
green_sc <- colorRampPalette(c("#dae7c7","#7fb069","#193c1e"))(12)
cal_fill <- function(limits=NULL,name="type-I\nerror")
  scale_fill_gradient2(low="#1565C0", mid="#9E9E9E", high="#C9265E",
                       midpoint=NOMINAL, limits=limits, oob=scales::squish, name=name)
cal_color <- function(limits=NULL,name="type-I\nerror")
  scale_color_gradient2(low="#1565C0", mid="#9E9E9E", high="#C9265E",
                        midpoint=NOMINAL, limits=limits, oob=scales::squish, name=name)
sym_lim <- function(v){d<-max(abs(v-NOMINAL),na.rm=TRUE); NOMINAL+c(-1,1)*(if(is.finite(d)&&d>0)d else 0.02)}
prov <- function(...) NULL   # on-panel Source captions removed (provenance lives in the README)
written <- character(0); mark <- function(p){written<<-c(written,p); message("  wrote ",basename(p))}

rd <- function(...) { f <- file.path(...); if (file.exists(f)) utils::read.csv(f, check.names=FALSE, stringsAsFactors=FALSE) else NULL }

# ===========================================================================
# G. Type-I calibration -- parametric (NB sim) vs permutation (real-data) null
#    Merges the two look-alike type-I bars into one cross-check scatter.
# G2. Parametric type-I vs n per method (moved here from the retired degx_robustness/).
# ===========================================================================
perm <- rd(DEGX_RUNS, CONTRAST, paste0("resample_permutation_summary_",CONTRAST,".csv"))
simn <- rd(DEGX_RUNS, "sim", "sim_null", "batched_summary.csv")
if (!is.null(simn)) simn <- simn[!simn$estimator %in% c("metafor_fe","metafor_hk"), , drop=FALSE]
if (!is.null(perm) && !is.null(simn)) {
  message("[render] G calibration: parametric vs permutation null")
  para <- simn %>% group_by(estimator) %>%
    summarise(parametric=mean(type_i_mean, na.rm=TRUE), .groups="drop")
  cal <- dplyr::inner_join(para, perm %>% transmute(estimator, permutation=type_i_mean), by="estimator") %>%
    mutate(fam=fam_resample(estimator), lab=pm(estimator))
  rng <- range(c(cal$parametric, cal$permutation, CAL_LO, 0.072), na.rm=TRUE)
  lab_layer <- if (requireNamespace("ggrepel", quietly=TRUE))
      ggrepel::geom_text_repel(aes(label=lab), size=PUB_GEOM_TEXT, color="grey25",
                               max.overlaps=20, seed=42, segment.size=0.2, min.segment.length=0)
    else geom_text(aes(label=lab), size=PUB_GEOM_TEXT, color="grey25", vjust=-0.7)
  pG <- ggplot(cal, aes(parametric, permutation, color=fam)) +
    annotate("rect", xmin=CAL_LO,xmax=CAL_HI, ymin=-Inf,ymax=Inf, fill="#9E9E9E", alpha=0.10) +
    annotate("rect", ymin=CAL_LO,ymax=CAL_HI, xmin=-Inf,xmax=Inf, fill="#9E9E9E", alpha=0.10) +
    geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=0.3, color="grey55") +
    geom_vline(xintercept=NOMINAL, linetype="dotted", linewidth=0.25, color="grey45") +
    geom_hline(yintercept=NOMINAL, linetype="dotted", linewidth=0.25, color="grey45") +
    geom_point(size=2) + lab_layer +
    scale_color_manual(values=fam_cols, name="family") +
    coord_fixed(xlim=rng, ylim=rng) +
    labs(x="parametric (simulation) type-I error", y="permutation (real-data) type-I error",
         caption=NULL) +
    theme_masld()+theme_pub()+theme(legend.position="right")
  save_fig(pG, file.path(OUT,"calibration_two_nulls.pdf"), width=fig_col_width, height=fig_col_width*0.92)
  mark(file.path(OUT,"calibration_two_nulls.pdf"))

  # G2 -- parametric type-I vs samples per group, small multiples
  ord2 <- para %>% arrange(parametric) %>% pull(estimator)
  b2g <- simn %>% transmute(lab=factor(pm(estimator), levels=pm(ord2)),
                            n_per_group=factor(n_per_group, levels=sort(unique(n_per_group))),
                            type_i_mean)
  pG2 <- ggplot(b2g, aes(n_per_group, type_i_mean, group=lab)) +
    annotate("rect", xmin=-Inf,xmax=Inf,ymin=CAL_LO,ymax=CAL_HI, fill="#9E9E9E", alpha=0.12) +
    geom_hline(yintercept=NOMINAL, linetype="dashed", linewidth=0.25, color="grey40") +
    geom_line(linewidth=0.4, color="grey60") + geom_point(aes(color=type_i_mean), size=1.4) +
    cal_color(limits=sym_lim(b2g$type_i_mean)) +
    facet_wrap(~lab, ncol=5) +
    labs(x="samples per group", y="type-I error", caption=prov("sim_null 2,000 reps")) +
    theme_masld()+theme_pub()+
    theme(legend.position="right", legend.key.width=unit(0.18,"cm"), panel.spacing=unit(0.45,"lines"))
  save_fig(pG2, file.path(OUT,"typeI_vs_n.pdf"), width=fig_full_width, height=4.0)
  mark(file.path(OUT,"typeI_vs_n.pdf"))
}

# ===========================================================================
# G1. Type-I calibration -- PERMUTATION (real-data) null ONLY.
#     The two-nulls scatter with the parametric axis dropped: per-method
#     permutation type-I vs the nominal 0.05. Renders on `perm` alone (no sim).
# ===========================================================================
if (!is.null(perm)) {
  message("[render] G1 permutation-only calibration")
  pc <- perm %>% transmute(estimator, ti=type_i_mean, mcse=type_i_mcse,
                           fam=fam_resample(estimator), lab=pm(estimator)) %>%
    arrange(ti) %>% mutate(lab=factor(lab, levels=lab))
  n_perm <- format(max(perm$n_resample, na.rm=TRUE), big.mark=",")
  pG1 <- ggplot(pc, aes(lab, ti, color=ti)) +
    annotate("rect", xmin=-Inf,xmax=Inf, ymin=CAL_LO,ymax=CAL_HI, fill="#9E9E9E", alpha=0.12) +
    geom_hline(yintercept=NOMINAL, linetype="dashed", linewidth=0.3, color="grey40") +
    geom_errorbar(aes(ymin=ti-1.96*mcse, ymax=ti+1.96*mcse), width=0.25, linewidth=0.35, color="grey60") +
    geom_point(size=2.2) +
    cal_color(limits=sym_lim(pc$ti)) +
    scale_y_continuous(expand=expansion(mult=c(0.05,0.07))) +
    labs(x=NULL, y="permutation type-I error", caption=NULL) +
    theme_masld()+theme_pub()+
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))
  save_fig(pG1, file.path(OUT,"calibration_permutation.pdf"), width=fig_col_width, height=3.8)
  mark(file.path(OUT,"calibration_permutation.pdf"))
}

# ===========================================================================
# G3. Permutation-null type-I vs n -- real-data analog of G2 (perm_vs_n GPU job).
#     n = balanced samples PER GROUP, subsampled from real data (Control-limited,
#     max 156); labels permuted within cohort at each n.
# ===========================================================================
pvn <- rd(DEGX_RUNS, CONTRAST, "perm_vs_n", "permutation_vs_n_summary.csv")
if (!is.null(pvn)) {
  message("[render] G3 permutation type-I vs n")
  ordp <- pvn %>% group_by(estimator) %>% summarise(m=mean(type_i_mean), .groups="drop") %>%
    arrange(m) %>% pull(estimator)
  d3 <- pvn %>% mutate(lab=factor(pm(estimator), levels=pm(ordp)),
                       nf=factor(n_per_group, levels=sort(unique(n_per_group))))
  pG3 <- ggplot(d3, aes(nf, type_i_mean, group=lab)) +
    annotate("rect", xmin=-Inf,xmax=Inf,ymin=CAL_LO,ymax=CAL_HI, fill="#9E9E9E", alpha=0.12) +
    geom_hline(yintercept=NOMINAL, linetype="dashed", linewidth=0.25, color="grey40") +
    geom_errorbar(aes(ymin=type_i_mean-1.96*type_i_mcse, ymax=type_i_mean+1.96*type_i_mcse),
                  width=0.2, linewidth=0.25, color="grey65") +
    geom_line(linewidth=0.4, color="grey60") + geom_point(aes(color=type_i_mean), size=1.4) +
    cal_color(limits=sym_lim(d3$type_i_mean)) +
    facet_wrap(~lab, ncol=4) +
    labs(x="samples per group (balanced)", y="permutation type-I error",
         caption=prov("1,000 permutations per n")) +
    theme_masld()+theme_pub()+
    theme(legend.position="right", legend.key.width=unit(0.18,"cm"), panel.spacing=unit(0.45,"lines"))
  save_fig(pG3, file.path(OUT,"perm_typeI_vs_n.pdf"), width=fig_full_width, height=4.0)
  mark(file.path(OUT,"perm_typeI_vs_n.pdf"))
}

# ===========================================================================
# H. Resampling stability (Nogueira phi): CPSS subsampling + k-fold
# ===========================================================================
cpss  <- rd(DEGX_RUNS, CONTRAST, paste0("resample_cpss_summary_",CONTRAST,".csv"))
kfold <- rd(DEGX_RUNS, CONTRAST, paste0("resample_kfold_summary_",CONTRAST,".csv"))
# metafor two-stage RE arm (CPU side-run; CPSS only) -- spliced in if present
cpss_mf <- rd(DEGX_RUNS, CONTRAST, paste0("resample_cpss_summary_",CONTRAST,"_metafor.csv"))
if (!is.null(cpss) && !is.null(cpss_mf)) cpss <- dplyr::bind_rows(cpss, cpss_mf)
# CPSS and k-fold are rendered as SEPARATE, VERTICAL panels (phi on the y-axis).
# Their phi live on very different scales (CPSS ~0.59-0.76 vs k-fold ~0.96), so a
# shared axis compresses one protocol -- a per-protocol y-axis reads cleanly.
render_stability_vert <- function(df, protocol_lab, cap, outfile, y_lo) {
  st <- df %>% transmute(estimator, phi=stability_phi, lo=phi_ci_lo, hi=phi_ci_hi) %>%
    mutate(fam=fam_resample(estimator), lab=pm(estimator)) %>% arrange(phi)
  st$lab <- factor(st$lab, levels=st$lab)   # ascending phi -> lowest at left
  # phi is a [0,1]-bounded index (1 = perfect). y_lo = per-protocol y-axis floor for a
  # tight, context-rich view: CPSS 0.5 (keeps metafor ~0.59 in frame); k-fold (50-fold) 0.7.
  p <- ggplot(st, aes(lab, phi, color=fam)) +
    geom_errorbar(aes(ymin=lo, ymax=hi), width=0.25, linewidth=0.35, color="grey55") +
    geom_point(size=1.9) +
    scale_color_manual(values=fam_cols, name="family") +
    scale_y_continuous(breaks=seq(0,1,0.1), expand=expansion(mult=c(0.01,0.02))) +
    coord_cartesian(ylim=c(y_lo, 1)) +
    labs(x=NULL, y="stability phi  (-> more reproducible)", caption=prov(cap)) +
    theme_masld()+theme_pub()+
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))
  save_fig(p, file.path(OUT, outfile), width=fig_col_width, height=3.8)
  mark(file.path(OUT, outfile))
}
if (!is.null(cpss)) {
  message("[render] H stability phi -- CPSS", if(!is.null(cpss_mf)) " (+metafor)" else "")
  render_stability_vert(cpss, "CPSS (50% subsampling)", "CPSS 100 pairs", "stability_phi_cpss.pdf", y_lo=0.5)
}
if (!is.null(kfold)) {
  message("[render] H stability phi -- k-fold")
  render_stability_vert(kfold, "k-fold (50-fold)", "k-fold 50 folds", "stability_phi_kfold.pdf", y_lo=0.7)
}

# ===========================================================================
# C2. Bootstrap stable-gene counts (resample bootstrap summary)
# ===========================================================================
boot <- rd(DEGX_RUNS, CONTRAST, paste0("resample_bootstrap_summary_",CONTRAST,".csv"))
boot_mf <- rd(DEGX_RUNS, CONTRAST, paste0("resample_bootstrap_summary_",CONTRAST,"_metafor.csv"))
if (!is.null(boot) && !is.null(boot_mf)) boot <- dplyr::bind_rows(boot, boot_mf)
if (!is.null(boot)) {
  message("[render] C2 stable-gene counts", if(!is.null(boot_mf)) " (+metafor)" else "")
  bs <- boot %>% mutate(fam=fam_resample(estimator), lab=pm(estimator)) %>%
    arrange(n_freq_ge_0p6) %>% mutate(lab=factor(lab, levels=lab))
  pC2 <- ggplot(bs, aes(n_freq_ge_0p6, lab, fill=fam)) +
    geom_col(width=0.72, color="white", linewidth=0.25) +
    geom_point(aes(x=n_ci_excl_0), shape=21, fill="white", color="grey25", size=1.6, stroke=0.4) +
    scale_fill_manual(values=fam_cols, name="family") +
    scale_x_continuous(expand=expansion(mult=c(0,0.05)), labels=scales::comma) +
    labs(x="number of stably-selected genes", y=NULL, caption=prov("1,000 stratified bootstraps")) +
    theme_masld()+theme_pub()+theme(legend.position="right")
  save_fig(pC2, file.path(OUT,"stability_counts.pdf"), width=fig_col_width, height=3.6)
  mark(file.path(OUT,"stability_counts.pdf"))
}

# ===========================================================================
# C. Selection-frequency concordance (resample bootstrap per-gene) -- matrix
# ===========================================================================
rho_csv <- file.path(DEGX_RUNS, CONTRAST, "selfreq_rho.csv")  # precomputed sidecar (see note)
pq <- file.path(DEGX_RUNS, CONTRAST, paste0("resample_bootstrap_pergene_",CONTRAST,".parquet"))
if (file.exists(rho_csv) || file.exists(pq)) {
  message("[render] C selection-freq concordance")
  if (file.exists(rho_csv)) {
    rho <- as.matrix(read.csv(rho_csv, row.names=1, check.names=FALSE))
  } else if (requireNamespace("arrow", quietly=TRUE)) {
    w <- arrow::read_parquet(pq) %>% select(gene, estimator, boot_freq) %>%
      tidyr::pivot_wider(names_from=estimator, values_from=boot_freq)
    m <- as.matrix(w[,-1]); rho <- cor(m, method="spearman", use="pairwise.complete.obs")
  } else rho <- NULL
  if (!is.null(rho)) {
    hc <- hclust(as.dist(1-rho), method="average"); ordm <- rownames(rho)[hc$order]
    long <- as.data.frame(as.table(rho)); names(long) <- c("mi","mj","rho")
    long <- long %>% mutate(mi=factor(pm(mi),levels=pm(ordm)), mj=factor(pm(mj),levels=pm(ordm)))
    pC <- ggplot(long, aes(mi,mj,fill=rho)) +
      geom_tile(color="white", linewidth=0.4) +
      scale_fill_gradientn(colours=green_sc, limits=c(min(rho),1), oob=scales::squish,
                           name="Spearman\nrho") +
      coord_fixed() +
      labs(x=NULL, y=NULL, caption=prov(sprintf("1,000 bootstraps; %d methods, clustered", nrow(rho)))) +
      theme_masld()+theme_pub()+
      theme(axis.text.x=element_text(angle=40,hjust=1), legend.position="right",
            legend.key.width=unit(0.18,"cm"))
    # panelC_selfreq_concordance.pdf retired per user request 2026-06-05 — no longer generated.
    # save_fig(pC, file.path(OUT,"selfreq_concordance.pdf"), width=fig_col_width, height=fig_col_width)
    # mark(file.path(OUT,"selfreq_concordance.pdf"))
  }
}

# ===========================================================================
# I0/I1/I2. R-exact method agreement (real R packages on real data)
# ===========================================================================
manifest <- rd(REX_DIR, paste0("Rexact_",CONTRAST,"_manifest.csv"))
rex_files <- list.files(REX_DIR, pattern=paste0("^Rexact_",CONTRAST,"_.*\\.csv$"), full.names=TRUE)
rex_files <- rex_files[!grepl("manifest", rex_files)]
if (length(rex_files) > 0) {
  message("[render] I0/I1/I2 R-exact agreement (", length(rex_files), " methods)")
  meth <- sub(paste0("^Rexact_",CONTRAST,"_"), "", sub("\\.csv$","", basename(rex_files)))
  # Drop invalid/uninformative metafor variants from the comparison panels:
  #   HK -> near-empty (0/110 DEGs); FE -> assumes a common true effect, invalid
  #   under the high between-cohort heterogeneity seen here (RE is the valid model).
  MF_DROP <- c("metafor_deseq2_hk", "metafor_voom_hk", "metafor_deseq2_fe", "metafor_voom_fe")
  tabs <- setNames(lapply(rex_files, function(f) utils::read.csv(f, stringsAsFactors=FALSE)), meth)
  degset <- lapply(tabs, function(d) if(all(c("gene","padj") %in% names(d))) d$gene[!is.na(d$padj) & d$padj < 0.05] else character(0))

  # I0 -- DEG counts by method (dramatic divergence)
  if (!is.null(manifest)) {
    cnt <- manifest %>% filter(status=="ok", !method %in% MF_DROP) %>%
      mutate(fam=fam_rex(method), lab=pmr(method)) %>%
      arrange(n_sig_padj05) %>% mutate(lab=factor(lab, levels=lab))
    pI0 <- ggplot(cnt, aes(n_sig_padj05, lab, fill=fam)) +
      geom_col(width=0.74, color="white", linewidth=0.25) +
      scale_fill_manual(values=fam_cols_rex, name="family") +
      scale_x_continuous(expand=expansion(mult=c(0,0.06)), labels=scales::comma) +
      labs(x="number of DEGs", y=NULL,
           caption=prov("R/Bioconductor exact fits; metafor FE/HK omitted")) +
      theme_masld()+theme_pub()+theme(legend.position="right")
    save_fig(pI0, file.path(OUT,"rexact_deg_counts.pdf"), width=fig_col_width, height=4.4)
    mark(file.path(OUT,"rexact_deg_counts.pdf"))
  }

  # I1 -- pairwise Jaccard of DEG sets, clustered heatmap (metafor FE/HK dropped)
  meth_j <- setdiff(meth, MF_DROP); ds_j <- degset[meth_j]
  jac <- function(a,b){u<-length(union(a,b)); if(u==0) NA_real_ else length(intersect(a,b))/u}
  M <- outer(seq_along(ds_j), seq_along(ds_j), Vectorize(function(i,j) jac(ds_j[[i]], ds_j[[j]])))
  dimnames(M) <- list(meth_j, meth_j)
  hcj <- hclust(as.dist(1 - M), method="average"); ordj <- rownames(M)[hcj$order]
  lj <- as.data.frame(as.table(M)); names(lj) <- c("mi","mj","jac")
  lj <- lj %>% mutate(mi=factor(pmr(mi),levels=pmr(ordj)), mj=factor(pmr(mj),levels=pmr(ordj)))
  pI1 <- ggplot(lj, aes(mi,mj,fill=jac)) +
    geom_tile(color="white", linewidth=0.3) +
    scale_fill_gradientn(colours=green_sc, limits=c(0,1), na.value="grey90", name="Jaccard") +
    coord_fixed() +
    labs(x=NULL, y=NULL, caption=prov("R/Bioconductor exact fits, incl. ComBat-seq / sva; metafor FE/HK omitted")) +
    theme_masld()+theme_pub()+
    theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="right",
          legend.key.width=unit(0.18,"cm"))
  save_fig(pI1, file.path(OUT,"rexact_jaccard.pdf"), width=6.6, height=6.2)
  mark(file.path(OUT,"rexact_jaccard.pdf"))

  # I2 -- concordance at the top (CAT), two panels: vs dream | vs metafor-RE
  rank_by_p <- function(d) if(all(c("gene","padj") %in% names(d))) d$gene[order(d$padj)] else character(0)
  ranks <- lapply(tabs, rank_by_p)
  ks <- c(50,100,250,500,1000,2000,4000)
  REFS <- c("vs dream"="dream", "vs metafor (RE)"="metafor_voom_re")
  REFS <- REFS[REFS %in% meth]
  if (length(REFS)) {
    cat_df <- do.call(rbind, lapply(names(REFS), function(rn){
      rr <- ranks[[REFS[[rn]]]]; if(!length(rr)) return(NULL)
      do.call(rbind, lapply(setdiff(meth, c(REFS[[rn]], MF_DROP)), function(m){
        r <- ranks[[m]]; if(!length(r)) return(NULL)
        data.frame(ref=rn, method=m, fam=fam_rex(m), k=ks,
                   cat=vapply(ks, function(k) length(intersect(head(r,k), head(rr,k)))/k, 0.0))
      }))
    }))
    pI2 <- ggplot(cat_df, aes(k, cat, group=method, color=fam)) +
      geom_line(linewidth=0.45, alpha=0.85) +
      scale_color_manual(values=fam_cols_rex, name="family") +
      scale_x_log10(labels=scales::comma) + ylim(0,1) +
      facet_wrap(~ref, nrow=1) +
      labs(x="top-k genes (log scale)", y="concordance at the top") +
      theme_masld()+theme_pub()+theme(legend.position="right")
    save_fig(pI2, file.path(OUT,"rexact_cat.pdf"), width=fig_full_width, height=3.4)
    mark(file.path(OUT,"rexact_cat.pdf"))
  }

  # F -- metafor stage-1 engine sensitivity (voom RE vs DESeq2 RE), logFC concordance
  if (all(c("metafor_voom_re","metafor_deseq2_re") %in% meth)) {
    a <- tabs[["metafor_voom_re"]]; b <- tabs[["metafor_deseq2_re"]]
    mg <- merge(a[,c("gene","logFC")], b[,c("gene","logFC")], by="gene", suffixes=c("_voom","_deseq2"))
    rr <- suppressWarnings(cor(mg$logFC_voom, mg$logFC_deseq2, method="spearman", use="complete.obs"))
    lim <- quantile(abs(c(mg$logFC_voom, mg$logFC_deseq2)), 0.999, na.rm=TRUE)
    pF <- ggplot(mg, aes(logFC_voom, logFC_deseq2)) +
      rasterize_layer(geom_point(size=0.25, alpha=0.25, color="#518dc9")) +
      geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=0.3, color="grey45") +
      coord_fixed(xlim=c(-lim,lim), ylim=c(-lim,lim)) +
      annotate("text", x=-lim*0.95, y=lim*0.9, hjust=0, size=PUB_GEOM_TEXT+0.3,
               label=sprintf("Spearman rho = %.3f\nn = %s genes", rr, format(nrow(mg),big.mark=","))) +
      labs(x="meta log2FC (voom stage-1)", y="meta log2FC (DESeq2 stage-1)",
           caption=prov("R-exact metafor RE")) +
      theme_masld()+theme_pub()
    save_fig(pF, file.path(OUT,"metafor_engine_sensitivity.pdf"), width=fig_half_width, height=fig_half_width)
    mark(file.path(OUT,"metafor_engine_sensitivity.pdf"))
  }
}

# ===========================================================================
# B / B2 / D. Power-regime panels (sim_power) -- render when the grid merges
# ===========================================================================
powf <- file.path(DEGX_RUNS, "sim", "sim_power", "batched_summary.csv")
if (!file.exists(powf)) {
  message("[skip] sim_power not merged yet -- panels B (power), B2 (obs-FDR), D (AUROC) pending")
} else {
  message("[render] B/B2/D power-regime panels")
  pw <- utils::read.csv(powf, check.names=FALSE)
  pw <- pw[!pw$estimator %in% c("metafor_fe","metafor_hk"), , drop=FALSE]  # FE invalid / HK degenerate -> RE only
  capg <- sprintf("%s reps/cell x %d methods | 48-cell grid", format(max(pw$n_reps),big.mark=","), dplyr::n_distinct(pw$estimator))

  # B -- power vs effect size, facet by method, line per n (tau2=0.04 slice)
  b <- pw %>% filter(abs(tau2-0.04)<1e-9) %>%
    transmute(lab=pm(estimator), true_lfc, n=factor(n_per_group, levels=sort(unique(n_per_group))), power=power_mean)
  npal <- colorRampPalette(c("#90CAF9","#1565C0","#0D47A1"))(nlevels(b$n))
  pB <- ggplot(b, aes(true_lfc, power, color=n, group=n)) +
    geom_hline(yintercept=0.8, linetype="dotted", linewidth=0.25, color="grey50")+
    geom_line(linewidth=0.5)+geom_point(size=0.9)+scale_color_manual(values=npal,name="n / group")+
    ylim(0,1)+facet_wrap(~lab, ncol=5)+
    labs(x="true log2 fold change", y="power (TPR at FDR < 0.05)", caption=prov(capg))+
    theme_masld()+theme_pub()+theme(legend.position="top", panel.spacing=unit(0.45,"lines"))
  save_fig(pB, file.path(OUT,"power_curves.pdf"), width=fig_full_width, height=4.4); mark(file.path(OUT,"power_curves.pdf"))

  # B2 -- observed FDR per method, diverging calibration bars
  b2 <- pw %>% filter(!is.na(obs_fdr_mean)) %>% group_by(estimator) %>%
    summarise(fdr=mean(obs_fdr_mean), .groups="drop") %>%
    mutate(lab=pm(estimator)) %>% arrange(fdr) %>% mutate(lab=factor(lab,levels=lab))
  pB2 <- ggplot(b2, aes(lab, fdr, fill=fdr)) +
    annotate("rect", xmin=-Inf,xmax=Inf, ymin=NOMINAL, ymax=Inf, fill=masld_colors$up, alpha=0.07) +
    geom_col(width=0.72,color="white",linewidth=0.25)+
    geom_hline(yintercept=NOMINAL, linetype="dashed", linewidth=0.3, color="grey35")+
    cal_fill(limits=sym_lim(b2$fdr), name="observed\nFDR")+
    scale_y_continuous(expand=expansion(mult=c(0,0.05)))+coord_flip()+
    labs(x=NULL, y="observed FDR", caption=prov(capg))+
    theme_masld()+theme_pub()+theme(legend.position="right", legend.key.width=unit(0.18,"cm"))
  save_fig(pB2, file.path(OUT,"fdr_control.pdf"), width=fig_half_width, height=3.6); mark(file.path(OUT,"fdr_control.pdf"))

  # D -- AUROC heatmap (method x effect-size grid, facet by n)
  d <- pw %>% filter(abs(tau2-0.04)<1e-9) %>%
    mutate(n=factor(n_per_group,levels=sort(unique(n_per_group))), lfc=factor(true_lfc,levels=sort(unique(true_lfc))))
  mo <- d %>% group_by(estimator) %>% summarise(a=mean(auroc_mean),.groups="drop") %>% arrange(a) %>% pull(estimator)
  d$lab <- factor(pm(d$estimator), levels=pm(mo))
  pD <- ggplot(d, aes(lfc, lab, fill=auroc_mean))+geom_tile(color="white",linewidth=0.3)+
    scale_fill_gradientn(colours=colorRampPalette(c("#dae7c7","#193c1e"))(12), limits=c(0.5,1), oob=scales::squish, name="AUROC")+
    facet_wrap(~n, nrow=1, labeller=labeller(n=function(x) paste0("n = ",x)))+
    labs(x="true log2 fold change", y=NULL, caption=prov(capg))+
    theme_masld()+theme_pub()+theme(legend.position="right")
  save_fig(pD, file.path(OUT,"auroc_heatmap.pdf"), width=fig_full_width, height=4.0); mark(file.path(OUT,"auroc_heatmap.pdf"))
}

message("\n[degx multimethod] wrote ", length(written), " panel(s) to:\n  ", OUT)
invisible(lapply(written, function(p) cat("   -", basename(p), "\n")))
