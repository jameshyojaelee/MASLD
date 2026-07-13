#!/usr/bin/env Rscript
# Curated SUPPLEMENTARY isoform figures (post-review, per-cohort, no dream/pooling,
# no LIDPAD). Each panel has a KEY MESSAGE. Guideline-compliant: PDF, control=gray,
# bars->lollipops, horizontal sorted bars, named palettes. Run in rnaseq env.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
grDevices::pdf.options(useDingbats = FALSE)
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(PROJ, "scripts/figures/publication_theme.R"))
HRES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/human")
MRES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/mouse")
OUT  <- file.path(PROJ, "figures/supplementary/figS_isoform_diversity")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
pp  <- function(n) file.path(OUT, sprintf("figS_isoform_%s.pdf", n))
GRAY <- "#9E9E9E"; HIL <- palette1[1]
panel <- function(name, expr) tryCatch({ expr; cat("  [ok]", name, "\n") },
                                       error=function(e) cat("  [SKIP]", name, ":", conditionMessage(e), "\n"))

# ============================================================================
# A. Dominant-isoform usage rises as DE fold-change rises (per-cohort, reproducible)
#    KEY MESSAGE: strongly differentially-expressed genes are single-isoform-dominated.
# ============================================================================
panel("A_dominant_iso_vs_lfc", {
  d <- fread(file.path(HRES, "percohort_diversity_vs_lfc.tsv"))
  d[, absLFC := abs(lfc)]; d[, DIFpct := DIF_dis*100]        # % of gene that is its top isoform
  d[, bin := cut(absLFC, c(0,0.25,0.5,1,2,Inf), c("0–0.25","0.25–0.5","0.5–1","1–2",">2"), right=FALSE)]
  b <- d[!is.na(bin), .(m=mean(DIFpct,na.rm=TRUE), se=sd(DIFpct,na.rm=TRUE)/sqrt(.N)), by=.(cohort,bin)]
  cols <- c(GSE213621=HIL, GSE135251="#4baeef", GSE130970="#30d796")
  p <- ggplot(b, aes(bin, m, colour=cohort, group=cohort)) +
    geom_ribbon(aes(ymin=m-se, ymax=m+se, fill=cohort), alpha=.15, colour=NA) +
    geom_line(linewidth=.7) + geom_point(size=1.8) +
    scale_colour_manual(values=cols, name=NULL) + scale_fill_manual(values=cols, guide="none") +
    labs(x="|log2 fold-change| (per-cohort DE)", y="dominant isoform (% of gene's transcripts)",
         title="Bigger DE = more single-isoform-dominated") +
    theme_masld() + theme(legend.position=c(.72,.25), axis.text.x=element_text(angle=30,hjust=1))
  save_fig(p, pp("A_dominant_iso_vs_lfc"), width=3.4, height=3)
})

# ============================================================================
# B. Per-cohort DTU power/calibration — only GSE213621 is usable
#    KEY MESSAGE: of 3 human cohorts only GSE213621 supports a credible switch test.
# ============================================================================
panel("B_percohort_dtu", {
  dg <- fread(file.path(HRES, "dtu_null_diagnostics.tsv"))
  cd <- dg[grepl("cohort-GSE", stratum)]
  cd[, cohort := sub(".*cohort-", "", stratum)]
  # only well-calibrated (null_status=="ok") strata yield trustworthy switch counts;
  # non-calibrated strata set to 0 (their raw counts are empirical-null-collapse artifacts)
  cd[, switches := fifelse(null_status == "ok", emp_effect_genes, 0L)]
  p <- ggplot(cd, aes(reorder(cohort, switches), switches)) +
    geom_col(fill=HIL, width=.62) +
    coord_flip() + scale_y_continuous(expand=expansion(mult=c(0,.08))) +
    labs(x=NULL, y="effect-gated isoform switches (|Δprop|≥0.1)",
         title="Per-cohort isoform-switch detection") +
    theme_masld()
  save_fig(p, pp("B_percohort_dtu"), width=3.8, height=2.6)
})

# ============================================================================
# C. GSE213621 switches: top by effect size, coloured by MANE-switch
#    KEY MESSAGE: 91 switches; 67 (74%) shift away from the canonical MANE isoform.
# ============================================================================
panel("C_gse213621_switches", {
  t <- fread(file.path(HRES, "gse213621_dtu_targets.tsv"))
  t <- t[order(-disease_dprop)][1:20]
  t[, mane := fifelse(switch_away_from_mane, "switch away from MANE", "MANE retained")]
  p <- ggplot(t, aes(reorder(symbol, disease_dprop), disease_dprop, colour=mane)) +
    geom_segment(aes(xend=symbol, y=0, yend=disease_dprop), linewidth=.6) +
    geom_point(size=2) + coord_flip() +
    scale_colour_manual(values=c(`switch away from MANE`=HIL, `MANE retained`=GRAY), name=NULL) +
    labs(x=NULL, y="Δ isoform proportion (disease − control)",
         title="GSE213621 isoform switches (top 20 of 91)",
         subtitle="67/91 (74%) switch away from the MANE canonical isoform") +
    theme_masld() + theme(legend.position="top", axis.text.y=element_text(size=6))
  save_fig(p, pp("C_gse213621_switches"), width=3.6, height=4)
})

# ============================================================================
# D. Mouse targetability funnel — human signal -> Cas13-actionable mouse target
#    KEY MESSAGE: only 8 of 67 human switch genes reach isoform-selective mouse targetability.
# ============================================================================
panel("D_cas13_target_isoforms", {
  m <- fread(file.path(MRES, "gse213621_targets_mouse_isoform_structure.tsv"))
  lv <- c("Isoform switches (human)",
          "Mouse orthologs",
          "Expressed in mouse liver (≥1 TPM)",
          "≥2 expressed mouse isoforms (≥0.5 TPM)")
  steps <- data.table(
    stage = factor(lv, levels=rev(lv)),
    n = c(nrow(m), sum(m$ortholog!="none"),
          sum(!is.na(m$mouse_n_expr_iso) & m$mouse_n_expr_iso>0),
          sum(m$isoform_selective_possible==TRUE, na.rm=TRUE)))
  steps[, hl := stage==lv[4]]
  p <- ggplot(steps, aes(stage, n, fill=hl)) +
    geom_col(width=.7) + geom_text(aes(label=n), hjust=-0.3, size=3) + coord_flip() +
    scale_fill_manual(values=c(`TRUE`=HIL, `FALSE`=GRAY), guide="none") +
    scale_y_continuous(expand=expansion(mult=c(0,.15))) +
    labs(x=NULL, y="genes", title="Cas13 target isoforms") +
    theme_masld()
  save_fig(p, pp("D_cas13_target_isoforms"), width=5.2, height=2.4)
})

# ============================================================================
# E. Cross-species: no conserved isoform switching (honest null)
#    KEY MESSAGE: isoform switching is species-specific (observed ≈ expected, NS).
# ============================================================================
panel("E_crossspecies_null", {
  st <- fread(file.path(PROJ,"RNA-seq/results/isoform_diversity/cross_species_concordance_stats.tsv"))
  u <- st[mapping=="strict_1to1" & test=="union"]
  d <- data.table(cat=factor(c("expected by chance","observed"), levels=c("expected by chance","observed")),
                  n=c(u$expected, u$k_obs))
  p <- ggplot(d, aes(cat, n, fill=cat)) + geom_col(width=.55) +
    geom_text(aes(label=round(n,1)), vjust=-.3, size=3) +
    scale_fill_manual(values=c(observed=HIL, `expected by chance`=GRAY), guide="none") +
    scale_y_continuous(expand=expansion(mult=c(0,.2))) +
    labs(x=NULL, y="mouse switches with human-DTU 1:1 ortholog",
         title="No conserved cross-species switching",
         subtitle=sprintf("%d obs vs %.1f exp; %.2fx, p=%.2f (NS)", u$k_obs,u$expected,u$fold,u$hyper_p)) +
    theme_masld()
  save_fig(p, pp("E_crossspecies_null"), width=2.8, height=3)
})
cat("[supp] curated panels ->", OUT, "\n")
