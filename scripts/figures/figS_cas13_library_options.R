#!/usr/bin/env Rscript
# =============================================================================
# Figure S_lib_6 -- Library composition options (v8 strategy)
# KEY MESSAGE: under the v8 target strategy the library is assembled from
#   CORE (integrated disease-vs-control DEGs, canonical limma-voom QW C2)
#   UNION COHORT-REPLICATED (DE in >=N of 5 cohorts)
#   UNION MOUSE-CONFIRMED (cross-diet >=3, human-concordant)
#   UNION COLOC (canonical SuSiE PP.H4>0.5)
#   UNION POSITIVE CONTROLS.
# miRNA and the MASH tier are retired. v8 PC gate = mouse-HEPATOCYTE scRNA
# substrate (drop 'absent' or 'ambient_suspect' with hep/other ratio<0.1),
# replacing the v7 bulk-MCD-TPM gate; lncRNA ungated (exploratory arm); positive
# controls + high-COLOC (SuSiE>=0.9) exempt. This figure enumerates the candidate
# scenarios (cohort threshold 2+ vs 3+; human-only vs +mouse +COLOC +controls)
# with gene counts + experimental feasibility, and marks the chosen design inside
# the ~2,000+-500 target band.
#
# Panels:
#   A  Options comparison table (per-scenario tier + biotype counts)
#   B  Composition stacked bars (human-only vs FULL, at 2+ and 3+; band shaded)
#   C  Coverage-vs-mice curves per scenario (10 gRNA/target)
# Output: Cas13_Library_Design/figures/11{a,b,c}_*.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(grid); library(gridExtra); library(scales); library(gtable)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
OUT_DIR <- FIGS_CAS13LIB_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
strip_v <- function(x) sub("[.][0-9]+$", "", x)

# -- Sources (mirror rebuild_cas13_library.R v8) ------------------------------
CANONICAL <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
PERSTUDY  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study")
PERDIET   <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet_cas13")  # Cas13 library Western pool; decoupled from paper 4-model per_diet (2026-06-16)
COLOCFILE <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
POSCTRL   <- file.path(BASE, "results/library/positive_control.csv")
MCDTPM    <- file.path(BASE, "RNA-seq/Mouse/InHouse_MCD/results/mean_tpm_mcd.csv")
ORTHO     <- file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")
META      <- file.path(BASE, "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv")
LIBRARY_CSV <- file.path(BASE, "Cas13_Library_Design/data/cas13_library_v3.0.csv")
COHORTS   <- c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621")
DIETS     <- c("MCD","CDAHFD","Western","HFD")
PADJ <- 0.05; LFC <- 0.5; COLOC_PP4 <- 0.5; MIN_DIETS <- 3L
HEP_RATIO_MIN <- 0.1; COLOC_EXEMPT <- 0.9   # v8 mouse-hep gate + COLOC gate-exemption
MOUSEHEP  <- file.path(BASE, "Cas13_Library_Design/data/mouse_hep_specificity.csv")
CHOSEN_COHORT <- 2L
N_SGRNA <- 10L; N_NT <- 500L; N_EXTRA_CTRL <- 30L
POSCTRL_ALIASES <- c("SCD1" = "SCD"); POSCTRL_EXCLUDE <- c("GIPR", "GLP1R", "FAP")

# Tier colors (controls always grey per FIGURE_GUIDELINES)
tier_colors <- c(`Core` = "#e35070", `Cohort-replicated` = "#f0a050",
                 `Mouse-confirmed` = "#4baeef", `COLOC` = "#8e6fc7",
                 `Positive control` = "#9E9E9E")
biotype_colors <- c(protein_coding = "#5480C2", lncRNA = "#D4834A")

# =============================================================================
# 1. Build the mouse-gene tier sets (same logic as the rebuild)
# =============================================================================
message("Loading sources ...")
mm <- fread(META)
setnames(mm, c("mouse_ensembl_base","mouse_biotype"), c("gb","bt0"), skip_absent = TRUE)
mm[grepl("protein_coding", bt0), bt := "protein_coding"]
mm[grepl("lncRNA|lincRNA", bt0), bt := "lncRNA"]
pc_ids <- mm[bt == "protein_coding", gb]; lnc_ids <- mm[bt == "lncRNA", gb]
pclnc  <- c(pc_ids, lnc_ids)

ortho <- fread(cmd = paste0("zcat ", ORTHO),
               select = c("mouse_ensembl","human_ensembl","human_symbol","confidence_tier","is_one2one"))
ortho <- ortho[confidence_tier %in% c("H","M")]
ortho[, `:=`(hb = strip_v(human_ensembl), mb = strip_v(mouse_ensembl), hsym = toupper(human_symbol))]
ortho[, trank := match(confidence_tier, c("H","M"))]
ortho[, o2o := ifelse(is_one2one %in% c(TRUE,"True","TRUE","true"), 0L, 1L)]
oh <- copy(ortho); setorder(oh, hb, trank, o2o, mb); oh <- unique(oh, by = "hb")
om <- copy(ortho); setorder(om, mb, trank, o2o, hb); om <- unique(om, by = "mb")
osym <- ortho[hsym != ""]; setorder(osym, hsym, trank, o2o, mb); osym <- unique(osym, by = "hsym")
h2m  <- function(hv){ hv <- hv[hv != ""]; intersect(unique(oh[hb %in% hv, mb]), pclnc) }
sym2m<- function(sv){ intersect(unique(osym[hsym %in% toupper(sv), mb]), pclnc) }

can <- fread(CANONICAL, select = c("gene","logFC","padj")); can[, hb := strip_v(gene)]
hlfc <- setNames(can$logFC, can$hb)
core_h <- unique(can[!is.na(padj) & !is.na(logFC) & padj < PADJ & logFC > LFC, hb])
core_mouse <- h2m(core_h)

cohort_lists <- lapply(COHORTS, function(c){
  d <- fread(file.path(PERSTUDY, paste0(c,"_de_results.csv"))); d[, hb := strip_v(gene)]
  unique(d[!is.na(adj.P.Val) & !is.na(logFC) & adj.P.Val < PADJ & logFC > LFC, hb]) })
cohort_n <- table(unlist(cohort_lists))
cohort_mouse <- function(nmin) h2m(names(cohort_n)[cohort_n >= nmin])

de <- lapply(DIETS, function(d) fread(file.path(PERDIET, paste0(d,"_de_results.csv"))))
upd <- lapply(de, function(d){ d[, gb := strip_v(gene)]; d[lfsr < 0.05 & shrunk_logFC > 0.5, gb] })
ndiet <- table(unlist(upd)); mouse3 <- names(ndiet)[ndiet >= MIN_DIETS]
mcm <- om[mb %in% mouse3]; mcm[, hl := hlfc[hb]]
mouse_conf <- intersect(mcm[!is.na(hl) & hl > 0, mb], pclnc)

cl <- fread(COLOCFILE); cl[, hb := strip_v(ensembl)]
coloc_mouse <- h2m(unique(cl[hb != "" & !is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 > COLOC_PP4, hb]))

pcv <- fread(POSCTRL)[["Gene symbol"]]; pcv <- toupper(pcv)
pcv[pcv %in% names(POSCTRL_ALIASES)] <- POSCTRL_ALIASES[pcv[pcv %in% names(POSCTRL_ALIASES)]]
pcv <- setdiff(pcv, POSCTRL_EXCLUDE)                          # drop systemic-mechanism controls
posctrl_mouse <- sym2m(pcv)
scd1_id <- mm[toupper(mouse_symbol_gtf) == "SCD1", gb][1]     # SCD mis-maps to Scd3 -> force Scd1
posctrl_mouse <- unique(c(setdiff(posctrl_mouse, sym2m("SCD")), scd1_id))

# v8 PC screen-expression gate: mouse-HEPATOCYTE scRNA substrate. Drop PC genes that
# are 'absent' or 'ambient_suspect' with hep/other ratio < HEP_RATIO_MIN (ambient
# spillover). lncRNA ungated; positive controls + high-COLOC (SuSiE>=COLOC_EXEMPT) exempt.
mhs <- fread(MOUSEHEP)
id2sym  <- setNames(mm$mouse_symbol_gtf, mm$gb)
sym2sub <- setNames(mhs$mouse_hep_substrate, mhs$gene_symbol)
sym2rat <- setNames(mhs$mouse_hep_ratio, mhs$gene_symbol)
hep_fail <- function(ids) { s <- unname(sym2sub[id2sym[ids]]); r <- unname(sym2rat[id2sym[ids]])
  s[is.na(s)] <- "absent"
  s == "absent" | (s == "ambient_suspect" & (is.na(r) | r < HEP_RATIO_MIN)) }
clp <- fread(COLOCFILE); clp[, hb := strip_v(ensembl)]
clmap <- merge(clp[hb != "" & !is.na(coloc_best_susie_pp4), .(hb, pp4 = coloc_best_susie_pp4)],
               oh[, .(hb, mb)], by = "hb")
cpp <- clmap[, .(pp4 = max(pp4)), by = mb]; coloc_pp4_of <- setNames(cpp$pp4, cpp$mb)
gate_pc <- function(g) {
  pc <- intersect(g, pc_ids); lnc <- intersect(g, lnc_ids)
  exempt <- (pc %in% posctrl_mouse) | (!is.na(coloc_pp4_of[pc]) & coloc_pp4_of[pc] >= COLOC_EXEMPT)
  pc_ok <- pc[!hep_fail(pc) | exempt]
  union(pc_ok, lnc)
}

# priority-assigned tier composition for any scenario (core>cohort>mouse>coloc>control)
scenario <- function(label, nmin, include, chosen = FALSE) {
  coh <- if ("cohort" %in% include) cohort_mouse(nmin) else character(0)
  mo  <- if ("mouse"  %in% include) mouse_conf else character(0)
  co  <- if ("coloc"  %in% include) coloc_mouse else character(0)
  ct  <- if ("control"%in% include) posctrl_mouse else character(0)
  core_g <- core_mouse
  coh_g  <- setdiff(coh, core_g)
  mo_g   <- setdiff(mo,  union(core_g, coh_g))
  co_g   <- setdiff(co,  Reduce(union, list(core_g, coh_g, mo_g)))
  ct_g   <- setdiff(ct,  Reduce(union, list(core_g, coh_g, mo_g, co_g)))
  ungated <- Reduce(union, list(core_g, coh_g, mo_g, co_g, ct_g))
  allg <- gate_pc(ungated)                                   # PC mouse-hepatocyte gate (lncRNA ungated)
  ing <- function(s) length(intersect(s, allg))
  data.table(scenario = label, cohort_min = nmin,
             Core = ing(core_g), `Cohort-replicated` = ing(coh_g),
             `Mouse-confirmed` = ing(mo_g), COLOC = ing(co_g), `Positive control` = ing(ct_g),
             n_total = length(allg), n_ungated = length(ungated),
             n_pc = length(intersect(allg, pc_ids)), n_lnc = length(intersect(allg, lnc_ids)),
             chosen = chosen, gene_set = list(allg))
}
opts <- rbindlist(list(
  scenario("Core only",            2L, c("core")),
  scenario("Core ∪ 2+ cohort",     2L, c("core","cohort")),
  scenario("+ Mouse-confirmed",    2L, c("core","cohort","mouse")),
  scenario("+ COLOC",              2L, c("core","cohort","mouse","coloc")),
  scenario("FULL (2+)",            2L, c("core","cohort","mouse","coloc","control"), chosen = TRUE),
  scenario("FULL (3+)",            3L, c("core","cohort","mouse","coloc","control"))
), fill = TRUE)
opts[, n_sgrna := n_total * N_SGRNA + N_NT + N_EXTRA_CTRL * N_SGRNA]

# feasibility (delivery/sort/survival; PI 2026-06-01)
HEP<-1e7; TRANSD<-0.30; CRE_EFF<-0.80; FACS_EFF<-0.60; IN_FRAC<-0.13; GATE<-0.15; MORTALITY<-1/8
eff_cells <- min(HEP*TRANSD*CRE_EFF*IN_FRAC, HEP*TRANSD*CRE_EFF*(1-IN_FRAC)*GATE*FACS_EFF)
eff_inject<- eff_cells*(1-MORTALITY)
opts[, surv_mice := ceiling(n_sgrna*500/eff_cells)]
opts[, min_mice  := ceiling(surv_mice/(1-MORTALITY))]

cat("\n=== v8 Library Options ===\n")
print(opts[, .(scenario, cohort_min, n_total, Core, `Cohort-replicated`,
               `Mouse-confirmed`, COLOC, `Positive control`, n_pc, n_lnc, n_sgrna, min_mice, chosen)])

# canonical check vs the rebuilt CSV
if (file.exists(LIBRARY_CSV)) {
  lib <- fread(LIBRARY_CSV); chosenN <- opts[chosen == TRUE, n_total]
  cat(sprintf("CANONICAL CHECK: chosen scenario total=%d vs cas13_library_v3.0.csv rows=%d (%s)\n",
              chosenN, nrow(lib), ifelse(chosenN == nrow(lib), "MATCH", "MISMATCH")))
}

# =============================================================================
# PANEL A: options table
# =============================================================================
tbl <- data.table(
  Composition = opts$scenario,
  `Cohort` = paste0(">=", opts$cohort_min),
  `Pre-gate` = comma(opts$n_ungated),
  Total    = comma(opts$n_total),
  Core     = comma(opts$Core),
  `Cohort\nrep.` = comma(opts$`Cohort-replicated`),
  `Mouse\nconf.` = comma(opts$`Mouse-confirmed`),
  COLOC    = comma(opts$COLOC),
  `Pos.\nctrl` = comma(opts$`Positive control`),
  PC = comma(opts$n_pc), lncRNA = comma(opts$n_lnc),
  sgRNAs = comma(opts$n_sgrna), `Mice\n(500x)` = opts$min_mice)
tt <- ttheme_minimal(base_size = 7, base_family = "Helvetica",
  core = list(bg_params = list(fill = ifelse(opts$chosen, "#FDECEF", c("gray97","white")),
                               col = "gray80", lwd = 0.4),
              fg_params = list(fontsize = 7, hjust = 0.5, x = 0.5)),
  colhead = list(bg_params = list(fill = "#E8E8E8", col = "gray60", lwd = 0.5),
                 fg_params = list(fontsize = 7, fontface = "bold", hjust = 0.5, x = 0.5)))
tg <- tableGrob(tbl, rows = NULL, theme = tt)
tg <- gtable_add_grob(tg, rectGrob(gp = gpar(lwd = 1.0, col = "gray40", fill = NA)),
                      t = 1, b = nrow(tg), l = 1, r = ncol(tg))
panel_A <- wrap_elements(full = tg) + ggtitle("A") +
  labs(caption = paste0("v8 strategy: Core (integrated C2, padj<0.05 & logFC>0.5, UP) ",
        "∪ Cohort-replicated (DE in >=N of 5 cohorts) ∪ Mouse-confirmed (>=3 diets, concordant) ",
        "∪ COLOC (SuSiE PP.H4>0.5) ∪ positive controls. PC+lncRNA, no miRNA. ",
        "Pre-gate = ungated union; Total = after PC mouse-hepatocyte gate (drop 'absent' or ",
        "ambient_suspect with hep/other ratio<", HEP_RATIO_MIN, "; lncRNA ungated; positive controls + COLOC>=",
        COLOC_EXEMPT, " exempt). Library v8 = FULL (2+) = ", comma(opts[chosen == TRUE, n_total]),
        " genes. ", N_SGRNA, " gRNA/target; mice incl. 12.5% mortality.")) +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0),
        plot.caption = element_text(size = 5, color = "gray40", hjust = 0))

# =============================================================================
# PANEL B: composition stacked bars (human-only vs FULL, 2+ and 3+) + band
# =============================================================================
barset <- rbindlist(list(
  scenario("Human-only\n(2+)",    2L, c("core","cohort")),
  scenario("Human-only\n(3+)",    3L, c("core","cohort")),
  scenario("FULL (2+)",           2L, c("core","cohort","mouse","coloc","control"), chosen = TRUE),
  scenario("FULL\n(3+)",          3L, c("core","cohort","mouse","coloc","control"))
), fill = TRUE)
tier_lvls <- c("Positive control","COLOC","Mouse-confirmed","Cohort-replicated","Core")
comp <- melt(barset[, c("scenario", tier_lvls), with = FALSE], id.vars = "scenario",
             variable.name = "tier", value.name = "n")
comp[, scenario := factor(scenario, levels = barset$scenario)]
comp[, tier := factor(tier, levels = tier_lvls)]
tot <- barset[, .(scenario = factor(scenario, levels = barset$scenario), n_total, chosen)]
panel_B <- ggplot(comp, aes(scenario, n, fill = tier)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1500, ymax = 2500,
           fill = "#2E9A86", alpha = 0.10) +
  annotate("text", x = -Inf, y = 2500, label = "target ~2,000 +- 500", hjust = -0.05, vjust = 1.4,
           size = 2.0, color = "#2E9A86", fontface = "italic") +
  geom_col(width = 0.66, color = "white", linewidth = 0.2) +
  geom_text(data = tot, aes(scenario, n_total, label = comma(n_total), fill = NULL,
            fontface = ifelse(chosen, "bold", "plain")), vjust = -0.4, size = 2.5,
            inherit.aes = FALSE) +
  scale_fill_manual(values = tier_colors, name = NULL, breaks = rev(tier_lvls)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Library target genes", title = "B") +
  theme_masld() + theme_pub() +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0),
        legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6), axis.text.x = element_text(size = 6.2))

# =============================================================================
# PANEL C: coverage-vs-mice curves
# =============================================================================
COV_CAP <- 1000L; COV_LINES <- c(100L,250L,500L,750L)
opts[, mice_at_cap := ceiling(COV_CAP*n_sgrna/eff_inject)]
max_mice <- max(opts$mice_at_cap)+3; mice_range <- 5:max_mice
cc <- rbindlist(lapply(1:nrow(opts), function(i)
  data.table(scenario = opts$scenario[i], mice = mice_range,
             coverage = (mice_range*eff_inject)/opts$n_sgrna[i])))
cc[, scenario := factor(scenario, levels = opts$scenario)]
pts <- opts[, .(scenario = factor(scenario, levels = opts$scenario), min_mice, n_sgrna, chosen)]
pts[, cov_at_min := (min_mice*eff_inject)/n_sgrna]
setorder(pts, min_mice); pts[, label_y := 500 + 50 + (seq_len(.N)-1)*60]
opt_colors <- setNames(scales::hue_pal()(nrow(opts)), opts$scenario)
opt_colors[opts[chosen==TRUE, scenario]] <- "#C0143C"
panel_C <- ggplot(cc, aes(mice, coverage, color = scenario)) +
  geom_hline(yintercept = 500, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_text(data = data.table(y = COV_LINES, lab = paste0(COV_LINES,"x")), inherit.aes = FALSE,
            aes(x = 3, y = y, label = lab), size = 2.2, color = "gray45", hjust = 0, vjust = -0.3) +
  geom_line(linewidth = 0.7) +
  geom_point(data = pts, aes(min_mice, cov_at_min, color = scenario), size = 2, show.legend = FALSE) +
  geom_text(data = pts, aes(min_mice, label_y, color = scenario,
            label = paste0(min_mice, " mice")), size = 2.4, vjust = -0.3, show.legend = FALSE) +
  scale_color_manual(values = opt_colors, name = NULL) +
  scale_y_continuous(breaks = c(COV_LINES, COV_CAP), labels = comma) +
  coord_cartesian(ylim = c(0, COV_CAP)) +
  labs(x = "Mice to inject (for 500x coverage)", y = "Coverage per sgRNA", title = "C",
       caption = sprintf("%s hep/mouse x %.0f%% transduction x %.0f%% Cre = %s informative; binding pool ~%s cells/mouse; %d gRNA/target.",
                         comma(HEP), TRANSD*100, CRE_EFF*100, comma(HEP*TRANSD*CRE_EFF), comma(round(eff_cells)), N_SGRNA)) +
  theme_masld() + theme_pub() +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0),
        plot.caption = element_text(size = 4.8, color = "gray50", hjust = 0),
        legend.position = "right", legend.text = element_text(size = 5.5),
        legend.key.size = unit(0.3, "cm"))

# =============================================================================
# Save panels
# =============================================================================
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
panels <- list(`11a_library_options_table` = list(p = panel_A, w = 13.5, h = 4.2),
               `11b_library_composition`   = list(p = panel_B, w = fig_half_width + 0.5, h = 3.8),
               `11c_coverage_curves`       = list(p = panel_C, w = 9.0, h = 5.2))
for (nm in names(panels)) {
  ggsave(file.path(OUT_DIR, paste0(nm, ".pdf")), panels[[nm]]$p,
         width = panels[[nm]]$w, height = panels[[nm]]$h, device = pdf_device)
  message("Saved: ", nm)
}
cat("\nfigS_cas13_library_options.R complete.\n")
