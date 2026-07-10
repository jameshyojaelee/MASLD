#!/usr/bin/env Rscript
# =============================================================================
# Figure S_lib_6 -- Library composition options (v9 strategy)
# KEY MESSAGE: the v9 library is the union of CORE (integrated disease-vs-control
#   DEGs, canonical limma-voom QW C2), COHORT-REPLICATED (DE in >=N of 5 cohorts),
#   MOUSE-REPLICATED (cross-diet >=3, human-concordant), COLOC (SuSiE PP.H4>0.5),
#   and POSITIVE CONTROLS. miRNA and the MASH tier are retired.
# v9 expression gate = mouse-HEPATOCYTE scRNA pseudobulk substrate (CPM, NOT TPM): drop
#   PC <1.0 CPM and lncRNA <0.1 CPM in mouse hepatocytes; the hep/other ratio is not
#   used; positive controls + high-COLOC (SuSiE>=0.9) exempt. A final pool-guideability
#   gate drops genes with no designable Cas13 guide in the vM38 pool (hard constraint,
#   all tiers). PCGs use TREAT lfc=0.5 on EVERY axis (C1/C2/C3 + cohort), biotype-split;
#   lncRNAs use no LFC floor (TREAT lfc=0).
#     Option A = cohort 2+ (chosen); Option B = cohort 3+.
#   Table/bars count POSITIVE CONTROLS FIRST, then de-dup core > cohort > mouse > COLOC.
#
# Panels:
#   A  Options table (Option A/B; evidence-tier counts | biotype + experiment categories)
#   B  Composition stacked bars (Option A vs Option B, by tier)
#   C  Coverage-vs-mice curves per scenario (6 gRNA/target)
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

# -- Sources (mirror rebuild_cas13_library.R v9) ------------------------------
CANONICAL <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
MASHMASL  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/mash_vs_masl_dream.csv")
ADVFIB    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/adv_vs_early_fibrosis_dream.csv")
PERSTUDY  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study")
PERDIET   <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet_cas13")  # Cas13 library Western pool; decoupled from paper 4-model per_diet (2026-06-16)
COLOCFILE <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
POSCTRL   <- file.path(BASE, "results/library/positive_control.csv")
MCDTPM    <- file.path(BASE, "RNA-seq/Mouse/InHouse_MCD/results/mean_tpm_mcd.csv")
ORTHO     <- file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")
META      <- file.path(BASE, "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv")
LIBRARY_CSV <- file.path(BASE, "Cas13_Library_Design/data/cas13_library.csv")
COHORTS   <- c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621")
DIETS     <- c("MCD","CDAHFD","Western","HFD")
PADJ <- 0.05; LFC_PC <- 0.2; LFC_LNC <- 0.0; COLOC_PP4 <- 0.5; MIN_DIETS <- 3L   # PC LFC 0.2 (PI sizing 2026-06-29, ~2,177 lib)
COLOC_EXEMPT <- 0.9   # COLOC gate-exemption (>=this PP.H4 kept even if hep-absent)
HEP_CPM_PC <- 1.0; HEP_CPM_LNC <- 0.1   # mouse-hep pseudobulk CPM floors (PC / lncRNA)
MOUSEHEP  <- file.path(BASE, "Cas13_Library_Design/data/mouse_hep_specificity_vm38.csv")
CHOSEN_COHORT <- 2L
N_SGRNA <- 4L; N_NT <- 500L; N_EXTRA_CTRL <- 30L   # 4 guides/target (matches build_library_guides.py --n 4, 2026-06-29)
POSCTRL_ALIASES <- c("SCD1" = "SCD"); POSCTRL_EXCLUDE <- c("GIPR", "GLP1R", "FAP")
infer_df_total <- function(dt) {
  pr <- dt[is.finite(t) & is.finite(P.Value) & P.Value > 0 & P.Value < 1 & abs(t) > 1e-6]
  idx <- unique(round(seq(1, nrow(pr), length.out = min(nrow(pr), 12))))
  median(vapply(idx, function(i)
    uniroot(function(df) 2 * pt(-abs(pr$t[i]), df = df) - pr$P.Value[i], c(0.1, 1e6))$root,
    numeric(1)))
}
add_treat_fdr <- function(dt, lfc, se_col = NULL, df_col = NULL) {
  out <- copy(dt)
  se <- if (!is.null(se_col) && se_col %in% names(out)) out[[se_col]] else abs(out$logFC / out$t)
  se[!is.finite(se) | se <= 0] <- NA_real_
  df_use <- if (!is.null(df_col) && df_col %in% names(out)) out[[df_col]] else infer_df_total(out)
  out[, p_treat := pt((abs(logFC) - lfc) / se, df = df_use, lower.tail = FALSE) +
                   pt((abs(logFC) + lfc) / se, df = df_use, lower.tail = FALSE)]
  out[, fdr_treat := p.adjust(p_treat, method = "BH")]
  out
}
# Biotype-aware TREAT (mirrors rebuild_cas13_library.R): PC genes tested against the PC
# effect-size offset, lncRNAs against no floor. `dt` must carry a `bt` column.
add_treat_fdr_bt <- function(dt, se_col = NULL, df_col = NULL,
                             lfc_pc = LFC_PC, lfc_lnc = LFC_LNC) {
  pc  <- add_treat_fdr(dt[bt != "lncRNA"], lfc = lfc_pc,  se_col = se_col, df_col = df_col)
  lnc <- add_treat_fdr(dt[bt == "lncRNA"], lfc = lfc_lnc, se_col = se_col, df_col = df_col)
  rbindlist(list(pc, lnc), fill = TRUE)
}

# Tier colors (controls always grey per FIGURE_GUIDELINES)
tier_colors <- c(`Core` = "#e35070", `Cohort-replicated` = "#f0a050",
                 `Mouse-replicated` = "#4baeef", `COLOC` = "#8e6fc7",
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
bt_of <- setNames(mm$bt, mm$gb)

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

can <- fread(CANONICAL, select = c("gene","logFC","SE","t","P.Value","padj","symbol")); can[, hb := strip_v(gene)]
can_m <- merge(can, oh[, .(hb, mb)], by = "hb", all.x = TRUE)
can_m[, bt := bt_of[mb]]
can_m[is.na(bt), bt := "protein_coding"]

# UP-significant human genes for one axis under biotype-split TREAT (PC at LFC_PC, lnc at 0)
sel_axis <- function(dt, se_col = NULL, df_col = NULL, lnc_lfc = LFC_LNC) {
  x <- add_treat_fdr_bt(dt, se_col = se_col, df_col = df_col, lfc_lnc = lnc_lfc)
  unique(x[!is.na(logFC) & fdr_treat < PADJ & logFC > 0, .(hb, hlfc = logFC)])
}
# CORE = C1 integrated disease-vs-control DEGs ONLY (MASH-vs-MASL + fibrosis axes dropped
# 2026-06-29). The human arm = CORE integrated UNION COHORT-REPLICATED (cohort tier below).
core_mouse_for <- function(lnc_lfc = LFC_LNC) {
  s1 <- sel_axis(copy(can_m), se_col = "SE", lnc_lfc = lnc_lfc)
  h2m(unique(s1$hb))
}
core_mouse <- core_mouse_for(LFC_LNC)
human_support <- sel_axis(copy(can_m), se_col = "SE")
human_support <- human_support[, .(hlfc = max(hlfc, na.rm = TRUE)), by = hb]
hlfc <- setNames(human_support$hlfc, human_support$hb)

perstudy <- lapply(COHORTS, function(c){
  d <- fread(file.path(PERSTUDY, paste0(c,"_de_results.csv")))
  d[, hb := strip_v(gene)]
  dm <- merge(d, oh[, .(hb, mb)], by = "hb", all.x = TRUE)
  dm[, bt := bt_of[mb]]
  dm[is.na(bt), bt := "protein_coding"]
  add_treat_fdr_bt(dm, se_col = "SE", df_col = "df.total")   # PC at LFC_PC, lncRNA at 0
})
names(perstudy) <- COHORTS

cohort_n_for <- function(lnc_lfc = LFC_LNC) {
  cohort_lists <- lapply(perstudy, function(dm)
    unique(dm[!is.na(logFC) & !is.na(fdr_treat) & fdr_treat < PADJ & logFC > 0, hb]))
  table(unlist(cohort_lists))
}
cohort_mouse_for <- function(lnc_lfc = LFC_LNC, nmin) {
  cohort_n <- cohort_n_for(lnc_lfc)
  h2m(names(cohort_n)[cohort_n >= nmin])
}
cohort_mouse <- function(nmin) cohort_mouse_for(LFC_LNC, nmin)

de <- lapply(DIETS, function(d) fread(file.path(PERDIET, paste0(d,"_de_results.csv"))))
mouse_conf_for <- function(lnc_lfc = LFC_LNC) {
  upd <- lapply(de, function(d){
    d[, gb := strip_v(gene)]
    d[, bt := bt_of[gb]]
    d[is.na(bt), bt := "protein_coding"]
    d_pc <- add_treat_fdr(d[bt != "lncRNA"], lfc = LFC_PC, df_col = "df_total")
    d_lnc <- add_treat_fdr(d[bt == "lncRNA"], lfc = lnc_lfc, df_col = "df_total")
    rbindlist(list(d_pc, d_lnc), fill = TRUE)[fdr_treat < PADJ & logFC > 0, gb]
  })
  ndiet <- table(unlist(upd)); mouse3 <- names(ndiet)[ndiet >= MIN_DIETS]
  mcm <- om[mb %in% mouse3]; mcm[, hl := hlfc[hb]]
  intersect(mcm[!is.na(hl) & hl > 0, mb], pclnc)
}
mouse_conf <- mouse_conf_for(LFC_LNC)

cl <- fread(COLOCFILE); cl[, hb := strip_v(ensembl)]
coloc_mouse <- h2m(unique(cl[hb != "" & !is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 > COLOC_PP4, hb]))

pcv <- fread(POSCTRL)[["Gene symbol"]]; pcv <- toupper(pcv)
pcv[pcv %in% names(POSCTRL_ALIASES)] <- POSCTRL_ALIASES[pcv[pcv %in% names(POSCTRL_ALIASES)]]
pcv <- setdiff(pcv, POSCTRL_EXCLUDE)                          # drop systemic-mechanism controls
posctrl_mouse <- sym2m(pcv)
scd1_id <- mm[toupper(mouse_symbol_gtf) == "SCD1", gb][1]     # SCD mis-maps to Scd3 -> force Scd1
posctrl_mouse <- unique(c(setdiff(posctrl_mouse, sym2m("SCD")), scd1_id))

# Mouse-hepatocyte scRNA pseudobulk CPM = the v9 gate substrate (mirrors rebuild).
mhs <- fread(MOUSEHEP)
hep_cpm_of <- setNames(mhs$mouse_hep_cpm, strip_v(mhs$gene_id))

clp <- fread(COLOCFILE); clp[, hb := strip_v(ensembl)]
clmap <- merge(clp[hb != "" & !is.na(coloc_best_susie_pp4), .(hb, pp4 = coloc_best_susie_pp4)],
               oh[, .(hb, mb)], by = "hb")
cpp <- clmap[, .(pp4 = max(pp4)), by = mb]; coloc_pp4_of <- setNames(cpp$pp4, cpp$mb)

# Screen-expression gate: mouse-hepatocyte CPM (PC>=1.0, lncRNA>=0.1; CPM, NOT TPM);
# positive controls + high-COLOC (>=COLOC_EXEMPT) exempt. Mirrors rebuild_cas13_library.R.
gate_pc <- function(g) {
  pc <- intersect(g, pc_ids); lnc <- intersect(g, lnc_ids)
  pc_cpm <- unname(hep_cpm_of[pc]); pc_cpm[is.na(pc_cpm)] <- 0
  lnc_cpm <- unname(hep_cpm_of[lnc]); lnc_cpm[is.na(lnc_cpm)] <- 0

  exempt_g <- function(ids) {
    (ids %in% posctrl_mouse) | (!is.na(coloc_pp4_of[ids]) & coloc_pp4_of[ids] >= COLOC_EXEMPT)
  }

  pc_ok <- pc[pc_cpm >= HEP_CPM_PC | exempt_g(pc)]
  lnc_ok <- lnc[lnc_cpm >= HEP_CPM_LNC | exempt_g(lnc)]

  union(pc_ok, lnc_ok)
}

# pool-guideability gate (mirrors rebuild_cas13_library.R 2026-06-24): drop genes with
# NO designable Cas13 guide in the upstream vM38 pool (single-exon / paralog-family /
# mitochondrial / wrong-biotype / no-CCDS). HARD constraint -- applies to every tier
# including pos-ctrl and high-COLOC. Subtracted from every scenario below.
UNGUIDEABLE <- file.path(BASE, "Cas13_Library_Design/data/guides/unguideable_vM38.csv")
ung_ids <- if (file.exists(UNGUIDEABLE)) strip_v(fread(UNGUIDEABLE)$gene_id_mouse) else character(0)

# priority-assigned tier composition for any scenario (core>cohort>mouse>coloc>control)
scenario <- function(label, nmin, include, chosen = FALSE, lnc_lfc = LFC_LNC) {
  core_set <- core_mouse_for(lnc_lfc)
  coh <- if ("cohort" %in% include) cohort_mouse_for(lnc_lfc, nmin) else character(0)
  mo  <- if ("mouse"  %in% include) mouse_conf_for(lnc_lfc) else character(0)
  co  <- if ("coloc"  %in% include) coloc_mouse else character(0)
  ct  <- if ("control"%in% include) posctrl_mouse else character(0)
  # dedup priority: Positive control FIRST, then Core > Cohort-rep > Mouse-rep > COLOC
  # (controls are a deliberate set, attributed first so the column shows the true count)
  ct_g   <- ct
  core_g <- setdiff(core_set, ct_g)
  coh_g  <- setdiff(coh, union(ct_g, core_g))
  mo_g   <- setdiff(mo,  Reduce(union, list(ct_g, core_g, coh_g)))
  co_g   <- setdiff(co,  Reduce(union, list(ct_g, core_g, coh_g, mo_g)))
  ungated <- Reduce(union, list(ct_g, core_g, coh_g, mo_g, co_g))
  allg <- setdiff(gate_pc(ungated), ung_ids)                 # PC mouse-hep gate (lncRNA ungated) + pool-guideability gate
  ing <- function(s) length(intersect(s, allg))
  data.table(scenario = label, cohort_min = nmin,
             lnc_lfc = lnc_lfc,
             Core = ing(core_g), `Cohort-replicated` = ing(coh_g),
             `Mouse-replicated` = ing(mo_g), COLOC = ing(co_g), `Positive control` = ing(ct_g),
             # NON-deduplicated: genes in the final library meeting EACH tier's criterion,
             # counted in every tier they qualify for (overlapping; sum > n_total).
             Core_all = ing(core_set), Cohort_all = ing(coh), Mouse_all = ing(mo),
             COLOC_all = ing(co), Control_all = ing(ct),
             n_total = length(allg), n_ungated = length(ungated),
             n_pc = length(intersect(allg, pc_ids)), n_lnc = length(intersect(allg, lnc_ids)),
             chosen = chosen, gene_set = list(allg))
}
opts <- rbindlist(list(
  scenario("Core only",            2L, c("core")),
  scenario("Core ∪ 2+ cohort",     2L, c("core","cohort")),
  scenario("+ Mouse-replicated",    2L, c("core","cohort","mouse")),
  scenario("+ COLOC",              2L, c("core","cohort","mouse","coloc")),
  scenario("Option A",            2L, c("core","cohort","mouse","coloc","control"), chosen = TRUE),
  scenario("Option B",            3L, c("core","cohort","mouse","coloc","control"))
), fill = TRUE)
opts[, n_sgrna := n_total * N_SGRNA + N_NT + N_EXTRA_CTRL * N_SGRNA]

# feasibility (delivery/sort/survival; PI 2026-06-01)
HEP<-1e7; TRANSD<-0.30; CRE_EFF<-0.80; FACS_EFF<-0.60; IN_FRAC<-0.13; GATE<-0.15; MORTALITY<-1/8
eff_cells <- min(HEP*TRANSD*CRE_EFF*IN_FRAC, HEP*TRANSD*CRE_EFF*(1-IN_FRAC)*GATE*FACS_EFF)
eff_inject<- eff_cells*(1-MORTALITY)
opts[, surv_mice := ceiling(n_sgrna*500/eff_cells)]
opts[, min_mice  := ceiling(surv_mice/(1-MORTALITY))]

cat("\n=== v9 Library Options ===\n")
print(opts[, .(scenario, cohort_min, n_total, Core, `Cohort-replicated`,
               `Mouse-replicated`, COLOC, `Positive control`, n_pc, n_lnc, n_sgrna, min_mice, chosen)])

# canonical check vs the rebuilt CSV
if (file.exists(LIBRARY_CSV)) {
  lib <- fread(LIBRARY_CSV); chosenN <- opts[chosen == TRUE, n_total]
  cat(sprintf("CANONICAL CHECK: chosen scenario total=%d vs cas13_library.csv rows=%d (%s)\n",
              chosenN, nrow(lib), ifelse(chosenN == nrow(lib), "MATCH", "MISMATCH")))
}

# =============================================================================
# PANEL A: options table
# =============================================================================
opts_tbl <- opts[scenario %in% c("Option A", "Option B")]
# Columns split into two CATEGORIES (divider drawn below):
#   [evidence tiers: dedup-counted, Pos.ctrl first, sum to Total] | [library composition + experiment]
tbl <- data.table(
  Composition = opts_tbl$scenario,
  Total    = comma(opts_tbl$n_total),
  `Pos.\nctrl` = comma(opts_tbl$`Positive control`),
  Core     = comma(opts_tbl$Core),
  `Cohort` = paste0(">=", opts_tbl$cohort_min),
  `PC\nLFC` = paste0(">=", LFC_PC),
  `lnc\nLFC` = "none",
  `PC hep\nCPM` = ">=1",
  `lnc hep\nCPM` = ">=0.1",
  `Cohort\nrep.` = comma(opts_tbl$`Cohort-replicated`),
  `Mouse\nrep.` = comma(opts_tbl$`Mouse-replicated`),
  COLOC    = comma(opts_tbl$COLOC),
  PC = comma(opts_tbl$n_pc), lncRNA = comma(opts_tbl$n_lnc),
  sgRNAs = comma(opts_tbl$n_sgrna), `Mice\n(500x)` = opts_tbl$min_mice)
N_TIER_COLS <- 12L   # Composition,Total,Pos.ctrl,Core,Cohort,PC/lnc LFC,PC/lnc CPM,Cohort rep.,Mouse rep.,COLOC ; PC.. = category 2
# header shading by category: evidence tiers grey; biotype+guides (PC/lncRNA/sgRNAs) blue-grey; Mice (output) green
hdr_fill <- c(rep("#E8E8E8", N_TIER_COLS), rep("#D6E3EA", ncol(tbl) - N_TIER_COLS - 1L), "#DCE8D8")
tt <- ttheme_minimal(base_size = 6, base_family = "Helvetica",
  core = list(bg_params = list(fill = "white",
                               col = "gray80", lwd = 0.4),
              fg_params = list(fontsize = 6, hjust = 0.5, x = 0.5)),
  colhead = list(bg_params = list(fill = hdr_fill, col = "gray60", lwd = 0.5),
                 fg_params = list(fontsize = 6, fontface = "plain", hjust = 0.5, x = 0.5)))
tg <- tableGrob(tbl, rows = NULL, theme = tt)
tg <- gtable_add_grob(tg, rectGrob(gp = gpar(lwd = 1.0, col = "gray40", fill = NA)),
                      t = 1, b = nrow(tg), l = 1, r = ncol(tg))
# thick vertical divider between the evidence-tier block and the composition/experiment block
tg <- gtable_add_grob(tg, segmentsGrob(x0 = 0, y0 = 0, x1 = 0, y1 = 1,
                      gp = gpar(lwd = 2.0, col = "gray30")),
                      t = 1, b = nrow(tg), l = N_TIER_COLS + 1, r = N_TIER_COLS + 1)
panel_A <- wrap_elements(full = tg) + ggtitle("A") +
  theme(plot.title = element_text(face = "plain", size = 6, hjust = 0))

# =============================================================================
# PANEL B: composition stacked bar (canonical library only)
# =============================================================================
chosen_dt <- scenario("Canonical library", CHOSEN_COHORT, c("core","cohort","mouse","coloc","control"), chosen = TRUE)
tier_lvls <- c("Positive control","COLOC","Mouse-replicated","Cohort-replicated","Core")
comp <- melt(chosen_dt[, c("scenario", tier_lvls), with = FALSE], id.vars = "scenario",
             variable.name = "tier", value.name = "n")
comp[, tier := factor(tier, levels = tier_lvls)]
panel_B <- ggplot(comp, aes("", n, fill = tier)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.2) +
  geom_text(data = data.table(n_total = chosen_dt$n_total),
            aes(x = "", y = n_total, label = comma(n_total)),
            vjust = -0.4, size = GEOM_TEXT_6PT, fontface = "plain", inherit.aes = FALSE) +
  scale_fill_manual(values = tier_colors, name = NULL, breaks = rev(tier_lvls)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Library target genes", title = "B") +
  theme_masld() + theme_pub() +
  theme(plot.title = element_text(face = "plain", size = 6, hjust = 0),
        legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# --- Non-deduplicated tier-count table (placed UNDER the Panel B bar) ----------
# The bar above is de-duplicated (priority core>cohort>mouse>COLOC>ctrl; sums to the
# library total). This table adds the NON-deduplicated count: how many library genes
# meet each tier's criterion regardless of which tier they were assigned to. Rows
# overlap (a gene counts in every tier it satisfies) so they sum to > the total --
# which is why an assigned bar segment (e.g. Mouse-replicated) can be far smaller than
# the number of genes that actually pass that tier's criterion.
nd_tbl <- data.table(
  ` `           = c("Assigned (bar)", "Meets criterion"),
  Core          = comma(c(chosen_dt$Core,                chosen_dt$Core_all)),
  `Cohort\nrep.`= comma(c(chosen_dt$`Cohort-replicated`, chosen_dt$Cohort_all)),
  `Mouse\nrep.` = comma(c(chosen_dt$`Mouse-replicated`,  chosen_dt$Mouse_all)),
  COLOC         = comma(c(chosen_dt$COLOC,               chosen_dt$COLOC_all)),
  `Pos.\nctrl`  = comma(c(chosen_dt$`Positive control`,  chosen_dt$Control_all)))
ttB <- ttheme_minimal(base_size = 6, base_family = "Helvetica",
  core = list(bg_params = list(fill = c("white", "#F2F2F2"), col = "gray80", lwd = 0.4),
              fg_params = list(fontsize = 6, hjust = 0.5, x = 0.5)),
  colhead = list(bg_params = list(fill = "#E8E8E8", col = "gray60", lwd = 0.5),
                 fg_params = list(fontsize = 6, fontface = "plain", hjust = 0.5, x = 0.5)))
ndg <- tableGrob(nd_tbl, rows = NULL, theme = ttB)
ndg <- gtable_add_grob(ndg, rectGrob(gp = gpar(lwd = 1.0, col = "gray40", fill = NA)),
                       t = 1, b = nrow(ndg), l = 1, r = ncol(ndg))
panel_B_full <- panel_B / wrap_elements(full = ndg) + plot_layout(heights = c(3.1, 1))
message(sprintf("Fig 11b: bar = de-dup assignment (sums to %s). Table 'Meets criterion' row = non-dedup, overlapping (e.g. Mouse-rep %s genes pass >=3-diet replication vs %s uniquely assigned).",
                comma(chosen_dt$n_total), comma(chosen_dt$Mouse_all), comma(chosen_dt$`Mouse-replicated`)))

# =============================================================================
# PANEL C: coverage-vs-mice curve (canonical library only)
# =============================================================================
COV_CAP <- 1000L; COV_LINES <- c(100L,250L,500L,750L)
chosen_row <- opts[chosen == TRUE]
chosen_sgrna <- chosen_row$n_sgrna
chosen_mice  <- chosen_row$min_mice
mice_range <- 5:(chosen_mice + 10)
cc <- data.table(mice = mice_range,
                 coverage = (mice_range * eff_inject) / chosen_sgrna)
cov_at_min <- (chosen_mice * eff_inject) / chosen_sgrna
panel_C <- ggplot(cc, aes(mice, coverage)) +
  geom_hline(yintercept = 500, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_text(data = data.table(y = COV_LINES, lab = paste0(COV_LINES,"x")), inherit.aes = FALSE,
            aes(x = min(mice_range), y = y, label = lab), size = GEOM_TEXT_6PT, color = "gray45", hjust = 0, vjust = -0.3) +
  geom_line(linewidth = 0.9, color = "#C0143C") +
  geom_point(data = data.table(x = chosen_mice, y = cov_at_min),
             aes(x = x, y = y), size = 3, color = "#C0143C") +
  annotate("text", x = chosen_mice, y = cov_at_min + 60,
           label = paste0(chosen_mice, " mice"), size = GEOM_TEXT_6PT, fontface = "plain") +
  scale_y_continuous(breaks = c(COV_LINES, COV_CAP), labels = comma) +
  coord_cartesian(ylim = c(0, COV_CAP)) +
  labs(x = "Mice to inject (for 500x coverage)", y = "Coverage per sgRNA", title = "C",
       caption = sprintf("%s hep/mouse x %.0f%% transduction x %.0f%% Cre = %s informative; binding pool ~%s cells/mouse; %d gRNA/target.",
                         comma(HEP), TRANSD*100, CRE_EFF*100, comma(HEP*TRANSD*CRE_EFF), comma(round(eff_cells)), N_SGRNA)) +
  theme_masld() + theme_pub() +
  theme(plot.title = element_text(face = "plain", size = 6, hjust = 0),
        plot.caption = element_text(size = 6, color = "gray50", hjust = 0))

# =============================================================================
# Save panels
# =============================================================================
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
panels <- list(`11a_library_options_table` = list(p = panel_A, w = 7.09, h = 1.85),
               `11b_library_composition`   = list(p = panel_B_full, w = 5.0, h = 4.7),
               `11c_coverage_curves`       = list(p = panel_C, w = 7.2, h = 4.4))
for (nm in names(panels)) {
  ggsave(file.path(OUT_DIR, paste0(nm, ".pdf")), panels[[nm]]$p,
         width = panels[[nm]]$w, height = panels[[nm]]$h, device = pdf_device)
  message("Saved: ", nm)
}
cat("\nfigS_cas13_library_options.R complete.\n")
