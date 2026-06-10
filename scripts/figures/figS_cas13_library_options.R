#!/usr/bin/env Rscript
# =============================================================================
# Figure S_lib_6 -- Library composition and experimental feasibility
# KEY MESSAGE: Library size trades off between human-evidence breadth (looser
# ashr threshold) and experimental feasibility (fewer mice). The human arm
# sweeps an ashr shrunk_logFC ladder (the SAME instrument used on the mouse
# arm). HUMAN-ANCHORED: human ashr spine UNION mouse cross-diet>=3 that is
# ALSO human-direction-concordant. Selected design = Human ashr>0.2 (PI 2026-06-01).
#
# v6 (2026-06-04): the human arm is now sourced from the limma-voom + metafor
# (REML) ashr DEGs (meta_results_ashr.csv), replacing dream. This is a human-arm
# THRESHOLD-EXPLORATION figure (PC+lncRNA, TPM-gated); the v6 library ALSO carries
# additive MASH-progression and miRNA-conserved tiers, whose full 4-tier x 3-biotype
# composition is shown in cas13_v6_tier_composition.pdf (not re-derived here).
# =============================================================================
# Panels:
#   A  Options comparison table (per-option 4-tier + biotype counts, via gridExtra)
#   B  Gene composition stacked bars (4 v6 tiers: human/MASH/mouse-conf/miRNA)
#   C  Coverage-vs-mice curves per option (efficiency- + mortality-aware)
#
# Output: Cas13_Library_Design/figures/11{a,b,c}_*.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(gridExtra)
  library(scales)
  library(gtable)
})

# -- Project paths & theme ----------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS_CAS13LIB_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -- Constants ----------------------------------------------------------------
DIETS_ALL <- c("MCD", "CDAHFD", "Western", "HFD")  # NASH+Western merged 2026-05-29 (4 groups)
PERDIET_DIR <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
ORTHO_TABLE <- file.path(BASE, "data/external/orthologs/master_ortholog_table.tsv.gz")
MOUSE_META  <- file.path(BASE, "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv")
ASHR_PATH   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                         "results/integration/meta_results_ashr.csv")  # v6: limma-voom+metafor ashr (was dream)
MASH_PATH   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                         "results/disease_signatures/nafl_vs_nash_meta_ashr.csv")  # v6: MASH-vs-MASL metafor ashr
LIBRARY_CSV <- file.path(BASE, "Cas13_Library_Design/data/cas13_library_v3.0.csv") # canonical v6 (for verification)

# Mouse UP-DEG definition: ashr-shrunk effect size (M02c), lfsr<0.05 & shrunk>0.5.
LFSR_THR       <- 0.05
SHRUNK_LFC_THR <- 0.5

# v6: the options ladder sweeps the human shrunk_logFC threshold (human spine).
# Each option = human spine(threshold) UNION the three FIXED v6 tiers:
#   mouse-confirmed (>=3 of 4 diets, human-concordant) + MASH-progression (MASH-vs-MASL
#   UP) + miRNA-conserved (tier-H ortholog). Tier priority on overlap:
#   human > mouse_confirmed > mash_progression > mirna_conserved. NO TPM gate, NO GWAS
#   tier (the v5 deployability/GWAS exploration is retired) -> option F == canonical v6.
HUMAN_LFSR_THR <- 0.05
HUMAN_THRS     <- c(1.00, 0.75, 0.50, 0.40, 0.30, 0.20, 0.15, 0.10)  # A..H (strict -> loose)
CHOSEN_THR     <- 0.20                                    # selected design (option F = canonical v6)
MASH_SHRUNK    <- 0.20                                    # MASH-vs-MASL UP threshold (fixed)
KEEP_BIOTYPES  <- c("protein_coding", "lncRNA", "miRNA")  # v6 adds miRNA
SUF     <- ""

# Experimental parameters
N_SGRNA_PC   <- 4L    # gRNA per protein-coding gene
N_SGRNA_LNC  <- 4L    # gRNA per lncRNA gene
N_SGRNA_MIR  <- 4L    # gRNA per miRNA gene (pri-miRNA / host transcript)
N_CTRL_GENES <- 100L  # control protein-coding genes
N_SGRNA_CTRL <- 4L    # gRNA per control gene
N_NT_GUIDES  <- 500L  # non-targeting (NT) guides

# Delivery / sort / survival efficiencies (PI 2026-06-01)
HEP_PER_MOUSE    <- 10000000L  # hepatocytes harvested per mouse
TRANSDUCTION_EFF <- 0.30       # lenti transduction (fraction carrying a guide)
CRE_EFF          <- 0.80       # Cre recombination -> only Cre+ cells have active Cas13 (informative)
FACS_EFF         <- 0.60       # FACS sorting efficiency (sorted arms only)
GATE_FRACTION    <- 0.15       # per FACS gate: top 15% AND bottom 15% (one sort, two bins)
COVERAGE_TARGET  <- 500L       # cells per gRNA per sequenced population
INPUT_FRACTION   <- 0.13       # transduced cells set aside PRE-SORT as unsorted input (no FACS loss)
MORTALITY        <- 1/8        # 1 in 8 mice die after viral injection

# Per-mouse INFORMATIVE (Cre+) cell yields feeding 3 sequenced pools: unsorted
# input + top-15% arm + bottom-15% arm. Sorted arms also lose FACS_EFF; mice are
# set by the smallest-yield pool (binding constraint), then scaled for mortality.
TRANSDUCED_PER_MOUSE  <- HEP_PER_MOUSE * TRANSDUCTION_EFF                 # 3.0e6 guide+
INFORMATIVE_PER_MOUSE <- TRANSDUCED_PER_MOUSE * CRE_EFF                   # 2.4e6 Cre+ (perturbed)
INPUT_CELLS_PER_MOUSE <- INFORMATIVE_PER_MOUSE * INPUT_FRACTION
ARM_CELLS_PER_MOUSE   <- INFORMATIVE_PER_MOUSE * (1 - INPUT_FRACTION) * GATE_FRACTION * FACS_EFF
EFF_CELLS_PER_MOUSE   <- min(ARM_CELLS_PER_MOUSE, INPUT_CELLS_PER_MOUSE) # ~187,920 (arm binds)

# Option color palette (colorblind-safe)
option_colors <- c(A = "#E64B35", B = "#4DBBD5", C = "#00A087", D = "#F39B7F",
                   E = "#3C5488", F = "#8491B4", G = "#7E6148", H = "#DC0000")
# v6 tier colors (match figS_cas13_library_core_vs_mash_overlap.R)
tier_colors <- c(`Human` = "#e35070", `MASH-progression` = "#b771e6",
                 `Mouse-confirmed` = "#4baeef", `miRNA-conserved` = "#30d796")

# =============================================================================
# 1. Load per-diet mouse DE, compute per-gene diet replication (mouse arm)
# =============================================================================
message("Loading per-diet mouse DE ...")
de_list <- lapply(DIETS_ALL, function(d) {
  f <- file.path(PERDIET_DIR, paste0(d, "_de_results.csv"))
  if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
  dt <- fread(f); dt[, gene_base := sub("\\..*", "", gene)]; dt
})
names(de_list) <- DIETS_ALL
up_genes_per_diet <- lapply(DIETS_ALL, function(d) {
  dt <- de_list[[d]]; if (is.null(dt)) return(character(0))
  dt[lfsr < LFSR_THR & shrunk_logFC > SHRUNK_LFC_THR, gene_base]
})
names(up_genes_per_diet) <- DIETS_ALL
all_up <- unique(unlist(up_genes_per_diet))
gene_diet_mat <- data.table(gene_base = all_up)
for (d in DIETS_ALL) gene_diet_mat[, (d) := gene_base %in% up_genes_per_diet[[d]]]
gene_diet_mat[, n_diets := rowSums(.SD), .SDcols = DIETS_ALL]
cat(sprintf("Mouse UP genes: total union = %d | >=3 diets: %d\n",
            nrow(gene_diet_mat), sum(gene_diet_mat$n_diets >= 3)))

# =============================================================================
# 2. Mouse gene metadata (biotype, normalized like rebuild_cas13_library.R)
# =============================================================================
mouse_meta <- fread(MOUSE_META)
mouse_meta[, gene_base := mouse_ensembl_base]
mouse_meta[grepl("protein_coding", mouse_biotype), bt2 := "protein_coding"]
mouse_meta[grepl("lncRNA|lincRNA", mouse_biotype), bt2 := "lncRNA"]
mouse_meta[grepl("miRNA", mouse_biotype), bt2 := "miRNA"]
pc_lnc_ids <- mouse_meta[bt2 %in% c("protein_coding", "lncRNA"), gene_base]  # human/mouse/MASH tiers
pc_ids     <- mouse_meta[bt2 == "protein_coding", gene_base]
lnc_ids    <- mouse_meta[bt2 == "lncRNA", gene_base]
mir_ids    <- mouse_meta[bt2 == "miRNA", gene_base]                          # miRNA tier

# v6: NO TPM scoreability gate (retired with the v5 deployability exploration). The
# canonical v6 library carries hepatocyte substrate as a SOFT annotation, not a hard
# filter -- so the option sizes here match rebuild_cas13_library.R exactly.

# =============================================================================
# 3. Ortholog bridge (deterministic, one2one-preferred) + human ashr DE
# =============================================================================
message("Loading ortholog bridge + human ashr DE ...")
ortho_raw <- fread(ORTHO_TABLE,
                   select = c("mouse_ensembl", "human_ensembl", "mouse_symbol",
                              "human_symbol", "confidence_tier", "is_one2one"))
ortho_raw[, mouse_ensembl := sub("\\..*", "", mouse_ensembl)]
ortho_raw[, human_ensembl := sub("\\..*", "", human_ensembl)]
ortho <- ortho_raw[confidence_tier %in% c("H", "M")]
ortho <- unique(ortho, by = c("mouse_ensembl", "human_ensembl"))
ortho[, trank := match(confidence_tier, c("H", "M"))]
ortho[, one2one_rank := ifelse(is_one2one %in% c(TRUE, "True", "TRUE", "true"), 0L, 1L)]
# best mouse per human (spine)
ortho_byhuman <- copy(ortho); setorder(ortho_byhuman, human_ensembl, trank, one2one_rank, mouse_ensembl)
ortho_byhuman <- unique(ortho_byhuman, by = "human_ensembl")
# best human per mouse (for mouse-tier concordance + control mapping)
ortho_m2h <- copy(ortho); setorder(ortho_m2h, mouse_ensembl, trank, one2one_rank, human_ensembl)
ortho_m2h <- unique(ortho_m2h, by = "mouse_ensembl")

ash <- fread(ASHR_PATH, select = c("gene", "logFC", "shrunk_logFC", "lfsr"))
ash[, hb := sub("\\..*", "", gene)]
human_logfc <- ash[, .(hb, hlfc = logFC)]
get_human_spine <- function(shrunk_thr) {
  hup <- unique(ash[!is.na(lfsr) & lfsr < HUMAN_LFSR_THR & shrunk_logFC > shrunk_thr, hb])
  unique(ortho_byhuman[human_ensembl %in% hup, mouse_ensembl])
}
mouse_3diets <- gene_diet_mat[n_diets >= 3, gene_base]
# mouse-confirmed tier: >=3 diets AND best human ortholog is up in human dream
mc_map <- merge(ortho_m2h[mouse_ensembl %in% mouse_3diets, .(mouse_ensembl, hb = human_ensembl)],
                human_logfc, by = "hb", all.x = TRUE)
mouse_conf <- mc_map[!is.na(hlfc) & hlfc > 0, mouse_ensembl]
# MASH-progression tier: mouse orthologs of MASH-vs-MASL UP human DEGs (lfsr<0.05 & shrunk>thr)
mash <- fread(MASH_PATH, select = c("gene", "shrunk_logFC", "lfsr"))
mash[, hb := sub("\\..*", "", gene)]
mash_up  <- unique(mash[!is.na(lfsr) & lfsr < HUMAN_LFSR_THR & shrunk_logFC > MASH_SHRUNK, hb])
mash_set <- unique(ortho_byhuman[human_ensembl %in% mash_up, mouse_ensembl])
# miRNA-conserved tier: super-confident ortholog-conserved mouse miRNAs (DE-independent):
# >=2 of {miRBase family-ID, MirGeneDB, biomaRt/Compara} (matches rebuild_cas13_library.R).
mir_raw <- fread(cmd = paste0("zcat ", ORTHO_TABLE),
                 select = c("mouse_ensembl", "mouse_biotype", "tier_H_mirbase", "tier_H_mirgenedb", "tier_H_biomart"))
mir_raw[, mouse_ensembl := sub("\\..*", "", mouse_ensembl)]
tier1 <- function(x) as.integer(x %in% c(1, "1", TRUE))
mir_raw[, n_ev := tier1(tier_H_mirbase) + tier1(tier_H_mirgenedb) + tier1(tier_H_biomart)]
mirna_set <- unique(mir_raw[mouse_biotype == "miRNA" & n_ev >= 2L, mouse_ensembl])
cat(sprintf("Human spine (mouse orthologs): >0.5=%d >0.3=%d >0.2=%d >0.15=%d | mouse-confirmed=%d | MASH=%d | miRNA=%d\n",
            length(get_human_spine(0.5)), length(get_human_spine(0.3)),
            length(get_human_spine(0.2)), length(get_human_spine(0.15)),
            length(mouse_conf), length(mash_set), length(mirna_set)))

# =============================================================================
# 4. Build option gene sets (human-anchored: spine UNION mouse-confirmed)
# =============================================================================
build_option <- function(human_thr, label) {
  h0  <- intersect(get_human_spine(human_thr), pc_lnc_ids)   # human spine (PC/lncRNA)
  mc0 <- intersect(mouse_conf, pc_lnc_ids)                   # mouse-confirmed (PC/lncRNA)
  ma0 <- intersect(mash_set,  pc_lnc_ids)                    # MASH-progression (PC/lncRNA)
  mi0 <- intersect(mirna_set, mir_ids)                       # miRNA-conserved
  # tier priority on overlap (matches rebuild_cas13_library.R):
  #   human > mouse_confirmed > mash_progression > mirna_conserved
  human   <- h0
  mouse_c <- setdiff(mc0, human)
  mash_p  <- setdiff(ma0, union(human, mouse_c))
  mirna_c <- setdiff(mi0, Reduce(union, list(human, mouse_c, mash_p)))
  gene_set <- Reduce(union, list(human, mouse_c, mash_p, mirna_c))
  n_pc  <- length(intersect(gene_set, pc_ids))
  n_lnc <- length(intersect(gene_set, lnc_ids))
  n_mir <- length(intersect(gene_set, mir_ids))
  n_sgrna <- n_pc * N_SGRNA_PC + n_lnc * N_SGRNA_LNC + n_mir * N_SGRNA_MIR +
             N_CTRL_GENES * N_SGRNA_CTRL + N_NT_GUIDES
  surv_mice <- ceiling(n_sgrna * COVERAGE_TARGET / EFF_CELLS_PER_MOUSE)  # survivors needed @500x
  inj_mice  <- ceiling(surv_mice / (1 - MORTALITY))                      # mice to INJECT
  data.table(
    option = label, human_thr = human_thr,
    n_total = length(gene_set), n_pc = n_pc, n_lnc = n_lnc, n_mir = n_mir,
    n_human = length(human), n_mouse = length(mouse_c),
    n_mash = length(mash_p), n_mirna = length(mirna_c),
    pct_human = round(100 * length(human) / max(length(gene_set), 1), 1),
    n_sgrna = n_sgrna, surv_mice = surv_mice, min_mice = inj_mice,
    chosen = isTRUE(all.equal(human_thr, CHOSEN_THR)),
    gene_set = list(gene_set)
  )
}
options_dt <- rbindlist(lapply(seq_along(HUMAN_THRS), function(i)
  build_option(HUMAN_THRS[i], LETTERS[i])), fill = TRUE)

cat("\n=== v6 Library Options (human spine ∪ mouse-confirmed ∪ MASH ∪ miRNA; priority-assigned) ===\n")
print(options_dt[, .(option, human_thr, n_total, n_human, n_mash, n_mouse, n_mirna,
                     n_pc, n_lnc, n_mir, pct_human, n_sgrna, min_mice, chosen)])

# verify option F reproduces the canonical v6 library tier counts
if (file.exists(LIBRARY_CSV)) {
  libv6 <- fread(LIBRARY_CSV); canon <- libv6[, .N, by = tier]; optF <- options_dt[chosen == TRUE]
  cv <- function(t) { v <- canon[tier == t, N]; if (length(v)) v else 0L }
  cat(sprintf("CANONICAL CHECK (option F vs cas13_library_v3.0.csv):\n  fig: human=%d mouse=%d mash=%d mirna=%d total=%d\n  csv: human=%d mouse_confirmed=%d mash_progression=%d mirna_conserved=%d total=%d\n",
      optF$n_human, optF$n_mouse, optF$n_mash, optF$n_mirna, optF$n_total,
      cv("human"), cv("mouse_confirmed"), cv("mash_progression"), cv("mirna_conserved"), nrow(libv6)))
}

# =============================================================================
# PANEL A: Options comparison table
# =============================================================================
message("Building Panel A (table) ...")
table_display <- data.table(
  Option                = options_dt$option,
  Definition            = sprintf("Human metafor-ashr>%.2f\n∪ mouse ∪ MASH ∪ miRNA", options_dt$human_thr),
  `Total`               = comma(options_dt$n_total),
  `Human`               = comma(options_dt$n_human),
  `MASH\nprog.`         = comma(options_dt$n_mash),
  `Mouse\nconf.`        = comma(options_dt$n_mouse),
  `miRNA`               = comma(options_dt$n_mirna),
  PC                    = comma(options_dt$n_pc),
  lncRNA                = comma(options_dt$n_lnc),
  `% Human`             = sprintf("%.1f%%", options_dt$pct_human),
  sgRNAs                = comma(options_dt$n_sgrna),
  `Mice inject\n(500x)` = options_dt$min_mice
)
tt_theme <- ttheme_minimal(
  base_size = 7, base_family = "Helvetica",
  core    = list(bg_params = list(fill = c("gray97", "white"), col = "gray80", lwd = 0.4),
                 fg_params = list(fontsize = 7, hjust = 0.5, x = 0.5)),
  colhead = list(bg_params = list(fill = "#E8E8E8", col = "gray60", lwd = 0.5),
                 fg_params = list(fontsize = 7, fontface = "bold", hjust = 0.5, x = 0.5)),
  rowhead = list(fg_params = list(fontsize = 7, fontface = "bold")))
tbl_grob <- tableGrob(table_display, rows = NULL, theme = tt_theme)
# (no chosen-row highlight per PI 2026-06-01)
tbl_grob <- gtable_add_grob(tbl_grob,
  grobs = rectGrob(gp = gpar(lwd = 1.0, col = "gray40", fill = NA)),
  t = 1, b = nrow(tbl_grob), l = 1, r = ncol(tbl_grob))
panel_A <- wrap_elements(full = tbl_grob) +
  ggtitle("A") +
  labs(caption = paste0("Each option = Human metafor-ashr>(threshold) spine ∪ mouse-confirmed (≥3 of 4 diets, ",
                        "human-concordant) ∪ MASH-progression (MASH-vs-MASL UP) ∪ miRNA-conserved (tier-H ortholog); ",
                        "tiers priority-assigned (human > mouse > MASH > miRNA). Option F (>0.20) = canonical v6 library. ",
                        "% human reported. Mice = injected (incl. 12.5% post-injection mortality).")) +
  theme(plot.title   = element_text(face = "bold", size = 9, hjust = 0),
        plot.caption = element_text(size = 5, color = "gray40", hjust = 0))

# =============================================================================
# PANEL B: Gene composition stacked bars
# =============================================================================
message("Building Panel B (composition bars) ...")
lfc_labels <- setNames(sprintf(">%.2f", options_dt$human_thr), options_dt$option)
comp_dt <- rbindlist(lapply(1:nrow(options_dt), function(i) {
  opt <- options_dt[i]
  rbindlist(list(
    data.table(option = opt$option, component = "Human",            n_genes = opt$n_human),
    data.table(option = opt$option, component = "MASH-progression", n_genes = opt$n_mash),
    data.table(option = opt$option, component = "Mouse-confirmed",  n_genes = opt$n_mouse),
    data.table(option = opt$option, component = "miRNA-conserved",  n_genes = opt$n_mirna)))
}))
comp_dt[, option := factor(option, levels = options_dt$option)]
comp_dt[, component := factor(component, levels = c("miRNA-conserved", "MASH-progression", "Mouse-confirmed", "Human"))]
totals <- options_dt[, .(option, n_total)]
totals[, option := factor(option, levels = options_dt$option)]
panel_B <- ggplot(comp_dt, aes(x = option, y = n_genes, fill = component)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.2) +
  geom_text(data = totals, aes(x = option, y = n_total, label = comma(n_total), fill = NULL),
            vjust = -0.4, size = 2.3, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = tier_colors, name = NULL) +
  scale_x_discrete(labels = lfc_labels) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Library option (human metafor-ashr threshold)", y = "Library target genes", title = "B") +
  theme_masld() + theme_pub() +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0),
        legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6), axis.text.x = element_text(size = 6))

# =============================================================================
# PANEL C: Coverage-vs-mice curves (x = surviving/contributing mice)
# =============================================================================
message("Building Panel C (coverage curves) ...")
effective_cells <- EFF_CELLS_PER_MOUSE
eff_inject <- effective_cells * (1 - MORTALITY)   # avg usable cells per INJECTED mouse (12.5% die)
COV_CAP   <- 1000L
COV_LINES <- c(100L, 250L, 500L, 750L)
options_dt[, mice_at_cap := ceiling(COV_CAP * n_sgrna / eff_inject)]
max_mice <- max(options_dt$mice_at_cap) + 3
mice_range <- 5:max_mice
coverage_curves <- rbindlist(lapply(1:nrow(options_dt), function(i) {
  opt <- options_dt[i]
  data.table(option = opt$option, mice = mice_range,
             coverage = (mice_range * eff_inject) / opt$n_sgrna)   # x = mice to inject
}))
# annotate at mice-to-inject needed for 500x; label = inject number only
min_mice_pts <- options_dt[, .(option, min_mice, n_sgrna)]
min_mice_pts[, coverage_at_min := (min_mice * eff_inject) / n_sgrna]
setorder(min_mice_pts, min_mice, option)
min_mice_pts[, label_y := COVERAGE_TARGET + 45 + (seq_len(.N) - 1) * 56]  # fits 8 options under COV_CAP
cov_ref <- data.table(y = COV_LINES, lab = paste0(COV_LINES, "x"))
panel_C <- ggplot(coverage_curves, aes(x = mice, y = coverage, color = option)) +
  geom_hline(data = cov_ref[y != COVERAGE_TARGET], inherit.aes = FALSE,
             aes(yintercept = y), linetype = "dotted", color = "gray80", linewidth = 0.3) +
  geom_hline(yintercept = COVERAGE_TARGET, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_text(data = cov_ref, inherit.aes = FALSE, aes(x = 3, y = y, label = lab),
            size = 2.2, color = "gray40", hjust = 0, vjust = -0.3) +
  geom_line(linewidth = 0.7) +
  geom_segment(data = min_mice_pts,
               aes(x = min_mice, xend = min_mice, y = coverage_at_min, yend = label_y, color = option),
               linewidth = 0.25, linetype = "dotted", show.legend = FALSE) +
  geom_point(data = min_mice_pts, aes(x = min_mice, y = coverage_at_min, color = option),
             size = 2.2, shape = 16, show.legend = FALSE) +
  geom_text(data = min_mice_pts,
            aes(x = min_mice, y = label_y, color = option, label = paste0(option, ": ", min_mice, " mice")),
            size = 2.8, fontface = "bold", vjust = -0.3, hjust = 0.5, show.legend = FALSE) +
  scale_color_manual(values = option_colors, name = "Option") +
  scale_x_continuous(breaks = seq(10, ceiling(max_mice / 10) * 10, by = 10),
                     limits = c(3, max_mice + 2)) +
  scale_y_continuous(breaks = c(COV_LINES, COV_CAP), labels = comma) +
  coord_cartesian(ylim = c(0, COV_CAP)) +
  labs(x = "Mice to inject (for 500x coverage)", y = "Coverage per sgRNA", title = "C",
       caption = paste0(
         comma(HEP_PER_MOUSE), " hep/mouse x ", TRANSDUCTION_EFF * 100, "% transduction x ",
         CRE_EFF * 100, "% Cre = ", comma(INFORMATIVE_PER_MOUSE), " informative cells/mouse.\n",
         INPUT_FRACTION * 100, "% unsorted input + top/bottom ", GATE_FRACTION * 100,
         "% gates at ", FACS_EFF * 100, "% FACS -> binding pool ~", comma(round(EFF_CELLS_PER_MOUSE)),
         " cells/mouse.\n",
         "x-axis = mice to inject (12.5% die post-injection).\n",
         N_SGRNA_PC, " gRNA/PC, ", N_SGRNA_LNC, " gRNA/lncRNA, ", N_SGRNA_MIR, " gRNA/miRNA, ",
         N_CTRL_GENES, " ctrl x ", N_SGRNA_CTRL, " + ", N_NT_GUIDES, " NT; 500x per population.")) +
  theme_masld() + theme_pub() +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0),
        plot.caption = element_text(size = 4.5, color = "gray50", hjust = 0),
        legend.position = c(0.87, 0.40), legend.key.width = unit(0.5, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA))

# (Panel D positive-control dot matrix removed from the library design per PI 2026-06-01)

# =============================================================================
# Save panels individually
# =============================================================================
message("Saving panels individually ...")
pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
panels <- list()
panels[[paste0("11a_library_options_table", SUF)]] <- list(plot = panel_A, w = 12.5, h = 4.8)
panels[[paste0("11b_library_composition",   SUF)]] <- list(plot = panel_B, w = fig_half_width, h = 3.5)
panels[[paste0("11c_coverage_curves",       SUF)]] <- list(plot = panel_C, w = 9.0, h = 5.5)
for (nm in names(panels)) {
  p <- panels[[nm]]
  ggsave(file.path(OUT_DIR, paste0(nm, ".pdf")), p$plot, width = p$w, height = p$h, device = pdf_device)
  message("Saved: ", nm)
}
cat("\n=== Final Option Statistics (v6 4-tier) ===\n")
for (i in 1:nrow(options_dt)) {
  opt <- options_dt[i]
  cat(sprintf("  Option %s%s: %s targets (%s PC, %s lnc, %s miRNA) | tiers H=%s/MASH=%s/mouse=%s/miR=%s | %s%% human | %s sgRNAs | %d surv / %d inject\n",
              opt$option, ifelse(opt$chosen, " *", ""), comma(opt$n_total),
              comma(opt$n_pc), comma(opt$n_lnc), comma(opt$n_mir),
              comma(opt$n_human), comma(opt$n_mash), comma(opt$n_mouse), comma(opt$n_mirna),
              opt$pct_human, comma(opt$n_sgrna), opt$surv_mice, opt$min_mice))
}
