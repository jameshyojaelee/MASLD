#!/usr/bin/env Rscript

# Compare canonical Option B (regular human LFC/FDR, ashr mouse) with an
# all-ashr sensitivity build. The all-ashr build uses lfsr < 0.05 and posterior
# mean LFC > 0.3 for integrated, per-cohort, and mouse DEG evidence.

suppressPackageStartupMessages({
  library(ashr)
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

# Load the exact Option B tier/gating implementation and its source objects.
source(file.path(BASE, "scripts/figures/figS_cas13_library_options.R"),
       local = environment())

LFC <- 0.3
LFSR <- 0.05

ashr_columns <- function(d) {
  fit_rows <- !is.na(d$logFC) & !is.na(d$SE) & d$SE > 0
  d[, `:=`(ashr_lfc = NA_real_, ashr_lfsr = NA_real_)]
  if (any(fit_rows)) {
    fit <- ash(d$logFC[fit_rows], d$SE[fit_rows], method = "fdr")
    d[fit_rows, `:=`(ashr_lfc = get_pm(fit), ashr_lfsr = get_lfsr(fit))]
  }
  d
}

# Integrated canonical file already carries the production ashr fit.
can_ash <- fread(CANONICAL,
                 select = c("gene", "shrunk_logFC", "lfsr"))
can_ash[, hb := strip_v(gene)]
core_h_ash <- can_ash[!is.na(lfsr) & lfsr < LFSR &
                        !is.na(shrunk_logFC) & shrunk_logFC > LFC, unique(hb)]
core_ash <- h2m(core_h_ash)
human_ashr_lfc <- setNames(can_ash$shrunk_logFC, can_ash$hb)

# Per-study files have limma estimates but no ashr columns, so fit each cohort
# independently from its moderated coefficient and standard error.
cohort_lists_ash <- lapply(COHORTS, function(cohort) {
  d <- fread(file.path(PERSTUDY, paste0(cohort, "_de_results.csv")),
             select = c("gene", "logFC", "SE"))
  d[, hb := strip_v(gene)]
  d <- ashr_columns(d)
  d[!is.na(ashr_lfsr) & ashr_lfsr < LFSR &
      !is.na(ashr_lfc) & ashr_lfc > LFC, unique(hb)]
})
cohort_n_ash <- table(unlist(cohort_lists_ash))
cohort_ash <- h2m(names(cohort_n_ash)[cohort_n_ash >= CHOSEN_COHORT])

# Mouse DEG membership is already ashr-based in the production implementation.
# For the all-ashr comparison, use the integrated human posterior mean for the
# directional-concordance check attached to mouse-confirmed genes.
upd <- lapply(de, function(d) {
  d[, gb := strip_v(gene)]
  d[!is.na(lfsr) & lfsr < LFSR &
      !is.na(shrunk_logFC) & shrunk_logFC > LFC, unique(gb)]
})
mouse_n_ash <- table(unlist(upd))
mouse3_ash <- names(mouse_n_ash)[mouse_n_ash >= MIN_DIETS]
mcm_ash <- copy(om[mb %in% mouse3_ash])
mcm_ash[, human_ashr_lfc := human_ashr_lfc[hb]]
mouse_ash <- intersect(mcm_ash[!is.na(human_ashr_lfc) & human_ashr_lfc > 0, mb], pclnc)

assemble <- function(core_set, cohort_set, mouse_set) {
  ct_g <- posctrl_mouse
  core_g <- setdiff(core_set, ct_g)
  cohort_g <- setdiff(cohort_set, union(ct_g, core_g))
  mouse_g <- setdiff(mouse_set, Reduce(union, list(ct_g, core_g, cohort_g)))
  coloc_g <- setdiff(coloc_mouse,
                     Reduce(union, list(ct_g, core_g, cohort_g, mouse_g)))
  all_g <- setdiff(gate_pc(Reduce(union,
                                  list(ct_g, core_g, cohort_g, mouse_g, coloc_g))),
                   ung_ids)
  tier <- data.table(
    gene_id_mouse = all_g,
    tier = fifelse(all_g %in% ct_g, "positive_control",
           fifelse(all_g %in% core_g, "core",
           fifelse(all_g %in% cohort_g, "cohort_replicated",
           fifelse(all_g %in% mouse_g, "mouse_confirmed", "coloc")))))
  tier
}

regular <- assemble(core_mouse_for(LFC_LNC), cohort_mouse_for(LFC_LNC, CHOSEN_COHORT),
                    mouse_conf_for(LFC_LNC))
all_ashr <- assemble(core_ash, cohort_ash, mouse_ash)

summarize_build <- function(x, label) {
  counts <- x[, .N, by = tier]
  data.table(
    method = label,
    total = nrow(x),
    protein_coding = sum(x$gene_id_mouse %in% pc_ids),
    lncRNA = sum(x$gene_id_mouse %in% lnc_ids),
    positive_control = counts[tier == "positive_control", N],
    core = counts[tier == "core", N],
    cohort_replicated = counts[tier == "cohort_replicated", N],
    mouse_confirmed = counts[tier == "mouse_confirmed", N],
    coloc = counts[tier == "coloc", N]
  )
}

summary <- rbind(
  summarize_build(regular, "regular human + ashr mouse"),
  summarize_build(all_ashr, "ashr human + ashr mouse")
)

delta <- merge(regular[, .(gene_id_mouse, regular_tier = tier)],
               all_ashr[, .(gene_id_mouse, ashr_tier = tier)],
               by = "gene_id_mouse", all = TRUE)
delta[, status := fifelse(is.na(regular_tier), "added_by_all_ashr",
                  fifelse(is.na(ashr_tier), "removed_by_all_ashr",
                  fifelse(regular_tier != ashr_tier, "tier_changed", "unchanged")))]
delta <- merge(delta, mm[, .(gene_id_mouse = gb,
                             gene_symbol_mouse = mouse_symbol_gtf,
                             biotype = bt)],
               by = "gene_id_mouse", all.x = TRUE)

out_summary <- file.path(BASE, "Cas13_Library_Design/data/all_ashr_library_comparison.csv")
out_delta <- file.path(BASE, "Cas13_Library_Design/data/all_ashr_library_gene_delta.csv")
fwrite(summary, out_summary)
fwrite(delta, out_delta)

print(summary)
cat("\nGene-level changes:\n")
print(delta[, .N, by = status][order(status)])
cat("\nWrote:\n", out_summary, "\n", out_delta, "\n")
