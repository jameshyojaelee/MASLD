#!/usr/bin/env Rscript
# sex_v3/13b_bootstrap_aggregate.R
# ---------------------------------------------------------------------------
# Pillar 4 aggregator — combine per-rep dream-bootstrap CSVs from
# intermediates/boot_v6/rep_*.csv into per-gene stability fractions across reps.
#
# Output: bootstrap_stability_v6.csv
# Validation: τ=0.8 retention of v5 Strong calls ≥ 80%.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SEXV3 <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "results/integration/sex_v3")
IDIR  <- file.path(SEXV3, "intermediates")

# Contrast routing (sex_v3_utils.R::contrast_paths) — overrides SEXV3 / IDIR
# when CONTRAST_NAME != "disease_vs_ctrl"; default preserves legacy layout.
if (!exists("contrast_paths")) source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT",
             "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3/sex_v3_utils.R"))
.cpaths <- contrast_paths()
SEXV3 <- .cpaths$sexv3
IDIR  <- .cpaths$idir
dir.create(IDIR, recursive = TRUE, showWarnings = FALSE)
BOOT_DIR <- file.path(IDIR, "boot_v6")
OUT_CSV  <- file.path(SEXV3, "bootstrap_stability_v6.csv")
V5_CSV   <- file.path(SEXV3, "interaction_classifier_v5.csv")

cat("=== sex_v6 P4 bootstrap aggregator ===\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

rep_files <- list.files(BOOT_DIR, pattern = "^rep_\\d+\\.csv$", full.names = TRUE)
n_files <- length(rep_files)
cat("rep_*.csv found:", n_files, "\n")
stopifnot(n_files >= 1)

rep_list <- list(); n_failed <- 0L
for (f in rep_files) {
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || all(is.na(d$gene))) { n_failed <- n_failed + 1L; next }
  if (!"class_v6_rep" %in% names(d)) { n_failed <- n_failed + 1L; next }
  rep_id <- as.integer(sub(".*rep_(\\d+)\\.csv$", "\\1", f))
  d[, rep_id := rep_id]
  rep_list[[length(rep_list) + 1L]] <- d[, .(gene, rep_id, class_v6_rep)]
}
cat("Valid reps:", length(rep_list), "   failed/NA reps:", n_failed, "\n")
stopifnot(length(rep_list) >= 1)

all_reps <- rbindlist(rep_list)
cat("Total gene × rep rows:", nrow(all_reps), "\n")

classes <- c("Female_biased", "Male_biased",
             "Female_biased_M_underpowered", "Male_biased_F_underpowered",
             "Divergent", "Divergent_one_sided",
             "Sex_modifier", "Concordant", "Not_DEG")
all_reps[, class_v6_rep := factor(class_v6_rep, levels = classes)]
wide <- dcast(all_reps[!is.na(class_v6_rep)],
              gene ~ class_v6_rep,
              fun.aggregate = length, value.var = "rep_id", drop = FALSE)
n_boot_used <- all_reps[!is.na(class_v6_rep),
                        .(n_boot_used = uniqueN(rep_id)), by = gene]
wide <- merge(wide, n_boot_used, by = "gene", all.x = TRUE)

for (cl in classes) {
  if (!cl %in% names(wide)) wide[, (cl) := 0L]
  wide[, paste0("stability_", cl) := get(cl) / pmax(n_boot_used, 1L)]
}
wide[, stability_F      := stability_Female_biased + stability_Female_biased_M_underpowered]
wide[, stability_M      := stability_Male_biased + stability_Male_biased_F_underpowered]
wide[, stability_div    := stability_Divergent + stability_Divergent_one_sided]
wide[, stability_conc   := stability_Concordant]
wide[, stability_notDEG := stability_Not_DEG]
wide[, stability_modif  := stability_Sex_modifier]

frac_mat <- as.matrix(wide[, .(stability_F, stability_M, stability_div,
                               stability_conc, stability_notDEG, stability_modif)])
mode_labels <- c("Female_biased", "Male_biased", "Divergent", "Concordant",
                 "Not_DEG", "Sex_modifier")
modal_idx <- max.col(frac_mat, ties.method = "first")
wide[, modal_class := mode_labels[modal_idx]]
wide[, stability_modal := frac_mat[cbind(seq_len(.N), modal_idx)]]
wide[, stability_pass_08 := stability_modal >= 0.8]
wide[, stability_pass_06 := stability_modal >= 0.6]

out_cols <- c("gene", "n_boot_used",
              paste0("stability_", classes),
              "stability_F", "stability_M", "stability_div",
              "stability_conc", "stability_notDEG", "stability_modif",
              "modal_class", "stability_modal",
              "stability_pass_08", "stability_pass_06")
out <- wide[, ..out_cols]
setorder(out, -stability_modal)
tmp <- paste0(OUT_CSV, ".tmp"); fwrite(out, tmp); file.rename(tmp, OUT_CSV)
cat("Wrote:", OUT_CSV, "  rows:", nrow(out), "\n")

cat("\nModal class distribution:\n"); print(table(out$modal_class))
cat("\nPass at τ=0.8 by modal class:\n")
print(out[, .(N = .N, pass = sum(stability_pass_08)), by = modal_class])
cat("\nPass at τ=0.6 by modal class:\n")
print(out[, .(N = .N, pass = sum(stability_pass_06)), by = modal_class])

# ---- Acceptance: τ=0.8 retention of v5 Strong calls ≥ 80% ----
# "v5 Strong" = interaction_classifier_v5.csv class ∈ {Female_biased, Male_biased,
# Divergent, Sex_modifier} with evidence_tier_A_B_C == "A"
if (file.exists(V5_CSV)) {
  v5 <- fread(V5_CSV, select = c("gene", "class_v5_interaction", "evidence_tier_A_B_C"))
  strong_classes <- c("Female_biased", "Male_biased", "Divergent",
                      "Divergent_one_sided", "Sex_modifier",
                      "Female_biased_M_underpowered", "Male_biased_F_underpowered")
  v5_strong <- v5[evidence_tier_A_B_C == "A" & class_v5_interaction %in% strong_classes, gene]
  cat("\nv5 Strong (Tier-A interaction calls):", length(v5_strong), "\n")
  if (length(v5_strong) > 0) {
    boot_pass <- out[stability_pass_08 == TRUE & modal_class %in%
                       c("Female_biased", "Male_biased", "Divergent", "Sex_modifier"), gene]
    retain <- length(intersect(v5_strong, boot_pass)) / length(v5_strong)
    cat("τ=0.8 retention of v5 Strong:", round(retain * 100, 1),
        "% (", length(intersect(v5_strong, boot_pass)), "/", length(v5_strong), ")\n")
    cat("Acceptance gate (≥ 80%):", if (retain >= 0.8) "PASS" else "FAIL", "\n")
  }
}
cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
