#!/usr/bin/env Rscript
# 93_kprogram_sanity.R — Phase 4 biological sanity checks after k-program refactor
#
# Checks:
#   1. Known MASLD genes localize to expected programs:
#      - PNPLA3, TM6SF2, HNF4A, HMGCS2 → Metabolic program top-100 loadings
#      - COL1A1, COL1A2, TGFB1, ACTA2 → Fibrotic program top-100
#      - IL6, TNF, CCL2 → Inflammatory program top-100 (if labeled)
#   2. Fibrotic program F-stage boundary test: dominant-sample frequency rises F0→F4
#      NOTE: "F2 switch" framing retired 2026-05-09; see paper_outline.md
#   3. Female-biased DEG enrichment for the Fibrotic program
#   4. No hardcoded S1/S2 remnants in refactored scripts

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

cat("=== k-Program Sanity Checks ===\n")
cat("Date:", format(Sys.time()), "\n\n")

labels <- fread(file.path(BASE, "RNA-seq/results/subtypes/program_labels.csv"))
markers <- fread(file.path(BASE, "RNA-seq/results/subtypes/subtype_markers.csv"))
nmf <- fread(file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv"))
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id","human_symbol"))
atlas[, ens_base := sub("\\..*", "", ensembl_id)]
markers[, ens_base := sub("\\..*", "", gene)]
markers <- merge(markers, atlas[, .(ens_base, human_symbol)], by = "ens_base", all.x = TRUE)

cat("k =", nrow(labels), "programs\n")
print(labels[, .(program_code, biological_label, label_category)], row.names = FALSE)

# ============================================================
# 1. Known gene localization
# ============================================================
cat("\n--- Check 1: Known MASLD genes in expected programs ---\n")
gene_expect <- list(
  Metabolic    = c("PNPLA3","TM6SF2","HNF4A","HMGCS2","CPT1A","CYP3A4","ALB","ACSL4"),
  Fibrotic     = c("COL1A1","COL1A2","COL3A1","TGFB1","ACTA2","TIMP1","LUM","SPARC"),
  Inflammatory = c("IL6","TNF","CCL2","NFKB1","STAT1","CD14"),
  "Sex-F"      = c("XIST","TSIX"),
  "Sex-M"      = c("KDM5D","USP9Y","UTY","RPS4Y1","DDX3Y"),
  "HCC-like"   = c("AKR1B10","GPC3","AFP","HKDC1")
)
check1 <- list()
for (lbl in names(gene_expect)) {
  # Get the program(s) matching this label_category
  matching_codes <- labels$program_code[labels$label_category == lbl]
  if (length(matching_codes) == 0) {
    check1[[lbl]] <- data.table(label_category = lbl, status = "NO_PROGRAM", hits = "-")
    next
  }
  # For each matching program, check top-100 upregulated markers (padj<0.05, logFC>0)
  hits_str <- vapply(matching_codes, function(pc) {
    top100 <- markers[program_code == pc & direction == "up" & padj < 0.05, human_symbol]
    top100 <- head(top100[order(-markers[program_code == pc & direction == "up" & padj < 0.05, t])], 100)
    found <- intersect(gene_expect[[lbl]], top100)
    paste0(pc, ": ", length(found), "/", length(gene_expect[[lbl]]),
           " (", paste(found, collapse=","), ")")
  }, character(1))
  check1[[lbl]] <- data.table(label_category = lbl, status = "CHECKED",
                              hits = paste(hits_str, collapse = " | "))
}
check1_dt <- rbindlist(check1)
print(check1_dt, row.names = FALSE)

# ============================================================
# 2. F-stage boundary — program dominance across fibrosis stages
# NOTE: "F2 switch" framing retired 2026-05-09; see paper_outline.md
# ============================================================
cat("\n--- Check 2: F-stage boundary per program ---\n")
meta <- fread(file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"))
nmf_m <- merge(nmf[, .(sample_id, dominant_program, dominant_program_code)],
               meta[, .(sample_id, fibrosis_stage)], by = "sample_id")
nmf_m <- nmf_m[!is.na(fibrosis_stage)]
f_dist <- nmf_m[, .N, by = .(dominant_program, fibrosis_stage)]
f_dist[, pct := N / sum(N) * 100, by = fibrosis_stage]
setorder(f_dist, dominant_program, fibrosis_stage)
print(dcast(f_dist, dominant_program ~ fibrosis_stage, value.var = "pct",
            fun.aggregate = function(x) round(mean(x), 1), fill = 0), row.names=FALSE)

# Per-program F-stage boundary ratio (F0-2 vs F3-4)
fib_code <- labels$program_code[labels$label_category == "Fibrogenic"]
stopifnot(length(fib_code) == 1)
if (length(fib_code) >= 1) {
  pc <- fib_code[1]
  lo <- nmf_m[fibrosis_stage <= 2 & dominant_program_code == pc, .N] /
        nmf_m[fibrosis_stage <= 2, .N]
  hi <- nmf_m[fibrosis_stage >= 3 & dominant_program_code == pc, .N] /
        nmf_m[fibrosis_stage >= 3, .N]
  cat(sprintf("  Fibrotic (%s) dominance: F0-2=%.1f%%, F3-4=%.1f%%, ratio=%.2f\n",
              pc, lo*100, hi*100, hi/lo))
}

# ============================================================
# 3. Female-biased DEG x Fibrotic-program overlap
# ============================================================
cat("\n--- Check 3: Female DEGs x Fibrotic program overlap ---\n")
sex_file <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_deg_classification.csv")
if (file.exists(sex_file)) {
  sex <- fread(sex_file)
  sex[, ens_base := sub("\\..*", "", gene)]
  sex <- merge(sex, atlas[, .(ens_base, human_symbol)], by = "ens_base", all.x = TRUE)
  female_cls <- if ("Female_biased" %in% sex$sex_class) "Female_biased" else "Female_specific"
  f_degs <- sex[sex_class == female_cls & !is.na(human_symbol), unique(human_symbol)]
  cat(sprintf("  %s genes: %d\n", female_cls, length(f_degs)))
  if (length(fib_code) >= 1) {
    fib_up <- markers[program_code == fib_code[1] & direction == "up" & padj < 0.05,
                       unique(human_symbol)]
    overlap <- length(intersect(f_degs, fib_up))
    universe <- unique(markers$human_symbol); universe <- universe[!is.na(universe)]
    n_fib <- length(intersect(fib_up, universe))
    n_fdeg <- length(intersect(f_degs, universe))
    a <- overlap; b <- n_fib - a; c_ <- n_fdeg - a; d <- length(universe) - a - b - c_
    ft <- fisher.test(matrix(c(a, b, max(c_,0), max(d,0)), nrow = 2))
    cat(sprintf("  Fibrotic-up ∩ female-biased: %d overlap (OR=%.2f, p=%.2e)\n",
                overlap, ft$estimate, ft$p.value))
  }
}

# ============================================================
# 4. Confounder-absence assertions (post comprehensive strip)
# ============================================================
cat("\n--- Check 4: confounder-absence in top-50 loadings per program ---\n")
confounder_report <- list()
for (pc in labels$program_code) {
  top_syms <- markers[program_code == pc & direction == "up" & padj < 0.05,
                      human_symbol][1:50]
  top_syms <- top_syms[!is.na(top_syms) & top_syms != ""]
  hits_y   <- top_syms[top_syms %in% c("BCORP1","CDY4P","ANOS2P","AGKP1",
                                        "KDM5D","USP9Y","UTY","RPS4Y1","DDX3Y")]
  hits_xi  <- top_syms[top_syms %in% c("XIST","TSIX")]
  hits_mt  <- top_syms[startsWith(top_syms, "MT-")]
  hits_rp  <- top_syms[grepl("^RPS[0-9]|^RPL[0-9]|^MRPS[0-9]|^MRPL[0-9]", top_syms)]
  hits_ig  <- top_syms[grepl("^IG[HKL][VJCD][0-9]", top_syms)]
  hits_hb  <- top_syms[top_syms %in% c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")]
  n_conf <- length(c(hits_y, hits_xi, hits_mt, hits_rp, hits_ig, hits_hb))
  confounder_report[[pc]] <- data.table(
    program = pc, n_confounders = n_conf,
    chrY = paste(hits_y, collapse=","),
    XIST = paste(hits_xi, collapse=","),
    MT = paste(hits_mt, collapse=","),
    RP = paste(hits_rp, collapse=","),
    IG = paste(hits_ig, collapse=","),
    HB = paste(hits_hb, collapse=","))
}
conf_dt <- rbindlist(confounder_report)
print(conf_dt[, .(program, n_confounders, chrY, XIST, MT, RP, IG, HB)], row.names = FALSE)
total_conf <- sum(conf_dt$n_confounders)
if (total_conf == 0) {
  cat("  PASS: zero confounder hits across all programs\n")
} else {
  cat(sprintf("  FAIL: %d confounder hits across programs\n", total_conf))
}

# ============================================================
# 5. Hardcoded S1/S2 remnants
# ============================================================
cat("\n--- Check 4: hardcoded S1/S2 remnants (grep audit) ---\n")
audit <- system(
  paste0("cd ", BASE, " && grep -RnE '\"S1\"|\"S2\"' RNA-seq/ scripts/figures/ --include='*.R' | ",
         "grep -vE 'legacy|deprecated|DEPRECATED|_legacy|# |aliased|LEGACY|nmf_subtype|S1 =|\\.R#|sex_class|S2_upregulated|S1_upregulated'"),
  intern = TRUE, ignore.stderr = TRUE)
if (length(audit) == 0) {
  cat("  CLEAN: no non-legacy S1/S2 strings found\n")
} else {
  cat(sprintf("  FOUND %d non-legacy hits (may be intentional or stragglers):\n", length(audit)))
  for (h in head(audit, 20)) cat("    ", h, "\n")
}

cat("\n=== Sanity checks complete ===\n")
