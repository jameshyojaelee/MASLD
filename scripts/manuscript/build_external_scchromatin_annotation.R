#!/usr/bin/env Rscript
# build_external_scchromatin_annotation.R
#
# Join the external Elison/Gaulton et al. 2025 single-cell-chromatin / cell-type
# regulatory calls to our FROZEN expression-QTL + disease-state evidence-class
# table (2026-07-15-r2) and quantify, per gene, how the single-cell chromatin
# layer relates to our genetic + disease-state classes. ADDITIVE / external
# annotation only -- does NOT mutate the frozen release and is NOT an input to the
# convergence/heuristic score. See data/external/elison2025_sc_chromatin/README.md
# and the number-verification separation.
#
# Analog of scripts/manuscript/build_external_pqtl_annotation.R (Gobeil pQTL), one
# regulatory layer over. Our per-gene values are pulled LIVE from the frozen table
# (never hardcoded); every Elison value comes from the curated, PDF-verified CSV.
#
# Out: RNA-seq/results/multi_evidence/external_scchromatin/scchromatin_external_annotation.tsv
#      RNA-seq/results/multi_evidence/external_scchromatin/scchromatin_external_summary.txt
#      RNA-seq/results/multi_evidence/external_scchromatin/figS2Q_source.tsv   (figure-ready)

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
elison_f <- file.path(BASE, "data/external/elison2025_sc_chromatin/elison2025_sc_chromatin_calls.csv")
frozen_f <- file.path(BASE, "RNA-seq/results/manuscript_release/2026-07-15-r2/evidence_class_table.tsv")
out_dir  <- file.path(BASE, "RNA-seq/results/multi_evidence/external_scchromatin")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(elison_f), file.exists(frozen_f))
x <- fread(elison_f)
e <- fread(frozen_f)

# --- pull our frozen per-gene expression + disease-state evidence ---------------
keep <- c("symbol", "primary_evidence_class", "sensitivity_evidence_class",
          "genetic_trait_scope", "bulk_logFC", "bulk_treat_fdr", "in_treat_deg",
          "max_susie_pp4_all", "max_abf_pp4_all",
          "max_susie_pp4_direct", "max_abf_pp4_direct",
          "max_susie_pp4_enzyme", "max_abf_pp4_enzyme")
e <- e[, ..keep]
setnames(e, keep,
         c("symbol", "our_primary_class", "our_sensitivity_class",
           "our_genetic_trait_scope", "our_bulk_logFC", "our_bulk_treat_fdr",
           "our_in_treat_deg", "our_susie_all", "our_abf_all",
           "our_susie_direct", "our_abf_direct",
           "our_susie_enzyme", "our_abf_enzyme"))

m <- merge(x, e, by.x = "gene", by.y = "symbol", all.x = TRUE, sort = FALSE)

num <- function(v) suppressWarnings(as.numeric(v))
GATE <- 0.5   # our SuSiE/ABF "genetic" gate (PP.H4 > 0.5)

# Coalesce missing PP4 to 0 (not tested / not colocalizing = below-gate, NOT NA)
sA <- fcoalesce(num(m$our_susie_all), 0); aA <- fcoalesce(num(m$our_abf_all), 0)
sD <- fcoalesce(num(m$our_susie_direct), 0); aD <- fcoalesce(num(m$our_abf_direct), 0)
m[, our_genetic_any    := (sA > GATE) | (aA > GATE)]          # any-trait (incl enzyme) coloc
m[, our_genetic_direct := (sD > GATE) | (aD > GATE)]          # direct-MASLD coloc (the honest anchor)
m[, our_is_deg         := our_in_treat_deg == TRUE]

# --- relationship call (uses our frozen values + their finding type) -----------
classify <- function(r) {
  ft <- r$elison_finding_type
  if (isTRUE(ft == "boundary_absent"))  return("boundary_lipid_canon")            # absent both maps
  if (isTRUE(ft == "cell_of_action")) {
    if (isTRUE(r$our_genetic_any)) return("gap_fill_celltype_resolved")           # we colocalize; they assign cell type
    return("gap_fill_they_stronger")                                              # our coloc null; their caQTL detects it
  }
  if (isTRUE(ft == "colocalized_locus")) return("co_localized_both_maps")         # EFHD1 (direction caveat)
  if (isTRUE(ft == "featured_target")) {
    if (isTRUE(r$our_genetic_any)) return("we_supply_coloc")                      # KRT8/PPP1R3B: loop-only for them, coloc for us
    return("loop_only_we_null")                                                   # CEBPA/XBP1
  }
  if (grepl("^grn_tf", ft)) {                                                     # SCENIC+ TF-activity
    if (isTRUE(r$our_is_deg)) return("corroborates_disease_state")               # our DEG agrees
    if (isTRUE(r$our_genetic_any)) return("corroborates_genetic")                # our coloc agrees
    return("corroborates_regulon_only")                                          # TF-activity-only agreement (e.g. THRB/HNF4A)
  }
  return("unclassified")
}
m[, relationship := vapply(seq_len(.N), function(i) classify(m[i]), character(1))]

# axis / trait-anchor caveat (honest reporting)
m[, caveat := fifelse(gene == "EFHD1", "direction/axis differ: their hep-eQTL risk-UP vs our bulk DE DOWN -> co-localized locus, NOT sign-match",
              fifelse(gene == "RORA",  "co-nomination NOT direction-match (their GRN-activity UP vs our disease DE DOWN)",
              fifelse(gene == "THRB",  "GGT/enzyme-anchored; direct-MASLD coloc null; TF-activity not a MASH genetic hit",
              fifelse(our_susie_direct < GATE & our_abf_direct < GATE & our_genetic_any == TRUE,
                      "headline coloc is enzyme/off-liver-injury-trait anchored; direct-MASLD sub-gate",
                      ""))))]

ordv <- c("gene", "elison_locus", "elison_finding_type", "elison_layer", "elison_celltype",
          "elison_lead_variant", "elison_direction", "our_rescue_layer",
          "our_primary_class", "our_sensitivity_class", "our_genetic_trait_scope",
          "our_susie_all", "our_abf_all", "our_susie_direct", "our_abf_direct",
          "our_susie_enzyme", "our_abf_enzyme", "our_is_deg", "our_bulk_logFC",
          "our_bulk_treat_fdr", "our_genetic_any", "our_genetic_direct",
          "relationship", "caveat", "elison_detail", "source_note")
m <- m[, intersect(ordv, names(m)), with = FALSE]
out_tsv <- file.path(out_dir, "scchromatin_external_annotation.tsv")
fwrite(m, out_tsv, sep = "\t")

# figure-ready long table (FigS2Q): corroboration set + boundary set
fig <- m[, .(gene, elison_finding_type, elison_celltype, relationship,
             our_primary_class, our_susie_direct, our_abf_direct,
             our_susie_enzyme, our_abf_enzyme, our_susie_all, our_abf_all,
             our_is_deg, our_bulk_treat_fdr, our_bulk_logFC, our_rescue_layer, caveat)]
fwrite(fig, file.path(out_dir, "figS2Q_source.tsv"), sep = "\t")

# --- headline counts -----------------------------------------------------------
sink(file.path(out_dir, "scchromatin_external_summary.txt"))
cat("External single-cell-chromatin (Elison/Gaulton 2025) vs our frozen maps\n")
cat("======================================================================\n\n")
cat("Rows:", nrow(m), " genes\n\n")
cat("Per-gene relationship to our expression-genetic + disease-state maps:\n")
print(m[, .(gene, elison_finding_type, our_primary_class,
            susieDir = round(fcoalesce(num(our_susie_direct),0),3),
            abfDir  = round(fcoalesce(num(our_abf_direct),0),3),
            deg = our_is_deg, relationship)])
cat("\nHeadline counts:\n")
tallies <- m[, .N, by = relationship][order(-N)]
print(tallies)
cat("\nCorroborated at our GENETIC map (direct-MASLD coloc > 0.5):\n  ",
    paste(m[our_genetic_direct == TRUE]$gene, collapse = ", "), "\n", sep = "")
cat("Corroborated at our DISEASE-STATE map (canonical DEG):\n  ",
    paste(m[our_is_deg == TRUE]$gene, collapse = ", "), "\n", sep = "")
cat("Cell-of-action gap-fills (they assign cell type / detect where we are null):\n  ",
    paste(m[grepl("gap_fill", relationship)]$gene, collapse = ", "), "\n", sep = "")
cat("Mutual boundary (lipid canon absent in Elison AND null on our expression-COLOC):\n  ",
    paste(m[relationship == "boundary_lipid_canon"]$gene, collapse = ", "), "\n", sep = "")
sink()

cat("Wrote:", out_tsv, "\n")
cat(readLines(file.path(out_dir, "scchromatin_external_summary.txt")), sep = "\n")
