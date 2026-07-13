#!/usr/bin/env Rscript
# build_external_pqtl_annotation.R
#
# Join the external Gobeil et al. 2026 liver-pQTL causal calls to our FROZEN
# expression-QTL evidence-class table (2026-07-10-r1) and quantify, per gene,
# how the protein-genetic layer relates to our expression-genetic + disease-state
# classes. ADDITIVE / external annotation only -- this does NOT mutate the frozen
# release and is NOT an input to the convergence/heuristic score. See
# data/external/gobeil2026_liver_pqtl/README.md and the number-verification firewall.
#
# Out: RNA-seq/results/multi_evidence/external_pqtl/pqtl_external_annotation.tsv
#      RNA-seq/results/multi_evidence/external_pqtl/pqtl_external_summary.txt

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
gobeil_f <- file.path(BASE, "data/external/gobeil2026_liver_pqtl/gobeil2026_liver_pqtl_causal.csv")
frozen_f <- file.path(BASE, "RNA-seq/results/manuscript_release/2026-07-10-r1/evidence_class_table.tsv")
out_dir  <- file.path(BASE, "RNA-seq/results/multi_evidence/external_pqtl")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(gobeil_f), file.exists(frozen_f))
g <- fread(gobeil_f)
e <- fread(frozen_f)

# --- pull our frozen per-gene expression evidence (symbol = col 'symbol') ------
keep <- c("symbol", "primary_evidence_class", "sensitivity_evidence_class",
          "genetic_confidence", "bulk_logFC", "bulk_treat_fdr", "in_treat_deg",
          "max_susie_pp4_all", "max_abf_pp4_all")
e <- e[, ..keep]
setnames(e, keep,
         c("symbol", "our_primary_class", "our_sensitivity_class",
           "our_genetic_confidence", "our_bulk_logFC", "our_bulk_treat_fdr",
           "our_in_treat_deg", "our_max_susie_pp4", "our_max_abf_pp4"))

m <- merge(g, e, by.x = "gene", by.y = "symbol", all.x = TRUE, sort = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))
GATE <- 0.5   # our SuSiE/ABF "genetic" gate (outline: PP.H4 > 0.5)
PQTL_GATE <- 0.80  # Gobeil's colocalization prioritization threshold

# Coalesce missing PP4 to 0 (a gene not tested / not colocalizing is below-gate,
# NOT NA) so the boolean gates and their tallies are not poisoned by NA.
susie0 <- fcoalesce(num(m$our_max_susie_pp4), 0)
abf0   <- fcoalesce(num(m$our_max_abf_pp4), 0)
pq0    <- fcoalesce(num(m$gobeil_pqtl_coloc_pph4), 0)
m[, our_genetic_primary     := susie0 > GATE]
m[, our_genetic_sensitivity := (susie0 > GATE) | (abf0 > GATE)]
m[, gobeil_pqtl_hit         := pq0 >= PQTL_GATE]

# --- concordance call ----------------------------------------------------------
classify <- function(r) {
  if (isTRUE(r$layer_relationship == "direction_discordant")) return("direction_discordant")
  if (isTRUE(r$our_primary_class == "disease_state_only"))    return("reactive_agreement")
  if (isTRUE(r$our_genetic_primary))                          return("corroborated_primary_susie")
  if (isTRUE(r$our_genetic_sensitivity))                      return("corroborated_sensitivity_abf")
  if (isTRUE(r$gobeil_pqtl_hit))                              return("expression_miss")   # pQTL causal, we miss on eQTL
  return("subthreshold_both")
}
m[, concordance_call := vapply(seq_len(.N), function(i) classify(m[i]), character(1))]

# reclassification if a protein-QTL channel (Gobeil pQTL coloc >= 0.80) were added
m[, reclass_if_pqtl_added := fifelse(
  our_primary_class == "neither" & gobeil_pqtl_hit,
  "neither -> genetic_only(pQTL)",
  fifelse(our_primary_class == "genetic_only" & gobeil_pqtl_hit,
          "genetic_only (protein-QTL confirmed)",
  fifelse(our_primary_class == "disease_state_only" & gobeil_pqtl_hit,
          "disease_state_only + protein-QTL (mixed)",
          as.character(our_primary_class))))]

ordv <- c("gene", "gobeil_coloc_prioritized", "gobeil_concordant_set",
          "gobeil_pqtl_coloc_pph4", "gobeil_pqtl_mr_beta", "gobeil_pqtl_risk_direction",
          "gobeil_eqtl_coloc_pph4", "layer_relationship",
          "our_primary_class", "our_sensitivity_class", "our_max_susie_pp4",
          "our_max_abf_pp4", "our_in_treat_deg", "our_bulk_logFC",
          "our_genetic_primary", "our_genetic_sensitivity", "gobeil_pqtl_hit",
          "concordance_call", "reclass_if_pqtl_added")
m <- m[, intersect(ordv, names(m)), with = FALSE]
out_tsv <- file.path(out_dir, "pqtl_external_annotation.tsv")
fwrite(m, out_tsv, sep = "\t")

# --- headline counts -----------------------------------------------------------
prio <- m[gobeil_coloc_prioritized == TRUE]           # Gobeil's 5 coloc>=0.80 genes
sink(file.path(out_dir, "pqtl_external_summary.txt"))
cat("External liver-pQTL (Gobeil 2026) vs our frozen expression-QTL map\n")
cat("=================================================================\n\n")
cat(sprintf("Gobeil coloc-prioritized causal genes (pQTL PPH4 >= %.2f): %d\n",
            PQTL_GATE, nrow(prio)))
cat("  ", paste(prio$gene, collapse = ", "), "\n\n", sep = "")
cat("Per-gene concordance with our expression map:\n")
print(m[, .(gene, gobeil_pqtl_coloc_pph4, our_primary_class,
            our_max_susie_pp4, our_max_abf_pp4, concordance_call,
            reclass_if_pqtl_added)])
cat("\nHeadline counts (over Gobeil coloc-prioritized genes):\n")
cat(sprintf("  corroborated by our PRIMARY SuSiE (>%.1f): %d  [%s]\n", GATE,
            prio[our_genetic_primary == TRUE, .N],
            paste(prio[our_genetic_primary == TRUE]$gene, collapse=", ")))
cat(sprintf("  corroborated only in SENSITIVITY (SuSiE|ABF>%.1f): %d  [%s]\n", GATE,
            prio[our_genetic_sensitivity == TRUE & our_genetic_primary == FALSE, .N],
            paste(prio[our_genetic_sensitivity == TRUE & our_genetic_primary == FALSE]$gene, collapse=", ")))
cat(sprintf("  reactive in BOTH maps (we: disease_state_only): %d  [%s]\n",
            prio[our_primary_class == "disease_state_only", .N],
            paste(prio[our_primary_class == "disease_state_only"]$gene, collapse=", ")))
cat(sprintf("  direction-discordant (eQTL vs pQTL sign flip): %d  [%s]\n",
            prio[layer_relationship == "direction_discordant", .N],
            paste(prio[layer_relationship == "direction_discordant"]$gene, collapse=", ")))
cat(sprintf("  MISSED by our expression map (not genetic even in sensitivity, not DEG): %d  [%s]\n",
            prio[our_genetic_sensitivity == FALSE & our_primary_class != "disease_state_only", .N],
            paste(prio[our_genetic_sensitivity == FALSE & our_primary_class != "disease_state_only"]$gene, collapse=", ")))
cat("\nReclassification the protein-QTL channel would induce (primary class):\n")
print(m[reclass_if_pqtl_added != our_primary_class,
        .(gene, our_primary_class, reclass_if_pqtl_added)])
sink()

cat("Wrote:", out_tsv, "\n")
cat(readLines(file.path(out_dir, "pqtl_external_summary.txt")), sep = "\n")
