#!/usr/bin/env Rscript
# ============================================================================
# q2_threshold_free_fgsea.R   (Vacca benchmark — Q2)
#
# THRESHOLD-FREE replacement for the cutoff-sensitive 172/1323 overlap.
# Rank ALL genes in concordance_atlas_unified.csv by a CONTINUOUS cross-species
# concordance statistic, and fgsea-test whether Vacca's 951-gene DSEA signature
# (MOESM4 'Gene_used_in_DSEA'==1, human_symbol) is enriched at the top.
#
# Rankings tested:
#   (1) n_concordant            (0-4 integer, # diets concordant)
#   (2) translatability_score   (continuous 0-0.764)
#   (3) signed_nconc = sign(mean_h_lfc) * n_concordant
#   (+) mean_h_lfc alone (continuous human effect, as a reference contrast)
#
# fgsea needs a UNIQUE-named numeric ranking vector; ties jittered deterministically.
# We also report a hypergeometric/Fisher sanity at the matched gene universe.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl); library(fgsea) })
set.seed(42)
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(CS, "vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── 1. Vacca 951-gene DSEA signature (human symbols) ─────────────────────────
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"),
                               sheet = "Table S4"))
dsea_col <- grep("Gene_used_in_DSEA", names(s4), value = TRUE)[1]
sym_col  <- grep("^GeneSymbol$", names(s4), value = TRUE)[1]
if (is.na(sym_col)) sym_col <- grep("GeneSymbol", names(s4), value = TRUE)[1]
cat(sprintf("MOESM4: %d rows; DSEA col='%s'; symbol col='%s'\n",
            nrow(s4), dsea_col, sym_col))
s4[, dsea := suppressWarnings(as.numeric(get(dsea_col)))]
vacca_sig <- unique(toupper(trimws(s4[dsea == 1][[sym_col]])))
vacca_sig <- vacca_sig[!is.na(vacca_sig) & vacca_sig != ""]
cat(sprintf("Vacca DSEA signature genes (in_dsea==1, unique upper symbol): %d\n",
            length(vacca_sig)))

# ── 2. Our concordance atlas — continuous ranking statistics ─────────────────
a <- fread(file.path(CS, "concordance_atlas_unified.csv"))
a[, sym := toupper(trimws(human_symbol))]
a <- a[!is.na(sym) & sym != ""]
# collapse duplicate symbols: keep the row with the strongest concordance
setorder(a, -n_concordant, -translatability_score)
a <- a[!duplicated(sym)]
cat(sprintf("Atlas genes (unique upper symbol): %d\n", nrow(a)))

# signed statistic: human direction × # concordant diets
a[, signed_nconc := sign(mean_h_lfc) * n_concordant]

# overlap of Vacca sig with our atlas universe (this is the testable set)
in_universe <- intersect(vacca_sig, a$sym)
cat(sprintf("Vacca sig genes present in our atlas universe: %d / %d (%.1f%%)\n",
            length(in_universe), length(vacca_sig),
            100 * length(in_universe) / length(vacca_sig)))

# ── 3. fgsea over each continuous ranking ────────────────────────────────────
# fgsea requires a uniquely-named numeric vector. Add a tiny deterministic
# jitter to break ties so ranking is well-defined (seeded; magnitude << min gap).
make_rank <- function(dt, statcol) {
  v <- dt[[statcol]]
  # deterministic jitter proportional to a hashed gene order, tiny magnitude
  jit <- (seq_len(nrow(dt)) / nrow(dt)) * 1e-6
  setNames(v + jit, dt$sym)
}
pathways <- list(Vacca_DSEA951 = in_universe)

rankings <- c("n_concordant", "translatability_score", "signed_nconc", "mean_h_lfc")
res <- rbindlist(lapply(rankings, function(rk) {
  stats <- sort(make_rank(a, rk), decreasing = TRUE)
  fg <- fgsea(pathways = pathways, stats = stats,
              minSize = 5, maxSize = length(in_universe) + 10,
              scoreType = "std", nPermSimple = 10000, eps = 0)
  data.table(ranking = rk, NES = fg$NES, pval = fg$pval, padj = fg$padj,
             ES = fg$ES, size = fg$size,
             leadingEdge_n = vapply(fg$leadingEdge, length, integer(1)))
}))

# Also a one-sided "positive" scoreType for the unsigned magnitude rankings
# (n_concordant / translatability are non-negative → top = most concordant;
#  scoreType="pos" asks specifically: is the sig enriched at the HIGH end?)
res_pos <- rbindlist(lapply(c("n_concordant", "translatability_score"), function(rk) {
  stats <- sort(make_rank(a, rk), decreasing = TRUE)
  fg <- fgsea(pathways = pathways, stats = stats,
              minSize = 5, maxSize = length(in_universe) + 10,
              scoreType = "pos", nPermSimple = 10000, eps = 0)
  data.table(ranking = paste0(rk, "_posScore"), NES = fg$NES, pval = fg$pval,
             padj = fg$padj, ES = fg$ES, size = fg$size,
             leadingEdge_n = vapply(fg$leadingEdge, length, integer(1)))
}))
res <- rbind(res, res_pos, fill = TRUE)

# ── 4. Fisher sanity at the matched universe (n_concordant>=1 as "concordant") ─
a[, our_conc := n_concordant >= 1]
a[, is_vacca := sym %in% vacca_sig]
ct <- table(our_conc = a$our_conc, is_vacca = a$is_vacca)
ft <- fisher.test(ct)
cat("\n=== Fisher sanity (n_concordant>=1 vs Vacca-sig membership) ===\n")
print(ct)
cat(sprintf("Fisher OR = %.2f, p = %.2e\n", ft$estimate, ft$p.value))

cat("\n=== THRESHOLD-FREE fgsea: Vacca DSEA-951 enrichment at top of our concordance rankings ===\n")
print(res[order(-NES)])

fwrite(res, file.path(OUT, "q2_threshold_free_fgsea.csv"))
fisher_dt <- data.table(test = "fisher_nconc1_vs_vaccasig",
                        OR = as.numeric(ft$estimate), p = ft$p.value,
                        n_both = ct["TRUE","TRUE"], n_vacca_universe = length(in_universe))
fwrite(fisher_dt, file.path(OUT, "q2_threshold_free_fisher.csv"))

cat(sprintf("\nWrote: %s\n", file.path(OUT, "q2_threshold_free_fgsea.csv")))
cat(strrep("=", 78), "\n")
cat("VERDICT: positive NES + padj<0.05 across rankings ⇒ Vacca's signature is\n")
cat("  enriched among OUR most cross-species-concordant genes WITHOUT any cutoff\n")
cat("  ⇒ the 172-overlap was real, not a cutoff artifact.\n")
cat(strrep("=", 78), "\n")
