#!/usr/bin/env Rscript
# ============================================================================
# q1_twoarm_pep_dhps.R   (Vacca follow-up, Q1 KEY)
#
# Vacca's DHPS is TWO-ARMED: a METABOLIC arm (MOESM9 col5) and a FIBROTIC arm
# (MOESM9 col9). Question: does a two-arm transcriptomic PEP (pathway-NES
# correlation, on Vacca's OWN ~37 mouse models, MOESM6) reproduce BOTH arms?
#
#   metabolic proximity = PEP corr( model pathway-NES , human 'Mild vs Control' )
#   fibrotic  proximity = PEP corr( model pathway-NES , human mean Severe-vs-Mild )
#                         where mean Severe-vs-Mild = rowMeans(UCAM Sev, EPoS Sev)
#
# Then Spearman:
#   metabolic proximity   vs published METABOLIC-arm DHPS (MOESM9 col5)
#   fibrotic proximity    vs published FIBROTIC-arm  DHPS (MOESM9 col9)
# CROSS checks included (does fibrotic PEP also track metabolic DHPS, etc.).
# Drop rat models (^R-). Env: rnaseq.
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
OUT   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
RES   <- file.path(BASE, "results/vacca_benchmark"); dir.create(RES, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
nrm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
sp  <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
sp.p <- function(a, b) suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE)$p.value)

# ── Published two-arm DHPS (MOESM9) ──────────────────────────────────────────
m9 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM9_ESM.xlsx"),
                               col_names = FALSE, skip = 2))
setnames(m9, 1:10, c("model","diet_group","mP","mH","m_DHPS","mM","fP","fH","f_DHPS","fM"))
m9 <- m9[!is.na(model)]
m9[, `:=`(key = nrm(model),
          m_DHPS = as.numeric(m_DHPS),   # METABOLIC-arm DHPS (col5)
          f_DHPS = as.numeric(f_DHPS))]  # FIBROTIC-arm  DHPS (col9)

# ── MOESM6 pathway NES: human reference + per-model PEPs ──────────────────────
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"), sheet = "NES"))
setnames(p6, 1, "pathway")
human_cols <- grep("UCAM|VCU|EPoS", names(p6), value = TRUE)
model_cols <- setdiff(names(p6)[-1], human_cols)
for (cc in c(human_cols, model_cols)) suppressWarnings(p6[[cc]] <- as.numeric(p6[[cc]]))

href_metab <- p6[["UCAM/VCU: Mild vs Control"]]                  # METABOLIC reference
href_fibro <- rowMeans(cbind(p6[["UCAM/VCU: Severe vs Mild"]],   # FIBROTIC reference
                             p6[["EPoS: Severe vs Mild"]]), na.rm = TRUE)

cat(sprintf("MOESM6 pathways: %d ; model cols: %d\n", nrow(p6), length(model_cols)))
cat(sprintf("Human metabolic ref finite: %d ; fibrotic ref finite: %d\n",
            sum(is.finite(href_metab)), sum(is.finite(href_fibro))))

# per-model two-arm PEP proximity (Pearson primary, Spearman robustness)
prox <- rbindlist(lapply(model_cols, function(mc) {
  mv <- p6[[mc]]
  dt <- data.table(model = mc,
    metab_pep_P = suppressWarnings(cor(mv, href_metab, method = "pearson",  use = "complete.obs")),
    fibro_pep_P = suppressWarnings(cor(mv, href_fibro, method = "pearson",  use = "complete.obs")),
    metab_pep_S = suppressWarnings(cor(mv, href_metab, method = "spearman", use = "complete.obs")),
    fibro_pep_S = suppressWarnings(cor(mv, href_fibro, method = "spearman", use = "complete.obs")))
  dt[, key := nrm(mc)]
  dt
}))

# ── Merge published DHPS, mouse models only ──────────────────────────────────
mg <- merge(prox, m9[, .(key, diet_group, m_DHPS, f_DHPS)], by = "key")
mg <- mg[!grepl("^R[-.]", model)]   # drop rat
cat(sprintf("\nMouse models matched (MOESM6 PEP x MOESM9 DHPS): %d\n", nrow(mg)))

# ── Two-arm matched + cross Spearmans ────────────────────────────────────────
out <- data.table(
  test = c("MATCHED: metabolic-PEP(Pearson) vs metabolic-DHPS",
           "MATCHED: fibrotic-PEP(Pearson)  vs fibrotic-DHPS",
           "MATCHED: metabolic-PEP(Spearman) vs metabolic-DHPS",
           "MATCHED: fibrotic-PEP(Spearman)  vs fibrotic-DHPS",
           "CROSS:   metabolic-PEP(Pearson) vs fibrotic-DHPS",
           "CROSS:   fibrotic-PEP(Pearson)  vs metabolic-DHPS",
           "ref:     PEP-proximity arms inter-corr (metab vs fibro PEP)",
           "ref:     published DHPS arms inter-corr (m_DHPS vs f_DHPS)"),
  spearman = c(
    sp(mg$metab_pep_P, mg$m_DHPS),
    sp(mg$fibro_pep_P, mg$f_DHPS),
    sp(mg$metab_pep_S, mg$m_DHPS),
    sp(mg$fibro_pep_S, mg$f_DHPS),
    sp(mg$metab_pep_P, mg$f_DHPS),
    sp(mg$fibro_pep_P, mg$m_DHPS),
    sp(mg$metab_pep_P, mg$fibro_pep_P),
    sp(mg$m_DHPS,      mg$f_DHPS)),
  p_value = c(
    sp.p(mg$metab_pep_P, mg$m_DHPS),
    sp.p(mg$fibro_pep_P, mg$f_DHPS),
    sp.p(mg$metab_pep_S, mg$m_DHPS),
    sp.p(mg$fibro_pep_S, mg$f_DHPS),
    sp.p(mg$metab_pep_P, mg$f_DHPS),
    sp.p(mg$fibro_pep_P, mg$m_DHPS),
    sp.p(mg$metab_pep_P, mg$fibro_pep_P),
    sp.p(mg$m_DHPS,      mg$f_DHPS)),
  n = nrow(mg))

cat("\n=== Q1: TWO-ARM transcriptomic PEP vs published two-arm DHPS ===\n")
print(out, digits = 3)

fwrite(out, file.path(RES, "q1_twoarm_pep_dhps.csv"))
fwrite(mg[order(-f_DHPS)], file.path(RES, "q1_twoarm_pep_per_model.csv"))

mP <- out[test=="MATCHED: metabolic-PEP(Pearson) vs metabolic-DHPS"]
fP <- out[test=="MATCHED: fibrotic-PEP(Pearson)  vs fibrotic-DHPS"]
cat(sprintf("\n>> METABOLIC arm reproduced: Spearman = %.3f (p=%.2g, n=%d)\n", mP$spearman, mP$p_value, mP$n))
cat(sprintf(">> FIBROTIC  arm reproduced: Spearman = %.3f (p=%.2g, n=%d)\n", fP$spearman, fP$p_value, fP$n))
cat("\nWrote:", file.path(RES, "q1_twoarm_pep_dhps.csv"), "\n")
cat(strrep("=", 78), "\n")
