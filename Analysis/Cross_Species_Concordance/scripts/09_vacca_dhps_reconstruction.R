#!/usr/bin/env Rscript
# ============================================================================
# 09_vacca_dhps_reconstruction.R
#
# The gene-level fgsea arm does NOT reproduce Vacca's DHPS (script 08: Spearman
# 0.06 even on their own data). DHPS is a DSEA/gep2pep construct over BOTH genes
# AND KEGG pathways in pathway-expression-profile (PEP) space. Here we REVERSE-
# ENGINEER the DHPS from Vacca's OWN provided intermediates (no gep2pep rebuild):
#   - per-model pathway NES (MOESM6, the DRP arm / model PEP)
#   - per-model %-agreement with human (MOESM7)
#   - per-model gene-arm fgsea (from script 08 output)
# and test which reconstruction reproduces published DHPS (MOESM9, ~41 models).
# A candidate with Spearman >= ~0.6 vs published DHPS = the DHPS construction we
# can then apply to OUR diets.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl); library(fgsea) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
OUT   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
set.seed(42)
nrm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
sp  <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))

# ── Published DHPS (MOESM9) ──────────────────────────────────────────────────
m9 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM9_ESM.xlsx"),
                               col_names = FALSE, skip = 2))
setnames(m9, 1:10, c("model","diet_group","mP","mH","m_DHPS","mM","fP","fH","f_DHPS","fM"))
m9 <- m9[!is.na(model)]; for (cc in c("m_DHPS","f_DHPS")) m9[[cc]] <- as.numeric(m9[[cc]])
m9[, `:=`(key = nrm(model), dhps_mean = (as.numeric(m_DHPS)+as.numeric(f_DHPS))/2)]

# ── MOESM6 pathway NES: human reference + per-model PEPs ──────────────────────
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"), sheet = "NES"))
setnames(p6, 1, "pathway")
human_cols <- grep("UCAM|EPoS", names(p6), value = TRUE)
model_cols <- setdiff(names(p6)[-1], human_cols)
for (cc in c(human_cols, model_cols)) suppressWarnings(p6[[cc]] <- as.numeric(p6[[cc]]))
# human reference PEPs (test several contrasts)
href <- list(
  UCAM_Sev = p6[["UCAM/VCU: Severe vs Mild"]],
  EPoS_Sev = p6[["EPoS: Severe vs Mild"]],
  UCAM_Mod = p6[["UCAM/VCU: Moderate vs Mild"]],
  Sev_mean = rowMeans(cbind(p6[["UCAM/VCU: Severe vs Mild"]], p6[["EPoS: Severe vs Mild"]]), na.rm = TRUE),
  Mild_ctrl = p6[["UCAM/VCU: Mild vs Control"]])

# per-model PEP-space proximity = correlation(model pathway-NES, human pathway-NES)
recon <- rbindlist(lapply(model_cols, function(mc) {
  mv <- p6[[mc]]
  row <- data.table(model = mc); row[, key := nrm(mc)]
  for (h in names(href)) {
    row[[paste0("pepP_", h)]] <- suppressWarnings(cor(mv, href[[h]], method = "pearson", use = "complete.obs"))
    row[[paste0("pepS_", h)]] <- suppressWarnings(cor(mv, href[[h]], method = "spearman", use = "complete.obs"))
  }
  # pathway-fgsea: human up/down pathways (by Sev_mean) enriched in model NES ranking
  hsev <- href$Sev_mean
  up_p   <- p6$pathway[is.finite(hsev) & hsev >  1]
  down_p <- p6$pathway[is.finite(hsev) & hsev < -1]
  st <- setNames(mv, p6$pathway); st <- st[is.finite(st)]
  if (length(st) > 20 && length(up_p) > 2 && length(down_p) > 2) {
    fg <- as.data.table(suppressWarnings(fgsea(list(up = up_p, down = down_p), st,
                                               scoreType = "std", nPermSimple = 5000)))
    nu <- fg[pathway=="up", NES]; nd <- fg[pathway=="down", NES]
    row[["pathArm"]] <- (ifelse(length(nu),nu,0) - ifelse(length(nd),nd,0))/2
  } else row[["pathArm"]] <- NA_real_
  row
}))

# ── MOESM7 %-agreement (DEGs_All) ────────────────────────────────────────────
p7 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM7_ESM.xlsx"), sheet = "% of agreement"))
mcol7 <- names(p7)[1]; setnames(p7, mcol7, "model")
agg_col <- grep("DEGs__All|DEGs_All|All", names(p7), value = TRUE)[1]
p7[, key := nrm(model)]
recon <- merge(recon, p7[, .(key, pct_agree_all = as.numeric(get(agg_col)))], by = "key", all.x = TRUE)

# ── gene-arm from script 08 (if present) ─────────────────────────────────────
d08 <- file.path(OUT, "diagnostic_D2D3_per_model.csv")
if (file.exists(d08)) {
  g <- fread(d08)[, .(key = nrm(model), geneArm = d3_vacca_sig)]
  recon <- merge(recon, g, by = "key", all.x = TRUE)
}

# ── Merge published DHPS + score every candidate ─────────────────────────────
mg <- merge(recon, m9[, .(key, diet_group, m_DHPS, f_DHPS, dhps_mean)], by = "key")
mg <- mg[!grepl("^R[-.]", model)]   # mouse only (drop rat CDAA)
cat(sprintf("Models matched (mouse): %d\n\n", nrow(mg)))

cand <- setdiff(names(mg), c("key","model","diet_group","m_DHPS","f_DHPS","dhps_mean"))
score <- rbindlist(lapply(cand, function(cc) data.table(
  candidate = cc,
  spearman_metab = sp(mg[[cc]], mg$m_DHPS),
  spearman_fibro = sp(mg[[cc]], mg$f_DHPS),
  spearman_mean  = sp(mg[[cc]], mg$dhps_mean),
  n = sum(is.finite(mg[[cc]]))) ))
setorder(score, -spearman_mean)
cat("=== Reconstruction candidate vs PUBLISHED DHPS (Spearman across mouse models) ===\n")
print(score)
cat("\n  >> Any candidate with |Spearman| >= 0.6 = the DHPS construction to adopt.\n")
fwrite(score, file.path(OUT, "dhps_reconstruction_scores.csv"))
fwrite(mg[order(-dhps_mean)], file.path(OUT, "dhps_reconstruction_per_model.csv"))

best <- score[which.max(abs(spearman_mean))]
cat(sprintf("\n  BEST candidate: %s  (Spearman vs DHPS_mean = %.2f)\n",
            best$candidate, best$spearman_mean))
cat(strrep("=", 78), "\n")
