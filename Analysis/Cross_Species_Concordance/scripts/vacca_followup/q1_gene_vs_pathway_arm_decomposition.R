#!/usr/bin/env Rscript
# ============================================================================
# q1_gene_vs_pathway_arm_decomposition.R
#
# Q1: WHY does gene-level reproduce Vacca's DHPS at only Spearman ~0.06 (script 07/08)
# but pathway-level at ~0.73 (scripts 09/10)?
#
# We compute, for every Vacca mouse model (~33), BOTH arms from source:
#   GENE-ARM  = DSEA-analog fgsea of the 951-gene human progression signature in
#               each model's per-gene L2FC ranking (MOESM4).  (NES_up - NES_down)/2.
#   PATH-ARM  = PEP correlation: Pearson(model pathway-NES vector, human pathway-NES
#               vector) over shared KEGG pathways (MOESM6, Severe-vs-Mild ref).
# Then:
#   (a) correlate gene-arm vs path-arm across models  -> are they even measuring
#       the same thing?
#   (b) per-model divergence = rank(path-arm) - rank(gene-arm); list worst divergers
#       (e.g. 6N-SURWIT-F59, 6J-AMLN: drastic gene perturbation, low coherent program).
#   (c) MECHANISM test: does the gene-arm track raw PERTURBATION MAGNITUDE
#       (mean |L2FC| of the signature genes, and fraction of strongly-perturbed genes)
#       rather than DHPS?  If gene-arm ~ magnitude but path-arm ~ DHPS, that is the
#       mechanism: a few violently-perturbed models hijack the gene-arm regardless of
#       whether the perturbation is the COHERENT disease program.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl); library(fgsea) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
OUT   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
strip <- function(x) gsub("\\..*", "", x)
nrm   <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
sp    <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
pr    <- function(a, b) suppressWarnings(cor(a, b, method = "pearson",  use = "complete.obs"))

# ── Published DHPS (MOESM9) ──────────────────────────────────────────────────
m9 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM9_ESM.xlsx"),
                               col_names = FALSE, skip = 2))
setnames(m9, 1:10, c("model","diet_group","mP","mH","m_DHPS","mM","fP","fH","f_DHPS","fM"))
m9 <- m9[!is.na(model)]; for (cc in c("m_DHPS","f_DHPS")) m9[[cc]] <- as.numeric(m9[[cc]])
m9[, `:=`(key = nrm(model), dhps_mean = (as.numeric(m_DHPS)+as.numeric(f_DHPS))/2)]

# ── Vacca S4: 951 progression signature (human) + per-model per-gene L2FC ─────
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"),
                               sheet = "Table S4"))
setnames(s4, c("Gene_used_in_DSEA(0:No/1:Yes)", "ENSGENEID_MOUSE"),
         c("in_dsea", "mouse_id"), skip_absent = TRUE)
s4[, mouse_id := strip(mouse_id)]
l2fc_cols  <- grep("^L2FC_", names(s4), value = TRUE)
human_cols <- grep("UCAM|EPoS", l2fc_cols, value = TRUE)
model_cols <- setdiff(l2fc_cols, human_cols)
for (cc in l2fc_cols) suppressWarnings(s4[, (cc) := as.numeric(get(cc))])
# human progression signature direction = mean(UCAM + EPoS Severe-vs-Mild)
prog_cols <- c("L2FC_UCAM/VCU: Severe vs Mild", "L2FC_EPoS: Severe vs Mild")
s4[, vacca_prog_lfc := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = prog_cols]
vsig <- s4[in_dsea == 1 & !is.na(mouse_id) & mouse_id != "" & is.finite(vacca_prog_lfc)]
vacca_sets <- list(up = unique(vsig[vacca_prog_lfc > 0]$mouse_id),
                   down = unique(vsig[vacca_prog_lfc < 0]$mouse_id))
sig_ids <- unique(c(vacca_sets$up, vacca_sets$down))
cat(sprintf("Vacca 951 progression sig -> mouse orthologs: %d up / %d down (%d total)\n",
            length(vacca_sets$up), length(vacca_sets$down), length(sig_ids)))

dsea_analog <- function(stat, sets) {
  stat <- stat[is.finite(stat)]; stat <- stat[!duplicated(names(stat))]
  if (length(stat) < 200) return(NA_real_)
  fg <- as.data.table(suppressWarnings(fgsea(sets, stat, scoreType = "std", nPermSimple = 5000)))
  nu <- fg[pathway == "up", NES]; nd <- fg[pathway == "down", NES]
  if (length(nu) == 0) nu <- 0; if (length(nd) == 0) nd <- 0
  (nu - nd) / 2
}

# ── Vacca pathway NES: human ref + per-model PEP (MOESM6) ─────────────────────
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"), sheet = "NES"))
setnames(p6, 1, "pathway")
p6_human <- grep("UCAM|EPoS", names(p6), value = TRUE)
p6_model <- setdiff(names(p6)[-1], p6_human)
for (cc in c(p6_human, p6_model)) suppressWarnings(p6[[cc]] <- as.numeric(p6[[cc]]))
href_path <- rowMeans(cbind(p6[["UCAM/VCU: Severe vs Mild"]], p6[["EPoS: Severe vs Mild"]]), na.rm = TRUE)

# ════════════════════════════════════════════════════════════════════════════
# Per-model: GENE-ARM, PATH-ARM, perturbation magnitude
# ════════════════════════════════════════════════════════════════════════════
# index path models by normalized key
p6_model_key <- setNames(p6_model, nrm(p6_model))

res <- rbindlist(lapply(model_cols, function(col) {
  mdl <- sub("^L2FC_", "", col)
  k   <- nrm(mdl)
  # gene-arm
  v  <- s4[is.finite(get(col)) & mouse_id != "", .(mouse_id, lfc = get(col))][!duplicated(mouse_id)]
  st <- setNames(v$lfc, v$mouse_id)
  gene_arm <- dsea_analog(st, vacca_sets)
  # perturbation magnitude over signature genes present in this model
  sigv <- st[names(st) %in% sig_ids]
  mag_mean_abs <- mean(abs(sigv), na.rm = TRUE)             # mean |L2FC| of sig genes
  mag_q90      <- as.numeric(quantile(abs(sigv), 0.90, na.rm = TRUE))
  frac_strong  <- mean(abs(sigv) > 1, na.rm = TRUE)         # fraction sig genes |L2FC|>1
  # GLOBAL magnitude (all genes in the model, not just sig) -> "how drastic overall"
  mag_global   <- mean(abs(st), na.rm = TRUE)
  # path-arm (PEP correlation), only if this model exists in MOESM6
  path_arm <- NA_real_
  if (!is.na(p6_model_key[k])) {
    mv <- p6[[ p6_model_key[k] ]]
    path_arm <- pr(mv, href_path)
  }
  dt <- data.table(model = mdl,
             gene_arm = gene_arm, path_arm = path_arm,
             mag_mean_abs = mag_mean_abs, mag_q90 = mag_q90,
             frac_strong = frac_strong, mag_global = mag_global,
             n_sig = sum(is.finite(sigv)))
  dt[, key := k][]
}))

# ── Merge published DHPS, drop rat (^R-) ─────────────────────────────────────
mg <- merge(res, m9[, .(key, diet_group, m_DHPS, f_DHPS, dhps_mean)], by = "key")
mg <- mg[!grepl("^R[-.]", model)]
mg <- mg[is.finite(gene_arm) & is.finite(path_arm)]
cat(sprintf("\nMouse models with BOTH arms + published DHPS: %d\n", nrow(mg)))

# ── (a) gene-arm vs path-arm correlation ─────────────────────────────────────
ga_pa_sp <- sp(mg$gene_arm, mg$path_arm)
ga_pa_pr <- pr(mg$gene_arm, mg$path_arm)
g_dhps_sp <- sp(mg$gene_arm, mg$dhps_mean)
p_dhps_sp <- sp(mg$path_arm, mg$dhps_mean)
cat(sprintf("\n(a) GENE-ARM vs PATH-ARM across models:  Spearman %.3f  Pearson %.3f\n",
            ga_pa_sp, ga_pa_pr))
cat(sprintf("    GENE-ARM vs published DHPS_mean:    Spearman %.3f\n", g_dhps_sp))
cat(sprintf("    PATH-ARM vs published DHPS_mean:    Spearman %.3f\n", p_dhps_sp))

# ── (b) per-model divergence = rank(path) - rank(gene) ───────────────────────
mg[, rank_gene := frank(-gene_arm, ties.method = "average")]
mg[, rank_path := frank(-path_arm, ties.method = "average")]
mg[, rank_dhps := frank(-dhps_mean, ties.method = "average")]
mg[, divergence := rank_path - rank_gene]   # +ve = path ranks it LOW, gene ranks it HIGH (gene-inflated)
setorder(mg, divergence)
cat("\n(b) WORST DIVERGERS (gene-arm inflates vs pathway-arm; rank_path - rank_gene most negative):\n")
print(mg[1:6, .(model, diet_group, gene_arm = round(gene_arm,2), path_arm = round(path_arm,2),
                dhps = round(dhps_mean,2), rank_gene, rank_path, rank_dhps,
                mag_mean_abs = round(mag_mean_abs,2), frac_strong = round(frac_strong,2))])
cat("\n   WORST DIVERGERS (other direction; pathway-high/gene-low):\n")
print(mg[(.N-2):.N, .(model, diet_group, gene_arm = round(gene_arm,2), path_arm = round(path_arm,2),
                dhps = round(dhps_mean,2), rank_gene, rank_path, rank_dhps,
                mag_mean_abs = round(mag_mean_abs,2), frac_strong = round(frac_strong,2))])

# ── (c) MECHANISM: gene-arm tracks PERTURBATION MAGNITUDE, not DHPS ──────────
cat("\n(c) MECHANISM — does each arm track raw gene-perturbation magnitude vs DHPS?\n")
mech <- data.table(
  predictor = c("mag_mean_abs (sig genes mean|L2FC|)", "mag_q90 (sig genes |L2FC| q90)",
                "frac_strong (sig genes |L2FC|>1)", "mag_global (all genes mean|L2FC|)"),
  cor_with_GENEarm = c(sp(mg$gene_arm, mg$mag_mean_abs), sp(mg$gene_arm, mg$mag_q90),
                       sp(mg$gene_arm, mg$frac_strong),  sp(mg$gene_arm, mg$mag_global)),
  cor_with_PATHarm = c(sp(mg$path_arm, mg$mag_mean_abs), sp(mg$path_arm, mg$mag_q90),
                       sp(mg$path_arm, mg$frac_strong),  sp(mg$path_arm, mg$mag_global)),
  cor_with_DHPS    = c(sp(mg$dhps_mean, mg$mag_mean_abs), sp(mg$dhps_mean, mg$mag_q90),
                       sp(mg$dhps_mean, mg$frac_strong),  sp(mg$dhps_mean, mg$mag_global)))
print(mech)

# headline numbers
hdr_gene_mag <- sp(mg$gene_arm, mg$mag_mean_abs)
hdr_path_mag <- sp(mg$path_arm, mg$mag_mean_abs)
cat(sprintf("\n   GENE-ARM tracks signature perturbation magnitude at Spearman %.3f\n", hdr_gene_mag))
cat(sprintf("   PATH-ARM tracks signature perturbation magnitude at Spearman %.3f\n", hdr_path_mag))
cat(sprintf("   (DHPS tracks magnitude at %.3f -- magnitude is NOT DHPS)\n",
            sp(mg$dhps_mean, mg$mag_mean_abs)))

# ── write outputs ─────────────────────────────────────────────────────────────
setorder(mg, -dhps_mean)
fwrite(mg, file.path(OUT, "q1_gene_vs_pathway_arm_decomposition.csv"))

summ <- data.table(
  metric = c("n_models",
             "geneArm_vs_pathArm_spearman", "geneArm_vs_pathArm_pearson",
             "geneArm_vs_DHPS_spearman", "pathArm_vs_DHPS_spearman",
             "geneArm_vs_magnitude_spearman", "pathArm_vs_magnitude_spearman",
             "DHPS_vs_magnitude_spearman",
             "worst_diverger_model", "worst_diverger_geneArm", "worst_diverger_pathArm",
             "worst_diverger_DHPS", "worst_diverger_mag_mean_abs"),
  value  = c(nrow(mg),
             round(ga_pa_sp,3), round(ga_pa_pr,3),
             round(g_dhps_sp,3), round(p_dhps_sp,3),
             round(hdr_gene_mag,3), round(hdr_path_mag,3),
             round(sp(mg$dhps_mean, mg$mag_mean_abs),3),
             mg[which.min(divergence)]$model,
             round(mg[which.min(divergence)]$gene_arm,3),
             round(mg[which.min(divergence)]$path_arm,3),
             round(mg[which.min(divergence)]$dhps_mean,3),
             round(mg[which.min(divergence)]$mag_mean_abs,3)))
fwrite(summ, file.path(OUT, "q1_gene_vs_pathway_arm_summary.csv"))
cat("\nWrote q1_gene_vs_pathway_arm_decomposition.csv + q1_gene_vs_pathway_arm_summary.csv\n")
cat(strrep("=", 78), "\n")
cat("VERDICT: if geneArm~magnitude high & pathArm~magnitude low & geneArm~DHPS<<pathArm~DHPS\n")
cat("  => gene-arm is hijacked by a few violently-perturbed models (drastic but INCOHERENT\n")
cat("     perturbation); pathway-arm averages over a KEGG program => robust to single-gene\n")
cat("     blow-ups => reproduces DHPS. THIS is why 0.06 (gene) vs 0.73 (pathway).\n")
cat(strrep("=", 78), "\n")
