#!/usr/bin/env Rscript
# ============================================================================
# 08_vacca_reproduction_diagnostic.R
#
# WHY: our cross-species per-diet ranking ANTI-correlates with Vacca 2024's DHPS
# (Spearman -1). Before publishing any model-ranking we must know WHY. This
# diagnostic localizes the cause with a controlled fgsea DSEA-analog ladder,
# reusing OUR fgsea + VACCA's own provided per-model L2FC (MOESM4) + published
# DHPS (MOESM9). (Vacca themselves used ctlab/fgsea.)
#
#   DSEA-analog per model/diet = (NES_up - NES_down)/2  [up genes up AND down
#   genes down → high; matches Vacca's "downregulated negated, averaged"].
#
#   D2: OUR disease-vs-ctrl signature, ranked by VACCA per-model L2FC.
#   D3: VACCA 951 progression signature, ranked by VACCA per-model L2FC.  ← KEY
#       (does fgsea + their signature reproduce their published DHPS on THEIR data?)
#   D0: OUR signature, OUR per-diet DE (the current, inverting baseline).
#   D1: VACCA signature, OUR per-diet DE (does their signature fix our ranking?).
#   D1b: OUR progression signature (nafl_vs_nash_lvqw), OUR per-diet DE.
#
# Decision: D3 Spearman vs published DHPS >= ~0.6 ⇒ method adequate (inversion is
# signature/data); D3 < ~0.6 ⇒ the DSEA pathway/[0,1] machinery is essential.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(readxl); library(fgsea)
})
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
ANNOT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation")
INT_I <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
DSIG  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures")
MPD   <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
OUT   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
strip <- function(x) gsub("\\..*", "", x)
nrm   <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

# DSEA-analog via fgsea: (NES_up - NES_down)/2
dsea_analog <- function(stat, sets) {
  stat <- stat[is.finite(stat)]; stat <- stat[!duplicated(names(stat))]
  if (length(stat) < 200) return(NA_real_)
  fg <- as.data.table(suppressWarnings(fgsea(sets, stat, scoreType = "std", nPermSimple = 5000)))
  nu <- fg[pathway == "up", NES]; nd <- fg[pathway == "down", NES]
  if (length(nu) == 0) nu <- 0; if (length(nd) == 0) nd <- 0
  (nu - nd) / 2
}

# ── Ortholog map (human↔mouse) ───────────────────────────────────────────────
ortho <- fread(file.path(ANNOT, "ortholog_mapping.tsv"))
ortho[, `:=`(hbase = strip(human_gene_id), mbase = strip(mouse_gene_id))]

# ── Vacca S4: signature + per-model L2FC ─────────────────────────────────────
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"),
                               sheet = "Table S4"))
setnames(s4, c("Gene_used_in_DSEA(0:No/1:Yes)", "ENSGENEID_MOUSE"),
         c("in_dsea", "mouse_id"), skip_absent = TRUE)
s4[, mouse_id := strip(mouse_id)]
l2fc_cols  <- grep("^L2FC_", names(s4), value = TRUE)
human_cols <- grep("UCAM|EPoS", l2fc_cols, value = TRUE)
model_cols <- setdiff(l2fc_cols, human_cols)                       # ~36 mouse models
for (cc in l2fc_cols) suppressWarnings(s4[, (cc) := as.numeric(get(cc))])
# Vacca human PROGRESSION logFC = mean(UCAM/VCU + EPoS Severe-vs-Mild)
prog_cols <- c("L2FC_UCAM/VCU: Severe vs Mild", "L2FC_EPoS: Severe vs Mild")
s4[, vacca_prog_lfc := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = prog_cols]
vsig <- s4[in_dsea == 1 & !is.na(mouse_id) & mouse_id != "" & is.finite(vacca_prog_lfc)]
vacca_sets <- list(up = unique(vsig[vacca_prog_lfc > 0]$mouse_id),
                   down = unique(vsig[vacca_prog_lfc < 0]$mouse_id))
cat(sprintf("Vacca progression signature: %d up / %d down (mouse orthologs)\n",
            length(vacca_sets$up), length(vacca_sets$down)))

# ── OUR disease-vs-control signature → mouse sets ────────────────────────────
hum <- fread(file.path(INT_I, "canonical_deg_results.csv"))
hum[, hbase := strip(gene)]
hm  <- merge(hum[, .(hbase, h_lfc = logFC, h_padj = padj)], ortho[, .(hbase, mbase)], by = "hbase")
our_sets <- list(up = unique(hm[h_padj < 0.05 & h_lfc > 0.5, mbase]),
                 down = unique(hm[h_padj < 0.05 & h_lfc < -0.5, mbase]))

# ── OUR progression signature (nafl_vs_nash LVQW) → mouse sets ───────────────
nn <- fread(file.path(DSIG, "nafl_vs_nash_lvqw.csv"))
nn[, hbase := strip(gene)]
nm <- merge(nn[, .(hbase, h_lfc = logFC, h_padj = padj)], ortho[, .(hbase, mbase)], by = "hbase")
ourprog_sets <- list(up = unique(nm[h_padj < 0.05 & h_lfc > 0.5, mbase]),
                     down = unique(nm[h_padj < 0.05 & h_lfc < -0.5, mbase]))
cat(sprintf("Our disease-vs-ctrl sig: %d up / %d down ; our progression sig: %d up / %d down\n",
            length(our_sets$up), length(our_sets$down),
            length(ourprog_sets$up), length(ourprog_sets$down)))

# ── Published DHPS (MOESM9) ──────────────────────────────────────────────────
m9 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM9_ESM.xlsx"),
                               col_names = FALSE, skip = 2))
setnames(m9, 1:10, c("model","diet_group","mP","mH","m_DHPS","mM","fP","fH","f_DHPS","fM"))
m9 <- m9[!is.na(model)]; for (cc in c("m_DHPS","f_DHPS")) m9[[cc]] <- as.numeric(m9[[cc]])
m9[, key := nrm(model)]

# ════════════════════════════════════════════════════════════════════════════
# D2 / D3 : reproduce on VACCA's per-model L2FC
# ════════════════════════════════════════════════════════════════════════════
cat("\n=== D2/D3: our fgsea on VACCA per-model L2FC vs published DHPS ===\n")
res_models <- rbindlist(lapply(model_cols, function(col) {
  v <- s4[is.finite(get(col)) & mouse_id != "", .(mouse_id, lfc = get(col))][!duplicated(mouse_id)]
  st <- setNames(v$lfc, v$mouse_id)
  data.table(model = sub("^L2FC_", "", col),
             d2_our_sig   = dsea_analog(st, our_sets),
             d3_vacca_sig = dsea_analog(st, vacca_sets),
             n_genes = length(st))
}))
res_models[, key := nrm(model)]
mg <- merge(res_models, m9[, .(key, diet_group, m_DHPS, f_DHPS)], by = "key")
mg[, dhps_mean := (m_DHPS + f_DHPS) / 2]
sp <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
cat(sprintf("  models matched to published DHPS: %d\n", nrow(mg)))
cat(sprintf("  D2 (OUR sig)   vs published DHPS: Spearman metab %.2f / fibro %.2f / mean %.2f\n",
            sp(mg$d2_our_sig, mg$m_DHPS), sp(mg$d2_our_sig, mg$f_DHPS), sp(mg$d2_our_sig, mg$dhps_mean)))
cat(sprintf("  D3 (VACCA sig) vs published DHPS: Spearman metab %.2f / fibro %.2f / mean %.2f   <-- KEY\n",
            sp(mg$d3_vacca_sig, mg$m_DHPS), sp(mg$d3_vacca_sig, mg$f_DHPS), sp(mg$d3_vacca_sig, mg$dhps_mean)))
fwrite(mg[order(-dhps_mean)], file.path(OUT, "diagnostic_D2D3_per_model.csv"))
cat("\n  per-model (top/bottom by published DHPS):\n")
print(mg[order(-dhps_mean), .(model, diet_group, d2_our_sig = round(d2_our_sig,2),
        d3_vacca_sig = round(d3_vacca_sig,2), m_DHPS = round(m_DHPS,2), f_DHPS = round(f_DHPS,2))][c(1:5, (.N-4):.N)])

# ════════════════════════════════════════════════════════════════════════════
# D0 / D1 / D1b : our 4 diets (rank by t-stat)
# ════════════════════════════════════════════════════════════════════════════
cat("\n=== D0/D1: our 4 diets — DSEA-analog under each signature ===\n")
res_diets <- rbindlist(lapply(c("MCD","HFD","CDAHFD","FPC"), function(d) {
  md <- fread(file.path(MPD, paste0(d, "_de_results.csv")))
  md[, mbase := strip(gene)]; md <- md[is.finite(t)]
  st <- setNames(md$t, md$mbase)
  data.table(our_diet = d,
             D0_our_dvc   = dsea_analog(st, our_sets),
             D1_vacca_sig = dsea_analog(st, vacca_sets),
             D1b_our_prog = dsea_analog(st, ourprog_sets))
}))
# Vacca same-diet published DHPS (MCD←CDD, CDAHFD←CDHFD, HFD←HFD, FPC←WD*)
vgrp <- function(d) switch(d, MCD = m9[diet_group == "CDD" & !grepl("^R-", model)],
  CDAHFD = m9[grepl("^CDHFD$", diet_group)], HFD = m9[diet_group == "HFD"],
  FPC = m9[grepl("^WD", diet_group)])
res_diets[, vacca_DHPS_metab := sapply(our_diet, function(d) mean(vgrp(d)$m_DHPS, na.rm = TRUE))]
print(res_diets)
cat(sprintf("\n  Across our 4 diets — Spearman vs Vacca same-diet metabolic DHPS:\n    D0 (our dvc)=%.2f  D1 (vacca sig)=%.2f  D1b (our prog)=%.2f\n",
            sp(res_diets$D0_our_dvc, res_diets$vacca_DHPS_metab),
            sp(res_diets$D1_vacca_sig, res_diets$vacca_DHPS_metab),
            sp(res_diets$D1b_our_prog, res_diets$vacca_DHPS_metab)))
fwrite(res_diets, file.path(OUT, "diagnostic_D0D1_our_diets.csv"))

cat("\n", strrep("=", 78), "\nVERDICT GUIDE:\n",
    "  D3 high (>=0.6): fgsea+their signature reproduces their DHPS on their data\n",
    "    → method adequate; inversion = signature (compare D0 vs D1) and/or our mouse data (D2 vs D3).\n",
    "  D3 low: DSEA pathway/[0,1] machinery essential → reimplement DSEA or ADOPT their DHPS.\n",
    strrep("=", 78), "\n", sep = "")
