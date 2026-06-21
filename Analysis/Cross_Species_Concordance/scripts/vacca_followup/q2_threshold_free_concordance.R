#!/usr/bin/env Rscript
# ============================================================================
# q2_threshold_free_concordance.R  (Vacca followup Q2)
#
# The thresholded "Conserved_Core" overlap with Vacca's 951-gene signature
# (172/1323; Fisher OR 2.80x, p=2.0e-24) is bottlenecked by ARBITRARY human+
# mouse LFC/padj cutoffs on both sides. Here we build a THRESHOLD-FREE,
# continuous per-gene cross-species concordance score (NO hard DEG cutoffs):
#
#   For each gene with a 1:1-mappable human↔mouse ortholog:
#     h_t  = human signed significance (canonical_deg_results.csv `t`)
#     m_t  = mouse signed significance, mean over the 4 per_diet DE `t`
#     S_geom = mean_over_diets[ sign(h_t)*sign(m_t)*sqrt(|h_t|*|m_t|) ]   (primary)
#     S_rp   = rank-product of human & mouse signed -log10(padj)  (alt ranking)
#
# Then, against Vacca's 951 signature (mapped to our human genes):
#   (a) fgsea of the Vacca-951 set in the continuous concordance ranking
#       (ES/NES/padj) — does the concordance ranking ENRICH for Vacca genes?
#   (b) Spearman of our per-gene concordance score vs Vacca progression logFC
#       (mean of UCAM/VCU + EPoS Severe-vs-Mild L2FC) over shared genes.
#
# Headline comparison: is the continuous (threshold-free) signal STRONGER than
# the thresholded core's Fisher OR 2.80x? fgsea NES/padj and the Spearman give
# a cutoff-independent answer.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl); library(fgsea) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
ANNOT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation")
INT_I <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
MPD   <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
OUT   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
strip <- function(x) gsub("\\..*", "", x)
sp <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))

# ── Ortholog map (human↔mouse, 1:1 for clean mapping) ────────────────────────
ortho <- fread(file.path(ANNOT, "ortholog_mapping.tsv"))
ortho[, `:=`(hbase = strip(human_gene_id), mbase = strip(mouse_gene_id))]
# keep 1:1 to avoid one-to-many duplication of the continuous score
o1 <- ortho[orthology_type == "ortholog_one2one", .(hbase, mbase, human_symbol)]
o1 <- o1[!duplicated(hbase) & !duplicated(mbase)]
cat(sprintf("1:1 orthologs: %d\n", nrow(o1)))

# ── Human signed significance (t) and signed -log10(padj) ────────────────────
hum <- fread(file.path(INT_I, "canonical_deg_results.csv"))
hum[, hbase := strip(gene)]
hum <- hum[is.finite(t) & !is.na(padj)]
hum[, h_signlp := sign(t) * -log10(pmax(padj, 1e-300))]
hum <- hum[, .(hbase, h_t = t, h_signlp, h_symbol = symbol)][!duplicated(hbase)]

# ── Mouse signed significance: mean t and mean signed -log10(padj) over 4 diets
diets <- c("MCD", "HFD", "CDAHFD", "FPC")
mlist <- lapply(diets, function(d) {
  md <- fread(file.path(MPD, paste0(d, "_de_results.csv")))
  md[, mbase := strip(gene)]
  md <- md[is.finite(t) & !is.na(adj.P.Val)]
  md[, m_signlp := sign(t) * -log10(pmax(adj.P.Val, 1e-300))]
  md <- md[, .(mbase, t, m_signlp)][!duplicated(mbase)]
  setnames(md, c("t", "m_signlp"), c(paste0("t_", d), paste0("slp_", d)))
  md
})
mouse <- Reduce(function(a, b) merge(a, b, by = "mbase", all = TRUE), mlist)
t_cols   <- paste0("t_",   diets)
slp_cols <- paste0("slp_", diets)
mouse[, m_t      := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = t_cols]
mouse[, m_signlp := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = slp_cols]

# ── Per-diet geometric concordance, then mean over diets (PRIMARY score) ─────
# S_geom_d = sign(h_t)*sign(m_t_d)*sqrt(|h_t|*|m_t_d|); averaged over available diets.
g <- merge(o1, hum, by = "hbase")
g <- merge(g, mouse, by = "mbase")
cat(sprintf("Genes with human + mouse(>=1 diet) + 1:1 ortholog: %d\n", nrow(g)))
for (d in diets) {
  g[[paste0("Sgeom_", d)]] <- sign(g$h_t) * sign(g[[paste0("t_", d)]]) *
    sqrt(abs(g$h_t) * abs(g[[paste0("t_", d)]]))
}
sg_cols <- paste0("Sgeom_", diets)
g[, S_geom := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = sg_cols]

# ── Rank-product of signed -log10(padj) (ALT score, fully threshold-free) ────
# Concordant-up genes get high (h_signlp & m_signlp both high +); concordant-down
# get low (both very negative). Use a signed product of within-species signed ranks.
g[, h_rank := frank(h_signlp, ties.method = "average")]
g[, m_rank := frank(m_signlp, ties.method = "average")]
N <- nrow(g)
# center ranks to [-1,1] so sign carries direction; product rewards agreement
g[, h_rc := 2 * (h_rank - 0.5) / N - 1]
g[, m_rc := 2 * (m_rank - 0.5) / N - 1]
g[, S_rankprod := h_rc * m_rc * sign(h_rc * m_rc) * sqrt(abs(h_rc * m_rc))]
# simpler signed agreement product also kept for the fgsea ranking
g[, S_signlp := sign(h_signlp) * sign(m_signlp) * sqrt(abs(h_signlp) * abs(m_signlp))]

# ── Vacca 951 signature → our human genes (via human ENSEMBL id, strip ver) ──
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"),
                               sheet = "Table S4"))
setnames(s4, "Gene_used_in_DSEA(0:No/1:Yes)", "in_dsea", skip_absent = TRUE)
s4[, hENS := strip(ENSEMBL_GENE_ID_HUMAN)]
prog_cols <- c("L2FC_UCAM/VCU: Severe vs Mild", "L2FC_EPoS: Severe vs Mild")
for (cc in prog_cols) suppressWarnings(s4[, (cc) := as.numeric(get(cc))])
s4[, vacca_prog_lfc := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = prog_cols]
vacca951_ens <- unique(s4[in_dsea == 1 & !is.na(hENS) & hENS != "", hENS])
cat(sprintf("Vacca-951 mapped to human ENSEMBL ids: %d\n", length(vacca951_ens)))

# Vacca-951 set in OUR gene space (genes present in our threshold-free ranking)
vset <- intersect(vacca951_ens, g$hbase)
cat(sprintf("Vacca-951 ∩ our ranked genes: %d\n", length(vset)))

# ════════════════════════════════════════════════════════════════════════════
# (a) fgsea: Vacca-951 enrichment in each continuous concordance ranking
# ════════════════════════════════════════════════════════════════════════════
run_fgsea <- function(stat_col) {
  v <- g[is.finite(get(stat_col)), .(hbase, s = get(stat_col))][!duplicated(hbase)]
  st <- setNames(v$s, v$hbase)
  st <- st[order(-st)]
  fg <- as.data.table(suppressWarnings(fgsea(list(vacca951 = vset), st,
                                             scoreType = "std", nPermSimple = 20000, eps = 0)))
  data.table(ranking = stat_col, n_ranked = length(st),
             ES = fg$ES, NES = fg$NES, pval = fg$pval, padj = fg$padj,
             size = fg$size, leadingEdge_n = length(fg$leadingEdge[[1]]))
}
cat("\n=== (a) fgsea of Vacca-951 in threshold-free concordance rankings ===\n")
fa <- rbindlist(lapply(c("S_geom", "S_signlp", "S_rankprod"), run_fgsea))
print(fa)

# control: human-only ranking (signed -log10 padj) — to show cross-species adds signal
hv <- g[is.finite(h_signlp), .(hbase, s = h_signlp)][!duplicated(hbase)]
hst <- setNames(hv$s, hv$hbase); hst <- hst[order(-hst)]
fh <- as.data.table(suppressWarnings(fgsea(list(vacca951 = vset), hst,
                                           scoreType = "std", nPermSimple = 20000, eps = 0)))
cat(sprintf("  [control] human-only signed -log10(padj): NES %.2f padj %.2e (size %d)\n",
            fh$NES, fh$padj, fh$size))

# ════════════════════════════════════════════════════════════════════════════
# (b) Spearman: our continuous concordance vs Vacca progression logFC
# ════════════════════════════════════════════════════════════════════════════
cat("\n=== (b) Spearman: continuous concordance vs Vacca progression logFC ===\n")
vlfc <- s4[!is.na(hENS) & hENS != "" & is.finite(vacca_prog_lfc),
           .(hbase = hENS, vacca_prog_lfc)][!duplicated(hbase)]
gv <- merge(g, vlfc, by = "hbase")
cat(sprintf("  genes with both our score and Vacca progression logFC: %d\n", nrow(gv)))
spear <- data.table(
  comparison = c("S_geom_vs_vaccaProgLFC", "S_signlp_vs_vaccaProgLFC",
                 "S_rankprod_vs_vaccaProgLFC", "h_t_only_vs_vaccaProgLFC",
                 "m_t_only_vs_vaccaProgLFC"),
  spearman = c(sp(gv$S_geom, gv$vacca_prog_lfc), sp(gv$S_signlp, gv$vacca_prog_lfc),
               sp(gv$S_rankprod, gv$vacca_prog_lfc), sp(gv$h_t, gv$vacca_prog_lfc),
               sp(gv$m_t, gv$vacca_prog_lfc)),
  n = nrow(gv))
print(spear)

# Restricted to the Vacca-951 signature genes (the published signature itself)
gv951 <- gv[hbase %in% vset]
cat(sprintf("\n  within Vacca-951 (n=%d): S_geom vs progLFC Spearman = %.3f ; h_t = %.3f ; m_t = %.3f\n",
            nrow(gv951), sp(gv951$S_geom, gv951$vacca_prog_lfc),
            sp(gv951$h_t, gv951$vacca_prog_lfc), sp(gv951$m_t, gv951$vacca_prog_lfc)))

# ════════════════════════════════════════════════════════════════════════════
# Threshold-FREE effect-size comparison to the thresholded OR 2.80x
# AUROC: can the continuous score rank Vacca-951 genes above the rest? (no cutoff)
# ════════════════════════════════════════════════════════════════════════════
auroc <- function(score, is_pos) {
  ok <- is.finite(score); score <- score[ok]; is_pos <- is_pos[ok]
  r <- rank(score); n1 <- sum(is_pos); n0 <- sum(!is_pos)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[is_pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
g[, is_vacca := hbase %in% vset]
cat("\n=== threshold-free AUROC: concordance score ranks Vacca-951 above rest ===\n")
auc_tab <- data.table(
  ranking  = c("S_geom", "S_signlp", "S_rankprod", "h_signlp(human-only)", "m_signlp(mouse-only)"),
  AUROC    = c(auroc(g$S_geom, g$is_vacca), auroc(g$S_signlp, g$is_vacca),
               auroc(g$S_rankprod, g$is_vacca), auroc(g$h_signlp, g$is_vacca),
               auroc(g$m_signlp, g$is_vacca)),
  n_pos = sum(g$is_vacca), n_tot = nrow(g))
# AUROC uses |score| direction-agnostically too? Vacca-951 is a mix of up/down →
# use absolute concordance magnitude for a direction-free "is it concordant at all"
g[, S_geom_abs := abs(S_geom)]
auc_tab_abs <- data.table(ranking = "S_geom_abs(|concordance|)",
  AUROC = auroc(g$S_geom_abs, g$is_vacca), n_pos = sum(g$is_vacca), n_tot = nrow(g))
print(rbind(auc_tab, auc_tab_abs))

# ── Write outputs ────────────────────────────────────────────────────────────
fwrite(g[, .(hbase, mbase, h_symbol, h_t, m_t, h_signlp, m_signlp,
             S_geom, S_signlp, S_rankprod, is_vacca)],
       file.path(OUT, "q2_threshold_free_concordance.csv"))
summary_dt <- rbindlist(list(
  fa[, .(metric = paste0("fgsea_NES_", ranking), value = NES, aux = padj)],
  spear[, .(metric = paste0("spearman_", comparison), value = spearman, aux = n)],
  rbind(auc_tab, auc_tab_abs)[, .(metric = paste0("AUROC_", ranking), value = AUROC, aux = n_pos)]
), use.names = TRUE)
summary_dt <- rbind(summary_dt,
  data.table(metric = "thresholded_core_Fisher_OR", value = 2.7956, aux = 172),
  data.table(metric = "control_humanOnly_fgsea_NES", value = fh$NES, aux = fh$padj))
fwrite(summary_dt, file.path(OUT, "q2_threshold_free_summary.csv"))

cat("\n", strrep("=", 78), "\n")
cat("VERDICT: compare fgsea NES/padj (continuous) + AUROC vs thresholded Fisher OR 2.80x.\n")
cat("  Strong NES (>2, padj<<0.05) and AUROC>>0.5 = threshold-free signal is robust,\n")
cat("  not an artifact of the LFC/padj cutoffs; magnitude (NES, AUROC) is the\n")
cat("  cutoff-independent analogue of the OR.\n")
cat(strrep("=", 78), "\n")
