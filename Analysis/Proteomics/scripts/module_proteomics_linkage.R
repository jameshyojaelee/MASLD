#!/usr/bin/env Rscript
# module_proteomics_linkage.R
# ─────────────────────────────────────────────────────────────────────────────
# Cross-modal, stage-resolved linkage of the Fig3 Hotspot single-cell co-expression
# modules to the Fig4 liver proteomics (PXD051911 DIA-MS, n=58, documented histology).
# Produces the tables consumed by scripts/figures/crossmodal_module_proteomics.R:
#   (A) module_proteomics_enrichment.csv    — fgsea NES of each disease-significant
#       Hotspot module gene set over the protein moderated-t rank (PXD051911 liver;
#       PXD052937 plasma as a flagged, compartment-mismatched replication).
#   (B) module_protein_stage_slopes.csv     — per-patient module protein score
#       (mean of per-protein z) regressed on the documented liver stage axes
#       (saf_diagnosis ordinal 0/1/2, NAS 0-8, Kleiner fibrosis).
#   (B')crossmodal_module_concordance.csv   — beta_scRNA (fig3g 3-stage slope) vs
#       beta_protein (saf ordinal) per module; Spearman rho + sign concordance.
#   (C) module_cascade_protein.csv          — per-module mean protein score across
#       No_MASLD -> MASL -> MASH (the fig3e handoff, in liver protein).
#
# MODALITY DOCTRINE: the Hotspot module GENE LISTS are used as a-priori sets only.
# Aggregation + statistics are proteomics-native (log2 -> per-protein z across
# patients -> mean-z module score; fgsea over the moderated-t rank). The Hotspot
# autocorrelation algorithm itself is NOT ported to the protein modality.
#
# STAGE AXIS: saf_diagnosis (No_MASLD/MASL/MASH) maps 1:1 onto fig3g's coarse
# Healthy/Steatosis/Steatohepatitis, so beta_scRNA and beta_protein share one
# ordinal. Kleiner fibrosis in this bariatric cohort tops out at F3 (n=3, no F4),
# so it is reported but not used as the primary axis.
#
# COMPOSITION CAVEAT: beta_scRNA is within-cell-type (composition-controlled);
# a bulk-liver-protein module score is not. Concordance/discordance is therefore an
# interpretive axis (per-cell regulation vs composition), not a pure replication.
# No bulk-protein deconvolution is attempted (stated limitation).
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(data.table); library(fgsea) })
set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HS   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules")
PROT <- file.path(BASE, "Analysis/Proteomics/results")
FIG3 <- file.path(BASE, "figures/main/fig3_RNAseq/panels/data")
OUT  <- PROT
MIN_DETECTED <- 8L    # module must have >= this many member proteins detected
CELL_TYPES <- c("hepatocytes", "macrophages", "fibroblasts", "cholangiocytes", "tcells")

# ── 0. ENSG -> HGNC symbol map (module_genes.tsv mixes symbols and ENSG ids) ───
gmeta   <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
ens2sym <- gmeta[, setNames(gene_name, ensembl_base)]
to_symbol <- function(g) {
  out <- g
  is_ens <- grepl("^ENSG", g)
  out[is_ens] <- unname(ens2sym[sub("\\..*", "", g[is_ens])])
  out
}

# ── 1. disease-significant Hotspot modules + their gene sets (symbols) ─────────
catmod <- fread(file.path(HS, "all_modules.tsv"))
sig <- catmod[is.finite(disease_stage_q) & disease_stage_q < 0.05,
              .(cell_type, module, module_name, disease_stage_beta, disease_stage_q,
                F_stage_beta, F_stage_q, progression_module, activity_module)]
sig[, mid := paste0(cell_type, "__", module)]
cat(sprintf("[1] disease-significant modules: %d across %d cell types\n",
            nrow(sig), uniqueN(sig$cell_type)))

gene_sets <- list()
for (ct in CELL_TYPES) {
  mg <- fread(file.path(HS, ct, "module_genes.tsv"))
  mg[, sym := to_symbol(gene)]
  mg <- mg[!is.na(sym) & sym != ""]
  for (m in unique(mg$module)) {
    mid <- paste0(ct, "__", m)
    if (mid %in% sig$mid) gene_sets[[mid]] <- unique(mg[module == m, sym])
  }
}
cat(sprintf("[1] built %d module gene sets (median %d genes)\n",
            length(gene_sets), as.integer(median(lengths(gene_sets)))))

# ── 2. protein DE ranks (moderated t) ─────────────────────────────────────────
de_all <- fread(file.path(PROT, "protein_differential_results_v3.csv"))
rank_of <- function(ds) {
  d <- de_all[dataset == ds & is.finite(t) & gene != ""][order(-abs(t))][!duplicated(gene)]
  sort(setNames(d$t, d$gene), decreasing = TRUE)
}

# ── 3. (A) module enrichment in proteomics ────────────────────────────────────
run_enr <- function(ds, cohort, tissue) {
  rk   <- rank_of(ds)
  sets <- lapply(gene_sets, function(s) intersect(s, names(rk)))
  ndet <- lengths(sets)
  keep <- names(sets)[ndet >= MIN_DETECTED]
  cat(sprintf("[A] %-16s %d modules tested, %d dropped (<%d proteins detected)\n",
              cohort, length(keep), sum(ndet < MIN_DETECTED), MIN_DETECTED))
  fg <- as.data.table(fgsea(pathways = sets[keep], stats = rk,
                            minSize = MIN_DETECTED, maxSize = 2000, eps = 0))
  fg <- fg[, .(mid = pathway, NES = round(NES, 3), pval = signif(pval, 4),
               padj = signif(padj, 4), n_detected = size)]
  fg[, `:=`(cohort = cohort, tissue = tissue)]
  fg
}
enrA <- rbindlist(list(
  run_enr("PXD051911", "PXD051911_liver",  "liver"),
  run_enr("PXD052937", "PXD052937_plasma", "plasma")   # compartment-mismatched: replication only
), fill = TRUE)
enrA <- merge(enrA, sig[, .(mid, cell_type, module, module_name,
                            scrna_beta = disease_stage_beta)], by = "mid")
enrA[, scrna_dir := ifelse(scrna_beta > 0, "up", "down")]
setorder(enrA, cohort, -NES)
fwrite(enrA, file.path(OUT, "module_proteomics_enrichment.csv"))

# ── 4. liver abundance matrix -> per-protein z across patients (mean-z scoring) ─
ab <- fread(file.path(BASE, "data/PXD051911/liver_protein_quant.txt"))
ab <- ab[Genes != "" & !is.na(Genes)]
ab[, gene := tstrsplit(Genes, ";", fixed = TRUE)[[1]]]
ab <- ab[!duplicated(gene)]
scols <- setdiff(names(ab), c("ProteinAccessions", "Genes", "ProteinDescriptions", "gene"))
abmat <- as.matrix(ab[, ..scols]); rownames(abmat) <- ab$gene; mode(abmat) <- "numeric"
labmat <- log2(abmat)                                    # raw intensities -> log2
rowZ <- function(M) {                                    # z per protein (row), NA-safe
  mu <- rowMeans(M, na.rm = TRUE)
  sdv <- apply(M, 1, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- NA
  (M - mu) / sdv
}
zmat <- rowZ(labmat)                                     # protein x patient, NA preserved

# ── 5. patient histology (documented, one baseline liver sample per patient) ──
mh <- fread(file.path(BASE, "data/PXD051911/meta_data.txt"))
mh <- mh[liver_proteomics_filename != "" & !is.na(liver_proteomics_filename)]
mh[, `:=`(raw = liver_proteomics_filename,
          saf_ord  = c(No_MASLD = 0, MASL = 1, MASH = 2)[saf_diagnosis],
          NAS      = as.numeric(nafld_activity_score),
          Fibrosis = as.numeric(sub("F", "", kleiner_fibrosis_grade)))]
mh <- mh[raw %in% colnames(zmat)]
zmat <- zmat[, mh$raw, drop = FALSE]                     # align columns to patients
cat(sprintf("[5] liver patients matched: %d (saf: %s)\n", nrow(mh),
            paste(sprintf("%s=%d", names(table(mh$saf_diagnosis)), table(mh$saf_diagnosis)),
                  collapse = " ")))

module_score <- function(mid) {                          # mean-z over detected members
  syms <- intersect(gene_sets[[mid]], rownames(zmat))
  if (length(syms) < MIN_DETECTED) return(NULL)
  list(score = colMeans(zmat[syms, , drop = FALSE], na.rm = TRUE), n = length(syms))
}
slope_on <- function(y, x) {                             # tidy lm(y ~ x)
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) < 5 || length(unique(x[ok])) < 2) return(c(beta = NA, se = NA, p = NA))
  s <- summary(lm(y[ok] ~ x[ok]))$coefficients
  c(beta = s[2, 1], se = s[2, 2], p = s[2, 4])
}

# ── 6. (B) per-module protein score + stage slopes ────────────────────────────
slopes <- rbindlist(lapply(names(gene_sets), function(mid) {
  ms <- module_score(mid); if (is.null(ms)) return(NULL)
  y <- ms$score
  b_saf <- slope_on(y, mh$saf_ord); b_nas <- slope_on(y, mh$NAS); b_fib <- slope_on(y, mh$Fibrosis)
  data.table(mid = mid, n_proteins = ms$n, n_patients = sum(is.finite(y)),
             beta_protein_saf = b_saf["beta"], se_saf = b_saf["se"], p_saf = b_saf["p"],
             beta_protein_nas = b_nas["beta"], p_nas = b_nas["p"],
             beta_protein_fib = b_fib["beta"], p_fib = b_fib["p"],
             spearman_nas = suppressWarnings(cor(y, mh$NAS, method = "spearman", use = "complete.obs")))
}), fill = TRUE)
slopes <- merge(slopes, sig[, .(mid, cell_type, module, module_name,
                                scrna_beta = disease_stage_beta)], by = "mid")
fwrite(slopes, file.path(OUT, "module_protein_stage_slopes.csv"))

# ── 7. (B') cross-modal concordance: beta_scRNA (fig3g) vs beta_protein (saf) ──
f3 <- fread(file.path(FIG3, "fig3g_module_matrix.csv"))   # cell_type, module, name, beta, q, dir
f3[, mid := paste0(cell_type, "__", module)]
conc <- merge(slopes, f3[, .(mid, name, scrna_beta_fig3g = beta, scrna_q = q, scrna_dir = dir)], by = "mid")
conc <- merge(conc, enrA[cohort == "PXD051911_liver", .(mid, NES, nes_padj = padj)], by = "mid", all.x = TRUE)
conc[, `:=`(
  protein_dir     = ifelse(beta_protein_saf > 0, "up", "down"),
  sign_concordant = sign(scrna_beta_fig3g) == sign(beta_protein_saf))]
conc[, quadrant := fifelse(scrna_beta_fig3g > 0 & beta_protein_saf > 0, "propagated_up",
                    fifelse(scrna_beta_fig3g < 0 & beta_protein_saf < 0, "propagated_down",
                    fifelse(scrna_beta_fig3g > 0 & beta_protein_saf <= 0, "RNA_up_protein_flat",
                                                                          "RNA_down_protein_up")))]
setorder(conc, -scrna_beta_fig3g)
fwrite(conc, file.path(OUT, "crossmodal_module_concordance.csv"))

ct_test <- suppressWarnings(cor.test(conc$scrna_beta_fig3g, conc$beta_protein_saf, method = "spearman"))
cat(sprintf("\n[B'] cross-modal concordance over %d modules: Spearman rho=%.3f p=%.2g | sign-concordant %d/%d (%.0f%%)\n",
            nrow(conc), ct_test$estimate, ct_test$p.value,
            sum(conc$sign_concordant), nrow(conc), 100 * mean(conc$sign_concordant)))

# ── 8. (C) cascade-in-protein: per-module mean score across saf stages ────────
sc_long <- rbindlist(lapply(names(gene_sets), function(mid) {
  ms <- module_score(mid); if (is.null(ms)) return(NULL)
  data.table(mid = mid, raw = mh$raw, saf = mh$saf_diagnosis, saf_ord = mh$saf_ord, score = ms$score)
}))
sc_long <- merge(sc_long, sig[, .(mid, cell_type, module_name,
                                  scrna_dir = ifelse(disease_stage_beta > 0, "up", "down"))], by = "mid")
casc <- sc_long[, .(mean_score = mean(score, na.rm = TRUE),
                    se = sd(score, na.rm = TRUE) / sqrt(sum(is.finite(score))),
                    n = sum(is.finite(score))),
                by = .(mid, cell_type, module_name, scrna_dir, saf, saf_ord)]
setorder(casc, mid, saf_ord)
fwrite(casc, file.path(OUT, "module_cascade_protein.csv"))

# ── 9. built-in positive-control directionality check ─────────────────────────
cat("\n===== POSITIVE-CONTROL directionality (liver enrichment) =====\n")
pc <- enrA[cohort == "PXD051911_liver"]
show_ct <- function(lab, dt) if (nrow(dt)) for (i in seq_len(nrow(dt)))
  cat(sprintf("  %-14s %-40s scRNA=%s  protein_NES=%+.2f (padj %.1g)\n",
              lab, substr(dt$module_name[i], 1, 40), dt$scrna_dir[i], dt$NES[i], dt$padj[i]))
show_ct("[fibroblast]",   pc[cell_type == "fibroblasts"][order(-NES)][1:3])
show_ct("[hepatocyte]",   pc[cell_type == "hepatocytes" & scrna_dir == "down"][order(NES)][1:3])
cat("\nEXPECT: fibroblast ECM/stellate up-modules -> NES > 0; hepatocyte metabolic down-modules -> NES < 0.\n")
cat("\n[done] wrote 4 tables to", OUT, "\n")
