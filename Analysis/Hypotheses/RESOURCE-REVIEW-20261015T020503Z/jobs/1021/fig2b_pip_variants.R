#!/usr/bin/env Rscript
# Review R02, job 1021. Re-derive Figure 2B points from the COLOC r2 release and
# recompute the x axis under method-specific fine-mapping PIPs.
# Replicates the data steps of scripts/figures/fig2_resolution_plane.R (no plotting).
suppressPackageStartupMessages(library(data.table))
setDTthreads(4)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT  <- file.path(BASE, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1021")
COLOC_INPUT <- file.path(BASE, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv")
INSTALLED <- file.path(BASE, "figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc_source.csv")

coloc <- fread(COLOC_INPUT, select = c("gwas_name", "gene", "PP.H4.susie", "PP.H4.abf", "top_snp", "ancestry"))
MAIN <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[placement == "main", study_name]
coloc <- coloc[gwas_name %in% MAIN]
setnames(coloc, c("PP.H4.susie", "PP.H4.abf"), c("pp4_susie", "pp4_abf"))
coloc <- coloc[!is.na(gene) & gene != ""]
coloc[, pp4_susie := suppressWarnings(as.numeric(pp4_susie))]
coloc <- coloc[!is.na(top_snp) & top_snp != ""]
coloc[, c("snp_chr", "snp_pos") := tstrsplit(top_snp, ":", fixed = TRUE, type.convert = TRUE)]
keys <- unique(coloc[, .(gwas_name, snp_chr, snp_pos)])
cat("coloc rows (main, gene, top_snp):", nrow(coloc), " unique study-variant keys:", nrow(keys), "\n")

cols <- c("study", "chromosome", "position", "ancestry", "susie_pip", "carma_pip",
          "susiex_pip", "mesusie_pip", "recommended_pip")
t0 <- Sys.time()
fm1 <- fread(file.path(BASE, "GWAS/finemapping/results/combined_finemapping.csv"), select = cols)
cat("combined_finemapping rows:", nrow(fm1), " read secs:", round(as.numeric(Sys.time() - t0, units = "secs")), "\n")
fm1 <- fm1[keys, on = .(study = gwas_name, chromosome = snp_chr, position = snp_pos), nomatch = 0L]
fm1[, x_source := "canonical_finemapping"]
fm2 <- fread(file.path(BASE, "GWAS/finemapping/results/combined_finemapping_coloc_targeted.csv"),
             select = c("study", "chromosome", "position", "ancestry", "susie_pip", "carma_pip", "recommended_pip"))
fm2 <- fm2[keys, on = .(study = gwas_name, chromosome = snp_chr, position = snp_pos), nomatch = 0L]
fm2[, `:=`(susiex_pip = NA_real_, mesusie_pip = NA_real_, x_source = "coloc_targeted_finemapping")]
fm <- rbindlist(list(fm1, fm2), use.names = TRUE)
num <- c("susie_pip", "carma_pip", "susiex_pip", "mesusie_pip", "recommended_pip")
fm[, (num) := lapply(.SD, function(v) suppressWarnings(as.numeric(v))), .SDcols = num]
fm[, source_rank := fifelse(x_source == "canonical_finemapping", 1L, 2L)]
cat("matched fine-mapping rows:", nrow(fm), "\n")

pmax_na <- function(...) { v <- suppressWarnings(pmax(..., na.rm = TRUE)); v[!is.finite(v)] <- NA_real_; v }
# Variant-level x definitions
fm[, x_installed := pmax_na(recommended_pip, susie_pip)]          # producer: max over SuSiE, CARMA, SuSiEx, MESuSiE
fm[, x_within := pmax_na(susie_pip, carma_pip)]                     # within-ancestry methods only (still a max of two)
fm[, x_susie := fifelse(is.finite(susie_pip), susie_pip, NA_real_)] # SuSiE only
fm[, joint_sets_max := is.finite(x_installed) & (is.na(x_within) | x_installed > x_within + 1e-12)]

dedupe <- function(col) {
  d <- fm[is.finite(get(col)), .(study, chromosome, position, v = get(col), x_source, source_rank, joint_sets_max)]
  setorder(d, study, chromosome, position, -v, source_rank)
  unique(d, by = c("study", "chromosome", "position"))
}

max_or_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
label_at_max <- function(value, label) if (all(is.na(value))) NA_character_ else label[which.max(fifelse(is.na(value), -Inf, value))]
gene_y <- coloc[, .(pp4_susie = max_or_na(pp4_susie), susie_gwas = label_at_max(pp4_susie, gwas_name),
                    susie_snp = label_at_max(pp4_susie, top_snp)), by = gene]

roster <- c("HSD17B13", "TM6SF2", "MBOAT7", "GCKR", "RORA")
build <- function(col, matched_study = FALSE) {
  u <- dedupe(col)
  cm <- merge(coloc, u, by.x = c("gwas_name", "snp_chr", "snp_pos"), by.y = c("study", "chromosome", "position"))
  if (matched_study) {
    cm <- merge(cm, gene_y[, .(gene, susie_gwas, susie_snp)], by = "gene")
    cm <- cm[gwas_name == susie_gwas & top_snp == susie_snp]
  }
  gx <- cm[order(gene, -v, x_source), .SD[1], by = gene][, .(gene, x = v, x_snp = top_snp, x_gwas = gwas_name, joint_sets_max)]
  dt <- merge(gene_y, gx, by = "gene", all.x = TRUE)
  pl <- dt[!is.na(pp4_susie) & is.finite(x)]
  pl[, pp4 := pmax(pp4_susie, 0)]
  pl[, plotted := (pp4 > 0.05) | (x >= 0.99) | gene %in% roster]
  pl
}

A <- build("x_installed")
inst <- fread(INSTALLED)
chk <- merge(inst[, .(gene, causal_pip, pp4, plotted)], A[, .(gene, x, pp4r = pp4, plottedr = plotted)], by = "gene", all = TRUE)
cat(sprintf("REPRO: installed rows %d, rederived rows %d, matched %d; max|dx| %.4g; max|dpp4| %.4g; plotted installed %d rederived %d; plotted disagreements %d\n",
            nrow(inst), nrow(A), chk[!is.na(causal_pip) & !is.na(x), .N],
            chk[, max(abs(causal_pip - round(x, 4)), na.rm = TRUE)], chk[, max(abs(pp4 - round(pp4r, 4)), na.rm = TRUE)],
            inst[, sum(plotted)], A[, sum(plotted)], chk[, sum(plotted != plottedr, na.rm = TRUE)]))
cat(sprintf("Installed plotted points whose x comes from a cross-ancestry joint PIP (SuSiEx/MESuSiE) exceeding SuSiE and CARMA: %d of %d\n",
            A[plotted == TRUE, sum(joint_sets_max, na.rm = TRUE)], A[, sum(plotted)]))
cat(sprintf("Installed plotted points whose x and y come from different GWAS: %d of %d; different lead SNP: %d\n",
            A[plotted == TRUE, sum(x_gwas != susie_gwas)], A[, sum(plotted)], A[plotted == TRUE, sum(x_snp != susie_snp)]))

compare <- function(B, label) {
  m <- merge(A[, .(gene, xA = x, plA = plotted, pp4)], B[, .(gene, xB = x, plB = plotted)], by = "gene", all.x = TRUE)
  pa <- m[plA == TRUE]
  cat(sprintf("\n[%s] genes with finite x: %d (installed %d); plotted: %d (installed %d)\n", label, nrow(B), nrow(A), B[, sum(plotted)], A[, sum(plotted)]))
  cat(sprintf("  of %d installed plotted points: x missing under %s %d; |dx|>0.1 %d; |dx|>0.5 %d; any change >1e-4 %d\n",
              nrow(pa), label, pa[is.na(xB), .N], pa[!is.na(xB) & abs(xA - xB) > 0.1, .N], pa[!is.na(xB) & abs(xA - xB) > 0.5, .N],
              pa[is.na(xB) | abs(xA - xB) > 1e-4, .N]))
  cat(sprintf("  cross PIP 0.5: %d; cross 0.99 gate: %d; leave plot: %d; enter plot: %d\n",
              pa[!is.na(xB) & ((xA >= 0.5) != (xB >= 0.5)), .N], pa[!is.na(xB) & ((xA >= 0.99) != (xB >= 0.99)), .N],
              m[plA == TRUE & (is.na(plB) | plB == FALSE), .N], m[(plA == FALSE | is.na(plA)) & plB == TRUE, .N]))
  cat(sprintf("  plotted genes with x>=0.99: installed %d, %s %d; x>=0.5: installed %d, %s %d\n",
              A[plotted == TRUE & x >= 0.99, .N], label, B[plotted == TRUE & x >= 0.99, .N],
              A[plotted == TRUE & x >= 0.5, .N], label, B[plotted == TRUE & x >= 0.5, .N]))
  # high-PIP low-PP.H4 corner (the panel's coding-blind claim) and high-high corner
  cat(sprintf("  corner PIP>=0.9 & PP.H4<0.5: installed %d, %s %d; PIP>=0.9 & PP.H4>0.5: installed %d, %s %d\n",
              A[plotted == TRUE & x >= 0.9 & pp4 < 0.5, .N], label, B[plotted == TRUE & x >= 0.9 & pp4 < 0.5, .N],
              A[plotted == TRUE & x >= 0.9 & pp4 > 0.5, .N], label, B[plotted == TRUE & x >= 0.9 & pp4 > 0.5, .N]))
  m[, variant := label]
  m
}
res <- rbindlist(list(
  compare(build("x_within"), "within_ancestry_max_susie_carma"),
  compare(build("x_susie"), "susie_only"),
  compare(build("x_susie", matched_study = TRUE), "susie_only_same_gwas_and_lead_as_y")), fill = TRUE)
for (g in roster) {
  cat(sprintf("roster %s: installed x=%s y=%s | susie_only x=%s | matched x=%s\n", g,
              A[gene == g, signif(x, 3)], A[gene == g, signif(pp4, 3)],
              build("x_susie")[gene == g, signif(x, 3)], build("x_susie", TRUE)[gene == g, signif(x, 3)]))
}
fwrite(res, file.path(OUT, "fig2b_x_variants_by_gene.tsv"), sep = "\t")
fwrite(A, file.path(OUT, "fig2b_rederived_installed_definition.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("DONE\n")
