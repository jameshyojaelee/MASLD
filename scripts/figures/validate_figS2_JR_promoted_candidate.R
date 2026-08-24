#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

base <- Sys.getenv("MASLD_PROJECT_ROOT")
out <- Sys.getenv("FIG2_SUPP_OUT_DIR")
stopifnot(nzchar(base), nzchar(out), dir.exists(out))
source(file.path(base, "scripts/figures/fig2_promoted_coloc_context.R"))

close_num <- function(x, y, tol = 1e-12) {
  # Fail closed on a missing/renamed column: is.na(NULL) is logical(0) and
  # all(logical(0)) is vacuously TRUE, which would pass this check silently.
  if (length(x) == 0L || length(y) == 0L || length(x) != length(y)) return(FALSE)
  same_na <- is.na(x) == is.na(y)
  all(same_na & (is.na(x) | abs(x - y) <= tol))
}

metrics <- data.table(panel = character(), metric = character(), value = character())
add_metric <- function(panel, metric, value) {
  metrics <<- rbind(metrics, data.table(
    panel = panel, metric = metric, value = as.character(value)))
}

# FigS2J: independently reconstruct the full EUR-restricted versus model-shared
# table and its displayed canonical-locus denominator.
gs <- fread(file.path(base,
  "GWAS/finemapping/results/susiex_mvp/susiex_gene_summary_mvp.csv"),
  select = c("GeneSymbol", "susiex_max_pip", "locus_ids"))
mg <- fread(file.path(base,
  "GWAS/finemapping/results/mesusie_mvp/mesusie_gene_summary_mvp.csv"))
gs <- gs[GeneSymbol != "" & !is.na(GeneSymbol)]
mg <- mg[GeneSymbol != "" & !is.na(GeneSymbol)]
mp_cols <- grep("^max_pip_", names(mg), value = TRUE)
mp_tok <- vapply(sub("^max_pip_", "", mp_cols),
  function(x) length(strsplit(x, "_", fixed = TRUE)[[1]]), integer(1))
mp_multi <- mp_cols[mp_tok >= 2L]
j <- merge(gs, mg[, c("GeneSymbol", "mesusie_max_pip", mp_cols), with = FALSE],
  by = "GeneSymbol", all = TRUE)
j[, eur_pip := suppressWarnings(as.numeric(max_pip_EUR))]
shmat <- as.matrix(j[, ..mp_multi])
shmat[is.na(shmat)] <- -Inf
j[, joint_pip := apply(shmat, 1L, max)]
j[!is.finite(joint_pip), joint_pip := NA_real_]
j[, delta := joint_pip - eur_pip]
j <- j[!is.na(joint_pip) & !is.na(eur_pip)]
jgot <- fread(file.path(out, "FigS2J_crossancestry_pip_concentration_source.csv"))
jexp <- j[, .(gene = GeneSymbol, eur_pip = round(eur_pip, 4),
              joint_pip = round(joint_pip, 4), delta = round(delta, 4))]
jgot <- jgot[, .(gene, eur_pip, joint_pip, delta)]
setorder(jgot, gene, eur_pip, joint_pip, delta)
setorder(jexp, gene, eur_pip, joint_pip, delta)
stopifnot(nrow(jgot) == nrow(j), identical(jgot$gene, jexp$gene),
  close_num(jgot$eur_pip, jexp$eur_pip),
  close_num(jgot$joint_pip, jexp$joint_pip),
  close_num(jgot$delta, jexp$delta))
canon <- c("GCKR", "PNPLA3", "TM6SF2", "MBOAT7", "HSD17B13", "TRIB1", "GPAM",
  "HNF1A", "HNF1B", "GATAD2A", "MLXIPL", "ABO", "HFE", "RECQL4", "GGT1",
  "SERPINA1", "MARC1", "APOE", "PPP1R3B", "TOR1B")
jplot <- j[GeneSymbol %chin% canon]
if (anyDuplicated(jplot$GeneSymbol)) {
  jplot <- jplot[, .SD[which.max(joint_pip)], by = GeneSymbol]
}
add_metric("FigS2J", "full_gene_records", nrow(j))
add_metric("FigS2J", "displayed_canonical_loci", nrow(jplot))
add_metric("FigS2J", "joint_raised", sum(jplot$delta > 0.005))
add_metric("FigS2J", "joint_lowered", sum(jplot$delta < -0.005))
add_metric("FigS2J", "joint_flat", sum(abs(jplot$delta) <= 0.005))

# FigS2K and FigS2N must remain blocked rather than silently reusing stale data.
k <- fread(file.path(out, "FigS2K_cs_regulatory_composition_BLOCKED.tsv"))
n <- fread(file.path(out,
  "FigS2N_progression_driver_coloc_phenotype_class_heatmap_BLOCKED.tsv"))
kraw <- fread(file.path(base,
  "RNA-seq/results/coloc_variant_classes/lead_causal_annotation.csv"))
stopifnot(k$status == "blocked_retired", k$invalid_input_rows == nrow(kraw),
  n$status == "blocked_retired", n$promoted_registry_strata == 35L)
add_metric("FigS2K", "status", k$status)
add_metric("FigS2K", "invalid_historical_rows", k$invalid_input_rows)
add_metric("FigS2N", "status", n$status)
add_metric("FigS2N", "promoted_registry_strata", n$promoted_registry_strata)

# FigS2L: check the three displayed posterior values directly against the
# promoted 35-stratum table.
sc <- fread(file.path(base,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
lgot <- fread(file.path(out, "FigS2L_f2rl1_spotlight_source.csv"))
l_exp <- c(
  sc[gene == "F2RL1" & gwas_name == "UKBB_GGT", PP.H4.abf][1],
  sc[gene == "F2RL1" & gwas_name == "UKBB_GGT", PP.H4.susie][1],
  sc[gene == "F2RL1" & gwas_name == "PanUKBB_CSA_GGT", PP.H4.abf][1])
stopifnot(close_num(lgot[evidence == "coloc", value], round(l_exp, 4)))
add_metric("FigS2L", "EUR_GGT_single_signal_PPH4", l_exp[1])
add_metric("FigS2L", "EUR_GGT_multi_signal_PPH4", l_exp[2])
add_metric("FigS2L", "SAS_GGT_single_signal_PPH4", l_exp[3])

# FigS2M: reconstruct every credible-set category and denominator.
ml <- fread(file.path(base,
  "GWAS/finemapping/results/mesusie_mvp/mesusie_locus_summary_mvp.csv"))[converged == TRUE]
ncs <- grep("^n_cs_", names(ml), value = TRUE)
ntok <- vapply(sub("^n_cs_", "", ncs),
  function(x) length(strsplit(x, "_", fixed = TRUE)[[1]]), integer(1))
spec <- ncs[ntok == 1L]
shared <- ncs[ntok >= 2L]
sh <- sum(unlist(ml[, ..shared]), na.rm = TRUE)
spec_n <- vapply(spec, function(z) sum(ml[[z]], na.rm = TRUE), numeric(1))
mexp <- data.table(
  type = c("Model-shared (≥2 ancestries)", paste0(sub("^n_cs_", "", spec), "-restricted")),
  n = c(sh, spec_n))
mgot <- fread(file.path(out, "FigS2M_mesusie_shared_specific_source.csv"))
mcmp <- merge(mgot[, .(type, n_got = n)], mexp[, .(type, n_exp = n)], by = "type", all = TRUE)
stopifnot(nrow(mcmp) == nrow(mexp), all(mcmp$n_got == mcmp$n_exp))
add_metric("FigS2M", "credible_sets_total", sum(mexp$n))
for (i in seq_len(nrow(mexp))) add_metric("FigS2M", mexp$type[i], mexp$n[i])
add_metric("FigS2M", "converged_loci", nrow(ml))

# FigS2O/P/Q/R: every Resource posterior shown or deposited must match the
# promoted loader. Q/R also retain the current canonical bulk gate.
check_context <- function(file, panel, genes, mapping) {
  got <- fread(file.path(out, file))
  ctx <- load_fig2_promoted_context(base, genes)
  z <- merge(got, ctx, by = "gene", all = TRUE, suffixes = c("_got", "_exp"))
  stopifnot(nrow(z) == length(unique(genes)))
  for (got_col in names(mapping)) {
    exp_col <- mapping[[got_col]]
    stopifnot(close_num(z[[got_col]], z[[exp_col]]))
  }
  add_metric(panel, "displayed_or_deposited_genes", nrow(got))
  invisible(got)
}

o <- fread(file.path(out, "FigS2O_pqtl_external_complementation_source.csv"))
check_context("FigS2O_pqtl_external_complementation_source.csv", "FigS2O", o$gene,
  c(expr_qtl_coloc_ours = "all_best", multi_signal_pph4 = "all_susie",
    single_signal_pph4 = "all_abf", direct_multi_signal_pph4 = "direct_susie",
    direct_single_signal_pph4 = "direct_abf", enzyme_multi_signal_pph4 = "enzyme_susie",
    enzyme_single_signal_pph4 = "enzyme_abf"))
p <- fread(file.path(out, "FigS2P_pqtl_eqtl_crossref_source.csv"))
check_context("FigS2P_pqtl_eqtl_crossref_source.csv", "FigS2P", p$gene,
  c(expr_qtl_coloc_ours = "all_best", multi_signal_pph4 = "all_susie",
    single_signal_pph4 = "all_abf", direct_multi_signal_pph4 = "direct_susie",
    direct_single_signal_pph4 = "direct_abf", enzyme_multi_signal_pph4 = "enzyme_susie",
    enzyme_single_signal_pph4 = "enzyme_abf"))
q <- fread(file.path(out, "FigS2Q_scchromatin_external_corroboration_source.csv"))
check_context("FigS2Q_scchromatin_external_corroboration_source.csv", "FigS2Q", q$gene,
  c(our_genetic_direct_pph4 = "direct_best", our_genetic_enzyme_pph4 = "enzyme_best",
    direct_multi_signal_pph4 = "direct_susie",
    direct_single_signal_pph4 = "direct_abf", enzyme_multi_signal_pph4 = "enzyme_susie",
    enzyme_single_signal_pph4 = "enzyme_abf"))
r <- fread(file.path(out, "FigS2R_mpra_external_corroboration_source.csv"))
check_context("FigS2R_mpra_external_corroboration_source.csv", "FigS2R", r$gene,
  c(our_genetic_direct_pph4 = "direct_best", our_genetic_enzyme_pph4 = "enzyme_best",
    direct_multi_signal_pph4 = "direct_susie",
    direct_single_signal_pph4 = "direct_abf", enzyme_multi_signal_pph4 = "enzyme_susie",
    enzyme_single_signal_pph4 = "enzyme_abf"))

deg <- load_fig2_current_deg_by_symbol(base, unique(c(q$gene, r$gene)))
for (x in list(FigS2Q = q, FigS2R = r)) {
  z <- merge(x, deg[, .(gene = gene_symbol, is_deg = current_is_deg)],
    by = "gene", all.x = TRUE)
  called <- grepl("^DEG ", z$our_disease_state)
  stopifnot(all(called == fifelse(is.na(z$is_deg), FALSE, z$is_deg)))
}

fwrite(metrics, file.path(out, "validation_summary.tsv"), sep = "\t")
cat("PASS: FigS2J-R promoted-release source rederivation\n")
