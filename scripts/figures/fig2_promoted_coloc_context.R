# Shared promoted-COLOC loader for Figure 2 supplementary panels.
# Returns one row per gene with method- and phenotype-specific maxima. Missing
# posteriors remain missing; they are never converted to negative evidence.

suppressPackageStartupMessages(library(data.table))

# The canonical DEG gate has exactly one definition (load_figure_data.R::
# is_canonical_deg). Never re-derive it inline here: hardcoded cutoffs silently
# desync if CANONICAL_DEG_PADJ / CANONICAL_DEG_LFC ever move.
if (!exists("is_canonical_deg", mode = "function")) {
  source(file.path(
    Sys.getenv("MASLD_PROJECT_ROOT",
               "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
    "scripts/figures/load_figure_data.R"))
}

fig2_finite_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

fig2_support_state <- function(susie, abf, n_evaluated) {
  fifelse(!is.na(susie) & susie > 0.5, "multi-signal COLOC",
    fifelse(!is.na(abf) & abf > 0.5, "single-signal COLOC only",
      fifelse(n_evaluated > 0, "evaluated without PP.H4 > 0.5 support",
              "not evaluable")))
}

load_fig2_promoted_pairs <- function(base, genes = NULL) {
  coloc_file <- file.path(base,
    "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv")
  registry_file <- file.path(base,
    "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/genetics_context/phenotype_registry.tsv")
  stopifnot(file.exists(coloc_file), file.exists(registry_file))

  reg <- fread(registry_file)[placement == "main",
    .(gwas_name = study_name, trait, tier, phenotype_stratum, ancestry)]
  stopifnot(nrow(reg) == 35L, uniqueN(reg$gwas_name) == 35L,
            all(reg$tier %in% c(1L, 2L)))

  x <- fread(coloc_file,
    select = c("gwas_name", "gene", "ensembl", "PP.H4.abf", "PP.H4.susie",
               "method", "ancestry", "ld_panel", "top_snp", "top_snp_PP"))
  x <- x[!is.na(gene) & gene != ""]
  if (!is.null(genes)) x <- x[gene %chin% genes]
  x <- merge(x, reg, by = "gwas_name", all = FALSE, suffixes = c("_coloc", "_registry"))
  stopifnot(uniqueN(x$gwas_name) <= 35L)
  x[, phenotype_group := fifelse(
    phenotype_stratum %chin% c("direct_masld_mash_diagnosis",
                               "mri_pdff_or_histologic_steatosis"),
    "direct_or_liver_fat", "enzyme")]
  x[]
}

summarize_fig2_promoted_context <- function(pairs) {
  if (!nrow(pairs)) return(data.table())
  one_group <- function(d, prefix) {
    out <- d[, .(
      susie = fig2_finite_max(PP.H4.susie),
      abf = fig2_finite_max(PP.H4.abf),
      n_evaluated = sum(is.finite(suppressWarnings(as.numeric(PP.H4.susie))) |
                        is.finite(suppressWarnings(as.numeric(PP.H4.abf))))
    ), by = gene]
    setnames(out, c("susie", "abf", "n_evaluated"),
             paste0(prefix, c("_susie", "_abf", "_n_evaluated")))
    out
  }

  all <- one_group(pairs, "all")
  direct <- one_group(pairs[phenotype_group == "direct_or_liver_fat"], "direct")
  enzyme <- one_group(pairs[phenotype_group == "enzyme"], "enzyme")
  out <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE),
                list(all, direct, enzyme))
  for (prefix in c("all", "direct", "enzyme")) {
    ncol <- paste0(prefix, "_n_evaluated")
    out[is.na(get(ncol)), (ncol) := 0L]
    out[, (paste0(prefix, "_best")) := {
      s <- get(paste0(prefix, "_susie")); a <- get(paste0(prefix, "_abf"))
      fifelse(is.na(s), a, fifelse(is.na(a), s, pmax(s, a)))
    }]
    out[, (paste0(prefix, "_support")) := fig2_support_state(
      get(paste0(prefix, "_susie")), get(paste0(prefix, "_abf")), get(ncol))]
  }
  out[]
}

load_fig2_promoted_context <- function(base, genes = NULL) {
  out <- summarize_fig2_promoted_context(load_fig2_promoted_pairs(base, genes))
  if (!is.null(genes)) {
    out <- merge(data.table(gene = unique(as.character(genes))), out,
                 by = "gene", all.x = TRUE)
    for (prefix in c("all", "direct", "enzyme")) {
      ncol <- paste0(prefix, "_n_evaluated")
      scol <- paste0(prefix, "_support")
      out[is.na(get(ncol)), (ncol) := 0L]
      out[is.na(get(scol)), (scol) := "not evaluable"]
    }
  }
  out[]
}

load_fig2_current_deg_by_symbol <- function(base, genes = NULL) {
  deg <- fread(file.path(base,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
    select = c("gene", "symbol", "logFC", "padj"))
  deg <- deg[!is.na(symbol) & symbol != ""]
  if (!is.null(genes)) deg <- deg[symbol %chin% genes]
  deg[, current_is_deg := is_canonical_deg(deg)]
  # Multiple Ensembl records can map to one deposited symbol. Retain a record
  # that passes the canonical gate when any passes; otherwise retain the record
  # with the smallest adjusted P value. Keep the Ensembl ID for provenance.
  deg[, abs_logFC := abs(logFC)]
  setorder(deg, symbol, -current_is_deg, padj, -abs_logFC, gene, na.last = TRUE)
  deg <- deg[, .SD[1L], by = symbol]
  deg[, abs_logFC := NULL]
  setnames(deg, "symbol", "gene_symbol")
  deg[]
}
