#!/usr/bin/env Rscript

# Build the frozen manuscript evidence classes and the audits that support them.
# The primary classes use only canonical bulk TREAT DE and Tier-1/2 SuSiE-COLOC;
# spatial, single-cell, and proteomic columns are reserved for independent tests.

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- Sys.getenv("MANUSCRIPT_RELEASE_ID", "2026-07-10-r1")
PP4_THRESHOLDS <- c(0.5, 0.7, 0.9)

paths <- list(
  deg = file.path(BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
  coloc = file.path(BASE,
    "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
  trait_tier = file.path(BASE,
    "GWAS/finemapping/config/gwas_trait_tier.tsv"),
  atlas = file.path(BASE,
    "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  convergence = file.path(BASE,
    "RNA-seq/results/multi_evidence/convergence_evidence.csv")
)

missing_inputs <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_inputs)) {
  stop("Missing release input(s): ", paste(missing_inputs, collapse = ", "))
}

RESULT_DIR <- file.path(BASE, "RNA-seq/results/manuscript_release", RELEASE_ID)
DOC_DIR <- file.path(BASE, "docs/manuscript/release")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOC_DIR, recursive = TRUE, showWarnings = FALSE)

clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_neglog10 <- function(x) {
  -log10(pmax(suppressWarnings(as.numeric(x)), .Machine$double.xmin))
}

message("[release] Loading canonical inputs")
deg_raw <- fread(paths$deg)
coloc_raw <- fread(paths$coloc)
tiers <- fread(paths$trait_tier)
atlas <- fread(paths$atlas)

stopifnot(all(c("symbol", "treat_fdr", "treat_p", "t", "AveExpr") %in% names(deg_raw)))
stopifnot(all(c("gwas_name", "gene", "ensembl", "PP.H4.susie", "PP.H4.abf") %in% names(coloc_raw)))
stopifnot(all(c("study_name", "tier", "tier_label", "placement") %in% names(tiers)))
stopifnot(all(c("human_symbol", "gene_biotype") %in% names(atlas)))

# Canonical bulk table, one row per non-empty symbol. The most significant TREAT
# row wins if duplicated symbols occur through annotation aliases.
deg <- copy(deg_raw)
deg[, symbol := clean_symbol(symbol)]
deg <- deg[!is.na(symbol)]
setorder(deg, symbol, treat_fdr, treat_p)
deg <- deg[!duplicated(symbol)]
deg[, in_treat_deg := !is.na(treat_fdr) & treat_fdr < 0.05]
deg <- deg[, .(
  symbol,
  ensembl_bulk = sub("\\..*$", "", gene),
  bulk_logFC = logFC,
  bulk_t = t,
  bulk_AveExpr = AveExpr,
  bulk_treat_p = treat_p,
  bulk_treat_fdr = treat_fdr,
  in_treat_deg
)]

# Restrict genetics to the 35 prespecified main Tier-1/2 studies.
coloc <- merge(
  coloc_raw,
  tiers,
  by.x = "gwas_name",
  by.y = "study_name",
  all = FALSE
)
coloc[, gene := clean_symbol(gene)]
coloc <- coloc[placement == "main" & tier %in% c(1L, 2L) & !is.na(gene)]
if (uniqueN(coloc$gwas_name) != 35L) {
  stop("Expected 35 primary GWAS, observed ", uniqueN(coloc$gwas_name))
}

aggregate_scope <- function(d, suffix) {
  out <- d[, .(
    susie = safe_max(`PP.H4.susie`),
    abf = safe_max(`PP.H4.abf`),
    n_gwas_tested = uniqueN(gwas_name),
    n_ancestries_tested = uniqueN(ancestry)
  ), by = .(symbol = gene)]
  setnames(
    out,
    c("susie", "abf", "n_gwas_tested", "n_ancestries_tested"),
    paste0(c("max_susie_pp4_", "max_abf_pp4_", "n_gwas_tested_", "n_ancestries_tested_"), suffix)
  )
  out
}

g_all <- aggregate_scope(coloc, "all")
g_direct <- aggregate_scope(coloc[tier_label == "direct_MASLD"], "direct")
g_enzyme <- aggregate_scope(coloc[tier_label == "liver_enzyme"], "enzyme")
genetics <- Reduce(
  function(x, y) merge(x, y, by = "symbol", all = TRUE),
  list(g_all, g_direct, g_enzyme)
)

atlas_genes <- unique(atlas[, .(
  symbol = clean_symbol(human_symbol),
  ensembl_atlas = sub("\\..*$", "", ensembl_id),
  gene_biotype
)])[!is.na(symbol)]

release <- merge(atlas_genes, deg, by = "symbol", all = TRUE)
release <- merge(release, genetics, by = "symbol", all = TRUE)
release[, bulk_tested := !is.na(bulk_treat_p)]
release[, coloc_tested := !is.na(max_susie_pp4_all) | !is.na(max_abf_pp4_all)]
release[, joint_testable := bulk_tested & coloc_tested]
release[, in_susie_primary := !is.na(max_susie_pp4_all) & max_susie_pp4_all > 0.5]
release[, in_union_sensitivity := pmax(max_susie_pp4_all, max_abf_pp4_all, na.rm = TRUE) > 0.5]
release[!is.finite(pmax(max_susie_pp4_all, max_abf_pp4_all, na.rm = TRUE)),
        in_union_sensitivity := FALSE]

release[, primary_evidence_class := fifelse(
  !joint_testable, "not_jointly_testable",
  fifelse(in_susie_primary & in_treat_deg, "convergent",
  fifelse(in_susie_primary, "genetic_only",
  fifelse(in_treat_deg, "disease_state_only", "neither")))
)]
release[, sensitivity_evidence_class := fifelse(
  !joint_testable, "not_jointly_testable",
  fifelse(in_union_sensitivity & in_treat_deg, "convergent",
  fifelse(in_union_sensitivity, "genetic_only",
  fifelse(in_treat_deg, "disease_state_only", "neither")))
)]

direct_hit <- !is.na(release$max_susie_pp4_direct) & release$max_susie_pp4_direct > 0.5
enzyme_hit <- !is.na(release$max_susie_pp4_enzyme) & release$max_susie_pp4_enzyme > 0.5
release[, genetic_trait_scope := fifelse(
  direct_hit & enzyme_hit, "both",
  fifelse(direct_hit, "direct_disease", fifelse(enzyme_hit, "enzyme", "none"))
)]
release[, genetic_confidence := fifelse(
  in_susie_primary,
  "susie",
  fifelse(in_union_sensitivity, "abf_only", "none")
)]
release[, analysis_release_id := RELEASE_ID]
setcolorder(release, c(
  "analysis_release_id", "symbol", "ensembl_atlas", "ensembl_bulk",
  "gene_biotype", "joint_testable", "primary_evidence_class",
  "sensitivity_evidence_class", "genetic_trait_scope", "genetic_confidence"
))
setorder(release, symbol)
fwrite(release, file.path(RESULT_DIR, "evidence_class_table.tsv"), sep = "\t", na = "")

message("[release] Building orthogonality audit")
joint <- release[joint_testable == TRUE]

scope_columns <- list(
  all = c("max_susie_pp4_all", "max_abf_pp4_all"),
  direct_disease = c("max_susie_pp4_direct", "max_abf_pp4_direct"),
  enzyme = c("max_susie_pp4_enzyme", "max_abf_pp4_enzyme")
)

audit_rows <- list()
continuous_rows <- list()
idx <- 0L
for (scope in names(scope_columns)) {
  susie_col <- scope_columns[[scope]][1]
  abf_col <- scope_columns[[scope]][2]
  for (definition in c("susie", "susie_or_abf")) {
    score <- if (definition == "susie") {
      joint[[susie_col]]
    } else {
      pmax(joint[[susie_col]], joint[[abf_col]], na.rm = TRUE)
    }
    score[!is.finite(score)] <- NA_real_

    ok_cont <- is.finite(score) & is.finite(joint$bulk_t)
    rho_abs_t <- if (sum(ok_cont) >= 5L) {
      suppressWarnings(cor.test(score[ok_cont], abs(joint$bulk_t[ok_cont]), method = "spearman", exact = FALSE))
    } else NULL
    rho_sig <- if (sum(ok_cont) >= 5L) {
      suppressWarnings(cor.test(score[ok_cont], safe_neglog10(joint$bulk_treat_p[ok_cont]), method = "spearman", exact = FALSE))
    } else NULL
    continuous_rows[[length(continuous_rows) + 1L]] <- data.table(
      analysis_release_id = RELEASE_ID,
      trait_scope = scope,
      genetic_definition = definition,
      n_joint = sum(ok_cont),
      rho_abs_bulk_t = if (is.null(rho_abs_t)) NA_real_ else unname(rho_abs_t$estimate),
      p_abs_bulk_t = if (is.null(rho_abs_t)) NA_real_ else rho_abs_t$p.value,
      rho_bulk_significance = if (is.null(rho_sig)) NA_real_ else unname(rho_sig$estimate),
      p_bulk_significance = if (is.null(rho_sig)) NA_real_ else rho_sig$p.value
    )

    for (threshold in PP4_THRESHOLDS) {
      genetic <- !is.na(score) & score > threshold
      total_score <- if (definition == "susie") {
        release[[susie_col]]
      } else {
        pmax(release[[susie_col]], release[[abf_col]], na.rm = TRUE)
      }
      total_score[!is.finite(total_score)] <- NA_real_
      disease <- joint$in_treat_deg %in% TRUE
      a <- sum(genetic & disease)
      b <- sum(genetic & !disease)
      c <- sum(!genetic & disease)
      d <- sum(!genetic & !disease)
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
      idx <- idx + 1L
      audit_rows[[idx]] <- data.table(
        analysis_release_id = RELEASE_ID,
        trait_scope = scope,
        genetic_definition = definition,
        pp4_threshold = threshold,
        n_joint_testable = nrow(joint),
        n_genetic_total = sum(!is.na(total_score) & total_score > threshold),
        n_genetic_joint = a + b,
        n_treat_deg_total = sum(release$in_treat_deg %in% TRUE),
        n_treat_deg_joint = a + c,
        n_overlap = a,
        pct_genetic_not_deg = 100 * b / max(1L, a + b),
        pct_deg_not_genetic = 100 * c / max(1L, a + c),
        jaccard = a / max(1L, a + b + c),
        fisher_or = unname(ft$estimate),
        fisher_p = ft$p.value
      )
    }
  }
}
orthogonality <- rbindlist(audit_rows)
continuous <- rbindlist(continuous_rows)
fwrite(orthogonality, file.path(RESULT_DIR, "orthogonality_audit.tsv"), sep = "\t", na = "")
fwrite(continuous, file.path(RESULT_DIR, "orthogonality_continuous.tsv"), sep = "\t", na = "")

class_counts <- release[, .N, by = .(
  analysis_release_id,
  primary_evidence_class,
  genetic_trait_scope,
  genetic_confidence
)][order(primary_evidence_class, genetic_trait_scope, genetic_confidence)]
fwrite(class_counts, file.path(RESULT_DIR, "evidence_class_counts.tsv"), sep = "\t", na = "")

message("[release] Testing evidence classes in independent atlas modalities")
endpoint_cols <- c(
  "human_symbol", "n_prot_datasets", "best_protein_padj",
  "spatial_is_svg", "spatial_morans_i", "sc_best_padj", "sc_n_celltypes_sig"
)
validation <- atlas[, ..endpoint_cols]
setnames(validation, "human_symbol", "symbol")
validation[, symbol := clean_symbol(symbol)]
validation <- validation[!is.na(symbol)]
validation <- validation[!duplicated(symbol)]
validation <- merge(validation, release, by = "symbol", all.x = TRUE)

endpoint_specs <- list(
  proteomics = list(tested = !is.na(validation$best_protein_padj),
                    positive = validation$best_protein_padj < 0.05),
  spatial_svg = list(tested = !is.na(validation$spatial_is_svg),
                     positive = validation$spatial_is_svg %in% TRUE),
  single_cell = list(tested = !is.na(validation$sc_best_padj),
                     positive = validation$sc_best_padj < 0.05)
)

class_levels <- c("neither", "genetic_only", "disease_state_only", "convergent")
summary_rows <- list()
pair_rows <- list()
model_rows <- list()

for (endpoint in names(endpoint_specs)) {
  spec <- endpoint_specs[[endpoint]]
  d <- validation[spec$tested & joint_testable == TRUE]
  d[, endpoint_positive := spec$positive[spec$tested & validation$joint_testable == TRUE]]
  d <- d[primary_evidence_class %in% class_levels]
  d[, primary_evidence_class := factor(primary_evidence_class, levels = class_levels)]
  d[, protein_coding := gene_biotype == "protein_coding"]
  d[, `:=`(analysis_release_id = RELEASE_ID, endpoint_name = endpoint)]

  summary_rows[[length(summary_rows) + 1L]] <- d[, .(
    n_tested = .N,
    n_positive = sum(endpoint_positive, na.rm = TRUE),
    positive_rate = mean(endpoint_positive, na.rm = TRUE)
  ), by = .(analysis_release_id, endpoint = endpoint_name, primary_evidence_class)]

  for (comparison in c("genetic_only", "disease_state_only")) {
    pdat <- d[primary_evidence_class %in% c("convergent", comparison)]
    pdat[, primary_evidence_class := droplevels(primary_evidence_class)]
    tab <- table(pdat$primary_evidence_class, pdat$endpoint_positive)
    if (all(dim(tab) == c(2L, 2L))) {
      ft <- fisher.test(tab)
      pair_rows[[length(pair_rows) + 1L]] <- data.table(
        analysis_release_id = RELEASE_ID,
        endpoint,
        comparison = paste0("convergent_vs_", comparison),
        n_convergent = sum(pdat$primary_evidence_class == "convergent"),
        n_comparator = sum(pdat$primary_evidence_class == comparison),
        odds_ratio = unname(ft$estimate),
        p_value = ft$p.value
      )
    }
  }

  fit_data <- d[is.finite(bulk_AveExpr)]
  if (nrow(fit_data) >= 50L && uniqueN(fit_data$endpoint_positive) == 2L) {
    fit <- suppressWarnings(glm(
      endpoint_positive ~ primary_evidence_class + scale(bulk_AveExpr) + protein_coding,
      data = fit_data,
      family = binomial()
    ))
    co <- summary(fit)$coefficients
    keep <- grep("^primary_evidence_class", rownames(co))
    if (length(keep)) {
      model_rows[[length(model_rows) + 1L]] <- data.table(
        analysis_release_id = RELEASE_ID,
        endpoint,
        term = rownames(co)[keep],
        log_odds = co[keep, "Estimate"],
        standard_error = co[keep, "Std. Error"],
        odds_ratio = exp(co[keep, "Estimate"]),
        ci_low = exp(co[keep, "Estimate"] - 1.96 * co[keep, "Std. Error"]),
        ci_high = exp(co[keep, "Estimate"] + 1.96 * co[keep, "Std. Error"]),
        p_value = co[keep, "Pr(>|z|)"]
      )
    }
  }
}

validation_summary <- rbindlist(summary_rows, fill = TRUE)
validation_pairwise <- rbindlist(pair_rows, fill = TRUE)
validation_models <- rbindlist(model_rows, fill = TRUE)
if (nrow(validation_pairwise)) validation_pairwise[, q_value := p.adjust(p_value, method = "BH")]
if (nrow(validation_models)) validation_models[, q_value := p.adjust(p_value, method = "BH")]
fwrite(validation_summary, file.path(RESULT_DIR, "evidence_class_validation_summary.tsv"), sep = "\t", na = "")
fwrite(validation_pairwise, file.path(RESULT_DIR, "evidence_class_validation_pairwise.tsv"), sep = "\t", na = "")
fwrite(validation_models, file.path(RESULT_DIR, "evidence_class_validation_adjusted.tsv"), sep = "\t", na = "")

message("[release] Writing manifest, gates, and claim ledger")
file_rows <- lapply(names(paths), function(role) {
  f <- paths[[role]]
  info <- file.info(f)
  data.table(
    analysis_release_id = RELEASE_ID,
    role,
    path = sub(paste0("^", BASE, "/?"), "", f),
    md5 = unname(tools::md5sum(f)),
    bytes = info$size,
    modified = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z")
  )
})
manifest <- rbindlist(file_rows)
fwrite(manifest, file.path(DOC_DIR, "analysis_release_manifest.tsv"), sep = "\t", na = "")

primary_rows <- orthogonality[
  pp4_threshold == 0.5 &
  ((trait_scope == "all" & genetic_definition %in% c("susie", "susie_or_abf")) |
   (trait_scope == "direct_disease" & genetic_definition == "susie_or_abf"))
]
orthogonality_pass <- nrow(primary_rows) == 3L && all(primary_rows$pct_genetic_not_deg >= 85)

class_gate_endpoints <- character()
if (nrow(validation_pairwise)) {
  passed <- validation_pairwise[odds_ratio > 1 & q_value < 0.05]
  class_gate_endpoints <- passed[, .N, by = endpoint][N >= 2L, endpoint]
}
class_validation_pass <- length(class_gate_endpoints) >= 2L

gates <- data.table(
  analysis_release_id = RELEASE_ID,
  gate = c("orthogonality_core", "class_validation_article"),
  passed = c(orthogonality_pass, class_validation_pass),
  criterion = c(
    ">=85% genetic-not-DE for SuSiE, union, and direct-disease union at PP4>0.5",
    "convergent class enriched versus both single-map classes at BH q<0.05 in >=2 modalities"
  ),
  detail = c(
    paste(sprintf("%s/%s=%.1f%%", primary_rows$trait_scope,
                  primary_rows$genetic_definition, primary_rows$pct_genetic_not_deg), collapse = "; "),
    if (length(class_gate_endpoints)) paste(class_gate_endpoints, collapse = ",") else "none"
  )
)
fwrite(gates, file.path(RESULT_DIR, "acceptance_gates.tsv"), sep = "\t", na = "")

pick <- function(scope, definition) {
  orthogonality[
    trait_scope == scope & genetic_definition == definition & pp4_threshold == 0.5
  ][1]
}
p_susie <- pick("all", "susie")
p_union <- pick("all", "susie_or_abf")
p_direct <- pick("direct_disease", "susie_or_abf")

claims <- data.table(
  analysis_release_id = RELEASE_ID,
  claim_id = c(
    "C_ORTHO_SUSIE", "C_ORTHO_UNION", "C_ORTHO_DIRECT",
    "C_CLASS_FRAME", "C_SCORE_SCOPE", "C_DRUG_RETIRED", "C_GSMAP_RETIRED",
    "C_STAGE_SCOPE", "C_CODING_SCOPE"
  ),
  status = c("supported", "supported", "supported", "supported", "supported",
             "retired_pending_rebuild", "retired_pending_rebuild", "supported", "supported"),
  allowed_wording = c(
    sprintf("The SuSiE set contains %d genes overall; among %d jointly tested genes, %d overlap canonical TREAT DEGs and %.1f%% are not differentially expressed.",
            p_susie$n_genetic_total, p_susie$n_genetic_joint,
            p_susie$n_overlap, p_susie$pct_genetic_not_deg),
    sprintf("The SuSiE/ABF union contains %d genes overall; among %d jointly tested genes, %d overlap TREAT DEGs and %.1f%% are not differentially expressed.",
            p_union$n_genetic_total, p_union$n_genetic_joint, p_union$n_overlap, p_union$pct_genetic_not_deg),
    sprintf("The direct-disease/PDFF union contains %d genes overall; among %d jointly tested genes, %d overlap TREAT DEGs and %.1f%% are not differentially expressed.",
            p_direct$n_genetic_total, p_direct$n_genetic_joint, p_direct$n_overlap, p_direct$pct_genetic_not_deg),
    "Genetic and disease-state evidence define susceptibility-linked, disease-state-associated, and convergent candidate classes.",
    "The convergence score is an evidence-weighted heuristic and not a posterior probability or gene-level significance measure.",
    "Drug-target calibration is under strict-label rebuild and is not a headline claim in this release.",
    "gsMap arm enrichment is under same-universe replication and is not a headline claim in this release.",
    "Cross-sectional adjacent-stage contrasts describe stage-ordered remodeling.",
    "Variant annotations describe the analyzed colocalized loci, not all inherited MASLD risk."
  ),
  prohibited_wording = c(
    "genetics and expression are statistically independent",
    "1,162-gene legacy denominator as the primary set",
    "the broad genetic map is entirely direct MASLD susceptibility",
    "convergence is required for credible targets",
    "Bayesian target probability; gene-level significant convergence",
    "recovers both approved targets; full drug-development gradient; prospective validation",
    "inherited risk is tissue-real only at the intersection",
    "multicellular cascade; therapeutic windows",
    "approximately 97% of all inherited MASLD risk is noncoding"
  ),
  source_output = c(
    rep("orthogonality_audit.tsv", 3),
    "evidence_class_table.tsv", "convergence_evidence.csv", "drug benchmark pending",
    "gsMap replication pending", "canonical adjacent-stage contrasts", "credible-set audit pending"
  )
)
fwrite(claims, file.path(DOC_DIR, "claim_ledger.tsv"), sep = "\t", na = "")

message("[release] Complete: ", RESULT_DIR)
print(gates)
print(primary_rows[, .(
  trait_scope, genetic_definition, n_genetic_total, n_genetic_joint,
  n_treat_deg_total, n_treat_deg_joint, n_overlap,
  pct_genetic_not_deg, pct_deg_not_genetic, fisher_or, fisher_p
)])
