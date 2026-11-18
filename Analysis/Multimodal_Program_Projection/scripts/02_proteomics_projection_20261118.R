#!/usr/bin/env Rscript

# PXD051911 liver DIA-MS projection (20261118 candidate, N10/N18/N22).
#
# Key methodological updates over 02_proteomics_projection.R:
# 1. Deterministic GENCODE v49 primary-assembly symbol mapping
#    (data/gencode_v49_gene_metadata.tsv.gz; no HGNChelper runtime branching).
# 2. Explicit audit + exclusion of `;`-delimited multi-gene shared peptide groups
#    (grepl(";", Genes, fixed = TRUE)) from gene-specific DE/concordance/program scores.
# 3. Bulk DEG input pinned to the 774-unit primary-source-controls fit
#    (RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/
#    primary_source_controls/deg_results.csv), resolving alt-haplotype HLA-DQA1 cleanly
#    to its primary-assembly GENCODE v49 gene ENSG00000196735.14.
# 4. Separate within-MASLD (n=46) normalization track: raw log2 intensities are subset
#    to the 46 MASLD samples FIRST, filtered (<=50% missing within MASLD), quantile-
#    normalized across the 46 MASLD samples alone, collapsed by median, and z-scored
#    (z_masld) using 46-MASLD means/SDs only before scoring programs and computing
#    within-MASLD histology associations. The 58-sample quantile norm is retained ONLY
#    for the 12-control vs 46-MASLD DE contrast and 58-sample descriptive abundance matrix.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(limma)
  library(fgsea)
})
set.seed(42)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
ROOT <- file.path(BASE, "Analysis/Multimodal_Program_Projection")
OUT <- Sys.getenv("PROTEOMICS_OUT_DIR", file.path(ROOT, "results/proteomics_20261118"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

REGISTRY_FILE <- file.path(ROOT, "results/frozen_programs.tsv")
MEMBERSHIP_FILE <- file.path(ROOT, "results/frozen_program_membership.tsv")
RAW_FILE <- file.path(BASE, "data/PXD051911/liver_protein_quant.txt")
META_FILE <- file.path(BASE, "data/PXD051911/meta_data.txt")
GENCODE_FILE <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
CONTRACT_FILE <- file.path(ROOT, "config/fig4_input_contract.tsv")
FIXED_ROWS_FILE <- file.path(ROOT, "config/panel4c_fixed_rows.tsv")
BULK_ID_FILE <- Sys.getenv("FIG5C_BULK_ID_FILE", file.path(ROOT, "config/panel4c_bulk_gene_ids.tsv"))
BULK_FILE <- Sys.getenv(
  "BULK_DEG_FILE",
  file.path(
    BASE,
    "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
  )
)
MSIGDB_RDS <- Sys.getenv(
  "MSIGDB_RDS",
  "/gpfs/commons/home/jameslee/.cache/R/msigdbr/msigdb.2025.1.Hs.rds"
)

required <- c(
  REGISTRY_FILE, MEMBERSHIP_FILE, RAW_FILE, META_FILE, GENCODE_FILE,
  CONTRACT_FILE, FIXED_ROWS_FILE, BULK_ID_FILE, BULK_FILE, MSIGDB_RDS
)
if (any(!file.exists(required))) {
  stop("Missing proteomics input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

# Deterministic GENCODE v49 primary-assembly map (chr1-22, X, Y, M). No HGNChelper.
primary_chr <- paste0("chr", c(1:22, "X", "Y", "M"))
gc_all <- fread(GENCODE_FILE, select = c("gene_id", "gene_name", "chromosome", "gene_biotype", "ensembl_base"))
gc_prim <- gc_all[chromosome %in% primary_chr & !is.na(gene_name) & nzchar(trimws(gene_name))]
gc_prim[, gene_name := trimws(gene_name)]
gc_sym_counts <- gc_prim[, .(
  n_primary_ids = uniqueN(ensembl_base),
  primary_ensembl_base = paste(sort(unique(ensembl_base)), collapse = ";"),
  primary_gene_id = paste(sort(unique(gene_id)), collapse = ";")
), by = gene_name]

registry <- fread(REGISTRY_FILE)
membership <- fread(MEMBERSHIP_FILE)
contract <- fread(CONTRACT_FILE)
row_contract <- contract[role == "panel4c_rows" & relative_path == sub(paste0("^", BASE, "/?"), "", FIXED_ROWS_FILE)]
if (nrow(row_contract) != 1L) stop("Panel 4C row contract is not uniquely pinned")
observed_row_sha256 <- digest(FIXED_ROWS_FILE, algo = "sha256", file = TRUE, serialize = FALSE)
if (observed_row_sha256 != row_contract$expected_sha256) {
  stop("Panel 4C fixed-row contract drifted from config/fig4_input_contract.tsv")
}
fixed <- fread(FIXED_ROWS_FILE)
setorder(fixed, prog_order, gene_order)
if (nrow(fixed) != 25L || anyDuplicated(fixed$gene)) {
  stop("Protected Panel 4C row contract must contain 25 unique proteins")
}
fwrite(fixed, file.path(OUT, "fixed_panel4c_rows.tsv"), sep = "\t", quote = FALSE)
fwrite(
  data.table(
    relative_path = row_contract$relative_path,
    observed_sha256 = observed_row_sha256,
    expected_sha256 = row_contract$expected_sha256,
    contract_match = TRUE
  ),
  file.path(OUT, "panel4c_contract_audit.tsv"),
  sep = "\t",
  quote = FALSE
)

message("[protein_20261118] loading PXD051911")
raw_all <- fread(RAW_FILE)
meta <- fread(META_FILE)
sample_cols <- setdiff(
  names(raw_all), c("ProteinAccessions", "Genes", "ProteinDescriptions")
)
meta <- meta[
  !is.na(liver_proteomics_filename) & liver_proteomics_filename != "" &
    liver_proteomics_filename %in% sample_cols
]
meta <- unique(meta, by = "liver_proteomics_filename")
if (nrow(meta) != 58L) stop("Expected 58 liver DIA-MS samples; observed ", nrow(meta))

meta[, sample_id := liver_proteomics_filename]
meta[, group := factor(
  fifelse(saf_diagnosis == "No_MASLD", "Control", "MASLD"),
  levels = c("Control", "MASLD")
)]
meta[, acquisition_batch := factor(
  fifelse(grepl("^2019", sample_id), "2019", "2020"),
  levels = c("2019", "2020")
)]
meta[, `:=`(
  age = as.numeric(alder),
  bmi_num = as.numeric(bmi),
  sex = factor(gender),
  Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)),
  Steatosis = as.numeric(steatosis_score),
  Ballooning = as.numeric(hepatocellular_ballooning_score),
  Inflammation = as.numeric(lobular_inflammation_score),
  NAS = as.numeric(nafld_activity_score)
)]
meta[, `:=`(
  age_z = as.numeric(scale(age)),
  bmi_z = as.numeric(scale(bmi_num))
)]
meta[, sample_order := match(sample_id, sample_cols)]
setorder(meta, sample_order)
samples <- meta$sample_id
masld_samples <- meta[group == "MASLD", sample_id]
if (length(masld_samples) != 46L) stop("Expected 46 MASLD samples; observed ", length(masld_samples))

# Explicit PG accounting and exclusion of missing/shared `;`-delimited PGs.
raw_all[, genes_raw := trimws(as.character(Genes))]
raw_all[, pg_missing_gene := is.na(genes_raw) | genes_raw == ""]
raw_all[, pg_shared_semicolon := !pg_missing_gene & grepl(";", genes_raw, fixed = TRUE)]
raw_all[, pg_unique_single_gene := !pg_missing_gene & !pg_shared_semicolon]

excluded_pg <- raw_all[pg_unique_single_gene == FALSE, .(
  ProteinAccessions,
  Genes = genes_raw,
  exclusion_reason = fifelse(pg_missing_gene == TRUE, "missing_gene_symbol", "shared_semicolon_multi_gene_pg")
)]
fwrite(excluded_pg, file.path(OUT, "excluded_shared_pg_groups.tsv"), sep = "\t", quote = FALSE)

raw <- raw_all[pg_unique_single_gene == TRUE]
raw_mat <- as.matrix(raw[, ..samples])
storage.mode(raw_mat) <- "double"
raw_mat[!is.finite(raw_mat) | raw_mat <= 0] <- NA_real_
raw_log2 <- log2(raw_mat)
rownames(raw_log2) <- raw$genes_raw

collapse_by_median <- function(mat) {
  collapse_gene <- function(idx) {
    z <- mat[idx, , drop = FALSE]
    if (nrow(z) == 1L) return(as.numeric(z[1, ]))
    apply(z, 2, median, na.rm = TRUE)
  }
  idx_by_gene <- split(seq_len(nrow(mat)), rownames(mat))
  out <- t(vapply(idx_by_gene, collapse_gene, numeric(ncol(mat))))
  colnames(out) <- colnames(mat)
  out[!is.finite(out)] <- NA_real_
  out
}

# ----------------------------------------------------------------------------
# Track A: All 58 samples (12 Control vs 46 MASLD DE + descriptive abundance)
# ----------------------------------------------------------------------------
keep_58 <- rowMeans(is.na(raw_log2)) <= 0.50
expr_58 <- normalizeBetweenArrays(raw_log2[keep_58, , drop = FALSE], method = "quantile")
gene_expr_58 <- collapse_by_median(expr_58)
if (!all(fixed$gene %in% rownames(gene_expr_58))) {
  stop(
    "Protected Panel 4C proteins missing after 58-sample filtering: ",
    paste(setdiff(fixed$gene, rownames(gene_expr_58)), collapse = ", ")
  )
}

# ----------------------------------------------------------------------------
# Track B: Within-MASLD (n=46) normalization FIRST, then 46-MASLD z-scoring
# ----------------------------------------------------------------------------
raw_log2_masld <- raw_log2[, masld_samples, drop = FALSE]
keep_masld <- rowMeans(is.na(raw_log2_masld)) <= 0.50
expr_masld <- normalizeBetweenArrays(raw_log2_masld[keep_masld, , drop = FALSE], method = "quantile")
gene_expr_masld <- collapse_by_median(expr_masld)
z_masld <- t(scale(t(gene_expr_masld)))
z_masld[!is.finite(z_masld)] <- NA_real_

# Audit PG & GENCODE v49 primary-assembly counts clearly.
n_prim_58 <- sum(rownames(gene_expr_58) %in% gc_sym_counts$gene_name)
n_prim_46 <- sum(rownames(gene_expr_masld) %in% gc_sym_counts$gene_name)
pg_audit <- data.table(
  metric = c(
    "n_pg_total",
    "n_pg_missing_gene",
    "n_pg_shared_semicolon_excluded",
    "n_pg_unique_single_gene_raw",
    "n_unique_symbols_raw",
    "n_pg_pass50_58sample",
    "n_unique_proteins_58sample",
    "n_unique_proteins_58sample_gencode_v49_primary",
    "n_unique_proteins_58sample_non_primary_or_alias",
    "n_pg_pass50_46masld",
    "n_unique_proteins_46masld",
    "n_unique_proteins_46masld_gencode_v49_primary",
    "n_unique_proteins_46masld_non_primary_or_alias"
  ),
  value = c(
    nrow(raw_all),
    sum(raw_all$pg_missing_gene),
    sum(raw_all$pg_shared_semicolon),
    sum(raw_all$pg_unique_single_gene),
    uniqueN(raw$genes_raw),
    sum(keep_58),
    nrow(gene_expr_58),
    n_prim_58,
    nrow(gene_expr_58) - n_prim_58,
    sum(keep_masld),
    nrow(gene_expr_masld),
    n_prim_46,
    nrow(gene_expr_masld) - n_prim_46
  )
)
fwrite(pg_audit, file.path(OUT, "pg_mapping_audit.tsv"), sep = "\t", quote = FALSE)

# ----------------------------------------------------------------------------
# DE (12 Control vs 46 MASLD on 58-sample normalized matrix)
# ----------------------------------------------------------------------------
design <- model.matrix(~ group + acquisition_batch + age_z + bmi_z + sex, data = meta)
coef_name <- "groupMASLD"
if (!coef_name %in% colnames(design)) stop("MASLD coefficient missing from design")
fit <- eBayes(lmFit(gene_expr_58[, samples, drop = FALSE], design), robust = TRUE)
tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
de <- as.data.table(tt, keep.rownames = "gene")
setnames(de, c("P.Value", "adj.P.Val"), c("pvalue", "padj"))

batch_effect <- function(batch_id) {
  idx <- meta$acquisition_batch == batch_id
  d <- droplevels(meta[idx])
  rhs <- c("group", "age_z", "bmi_z")
  if (nlevels(d$sex) > 1L) rhs <- c(rhs, "sex")
  dm <- model.matrix(reformulate(rhs), data = d)
  f <- eBayes(lmFit(gene_expr_58[, d$sample_id, drop = FALSE], dm), robust = TRUE)
  z <- as.data.table(topTable(f, coef = "groupMASLD", number = Inf, sort.by = "none"), keep.rownames = "gene")
  z[, .(gene, effect = logFC, pvalue = P.Value)]
}
b19 <- batch_effect("2019")
b20 <- batch_effect("2020")
setnames(b19, c("effect", "pvalue"), c("effect_2019", "pvalue_2019"))
setnames(b20, c("effect", "pvalue"), c("effect_2020", "pvalue_2020"))
interaction_design <- model.matrix(~ group * acquisition_batch + age_z + bmi_z + sex, data = meta)
interaction_fit <- eBayes(lmFit(gene_expr_58[, samples, drop = FALSE], interaction_design), robust = TRUE)
interaction <- as.data.table(
  topTable(interaction_fit, coef = "groupMASLD:acquisition_batch2020", number = Inf, sort.by = "none"),
  keep.rownames = "gene"
)[, .(gene, interaction_effect = logFC, interaction_pvalue = P.Value)]
de <- Reduce(function(x, y) merge(x, y, by = "gene", all.x = TRUE, sort = FALSE), list(de, b19, b20, interaction))
de <- merge(
  de,
  gc_sym_counts[, .(gene = gene_name, n_gencode_v49_primary_ids = n_primary_ids, gencode_v49_primary_ensembl_base = primary_ensembl_base)],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
de[is.na(n_gencode_v49_primary_ids), n_gencode_v49_primary_ids := 0L]
de[, in_gencode_v49_primary := n_gencode_v49_primary_ids >= 1L]
setorder(de, gene)
fwrite(de, file.path(OUT, "protein_de_adjusted.tsv"), sep = "\t", quote = FALSE)

# Batch audit (58 samples)
imp <- gene_expr_58
for (i in seq_len(nrow(imp))) {
  med <- median(imp[i, ], na.rm = TRUE)
  imp[i, is.na(imp[i, ])] <- med
}
pca <- prcomp(t(imp), center = TRUE, scale. = TRUE)
pc1_var <- summary(pca)$importance[2, 1]
batch_audit <- meta[, .(
  n = .N,
  n_masld = sum(group == "MASLD"),
  median_sample_log2 = median(colMeans(gene_expr_58[, sample_id, drop = FALSE], na.rm = TRUE))
), by = acquisition_batch]
batch_audit[, `:=`(
  pc1_variance_fraction = pc1_var,
  pc1_mean = vapply(acquisition_batch, function(b) mean(pca$x[meta$acquisition_batch == b, 1]), numeric(1))
)]
fwrite(batch_audit, file.path(OUT, "batch_audit.tsv"), sep = "\t", quote = FALSE)

# 58-sample descriptive heatmap matrix (retained strictly for display of Panel 4C)
nuisance <- cbind(
  age_z = meta$age_z,
  bmi_z = meta$bmi_z,
  sex_male = as.numeric(meta$sex == "Male")
)
keep_design <- model.matrix(~ group, data = meta)
adjusted_expr <- removeBatchEffect(
  gene_expr_58[, samples, drop = FALSE],
  batch = meta$acquisition_batch,
  covariates = nuisance,
  design = keep_design
)
adjusted_expr <- adjusted_expr[fixed$gene, samples, drop = FALSE]
fwrite(
  as.data.table(adjusted_expr, keep.rownames = "gene"),
  file.path(OUT, "panel4c_adjusted_abundance.tsv"),
  sep = "\t",
  quote = FALSE
)

panel_meta <- meta[, .(
  sample_id, group, acquisition_batch, sex, age, bmi = bmi_num,
  Steatosis, Ballooning, Inflammation, Fibrosis, NAS
)]
fwrite(panel_meta, file.path(OUT, "panel4c_metadata.tsv"), sep = "\t", quote = FALSE)

# ----------------------------------------------------------------------------
# Bulk DEG comparison (774-unit primary_source_controls fit)
# ----------------------------------------------------------------------------
bulk_header <- names(fread(BULK_FILE, nrows = 0L))
pick_bulk_col <- function(role, candidates) {
  hit <- intersect(candidates, bulk_header)
  if (length(hit) != 1L) {
    stop(
      "Bulk DEG file ", BULK_FILE, " needs exactly one ", role, " column from {",
      paste(candidates, collapse = ", "), "}; found {", paste(hit, collapse = ", "), "}"
    )
  }
  hit
}
bulk_cols <- c(
  gene_id = pick_bulk_col("Ensembl gene ID", c("gene", "gene_id", "ensembl_id")),
  symbol = pick_bulk_col("gene symbol", c("symbol", "gene_symbol")),
  logFC = pick_bulk_col("log2 fold change", "logFC"),
  padj = pick_bulk_col("conventional BH-adjusted p-value", c("padj", "adj.P.Val")),
  treat_fdr = pick_bulk_col("TREAT FDR", "treat_fdr")
)
bulk <- fread(
  BULK_FILE,
  select = unname(bulk_cols),
  colClasses = list(character = unname(bulk_cols[c("gene_id", "symbol")]))
)
setnames(bulk, unname(bulk_cols), names(bulk_cols))
bulk[, gene_id_unversioned := sub("\\.[0-9]+$", "", gene_id)]

# Map fixed Panel 4C symbols deterministically against primary-assembly GENCODE v49 + bulk fit
bulk_ids <- fread(BULK_ID_FILE)
if (nrow(bulk_ids) != 25L || anyDuplicated(bulk_ids$gene) ||
    anyDuplicated(bulk_ids$bulk_gene_id) || !setequal(bulk_ids$gene, fixed$gene)) {
  stop("Panel 4C bulk gene-ID table must map the 25 fixed proteins to 25 unique Ensembl IDs")
}
bulk_ids[, pinned_id_unversioned := sub("\\.[0-9]+$", "", bulk_gene_id)]
res_list <- lapply(seq_len(nrow(bulk_ids)), function(idx) {
  g <- bulk_ids$gene[idx]
  pid_base <- bulk_ids$pinned_id_unversioned[idx]
  hit <- bulk[gene_id_unversioned == pid_base]
  if (nrow(hit) == 1L && !is.na(hit$symbol) && hit$symbol == g) {
    return(data.table(gene = g, resolved_id_unversioned = pid_base, resolution_mode = "pinned_ensembl_id"))
  }
  # Fallback strictly deterministic: unique GENCODE v49 primary-assembly ID matching symbol g & bulk fit
  gc_hit <- unique(gc_prim[gene_name == g, ensembl_base])
  b_hit <- bulk[gene_id_unversioned %in% gc_hit & symbol == g]
  if (nrow(b_hit) != 1L) {
    stop("Cannot uniquely resolve primary-assembly bulk row for Panel 4C gene ", g)
  }
  data.table(
    gene = g,
    resolved_id_unversioned = b_hit$gene_id_unversioned[1L],
    resolution_mode = paste0("gencode_v49_primary_resolved_from_", pid_base)
  )
})
bulk_resolved_map <- rbindlist(res_list)
fwrite(bulk_resolved_map, file.path(OUT, "panel4c_bulk_id_resolution.tsv"), sep = "\t", quote = FALSE)

bulk_sel <- merge(
  bulk_resolved_map[, .(gene, gene_id_unversioned = resolved_id_unversioned)],
  bulk,
  by = "gene_id_unversioned",
  all.x = TRUE,
  sort = FALSE
)
if (nrow(bulk_sel) != 25L || any(is.na(bulk_sel$gene_id)) || any(bulk_sel$symbol != bulk_sel$gene)) {
  stop("Bulk resolution failed for Panel 4C fixed rows")
}

panel_pairs <- merge(
  fixed[, .(gene, program, prog_order, gene_order)],
  de[, .(
    gene,
    protein_logFC = logFC,
    protein_t = t,
    protein_padj = padj,
    protein_effect_2019 = effect_2019,
    protein_effect_2020 = effect_2020,
    protein_batch_interaction_p = interaction_pvalue
  )],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
panel_pairs <- merge(
  panel_pairs,
  bulk_sel[, .(
    gene,
    bulk_gene_id = gene_id,
    bulk_logFC = logFC,
    bulk_padj = padj,
    bulk_treat_fdr = treat_fdr
  )],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
setorder(panel_pairs, prog_order, gene_order)
panel_pairs[, `:=`(
  bulk_sig = is.finite(bulk_padj) & bulk_padj < 0.05 & abs(bulk_logFC) > 0.5,
  protein_sig = is.finite(protein_padj) & protein_padj < 0.05,
  direction_concordant = sign(bulk_logFC) == sign(protein_logFC),
  selection_conditioned = TRUE,
  interpretation = "descriptive_same_cohort_reestimate"
)]
fwrite(panel_pairs, file.path(OUT, "panel4c_mrna_protein.tsv"), sep = "\t", quote = FALSE)

# Panel 4C partial Spearman histology associations (emit both 58-sample and 46-MASLD-only normalized)
run_panel4c_histology <- function(expr_mat, sub_meta, tag_cohort) {
  d_m <- copy(sub_meta)
  d_m[, `:=`(age_z = as.numeric(scale(age)), bmi_z = as.numeric(scale(bmi_num)))]
  samp_sub <- d_m$sample_id
  nuisance_vars <- c("acquisition_batch", "age_z", "bmi_z", "sex")
  covariates_ok <- complete.cases(d_m[, ..nuisance_vars])
  features <- c("Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS")
  res_h <- rbindlist(lapply(fixed$gene, function(g) {
    rbindlist(lapply(features, function(feature) {
      protein <- expr_mat[g, samp_sub]
      score <- d_m[[feature]]
      ok <- is.finite(protein) & is.finite(score) & covariates_ok
      d <- droplevels(d_m[ok])
      rhs <- c("age_z", "bmi_z")
      if (nlevels(d$acquisition_batch) > 1L) rhs <- c("acquisition_batch", rhs)
      if (nlevels(d$sex) > 1L) rhs <- c(rhs, "sex")
      x <- model.matrix(reformulate(rhs), data = d)
      design_rank <- qr(x)$rank
      out <- data.table(
        gene = g, feature = feature, n = sum(ok), rho = NA_real_, pvalue = NA_real_,
        n_nuisance_terms = design_rank - 1L, df = sum(ok) - design_rank - 1L,
        cohort = tag_cohort, test_status = "tested"
      )
      if (design_rank < ncol(x)) return(out[, test_status := "rank_deficient_nuisance_design"])
      if (out$n < 8L || out$df < 3L) return(out[, test_status := "too_few_participants"])
      protein_resid <- residuals(lm.fit(x, rank(protein[ok], ties.method = "average")))
      feature_resid <- residuals(lm.fit(x, rank(score[ok], ties.method = "average")))
      r <- suppressWarnings(cor(protein_resid, feature_resid))
      if (!is.finite(r)) return(out[, test_status := "no_residual_variance"])
      t_stat <- r * sqrt(out$df / (1 - r^2))
      out[, `:=`(rho = r, pvalue = 2 * pt(-abs(t_stat), out$df))]
    }))
  }))
  res_h[, padj := p.adjust(pvalue, method = "BH")]
  res_h[, correction_family := "25_fixed_proteins_x_5_histology_features"]
  setcolorder(res_h, c("gene", "feature", "n", "rho", "pvalue", "padj"))
  res_h
}
histology_58 <- run_panel4c_histology(gene_expr_58, meta, "all_58")
fwrite(histology_58, file.path(OUT, "panel4c_histology_partial.tsv"), sep = "\t", quote = FALSE)
histology_masld <- run_panel4c_histology(gene_expr_masld, meta[group == "MASLD"], "MASLD_46_only_qnorm")
fwrite(histology_masld, file.path(OUT, "panel4c_histology_partial_masld.tsv"), sep = "\t", quote = FALSE)

# Pathway enrichment on adjusted t-statistic (58-sample DE)
message("[protein_20261118] loading frozen MSigDB 2025.1 cache")
mdb <- as.data.table(readRDS(MSIGDB_RDS))
reactome <- mdb[gs_collection == "C2" & gs_subcollection == "CP:REACTOME", .(gs_name, gene_symbol = db_gene_symbol)]
hallmark <- mdb[gs_collection == "H", .(gs_name, gene_symbol = db_gene_symbol)]
msig <- rbindlist(list(
  reactome[, .(pathway = gs_name, gene = gene_symbol, family = "Reactome")],
  hallmark[, .(pathway = gs_name, gene = gene_symbol, family = "Hallmark")]
))
Rset <- function(name) unique(reactome[gs_name == name, gene_symbol])
lipid_droplet <- c(
  "PLIN1", "PLIN2", "PLIN3", "PLIN4", "PLIN5", "ABHD5", "CIDEA", "CIDEB",
  "CIDEC", "G0S2", "HILPDA", "ACLY", "FASN", "ACACA", "SCD", "DGAT1", "DGAT2", "LPIN1"
)
immune <- c(
  "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
  "CD74", "CD68", "CD163", "AIF1", "TYROBP", "FCER1G", "LGALS3", "CTSS",
  "LAPTM5", "MARCO", "VSIG4", "CAPG", "SPP1", "CASP1", "PYCARD", "S100A8",
  "S100A9", "PTPRC", "LCP1", "CORO1A", "ITGAL", "C1QA", "C1QB", "C1QC", "C3"
)
curated <- list(
  "CURATED_LIPID_DROPLET" = lipid_droplet,
  "CURATED_IMMUNE_INFLAMMATION" = immune,
  "CURATED_DETOX_PHASE_I_II" = unique(c(
    Rset("REACTOME_PHASE_I_FUNCTIONALIZATION_OF_COMPOUNDS"),
    Rset("REACTOME_PHASE_II_CONJUGATION_OF_COMPOUNDS")
  )),
  "CURATED_MITO_OXPHOS_TCA" = unique(c(
    Rset("REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"),
    Rset("REACTOME_COMPLEX_I_BIOGENESIS"),
    Rset("REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE")
  ))
)
msig <- rbind(msig, rbindlist(lapply(names(curated), function(p) data.table(pathway = p, gene = curated[[p]], family = "Curated"))))
pathways <- split(msig$gene, msig$pathway)
ranked <- de[is.finite(t) & !duplicated(gene)]
ranks <- sort(setNames(ranked$t, ranked$gene), decreasing = TRUE)
pathways <- lapply(pathways, intersect, y = names(ranks))
fg <- as.data.table(fgsea(
  pathways = pathways, stats = ranks, minSize = 8, maxSize = 600,
  eps = 0, nproc = 1, BPPARAM = BiocParallel::SerialParam()
))
fg[, family := msig$family[match(pathway, msig$pathway)]]
fg[, abs_nes := abs(NES)]
setorder(fg, padj, -abs_nes)
fwrite(fg[, .(pathway, family, size, NES, pval, padj, leadingEdge)], file.path(OUT, "protein_pathway_enrichment.tsv"), sep = "\t", quote = FALSE)

process_map <- data.table(
  program = c("Lipid droplet", "ECM / BM", "Immune / inflammation", "Detox", "AA catabolism", "Mito / OXPHOS"),
  pathway = c(
    "CURATED_LIPID_DROPLET", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
    "CURATED_IMMUNE_INFLAMMATION", "CURATED_DETOX_PHASE_I_II",
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES", "CURATED_MITO_OXPHOS_TCA"
  ),
  representable = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  display_order = 1:6
)
process_enrichment <- merge(process_map, fg, by = "pathway", all.x = TRUE, sort = FALSE)
setorder(process_enrichment, display_order)
fwrite(process_enrichment, file.path(OUT, "panel4c_process_enrichment.tsv"), sep = "\t", quote = FALSE)

# ----------------------------------------------------------------------------
# 22-module projection:
# - Case-vs-Control (12 Ctrl vs 46 MASLD) uses z_58
# - Within-MASLD (n=46) scores and histology use z_masld (46-MASLD qnorm & scale!)
# ----------------------------------------------------------------------------
z_58 <- t(scale(t(gene_expr_58[, samples, drop = FALSE])))
z_58[!is.finite(z_58)] <- NA_real_

score_program_matrix <- function(z_mat, samp_ids, meta_sub, m_df) {
  m_sub <- m_df[, .(original_l1_weight = sum(original_l1_weight)), by = gene_symbol]
  measured <- m_sub[gene_symbol %in% rownames(z_mat)]
  n_measured <- nrow(measured)
  retained_weight <- sum(measured$original_l1_weight)
  testable <- n_measured >= 8L && retained_weight >= 0.20
  if (!testable) return(list(n_measured = n_measured, retained_weight = retained_weight, testable = FALSE, dt = NULL))
  measured[, weight := original_l1_weight / sum(original_l1_weight)]
  zz <- z_mat[measured$gene_symbol, samp_ids, drop = FALSE]
  obs_w <- colSums(is.finite(zz) * measured$weight)
  w_score <- colSums(zz * measured$weight, na.rm = TRUE) / obs_w
  w_score[!is.finite(w_score)] <- NA_real_
  eq_score <- colMeans(zz, na.rm = TRUE)
  eq_score[!is.finite(eq_score)] <- NA_real_
  top_gene <- measured$gene_symbol[which.max(measured$weight)]
  leave <- measured[gene_symbol != top_gene]
  leave[, weight := original_l1_weight / sum(original_l1_weight)]
  leave_zz <- z_mat[leave$gene_symbol, samp_ids, drop = FALSE]
  leave_w <- colSums(is.finite(leave_zz) * leave$weight)
  l_score <- colSums(leave_zz * leave$weight, na.rm = TRUE) / leave_w
  l_score[!is.finite(l_score)] <- NA_real_
  list(
    n_measured = n_measured,
    retained_weight = retained_weight,
    testable = TRUE,
    dt = data.table(
      sample_id = samp_ids,
      group = as.character(meta_sub$group),
      acquisition_batch = as.character(meta_sub$acquisition_batch),
      score = w_score,
      equal_score = eq_score,
      leave_top_score = l_score
    )
  )
}

fit_score <- function(score, subset = rep(TRUE, length(score))) {
  d <- copy(meta[subset])
  d[, score := score[subset]]
  d <- d[is.finite(score)]
  d[, `:=`(group = droplevels(group), acquisition_batch = droplevels(acquisition_batch), sex = droplevels(sex))]
  if (nrow(d) <= 8L || length(unique(d$group)) < 2L) return(c(effect = NA_real_, se = NA_real_, pvalue = NA_real_))
  rhs <- c("group", "age_z", "bmi_z")
  if (nlevels(d$acquisition_batch) > 1L) rhs <- c(rhs, "acquisition_batch")
  if (nlevels(d$sex) > 1L) rhs <- c(rhs, "sex")
  mod <- lm(reformulate(rhs, response = "score"), data = d)
  co <- summary(mod)$coefficients
  if (!"groupMASLD" %in% rownames(co)) return(c(effect = NA_real_, se = NA_real_, pvalue = NA_real_))
  c(effect = co["groupMASLD", 1], se = co["groupMASLD", 2], pvalue = co["groupMASLD", 4])
}

meta_masld_sub <- copy(meta[group == "MASLD"])
meta_masld_sub[, sample_order := match(sample_id, masld_samples)]
setorder(meta_masld_sub, sample_order)
meta_masld_sub[, `:=`(age_z = as.numeric(scale(age)), bmi_z = as.numeric(scale(bmi_num)))]

score_rows_58 <- list()
score_rows_masld <- list()
score_results <- list()

for (i in seq_len(nrow(registry))) {
  pid <- registry$program_id[i]
  m_df <- membership[program_id == pid & mapped_symbol == TRUE & !is.na(gene_symbol)]
  # Restrict to deterministic primary-assembly confirmed symbols too (we checked v1 has 0 non-primary, still explicit)
  sc58 <- score_program_matrix(z_58, samples, meta, m_df)
  sc46 <- score_program_matrix(z_masld, masld_samples, meta_masld_sub, m_df)
  base <- data.table(
    program_id = pid,
    n_measured = sc58$n_measured,
    retained_l1_weight = sc58$retained_weight,
    testable = sc58$testable,
    n_measured_masld46 = sc46$n_measured,
    retained_l1_weight_masld46 = sc46$retained_weight,
    testable_masld46 = sc46$testable
  )
  if (!sc58$testable) {
    score_results[[pid]] <- base[, `:=`(
      effect = NA_real_, se = NA_real_, pvalue = NA_real_,
      equal_effect = NA_real_, leave_top_effect = NA_real_,
      effect_2019 = NA_real_, effect_2020 = NA_real_, interaction_p = NA_real_
    )]
  } else {
    primary <- fit_score(sc58$dt$score)
    equal <- fit_score(sc58$dt$equal_score)
    leave_top <- fit_score(sc58$dt$leave_top_score)
    b19 <- fit_score(sc58$dt$score, meta$acquisition_batch == "2019")
    b20 <- fit_score(sc58$dt$score, meta$acquisition_batch == "2020")
    int_dat <- copy(meta); int_dat[, score := sc58$dt$score]
    int_fit <- lm(score ~ group * acquisition_batch + age_z + bmi_z + sex, data = int_dat)
    int_name <- grep("groupMASLD:acquisition_batch", names(coef(int_fit)), value = TRUE)
    int_p <- if (length(int_name) == 1L) summary(int_fit)$coefficients[int_name, 4] else NA_real_
    score_results[[pid]] <- base[, `:=`(
      effect = unname(primary["effect"]), se = unname(primary["se"]), pvalue = unname(primary["pvalue"]),
      equal_effect = unname(equal["effect"]), leave_top_effect = unname(leave_top["effect"]),
      effect_2019 = unname(b19["effect"]), effect_2020 = unname(b20["effect"]), interaction_p = int_p
    )]
    score_rows_58[[pid]] <- cbind(program_id = pid, sc58$dt)
  }
  if (sc46$testable) {
    score_rows_masld[[pid]] <- cbind(program_id = pid, sc46$dt)
  }
}

module_results <- rbindlist(score_results, fill = TRUE)
module_results[testable == TRUE, padj := p.adjust(pvalue, method = "BH")]
module_results[, sensitivity_sign_agree := testable & is.finite(effect) &
  sign(effect) == sign(equal_effect) & sign(effect) == sign(leave_top_effect)]
module_results[, robust := testable & is.finite(padj) & padj < 0.05 & sensitivity_sign_agree]
module_results <- merge(registry[, .(program_id, display_order, cell_type, module, program_name)], module_results, by = "program_id", all.x = TRUE, sort = FALSE)
setorder(module_results, display_order)
fwrite(module_results, file.path(OUT, "module_protein_results.tsv"), sep = "\t", quote = FALSE)

if (length(score_rows_58)) {
  fwrite(rbindlist(score_rows_58), file.path(OUT, "module_protein_scores_58.tsv"), sep = "\t", quote = FALSE)
}
if (length(score_rows_masld)) {
  module_scores_masld <- rbindlist(score_rows_masld)
  fwrite(module_scores_masld, file.path(OUT, "module_protein_scores_masld.tsv"), sep = "\t", quote = FALSE)

  severity_features <- c("Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS")
  severity_data <- merge(
    module_scores_masld,
    meta_masld_sub[, c("sample_id", "sex", "age_z", "bmi_z", severity_features), with = FALSE],
    by = "sample_id", all.x = TRUE, sort = FALSE
  )
  rank_residualize_severity <- function(y, d, include_batch = TRUE) {
    d <- droplevels(copy(d))
    rhs <- c("age_z", "bmi_z")
    if (include_batch && uniqueN(d$acquisition_batch) > 1L) rhs <- c("acquisition_batch", rhs)
    if (uniqueN(d$sex[!is.na(d$sex)]) > 1L) rhs <- c(rhs, "sex")
    nuisance <- model.matrix(reformulate(rhs), data = d)
    ok <- is.finite(y) & apply(nuisance, 1L, function(z) all(is.finite(z)))
    out <- rep(NA_real_, length(y))
    if (sum(ok) > ncol(nuisance) + 2L) {
      out[ok] <- residuals(lm.fit(nuisance[ok, , drop = FALSE], rank(y[ok], ties.method = "average")))
    }
    out
  }
  severity_cor <- function(d, score_col, feature, include_batch = TRUE) {
    score_resid <- rank_residualize_severity(d[[score_col]], d, include_batch)
    hist_resid <- rank_residualize_severity(d[[feature]], d, include_batch)
    ok <- is.finite(score_resid) & is.finite(hist_resid)
    if (sum(ok) < 8L) return(c(n = sum(ok), rho = NA_real_, pvalue = NA_real_))
    tst <- suppressWarnings(cor.test(score_resid[ok], hist_resid[ok], method = "pearson"))
    c(n = sum(ok), rho = unname(tst$estimate), pvalue = tst$p.value)
  }
  severity_rows <- rbindlist(lapply(unique(module_scores_masld$program_id), function(pid) {
    d <- severity_data[program_id == pid]
    rbindlist(lapply(severity_features, function(feature) {
      primary <- severity_cor(d, "score", feature)
      equal <- severity_cor(d, "equal_score", feature)
      leave_top <- severity_cor(d, "leave_top_score", feature)
      batch_rho <- vapply(c("2019", "2020"), function(bid) {
        unname(severity_cor(d[acquisition_batch == bid], "score", feature, include_batch = FALSE)["rho"])
      }, numeric(1))
      data.table(
        program_id = pid, feature = feature, n = as.integer(primary["n"]),
        rho = unname(primary["rho"]), pvalue = unname(primary["pvalue"]),
        equal_rho = unname(equal["rho"]), leave_top_rho = unname(leave_top["rho"]),
        rho_2019 = batch_rho[[1]], rho_2020 = batch_rho[[2]]
      )
    }))
  }))
  severity_template <- merge(
    CJ(program_id = registry$program_id, feature = severity_features, unique = TRUE, sorted = FALSE),
    module_results[, .(program_id, display_order, cell_type, module, program_name,
                       n_measured = n_measured_masld46, retained_l1_weight = retained_l1_weight_masld46,
                       testable = testable_masld46)],
    by = "program_id", all.x = TRUE, sort = FALSE
  )
  module_histology <- merge(severity_template, severity_rows, by = c("program_id", "feature"), all.x = TRUE, sort = FALSE)
  module_histology[, feature_order := match(feature, severity_features)]
  module_histology[testable == TRUE & is.finite(pvalue), padj := p.adjust(pvalue, method = "BH")]
  module_histology[, sensitivity_sign_agree := testable & is.finite(rho) &
    is.finite(equal_rho) & is.finite(leave_top_rho) &
    sign(rho) == sign(equal_rho) & sign(rho) == sign(leave_top_rho)]
  module_histology[, batch_sign_agree := testable & is.finite(rho) &
    is.finite(rho_2019) & is.finite(rho_2020) &
    sign(rho) == sign(rho_2019) & sign(rho) == sign(rho_2020)]
  module_histology[, robust := testable & is.finite(padj) & padj < 0.05 & sensitivity_sign_agree & batch_sign_agree]
  module_histology[, `:=`(
    cohort = "MASLD_only_46_qnorm",
    correction_scope = paste0(sum(module_results$testable_masld46), "_testable_programs_x_5_histology_features"),
    interpretation = "independent_frozen_program_severity_association"
  )]
  setorder(module_histology, display_order, feature_order)
  fwrite(module_histology, file.path(OUT, "module_protein_histology_masld.tsv"), sep = "\t", quote = FALSE)
}

cat("[protein_20261118] PG accounting:\n")
print(pg_audit)
cat("[protein_20261118] fixed Panel 4C direction concordant:",
    sum(panel_pairs$direction_concordant, na.rm = TRUE), "/25\n")
cat("[protein_20261118] fixed Panel 4C protein padj<0.05:",
    sum(panel_pairs$protein_sig, na.rm = TRUE), "/25\n")
cat("[protein_20261118] 22-module DE testable:", sum(module_results$testable), "/22; robust:",
    sum(module_results$robust), "\n")
