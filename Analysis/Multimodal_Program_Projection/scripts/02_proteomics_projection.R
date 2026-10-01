#!/usr/bin/env Rscript

# PXD051911 liver DIA-MS projection for Figure 4.
#
# The assay-native 25-protein display was selected in an earlier exploratory
# analysis of PXD051911 and is now held fixed as a descriptive row contract.
# Protein-group intensities are log2 transformed and quantile normalized, and the primary
# MASLD-vs-control model adjusts for acquisition batch, age, BMI, and sex. The
# same normalized matrix is then used to project all 22 frozen Hotspot modules.

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
OUT <- Sys.getenv("PROTEOMICS_OUT_DIR", file.path(ROOT, "results/proteomics"))
# Never overwrite a previous Panel 4C build; a rerun writes to a new directory.
existing_outputs <- file.path(OUT, c("panel4c_mrna_protein.tsv", "panel4c_histology_partial.tsv"))
if (any(file.exists(existing_outputs))) {
  stop(
    "Refusing to overwrite existing output(s): ",
    paste(existing_outputs[file.exists(existing_outputs)], collapse = ", "),
    ". Set PROTEOMICS_OUT_DIR to a new directory."
  )
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

REGISTRY_FILE <- file.path(ROOT, "results/frozen_programs.tsv")
MEMBERSHIP_FILE <- file.path(ROOT, "results/frozen_program_membership.tsv")
RAW_FILE <- file.path(BASE, "data/PXD051911/liver_protein_quant.txt")
META_FILE <- file.path(BASE, "data/PXD051911/meta_data.txt")
CONTRACT_FILE <- file.path(ROOT, "config/fig4_input_contract.tsv")
FIXED_ROWS_FILE <- file.path(ROOT, "config/panel4c_fixed_rows.tsv")
# Stable Ensembl ID of the bulk row for each fixed display gene. HLA-DQA1 has
# four bulk rows; the display row is ENSG00000206305.
BULK_ID_FILE <- file.path(ROOT, "config/panel4c_bulk_gene_ids.tsv")
BULK_FILE <- Sys.getenv(
  "BULK_DEG_FILE",
  file.path(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"
  )
)
MSIGDB_RDS <- Sys.getenv(
  "MSIGDB_RDS",
  "/gpfs/commons/home/jameslee/.cache/R/msigdbr/msigdb.2025.1.Hs.rds"
)

required <- c(
  REGISTRY_FILE, MEMBERSHIP_FILE, RAW_FILE, META_FILE, CONTRACT_FILE, FIXED_ROWS_FILE,
  BULK_ID_FILE, BULK_FILE, MSIGDB_RDS
)
if (any(!file.exists(required))) {
  stop("Missing proteomics input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

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

message("[protein] loading PXD051911")
raw <- fread(RAW_FILE)
meta <- fread(META_FILE)
sample_cols <- setdiff(
  names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions")
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

# Restrict to unambiguous single-gene protein groups. Quantile normalization is
# performed before duplicate symbols are collapsed and without disease labels.
raw <- raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
expr <- as.matrix(raw[, ..samples])
storage.mode(expr) <- "double"
expr[!is.finite(expr) | expr <= 0] <- NA_real_
expr <- log2(expr)
rownames(expr) <- raw$Genes
keep <- rowMeans(is.na(expr)) <= 0.50
expr <- expr[keep, , drop = FALSE]
expr <- normalizeBetweenArrays(expr, method = "quantile")

collapse_gene <- function(idx) {
  z <- expr[idx, , drop = FALSE]
  if (nrow(z) == 1L) return(as.numeric(z[1, ]))
  apply(z, 2, median, na.rm = TRUE)
}
idx_by_gene <- split(seq_len(nrow(expr)), rownames(expr))
gene_expr <- t(vapply(idx_by_gene, collapse_gene, numeric(ncol(expr))))
colnames(gene_expr) <- colnames(expr)
gene_expr[!is.finite(gene_expr)] <- NA_real_
if (!all(fixed$gene %in% rownames(gene_expr))) {
  stop(
    "Protected Panel 4C proteins missing after source-compatible filtering: ",
    paste(setdiff(fixed$gene, rownames(gene_expr)), collapse = ", ")
  )
}

design <- model.matrix(
  ~ group + acquisition_batch + age_z + bmi_z + sex,
  data = meta
)
coef_name <- "groupMASLD"
if (!coef_name %in% colnames(design)) stop("MASLD coefficient missing from design")
fit <- eBayes(lmFit(gene_expr[, samples, drop = FALSE], design), robust = TRUE)
tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
de <- as.data.table(tt, keep.rownames = "gene")
setnames(de, c("P.Value", "adj.P.Val"), c("pvalue", "padj"))

# Cheapest refutation of residual batch confounding: within-batch adjusted
# effects plus the explicit disease-by-batch interaction. These are sensitivity
# columns; the prespecified main effect above remains primary.
batch_effect <- function(batch_id) {
  idx <- meta$acquisition_batch == batch_id
  d <- droplevels(meta[idx])
  rhs <- c("group", "age_z", "bmi_z")
  if (nlevels(d$sex) > 1L) rhs <- c(rhs, "sex")
  dm <- model.matrix(reformulate(rhs), data = d)
  f <- eBayes(lmFit(gene_expr[, d$sample_id, drop = FALSE], dm), robust = TRUE)
  z <- as.data.table(topTable(f, coef = "groupMASLD", number = Inf, sort.by = "none"), keep.rownames = "gene")
  z[, .(gene, effect = logFC, pvalue = P.Value)]
}
b19 <- batch_effect("2019")
b20 <- batch_effect("2020")
setnames(b19, c("effect", "pvalue"), c("effect_2019", "pvalue_2019"))
setnames(b20, c("effect", "pvalue"), c("effect_2020", "pvalue_2020"))
interaction_design <- model.matrix(
  ~ group * acquisition_batch + age_z + bmi_z + sex,
  data = meta
)
interaction_name <- "groupMASLD:acquisition_batch2020"
interaction_fit <- eBayes(
  lmFit(gene_expr[, samples, drop = FALSE], interaction_design),
  robust = TRUE
)
interaction <- as.data.table(
  topTable(interaction_fit, coef = interaction_name, number = Inf, sort.by = "none"),
  keep.rownames = "gene"
)[, .(gene, interaction_effect = logFC, interaction_pvalue = P.Value)]
de <- Reduce(function(x, y) merge(x, y, by = "gene", all.x = TRUE, sort = FALSE), list(de, b19, b20, interaction))
setorder(de, gene)
fwrite(de, file.path(OUT, "protein_de_adjusted.tsv"), sep = "\t", quote = FALSE)

# Batch audit: this is a diagnostic, not an adjustment surrogate.
imp <- gene_expr
for (i in seq_len(nrow(imp))) {
  med <- median(imp[i, ], na.rm = TRUE)
  imp[i, is.na(imp[i, ])] <- med
}
pca <- prcomp(t(imp), center = TRUE, scale. = TRUE)
pc1_var <- summary(pca)$importance[2, 1]
batch_audit <- meta[, .(
  n = .N,
  n_masld = sum(group == "MASLD"),
  median_sample_log2 = median(colMeans(gene_expr[, sample_id, drop = FALSE], na.rm = TRUE))
), by = acquisition_batch]
batch_audit[, `:=`(
  pc1_variance_fraction = pc1_var,
  pc1_mean = vapply(
    acquisition_batch,
    function(b) mean(pca$x[meta$acquisition_batch == b, 1]),
    numeric(1)
  )
)]
fwrite(batch_audit, file.path(OUT, "batch_audit.tsv"), sep = "\t", quote = FALSE)

# Remove nuisance effects while explicitly retaining the disease contrast for
# the patient heatmap. This matrix is display-only; inference comes from the
# full limma model above.
nuisance <- cbind(
  age_z = meta$age_z,
  bmi_z = meta$bmi_z,
  sex_male = as.numeric(meta$sex == "Male")
)
keep_design <- model.matrix(~ group, data = meta)
adjusted_expr <- removeBatchEffect(
  gene_expr[, samples, drop = FALSE],
  batch = meta$acquisition_batch,
  covariates = nuisance,
  design = keep_design
)
adjusted_expr <- adjusted_expr[fixed$gene, , drop = FALSE]
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

# The bulk source is either the live canonical table or a release
# deg_results.csv. Resolve each column role explicitly and stop if any is
# missing or ambiguous. The conventional BH p-value (padj) and the TREAT FDR
# are carried as separate fields; TREAT is the sensitivity arm only.
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

# Select the bulk row for each fixed display gene by its stable Ensembl ID
# (version-insensitive), never by the lowest FDR among rows sharing a symbol.
bulk_ids <- fread(BULK_ID_FILE)
if (nrow(bulk_ids) != 25L || anyDuplicated(bulk_ids$gene) ||
    anyDuplicated(bulk_ids$bulk_gene_id) || !setequal(bulk_ids$gene, fixed$gene)) {
  stop("Panel 4C bulk gene-ID table must map the 25 fixed proteins to 25 unique Ensembl IDs")
}
bulk_ids[, gene_id_unversioned := sub("\\.[0-9]+$", "", bulk_gene_id)]
bulk_sel <- merge(
  bulk_ids[, .(gene, gene_id_unversioned)],
  bulk,
  by = "gene_id_unversioned",
  all.x = TRUE,
  sort = FALSE
)
bulk_hits <- bulk_sel[, .(n_rows = sum(!is.na(gene_id))), by = gene]
if (any(bulk_hits$n_rows != 1L)) {
  stop(
    "Each fixed gene must match exactly one bulk row by Ensembl ID; offending: ",
    paste0(bulk_hits[n_rows != 1L, paste0(gene, "=", n_rows)], collapse = ", ")
  )
}
if (any(is.na(bulk_sel$symbol) | bulk_sel$symbol != bulk_sel$gene)) {
  stop(
    "Bulk symbol disagrees with the fixed display gene for: ",
    paste(bulk_sel[is.na(symbol) | symbol != gene, paste0(gene, "/", gene_id, "/", symbol)], collapse = ", ")
  )
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
# bulk_sig mirrors is_canonical_deg() in scripts/figures/load_figure_data.R:
# conventional padj < 0.05 and |log2FC| > 0.5.
panel_pairs[, `:=`(
  bulk_sig = is.finite(bulk_padj) & bulk_padj < 0.05 & abs(bulk_logFC) > 0.5,
  protein_sig = is.finite(protein_padj) & protein_padj < 0.05,
  direction_concordant = sign(bulk_logFC) == sign(protein_logFC),
  selection_conditioned = TRUE,
  interpretation = "descriptive_same_cohort_reestimate"
)]
fwrite(panel_pairs, file.path(OUT, "panel4c_mrna_protein.tsv"), sep = "\t", quote = FALSE)

# Partial Spearman associations for all 25 proteins x five prespecified
# histology features. Each pair uses one complete-case participant set in which
# the protein, the histology feature, and every nuisance covariate are observed.
# Protein and histology are ranked within that set, and both ranks are
# residualized on the identical nuisance design for the same participants.
# age_z and bmi_z keep their all-58 scaling: least-squares residuals from a
# design with an intercept do not change under an affine rescaling of a
# covariate, so re-z-scoring within the set would give the same residuals.
# The t approximation uses df = n - rank(nuisance design) - 1, i.e. n - 2 minus
# the non-intercept nuisance columns (n - 6 for batch, age, BMI, and sex).
nuisance_vars <- c("acquisition_batch", "age_z", "bmi_z", "sex")
covariates_ok <- complete.cases(meta[, ..nuisance_vars])
features <- c("Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS")
partial_spearman <- function(g, feature) {
  protein <- gene_expr[g, samples]
  score <- meta[[feature]]
  protein_ok <- is.finite(protein)
  feature_ok <- is.finite(score)
  ok <- protein_ok & feature_ok & covariates_ok
  excluded <- c(
    protein_not_quantified = sum(!protein_ok),
    histology_missing = sum(!feature_ok),
    covariate_missing = sum(!covariates_ok)
  )
  missing_reason <- if (any(excluded > 0L)) {
    paste0(names(excluded)[excluded > 0L], ":", excluded[excluded > 0L], collapse = ";")
  } else {
    "none"
  }
  d <- droplevels(meta[ok])
  rhs <- c("age_z", "bmi_z")
  if (nlevels(d$acquisition_batch) > 1L) rhs <- c("acquisition_batch", rhs)
  if (nlevels(d$sex) > 1L) rhs <- c(rhs, "sex")
  x <- model.matrix(reformulate(rhs), data = d)
  design_rank <- qr(x)$rank
  out <- data.table(
    gene = g,
    feature,
    n = sum(ok),
    rho = NA_real_,
    pvalue = NA_real_,
    n_nuisance_terms = design_rank - 1L,
    df = sum(ok) - design_rank - 1L,
    missing_reason,
    test_status = "tested"
  )
  if (design_rank < ncol(x)) return(out[, test_status := "rank_deficient_nuisance_design"])
  if (out$n < 8L || out$df < 3L) return(out[, test_status := "too_few_participants"])
  protein_resid <- residuals(lm.fit(x, rank(protein[ok], ties.method = "average")))
  feature_resid <- residuals(lm.fit(x, rank(score[ok], ties.method = "average")))
  r <- suppressWarnings(cor(protein_resid, feature_resid))
  if (!is.finite(r)) return(out[, test_status := "no_residual_variance"])
  t_stat <- r * sqrt(out$df / (1 - r^2))
  out[, `:=`(rho = r, pvalue = 2 * pt(-abs(t_stat), out$df))]
}
histology <- rbindlist(lapply(fixed$gene, function(g) {
  rbindlist(lapply(features, function(feature) partial_spearman(g, feature)))
}))
histology[, padj := p.adjust(pvalue, method = "BH")]
histology[, correction_family := "25_fixed_proteins_x_5_histology_features"]
setcolorder(histology, c("gene", "feature", "n", "rho", "pvalue", "padj"))
fwrite(histology, file.path(OUT, "panel4c_histology_partial.tsv"), sep = "\t", quote = FALSE)

# Full-family enrichment on the adjusted protein t statistic. The six process
# displays are extracted only after BH correction across the combined Hallmark,
# Reactome, and curated-set family.
message("[protein] loading frozen MSigDB 2025.1 cache")
mdb <- as.data.table(readRDS(MSIGDB_RDS))
stopifnot(all(c("db_gene_symbol", "gs_name", "gs_collection", "gs_subcollection") %in% names(mdb)))
reactome <- mdb[
  gs_collection == "C2" & gs_subcollection == "CP:REACTOME",
  .(gs_name, gene_symbol = db_gene_symbol)
]
hallmark <- mdb[
  gs_collection == "H",
  .(gs_name, gene_symbol = db_gene_symbol)
]
msig <- rbindlist(list(
  reactome[, .(pathway = gs_name, gene = gene_symbol, family = "Reactome")],
  hallmark[, .(pathway = gs_name, gene = gene_symbol, family = "Hallmark")]
))
Rset <- function(name) unique(reactome[gs_name == name, gene_symbol])
lipid_droplet <- c(
  "PLIN1", "PLIN2", "PLIN3", "PLIN4", "PLIN5", "ABHD5", "CIDEA", "CIDEB",
  "CIDEC", "G0S2", "HILPDA", "ACLY", "FASN", "ACACA", "SCD", "DGAT1",
  "DGAT2", "LPIN1"
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
msig <- rbind(
  msig,
  rbindlist(lapply(names(curated), function(p) {
    data.table(pathway = p, gene = curated[[p]], family = "Curated")
  }))
)
pathways <- split(msig$gene, msig$pathway)
ranked <- de[is.finite(t) & !duplicated(gene)]
ranks <- sort(setNames(ranked$t, ranked$gene), decreasing = TRUE)
pathways <- lapply(pathways, intersect, y = names(ranks))
fg <- as.data.table(fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = 8,
  maxSize = 600,
  eps = 0,
  nproc = 1,
  BPPARAM = BiocParallel::SerialParam()
))
fg[, family := msig$family[match(pathway, msig$pathway)]]
fg[, abs_nes := abs(NES)]
setorder(fg, padj, -abs_nes)
fwrite(
  fg[, .(pathway, family, size, NES, pval, padj, leadingEdge)],
  file.path(OUT, "protein_pathway_enrichment.tsv"),
  sep = "\t",
  quote = FALSE
)

process_map <- data.table(
  program = c(
    "Lipid droplet", "ECM / BM", "Immune / inflammation", "Detox",
    "AA catabolism", "Mito / OXPHOS"
  ),
  pathway = c(
    "CURATED_LIPID_DROPLET",
    "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
    "CURATED_IMMUNE_INFLAMMATION",
    "CURATED_DETOX_PHASE_I_II",
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
    "CURATED_MITO_OXPHOS_TCA"
  ),
  representable = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  display_order = 1:6
)
process_enrichment <- merge(process_map, fg, by = "pathway", all.x = TRUE, sort = FALSE)
setorder(process_enrichment, display_order)
fwrite(
  process_enrichment,
  file.path(OUT, "panel4c_process_enrichment.tsv"),
  sep = "\t",
  quote = FALSE
)

# Frozen Hotspot-module projection. A protein is measured only if it passed the
# group-blind source-compatible missingness filter above.
z_all <- t(scale(t(gene_expr[, samples, drop = FALSE])))
z_all[!is.finite(z_all)] <- NA_real_

fit_score <- function(score, subset = rep(TRUE, length(score))) {
  d <- copy(meta[subset])
  d[, score := score[subset]]
  d <- d[is.finite(score)]
  d[, `:=`(
    group = droplevels(group),
    acquisition_batch = droplevels(acquisition_batch),
    sex = droplevels(sex)
  )]
  if (nrow(d) <= 8L || length(unique(d$group)) < 2L) {
    return(c(effect = NA_real_, se = NA_real_, pvalue = NA_real_))
  }
  rhs <- c("group", "age_z", "bmi_z")
  if (nlevels(d$acquisition_batch) > 1L) rhs <- c(rhs, "acquisition_batch")
  if (nlevels(d$sex) > 1L) rhs <- c(rhs, "sex")
  mod <- lm(reformulate(rhs, response = "score"), data = d)
  co <- summary(mod)$coefficients
  if (!"groupMASLD" %in% rownames(co)) {
    return(c(effect = NA_real_, se = NA_real_, pvalue = NA_real_))
  }
  c(effect = co["groupMASLD", 1], se = co["groupMASLD", 2], pvalue = co["groupMASLD", 4])
}

score_rows <- list()
score_results <- list()
for (i in seq_len(nrow(registry))) {
  pid <- registry$program_id[i]
  m <- membership[program_id == pid & mapped_symbol == TRUE & !is.na(gene_symbol)]
  m <- m[, .(original_l1_weight = sum(original_l1_weight)), by = gene_symbol]
  measured <- m[gene_symbol %in% rownames(z_all)]
  n_measured <- nrow(measured)
  retained_weight <- sum(measured$original_l1_weight)
  testable <- n_measured >= 8L && retained_weight >= 0.20
  base <- data.table(
    program_id = pid,
    n_measured = n_measured,
    retained_l1_weight = retained_weight,
    testable = testable
  )
  if (!testable) {
    score_results[[pid]] <- base[, `:=`(
      effect = NA_real_, se = NA_real_, pvalue = NA_real_,
      equal_effect = NA_real_, leave_top_effect = NA_real_,
      effect_2019 = NA_real_, effect_2020 = NA_real_, interaction_p = NA_real_
    )]
    next
  }
  measured[, weight := original_l1_weight / sum(original_l1_weight)]
  zz <- z_all[measured$gene_symbol, samples, drop = FALSE]
  observed_weight <- colSums(is.finite(zz) * measured$weight)
  weighted_score <- colSums(zz * measured$weight, na.rm = TRUE) / observed_weight
  weighted_score[!is.finite(weighted_score)] <- NA_real_
  equal_score <- colMeans(zz, na.rm = TRUE)
  equal_score[!is.finite(equal_score)] <- NA_real_
  top_gene <- measured$gene_symbol[which.max(measured$weight)]
  leave <- measured[gene_symbol != top_gene]
  leave[, weight := original_l1_weight / sum(original_l1_weight)]
  leave_zz <- z_all[leave$gene_symbol, samples, drop = FALSE]
  leave_observed_weight <- colSums(is.finite(leave_zz) * leave$weight)
  leave_score <- colSums(leave_zz * leave$weight, na.rm = TRUE) / leave_observed_weight
  leave_score[!is.finite(leave_score)] <- NA_real_

  primary <- fit_score(weighted_score)
  equal <- fit_score(equal_score)
  leave_top <- fit_score(leave_score)
  b19 <- fit_score(weighted_score, meta$acquisition_batch == "2019")
  b20 <- fit_score(weighted_score, meta$acquisition_batch == "2020")
  int_dat <- copy(meta)
  int_dat[, score := weighted_score]
  int_fit <- lm(score ~ group * acquisition_batch + age_z + bmi_z + sex, data = int_dat)
  int_name <- grep("groupMASLD:acquisition_batch", names(coef(int_fit)), value = TRUE)
  int_p <- if (length(int_name) == 1L) {
    summary(int_fit)$coefficients[int_name, 4]
  } else {
    NA_real_
  }
  score_results[[pid]] <- base[, `:=`(
    effect = unname(primary["effect"]),
    se = unname(primary["se"]),
    pvalue = unname(primary["pvalue"]),
    equal_effect = unname(equal["effect"]),
    leave_top_effect = unname(leave_top["effect"]),
    effect_2019 = unname(b19["effect"]),
    effect_2020 = unname(b20["effect"]),
    interaction_p = int_p
  )]
  score_rows[[pid]] <- data.table(
    program_id = pid,
    sample_id = samples,
    group = as.character(meta$group),
    acquisition_batch = as.character(meta$acquisition_batch),
    score = weighted_score,
    equal_score = equal_score,
    leave_top_score = leave_score
  )
}
module_results <- rbindlist(score_results, fill = TRUE)
module_results[testable == TRUE, padj := p.adjust(pvalue, method = "BH")]
module_results[, sensitivity_sign_agree := testable & is.finite(effect) &
  sign(effect) == sign(equal_effect) & sign(effect) == sign(leave_top_effect)]
module_results[, robust := testable & is.finite(padj) & padj < 0.05 & sensitivity_sign_agree]
module_results <- merge(
  registry[, .(program_id, display_order, cell_type, module, program_name)],
  module_results,
  by = "program_id",
  all.x = TRUE,
  sort = FALSE
)
setorder(module_results, display_order)
fwrite(module_results, file.path(OUT, "module_protein_results.tsv"), sep = "\t", quote = FALSE)
if (length(score_rows)) {
  module_scores <- rbindlist(score_rows)
  fwrite(
    module_scores,
    file.path(OUT, "module_protein_scores.tsv"),
    sep = "\t",
    quote = FALSE
  )

  # Severity-resolved protein realization of the independently frozen Figure 3
  # modules. The primary cohort is MASLD-only so correlations cannot be driven
  # by control-versus-case separation. Both module score and histology are
  # residualized on the same nuisance design before Spearman correlation. BH
  # spans every testable program x five prespecified histology features.
  severity_features <- c("Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS")
  severity_meta <- meta[, c(
    "sample_id", "sex", "age_z", "bmi_z", severity_features
  ), with = FALSE]
  severity_data <- merge(
    module_scores,
    severity_meta,
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )[group == "MASLD"]

  rank_residualize_severity <- function(y, d, include_batch = TRUE) {
    rhs <- c("age_z", "bmi_z")
    if (include_batch && uniqueN(d$acquisition_batch[!is.na(d$acquisition_batch)]) > 1L) {
      rhs <- c("acquisition_batch", rhs)
    }
    if (uniqueN(d$sex[!is.na(d$sex)]) > 1L) rhs <- c(rhs, "sex")
    nuisance <- model.matrix(reformulate(rhs), data = d)
    ok <- is.finite(y) & apply(nuisance, 1, function(z) all(is.finite(z)))
    out <- rep(NA_real_, length(y))
    if (sum(ok) > ncol(nuisance) + 2L) {
      out[ok] <- residuals(lm.fit(
        nuisance[ok, , drop = FALSE],
        rank(y[ok], ties.method = "average")
      ))
    }
    out
  }

  severity_cor <- function(d, score_col, feature, include_batch = TRUE) {
    score_resid <- rank_residualize_severity(d[[score_col]], d, include_batch)
    hist_resid <- rank_residualize_severity(d[[feature]], d, include_batch)
    ok <- is.finite(score_resid) & is.finite(hist_resid)
    if (sum(ok) < 8L) {
      return(c(n = sum(ok), rho = NA_real_, pvalue = NA_real_))
    }
    tst <- suppressWarnings(cor.test(
      score_resid[ok], hist_resid[ok], method = "pearson"
    ))
    c(n = sum(ok), rho = unname(tst$estimate), pvalue = tst$p.value)
  }

  severity_rows <- rbindlist(lapply(unique(module_scores$program_id), function(pid) {
    d <- severity_data[program_id == pid]
    rbindlist(lapply(severity_features, function(feature) {
      primary <- severity_cor(d, "score", feature)
      equal <- severity_cor(d, "equal_score", feature)
      leave_top <- severity_cor(d, "leave_top_score", feature)
      batch_rho <- vapply(c("2019", "2020"), function(batch_id) {
        unname(severity_cor(
          d[acquisition_batch == batch_id], "score", feature, include_batch = FALSE
        )["rho"])
      }, numeric(1))
      data.table(
        program_id = pid,
        feature = feature,
        n = as.integer(primary["n"]),
        rho = unname(primary["rho"]),
        pvalue = unname(primary["pvalue"]),
        equal_rho = unname(equal["rho"]),
        leave_top_rho = unname(leave_top["rho"]),
        rho_2019 = batch_rho[[1]],
        rho_2020 = batch_rho[[2]]
      )
    }))
  }))

  severity_template <- merge(
    CJ(
      program_id = registry$program_id,
      feature = severity_features,
      unique = TRUE,
      sorted = FALSE
    ),
    module_results[, .(
      program_id, display_order, cell_type, module, program_name,
      n_measured, retained_l1_weight, testable
    )],
    by = "program_id",
    all.x = TRUE,
    sort = FALSE
  )
  module_histology <- merge(
    severity_template,
    severity_rows,
    by = c("program_id", "feature"),
    all.x = TRUE,
    sort = FALSE
  )
  module_histology[, feature_order := match(feature, severity_features)]
  module_histology[testable == TRUE & is.finite(pvalue),
                   padj := p.adjust(pvalue, method = "BH")]
  module_histology[, sensitivity_sign_agree := testable & is.finite(rho) &
    is.finite(equal_rho) & is.finite(leave_top_rho) &
    sign(rho) == sign(equal_rho) & sign(rho) == sign(leave_top_rho)]
  module_histology[, batch_sign_agree := testable & is.finite(rho) &
    is.finite(rho_2019) & is.finite(rho_2020) &
    sign(rho) == sign(rho_2019) & sign(rho) == sign(rho_2020)]
  module_histology[, robust := testable & is.finite(padj) & padj < 0.05 &
    sensitivity_sign_agree & batch_sign_agree]
  module_histology[, `:=`(
    cohort = "MASLD_only",
    correction_scope = "17_testable_programs_x_5_histology_features",
    interpretation = "independent_frozen_program_severity_association"
  )]
  setorder(module_histology, display_order, feature_order)
  fwrite(
    module_histology,
    file.path(OUT, "module_protein_histology_masld.tsv"),
    sep = "\t",
    quote = FALSE
  )
}

cat("[protein] batch counts:\n")
print(meta[, .N, by = .(acquisition_batch, group)])
cat("[protein] fixed Panel 4C direction retained:",
    sum(panel_pairs$direction_concordant, na.rm = TRUE), "/25\n")
cat("[protein] fixed Panel 4C adjusted q<0.05:",
    sum(panel_pairs$protein_sig, na.rm = TRUE), "/25\n")
cat("[protein] Panel 4C partial Spearman tested:", sum(histology$test_status == "tested"),
    "/125; BH q<0.05:", sum(histology$padj < 0.05, na.rm = TRUE), "\n")
cat("[protein] module testability:", sum(module_results$testable), "/22; robust:",
    sum(module_results$robust), "\n")
if (exists("module_histology")) {
  cat("[protein] MASLD-only module-histology robust:",
      sum(module_histology$robust), "/",
      sum(module_histology$testable), "testable program-feature pairs\n")
}
