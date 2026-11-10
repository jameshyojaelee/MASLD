# CONF1 C3 (PRESPEC_C3_C4.md): composition versus within-lineage partition of replicated stage genes.
suppressPackageStartupMessages({library(data.table); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
if (is.na(OUT) || dir.exists(OUT)) stop("output directory missing or already exists")
dir.create(OUT)
NCORE <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
HYP <- file.path(ROOT, "Analysis/Hypotheses")
SS <- file.path(HYP, "SS-20261025T135237Z")
V5 <- file.path(HYP, "V5-20261021T122723Z/outputs")
V3 <- file.path(HYP, "V3-20261008T015239Z/pseudobulk")
NF <- file.path(HYP, "NUCLEUS-FIX-20261107/unique_counts")
ATLAS_LIN <- c(Hepatocytes = "Hepatocytes", Cholangiocytes = "Cholangiocytes", Fibroblasts = "Fibroblasts",
               Macrophages = "Macrophages", T_cells = "T_cells", Endothelial = "Endothelial_cells")
V5_LIN <- c(Hepatocytes = "hep", Cholangiocytes = "cholangiocytes", Fibroblasts = "fibroblasts",
            Macrophages = "macrophages", T_cells = "tcells")
V3_LIN <- c(Hepatocytes = "hepatocyte", Cholangiocytes = "cholangiocyte", Fibroblasts = "stellate_fibroblast")
IN <- c(
  sidecar = file.path(SS, "out_s3b-20261028T035515Z/steatosite_catalog_sidecar.tsv.gz"),
  bulk_unadjusted = file.path(SS, "out_s3a-20261026T234211Z/gene_bin.tsv.gz"),
  bulk_fits = file.path(SS, "out_s3a-20261026T234211Z/observed_fits.rds"),
  v5_donors = file.path(V5, "donor_qc.tsv"), v5_cells = file.path(V5, "job3/cells_per_donor_set.tsv"),
  v3_donors = file.path(V3, "donor_metadata.tsv"), v3_cells = file.path(V3, "nuclei_per_donor_lineage.tsv"),
  atlas_donors = file.path(NF, "metadata/donor_metadata.tsv"), atlas_libraries = file.path(NF, "metadata/library_metadata.tsv"),
  atlas_sex = file.path(HYP, "RV1-20261031T155620Z/out-20261102T130905Z/donors.tsv"),
  annotation = file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
  setNames(file.path(V5, sprintf("job3/pseudobulk_%s.tsv.gz", V5_LIN)), paste0("v5_", names(V5_LIN))),
  setNames(file.path(V3, sprintf("pseudobulk_%s.tsv.gz", V3_LIN)), paste0("v3_", names(V3_LIN))),
  setNames(file.path(NF, sprintf("pseudobulk/%s_pseudobulk.csv", ATLAS_LIN)), paste0("atlas_", names(ATLAS_LIN))))
stopifnot(all(file.exists(IN)))
SEED <- 20261109; N_PERM <- 20L; MIN_CELLS <- 30; MIN_GROUP <- 3
wr <- function(x, f) fwrite(as.data.table(x), file.path(OUT, f), sep = "\t", na = "NA")
wr(data.table(input = names(IN), path = IN, sha256 = vapply(IN, function(f)
  sub(" .*", "", system2("sha256sum", shQuote(f), stdout = TRUE)), "")), "inputs_sha256.tsv")

# ---- bulk side ---------------------------------------------------------------------------
sc <- fread(IN["sidecar"])
sc <- sc[!steatosite_stage_state %in% c("not a discovery F3F4 hit", "not in release gene universe")]
stopifnot(nrow(sc) == 8542)
un <- fread(IN["bulk_unadjusted"])[coef == "binF3F4", .(gene, bulk_log2fc = logFC, bulk_se = SE, bulk_p = P)]
ad <- as.data.table(readRDS(IN["bulk_fits"])$comp$tab)
stopifnot(all(c("gene", "coef", "logFC", "SE", "P") %in% names(ad)))
ad <- ad[coef == "binF3F4", .(gene, adj_log2fc = logFC, adj_se = SE, adj_p = P)]
G <- merge(merge(sc[, .(ensembl_id, symbol, gene = gene_id_versioned, primary_evidence_class, steatosite_stage_state,
                        D0_F3F4_log2FC)], un, by = "gene", all.x = TRUE), ad, by = "gene", all.x = TRUE)
G[, replicated := steatosite_stage_state == "replicated"]
stopifnot(sum(G$replicated) == 7214, !anyNA(G[replicated == TRUE, bulk_log2fc]))
G[, retained := adj_log2fc / bulk_log2fc]
G[replicated == TRUE, adj_q := p.adjust(adj_p, "BH")]
dup_sym <- G$symbol[duplicated(G$symbol) | is.na(G$symbol) | G$symbol == ""]
TARGET <- setdiff(G[replicated == TRUE, symbol], dup_sym)

# ---- lineage pseudobulks -------------------------------------------------------------------
annot <- fread(IN["annotation"], select = c("ensembl_base", "gene_name"))
by_symbol <- function(counts, symbol) {
  ok <- !is.na(symbol) & symbol != "" & !(symbol %in% symbol[duplicated(symbol)])
  counts <- counts[ok, , drop = FALSE]; rownames(counts) <- symbol[ok]; counts
}
disease_of <- function(stage_ord) fifelse(stage_ord == 0, 0, fifelse(stage_ord >= 2, 1, NA_real_))
load_v5 <- function() {
  dq <- fread(IN["v5_donors"])[included == TRUE]; cells <- fread(IN["v5_cells"])
  meta <- data.table(donor = dq$gsm, female = as.integer(dq$sex_f), disease = as.numeric(dq$disease))
  lapply(V5_LIN, function(s) { pb <- fread(file.path(V5, sprintf("job3/pseudobulk_%s.tsv.gz", s)))
    list(meta = copy(meta)[, n_cells := cells[set == s][match(meta$donor, gsm), n_cells]],
         counts = by_symbol(as.matrix(pb[, meta$donor, with = FALSE]),
                            annot$gene_name[match(sub("\\..*$", "", pb$gene_id), annot$ensembl_base)])) })
}
load_v3 <- function() {
  m <- fread(IN["v3_donors"]); n3 <- fread(IN["v3_cells"])
  meta <- data.table(donor = m$donor, female = as.integer(m$sex == "female"), disease = disease_of(m$stage_ordinal))
  lapply(V3_LIN, function(s) { pb <- fread(file.path(V3, sprintf("pseudobulk_%s.tsv.gz", s)))
    list(meta = copy(meta)[, n_cells := n3[[s]][match(meta$donor, n3$donor)]],
         counts = by_symbol(as.matrix(pb[, meta$donor, with = FALSE]), pb$symbol)) })
}
load_atlas <- function(ds) {
  dm <- fread(IN["atlas_donors"])[dataset == ds & toupper(as.character(exclude)) == "FALSE"]
  lib <- fread(IN["atlas_libraries"])[donor %in% dm$donor]
  sex <- fread(IN["atlas_sex"], select = c("donor", "sex"))
  stopifnot(all(dm$donor %in% sex$donor))
  meta <- data.table(donor = dm$donor, female = as.integer(sex$sex[match(dm$donor, sex$donor)] == "Female"),
                     disease = disease_of(dm$stage_ordinal))
  lapply(ATLAS_LIN, function(s) { pb <- fread(file.path(NF, sprintf("pseudobulk/%s_pseudobulk.csv", s)))
    ncol_name <- paste0("n_", s); stopifnot(ncol_name %in% names(lib))
    libs <- lapply(setNames(dm$donor, dm$donor), function(d) intersect(lib[donor == d, sample], names(pb)))
    has <- lengths(libs) > 0
    counts <- vapply(libs[has], function(l) rowSums(as.matrix(pb[, l, with = FALSE])), numeric(nrow(pb)))
    mm <- copy(meta)[has][, n_cells := vapply(donor, function(d) sum(lib[donor == d][[ncol_name]]), numeric(1))]
    list(meta = mm, counts = by_symbol(counts, pb$gene)) })
}
CO <- list(GSE298719 = load_v5(), GSE256398 = load_v3(), GSE202379 = load_atlas("GSE202379"), GSE244832 = load_atlas("GSE244832"))
UNITS <- rbindlist(lapply(names(CO), function(k) rbindlist(lapply(names(CO[[k]]), function(l) {
  m <- CO[[k]][[l]]$meta[!is.na(disease) & !is.na(n_cells) & n_cells >= MIN_CELLS]
  data.table(cohort = k, lineage = l, n0 = sum(m$disease == 0), n1 = sum(m$disease == 1))
}))))
UNITS[, usable := n0 >= MIN_GROUP & n1 >= MIN_GROUP]
UNITS[, weight := sqrt(n0 * n1 / (n0 + n1))]
wr(UNITS, "lineage_cohorts.tsv")
TASKS <- UNITS[usable == TRUE]

# Donor-level disease labels per cohort: observed (column 1) and permuted.
set.seed(SEED)
LAB <- lapply(CO, function(x) { m <- x$Hepatocytes$meta[!is.na(disease), .(donor, disease)]
  P <- cbind(m$disease, vapply(seq_len(N_PERM), function(i) sample(m$disease), numeric(nrow(m))))
  rownames(P) <- m$donor; P })

fit_one <- function(k, l, r) {
  x <- CO[[k]][[l]]
  m <- x$meta[!is.na(disease) & !is.na(n_cells) & n_cells >= MIN_CELLS & donor %in% rownames(LAB[[k]])]
  dis <- LAB[[k]][m$donor, r]
  if (sum(dis == 0) < 2 || sum(dis == 1) < 2) return(NULL)
  X <- cbind(1, disease = dis, female = m$female)
  if (qr(X)$rank < ncol(X)) X <- X[, 1:2]
  y <- DGEList(x$counts[, m$donor, drop = FALSE])
  y <- calcNormFactors(y[filterByExpr(y, design = X), , keep.lib.sizes = FALSE])
  fit <- glmQLFit(estimateDisp(y, X), X, robust = TRUE)
  tt <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
  keep <- rownames(tt) %in% TARGET
  data.table(cohort = k, lineage = l, symbol = rownames(tt)[keep], logFC = tt$logFC[keep], P = tt$PValue[keep])
}
stouffer <- function(de) {
  de <- merge(de, TASKS[, .(cohort, lineage, weight)], by = c("cohort", "lineage"))
  de[, z := sign(logFC) * qnorm(pmax(P, 1e-300) / 2, lower.tail = FALSE)]
  s <- de[, .(n_cohorts = .N, Z = sum(weight * z) / sqrt(sum(weight^2)), n_pos = sum(logFC > 0), n_neg = sum(logFC < 0),
              mean_logFC = sum(weight^2 * logFC) / sum(weight^2)), by = .(lineage, symbol)][n_cohorts >= 2]
  s[, p := 2 * pnorm(-abs(Z))]; s[, q := p.adjust(p, "BH")]
  s[, bulk_sign := sign(G$bulk_log2fc[match(symbol, G$symbol)])]
  s[, call := q < 0.05 & sign(Z) == bulk_sign & fifelse(Z > 0, n_pos, n_neg) >= 2]
  s
}
run <- function(r) stouffer(rbindlist(lapply(seq_len(nrow(TASKS)), function(i) fit_one(TASKS$cohort[i], TASKS$lineage[i], r))))
call_counts <- function(s) data.table(hepatocyte_calls = s[lineage == "Hepatocytes", sum(call)],
  genes_any_lineage = uniqueN(s[call == TRUE, symbol]), pairs_tested = nrow(s), pairs_q05_any_sign = sum(s$q < 0.05))

de_obs <- rbindlist(lapply(seq_len(nrow(TASKS)), function(i) fit_one(TASKS$cohort[i], TASKS$lineage[i], 1L)))
fwrite(de_obs, file.path(OUT, "lineage_cohort_de.tsv.gz"), sep = "\t")
S <- stouffer(de_obs)
fwrite(S, file.path(OUT, "lineage_stouffer.tsv.gz"), sep = "\t")
perm <- rbindlist(mclapply(seq_len(N_PERM), function(i) cbind(run = i, call_counts(run(i + 1L))), mc.cores = NCORE, mc.preschedule = FALSE))
stopifnot(nrow(perm) == N_PERM)
obs <- call_counts(S)
wr(rbind(cbind(run = 0L, obs), perm)[, `:=`(label = fifelse(run == 0L, "observed", "permuted"))], "permutation_calls.tsv")
nc_pass <- obs$hepatocyte_calls > max(perm$hepatocyte_calls) && obs$genes_any_lineage > max(perm$genes_any_lineage)

# ---- classes ---------------------------------------------------------------------------------
hep <- S[lineage == "Hepatocytes", .(symbol, hep_Z = Z, hep_q = q, hep_logFC = mean_logFC, hep_n_cohorts = n_cohorts, hep_call = call)]
oth <- S[lineage != "Hepatocytes" & call == TRUE, .(other_lineages = paste(sort(lineage), collapse = ";")), by = symbol]
G <- merge(merge(G, hep, by = "symbol", all.x = TRUE), oth, by = "symbol", all.x = TRUE)
G[symbol %in% dup_sym, c("hep_Z", "hep_q", "hep_logFC", "hep_n_cohorts", "hep_call", "other_lineages") := NA]
G[, same_sign_adj := sign(adj_log2fc) == sign(bulk_log2fc)]
G[replicated == TRUE, class := fifelse(!is.na(hep_call) & hep_call, "within_hepatocyte",
  fifelse(!is.na(other_lineages), "within_other_lineage",
  fifelse(!is.na(adj_q) & adj_q >= 0.05 & (retained < 0.5 | !same_sign_adj), "composition_associated",
  fifelse(!is.na(adj_q) & adj_q < 0.05 & same_sign_adj & retained >= 0.5, "bulk_composition_independent", "unresolved"))))]
G[replicated == TRUE, bulk_only_class := fifelse(!is.na(adj_q) & adj_q >= 0.05 & (retained < 0.5 | !same_sign_adj), "composition_associated",
  fifelse(!is.na(adj_q) & adj_q < 0.05 & same_sign_adj & retained >= 0.5, "bulk_composition_independent", "unresolved"))]
G[, negative_control_passed := nc_pass]
fwrite(G, file.path(OUT, "gene_partition.tsv.gz"), sep = "\t", na = "NA")
wr(rbindlist(list(
  G[replicated == TRUE, .(scheme = "with_lineage", direction = fifelse(bulk_log2fc > 0, "up", "down"), class)][, .N, by = .(scheme, direction, class)],
  G[replicated == TRUE, .(scheme = "bulk_only", direction = fifelse(bulk_log2fc > 0, "up", "down"), class = bulk_only_class)][, .N, by = .(scheme, direction, class)],
  G[replicated == TRUE, .(scheme = "with_lineage_by_catalog_class", direction = primary_evidence_class, class)][, .N, by = .(scheme, direction, class)],
  data.table(scheme = "summary", direction = c("median_retained", "frac_retained_lt_0.5", "frac_adj_q_lt_0.05", "negative_control_passed"),
             class = NA_character_, N = c(G[replicated == TRUE, median(retained, na.rm = TRUE)], G[replicated == TRUE, mean(retained < 0.5, na.rm = TRUE)],
                                          G[replicated == TRUE, mean(adj_q < 0.05, na.rm = TRUE)], as.numeric(nc_pass)))), use.names = TRUE), "class_counts.tsv")
NAMED <- c("COL1A1", "COL1A2", "COL3A1", "LUM", "KRT7", "KRT19", "EPCAM", "SOX9", "SPP1", "GNMT", "MAT1A", "CYP2C19", "CCL20",
           "FCAMR", "TRIM31", "CXCL10", "IL32", "UBD", "AKR1B10", "HKDC1", "FADS2", "GAS6", "KRT8", "PKM", "TNFSF10", "IGFBP7", "THY1", "PDGFRB")
wr(G[symbol %in% NAMED, .(symbol, primary_evidence_class, steatosite_stage_state, bulk_log2fc, adj_log2fc, retained, adj_q, hep_logFC, hep_q,
                           hep_n_cohorts, other_lineages, class)][order(class, symbol)], "named_genes.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
