# DIG1 job 23 (exploratory; follows jobs 20-22). Which plasma proteins track the three liver protein
# axes in the same people? PXD051911, participants with liver protein and an initial plasma sample.
# Liver axes (fixed in jobs 21-22): CI_minus_hep_mito (Complex I relative to hepatocyte matrix
# enzymes), LD_coat (PLIN2, PLIN4, ABHD5), collagen (COL*A* chains). Plasma: deposited batch-corrected
# matrix, accessions mapped to genes as in job 20. Partial Spearman, ranks residualised on age, sex,
# BMI and plasma batch; BH over every plasma protein x 3 axes. Secondary: the same with the matching
# histology score added (ballooning for CI, steatosis for LD_coat, fibrosis for collagen), to ask
# whether the plasma protein reads the liver protein beyond the biopsy grade; and plasma vs histology
# in all initial participants for the proteins that pass. Calibration: same-gene liver-vs-plasma rho
# for every protein measured in both (liver tissue holds blood, so this is an upper reference).
# Usage: Rscript 23_plasma_liver_axes.R <new_out_dir>
suppressPackageStartupMessages({library(data.table)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
PX <- file.path(ROOT, "data/PXD051911"); RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
lv <- fread(file.path(PX, "liver_protein_quant.txt"))
amap <- unique(lv[!is.na(Genes) & Genes != "", .(acc = unlist(strsplit(ProteinAccessions, ";"))), by = Genes], by = "acc")
amap <- rbind(amap, data.table(Genes = c("IGFBP2", "IGF1"), acc = c("P18065", "P05019")))
pl <- fread(file.path(PX, "plasma_protein_quant_filtered_and_batch_corrected.txt"))
pl[, gene := vapply(strsplit(ProteinAccessions, ";"), function(a) { g <- amap$Genes[match(a, amap$acc)]; g <- g[!is.na(g)]; if (length(g)) g[1] else NA_character_ }, "")]
pl <- pl[!is.na(gene) & !grepl(";", gene, fixed = TRUE)]
pl <- pl[, lapply(.SD, median, na.rm = TRUE), by = gene, .SDcols = setdiff(names(pl), c("ProteinAccessions", "gene"))]
X <- as.matrix(pl[, -1]); rownames(X) <- pl$gene
md <- fread(file.path(PX, "meta_data.txt"))
pm <- md[unique_identifier %in% colnames(X) & sample_group == "initial_sample"]
pm[, `:=`(age = as.numeric(alder), bmi = as.numeric(bmi), sex = factor(gender), batch = factor(plasma_batch_effects),
          masld = as.integer(saf_diagnosis != "No_MASLD"), Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)),
          Steatosis = as.numeric(steatosis_score), Ballooning = as.numeric(hepatocellular_ballooning_score))]
pm <- unique(pm[complete.cases(pm[, .(age, bmi, sex, batch)])], by = "patient_name")
# liver axes, keyed by liver file -> patient
s22 <- fread(file.path(RUN, "out_mito_pattern/participant_scores.tsv"))
s21 <- fread(file.path(RUN, "out_protein_histology/participant_scores.tsv"))[, .(sample_id, collagen)]
la <- merge(s22[, .(sample_id, CI_minus_hep_mito, LD_coat)], s21, by = "sample_id")
lk <- unique(md[liver_proteomics_filename %in% la$sample_id, .(sample_id = liver_proteomics_filename, patient_name)], by = "sample_id")
la <- merge(la, lk, by = "sample_id")
D <- merge(la, pm[, .(patient_name, unique_identifier, age, bmi, sex, batch, masld, Fibrosis, Steatosis, Ballooning)], by = "patient_name")
stopifnot(!anyDuplicated(D$patient_name))
P <- X[, D$unique_identifier, drop = FALSE]
keep <- rowSums(!is.na(P)) >= 30; P <- P[keep, , drop = FALSE]
fwrite(data.table(n_participants = nrow(D), n_masld = sum(D$masld), n_plasma_proteins = nrow(P)), file.path(out, "design.tsv"), sep = "\t")
pc <- function(a, b, Z) {
  ok <- complete.cases(a, b, Z); if (sum(ok) < 20) return(c(NA, NA, sum(ok)))
  xx <- model.matrix(~ ., data = Z[ok, , drop = FALSE])
  rr <- qr.resid(qr(xx), cbind(rank(a[ok]), rank(b[ok]))); r <- cor(rr[, 1], rr[, 2]); df <- sum(ok) - ncol(xx) - 1
  c(r, 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df), sum(ok)) }
Zb <- as.data.frame(D[, .(age, sex, bmi, batch)])
AX <- c(CI_minus_hep_mito = "Ballooning", LD_coat = "Steatosis", collagen = "Fibrosis")
res <- rbindlist(lapply(names(AX), function(ax) rbindlist(lapply(rownames(P), function(p) {
  r0 <- pc(D[[ax]], P[p, ], Zb); r1 <- pc(D[[ax]], P[p, ], cbind(Zb, h = D[[AX[[ax]]]]))
  data.table(axis = ax, plasma = p, n = r0[3], rho = r0[1], p = r0[2], rho_plus_histology = r1[1], p_plus_histology = r1[2]) }))))
res[, bh := p.adjust(p, "BH")]
fwrite(res[order(p)], file.path(out, "plasma_vs_liver_axes.tsv"), sep = "\t")
# plasma vs histology in all initial participants, for proteins with bh < 0.10
hits <- unique(res[bh < 0.10, plasma])
Pa <- X[, pm$unique_identifier, drop = FALSE]; Za <- as.data.frame(pm[, .(age, sex, bmi, batch)])
hh <- rbindlist(lapply(hits, function(p) rbindlist(lapply(c("Ballooning", "Steatosis", "Fibrosis"), function(h) {
  r <- pc(pm[[h]], Pa[p, ], Za); data.table(plasma = p, histology = h, n = r[3], rho = r[1], p = r[2]) }))))
fwrite(hh, file.path(out, "hits_vs_histology_all_initial.tsv"), sep = "\t")
# calibration: same gene in liver tissue and plasma
raw <- lv[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
L <- as.matrix(raw[, D$sample_id, with = FALSE]); storage.mode(L) <- "double"; L[!is.finite(L) | L <= 0] <- NA; L <- log2(L)
rownames(L) <- raw$Genes; L <- L[!duplicated(rownames(L)), , drop = FALSE]
both <- intersect(rownames(L), rownames(P))
cal <- rbindlist(lapply(both, function(gg) { r <- pc(L[gg, ], P[gg, ], Zb); data.table(gene = gg, n = r[3], rho = r[1], p = r[2]) }))
fwrite(cal[order(-rho)], file.path(out, "same_gene_liver_vs_plasma.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(fread(file.path(out, "design.tsv"))); print(res[order(p)][1:30], digits = 2); print(summary(cal$rho))
