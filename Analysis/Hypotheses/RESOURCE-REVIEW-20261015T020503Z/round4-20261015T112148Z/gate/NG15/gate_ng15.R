#!/usr/bin/env Rscript
# NG15 (L06-04) feasibility counts. Does NOT estimate the sex effect on the ECM or ductular score.
# Outputs: residual SDs within stage, MDE, GTEx-predicted sex shift of program scores, sex-score positive control.
suppressPackageStartupMessages({ library(data.table); library(edgeR) })
set.seed(20260923)
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
gate <- file.path(repo, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG15")
out <- file.path(gate, "out", format(Sys.time(), "%Y%m%dT%H%M%S")); dir.create(out, recursive = TRUE)
sdir <- file.path(repo, "Analysis/Hypotheses/H3-20261009T180331Z/results/scores-20261009T181306Z")
ph <- fread(file.path(sdir, "person_phenotypes.tsv"))
pc <- readRDS(file.path(sdir, "person_counts.rds"))
memb <- fread(file.path(repo, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"))
ann <- fread(file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"))
ann[, symbol := toupper(gene_name)]
id2sym <- ann[symbol != "", .(symbol = symbol[1], chrom = chromosome[1]), by = ensembl_base]
PROG <- c(ecm = "hotspot_hepatocytes_f05c535ae5bbc0b9", duct = "hotspot_hepatocytes_48f39dd4d817a10e", upr = "hotspot_hepatocytes_b75cab32b41701a0")
mw <- memb[program_uid %in% PROG & !is.na(mapped_symbol) & mapped_symbol != "",
           .(weight = sum(as.numeric(original_l1_weight))), by = .(program_uid, symbol = toupper(mapped_symbol))]
es <- read.delim(file.path(repo, "data/gtex_sb_eqtl/GTEx_Analysis_v8_sbgenes/effect_size.tsv"), check.names = FALSE)
gtex_liver <- data.table(ensembl_base = sub("\\..*$", "", rownames(es)), beta = es[["Liver"]])
gtex_liver <- gtex_liver[!is.na(beta)][, .(beta = beta[1]), by = ensembl_base]
sig <- fread(file.path(repo, "data/gtex_sb_eqtl/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt"))[tissue == "Liver"]
sig[, ensembl_base := sub("\\..*$", "", gene)]
sig <- merge(sig, id2sym, by = "ensembl_base", all.x = TRUE)
cat("GTEx liver signif:", nrow(sig), " autosomal:", sig[!chrom %in% c("chrX", "chrY") & !is.na(chrom), .N], "\n")
sig_auto <- sig[!chrom %in% c("chrX", "chrY") & !is.na(chrom)]

cohorts <- unique(ph$cohort)
gt <- list(); pcrows <- list(); rs <- list(); recheck <- list()
for (co in cohorts) {
  idx <- which(pc$cohort == co)
  m <- pc$counts[, idx, drop = FALSE]
  lc <- cpm(DGEList(m), log = TRUE, prior.count = 1)
  lc <- lc[rowSums(m) > ncol(m), , drop = FALSE]
  sym <- id2sym$symbol[match(rownames(lc), id2sym$ensembl_base)]
  ok <- !is.na(sym)
  s <- rowsum(lc[ok, , drop = FALSE], sym[ok], reorder = FALSE)
  s <- s / as.numeric(table(sym[ok])[rownames(s)])
  mu <- rowMeans(s); sdv <- apply(s, 1, sd); keep <- is.finite(sdv) & sdv > 0
  z <- (s[keep, ] - mu[keep]) / sdv[keep]
  pco <- ph[match(colnames(m), person_id)]
  for (nm in names(PROG)) {
    mm <- mw[program_uid == PROG[[nm]] & symbol %in% rownames(z)]
    w <- mm$weight / sum(mm$weight)
    v <- as.numeric(crossprod(w, z[mm$symbol, , drop = FALSE]))
    sdv_score <- sd(v)
    sc <- (v - mean(v)) / sdv_score
    pco[[paste0("re_", nm)]] <- sc
    # GTEx-predicted female-minus-male shift in score SD units (healthy-liver sex effect, GTEx v8)
    ens <- id2sym$ensembl_base[match(mm$symbol, id2sym$symbol)]
    b <- gtex_liver$beta[match(ens, gtex_liver$ensembl_base)]
    okb <- !is.na(b)
    pred <- sum(w[okb] * b[okb] / sdv[mm$symbol[okb]]) / sdv_score
    gt[[length(gt) + 1]] <- data.table(cohort = co, program = nm, n_genes = nrow(mm), n_gtex = sum(okb),
                                       gtex_weight_frac = sum(w[okb]), pred_F_minus_M_sd = pred)
  }
  # positive control 2: GTEx autosomal female-biased minus male-biased score
  fb <- intersect(sig_auto[effsize > 0, symbol], rownames(z)); mb <- intersect(sig_auto[effsize < 0, symbol], rownames(z))
  pco[, sexscore := colMeans(z[fb, , drop = FALSE]) - colMeans(z[mb, , drop = FALSE])]
  recheck[[co]] <- data.table(cohort = co, cor_ecm = cor(pco$re_ecm, pco$score_ecm, use = "complete"),
                              cor_upr = cor(pco$re_upr, pco$score_upr, use = "complete"), n_fb = length(fb), n_mb = length(mb))
  d <- pco[group_binary == "Disease" & !is.na(stage) & sex %in% c("F", "M")]
  if (nrow(d) == 0) next
  pcrows[[co]] <- data.table(cohort = co, nF = d[sex == "F", .N], nM = d[sex == "M", .N],
                             d_sexscore = (mean(d[sex == "F", sexscore]) - mean(d[sex == "M", sexscore])) / sd(d$sexscore),
                             auc_sexscore = { r <- rank(d$sexscore); nf <- sum(d$sex == "F"); nm <- sum(d$sex == "M")
                               (sum(r[d$sex == "F"]) - nf * (nf + 1) / 2) / (nf * nm) })
  # residual SD within stage (sex NOT in the model; no sex effect estimated)
  d[, st := factor(stage)]
  tab <- table(d$st); if (any(tab < 3)) { lv <- as.numeric(names(tab)); for (l in names(tab)[tab < 3]) { tgt <- names(tab)[tab >= 3][which.min(abs(as.numeric(names(tab)[tab >= 3]) - as.numeric(l)))]; d[st == l, st := tgt] }; d[, st := droplevels(st)] }
  hasnas <- all(!is.na(d$nas_score))
  for (nm in names(PROG)) {
    y <- d[[paste0("re_", nm)]]
    f <- if (hasnas) lm(y ~ st + nas_score, data = d) else lm(y ~ st, data = d)
    rs[[length(rs) + 1]] <- data.table(cohort = co, program = nm, n = nrow(d), nF = sum(d$sex == "F"), nM = sum(d$sex == "M"),
                                       nas_in_model = hasnas, resid_sd = sigma(f), total_sd = sd(y))
  }
}
gt <- rbindlist(gt); pcr <- rbindlist(pcrows); rs <- rbindlist(rs); recheck <- rbindlist(recheck)
rs[, se_k := resid_sd * sqrt(1 / nF + 1 / nM)]
prim <- rs[cohort != "GSE193066"]
mde <- prim[, .(k = .N, n = sum(n), nF = sum(nF), nM = sum(nM), se_fe = sqrt(1 / sum(1 / se_k^2))), by = program]
mde[, `:=`(mde_1s_a025 = (qnorm(0.975) + qnorm(0.8)) * se_fe, mde_1s_a05 = (qnorm(0.95) + qnorm(0.8)) * se_fe,
           mde_2s_a025 = (qnorm(0.9875) + qnorm(0.8)) * se_fe,
           mde_1s_a025_I2_50 = (qnorm(0.975) + qnorm(0.8)) * se_fe * sqrt(2))]
fwrite(rs, file.path(out, "resid_sd_by_cohort.tsv"), sep = "\t")
fwrite(mde, file.path(out, "mde.tsv"), sep = "\t")
fwrite(gt, file.path(out, "gtex_predicted_sex_shift.tsv"), sep = "\t")
fwrite(pcr, file.path(out, "positive_control_sexscore.tsv"), sep = "\t")
fwrite(recheck, file.path(out, "score_recheck.tsv"), sep = "\t")
print(rs); print(mde); print(gt); print(pcr); print(recheck)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
