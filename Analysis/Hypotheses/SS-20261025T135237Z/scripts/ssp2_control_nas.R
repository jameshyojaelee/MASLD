# SS SSP2 (PRESPEC_SS_addendum_SSP2_control_nas.md): disease vs control and NAS adjacent-step replication in SteatoSITE,
# with control status and NAS components from the LiRNA atlas (unlabelled columns, names inferred).
# Usage:
#   Rscript ssp2_control_nas.R draws <out>          # units, kept genes, permutation and bootstrap draws
#   Rscript ssp2_control_nas.R fit <out> <i> <k>    # fit tasks whose index mod k == i (0-based)
#   Rscript ssp2_control_nas.R summarize <out>
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_common.R")
a <- commandArgs(trailingOnly = TRUE); mode <- a[1]; OUT <- a[2]
S2 <- file.path(SS, "out_s2_qc")
STAGE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv"
LIRNA <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/downloads_untrusted/metadata_search_20261028T005409Z/lirna/steatosite_LiRNA_meta_by_sample.tsv"
C3 <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv"
B_PERM <- 1000L; B_BOOT <- 200L
NAS_LV <- c("NAS0", "NAS1_2", "NAS3_4", "NAS5_8")
ADJ <- c(NAS1_2_vs_NAS0 = "NAS1_2", NAS3_4_vs_NAS1_2 = "NAS3_4", NAS5_8_vs_NAS3_4 = "NAS5_8")

adj_coefs <- function(cm) {  # adjacent NAS contrasts from NAS0-referenced coefficients
  g <- function(l) if (l == "NAS0") rep(0, nrow(cm)) else cm[, paste0("nas_bin", l)]
  cbind(NAS1_2_vs_NAS0 = g("NAS1_2"), NAS3_4_vs_NAS1_2 = g("NAS3_4") - g("NAS1_2"), NAS5_8_vs_NAS3_4 = g("NAS5_8") - g("NAS3_4"))
}

if (mode == "draws") {
  OUT <- new_outdir(OUT)
  units <- fread(file.path(S2, "units.tsv"))[include == TRUE & !is.na(inferred_sex)]
  st <- read_stage(STAGE); L <- fread(LIRNA, na.strings = c("NA", ""))
  u <- merge(merge(units, st, by = "sample_alias", all.x = TRUE), L[, .(sample_alias, ctrl = Info_2, steat = Info_7, infl = Info_4, ball = Info_3)],
             by = "sample_alias")
  setorder(u, sample_alias)
  u[, group := fifelse(ctrl == "yes", "control", fifelse(ctrl == "no" & !is.na(stage_code), "masld", NA_character_))]
  u[, nas := steat + infl + ball]
  u[ctrl == "yes", nas := 0L]
  u[, nas_bin := fifelse(nas == 0, "NAS0", fifelse(nas <= 2, "NAS1_2", fifelse(nas <= 4, "NAS3_4", "NAS5_8")))]
  u[!(ctrl %in% c("yes", "no")), nas_bin := NA_character_]
  UD <- u[!is.na(group)]; UD[, group := factor(group, levels = c("control", "masld"))]
  UN <- u[!is.na(nas_bin)]; UN[, nas_bin := factor(nas_bin, levels = NAS_LV)]
  dge <- readRDS(file.path(S2, "dge_included.rds"))
  if (!"group" %in% names(dge$samples)) dge$samples <- cbind(group = factor(rep(1L, ncol(dge))), dge$samples)
  KD <- rownames(dge)[filterByExpr(dge[, UD$sample_alias], group = UD$group)]
  KN <- rownames(dge)[filterByExpr(dge[, UN$sample_alias], group = UN$nas_bin)]
  set.seed(20261028L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  perm_within <- function(lab, strata) { for (s in sort(unique(strata))) { i <- which(strata == s); lab[i] <- lab[i][sample.int(length(i))] }; lab }
  sD <- paste(UD$batch, UD$inferred_sex); sN <- paste(UN$batch, UN$inferred_sex)
  PD <- t(vapply(seq_len(B_PERM), function(k) as.character(perm_within(as.character(UD$group), sD)), character(nrow(UD))))
  PN <- t(vapply(seq_len(B_PERM), function(k) as.character(perm_within(as.character(UN$nas_bin), sN)), character(nrow(UN))))
  set.seed(20261029L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  BD <- t(vapply(seq_len(B_BOOT), function(k) boot_draw(paste(UD$group, UD$inferred_sex)), integer(nrow(UD))))
  BN <- t(vapply(seq_len(B_BOOT), function(k) boot_draw(paste(UN$nas_bin, UN$inferred_sex)), integer(nrow(UN))))
  tasks <- c("obsD", "obsN", "qualD", "qualN", paste0("pD", seq_len(B_PERM)), paste0("pN", seq_len(B_PERM)),
             paste0("bD", seq_len(B_BOOT)), paste0("bN", seq_len(B_BOOT)))
  saveRDS(list(UD = UD, UN = UN, KD = KD, KN = KN, PD = PD, PN = PN, BD = BD, BN = BN, tasks = tasks), file.path(OUT, "draws.rds"))
  writeLines(c(paste("DvC units:", nrow(UD), "controls", sum(UD$group == "control"), "masld", sum(UD$group == "masld")),
               paste("NAS units:", nrow(UN), paste(names(table(UN$nas_bin)), table(UN$nas_bin), collapse = " ")),
               paste("kept genes DvC", length(KD), "NAS", length(KN)), paste("tasks", length(tasks))), file.path(OUT, "draw_notes.txt"))
  log_msg("draws written"); quit(save = "no")
}

Z <- readRDS(file.path(OUT, "draws.rds"))
if (mode == "fit") {
  i <- as.integer(a[3]); k <- as.integer(a[4])
  dge <- readRDS(file.path(S2, "dge_included.rds"))
  if (!"group" %in% names(dge$samples)) dge$samples <- cbind(group = factor(rep(1L, ncol(dge))), dge$samples)
  dD <- calcNormFactors(dge[Z$KD, Z$UD$sample_alias, keep.lib.sizes = TRUE], method = "RLE")
  dN <- calcNormFactors(dge[Z$KN, Z$UN$sample_alias, keep.lib.sizes = TRUE], method = "RLE")
  dir.create(file.path(OUT, "fits"), showWarnings = FALSE)
  todo <- Z$tasks[(seq_along(Z$tasks) - 1L) %% k == i]
  fitD <- function(info, d, extra = character(), eb = FALSE) fit_qw(d, design_of(info, "group", extra), eb)
  fitN <- function(info, d, extra = character(), eb = FALSE) fit_qw(d, design_of(info, "nas_bin", extra), eb)
  res <- mclapply(todo, function(t) {
    f <- file.path(OUT, "fits", paste0(t, ".rds")); if (file.exists(f)) return(TRUE)
    r <- if (t == "obsD") { fit <- fitD(Z$UD, dD, eb = TRUE); topTable(fit, coef = "groupmasld", number = Inf, sort.by = "none")
    } else if (t == "obsN") { fit <- fitN(Z$UN, dN, eb = TRUE)
      cn <- colnames(fit$design); cm <- matrix(0, length(cn), 3, dimnames = list(cn, names(ADJ)))
      cm["nas_binNAS1_2", 1] <- 1; cm[c("nas_binNAS3_4", "nas_binNAS1_2"), 2] <- c(1, -1); cm[c("nas_binNAS5_8", "nas_binNAS3_4"), 3] <- c(1, -1)
      f2 <- eBayes(contrasts.fit(fit, cm)); setNames(lapply(1:3, function(j) topTable(f2, coef = j, number = Inf, sort.by = "none")), names(ADJ))
    } else if (t == "qualD") fitD(Z$UD, dD, QUAL)$coefficients[, "groupmasld"]
    else if (t == "qualN") adj_coefs(fitN(Z$UN, dN, QUAL)$coefficients)
    else { j <- as.integer(sub("^[pb][DN]", "", t)); D <- substr(t, 2, 2) == "D"
      if (startsWith(t, "p")) { info <- if (D) copy(Z$UD) else copy(Z$UN)
        if (D) info[, group := factor(Z$PD[j, ], levels = c("control", "masld"))] else info[, nas_bin := factor(Z$PN[j, ], levels = NAS_LV)]
        if (D) fitD(info, dD)$coefficients[, "groupmasld"] else adj_coefs(fitN(info, dN)$coefficients)
      } else { idx <- if (D) Z$BD[j, ] else Z$BN[j, ]; info <- if (D) Z$UD[idx] else Z$UN[idx]
        d <- if (D) dD[, idx] else dN[, idx]; colnames(d) <- make.unique(colnames(d))
        if (D) fitD(info, d)$coefficients[, "groupmasld"] else adj_coefs(fitN(info, d)$coefficients) } }
    saveRDS(r, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f); TRUE }, mc.cores = NCORES, mc.preschedule = FALSE)
  if (!all(vapply(res, isTRUE, TRUE))) { print(res[!vapply(res, isTRUE, TRUE)][1]); stop("fit failed") }
  log_msg("worker", i, "of", k, "done:", length(todo)); quit(save = "no")
}

if (mode == "summarize") {
  rd <- function(t) readRDS(file.path(OUT, "fits", paste0(t, ".rds")))
  rel <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/frozen/release/deg_results.csv")
  canon <- rel[padj < 0.05 & abs(logFC) > 0.5]; stopifnot(nrow(canon) == 1854L)
  se3 <- fread(file.path(FROZEN, "stage/stage_extension_all_gene_results.tsv"))[axis == "NAS_adjacent"]
  sets <- c(list(DvC = list(g = intersect(canon$gene, Z$KD), x = function(g) canon$logFC[match(g, canon$gene)])),
            lapply(setNames(names(ADJ), names(ADJ)), function(cn) { h <- se3[contrast == cn & FDR < 0.05]
              list(g = intersect(h$gene_id_versioned, Z$KN), x = function(g) h$logFC[match(g, h$gene_id_versioned)]) }))
  yv <- function(obj, s, g) if (s == "DvC") obj[g] else obj[g, s]
  oD <- rd("obsD"); oN <- rd("obsN")
  yobs <- function(s, g) if (s == "DvC") oD[g, "logFC"] else oN[[s]][g, "logFC"]
  stat <- function(y, x) c(b = ols_slope(y, x), conc = conc(y, x))
  R <- rbindlist(lapply(names(sets), function(s) { g <- sets[[s]]$g; x <- sets[[s]]$x(g); D <- s == "DvC"
    if (length(g) < 2L) return(data.table(test = s, n_genes = length(g), testable = FALSE))
    ob <- stat(yobs(s, g), x)
    nul <- t(vapply(seq_len(B_PERM), function(j) stat(yv(rd(paste0(if (D) "pD" else "pN", j)), s, g), x), numeric(2)))
    bt <- vapply(seq_len(B_BOOT), function(j) ols_slope(yv(rd(paste0(if (D) "bD" else "bN", j)), s, g), x), 0)
    q <- ols_slope(yv(rd(if (D) "qualD" else "qualN"), s, g), x)
    up <- x > 0; yo <- yobs(s, g)
    data.table(test = s, n_genes = length(g), b = ob[["b"]], conc = ob[["conc"]], null_mean_b = mean(nul[, 1]), null_sd_b = sd(nul[, 1]),
               null_conc_mean = mean(nul[, 2]), p_one = (1 + sum(nul[, 1] >= ob[["b"]])) / (1 + B_PERM),
               testable = 2.80 * sd(nul[, 1]) <= 0.30 && length(g) >= 10L, boot_se = sd(bt),
               ci_lo = ob[["b"]] - 1.96 * sd(bt), ci_hi = ob[["b"]] + 1.96 * sd(bt), b_quality = q,
               b_up = if (sum(up) > 1) ols_slope(yo[up], x[up]) else NA_real_, b_down = if (sum(!up) > 1) ols_slope(yo[!up], x[!up]) else NA_real_) }), fill = TRUE)
  fam <- c("DvC", "NAS3_4_vs_NAS1_2", "NAS5_8_vs_NAS3_4")
  R[, family := fifelse(test %in% fam, "SSP2-F1", "secondary")]
  R[family == "SSP2-F1", p_holm := p.adjust(fifelse(testable, p_one, 1), "holm", n = 3L)]
  R[, nc_perm_ok := abs(null_mean_b) <= 0.05][, quality_sensitive := b_quality < 0.5 * b]
  top <- canon[gene %in% Z$KD][order(-abs(t))][1:10]
  PC <- data.table(gene = top$gene, symbol = top$symbol, disc = top$logFC, ss = oD[top$gene, "logFC"])[, agree := sign(disc) == sign(ss)]
  ssc <- data.table(gene = rownames(oD), logFC = oD$logFC, q = oD$adj.P.Val)[q < 0.05]
  c3 <- fread(C3)
  OV <- data.table(set = "pooled_1854", n_set = nrow(canon), n_ss_calls = nrow(ssc), shared = length(intersect(ssc$gene, canon$gene)),
                   sign_agree_shared = mean(sign(ssc$logFC[match(intersect(ssc$gene, canon$gene), ssc$gene)]) == sign(canon$logFC[match(intersect(ssc$gene, canon$gene), canon$gene)])))
  fwrite(R, file.path(OUT, "ssp2_family.tsv"), sep = "\t"); fwrite(PC, file.path(OUT, "pc_dvc.tsv"), sep = "\t")
  fwrite(ssc, file.path(OUT, "steatosite_dvc_calls.tsv"), sep = "\t"); fwrite(OV, file.path(OUT, "overlap_pooled.tsv"), sep = "\t")
  saveRDS(oD, file.path(OUT, "obs_dvc_toptable.rds")); saveRDS(oN, file.path(OUT, "obs_nas_toptables.rds"))
  writeLines(c(paste("PC-DvC agree:", sum(PC$agree), "of", nrow(PC), "pass:", sum(PC$agree) >= 8), paste("cohort table columns:", paste(names(c3), collapse = ","))),
             file.path(OUT, "summary_notes.txt"))
  writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  print(R); print(PC); print(OV); log_msg("summary written")
}
