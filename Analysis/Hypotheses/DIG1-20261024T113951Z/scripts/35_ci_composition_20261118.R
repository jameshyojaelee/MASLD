#!/usr/bin/env Rscript
# DIG1 job 35 (20261118 candidate, N10/N18):
# Updates over 35_ci_composition.R:
# 1. For within-MASLD (n=46) PXD051911 analysis (K1 & ballooning contrasts),
#    subset raw PXD051911 intensities (single-gene PGs) to the 46 MASLD cases FIRST,
#    filter <=50% missingness, run limma::normalizeBetweenArrays(..., method="quantile")
#    on those 46 MASLD samples alone, collapse by median, and z-score (46-MASLD means/SDs)
#    before scoring ER_TM_other, PM_TM, other_TM (MF), cytosol_nonTM (CF), composition
#    markers (COMP), CI, CII, CI_minus_CII, OXPHOS_mitoribo, and CCT.
#    (Retains 58-sample PXD051911 normalization ONLY for the disease-vs-control contrast.)
# 2. Studentize Freedman-Lane permutations in fl() by HC3 SE: T* = |b1* / SE_HC3*|.

suppressPackageStartupMessages({library(data.table); library(limma); library(org.Hs.eg.db); library(AnnotationDbi)})
args <- commandArgs(TRUE)
out <- if (length(args) >= 1L && nzchar(args[1L])) {
  args[1L]
} else {
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/DIG1-20261024T113951Z/out_35_20261118"
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261118)

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
EXT <- file.path(ROOT, "data/external")
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
COMP <- c("immune", "macrophage", "stroma", "ductular", "endothelial", "hepatocyte")

# Annotations matching 31_membrane_topology.R & 33_cct_composition.R
hpa <- fread(file.path(EXT, "hpa/proteinatlas.tsv"), select = c("Gene", "Protein class"))
tm <- unique(hpa[grepl("Predicted membrane proteins", `Protein class`, fixed = TRUE), Gene])
GO <- c(ERM = "GO:0005789", ERL = "GO:0005788", MIM = "GO:0005743", MAT = "GO:0005759", PM = "GO:0005886", CYT = "GO:0005829")
go <- as.data.table(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unname(GO), keytype = "GOALL", columns = "SYMBOL")))
go <- unique(go[!is.na(SYMBOL), .(go = GOALL, gene = SYMBOL)])
inGO <- function(k) unique(go[go == GO[[k]], gene])
ER_INS <- c("SEC61A1", "SEC61A2", "SEC61B", "SEC61G", "SEC62", "SEC63", "SSR1", "SSR2", "SSR3", "SSR4", "STT3A", "STT3B", "RPN1", "RPN2",
            "DDOST", "DAD1", "OSTC", "TMEM258", "MAGT1", "TUSC3", "KRTCAP2", "OST4", "SPCS1", "SPCS2", "SPCS3", "SEC11A", "SEC11C",
            paste0("EMC", 1:10), "MMGT1", "SRPRA", "SRPRB", "TRAM1", "GET1", "CAMLG")
is_ox <- function(g) (grepl("^NDUF[ABSCV][0-9]", g) | g %in% c("SDHA", "SDHB", "SDHC", "SDHD") | grepl("^UQCR|^CYC1$", g) |
                      grepl("^COX[4-8]", g) | grepl("^ATP5", g) | grepl("^MRP[LS][0-9]", g))
classify <- function(g) {
  t <- g %in% tm
  fcase(g %in% ER_INS, "ER_insertion", is_ox(g), "OXPHOS_mitoribo",
        t & g %in% inGO("ERM"), "ER_TM_other", t & g %in% inGO("MIM"), "mitoIM_TM_other", t & g %in% inGO("PM"), "PM_TM", t, "other_TM",
        g %in% inGO("ERL"), "ER_lumen_nonTM", g %in% inGO("MAT"), "matrix_nonTM", g %in% inGO("CYT"), "cytosol_nonTM", default = "other_nonTM")
}
MK <- list(
  immune = c("PTPRC", "LCP1", "CORO1A", "LSP1", "CD74", "HLA-DRA"),
  macrophage = c("CD68", "CD163", "MARCO", "LYZ"),
  stroma = c("COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "LUM", "DCN", "TAGLN", "ACTA2", "MYL9"),
  ductular = c("KRT7", "KRT19", "EPCAM", "ANXA4"),
  endothelial = c("PECAM1", "VWF", "STAB1", "STAB2", "LYVE1", "FCGR2B", "CLEC4G"),
  hepatocyte = c("ALB", "APOA1", "APOC3", "SERPINA1", "TF", "HP", "AHSG", "FGB")
)
CCT_SET <- c("TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8")

build_pxd051911_scores <- function(x_mat, gsym, meta_sub) {
  ok <- !is.na(gsym) & trimws(gsym) != "" & !grepl(";", gsym, fixed = TRUE)
  x_mat <- x_mat[ok, , drop = FALSE]; gsym <- trimws(gsym[ok]); x_mat[!is.finite(x_mat)] <- NA_real_
  keep <- rowMeans(is.na(x_mat)) <= 0.50
  x_mat <- normalizeBetweenArrays(x_mat[keep, , drop = FALSE], method = "quantile"); gsym <- gsym[keep]
  idx <- split(seq_len(nrow(x_mat)), gsym)
  G <- t(vapply(idx, function(i) if (length(i) == 1L) x_mat[i, ] else apply(x_mat[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x_mat))))
  colnames(G) <- colnames(x_mat)
  Z <- t(scale(t(G)))
  cl <- classify(rownames(G))
  S <- c(
    split(rownames(G), cl),
    MK,
    list(
      CCT = intersect(CCT_SET, rownames(G)),
      CI = rownames(G)[grepl("^NDUF[ABSCV][0-9]", rownames(G)) & rownames(G) != "NDUFA4"],
      CII = intersect(c("SDHA", "SDHB", "SDHC", "SDHD"), rownames(G))
    )
  )
  S <- lapply(S, intersect, rownames(G))
  sc <- as.data.table(lapply(S, function(s) colMeans(Z[s, , drop = FALSE], na.rm = TRUE)))
  sc[, `:=`(
    MF = (ER_TM_other + PM_TM + other_TM) / 3,
    CF = cytosol_nonTM,
    CI_minus_CII = CI - CII
  )]
  cbind(meta_sub, sc)
}

raw <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt"))
md <- fread(file.path(ROOT, "data/PXD051911/meta_data.txt"))
smp <- setdiff(names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% smp], by = "liver_proteomics_filename"); setorder(md, liver_proteomics_filename)
x_all <- as.matrix(raw[, md$liver_proteomics_filename, with = FALSE]); storage.mode(x_all) <- "double"
x_all[!is.finite(x_all) | x_all <= 0] <- NA_real_; x_all <- log2(x_all)
a1_all <- md[, .(
  sample = liver_proteomics_filename,
  grp = fifelse(saf_diagnosis != "No_MASLD", "disease", "control"),
  batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")),
  age = as.numeric(alder), bmi = as.numeric(bmi), sex = factor(gender),
  Ballooning = as.numeric(hepatocellular_ballooning_score),
  Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade))
)]
# 58-sample matrix for case-vs-control contrast only
D1 <- build_pxd051911_scores(x_all, raw$Genes, a1_all)

# 46-MASLD-only quantile-normalized + z-scored matrix for within-MASLD (K1 & ballooning)
is_m <- a1_all$grp == "disease"
M1 <- build_pxd051911_scores(x_all[, is_m, drop = FALSE], raw$Genes, droplevels(a1_all[is_m]))
stopifnot(nrow(M1) == 46L)
fwrite(M1, file.path(out, "scores_PXD051911_MASLD46_only_qnorm.tsv"), sep = "\t")

# External validation cohorts D2 & D3
rd_ext <- function(k) {
  a <- fread(file.path(RUN, "out_cct_composition", paste0("scores_", k, ".tsv")))
  b <- fread(file.path(RUN, "out_membrane_topology", paste0("scores_", k, ".tsv")))
  stopifnot(identical(a$sample, b$sample))
  cbind(a, b[, .(CI, CII, CI_minus_CII, OXPHOS_mitoribo)])
}
D2 <- rd_ext("PXD027722")
D3 <- rd_ext("PXD065892")

pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]])))
  r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1L
  data.table(n = nrow(z), df = df, estimate = r, p = 2 * pt(-abs(r * sqrt(df / pmax(1e-15, 1 - r^2))), df))
}

# HC3-studentized Freedman-Lane permutation test
fl_hc3 <- function(z, y, g, vars, nperm) {
  z <- z[complete.cases(z[, c(y, vars), with = FALSE])]
  lab <- as.numeric(z$grp == g)
  X0 <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z)
  X1 <- cbind(lab = lab, X0)
  yy <- z[[y]]
  XtXi <- solve(crossprod(X1)); C <- XtXi %*% t(X1)
  h <- rowSums((X1 %*% XtXi) * X1)
  w_hc3 <- C[1L, ]^2 / (1 - h)^2
  b_obs <- as.numeric(C %*% yy)
  e_obs <- yy - as.numeric(X1 %*% b_obs)
  se_obs <- sqrt(sum(w_hc3 * e_obs^2))
  T_obs <- abs(b_obs[1L] / se_obs)

  f0 <- qr.fitted(qr(X0), yy)
  e0 <- yy - f0
  chunk_size <- 5000L
  exceed <- 0L
  for (s_idx in seq(1L, nperm, by = chunk_size)) {
    B_cur <- min(chunk_size, nperm - s_idx + 1L)
    Perm <- replicate(B_cur, sample.int(length(e0)))
    Y_star <- matrix(f0, nrow = length(e0), ncol = B_cur) + matrix(e0[Perm], nrow = length(e0), ncol = B_cur)
    B_mat <- C %*% Y_star
    E_mat <- Y_star - X1 %*% B_mat
    se_star <- sqrt(as.numeric(crossprod(w_hc3, E_mat^2)))
    T_star <- abs(B_mat[1L, ] / se_star)
    exceed <- exceed + sum(T_star >= T_obs - 1e-12, na.rm = TRUE)
  }
  df_res <- nrow(X1) - ncol(X1)
  data.table(
    n = nrow(z), n1 = sum(lab), df = df_res,
    estimate = unname(b_obs[1L]), se_hc3 = se_obs, t_hc3 = abs(b_obs[1L] / se_obs),
    p_hc3_asym = 2 * pt(-abs(b_obs[1L] / se_obs), df = df_res),
    sd_y = sd(yy), p = (1 + exceed) / (nperm + 1)
  )
}

BASE1 <- c("batch", "age", "bmi", "sex")
ADJ <- list(none = character(), MF_CF = c("MF", "CF"), composition = c("MF", "CF", COMP), composition_fibrosis = c("MF", "CF", COMP, "Fibrosis"))
prim <- rbind(
  cbind(id = "K1", what = "PXD051911 MASLD (46-only qnorm): CI-CII vs ballooning | cov + MF + CF + composition",
        pcor(M1, "CI_minus_CII", "Ballooning", c(BASE1, ADJ$composition))),
  cbind(id = "K2", what = "PXD065892 NAFLD vs healthy: CI-CII | MF + CF + composition (HC3 FL)",
        fl_hc3(D3[grp %in% c("disease", "control")], "CI_minus_CII", "disease", ADJ$composition, 20000L)),
  cbind(id = "K3", what = "PXD065892 AIH vs healthy: CI-CII | MF + CF + composition (HC3 FL)",
        fl_hc3(D3[grp %in% c("AIH", "control")], "CI_minus_CII", "AIH", ADJ$composition, 20000L)),
  fill = TRUE
)
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")

Y <- c("CI_minus_CII", "CI", "CII", "OXPHOS_mitoribo")
CON <- list(
  list("PXD051911", D1, "disease", BASE1),
  list("PXD027722", D2, "disease", character()),
  list("PXD027722", D2, "cirrhosis", character()),
  list("PXD065892", D3, "disease", character()),
  list("PXD065892", D3, "AIH", character())
)
desc <- rbindlist(c(
  lapply(names(ADJ), function(a) rbindlist(lapply(Y, function(y)
    cbind(analysis = "PXD051911 MASLD (46-only qnorm) vs ballooning", y = y, adjust = a, pcor(M1, y, "Ballooning", c(BASE1, ADJ[[a]])))))),
  lapply(CON, function(k) rbindlist(lapply(names(ADJ)[1:3], function(a) rbindlist(lapply(Y, function(y)
    cbind(analysis = paste(k[[1]], k[[3]], "vs control"), y = y, adjust = a, fl_hc3(k[[2]][grp %in% c(k[[3]], "control")], y, k[[3]], c(k[[4]], ADJ[[a]]), 5000L)))))))
), fill = TRUE)
desc[, bh := p.adjust(p, "BH")][, sd_units := estimate / sd_y]
fwrite(desc, file.path(out, "adjusted_contrasts.tsv"), sep = "\t")

cc <- rbindlist(lapply(c(COMP, "MF", "CF"), function(y) cbind(y = y, pcor(M1, "CI_minus_CII", y, BASE1))))
fwrite(cc, file.path(out, "PXD051911_MASLD_CIminusCII_vs_composition.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220)
print(prim, digits = 3)
cat("DIG1_35_20261118 DONE\n")
