# DIG1 job 35 (exploratory; follows jobs 32-33). Job 33 showed that most of the CCT rise goes with
# cell-composition markers. Does the Complex I deficit (CI - CII) survive the same adjustment?
# Inputs: job 31 scores (CI, CII, CI_minus_CII, OXPHOS_mitoribo) and job 33 scores (MF, CF, six
# composition scores: immune, macrophage, stroma, ductular, endothelial, hepatocyte).
# Primary family (Holm over 3), adjustment = MF + CF + six composition scores (+ batch, age, BMI, sex
# in PXD051911):
#   K1  PXD051911 MASLD (n 46): CI - CII vs ballooning, partial Spearman.
#   K2  PXD065892 NAFLD vs healthy: CI - CII, Freedman-Lane, 20,000.
#   K3  PXD065892 AIH vs healthy: same.
# Descriptive (BH): CI, CII, CI - CII, OXPHOS_mitoribo under none / MF_CF / composition /
# composition + fibrosis (ballooning only); all five group contrasts. Seed 20261024.
# Usage: Rscript 35_ci_composition.R <new_out_dir>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z"
COMP <- c("immune", "macrophage", "stroma", "ductular", "endothelial", "hepatocyte")
rd <- function(k) { a <- fread(file.path(RUN, "out_cct_composition", paste0("scores_", k, ".tsv")))
  b <- fread(file.path(RUN, "out_membrane_topology", paste0("scores_", k, ".tsv")))
  stopifnot(identical(a$sample, b$sample)); cbind(a, b[, .(CI, CII, CI_minus_CII, OXPHOS_mitoribo)]) }
D1 <- rd("PXD051911"); D2 <- rd("PXD027722"); D3 <- rd("PXD065892")
D1[, sex := factor(sex)][, batch := factor(batch)]
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(n = nrow(z), df = df, estimate = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
fl <- function(z, y, g, vars, nperm) {
  z <- z[complete.cases(z[, c(y, vars), with = FALSE])]; lab <- as.numeric(z$grp == g)
  X0 <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z); X1 <- cbind(lab, X0)
  b1 <- function(yy) qr.coef(qr(X1), yy)[1]
  f0 <- qr.fitted(qr(X0), z[[y]]); e0 <- z[[y]] - f0; obs <- b1(z[[y]])
  nd <- replicate(nperm, b1(f0 + e0[sample.int(length(e0))]))
  data.table(n = nrow(z), n1 = sum(lab), df = nrow(z) - ncol(X1), estimate = unname(obs), sd_y = sd(z[[y]]), p = (1 + sum(abs(nd) >= abs(obs) - 1e-12)) / (nperm + 1)) }
BASE1 <- c("batch", "age", "bmi", "sex")
ADJ <- list(none = character(), MF_CF = c("MF", "CF"), composition = c("MF", "CF", COMP), composition_fibrosis = c("MF", "CF", COMP, "Fibrosis"))
M1 <- D1[grp == "disease"]
prim <- rbind(cbind(id = "K1", what = "PXD051911 MASLD: CI-CII vs ballooning | cov + MF + CF + composition", pcor(M1, "CI_minus_CII", "Ballooning", c(BASE1, ADJ$composition))),
              cbind(id = "K2", what = "PXD065892 NAFLD vs healthy: CI-CII | MF + CF + composition", fl(D3[grp %in% c("disease", "control")], "CI_minus_CII", "disease", ADJ$composition, 20000)),
              cbind(id = "K3", what = "PXD065892 AIH vs healthy: CI-CII | MF + CF + composition", fl(D3[grp %in% c("AIH", "control")], "CI_minus_CII", "AIH", ADJ$composition, 20000)), fill = TRUE)
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
Y <- c("CI_minus_CII", "CI", "CII", "OXPHOS_mitoribo")
CON <- list(list("PXD051911", D1, "disease", BASE1), list("PXD027722", D2, "disease", character()), list("PXD027722", D2, "cirrhosis", character()),
            list("PXD065892", D3, "disease", character()), list("PXD065892", D3, "AIH", character()))
desc <- rbindlist(c(
  lapply(names(ADJ), function(a) rbindlist(lapply(Y, function(y)
    cbind(analysis = "PXD051911 MASLD vs ballooning", y = y, adjust = a, pcor(M1, y, "Ballooning", c(BASE1, ADJ[[a]])))))),
  lapply(CON, function(k) rbindlist(lapply(names(ADJ)[1:3], function(a) rbindlist(lapply(Y, function(y)
    cbind(analysis = paste(k[[1]], k[[3]], "vs control"), y = y, adjust = a, fl(k[[2]][grp %in% c(k[[3]], "control")], y, k[[3]], c(k[[4]], ADJ[[a]]), 5000)))))))), fill = TRUE)
desc[, bh := p.adjust(p, "BH")][, sd_units := estimate / sd_y]
fwrite(desc, file.path(out, "adjusted_contrasts.tsv"), sep = "\t")
cc <- rbindlist(lapply(c(COMP, "MF", "CF"), function(y) cbind(y = y, pcor(M1, "CI_minus_CII", y, BASE1))))
fwrite(cc, file.path(out, "PXD051911_MASLD_CIminusCII_vs_composition.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220)
print(prim, digits = 3)
print(dcast(desc, analysis + y ~ adjust, value.var = "estimate"), digits = 2)
print(dcast(desc, analysis + y ~ adjust, value.var = "p"), digits = 2)
print(dcast(desc[!grepl("ballooning", analysis)], analysis + y ~ adjust, value.var = "sd_units"), digits = 2)
print(cc, digits = 2)
