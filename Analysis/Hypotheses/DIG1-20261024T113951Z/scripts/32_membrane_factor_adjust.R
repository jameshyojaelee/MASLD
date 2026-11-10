# DIG1 job 32 (exploratory; follows job 31). Job 31 found a strong sample-level factor shared by
# non-mitochondrial integral membrane proteins: within 46 MASLD (PXD051911) the ER insertion, other
# ER TM and plasma-membrane TM scores correlate 0.85-0.92, and the disease shift of these scores
# changes sign between cohorts. Do the earlier protein findings survive adjustment for it?
# Membrane factor MF = mean of the ER_TM_other, PM_TM and other_TM class scores (job 31, mean z;
# no mitochondrial, OXPHOS or ER-insertion proteins). Cytosol factor CF = cytosol_nonTM score.
# Inputs: out_membrane_topology/scores_*.tsv (job 31).
# Primary family (Holm over 3), the three earlier positive CI - CII results, now + MF:
#   A1  PXD051911 MASLD (n 46): CI - CII vs ballooning, partial Spearman; batch, age, BMI, sex, MF.
#   A2  PXD065892 NAFLD vs healthy: CI - CII, lm(score ~ group + MF), Freedman-Lane, 20,000.
#   A3  PXD065892 AIH vs healthy: same.
# Descriptive (BH over the table): CI, CII, CI - CII, ER_insertion, CCT, CCT - CF, ER_chaperones,
# OXPHOS_mitoribo for all five contrasts and for ballooning, each unadjusted, + MF, and + MF + CF.
# PXD051911 contrasts keep batch, age, BMI, sex. Seed 20261024.
# Usage: Rscript 32_membrane_factor_adjust.R <new_out_dir>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z"
rd <- function(k) { d <- fread(file.path(RUN, "out_membrane_topology", paste0("scores_", k, ".tsv")))
  d[, MF := (ER_TM_other + PM_TM + other_TM) / 3][, CF := cytosol_nonTM][, CCT_minus_CF := CCT - CF][] }
D1 <- rd("PXD051911"); D2 <- rd("PXD027722"); D3 <- rd("PXD065892")
D1[, sex := factor(sex)][, batch := factor(batch)]
Y <- c("CI_minus_CII", "CI", "CII", "OXPHOS_mitoribo", "ER_insertion", "CCT", "CCT_minus_CF", "ER_chaperones", "MF", "CF")
ADJ <- list(none = character(), MF = "MF", MF_CF = c("MF", "CF"))
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(n = nrow(z), estimate = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
fl <- function(z, y, g, vars, nperm) {  # Freedman-Lane permutation for the group coefficient
  z <- z[complete.cases(z[, c(y, vars), with = FALSE])]; lab <- as.numeric(z$grp == g)
  X0 <- model.matrix(reformulate(if (length(vars)) vars else "1"), data = z); X1 <- cbind(lab, X0)
  b1 <- function(yy) qr.coef(qr(X1), yy)[1]
  f0 <- qr.fitted(qr(X0), z[[y]]); e0 <- z[[y]] - f0; obs <- b1(z[[y]])
  nd <- replicate(nperm, b1(f0 + e0[sample.int(length(e0))]))
  data.table(n = nrow(z), n1 = sum(lab), estimate = unname(obs), sd_y = sd(z[[y]]), p = (1 + sum(abs(nd) >= abs(obs) - 1e-12)) / (nperm + 1)) }
BASE1 <- c("batch", "age", "bmi", "sex")
M1 <- D1[grp == "disease"]
CON <- list(list("PXD051911", D1, "disease", BASE1), list("PXD027722", D2, "disease", character()), list("PXD027722", D2, "cirrhosis", character()),
            list("PXD065892", D3, "disease", character()), list("PXD065892", D3, "AIH", character()))
prim <- rbind(cbind(id = "A1", what = "PXD051911 MASLD: CI-CII vs ballooning | cov + MF", pcor(M1, "CI_minus_CII", "Ballooning", c(BASE1, "MF")), n1 = NA_integer_, sd_y = NA_real_),
              cbind(id = "A2", what = "PXD065892 NAFLD vs healthy: CI-CII | MF", fl(D3[grp %in% c("disease", "control")], "CI_minus_CII", "disease", "MF", 20000)),
              cbind(id = "A3", what = "PXD065892 AIH vs healthy: CI-CII | MF", fl(D3[grp %in% c("AIH", "control")], "CI_minus_CII", "AIH", "MF", 20000)), fill = TRUE)
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
desc <- rbindlist(c(
  lapply(names(ADJ), function(a) rbindlist(lapply(setdiff(Y, ADJ[[a]]), function(y)
    cbind(analysis = "PXD051911 MASLD vs ballooning (partial Spearman)", y = y, adjust = a, pcor(M1, y, "Ballooning", c(BASE1, ADJ[[a]])))))),
  lapply(CON, function(k) rbindlist(lapply(names(ADJ), function(a) rbindlist(lapply(setdiff(Y, ADJ[[a]]), function(y)
    cbind(analysis = paste(k[[1]], k[[3]], "vs control (coef, score units)"), y = y, adjust = a,
          fl(k[[2]][grp %in% c(k[[3]], "control")], y, k[[3]], c(k[[4]], ADJ[[a]]), 5000)))))))), fill = TRUE)
desc[, bh := p.adjust(p, "BH")][, sd_units := estimate / sd_y]
fwrite(desc, file.path(out, "adjusted_contrasts.tsv"), sep = "\t")
fwrite(rbindlist(lapply(list(PXD051911 = D1, PXD027722 = D2, PXD065892 = D3), function(d)
  data.table(pair = c("MF~CI", "MF~CII", "MF~CF", "CI~CF", "CCT~CF"),
             rho = c(cor(d$MF, d$CI, method = "spearman"), cor(d$MF, d$CII, method = "spearman"), cor(d$MF, d$CF, method = "spearman"),
                     cor(d$CI, d$CF, method = "spearman"), cor(d$CCT, d$CF, method = "spearman")))), idcol = "cohort"),
  file.path(out, "factor_correlations_all_samples.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220)
print(prim, digits = 3)
print(dcast(desc[y %in% c("CI_minus_CII", "CI", "CII", "ER_insertion", "CCT", "CCT_minus_CF", "MF")], analysis + y ~ adjust, value.var = "estimate"), digits = 2)
print(dcast(desc[y %in% c("CI_minus_CII", "CI", "CII", "ER_insertion", "CCT", "CCT_minus_CF", "MF")], analysis + y ~ adjust, value.var = "p"), digits = 2)
print(fread(file.path(out, "factor_correlations_all_samples.tsv")), digits = 2)
