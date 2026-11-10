# DIG1 job 8 (exploratory; follows sections 1-3 of RESULT_DIG1.md). Do the single-cell hepatocyte
# patterns appear at protein level in an independent cohort (PXD051911 liver DIA-MS, 46 MASLD and
# 12 control participants)? Sets are fixed from pathway membership before reading protein values:
#   sterol  = 19 cholesterol-biosynthesis enzymes (HMGCS1 ... DHCR24), mevalonate + post-squalene
#   methyl  = GNMT, MAT1A, BHMT (lower per hepatocyte in both single-cell cohorts)
#   biliary = KRT7, KRT19, ANXA4 (proteins detected in >= 90% of samples among the plastic markers)
# Primary family (Holm over 3): T1 sterol score MASLD vs control; T2 within MASLD partial Spearman
# GNMT vs sterol score; T3 within MASLD partial Spearman GNMT vs biliary score. Adjusted for
# acquisition batch, age, BMI, sex (as the existing audit). Participant is the unit.
# Positive control: sterol enzymes must covary with each other (median pairwise partial rho > 0.2).
# Usage: Rscript 08_protein_links.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
STEROL <- c("HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","TM7SF2",
            "MSMO1","NSDHL","HSD17B7","EBP","SC5D","DHCR7","DHCR24")
METHYL <- c("GNMT","MAT1A","BHMT"); BILIARY <- c("KRT7","KRT19","ANXA4")
EXTRA <- c("AHCY","GAMT","PEMT","NNMT","CYP2E1","CYP1A2","AKR1B10","KRT8","KRT18","KRT23","EPCAM","SCD","FASN","LDLR","FMO1")
raw <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt"))
md <- fread(file.path(ROOT, "data/PXD051911/meta_data.txt"))
samples <- setdiff(names(raw), c("ProteinAccessions","Genes","ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% samples], by = "liver_proteomics_filename")
stopifnot(nrow(md) == 58L)
md[, `:=`(sample_id = liver_proteomics_filename, masld = as.integer(saf_diagnosis != "No_MASLD"),
          batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")),
          age = as.numeric(alder), bmi_num = as.numeric(bmi), sex = factor(gender),
          Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)), NAS = as.numeric(nafld_activity_score))]
setorder(md, sample_id)
raw <- raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
x <- as.matrix(raw[, md$sample_id, with = FALSE]); storage.mode(x) <- "double"
x[!is.finite(x) | x <= 0] <- NA_real_; x <- log2(x); rownames(x) <- raw$Genes
x <- x[rowMeans(is.na(x)) <= .5, , drop = FALSE]
x <- normalizeBetweenArrays(x, method = "quantile")
idx <- split(seq_len(nrow(x)), rownames(x))
gene <- t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x))))
colnames(gene) <- md$sample_id
present <- function(gs) intersect(gs, rownames(gene))
fwrite(data.table(set = c(rep("sterol", length(STEROL)), rep("methyl", 3), rep("biliary", 3)), gene = c(STEROL, METHYL, BILIARY),
                  quantified = c(STEROL, METHYL, BILIARY) %in% rownames(gene)), file.path(out, "sets_used.tsv"), sep = "\t")
score <- function(gs) { z <- t(scale(t(gene[present(gs), , drop = FALSE]))); colMeans(z, na.rm = TRUE) }
md[, `:=`(sterol = score(STEROL), methyl = score(METHYL), biliary = score(BILIARY))]
for (g in present(c(STEROL, METHYL, BILIARY, EXTRA))) md[, (g) := gene[g, sample_id]]
base <- c("batch", "age", "bmi_num", "sex")
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(vars), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]])))
  r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  data.table(a = a, b = b, adjust = paste(vars, collapse = "+"), n = nrow(z), rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df))
}
grp <- function(y, z = md) {
  z <- z[complete.cases(z[, c(y, "masld", base), with = FALSE])]
  f <- lm(reformulate(c("masld", base), y), data = z); s <- summary(f)$coefficients["masld", ]
  data.table(outcome = y, n_masld = sum(z$masld), n_control = sum(1 - z$masld), beta = s[1], se = s[2], p = s[4])
}
M <- md[masld == 1]
prim <- rbind(cbind(test = "T1_sterol_masld_vs_control", grp("sterol")[, .(n = n_masld + n_control, estimate = beta, p)]),
              cbind(test = "T2_GNMT_vs_sterol_within_MASLD", pcor(M, "GNMT", "sterol", base)[, .(n, estimate = rho, p)]),
              cbind(test = "T3_GNMT_vs_biliary_within_MASLD", pcor(M, "GNMT", "biliary", base)[, .(n, estimate = rho, p)]))
prim[, p_holm := p.adjust(p, "holm")]
fwrite(prim, file.path(out, "primary_tests.tsv"), sep = "\t")
# positive control: sterol enzyme coherence within MASLD
st <- present(STEROL); pc <- rbindlist(combn(st, 2, function(v) pcor(M, v[1], v[2], base), simplify = FALSE))
fwrite(data.table(n_pairs = nrow(pc), median_rho = median(pc$rho), frac_positive = mean(pc$rho > 0)), file.path(out, "positive_control_sterol_coherence.tsv"), sep = "\t")
# secondary (descriptive): per protein group effects, fibrosis and NAS within MASLD, fibrosis-adjusted pairs
prots <- present(c(STEROL, METHYL, BILIARY, EXTRA))
sec_grp <- rbindlist(lapply(c("sterol", "methyl", "biliary", prots), grp)); sec_grp[, bh := p.adjust(p, "BH")]
fwrite(sec_grp, file.path(out, "masld_vs_control_by_protein.tsv"), sep = "\t")
sev <- rbindlist(lapply(c("sterol", "methyl", "biliary", prots), function(y) rbind(pcor(M, y, "Fibrosis", base), pcor(M, y, "NAS", base))))
fwrite(sev, file.path(out, "severity_within_MASLD.tsv"), sep = "\t")
pairs2 <- rbindlist(lapply(c("sterol", "biliary", "AKR1B10", "CYP2E1", "KRT8", "NNMT", "PEMT", st), function(y)
  rbind(pcor(M, "GNMT", y, base), pcor(M, "GNMT", y, c(base, "Fibrosis")), pcor(M, "GNMT", y, c(base, "NAS")),
        pcor(M, "methyl", y, base))))
fwrite(pairs2, file.path(out, "gnmt_pairs_within_MASLD.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(prim); print(fread(file.path(out, "positive_control_sterol_coherence.tsv")))
