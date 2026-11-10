# DIG1 job 68b (post hoc to job 68, whose results were read: the state protein score followed ballooning, but the top
# proteins were PALLD, ANXA2, PLOD2, S100A10, CAP2, MSN; several are made by stromal cells too). Same tests, restricted
# to hepatocyte-made proteins: genes whose median hepatocyte UMI share over MASLD donors (job 68a) is >= 0.85 in every
# cohort that measures them (Flex, GSE202379, GSE244832). The matched null pool is restricted the same way.
# Job 68 header follows.
# DIG1 job 68 (exploratory; post hoc to sections 34-41). Is the inflamed-hepatocyte (CCL20, NF-kB) state visible
# in liver protein, and does it follow ballooning there as the bulk RNA score does (section 38)?
# Read before writing: only which state genes are present in the PXD051911 protein matrix (KRT23, NFKB2, SOD2,
# TNFAIP2, TNIP1, ICAM1, LCN2 present; CCL20, CXCL10, IL32, BIRC3, RELB, TRIM31, FCAMR absent). No protein value read.
# Set: CCL20_STATE_UP of job 67b (rises > 0.25 in CCL20+ hepatocytes in all four single-cell cohorts, depth mode).
# Data: PXD051911 (46 MASLD, 12 controls), processing as job 21 (single-gene groups, <= 50% missing, quantile
# normalised, median per gene). Score: z within protein across the 58 participants, mean over set proteins.
# Null: 2,000 random protein sets of the same size, matched on mean abundance (10 bins) x missing fraction (3 bins),
# drawn from proteins outside the set; seed 20261025. Empirical two-sided p of the statistic against the null,
# (1 + #|null - mean| >= |obs - mean|) / (1 + B).
# Primary (one test, named here): within MASLD, partial Spearman of the state protein score with ballooning grade,
#   adjusted for batch, age, BMI and sex, against the matched-set null.
# Secondary (descriptive, BH over the table): steatosis, lobular inflammation, fibrosis, NAS; MASLD vs control
#   (OLS beta with the same covariates); the ballooning term adjusted for an immune protein score (PTPRC, CD74,
#   HLA-DRA, LCP1, CORO1A, LYZ, CD14, CD163, CD68, MNDA) and for a hepatocyte-mass score (ALB, APOA1, CYP2E1, HP);
#   single proteins KRT23, NFKB2, TNFAIP2, TNIP1, SOD2, ICAM1, LCN2, SAA1, SAA2.
# Replication (descriptive): PXD027722 NASH, cirrhosis vs control; PXD065892 NAFLD, AIH vs healthy (job 27
#   processing; no covariates); state score difference with 5,000 label shuffles and the matched-set excess.
# Unit: one liver sample per participant. Cross-sectional.
# Usage: Rscript 68b_state_protein_hep.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261025)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
EXT <- file.path(ROOT, "data/external")
# state set
g53 <- fread(file.path(RUN, "out_chemokine_hepatocyte_state/gene_lfc_matched.tsv.gz"))[donors == "masld"][, donors := NULL]
g56 <- fread(file.path(RUN, "out_atlas_state_replication/gene_lfc_matched.tsv.gz"))
SC <- rbind(g53, g56, use.names = TRUE)[mode == "depth" & is.finite(mean_excess) & anchor == "CCL20"]
xs <- SC[, .(k = .N, mn = min(mean_excess)), by = gene][k == 4]
STATE <- setdiff(xs[mn > 0.25, gene], c("FCAMR", "TRIM31", "UBD", "CXCL10", "IL32", "CCL20", "HSD17B13", "PNPLA3"))
# PXD051911, as job 21
raw <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt"))
md <- fread(file.path(ROOT, "data/PXD051911/meta_data.txt"))
samples <- setdiff(names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% samples], by = "liver_proteomics_filename"); stopifnot(nrow(md) == 58L)
md[, `:=`(sample_id = liver_proteomics_filename, masld = as.integer(saf_diagnosis != "No_MASLD"),
          batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")),
          age = as.numeric(alder), bmi_num = as.numeric(bmi), sex = factor(gender),
          Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)), NAS = as.numeric(nafld_activity_score),
          Steatosis = as.numeric(steatosis_score), Ballooning = as.numeric(hepatocellular_ballooning_score),
          Inflammation = as.numeric(lobular_inflammation_score))]
setorder(md, sample_id)
raw <- raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
x <- as.matrix(raw[, md$sample_id, with = FALSE]); storage.mode(x) <- "double"
x[!is.finite(x) | x <= 0] <- NA_real_; x <- log2(x); rownames(x) <- raw$Genes
x <- x[rowMeans(is.na(x)) <= .5, , drop = FALSE]
x <- normalizeBetweenArrays(x, method = "quantile")
idx <- split(seq_len(nrow(x)), rownames(x))
gene <- t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x))))
colnames(gene) <- md$sample_id
Z <- t(scale(t(gene)))
HS <- fread(file.path(RUN, "out_hep_share/gene_hep_share.tsv.gz"))[is.finite(median_hep_share)]
HEPOK <- HS[, .(mn = min(median_hep_share), k = .N), by = gene][mn >= 0.85, gene]
S_all <- intersect(STATE, rownames(gene)); S_in <- intersect(S_all, HEPOK)
fwrite(data.table(protein = S_all, hep_made = S_all %in% HEPOK), file.path(out, "state_proteins_hep_filter.tsv"), sep = "\t")
IMM <- intersect(c("PTPRC", "CD74", "HLA-DRA", "LCP1", "CORO1A", "LYZ", "CD14", "CD163", "CD68", "MNDA"), rownames(gene))
HEP <- intersect(c("ALB", "APOA1", "CYP2E1", "HP"), rownames(gene))
SINGLE <- intersect(c("KRT23", "NFKB2", "TNFAIP2", "TNIP1", "SOD2", "ICAM1", "LCN2", "SAA1", "SAA2"), rownames(gene))
writeLines(c(paste("state set:", length(STATE), "genes;", length(S_all), "proteins measured;", length(S_in), "hepatocyte-made; pool", length(setdiff(intersect(rownames(gene), HEPOK), c(S_all, IMM, HEP)))), paste(sort(S_in), collapse = ","),
             paste("immune:", paste(IMM, collapse = ",")), paste("hepmass:", paste(HEP, collapse = ","))), file.path(out, "sets_used.txt"))
sc <- function(gs) colMeans(Z[gs, , drop = FALSE], na.rm = TRUE)
md[, STATE := sc(S_in)][, IMMUNE := sc(IMM)][, HEPMASS := sc(HEP)]
for (p in SINGLE) md[, (p) := Z[p, ]]
base <- c("batch", "age", "bmi_num", "sex")
pcor <- function(z, a, b, vars) {
  z <- z[complete.cases(z[, c(a, b, vars), with = FALSE])]
  xx <- model.matrix(reformulate(vars), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z[[a]]), rank(z[[b]]))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  c(n = nrow(z), rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
# matched null sets
ab <- rowMeans(gene, na.rm = TRUE); mi <- rowMeans(is.na(gene))
strat <- paste(cut(ab, unique(quantile(ab, seq(0, 1, 0.1))), include.lowest = TRUE, labels = FALSE),
               cut(mi, c(-1, 0, 0.2, 1), labels = FALSE))
names(strat) <- rownames(gene)
pool <- setdiff(intersect(rownames(gene), HEPOK), c(S_all, IMM, HEP))
draw <- function() unlist(lapply(split(S_in, strat[S_in]), function(v) { cand <- pool[strat[pool] == strat[v[1]]]; cand[sample.int(length(cand), length(v), replace = length(cand) < length(v))] }))
B <- 2000
M <- md[masld == 1]
HIST <- c("Ballooning", "Steatosis", "Inflammation", "Fibrosis", "NAS")
obs <- sapply(HIST, function(h) pcor(M, "STATE", h, base)["rho"])
obs_grp <- coef(summary(lm(reformulate(c("masld", base), "STATE"), data = md)))["masld", 1]
nul <- t(replicate(B, {
  s <- draw(); v <- colMeans(Z[s, , drop = FALSE], na.rm = TRUE); m2 <- copy(md); m2[, R := v]
  c(sapply(HIST, function(h) pcor(m2[masld == 1], "R", h, base)["rho"]), grp = coef(summary(lm(reformulate(c("masld", base), "R"), data = m2)))["masld", 1])
}))
colnames(nul) <- c(HIST, "MASLD_vs_control")
ob <- c(obs, obs_grp); names(ob) <- colnames(nul)
emp <- function(o, nv) (1 + sum(abs(nv - mean(nv)) >= abs(o - mean(nv)))) / (1 + length(nv))
NULLT <- data.table(stat = colnames(nul), obs = ob, null_mean = colMeans(nul), null_sd = apply(nul, 2, sd),
                    excess = ob - colMeans(nul), p_matched = sapply(seq_along(ob), function(k) emp(ob[k], nul[, k])))
NULLT[stat == "Ballooning", role := "PRIMARY"]
fwrite(NULLT, file.path(out, "state_score_matched_null.tsv"), sep = "\t")
saveRDS(nul, file.path(out, "matched_null.rds"))
# descriptive table
rows <- list()
for (y in c("STATE", "IMMUNE", "HEPMASS", SINGLE)) {
  for (h in HIST) { r <- pcor(M, y, h, base); rows[[length(rows) + 1]] <- data.table(y = y, test = paste0(h, "_within_MASLD"), n = r["n"], est = r["rho"], p = r["p"]) }
  z <- md[complete.cases(md[, c(y, "masld", base), with = FALSE])]; s <- coef(summary(lm(reformulate(c("masld", base), y), data = z)))["masld", ]
  rows[[length(rows) + 1]] <- data.table(y = y, test = "MASLD_vs_control_beta", n = nrow(z), est = s[1], p = s[4])
}
for (adj in list(c("IMMUNE"), c("HEPMASS"), c("IMMUNE", "HEPMASS"))) {
  r <- pcor(M, "STATE", "Ballooning", c(base, adj))
  rows[[length(rows) + 1]] <- data.table(y = "STATE", test = paste0("Ballooning_within_MASLD_adj_", paste(adj, collapse = "+")), n = r["n"], est = r["rho"], p = r["p"])
}
r <- pcor(M, "STATE", "Ballooning", c(base, "Steatosis", "Inflammation", "Fibrosis"))
rows[[length(rows) + 1]] <- data.table(y = "STATE", test = "Ballooning_within_MASLD_adj_other_histology", n = r["n"], est = r["rho"], p = r["p"])
DT <- rbindlist(rows); DT[, bh := p.adjust(p, "BH")]
fwrite(DT, file.path(out, "pxd051911_descriptive.tsv"), sep = "\t")
fwrite(md[, c("sample_id", "masld", HIST, "STATE", "IMMUNE", "HEPMASS", SINGLE), with = FALSE], file.path(out, "pxd051911_participant_scores.tsv"), sep = "\t")
# per-protein ballooning rho for the set
pp <- rbindlist(lapply(S_in, function(g) { m2 <- copy(M); m2[, R := Z[g, M$sample_id]]; r <- pcor(m2, "R", "Ballooning", base)
  data.table(protein = g, n = r["n"], rho = r["rho"], p = r["p"]) }))
fwrite(pp[order(p)], file.path(out, "pxd051911_state_proteins_ballooning.tsv"), sep = "\t")
# replication cohorts (job 27 processing)
to_gene <- function(x, gsym) {
  ok <- !is.na(gsym) & gsym != "" & !grepl(";", gsym, fixed = TRUE)
  x <- x[ok, , drop = FALSE]; gsym <- gsym[ok]; x[!is.finite(x)] <- NA_real_
  x <- normalizeBetweenArrays(x, method = "quantile"); idx <- split(seq_len(nrow(x)), gsym)
  t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x)))) }
pa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_proteome.csv"))
aa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_annotation_file.csv"), select = 1:6)
setnames(aa, c("patient", "sample", "group", "bmi", "age", "sex"))
sa <- intersect(aa$sample, names(pa)); aa <- aa[match(sa, sample)]
XA <- as.matrix(pa[, ..sa]); storage.mode(XA) <- "double"; GA <- to_gene(XA, pa[["Gene name"]]); colnames(GA) <- sa
aa[, grp := fcase(group == "NASH", "disease", group == "control", "control", group == "Cirrhosis", "cirrhosis")]
pb <- fread(file.path(EXT, "PXD065892/liver_20261017T003606Z/liver_sn_minprob_processed.tsv"))
XB <- t(as.matrix(pb[, -(1:3)])); storage.mode(XB) <- "double"
GB <- to_gene(XB, sub("~.*$", "", rownames(XB))); colnames(GB) <- pb$sample
ab2 <- data.table(sample = pb$sample, subject = pb$subject, group = pb$group)
ab2[, grp := fcase(group == "NAFLD", "disease", group == "Healthy", "control", group == "AIH", "AIH")]
rep_rows <- list()
for (co in list(list("PXD027722", GA, aa, c("disease", "cirrhosis")), list("PXD065892", GB, ab2, c("disease", "AIH")))) {
  G <- co[[2]]; ann <- co[[3]]; Zc <- t(scale(t(G))); s_all <- intersect(STATE, rownames(G)); s_in <- intersect(s_all, HEPOK)
  ab_c <- rowMeans(G, na.rm = TRUE); st_c <- cut(ab_c, unique(quantile(ab_c, seq(0, 1, 0.1))), include.lowest = TRUE, labels = FALSE); names(st_c) <- rownames(G)
  pool_c <- setdiff(intersect(rownames(G), HEPOK), s_all)
  for (cmp in co[[4]]) {
    keep <- ann$grp %in% c("control", cmp); lab <- ann$grp[keep] == cmp
    stat_of <- function(v) mean(v[keep][lab]) - mean(v[keep][!lab])
    v <- colMeans(Zc[s_in, , drop = FALSE], na.rm = TRUE); o <- stat_of(v)
    perm <- replicate(5000, { l2 <- sample(lab); mean(v[keep][l2]) - mean(v[keep][!l2]) })
    nm <- replicate(1000, { s <- unlist(lapply(split(s_in, st_c[s_in]), function(q) { cd <- pool_c[st_c[pool_c] == st_c[q[1]]]; cd[sample.int(length(cd), length(q), replace = length(cd) < length(q))] })); stat_of(colMeans(Zc[s, , drop = FALSE], na.rm = TRUE)) })
    rep_rows[[length(rep_rows) + 1]] <- data.table(cohort = co[[1]], contrast = paste(cmp, "vs control"), n_case = sum(lab), n_control = sum(!lab),
      n_proteins = length(s_in), diff_sd = o, p_label_perm = (1 + sum(abs(perm) >= abs(o))) / 5001, matched_mean = mean(nm), excess = o - mean(nm),
      p_matched = emp(o, nm))
    for (p in intersect(SINGLE, rownames(G))) {
      w <- Zc[p, ]; rep_rows[[length(rep_rows) + 1]] <- data.table(cohort = co[[1]], contrast = paste(cmp, "vs control"), n_case = sum(lab), n_control = sum(!lab),
        n_proteins = 1, protein = p, diff_sd = stat_of(w), p_label_perm = wilcox.test(w[keep][lab], w[keep][!lab])$p.value)
    }
  }
}
REP <- rbindlist(rep_rows, fill = TRUE); fwrite(REP, file.path(out, "replication_cohorts.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220)
print(readLines(file.path(out, "sets_used.txt"))); print(NULLT, digits = 3)
print(DT, digits = 3); print(head(pp[order(p)], 25), digits = 3); print(REP, digits = 3)
