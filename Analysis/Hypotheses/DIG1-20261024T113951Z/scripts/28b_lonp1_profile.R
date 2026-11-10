# DIG1 job 28b: job 28 code plus per-group per-protein tables for every protease, and a profile check:
# across mitochondrial proteins, does the protease-linked profile from the three cohorts other than
# PXD051911 match the PXD051911 within-MASLD ballooning profile (job 22b)? Descriptive only;
# proteins are not independent, so the other proteases serve as the comparison.
# Original job 28 header follows.
# DIG1 job 28 (exploratory; follows jobs 26 and 27). Within diseased liver, LONP1 protein is higher
# where Complex I is lower in all four disease groups tested so far. In mouse heart, LonP1
# overexpression lowers Complex I content (Venkatesh 2019, doi:10.1016/j.yjmcc.2018.12.017). Is
# LONP1 abundance tied to less Complex I *relative to mitochondrial content*, is it specific to
# LONP1 among mitochondrial proteases, and does LONP1 also go with less of its known substrates?
# Cohorts and groups (donor = unit):
#   PXD051911 MASLD 46, control 12 (matrix and covariates as job 22: batch, age, BMI, sex)
#   PXD027722 NASH + cirrhosis 30, control 15; PXD065892 NAFLD 15, AIH 19, healthy 9 (as job 27)
# Mitochondrial proteins: HPA main or additional location contains "Mitochondria", plus all
# OXPHOS and mitoribosome subunits. Mito content = mean z of mitochondrial proteins not in any
# tested set and not a protease. Scores: mean z over members.
# Statistic: Spearman of LONP1 with a set score, both ranks residualised on mito content (+ the
# cohort covariates). Combined over the four disease groups by Fisher z, weights n - 3 - k.
# Fixed substrate set (literature, checked only for DHODH here; others [Unverified] recall):
#   DHODH (Xu 2025, doi:10.1016/j.jhep.2025.07.013), TFAM, ACO2, ALAS1, COX4I1.
# Primary family (Holm over 2), disease groups combined:
#   L1  LONP1 vs CI (given mito content) < 0.
#   L2  LONP1 vs substrate set (given mito content) < 0.
# Descriptive (BH over the table): every set x protease (LONP1, CLPP, YME1L1, AFG3L2, SPG7, and
# chaperone HSPD1 as a non-protease comparator) in disease and in control groups; per mitochondrial
# protein, combined disease-group rho with LONP1 given mito content.
# Usage: Rscript 28_lonp1_axis.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; EXT <- file.path(ROOT, "data/external")
to_gene <- function(x, gsym, maxmiss = 1) {
  ok <- !is.na(gsym) & gsym != "" & !grepl(";", gsym, fixed = TRUE)
  x <- x[ok, , drop = FALSE]; gsym <- gsym[ok]; x[!is.finite(x)] <- NA_real_
  keep <- rowMeans(is.na(x)) <= maxmiss; x <- x[keep, , drop = FALSE]; gsym <- gsym[keep]
  x <- normalizeBetweenArrays(x, method = "quantile")
  idx <- split(seq_len(nrow(x)), gsym)
  t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x)))) }
# PXD051911
raw <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt")); md <- fread(file.path(ROOT, "data/PXD051911/meta_data.txt"))
smp <- setdiff(names(raw), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% smp], by = "liver_proteomics_filename"); setorder(md, liver_proteomics_filename)
x1 <- as.matrix(raw[, md$liver_proteomics_filename, with = FALSE]); storage.mode(x1) <- "double"; x1[x1 <= 0] <- NA; x1 <- log2(x1)
G1 <- to_gene(x1, raw$Genes, 0.5); colnames(G1) <- md$liver_proteomics_filename
a1 <- md[, .(sample = liver_proteomics_filename, grp = fifelse(saf_diagnosis != "No_MASLD", "disease", "control"),
             batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")), age = as.numeric(alder),
             bmi = as.numeric(bmi), sex = factor(gender))]
# PXD027722
pa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_proteome.csv"))
aa <- fread(file.path(EXT, "PXD027722/download_20261017T012300Z/PatientSamples_annotation_file.csv"), select = 1:3)
setnames(aa, c("patient", "sample", "group")); aa <- aa[sample %in% names(pa)]
x2 <- as.matrix(pa[, aa$sample, with = FALSE]); storage.mode(x2) <- "double"
G2 <- to_gene(x2, pa[["Gene name"]]); colnames(G2) <- aa$sample
a2 <- aa[, .(sample, grp = fifelse(group == "control", "control", "disease"))]
# PXD065892
pb <- fread(file.path(EXT, "PXD065892/liver_20261017T003606Z/liver_sn_minprob_processed.tsv"))
x3 <- t(as.matrix(pb[, -(1:3)])); storage.mode(x3) <- "double"
G3 <- to_gene(x3, sub("~.*$", "", rownames(x3))); colnames(G3) <- pb$sample
a3 <- data.table(sample = pb$sample, grp = fcase(pb$group == "Healthy", "control", pb$group == "NAFLD", "disease", pb$group == "AIH", "AIH"))
GROUPS <- list(
  list(id = "PXD051911_MASLD", G = G1, a = a1, g = "disease", cov = c("batch", "age", "bmi", "sex")),
  list(id = "PXD027722_NASH_cirrhosis", G = G2, a = a2, g = "disease", cov = character()),
  list(id = "PXD065892_NAFLD", G = G3, a = a3, g = "disease", cov = character()),
  list(id = "PXD065892_AIH", G = G3, a = a3, g = "AIH", cov = character()),
  list(id = "PXD051911_control", G = G1, a = a1, g = "control", cov = c("batch", "age", "bmi", "sex")),
  list(id = "PXD027722_control", G = G2, a = a2, g = "control", cov = character()),
  list(id = "PXD065892_healthy", G = G3, a = a3, g = "control", cov = character()))
hpa <- fread(file.path(EXT, "hpa/proteinatlas.tsv"), select = c("Gene", "Subcellular main location", "Subcellular additional location"))
hmito <- hpa[grepl("Mitochondria", `Subcellular main location`) | grepl("Mitochondria", `Subcellular additional location`), unique(Gene)]
PROT <- c("LONP1", "CLPP", "YME1L1", "AFG3L2", "SPG7", "HSPD1")
SUBS <- c("DHODH", "TFAM", "ACO2", "ALAS1", "COX4I1")
sets_for <- function(g) {
  S <- list(CI = setdiff(g[grepl("^NDUF[ABSCV][0-9]", g)], "NDUFA4"), CII = intersect(c("SDHA", "SDHB", "SDHC", "SDHD"), g),
            CIII = g[grepl("^UQCR|^CYC1$", g)], CIV = setdiff(union(g[grepl("^COX[4-8]", g)], intersect("NDUFA4", g)), SUBS),
            CV = g[grepl("^ATP5", g)], mt_ribosome = g[grepl("^MRP[LS][0-9]", g)], substrates = intersect(SUBS, g),
            hep_mito = intersect(c("CPS1", "OTC", "HMGCS2", "GLS2"), g),
            FAO = intersect(c("CPT2", "ACADM", "ACADL", "ACADVL", "ACADS", "HADHA", "HADHB", "ECHS1", "HADH", "ETFA", "ETFB", "ETFDH"), g))
  S[lengths(S) > 0] }
pres <- function(v, Z) { R <- cbind(rank(v)); qr.resid(qr(Z), R)[, 1] }
res <- list(); per <- list(); info <- list()
for (k in GROUPS) {
  an <- k$a[grp == k$g]; if (length(k$cov)) an <- an[complete.cases(an[, k$cov, with = FALSE])]
  s <- an$sample; G <- k$G[, s, drop = FALSE]
  G <- G[rowSums(is.na(G)) == 0, , drop = FALSE]  # complete proteins only within the group, so ranks are valid
  Zs <- t(scale(t(G))); Zs <- Zs[rowSums(!is.finite(Zs)) == 0, , drop = FALSE]
  g <- rownames(Zs); S <- sets_for(g); G <- G[g, , drop = FALSE]
  mito <- union(intersect(hmito, g), unlist(S[c("CI", "CII", "CIII", "CIV", "CV", "mt_ribosome")]))
  bg <- setdiff(mito, c(unlist(S), PROT))
  content <- colMeans(Zs[bg, , drop = FALSE], na.rm = TRUE)
  Xc <- model.matrix(reformulate(c("content", k$cov)), data = data.table(an, content = content))
  stopifnot(nrow(Xc) == length(s))
  df <- length(s) - ncol(Xc) - 1
  info[[k$id]] <- data.table(group = k$id, n = length(s), n_mito = length(mito), n_background = length(bg), df = df)
  for (p in intersect(PROT, g)) {
    rp <- pres(G[p, ], Xc)
    for (sn in names(S)) { sc <- colMeans(Zs[S[[sn]], , drop = FALSE], na.rm = TRUE); r <- cor(rp, pres(sc, Xc))
      res[[length(res) + 1]] <- data.table(group = k$id, protease = p, set = sn, n = length(s), k = ncol(Xc), rho = r) }
    mm <- setdiff(mito, p); RM <- qr.resid(qr(Xc), apply(G[mm, , drop = FALSE], 1, rank))
    per[[paste(k$id, p)]] <- data.table(group = k$id, protease = p, gene = mm, rho = as.vector(cor(rp, RM)), n = length(s), k = ncol(Xc))
  }
}
res <- rbindlist(res); per <- rbindlist(per); info <- rbindlist(info)
fwrite(info, file.path(out, "design.tsv"), sep = "\t"); fwrite(res, file.path(out, "per_group_rho.tsv"), sep = "\t")
comb <- function(d) { w <- d$n - 3 - (d$k - 1); z <- atanh(d$rho); zc <- sum(w * z) / sum(w); se <- 1 / sqrt(sum(w))
  data.table(n_groups = nrow(d), n_total = sum(d$n), rho_comb = tanh(zc), z = zc / se, p = 2 * pnorm(-abs(zc / se)), n_neg = sum(d$rho < 0)) }
res[, kind := fifelse(grepl("control|healthy", group), "control", "disease")]
cm <- res[, comb(.SD), by = .(kind, protease, set)]
prim <- cm[kind == "disease" & protease == "LONP1" & set %in% c("CI", "substrates")]
prim[, p_one_sided := pnorm(z)][, p_holm := p.adjust(p_one_sided, "holm")]
fwrite(prim, file.path(out, "primary_holm.tsv"), sep = "\t")
cm[, bh := p.adjust(p, "BH")]
fwrite(cm[order(kind, protease, p)], file.path(out, "combined_sets.tsv"), sep = "\t")
fwrite(per, file.path(out, "per_group_per_protein.tsv.gz"), sep = "\t")
pc <- per[protease == "LONP1" & !grepl("control|healthy", group)][, comb(.SD), by = gene]
pc[, bh := p.adjust(p, "BH")]
lab <- function(x) fcase(grepl("^NDUF[ABSCV][0-9]", x) & x != "NDUFA4", "CI", x %in% c("SDHA", "SDHB", "SDHC", "SDHD"), "CII",
                         grepl("^UQCR|^CYC1$", x), "CIII", grepl("^COX[4-8]", x) | x == "NDUFA4", "CIV", grepl("^ATP5", x), "CV",
                         grepl("^MRP[LS][0-9]", x), "mt_ribosome", x %in% SUBS, "substrate", default = "other_mito")
pc[, class := lab(gene)]
fwrite(pc[order(p)], file.path(out, "per_mito_protein_lonp1_disease_combined.tsv"), sep = "\t")
fwrite(pc[n_groups == 4, .(n = .N, mean_rho = mean(rho_comb), frac_neg = mean(rho_comb < 0)), by = class][order(mean_rho)], file.path(out, "per_protein_class_summary.tsv"), sep = "\t")
bal <- fread(file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z/out_mtdna_link/per_protein_effects.tsv"))[, .(gene, rho_ballooning)]
prof <- rbindlist(lapply(PROT, function(pr) { d <- per[protease == pr & !grepl("control|healthy|PXD051911", group)][, comb(.SD), by = gene]
  d <- merge(d[n_groups == 3], bal, by = "gene"); d[, class := lab(gene)]
  ct <- cor.test(d$rho_comb, d$rho_ballooning, method = "spearman", exact = FALSE)
  rbind(data.table(protease = pr, proteins = "all_mito", n = nrow(d), spearman = unname(ct$estimate), p = ct$p.value),
        d[class == "other_mito", { c2 <- cor.test(rho_comb, rho_ballooning, method = "spearman", exact = FALSE)
          data.table(protease = pr, proteins = "other_mito_only", n = .N, spearman = unname(c2$estimate), p = c2$p.value) }]) }))
fwrite(prof, file.path(out, "profile_vs_PXD051911_ballooning.tsv"), sep = "\t")
print(prof, digits = 3)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(info); print(prim, digits = 3); print(cm[protease == "LONP1"][order(kind, p)], digits = 2)
print(cm[set == "CI"][order(kind, p)], digits = 2); print(fread(file.path(out, "per_protein_class_summary.tsv")), digits = 2); print(pc[order(p)][1:30], digits = 2)
