# DIG1 job 66 (post hoc, after the job 63 results were read). Checks and extensions of the genotype -> state result.
# Read before this script: job 63 output (PNPLA3 148M +0.18 SD/allele, flat in GSE135251, removed by NAS+fibrosis
# adjustment, not specific: FIBRO and DUCTULAR rise as much; HSD17B13 junction allele -0.14, Holm 0.013, carried by
# GSE213621 and GSE126848, specific to the hepatocyte genes: FIBRO/IMMUNE/DUCTULAR flat).
# A. Per-cohort call-validity control: NAS, fibrosis, disease group ~ allele + sex, each cohort, both variants.
# B. Cis check for the junction call: HSD17B13 log2 CPM (within-cohort z) ~ junction allele + dataset + sex.
#    Literature: HSD17B13 expression falls with allele dosage (Pirola 2018, doi 10.1194/jlr.P089953).
# C. GSE213621 (where the HSD17B13 effect sits): INFLAMED_HEP ~ allele + sex [+ fibrosis] [+ group].
# D. Reverse-direction check: in junction heterozygotes, continuous skip ratio vs INFLAMED_HEP (+ dataset + sex).
#    If the inflamed state changed exon 6 splicing, the ratio would track the score within a genotype class.
# E. Transcriptome-wide genotype scan (bulk counts; limma-voom with quality weights is native to this modality):
#    expression ~ allele + dataset + sex, genes kept by filterByExpr(group = dataset), robust eBayes, BH over genes.
#    Then cameraPR on the allele t (inter-gene correlation 0.01; rank version as sensitivity):
#    family H = 50 Hallmarks (BH); family T (named here, Holm) = single-cell state sets from jobs 53/56, CCL20 anchor,
#    depth mode: STATE_UP = mean excess > +0.25 in all four cohorts, STATE_DN = < -0.25 in all four; same for the
#    CXCL10 anchor (comparator). Score genes (FCAMR, TRIM31, UBD, CXCL10, IL32, CCL20), HSD17B13 and PNPLA3 removed.
# Unit: one bulk RNA-seq sample per participant (774). No ancestry covariates. Cross-sectional, disease-selected.
# Usage: Rscript 66_genotype_posthoc.R <job63_out_dir> <new_out_dir>
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table)})
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 2)
g63 <- args[1]; out <- args[2]; stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
rel <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
GMT <- file.path(ROOT, "Analysis/downstream_analysis/pathway_analysis/data/msigdb.v2025.1.Hs.symbols.gmt")
D <- read.delim(file.path(g63, "sample_scores_genotypes.tsv"))
d <- readRDS(file.path(rel, "merged_dge.rds"))
stopifnot(identical(colnames(d$counts), as.character(D$sample_id)))
deg <- read.csv(file.path(rel, "deg_results.csv")); sym <- setNames(deg$symbol, deg$gene)
rn <- rownames(d$counts); g_sym <- ifelse(rn %in% names(sym) & !is.na(sym[rn]) & sym[rn] != "", sym[rn], rn)
D$disease <- as.numeric(D$group_binary == "Disease")
fit <- function(f, dat) {
  v <- all.vars(f); cc <- dat[complete.cases(dat[, v, drop = FALSE]), , drop = FALSE]
  rhs <- attr(terms(f), "term.labels"); rhs <- rhs[!(rhs %in% c("dataset", "sex_final") & sapply(rhs, function(t) t %in% names(cc) && length(unique(cc[[t]])) < 2))]
  f <- reformulate(rhs, response = v[1])
  mod <- lm(f, data = cc); X <- model.matrix(mod); e <- residuals(mod); h <- hatvalues(mod)
  XtXi <- solve(crossprod(X)); V <- XtXi %*% crossprod(X * (e / (1 - h))) %*% XtXi
  b <- coef(mod); se <- sqrt(diag(V))
  data.frame(term = names(b), estimate = b, se = se, p = 2 * pt(-abs(b / se), df.residual(mod)), n = nobs(mod), row.names = NULL)
}
VAR <- c("PNPLA3_148M", "HSD17B13_junction")
rows <- list(); add <- function(r, ...) { a <- list(...); for (k in names(a)) r[[k]] <- a[[k]]; rows[[length(rows) + 1]] <<- r }
# A
for (v in VAR) for (y in c("nas_score", "fibrosis_stage", "disease")) for (co in unique(D$dataset)) {
  dd <- D[D$dataset == co & !is.na(D[[y]]) & !is.na(D[[v]]), ]
  if (nrow(dd) < 20 || length(unique(dd[[y]])) < 2 || length(unique(dd[[v]])) < 2) next
  r <- fit(as.formula(paste(y, "~", v, "+ sex_final")), dd); add(r[r$term == v, ], part = "A_control_per_cohort", outcome = y, variant = v, cohort = co)
}
for (y in c("nas_score", "fibrosis_stage", "disease")) {
  dd <- D[!is.na(D[[y]]) & !is.na(D$HSD17B13_junction), ]
  r <- fit(as.formula(paste(y, "~ HSD17B13_junction + dataset + sex_final")), dd); add(r[r$term == "HSD17B13_junction", ], part = "A_control_pooled", outcome = y, variant = "HSD17B13_junction", cohort = "pooled")
}
# B
d <- calcNormFactors(d); lc <- cpm(d, log = TRUE, prior.count = 2)
zc <- function(x) { v <- rep(NA_real_, length(x)); for (co in unique(D$dataset)) { j <- D$dataset == co; v[j] <- as.numeric(scale(x[j])) }; v }
for (g in c("HSD17B13", "PNPLA3")) {
  D[[paste0("cis_", g)]] <- zc(lc[which(g_sym == g)[1], ])
  for (v in VAR) { r <- fit(as.formula(paste0("cis_", g, " ~ ", v, " + dataset + sex_final")), D); add(r[r$term == v, ], part = "B_cis", outcome = paste0("cis_", g), variant = v, cohort = "pooled") }
}
# C
G2 <- D[D$dataset == "GSE213621", ]
for (f in c("INFLAMED_HEP ~ HSD17B13_junction + sex_final", "INFLAMED_HEP ~ HSD17B13_junction + sex_final + fibrosis_stage",
            "INFLAMED_HEP ~ HSD17B13_junction + sex_final + group_binary", "INFLAMED_HEP ~ HSD17B13_junction + sex_final + fibrosis_stage + group_binary",
            "fibrosis_stage ~ HSD17B13_junction + sex_final")) {
  r <- fit(as.formula(f), G2); add(r[r$term == "HSD17B13_junction", ], part = "C_GSE213621", outcome = f, variant = "HSD17B13_junction", cohort = "GSE213621")
}
tab213 <- with(G2, table(group_binary, fibrosis_stage, useNA = "ifany"))
# D
J <- read.delim(file.path(RUN, "out_hsd17b13_junctions/hsd17b13_skip6.tsv"))
D$skip_ratio <- J$skip_ratio[match(D$sample_id, J$sample_id)]
het <- D[D$HSD17B13_junction %in% 1, ]
r <- fit(INFLAMED_HEP ~ skip_ratio + dataset + sex_final, het); add(r[r$term == "skip_ratio", ], part = "D_within_het", outcome = "INFLAMED_HEP", variant = "skip_ratio", cohort = "pooled")
r <- fit(nas_score ~ skip_ratio + dataset + sex_final, het); add(r[r$term == "skip_ratio", ], part = "D_within_het", outcome = "nas_score", variant = "skip_ratio", cohort = "pooled")
R <- do.call(rbind, rows)
fwrite(R, file.path(out, "posthoc_models.tsv"), sep = "\t")
# E
gl <- strsplit(readLines(GMT), "\t"); all_sets <- setNames(lapply(gl, function(x) unique(x[-(1:2)])), sapply(gl, `[`, 1))
HALL <- grep("^HALLMARK_", names(all_sets), value = TRUE); stopifnot(length(HALL) == 50)
g53 <- fread(file.path(RUN, "out_chemokine_hepatocyte_state/gene_lfc_matched.tsv.gz"))[donors == "masld"][, donors := NULL]
g56 <- fread(file.path(RUN, "out_atlas_state_replication/gene_lfc_matched.tsv.gz"))
SC <- rbind(g53, g56, use.names = TRUE)[mode == "depth" & is.finite(mean_excess)]
DROP <- c("FCAMR", "TRIM31", "UBD", "CXCL10", "IL32", "CCL20", "HSD17B13", "PNPLA3")
state_set <- function(a, dir) {
  x <- SC[anchor == a, .(k = .N, mn = min(mean_excess), mx = max(mean_excess)), by = gene][k == 4]
  setdiff(if (dir > 0) x[mn > 0.25, gene] else x[mx < -0.25, gene], DROP)
}
TSETS <- list(CCL20_STATE_UP = state_set("CCL20", 1), CCL20_STATE_DN = state_set("CCL20", -1),
              CXCL10_STATE_UP = state_set("CXCL10", 1), CXCL10_STATE_DN = state_set("CXCL10", -1))
E <- list(); CAM <- list()
for (v in VAR) {
  keep_s <- !is.na(D[[v]]); dv <- d[, keep_s]; Dv <- D[keep_s, ]
  keep_g <- filterByExpr(dv, group = Dv$dataset); dv <- dv[keep_g, , keep.lib.sizes = FALSE]; dv <- calcNormFactors(dv)
  X <- model.matrix(as.formula(paste("~", v, "+ dataset + sex_final")), data = Dv)
  vw <- voomWithQualityWeights(dv, X)
  fm <- eBayes(lmFit(vw, X), robust = TRUE)
  tt <- topTable(fm, coef = v, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt); tt$symbol <- g_sym[match(tt$gene, rn)]; tt$variant <- v
  E[[v]] <- tt
  st <- tapply(tt$t, tt$symbol, function(z) z[which.max(abs(z))])
  st <- st[!(names(st) %in% DROP)]
  for (fam in c("H", "T")) {
    sets <- if (fam == "H") all_sets[HALL] else TSETS
    idx <- ids2indices(sets, names(st)); idx <- idx[lengths(idx) >= 10]
    a <- cameraPR(st, idx, inter.gene.cor = 0.01, sort = FALSE); rk <- cameraPR(st, idx, use.ranks = TRUE, sort = FALSE)
    x <- data.table(variant = v, family = fam, set = names(idx), n_genes = lengths(idx), mean_t = sapply(idx, function(i) mean(st[i])),
                    direction = a$Direction, p = a$PValue, dir_rank = rk$Direction, p_rank = rk$PValue)
    x[, p_adj := p.adjust(p, if (fam == "H") "BH" else "holm")]
    CAM[[paste(v, fam)]] <- x
  }
}
EE <- rbindlist(lapply(E, as.data.table)); fwrite(EE, file.path(out, "genotype_transcriptome_scan.tsv.gz"), sep = "\t")
CC <- rbindlist(CAM); fwrite(CC, file.path(out, "genotype_scan_camera.tsv"), sep = "\t")
writeLines(sprintf("%s\t%d", names(TSETS), lengths(TSETS)), file.path(out, "state_set_sizes.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 220)
print(R, digits = 3); print(tab213)
print(readLines(file.path(out, "state_set_sizes.tsv")))
print(CC[family == "T"], digits = 3)
print(CC[family == "H" & p_adj < 0.05][order(variant, p)], digits = 3)
for (v in VAR) { x <- EE[variant == v][order(P.Value)][1:30, .(symbol, logFC = round(logFC, 3), t = round(t, 2), P.Value = signif(P.Value, 2), adj.P.Val = signif(adj.P.Val, 2))]; cat(v, "n BH<0.05:", sum(EE[variant == v]$adj.P.Val < 0.05), "\n"); print(x) }
