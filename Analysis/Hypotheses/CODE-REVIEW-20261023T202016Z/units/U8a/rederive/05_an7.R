# U8a: AN7 sensitivities under the H59 model with the independent GLS and permutation code.
# S1: donor means over QC-passing cells (AN7 donor table), all 117 programs.
# S2: T-cell donors with < 20 T cells removed (cell counts from 04_donor_scores.py), 11 T-cell programs;
#     the other 106 programs keep their primary permutation p (from 02_null.R, seed 20260923).
# Usage: Rscript --vanilla 05_an7.R <seed> <B> <outdir> <ncores>
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1]); B <- as.integer(args[2]); out <- args[3]; ncores <- as.integer(args[4])
here <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U8a/rederive"
source(file.path(here, "u8a_lib.R"))
suppressPackageStartupMessages(library(parallel))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
AN7 <- file.path(ROOT, "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN7/20261021T124200Z/out")

null_family <- function(fit, progs, K) {
  one <- function(p) {
    d <- fit[program == p]
    g <- gls_fit(d$score, d$stage, d$dataset)
    tn <- vapply(seq_len(nrow(K)), function(k)
      tryCatch(gls_fit(d$score, permute_stage(d$donor, d$stage, d$dataset, K[k, ]), d$dataset)$t,
               error = function(e) NA_real_), 0)
    data.table(program = p, n_fit = nrow(d), beta = g$beta, t = g$t, perm_p = perm_p(g$t, tn), n_failed = sum(!is.finite(tn)))
  }
  rbindlist(mclapply(progs, one, mc.cores = ncores, mc.preschedule = FALSE))
}

prim <- load_donor_scores()
donors62 <- sort(unique(prim$donor), method = "radix")
K <- make_keys(seed, B, donors62)
reg <- fread(REGISTRY_FILE); reg[, program := paste0(cell_type, "_", module)]
h59 <- fread(file.path(here, sprintf("out/null_seed_%d/perm_results_117.tsv", seed)))

# S1
qc <- load_donor_scores(file.path(AN7, "donor_scores_v2qc_cells.tsv"))
qcfit <- qc[dataset %in% FIT_DS]
progs <- sort(unique(qcfit$program))
s1 <- null_family(qcfit, progs, K)
s1 <- merge(s1, reg[, .(program, module_name, stability_median)], by = "program")
s1[, perm_q := p.adjust(perm_p, "BH", n = 117)]
s1[, selected := perm_q < 0.05 & stability_median >= 0.5]
cmp <- merge(s1, h59[, .(program, beta_h59 = beta)], by = "program")
say("S1: programs %d; selected %d; smallest q %.4f (%d programs at it); Spearman beta vs H59 %.4f; sign changes %d (max |beta| among them %.4f)",
    nrow(s1), sum(s1$selected), min(s1$perm_q), sum(abs(s1$perm_q - min(s1$perm_q)) < 1e-12),
    cor(cmp$beta, cmp$beta_h59, method = "spearman"), sum(sign(cmp$beta) != sign(cmp$beta_h59)),
    max(pmax(abs(cmp$beta), abs(cmp$beta_h59))[sign(cmp$beta) != sign(cmp$beta_h59)]))
print(s1[program %in% c("hepatocytes_20", "hepatocytes_8")])
say("S1 T-cell fitted donors: %d", qcfit[cell_type == "tcells", uniqueN(donor)])
fwrite(s1, file.path(out, "S1_qc_cells_117.tsv"), sep = "\t")
an7s1 <- fread(file.path(AN7, "S1_v2qc_results_117.tsv"))
say("AN7 S1 columns: %s", paste(names(an7s1), collapse = ","))

# S2
cells <- fread(file.path(here, "out/donor_scores/cells_per_donor_lineage.tsv"))
fit <- prim[dataset %in% FIT_DS & dataset != "GSE189600"]
tc <- merge(unique(fit[cell_type == "tcells", .(donor, dataset)]), cells[cell_type == "tcells", .(donor, n_cells)], by = "donor", all.x = TRUE)
say("S2: fitted T-cell donors %d; missing cell counts %d; < 20 T cells %d (%s)", nrow(tc), sum(is.na(tc$n_cells)),
    sum(tc$n_cells < 20), paste(tc[n_cells < 20, .N, by = dataset][, paste(dataset, N)], collapse = "; "))
drop <- tc[n_cells < 20, donor]
tfit <- fit[cell_type == "tcells" & !donor %in% drop]
s2 <- null_family(tfit, sort(unique(tfit$program)), K)
pall <- rbind(h59[!startsWith(program, "tcells_"), .(program, perm_p)], s2[, .(program, perm_p)])
stopifnot(nrow(pall) == 117L)
pall[, q := p.adjust(perm_p, "BH", n = 117)]
pall <- merge(pall, reg[, .(program, stability_median)], by = "program")
say("S2: T-cell n_fit %s; T-cell q<0.05 %d; selected over 117 %d; smallest q over 117 %.5f (%s); smallest T-cell q %.5f",
    paste(unique(s2$n_fit), collapse = ","), sum(pall[startsWith(program, "tcells_")]$q < 0.05),
    sum(pall$q < 0.05 & pall$stability_median >= 0.5), min(pall$q), pall$program[which.min(pall$q)],
    min(pall[startsWith(program, "tcells_")]$q))
fwrite(s2, file.path(out, "S2_tcell_gate20.tsv"), sep = "\t")
fwrite(pall, file.path(out, "S2_family117.tsv"), sep = "\t")
writeLines(c(sprintf("seed\t%d", seed), sprintf("B\t%d", B), "family\t117 (BH n=117)"), file.path(out, "run_params.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
system2("sha256sum", c(SCORE_FILE, META_FILE, REGISTRY_FILE, file.path(AN7, "donor_scores_v2qc_cells.tsv")),
        stdout = file.path(out, "inputs.sha256"))
