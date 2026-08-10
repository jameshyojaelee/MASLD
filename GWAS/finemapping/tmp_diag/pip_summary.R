suppressPackageStartupMessages(library(data.table))
D <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/pip_paired"
d <- rbindlist(lapply(list.files(D, pattern = "[.]tsv$", full.names = TRUE), fread), fill = TRUE)
cat(sprintf("=== TASK #15: do the FITS differ? %d loci, %d studies, anc %s ===\n\n",
    nrow(d), uniqueN(d$study), paste(sort(unique(d$ancestry)), collapse = "/")))
cat("--- the shipped arm really was indefinite (else this tests nothing) ---\n")
print(d[, .(n = .N, med_min_eig_shipped = round(median(min_eig_shipped), 5),
            med_min_eig_gram = round(median(min_eig_gram), 7),
            worst_shipped = round(min(min_eig_shipped), 4)), by = ancestry])
cat("\n--- PIP agreement ---\n")
print(d[, .(n = .N, med_pip_cor = round(median(pip_cor), 6), min_pip_cor = round(min(pip_cor), 6),
            med_maxdiff = round(median(pip_max_abs_diff), 5),
            worst_maxdiff = round(max(pip_max_abs_diff), 4),
            lead_same = sum(lead_same)), by = ancestry])
cat("\n--- credible sets ---\n")
cat(sprintf("  loci where n_CS differs      : %d / %d\n", sum(d$n_cs_shipped != d$n_cs_gram), nrow(d)))
cat(sprintf("  loci where LEAD variant moves: %d / %d\n", sum(!d$lead_same), nrow(d)))
sub <- d[!is.na(cs_jaccard) & (n_cs_shipped > 0 | n_cs_gram > 0)]
cat(sprintf("  loci with any credible set   : %d ; median CS Jaccard %.4f\n",
            nrow(sub), if (nrow(sub)) median(sub$cs_jaccard, na.rm = TRUE) else NA))
cat(sprintf("  n PIP>0.5 shipped %d vs gram %d | PIP>0.1 shipped %d vs gram %d\n",
    sum(d$n_pip_gt05_shipped), sum(d$n_pip_gt05_gram),
    sum(d$n_pip_gt01_shipped), sum(d$n_pip_gt01_gram)))
cat("\n--- worst 6 loci by max |PIP difference| ---\n")
print(d[order(-pip_max_abs_diff)][1:6, .(study, locus, p,
      min_eig = round(min_eig_shipped, 4), pip_cor = round(pip_cor, 6),
      maxdiff = round(pip_max_abs_diff, 4), lead_same,
      cs = paste0(n_cs_shipped, "/", n_cs_gram))])
