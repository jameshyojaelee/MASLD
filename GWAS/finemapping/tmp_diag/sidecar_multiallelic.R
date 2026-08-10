# Is the panel AF sidecar wrong because MULTIPLE records share a position and
# the wrong one is being used? 00a asserts nrow(afreq)==nrow(bim), so the
# sidecar has one row PER BIM LINE -- and a multi-allelic site occupies several
# BIM lines at the SAME position. Any consumer keying on position alone then
# picks arbitrarily.
suppressPackageStartupMessages(library(data.table))
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
h <- unique(fread(file.path(FM,"tmp_diag/reoriented_highpip.tsv")),
            by = c("chromosome","position"))
cat(sprintf("checking %d distinct positions from the re-oriented set\n\n", nrow(h)))
res <- rbindlist(lapply(unique(h$chromosome), function(cc) {
  p <- fread(cmd = sprintf("zcat %s/data/ld_ref/panel_af/eur/chr%d.af.tsv.gz", FM, cc),
             select = c("position","bim_a1","bim_a2","af_a1"), showProgress = FALSE)
  # how many sidecar rows exist at each of our positions?
  p[position %in% h[chromosome==cc]$position][, chromosome := cc][]
}), fill = TRUE)
n_per <- res[, .(n_records = .N,
                 alleles = paste(sprintf("%s/%s=%.4f", bim_a1, bim_a2, af_a1), collapse = "  ")),
             by = .(chromosome, position)]
setorder(n_per, -n_records)
print(n_per, nrows = 30)
cat(sprintf("\npositions with >1 sidecar record: %d of %d\n",
            sum(n_per$n_records > 1), nrow(n_per)))
cat(sprintf("total sidecar rows at these positions: %d\n", nrow(res)))
cat("\n=== GLOBAL: how common are multi-record positions in the EUR sidecar? ===\n")
for (cc in c(19, 11, 1)) {
  p <- fread(cmd = sprintf("zcat %s/data/ld_ref/panel_af/eur/chr%d.af.tsv.gz", FM, cc),
             select = c("position","bim_a1","bim_a2"), showProgress = FALSE)
  d <- p[, .N, by = position]
  pal <- p[(bim_a1 %in% c("A","T") & bim_a2 %in% c("A","T")) |
           (bim_a1 %in% c("C","G") & bim_a2 %in% c("C","G"))]
  dp <- pal[, .N, by = position]
  cat(sprintf("  chr%-3d %8d rows | %8d distinct pos | multi-record pos %6d (%.3f%%) | palindromic rows %7d, multi-record %5d (%.3f%%)\n",
      cc, nrow(p), nrow(d), sum(d$N>1), 100*mean(d$N>1),
      nrow(pal), sum(dp$N>1), 100*mean(dp$N>1)))
}
