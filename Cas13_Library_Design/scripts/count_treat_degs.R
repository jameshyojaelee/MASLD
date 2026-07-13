#!/usr/bin/env Rscript

# Count canonical integrated human DEGs using analytical limma TREAT p-values
# reconstructed from canonical_deg_results.csv. This avoids refitting the full
# limma-voom model when we only need DEG counts across an LFC cutoff sweep.

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RDIR <- file.path(
  ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration"
)

cuts <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0)
d <- fread(file.path(RDIR, "canonical_deg_results.csv"))

probe <- d[is.finite(t) & is.finite(P.Value) & P.Value > 0 & P.Value < 1 & abs(t) > 1e-6]
probe_idx <- unique(round(seq(1, nrow(probe), length.out = min(nrow(probe), 12))))
df_total <- median(vapply(
  probe_idx,
  function(i) {
    uniroot(function(df) 2 * pt(-abs(probe$t[i]), df = df) - probe$P.Value[i],
            c(0.1, 1e6))$root
  },
  numeric(1)
))

counts <- rbindlist(lapply(cuts, function(cut) {
  p_treat <- pt((abs(d$logFC) - cut) / d$SE, df = df_total, lower.tail = FALSE) +
    pt((abs(d$logFC) + cut) / d$SE, df = df_total, lower.tail = FALSE)
  fdr_treat <- p.adjust(p_treat, method = "BH")
  sig <- fdr_treat < 0.05
  data.table(
    treat_lfc = cut,
    fdr = 0.05,
    total = sum(sig),
    up = sum(sig & d$logFC > 0),
    down = sum(sig & d$logFC < 0),
    min_abs_logFC = if (any(sig)) min(abs(d$logFC[sig])) else NA_real_,
    inferred_df_total = df_total
  )
}))

out <- file.path(ROOT, "Cas13_Library_Design/data/treat_deg_counts.csv")
fwrite(counts, out)
print(counts)
cat("Wrote", out, "\n")
