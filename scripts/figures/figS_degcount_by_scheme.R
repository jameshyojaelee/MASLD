#!/usr/bin/env Rscript
# figS_degcount_by_scheme.R
# DEG count vs confound-correction scheme, from the degx_factorial DE benchmark
# (286 method x correction x k cells on the canonical 5-cohort disease_vs_control
# data; DEGs = padj < 0.05 counted per-gene from cells/cell_*.csv). Each point is
# one (method, k) cell; boxes summarize across methods/k per scheme.
#
# Story: correction choice swings the DEG list >2x. dream random-intercept (C8)
# and dataset fixed-effect (C1) are the principled middle; SVA-replace (C4) and
# RUVr (C6r) over-call. Source: degx_factorial. Publication theme; cairo_pdf; ASCII.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE,"scripts/figures/publication_theme.R"))
source(file.path(BASE,"scripts/figures/method_correction_labels.R"))
CSV <- "/gpfs/commons/home/jameslee/degx/runs/robustness/degcount_by_scheme.csv"
OUT <- file.path(BASE,"figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

d <- read.csv(CSV, stringsAsFactors=FALSE)
d <- d[d$scheme != "C10", ]                       # C10 = single degenerate cell (438 DEGs), excluded

# Readable correction labels from the shared lookup (compact "code: scheme" form);
# correction FAMILY groups schemes for the colour axis.
corr_family <- function(s) {
  s <- as.character(s)
  ifelse(s %in% c("C0"), "None",
  ifelse(s %in% c("C1","C2","C3"), "Fixed-effect",
  ifelse(s %in% c("C4","C5","C6hk"), "SVA",
  ifelse(s %in% c("C6g","C6s","C6r"), "RUV",
  ifelse(s %in% c("C8","C9","C10"), "dream", "Other")))))
}
famcol <- c(None="#9E9E9E", `Fixed-effect`="#1565C0", SVA="#C9265E", RUV="#E8A33D",
            dream="#2E7D32", Other="#7E57C2")

scheme_lab <- setNames(correction_label(unique(d$scheme), short=TRUE), unique(d$scheme))
d$label  <- factor(scheme_lab[d$scheme],
                   levels=scheme_lab[names(sort(tapply(d$n_deg_05, d$scheme, median)))])
d$family <- factor(corr_family(d$scheme), levels=names(famcol))

p <- ggplot(d, aes(label, n_deg_05)) +
  geom_boxplot(aes(color=family), fill=NA, outlier.shape=NA, width=0.62, linewidth=0.4) +
  geom_jitter(aes(color=family), width=0.16, height=0, size=0.9, alpha=0.65) +
  scale_color_manual(values=famcol, name="Correction family") +
  scale_y_continuous(breaks=seq(0,20000,4000)) +
  coord_flip() +
  labs(title="DEG count depends strongly on confound-correction scheme",
       subtitle="Disease vs Control (846 samples, 5 cohorts); each point = one (method, k) cell; DEGs at padj < 0.05",
       x=NULL, y="# DEGs (padj < 0.05)") +
  theme_masld() + theme_pub() +
  theme(legend.position="right", panel.grid.major.y=element_blank())

save_fig(p, file.path(OUT,"DEGcount_by_scheme.pdf"), width=7.0, height=4.4)
cat("wrote", file.path(OUT,"DEGcount_by_scheme.pdf"), "\n")

# companion: median DEGs per scheme table
s <- d %>% group_by(scheme,label,family) %>%
  summarise(median_deg=median(n_deg_05), min=min(n_deg_05), max=max(n_deg_05), n_cells=n(), .groups="drop") %>%
  arrange(median_deg)
write.csv(s, file.path(OUT,"DEGcount_by_scheme_summary.csv"), row.names=FALSE)
cat("wrote summary csv (", nrow(s), "schemes)\n")
