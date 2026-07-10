#!/usr/bin/env Rscript
# composite_gwas.R — CANDIDATE composite (⑤ GWAS genetic-architecture)
# LEFT : LDSC cross-trait genetic correlation (rg) matrix among EUR liver-disease
#        GWAS (pie glyphs), grouped by trait category. RIGHT: observed vs liability
#        SNP-heritability lollipop per trait (the total-vs-scale gap; the direct
#        analog of the example's twin-vs-SNP heritability).
# Statistical genetics x functional genomics. Output: FIG2_DIR/panels/composite_gwas.pdf
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
source(file.path(BASE, "scripts/figures/composite_helpers.R"))
PANEL_DIR <- file.path(FIG2_DIR, "panels"); DATA_DIR <- file.path(PANEL_DIR, "data")
RG <- file.path(BASE, "GWAS/ldsc/rg")

# trait metadata: ONE well-powered GWAS per DISTINCT trait (avoids same-trait cross-cohort
# rg>1 caps + halves NA). Redundant cohorts (MVP/FinnGen replicates), PDFF (sign anomaly),
# and low-power FinnGen diseases dropped. K = pop prev (NA=continuous), P = sample prev.
meta <- fread(text = "key,display,group,K,P
ghodsian_nafld,NAFLD,Liver disease,0.25,0.0108
ghouse_cirr,Cirrhosis,Liver disease,0.003,0.02
ukbb_alt,ALT,Liver enzymes,NA,NA
ukbb_ast,AST,Liver enzymes,NA,NA
finngen_obesity,Obesity,Metabolic,0.4,0.089")

# ── parse LDSC h2 (observed) ──
h2_obs <- function(key) {
  f <- file.path(RG, paste0("h2_", key, ".log")); if (!file.exists(f)) return(NA_real_)
  l <- grep("Total Observed scale h2", readLines(f), value = TRUE)
  if (!length(l)) return(NA_real_); as.numeric(sub(".*h2: *([0-9.eE-]+).*", "\\1", l[1]))
}
# analytic observed -> liability conversion (Lee 2011)
lia <- function(h2o, K, P) {
  if (is.na(K) || is.na(P)) return(h2o)              # continuous: no conversion
  t <- qnorm(1 - K); z <- dnorm(t)
  h2o * K^2 * (1 - K)^2 / (P * (1 - P) * z^2)
}
meta[, h2_observed := sapply(key, h2_obs)]
meta[, h2_liability := mapply(lia, h2_observed, K, P)]
meta <- meta[!is.na(h2_observed)]
avail <- meta$key

# ── parse rg from all wp_*.log into a symmetric matrix ──
rg <- matrix(NA_real_, length(avail), length(avail), dimnames = list(avail, avail))
diag(rg) <- 1
for (f in list.files(RG, pattern = "^ma_.*\\.log$", full.names = TRUE)) {
  ln <- readLines(f); s <- grep("Summary of Genetic Correlation", ln)
  if (!length(s)) next
  for (row in ln[(s + 2):length(ln)]) {
    p <- strsplit(trimws(row), " +")[[1]]; if (length(p) < 3 || !grepl("sumstats", p[1])) next
    a <- sub(".*/", "", sub(".sumstats.gz", "", p[1])); b <- sub(".*/", "", sub(".sumstats.gz", "", p[2]))
    v <- suppressWarnings(as.numeric(p[3]))
    if (a %in% avail && b %in% avail && !is.na(v)) { rg[a, b] <- rg[b, a] <- max(min(v, 1), -1) }
  }
}
cat(sprintf("[gwas] %d traits; rg non-NA off-diag: %d/%d\n", length(avail),
            sum(!is.na(rg[lower.tri(rg)])), sum(lower.tri(rg))))

# ── build MATRIX-ONLY genetic-correlation figure ──
# (obs->liability h2 lollipop dropped: deterministic scale conversion, not an informative
#  "captured" gap. The cross-trait rg matrix is the legitimate, example-faithful content.)
gl <- c("Liver disease", "Liver enzymes", "Metabolic")
gl <- gl[gl %in% meta$group]
group_of <- setNames(meta$group, meta$display)
rgd <- rg[avail, avail]; dimnames(rgd) <- list(meta$display, meta$display)
ord_keys <- order_by_group_then_clust(rgd, group_of, gl)
mat <- pie_glyph_matrix(rgd, ord_keys, group_of, group_levels = gl, italic_items = FALSE)

n_est <- sum(!is.na(rgd[lower.tri(rgd)])); n_tot <- sum(lower.tri(rgd))
# standalone matrix: move legends into the empty upper-right triangle + tighten margins
p_gwas <- mat$plot + theme(legend.position = c(0.82, 0.72), legend.background = element_blank(),
                           plot.margin = margin(20, 6, 4, 30))
save_fig(p_gwas, file.path(PANEL_DIR, "composite_gwas.pdf"), width = 4.4, height = 4.0)
cat(sprintf("[gwas] matrix-only saved; %d/%d off-diag estimable\n", n_est, n_tot))
message(sprintf(paste0("CAPTION (GWAS genetic-architecture): LDSC cross-trait genetic correlation (rg) ",
  "among %d distinct EUR liver-disease GWAS (one well-powered GWAS per trait; pie fill proportional to ",
  "|rg|, blue +/red -), grouped by category; %d/%d off-diagonal pairs estimable. All sumstats munged with ",
  "a common HapMap3 --merge-alleles reference. Observed-scale SNP h2 (sidecar): NAFLD %.3f / cirrhosis ",
  "%.3f / ALT %.3f / AST %.3f / obesity %.3f (biomarker traits far exceed rare EHR-defined diseases)."),
  length(avail), n_est, n_tot, meta[key=="ghodsian_nafld"]$h2_observed, meta[key=="ghouse_cirr"]$h2_observed,
  meta[key=="ukbb_alt"]$h2_observed, meta[key=="ukbb_ast"]$h2_observed, meta[key=="finngen_obesity"]$h2_observed))
fwrite(meta[, .(key, display, group, h2_observed = round(h2_observed, 4))],
       file.path(DATA_DIR, "composite_gwas_h2.csv"))
fwrite(as.data.table(rgd, keep.rownames = "trait"), file.path(DATA_DIR, "composite_gwas_rg.csv"))
