#!/usr/bin/env Rscript
# fig2_finemap_headtohead.R  (2026-06-23)
# Companion to fig2_finemap_resolution.R (the own-denominator 184/132/139/119 bars).
# Here every tool is restricted to the COMMON locus set it CAN actually run on, so
# all four bars are the same length and only the PIP-resolution colors change — the
# honest head-to-head the own-denominator figure could not give.
#
# Common set = physical regions where all four fine-mappers produced a result.
# It is bounded by the cross-ancestry tools (SuSiEx/meSuSiE need >=2 ancestries, so
# only the multi-ancestry MVP loci are eligible). The cross-ancestry arms are the
# within-MVP N-way run (SuSiEx + meSuSiE; EUR/AFR/AMR/EAS; 7 MVP traits, 411 shared
# loci) that REPLACES the retired enzyme-only EUR/EAS 2-way run. Computed live from
# on-disk fine-mapping outputs (no hard-coded counts):
#   merged_loci_map.csv  carma_all_results.csv
#   susiex_mvp/susiex_cs_mvp.csv   mesusie_mvp/mesusie_locus_summary_mvp.csv
#
# Out: figures/main/fig2_genetics/panels/FigS2C_finemap_method_headtohead.pdf
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")
FM        <- file.path(BASE, "GWAS/finemapping/results")
TOL       <- 500000L   # physical-region match window

# ---- within-ancestry independent loci (+ CARMA coverage) -------------------
lm  <- fread(file.path(FM, "merged_loci_map.csv"))
ind <- lm[kept == TRUE, .(chr = as.integer(chr), pos = as.integer(bp),
                          susie_pip = as.numeric(max_pip), merged_locus)]
ca  <- fread(file.path(FM, "carma_all_results.csv"))
key <- lm[, .(study, original_locus, merged_locus)]
ca  <- merge(ca, key, by.x = c("study","locus"), by.y = c("study","original_locus"), all.x = TRUE)
carma_cov <- ca[merged_locus %in% ind$merged_locus, .(carma_pip = max(PIP, na.rm = TRUE)), by = merged_locus]
ind <- merge(ind, carma_cov, by = "merged_locus", all.x = TRUE)

# ---- cross-ancestry scaffold (the limiting tools) --------------------------
sx  <- fread(file.path(FM, "susiex_mvp/susiex_cs_mvp.csv"))
sxl <- sx[, .(pip = max(OVRL_PIP, na.rm = TRUE)), by = locus_id]
# MVP locus_ids are "locus_MVP_<trait>_chr<C>_<POS>" (two tokens before chr) — parse
# chr/pos robustly with regexec (a fixed anchored sub() silently fails on this form).
sxm <- regmatches(sxl$locus_id, regexec("chr([0-9XY]+)_([0-9]+)$", sxl$locus_id))
sxl[, `:=`(chr = as.integer(sapply(sxm, function(x) if (length(x) == 3) x[2] else NA)),
           pos = as.integer(sapply(sxm, function(x) if (length(x) == 3) x[3] else NA)),
           tool = "SuSiEx")]
me  <- fread(file.path(FM, "mesusie_mvp/mesusie_locus_summary_mvp.csv"))[converged == TRUE]
mel <- me[, .(pip = as.numeric(max_pip), chr = as.integer(chr),
              pos = as.integer((window_start + window_end)/2), tool = "meSuSiE")]
xa  <- rbindlist(list(sxl[,.(chr,pos,pip,tool)], mel[,.(chr,pos,pip,tool)]))[order(chr,pos)]
xa[, region := { r<-integer(.N); c<-0L; last<--Inf
  for(i in seq_len(.N)){ if(pos[i]-last>TOL) c<-c+1L; r[i]<-c; last<-pos[i] }; r }, by=chr]
xa[, region := paste0(chr,"_",region)]
reg <- xa[, .(chr = chr[1], pos_lo = min(pos), pos_hi = max(pos),
              susiex_pip  = max(c(-Inf, pip[tool=="SuSiEx"])),
              mesusie_pip = max(c(-Inf, pip[tool=="meSuSiE"]))), by = region]
reg[susiex_pip  == -Inf, susiex_pip  := NA_real_]
reg[mesusie_pip == -Inf, mesusie_pip := NA_real_]

# ---- match within-ancestry loci to each cross-ancestry region --------------
reg[, `:=`(susie_pip = NA_real_, carma_pip = NA_real_)]
for (i in seq_len(nrow(reg))) {
  cand <- ind[chr == reg$chr[i] & pos >= reg$pos_lo[i]-TOL & pos <= reg$pos_hi[i]+TOL]
  if (nrow(cand)) {
    cand <- cand[order(abs(pos - (reg$pos_lo[i]+reg$pos_hi[i])/2))][1]
    reg$susie_pip[i] <- cand$susie_pip; reg$carma_pip[i] <- cand$carma_pip
  }
}
tcols  <- c("susie_pip","carma_pip","susiex_pip","mesusie_pip")
reg[, n_tools := rowSums(!is.na(.SD)), .SDcols = tcols]
common <- reg[n_tools == 4]
NTOT <- nrow(common)

# ---- per-tool resolution bands on the common set ---------------------------
band <- function(p) factor(fifelse(p>=0.9,"PIP ≥ 0.9", fifelse(p>=0.5,"PIP 0.5–0.9","PIP < 0.5")),
                           levels = c("PIP ≥ 0.9","PIP 0.5–0.9","PIP < 0.5"))
lab  <- c(susie_pip="SuSiE", carma_pip="CARMA", susiex_pip="SuSiEx", mesusie_pip="meSuSiE")
dB <- rbindlist(lapply(tcols, function(t)
  data.table(tool = lab[t], band = band(common[[t]]))))[, .N, by = .(tool, band)]
dB <- merge(CJ(tool = unname(lab), band = levels(band(0))), dB, by = c("tool","band"), all.x = TRUE)
dB[is.na(N), N := 0L]
dB[, tool := factor(tool, levels = rev(c("SuSiE","CARMA","SuSiEx","meSuSiE")))]   # SuSiE on top
dB[, band := factor(band, levels = c("PIP < 0.5","PIP 0.5–0.9","PIP ≥ 0.9"))]   # low PIP at left, >=0.9 at right
pct <- dB[band == "PIP ≥ 0.9", .(tool, x = NTOT, lab = paste0(round(100*N/NTOT),"%"))]  # % resolved per bar end

colB <- c("PIP ≥ 0.9" = "#0D3B66", "PIP 0.5–0.9" = "#7FB3D5", "PIP < 0.5" = "#D9DCE0")

p <- ggplot(dB, aes(N, tool, fill = band)) +
  geom_col(width = 0.66, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(N > 0, N, "")), position = position_stack(vjust = 0.5, reverse = TRUE),
            size = GEOM_TEXT_6PT, color = "black") +
  geom_text(data = pct, aes(x = x, y = tool, label = lab), inherit.aes = FALSE,
            hjust = -0.25, size = GEOM_TEXT_6PT, fontface = "plain", color = "black") +
  scale_fill_manual(values = colB, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16)), limits = c(0, NTOT),
                     breaks = seq(0, NTOT, 10)) +
  labs(x = "Loci", y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.position = "bottom", legend.text = element_text(size = 6, color = "black"),
        legend.key.size = unit(0.32, "cm"))

save_fig(p, file.path(PANEL_DIR, "FigS2C_finemap_method_headtohead.pdf"),
         width = fig_col_width * 1.15, height = 2.6)

# ---- caption + numbers to stdout (NOT on the plot) -------------------------
cat(sprintf("[finemap_headtohead] wrote panel — common set N = %d physical loci\n", NTOT))
cat("Per-tool best-PIP resolution on the common set:\n")
print(dcast(dB, tool ~ band, value.var = "N")[order(match(tool, c("SuSiE","CARMA","SuSiEx","meSuSiE")))])
cat(sprintf(
"\nCAPTION: Fine-mapping resolution head-to-head on the %d physical loci that all four\nmethods can analyse. SuSiE and CARMA are within-ancestry; SuSiEx and meSuSiE are\nthe within-MVP N-way cross-ancestry run (EUR/AFR/AMR/EAS, 7 MVP traits) and can only\nrun where >=2 ancestries have data (the multi-ancestry MVP loci), which bounds the\ncommon set. Each bar holds the same %d loci; colour gives the best variant posterior\ninclusion probability (PIP) in the locus. The bold value is the share resolved to a\nsingle high-confidence variant (PIP>=0.9). Cross-ancestry fine-mapping raises\nhigh-confidence resolution on the same loci (SuSiEx %s / meSuSiE %s vs SuSiE %s /\nCARMA %s).\n",
  NTOT, NTOT,
  pct[tool=="SuSiEx", lab], pct[tool=="meSuSiE", lab],
  pct[tool=="SuSiE", lab], pct[tool=="CARMA", lab]))
