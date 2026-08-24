#!/usr/bin/env Rscript
# figS_mesusie_shared.R  (2026-07-05)  — Fig 2 supplement (defends para 4)
# Multi-ancestry joint fine-mapping (meSuSiE) assigns credible sets to
# model-shared or ancestry-restricted components. These are model assignments,
# not proof of ancestry-specific biology. All numbers are read from disk.
#
# SOURCE (2026-07-05): repointed from the retired enzyme-only EUR/EAS 2-way run
# (mesusie/mesusie_locus_summary.csv) to the within-MVP N-way cross-ancestry run
# (mesusie_mvp/mesusie_locus_summary_mvp.csv; EUR/AFR/AMR/EAS, 7 MVP traits). The
# per-credible-set ancestry composition is now encoded as one n_cs_<combo> column
# per ancestry combination (e.g. n_cs_EUR, n_cs_AFR, n_cs_AMR = single-ancestry;
# n_cs_EUR_AFR, n_cs_EUR_AFR_AMR, ... = multi-ancestry "shared"). We split them
# programmatically by token count so the panel is robust to which combos appear.
#
# Out: figures/main/fig2_genetics/panels/FigS2M_mesusie_shared_specific.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- Sys.getenv("FIG2_SUPP_OUT_DIR", unset = file.path(FIG3_DIR, "panels"))
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

ml <- fread(file.path(BASE, "GWAS/finemapping/results/mesusie_mvp/mesusie_locus_summary_mvp.csv"))
ml <- ml[converged == TRUE]

# --- split n_cs_<combo> columns into ancestry-restricted vs model-shared
ncs_cols  <- grep("^n_cs_", names(ml), value = TRUE)
n_token   <- vapply(sub("^n_cs_", "", ncs_cols),
                    function(x) length(strsplit(x, "_")[[1]]), integer(1))
spec_cols <- ncs_cols[n_token == 1]        # e.g. n_cs_EUR / n_cs_AFR / n_cs_AMR
shar_cols <- ncs_cols[n_token >= 2]        # any credible set spanning >= 2 ancestries

# per-locus shared/specific CS counts (for the loci-level summary)
ml[, shared_n := rowSums(as.matrix(.SD), na.rm = TRUE), .SDcols = shar_cols]
ml[, spec_n   := rowSums(as.matrix(.SD), na.rm = TRUE), .SDcols = spec_cols]
nloci        <- nrow(ml)
nshared_loci <- sum(ml$shared_n > 0)
amr_cols     <- ncs_cols[grepl("AMR", ncs_cols)]
ml[, amr_n := rowSums(as.matrix(.SD), na.rm = TRUE), .SDcols = amr_cols]
n_amr_loci   <- sum(ml$amr_n > 0)

# credible-set totals: one aggregate "Shared" bar + one bar per single ancestry
sh        <- sum(unlist(ml[, ..shar_cols]), na.rm = TRUE)
spec_tot  <- vapply(spec_cols, function(c) sum(ml[[c]], na.rm = TRUE), numeric(1))
spec_anc  <- sub("^n_cs_", "", spec_cols)                       # EUR / AFR / AMR
tot       <- sh + sum(spec_tot)

dec <- rbindlist(c(
  list(data.table(type = "Model-shared (≥2 ancestries)", n = sh, key_anc = "Shared")),
  lapply(seq_along(spec_cols), function(i)
    data.table(type = paste0(spec_anc[i], "-restricted"), n = spec_tot[i], key_anc = spec_anc[i]))
))
dec[, pct := 100 * n / tot]
# order: single-ancestry bars first (ascending n), Shared last -> plotted at the
# TOP (ggplot maps the last factor level to the top of a horizontal bar).
dec[, is_shared := as.integer(key_anc == "Shared")]
setorder(dec, is_shared, n)
dec[, type := factor(type, levels = type)]
dec[, is_shared := NULL]

# colours: shared = teal; single-ancestry bars = ANCESTRY_COLORS (EUR/AFR/AMR/EAS)
dec_cols <- c("Shared" = "#00695C", ANCESTRY_COLORS)
names(dec_cols)[1] <- "Shared"
fill_map <- setNames(dec_cols[dec$key_anc], dec$type)

p <- ggplot(dec, aes(n, type, fill = type)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%d  (%.0f%%)", n, pct)), hjust = -0.12, size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_manual(values = fill_map, guide = "none") +
  scale_x_continuous(limits = c(0, 1.18 * max(dec$n)), expand = expansion(mult = c(0, 0))) +
  labs(x = "meSuSiE credible sets", y = NULL) +
  theme_masld(base_size = 9) +
  theme(axis.text.y = element_text(size = 6))

save_fig(p, file.path(PANEL_DIR, "FigS2M_mesusie_shared_specific.pdf"),
         width = fig_col_width * 1.05, height = 2.2)

fwrite(dec[order(-n), .(type, n, pct = round(pct, 1))],
       file.path(PANEL_DIR, "FigS2M_mesusie_shared_specific_source.csv"))
message("[caption] meSuSiE model-assigned shared and ancestry-restricted credible sets; these categories do not prove ancestry-specific biology.")
cat(sprintf("[figS meSuSiE MVP] shared CS=%d (%.0f%%); specific CS: %s\n",
            sh, 100 * sh / tot,
            paste(sprintf("%s=%d", spec_anc, spec_tot), collapse = " ")))
cat(sprintf("[figS meSuSiE MVP] loci: %d total; %d with >=1 shared CS; %d only-specific; %d involve AMR\n",
            nloci, nshared_loci, sum(ml$shared_n == 0 & ml$spec_n > 0), n_amr_loci))
