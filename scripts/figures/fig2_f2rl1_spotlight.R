#!/usr/bin/env Rscript
# fig2_f2rl1_spotlight.R — Fig 2 supplementary method-resolution example.
# F2RL1 / PAR2 colocalizes with serum GGT and is induced in the disease-state
# map. The comparison is multi-signal versus single-signal COLOC; neither model
# alone establishes a causal gene.
# All values read from disk.
#
# Out: figures/main/fig2_genetics/panels/FigS2L_f2rl1_spotlight.pdf (+ source CSV)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- Sys.getenv("FIG2_SUPP_OUT_DIR", unset = file.path(FIG3_DIR, "panels"))
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

sc  <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
deg <- fread(file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))

f <- sc[gene == "F2RL1"]
get1 <- function(gw, col) f[gwas_name == gw, get(col)][1]
co <- data.table(
  label  = c("EUR GGT\nsingle-signal", "EUR GGT\nmulti-signal", "SAS GGT\nsingle-signal"),
  pp4    = c(get1("UKBB_GGT", "PP.H4.abf"), get1("UKBB_GGT", "PP.H4.susie"),
             get1("PanUKBB_CSA_GGT", "PP.H4.abf")),
  method = c("single-signal COLOC", "multi-signal COLOC", "single-signal COLOC"))
co[, label := factor(label, levels = rev(label))]
mcols <- c("single-signal COLOC" = "#90A4AE", "multi-signal COLOC" = "#1565C0")

pA <- ggplot(co, aes(pp4, label, fill = method)) +
  geom_col(width = 0.66) +
  geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.3, color = "grey45") +
  geom_text(aes(label = sprintf("%.3f", pp4)), hjust = -0.15, size = GEOM_TEXT_6PT, color = "grey15") +
  scale_fill_manual(values = mcols, name = NULL, guide = "none") +
  scale_x_continuous(limits = c(0, 1.12), breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "PP.H4", y = NULL) +
  theme_masld(base_size = 9) +
  theme(axis.text.y = element_text(size = 6, lineheight = 0.8))

d <- deg[symbol == "F2RL1"]
ci <- 1.96 * d$SE[1]
ed <- data.table(lab = "F2RL1", lfc = d$logFC[1], lo = d$logFC[1] - ci, hi = d$logFC[1] + ci,
                 padj = d$padj[1])
ed[, lab_x := hi + 0.07]
pB <- ggplot(ed, aes(lfc, lab)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.10, color = "#C9265E", linewidth = 0.5) +
  geom_point(size = 3, color = "#C9265E") +
  geom_text(aes(x = lab_x, label = sprintf("+%.2f\npadj %.0e", lfc, padj)), hjust = 0,
            size = GEOM_TEXT_6PT, color = "grey15", lineheight = 0.85) +
  scale_x_continuous(limits = c(0, 1.5), breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression("Disease log"[2]*"FC"), y = NULL) +
  theme_masld(base_size = 9) +
  theme(axis.text.y = element_text(size = 6, face = "italic"))

p <- (pA | pB) + plot_layout(widths = c(1.85, 1))
message("[caption] F2RL1 / PAR2 -- left: colocalization (GGT); right: disease expression")

save_fig(p, file.path(PANEL_DIR, "FigS2L_f2rl1_spotlight.pdf"),
         width = fig_col_width * 1.25, height = 2.5)

fwrite(rbind(
  co[, .(evidence = "coloc", detail = gsub("\n", " ", label),
         signal_model = method, value = round(pp4, 4))],
  data.table(evidence = "DEG", detail = "disease vs control log2FC",
             signal_model = NA_character_, value = round(d$logFC[1], 3))),
  file.path(PANEL_DIR, "FigS2L_f2rl1_spotlight_source.csv"))
cat(sprintf("[fig2 F2RL1] wrote FigS2L_f2rl1_spotlight.pdf | EUR abf=%.3f susie=%.3f, SAS abf=%.3f | logFC=%.2f padj=%.1e\n",
            co$pp4[1], co$pp4[2], co$pp4[3], d$logFC[1], d$padj[1]))
message(sprintf("CAPTION: F2RL1/PAR2 is a GGT-colocalized candidate that is induced in the current disease-state map. Multi-signal COLOC resolves the EUR GGT signal to PP.H4 %.3f, whereas the single-signal model gives %.3f; disease log2FC +%.2f (padj %.0e). COLOC is shared-signal evidence, not causal-gene proof.",
                co$pp4[2], co$pp4[1], d$logFC[1], d$padj[1]))
