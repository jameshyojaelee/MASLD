# figS05 A: LD-corrected vs uncorrected enrichment of fine-mapped GWAS variants in cell-type ATAC peaks
# Side-by-side bars sorted by LD-corrected FE descending.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

LD_F     <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/enrichment_ld_null.csv")
FISHER_F <- file.path(BASE, "GWAS/finemapping/results/gwas_atac/enrichment_statistics.csv")

ld     <- fread(LD_F)
fisher <- fread(FISHER_F)

# canonical display names (NK and Plasma excluded)
ct_map <- c(
  "Hepatocytes"       = "Hep",
  "Macrophages"       = "Mac",
  "Fibroblasts"       = "Stellate",
  "Endothelial_cells" = "Endo",
  "Cholangiocytes"    = "Chol",
  "T_cells"           = "T"
)

ld[, ct := ct_map[cell_type]]
ld <- ld[!is.na(ct)]
ld[, fe    := 2 ^ log2FE]
ld[, fe_lo := 2 ^ ci_lo]
ld[, fe_hi := 2 ^ ci_hi]

fisher[, ct := ct_map[cell_type]]
fisher <- fisher[!is.na(ct)]
fisher[, fe_fisher := fold_enrichment]

# sort cell types by LD-corrected FE descending
setorder(ld, -fe)
ct_order <- ld$ct

merged <- merge(
  ld[, .(ct, fe_ld = fe, fe_lo, fe_hi, p_emp_bh)],
  fisher[, .(ct, fe_fisher, fisher_padj)],
  by = "ct"
)

# long format
long <- rbind(
  merged[, .(ct, method = "Uncorrected",  fe = fe_fisher,
             lo = NA_real_, hi = NA_real_, padj = fisher_padj)],
  merged[, .(ct, method = "LD-corrected", fe = fe_ld,
             lo = fe_lo,    hi = fe_hi,   padj = p_emp_bh)]
)
long <- long[!is.na(fe)]
long[, ct     := factor(ct, levels = ct_order)]
long[, method := factor(method, levels = c("Uncorrected", "LD-corrected"))]

bar_colors <- c("Uncorrected" = "#1565C0", "LD-corrected" = "#C9265E")

p <- ggplot(long, aes(x = ct, y = fe, fill = method)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#9E9E9E", linewidth = 0.35) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65,
           color = "black", linewidth = 0.25) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(width = 0.72),
                width = 0.18, linewidth = 0.3, color = "black",
                na.rm = TRUE) +
  scale_fill_manual(values = bar_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL,
       y = "Fold enrichment\n(fine-mapped variants in cell-type ATAC peaks)") +
  theme_masld(base_size = 9) +
  theme(
    axis.title.y    = element_text(size = 6, face = "plain", color = "black",
                                   margin = margin(r = 4)),
    axis.text.x     = element_text(size = 6, angle = 30, hjust = 1, color = "black"),
    axis.text.y     = element_text(size = 6, color = "black"),
    legend.text     = element_text(size = 6),
    legend.position = c(0.85, 0.92),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.35, "cm")
  )
message("[caption] GWAS variant enrichment in cell-type ATAC peaks")

out_pdf <- file.path(FIGS05_DIR, "figS05_a_ld_aware_null_comparison.pdf")
dir.create(FIGS05_DIR, showWarnings = FALSE, recursive = TRUE)
ggsave(out_pdf, p,
       width  = 120 / 25.4,
       height =  78 / 25.4,
       units  = "in",
       device = cairo_pdf)
message("Wrote: ", out_pdf)

out_csv <- sub("\\.pdf$", ".csv", out_pdf)
fwrite(long, out_csv)
message("Wrote: ", out_csv)
