#!/usr/bin/env Rscript
# KEY MESSAGE (Fig 4, Panel D — DESCRIPTIVE / EXPLORATORY COMPANION):
# Donor-level (n=18) hepatocyte chromVAR motif-accessibility effect sizes for the
# canonical metabolic nuclear receptors (FXR/NR1H4, LXR/NR1H3+NR1H2, PPARA) and
# the lipogenic + hypoxia TFs (ChREBP/MLXIPL, HIF1A/HIF3A). This panel pairs with
# the well-powered bulk (Fig 3, n=846) and GWAS/COLOC (Fig 2) evidence; the
# single-cohort multiome is UNDERPOWERED, so we show ONLY the effect-size
# direction (logFC of motif accessibility, disease vs control), with NO
# significance stars. An explicit on-figure note states 0 motifs pass FDR<0.05.
# This avoids any pseudoreplicated / overclaimed n=18 disease-DE p-value claim.
#
# Source (READ, never hardcoded): donor-level limma chromVAR per cell type
#   Analysis/ATAC/Human_Multiome/results/chromvar_v2/chromvar_limma_per_ct.csv
#   (cell_type, TF, logFC, P.Value, adj.P.Val, n_donors; hepatocyte rows only).
#
# Exact hepatocyte values pulled from the CSV (logFC / adj.P.Val), for the
# JASPAR monomer motifs of each requested factor:
#   Nr1H4 (FXR)        logFC = -0.225853  adj.P.Val = 0.929998   (closing)
#   Nr1h3 (LXRa)       logFC = -0.246356  adj.P.Val = 0.929998   (closing)
#   Nr1H2 (LXRb)       logFC = -0.225853  adj.P.Val = 0.929998   (closing)
#   Ppara (PPARA)      logFC = -0.183877  adj.P.Val = 0.929998   (closing)
#   MLXIPL (ChREBP)    logFC = +0.179249  adj.P.Val = 0.929998   (opening)
#   HIF1A              logFC = +0.077981  adj.P.Val = 0.929998   (opening)
#   HIF3A              logFC = +0.180179  adj.P.Val = 0.929998   (opening)
# Minimum adj.P.Val across ALL 993 hepatocyte motifs = 0.93 -> 0 pass FDR<0.05.
#
# Output: figures/main/fig4_validation/fig4d_chromvar_metabolic_nr.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data: donor-level hepatocyte chromVAR (READ from canonical CSV) ───────────
cv <- read.csv(file.path(BASE,
  "Analysis/ATAC/Human_Multiome/results/chromvar_v2/chromvar_limma_per_ct.csv"),
  stringsAsFactors = FALSE)

hep <- cv[cv$cell_type == "Hepatocytes", ]
stopifnot(nrow(hep) > 0)

# Canonical metabolic-NR + HIF motifs. `motif` = exact JASPAR id as it appears in
# the CSV (case-sensitive); `label` = display name; `factor` = pathway grouping.
sel <- data.frame(
  motif = c("Nr1H4", "Nr1h3", "Nr1H2", "Ppara", "MLXIPL", "HIF1A", "HIF3A"),
  label = c("NR1H4 (FXR)", "NR1H3 (LXRa)", "NR1H2 (LXRb)", "PPARA",
            "MLXIPL (ChREBP)", "HIF1A", "HIF3A"),
  group = c("Nuclear receptor", "Nuclear receptor", "Nuclear receptor",
            "Nuclear receptor", "Lipogenic /\nhypoxia", "Lipogenic /\nhypoxia",
            "Lipogenic /\nhypoxia"),
  stringsAsFactors = FALSE
)

df <- merge(sel, hep[, c("TF", "logFC", "adj.P.Val", "n_donors")],
            by.x = "motif", by.y = "TF", all.x = TRUE, sort = FALSE)
# Guard: every requested motif must resolve to exactly one hepatocyte row.
if (any(is.na(df$logFC)) || nrow(df) != nrow(sel))
  stop("Motif lookup failed; check JASPAR ids vs chromVAR CSV.")

n_don <- unique(df$n_donors)
stopifnot(length(n_don) == 1L)             # all rows share donor count (=18)
min_fdr <- min(hep$adj.P.Val, na.rm = TRUE)  # smallest FDR across all hep motifs
n_pass  <- sum(hep$adj.P.Val < 0.05, na.rm = TRUE)
stopifnot(n_pass == 0)                       # honesty guard: nothing significant

# Direction = sign of the accessibility effect (disease vs control), NOT a
# significance call. Closing (down) = blue, opening (up) = magenta; both muted.
df$direction <- ifelse(df$logFC < 0, "Closing (down)", "Opening (up)")
dir_cols <- c("Closing (down)" = "#518dc9",   # Liang medium blue
              "Opening (up)"   = "#C9265E")    # Liang deep magenta

# Order: nuclear receptors first, then lipogenic/hypoxia; within group by logFC.
df$group <- factor(df$group, levels = c("Nuclear receptor", "Lipogenic /\nhypoxia"))
df <- df[order(df$group, df$logFC), ]
df$label <- factor(df$label, levels = df$label)

# Value labels placed just OUTSIDE the bar end, away from 0, all-black text.
lab_pad <- 0.012
df$hjust <- ifelse(df$logFC < 0, 1.12, -0.12)
df$vlab  <- sprintf("%+.2f", df$logFC)

xr <- max(abs(df$logFC)) + 0.10

note <- sprintf(paste0("n=%d donors, donor-level, exploratory; ",
                       "0/%d motifs pass FDR<0.05 (min FDR=%.2f)"),
                n_don, nrow(hep), min_fdr)

p <- ggplot(df, aes(x = logFC, y = label, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_col(width = 0.66) +
  geom_text(aes(label = vlab, hjust = hjust), size = PUB_GEOM_TEXT,
            color = "black") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = dir_cols, name = NULL) +
  scale_x_continuous(limits = c(-xr, xr),
                     breaks = c(-0.2, 0, 0.2),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "chromVAR motif accessibility log2FC (MASLD vs control)",
       y = NULL, title = "Hepatocyte metabolic-TF motif accessibility",
       caption = note) +
  theme_masld() + theme_pub() +
  theme(
    axis.text.y  = element_text(face = "italic", color = "black"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, face = "bold", color = "black"),
    legend.position = "top",
    legend.key.size = PUB_LEGEND_KEY,
    plot.margin  = margin(3, 6, 3, 4),
    plot.caption = element_text(size = PUB_SUBTITLE - 1, color = "gray35",
                                hjust = 0, margin = margin(t = 4))
  )

out <- file.path(FIG4_DIR, "fig4d_chromvar_metabolic_nr.pdf")
pdf(out, width = fig_half_width, height = 2.35, useDingbats = FALSE)
print(p)
invisible(dev.off())
message("Saved: ", out)
message(sprintf("  %d motifs | donors=%d | min hep FDR=%.4f | n pass FDR<0.05=%d",
                nrow(df), n_don, min_fdr, n_pass))
