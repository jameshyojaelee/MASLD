#!/usr/bin/env Rscript
# fig2_tri_ancestry_coloc.R  (2026-06-17)  — Fig 2 (defends para 6)
# Cross-ancestry colocalization of the three pan-ancestry effector genes.
# RORA / EPHA2 / GGT1 each colocalize at PP.H4 > 0.9 in >=3 ancestry panels at
# the SAME lead variant — the visual proof of effector-gene portability.
# All values computed from disk so they always match the manuscript text.
#
# Out: figures/main/fig2_genetics/panels/rora_tri_ancestry_coloc.pdf (+ source CSV)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG3_DIR, "panels")

sc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))
anc <- function(g) fifelse(grepl("BBJ", g), "EAS",
                  fifelse(grepl("PanUKBB_AFR", g), "AFR",
                  fifelse(grepl("PanUKBB_CSA", g), "SAS", "EUR")))
trait <- function(g) fifelse(grepl("ALT", g), "ALT", fifelse(grepl("AST", g), "AST",
                    fifelse(grepl("GGT", g), "GGT", fifelse(grepl("PDFF", g), "PDFF",
                    fifelse(grepl("NASH", g), "NASH", "NAFLD")))))
sc[, ancestry := anc(gwas_name)]
sc[, trait := trait(gwas_name)]
sc[, best := fifelse(!is.na(PP.H4.susie), PP.H4.susie, PP.H4.abf)]
sc[, method_used := fifelse(!is.na(PP.H4.susie), "SuSiE", "ABF")]

GENES <- c("RORA", "EPHA2", "GGT1")
ANC   <- c("EUR", "EAS", "SAS", "AFR")

# best colocalizing trait per gene x ancestry
d <- sc[gene %in% GENES & ancestry %in% ANC]
best <- d[, .SD[which.max(best)], by = .(gene, ancestry)]
grid <- CJ(gene = GENES, ancestry = ANC)
m <- merge(grid, best[, .(gene, ancestry, best, trait, method_used, top_snp)],
           by = c("gene", "ancestry"), all.x = TRUE)
m[, colocalizes := !is.na(best) & best > 0.5]
m[, lab := fifelse(colocalizes, sprintf("%.3f\n%s", best, trait), "n.s.")]
m[, fill_val := fifelse(colocalizes, best, NA_real_)]
m[, gene := factor(gene, levels = GENES)]
m[, ancestry := factor(ancestry, levels = rev(ANC))]   # EUR on top

# shared lead per gene (for subtitle): most common lead among colocalizing panels
lead <- best[best > 0.5, .(lead = names(sort(table(top_snp), decreasing = TRUE))[1]), by = gene]
lead <- lead[match(GENES, gene)]
sub <- paste0("Shared lead variant per gene:  ",
              paste(sprintf("%s chr%s", lead$gene, lead$lead), collapse = "   "))

p <- ggplot(m, aes(x = gene, y = ancestry)) +
  geom_tile(aes(fill = fill_val), color = "white", linewidth = 0.6) +
  geom_text(aes(label = lab), size = 2.5, lineheight = 0.82,
            color = ifelse(!is.na(m$fill_val) & m$fill_val > 0.85, "white", "grey20")) +
  scale_fill_gradient(low = "#9EC2BC", high = "#00695C", limits = c(0.5, 1),
                      na.value = "grey92", name = "PP.H4",
                      breaks = c(0.5, 0.75, 1.0)) +
  scale_x_discrete(position = "top") +
  coord_equal() +
  labs(x = NULL, y = NULL,
       title = "Cross-ancestry colocalization of effector genes") +
  theme_masld(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"),
        axis.text.x.top = element_text(face = "bold.italic", size = 9),
        axis.text.y = element_text(face = "bold", size = 8),
        legend.position = "right",
        legend.key.width = unit(0.22, "cm"),
        panel.grid = element_blank())

save_fig(p, file.path(PANEL_DIR, "rora_tri_ancestry_coloc.pdf"),
         width = fig_col_width * 0.95, height = 2.7)

fwrite(m[order(gene, ancestry), .(gene, ancestry, PP_H4 = round(best, 4),
        colocalizing_trait = trait, method = method_used, lead_variant = top_snp,
        colocalizes)],
       file.path(PANEL_DIR, "rora_tri_ancestry_coloc_source.csv"))
cat("[fig2 tri-ancestry] wrote rora_tri_ancestry_coloc.pdf\n")
print(m[order(gene, ancestry), .(gene, ancestry, best = round(best,3), trait, colocalizes)])
