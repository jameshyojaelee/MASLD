#!/usr/bin/env Rscript
# F2abce candidate for Fig 2C (finding R02-07). Derived from
# scripts/figures/fig2_coding_noncoding.R (copy: fig2_coding_noncoding.orig.R).
# Membership, lead variant and fine_class are re-derived exactly as the r2
# producer does (r2 COLOC, main strata, best of SuSiE/ABF PP.H4 > 0.5, class
# from the r2 variant_classification.csv). New: the coding segment is split by
# whether the lead variant lies in the CDS or splice site of the listed gene
# itself or only in another gene's (TxDb hg19 knownGene GENEID).
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(20260923)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
WS   <- file.path(BASE, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F2abce")
OUTF <- file.path(BASE, "figures/candidates/review-r3-fixes-20261015T112148Z/F2abce")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/fig2_consequence_palette.R"))
source(file.path(WS, "scripts/lead_coding_owner.R"))
FM <- file.path(BASE, "GWAS/finemapping/results")
COLOC_INPUT <- file.path(FM, "susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv")
VC <- file.path(BASE, "figures/candidates/fig2-coloc-r2-20260924T053244Z/variant_class/variant_classification.csv")

MAIN <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[placement == "main", study_name]
av <- fread(COLOC_INPUT, select = c("gwas_name", "gene", "PP.H4.susie", "PP.H4.abf", "top_snp"))
av <- av[gwas_name %in% MAIN & !is.na(gene) & gene != ""]
av[, pp4_best := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
av[!is.finite(pp4_best), pp4_best := NA_real_]
setorder(av, gene, -pp4_best, gwas_name)
best <- av[, .SD[1], by = gene][is.finite(pp4_best) & pp4_best > 0.5]
vc <- unique(fread(VC, select = c("variant_key", "fine_class", "nearest_gene_symbol")), by = "variant_key")
best <- merge(best, vc, by.x = "top_snp", by.y = "variant_key", all.x = TRUE)
stopifnot(!anyNA(best$fine_class))
best[, fc := fine_class]
best[fc == "spliceSite", fc := "coding"]
best[fc %in% c("fiveUTR", "threeUTR"), fc := "UTR"]
cat("genes:", nrow(best), "\n"); print(best[, .N, by = fc][order(-N)])

ann <- annotate_leads(best$gene, best$top_snp)
best <- merge(best, ann, by.x = c("gene", "top_snp"), by.y = c("gene", "variant_key"), all.x = TRUE)
cat("fine_class vs fresh TxDb class disagreements:", best[fine_class != txdb_class, .N], "\n")
best[, fc2 := fifelse(fc == "coding" & own_gene_coding, "coding_this_gene",
               fifelse(fc == "coding", "coding_other_gene", fc))]
print(best[, .N, by = fc2][order(-N)])
cod <- best[fc == "coding", .(gene, gwas_name, top_snp, fine_class, pp4_best = signif(pp4_best, 4),
                              coding_genes_at_lead = coding_symbols, own_gene_coding,
                              nearest_gene_symbol)][order(own_gene_coding, top_snp)]
print(cod)
cat("distinct coding lead variants:", uniqueN(cod$top_snp), "\n")
fwrite(cod, file.path(OUTF, "Fig2C_coding_lead_gene_ownership.tsv"), sep = "\t")

counts <- best[, .(n_genes = .N), by = fc2]
stack_order <- c("intergenic", "intron", "promoter", "UTR", "coding_other_gene", "coding_this_gene")
counts <- counts[match(stack_order, fc2)][!is.na(fc2)]
counts[, ymax := cumsum(n_genes)][, ymin := c(0, head(ymax, -1))][, ymid := (ymin + ymax) / 2]
pal <- c(intergenic = FIG2_CONSEQUENCE_PALETTE[["intergenic"]], intron = FIG2_CONSEQUENCE_PALETTE[["intron"]],
         promoter = FIG2_CONSEQUENCE_PALETTE[["promoter"]], UTR = FIG2_CONSEQUENCE_PALETTE[["UTR"]],
         coding_other_gene = "#E8A3BC", coding_this_gene = FIG2_CONSEQUENCE_PALETTE[["coding"]])
nm <- c(intergenic = "intergenic", intron = "intron", promoter = "promoter", UTR = "UTR",
        coding_other_gene = "coding, other gene", coding_this_gene = "coding, this gene")
counts[, fill_col := unname(pal[fc2])]
counts[, label := paste0(nm[fc2], " (", n_genes, ")")]
ntot <- sum(counts$n_genes)
label_y <- c(intergenic = 171, intron = 570, promoter = 820, UTR = 890, coding_other_gene = 955, coding_this_gene = 1020)
counts[, label_y := unname(label_y[fc2])]
BAR_LEFT <- 0.06; BAR_RIGHT <- 0.26; LABEL_X <- 0.36
p <- ggplot() +
  geom_rect(data = counts, aes(xmin = BAR_LEFT, xmax = BAR_RIGHT, ymin = ymin, ymax = ymax, fill = fill_col),
            color = "white", linewidth = 0.18) +
  geom_segment(data = counts, aes(x = BAR_RIGHT, xend = LABEL_X - 0.02, y = ymid, yend = label_y),
               color = "#666666", linewidth = 0.08) +
  geom_text(data = counts, aes(x = LABEL_X, y = label_y, label = label), hjust = 0, vjust = 0.5,
            size = GEOM_TEXT_6PT, family = "Helvetica", color = "black") +
  annotate("text", x = BAR_LEFT, y = 1090, label = paste0(format(ntot, big.mark = ","), " colocalizing genes"),
           hjust = 0, vjust = 0.5, size = GEOM_TEXT_6PT, family = "Helvetica", color = "black") +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, 1120), expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") + theme_void(base_size = 6) +
  theme(text = element_text(family = "Helvetica", size = 6, color = "black"),
        legend.position = "none", plot.margin = margin(3, 2, 3, 2, "pt"))
ggsave(file.path(OUTF, "Fig2C_coding_noncoding_split_ownership.pdf"), p,
       width = 1.30, height = 2.35, units = "in", device = cairo_pdf)
fwrite(counts[, .(consequence = fc2, n_genes, method_total = ntot)],
       file.path(OUTF, "Fig2C_coding_noncoding_ownership_source.csv"))
fwrite(best[, .(gene, gwas_name, top_snp, pp4_best, fine_class, txdb_class, consequence_panel = fc2,
                coding_genes_at_lead = coding_symbols, own_gene_coding, nearest_gene_symbol)],
       file.path(OUTF, "Fig2C_gene_level_lead_annotation.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(WS, "sessionInfo_fig2c.txt"))
cat("done\n")
