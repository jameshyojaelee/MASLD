# H86 audit: recompute writer remap counts (symbols only; no evidence class read).
set.seed(20260923)
suppressPackageStartupMessages({library(data.table); library(HGNChelper); library(org.Hs.eg.db); library(limma)})
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
A <- file.path(root, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H86")
g <- fread(file.path(A, "tmp/catalog_symbols_noclass.tsv"))
p <- fread(file.path(root, "Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv"), na.strings = "")
cat("BH recompute max abs diff (R)", max(abs(p.adjust(p$pvalue, "BH") - p$padj)), "\n")
cat("HGNChelper", as.character(packageVersion("HGNChelper")), "org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db")), "\n")
um <- setdiff(p$gene, g$symbol); cat("unmatched", length(um), "\n")
h <- as.data.table(suppressWarnings(checkGeneSymbols(um, unmapped.as.na = TRUE)))
h[, n_sugg := ifelse(is.na(Suggested.Symbol), 0L, lengths(strsplit(Suggested.Symbol, " /// ", fixed = TRUE)))]
h[, in_cat := n_sugg == 1L & Suggested.Symbol %in% g$symbol]
h[, target_direct := in_cat & Suggested.Symbol %in% p$gene]
cat("single", sum(h$n_sugg == 1), "ambiguous", sum(h$n_sugg > 1), "none", sum(h$n_sugg == 0),
    "single_in_catalog", sum(h$in_cat), "target_already_direct", sum(h$target_direct),
    "dup_targets", sum(duplicated(h[in_cat == TRUE]$Suggested.Symbol)),
    "single_is_self", sum(h$n_sugg == 1 & h$Suggested.Symbol == h$x), "\n")
cat("ambiguous:", h[n_sugg > 1]$x, "\n")
a <- data.table(x = um, org = alias2SymbolTable(um, species = "Hs"))
m <- merge(h, a, by = "x")
both <- m[n_sugg == 1 & !is.na(org)]
cat("agree", sum(both$Suggested.Symbol == both$org), "/", nrow(both), "\n")
print(both[Suggested.Symbol != org, .(x, hgnc = Suggested.Symbol, org)])
print(h[x %in% c("MARC1","MARC2","AARS","CD97","CERT","ATP5MD"), .(x, Suggested.Symbol, in_cat)])
# NC-2 ceiling: remap can add at most sum(in_cat) genes to any class
cat("max untested pct after remap", round(100 * (409 + sum(h$in_cat)) / 13551, 2), "\n")
# HGNChelper map date
md <- grep("updated", capture.output(print(attributes(HGNChelper::hgnc.table))), value = TRUE)
cat(head(md, 2), sep = "\n")
# org.Hs.eg.db: alias -> Ensembl available for the remap targets (feasibility only; no class read)
tg <- h[in_cat == TRUE]
ens <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = tg$x, keytype = "ALIAS", columns = "ENSEMBL"))
setDT(ens)
cat("remapped aliases with >=1 org.Hs.eg.db Ensembl ID", ens[!is.na(ENSEMBL), uniqueN(ALIAS)], "/", nrow(tg), "\n")
