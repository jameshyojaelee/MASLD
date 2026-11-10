# DIG1 job 26 (exploratory; follows jobs 22, 22b, 24, 25). OXPHOS protein falls with ballooning
# within MASLD while its RNA does not. Which liver proteins rise or fall with the Complex I deficit,
# and do they point to a removal route (proteasome, immunoproteasome, mitophagy, lysosome) or to
# a cell source (macrophages, neutrophils)? PXD051911, 46 MASLD, same matrix and covariates as job 22.
# Score: CI - CII (job 22 participant_scores). Per protein, within MASLD, partial Spearman with the
# score (batch + age + BMI + sex), and again with ballooning added. OXPHOS, mitoribosome and nucleoid
# proteins are excluded from the partner table (they define or co-define the score).
# Set tests: for each fixed set and each Hallmark set, Wilcoxon of member rho vs all other partner
# proteins (sets with >= 5 measured members); BH over all set tests per adjustment.
# Usage: Rscript 26_ci_loss_partners.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(limma)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z")
raw <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt"))
md <- fread(file.path(ROOT, "data/PXD051911/meta_data.txt"))
samples <- setdiff(names(raw), c("ProteinAccessions","Genes","ProteinDescriptions"))
md <- unique(md[liver_proteomics_filename %in% samples], by = "liver_proteomics_filename")
stopifnot(nrow(md) == 58L)
md[, `:=`(sample_id = liver_proteomics_filename, masld = as.integer(saf_diagnosis != "No_MASLD"),
          batch = factor(ifelse(grepl("^2019", liver_proteomics_filename), "2019", "2020")),
          age = as.numeric(alder), bmi_num = as.numeric(bmi), sex = factor(gender),
          Fibrosis = as.numeric(sub("^F", "", kleiner_fibrosis_grade)), NAS = as.numeric(nafld_activity_score),
          Steatosis = as.numeric(steatosis_score), Ballooning = as.numeric(hepatocellular_ballooning_score),
          Inflammation = as.numeric(lobular_inflammation_score))]
setorder(md, sample_id)
raw <- raw[!is.na(Genes) & Genes != "" & !grepl(";", Genes, fixed = TRUE)]
x <- as.matrix(raw[, md$sample_id, with = FALSE]); storage.mode(x) <- "double"
x[!is.finite(x) | x <= 0] <- NA_real_; x <- log2(x); rownames(x) <- raw$Genes
x <- x[rowMeans(is.na(x)) <= .5, , drop = FALSE]
x <- normalizeBetweenArrays(x, method = "quantile")
idx <- split(seq_len(nrow(x)), rownames(x))
gene <- t(vapply(idx, function(i) if (length(i) == 1) x[i, ] else apply(x[i, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(x))))
colnames(gene) <- md$sample_id
g <- rownames(gene)
sc <- fread(file.path(RUN, "out_mito_pattern/participant_scores.tsv"))[, .(sample_id, CI_minus_CII)]
md <- merge(md, sc, by = "sample_id"); setorder(md, sample_id); stopifnot(identical(md$sample_id, colnames(gene)))
excl <- g[grepl("^NDUF|^SDH[ABCD]$|^UQCR|^CYC1$|^COX[4-8]|^ATP5|^MT-|^MRP[LS][0-9]", g) | g %in% c("TFAM", "SSBP1", "POLG", "LRPPRC")]
M <- md[masld == 1]; GM <- gene[setdiff(g, excl), M$sample_id]
pc <- function(v, vars) { z <- data.table(M, v = v); z <- z[complete.cases(z[, c("v", "CI_minus_CII", vars), with = FALSE])]
  xx <- model.matrix(reformulate(vars), data = z)
  rr <- qr.resid(qr(xx), cbind(rank(z$v), rank(z$CI_minus_CII))); r <- cor(rr[, 1], rr[, 2]); df <- nrow(xx) - ncol(xx) - 1
  c(r, 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df), nrow(z)) }
base <- c("batch", "age", "bmi_num", "sex")
tab <- rbindlist(lapply(rownames(GM), function(p) { a <- pc(GM[p, ], base); b <- pc(GM[p, ], c(base, "Ballooning"))
  data.table(gene = p, n = a[3], rho = a[1], p = a[2], rho_adj_ballooning = b[1], p_adj_ballooning = b[2]) }))
tab[, `:=`(bh = p.adjust(p, "BH"), bh_adj_ballooning = p.adjust(p_adj_ballooning, "BH"))]
fwrite(tab[order(p)], file.path(out, "partners_of_CI_minus_CII.tsv"), sep = "\t")
FIX <- list(
  immunoproteasome = c("PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2"),
  proteasome_core = c(paste0("PSMA", 1:7), paste0("PSMB", 1:7)),
  proteasome_19S = c(paste0("PSMC", 1:6), paste0("PSMD", 1:14)),
  mitophagy = c("PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1", "OPTN", "TAX1BP1", "CALCOCO2", "SQSTM1", "FKBP8", "PHB2"),
  lysosome = c("LAMP1", "LAMP2", "CTSB", "CTSD", "CTSZ", "CTSS", "HEXA", "HEXB", "GBA1", "NPC2", "ATP6V1A", "ATP6V1B2", "ATP6V0D1", "LIPA", "ASAH1"),
  macrophage = c("CD68", "CD163", "MSR1", "MRC1", "MARCO", "LYZ", "CD14", "FCGR3A", "C1QA", "C1QB", "C1QC", "VSIG4"),
  neutrophil = c("MPO", "ELANE", "S100A8", "S100A9", "LCN2", "LTF", "CAMP", "PRTN3", "AZU1"),
  interferon_ISG = c("ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2", "OAS3", "STAT1", "IFI35", "GBP1", "IRF9", "BST2"),
  MHC_II = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74"),
  nitrosative = c("NOS2", "SOD2", "PRDX5", "TXNRD1", "GPX1", "CAT"),
  mt_proteases = c("LONP1", "CLPP", "CLPX", "YME1L1", "AFG3L2", "SPG7", "OMA1", "PARL", "HTRA2"),
  MDB_keratins = c("KRT8", "KRT18", "SQSTM1", "UBB", "HSPB1", "DNAJB1"))
hm <- fread(file.path(ROOT, "data/reference_signatures/hallmark_v2025_1.tsv"))
SETS <- c(FIX, split(hm$gene, hm$program))
st <- rbindlist(lapply(names(SETS), function(s) { mem <- intersect(SETS[[s]], tab$gene); if (length(mem) < 5) return(NULL)
  rbindlist(lapply(c("rho", "rho_adj_ballooning"), function(col) { a <- tab[[col]][tab$gene %in% mem]; o <- tab[[col]][!tab$gene %in% mem]
    data.table(set = s, kind = if (s %in% names(FIX)) "fixed" else "hallmark", adjust = col, n = length(mem), mean_rho = mean(a),
               mean_rho_other = mean(o), p = wilcox.test(a, o)$p.value, members = paste(mem, collapse = ",")) })) }))
st[, bh := p.adjust(p, "BH"), by = adjust]
fwrite(st[order(adjust, p)], file.path(out, "set_tests.tsv"), sep = "\t")
fwrite(data.table(n_masld = nrow(M), n_partner_proteins = nrow(tab), n_excluded = length(excl)), file.path(out, "design.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(tab[order(p)][1:40], digits = 2); print(st[order(adjust, p)][, .(set, adjust, n, mean_rho, p, bh)][, head(.SD, 25), by = adjust], digits = 2)
