# DIG1 job 30 (exploratory; follows jobs 19, 22b, 25). Unbiased scan: which protein complexes lose
# protein with ballooning (PXD051911, within MASLD) without losing RNA (GSE130970, within disease),
# and which differ between MASLD protein and bulk RNA disease effects? OXPHOS and the mitoribosome
# are the known positives; the question is what else shows the same pattern.
# Complexes: GO cellular-component terms that are descendants of protein-containing complex
# (GO:0032991), org.Hs.eg.db / GO.db as installed, 5-200 members measured in the gene table.
# Gene tables (fixed upstream): out_bulk_rna_ballooning/per_gene_rna_vs_protein_ballooning.tsv.gz
# (rho_prot, rho_rna) and out_prot_rna_discord/gene_discordance_disease.tsv.gz (logFC_prot,
# logFC_rna, d = discordance z).
# Per complex: mean rho_prot, mean rho_rna, gap = mean(rho_prot - rho_rna); mean d. Null: 10,000
# random gene sets of the same size from the measured genes (seed 20261024); two-sided empirical p;
# BH over complexes per statistic. Genes in a complex are co-regulated, so this null is
# anti-conservative; the ranking, not the p, is the result. Terms overlap; the output keeps the
# member list so overlap can be read.
# Usage: Rscript 30_complex_scan.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(org.Hs.eg.db); library(GO.db); library(AnnotationDbi)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z"
bal <- fread(file.path(RUN, "out_bulk_rna_ballooning/per_gene_rna_vs_protein_ballooning.tsv.gz"))[, .(gene, rho_prot, rho_rna, gap = rho_prot - rho_rna)]
dis <- fread(file.path(RUN, "out_prot_rna_discord/gene_discordance_disease.tsv.gz"))[, .(gene, logFC_prot, logFC_rna, d)]
bal <- unique(bal[complete.cases(bal)], by = "gene"); dis <- unique(dis[complete.cases(dis)], by = "gene")
off <- unique(c("GO:0032991", get("GO:0032991", GOCCOFFSPRING)))
m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = off, keytype = "GOALL", columns = "SYMBOL"))
m <- unique(as.data.table(m)[ONTOLOGYALL == "CC" & !is.na(SYMBOL), .(go = GOALL, gene = SYMBOL)])
m[, term := Term(go)]
fwrite(data.table(n_terms = uniqueN(m$go), n_genes = uniqueN(m$gene), GO_db = as.character(packageVersion("GO.db")),
                  org_db = as.character(packageVersion("org.Hs.eg.db"))), file.path(out, "go_source.tsv"), sep = "\t")
scan <- function(tab, stats, label) {
  mm <- m[gene %in% tab$gene]; sz <- mm[, .N, by = go][N >= 5 & N <= 200]
  mm <- mm[go %in% sz$go]
  vals <- setNames(lapply(stats, function(s) setNames(tab[[s]], tab$gene)), stats)
  obs <- mm[, c(list(term = term[1], n = .N, members = paste(sort(gene), collapse = ",")),
               setNames(lapply(stats, function(s) mean(vals[[s]][gene])), paste0("mean_", stats))), by = go]
  sizes <- sort(unique(obs$n)); pool <- tab$gene
  nulls <- lapply(sizes, function(k) { idx <- replicate(10000, sample.int(length(pool), k))
    if (k == 1) idx <- matrix(idx, nrow = 1)
    sapply(stats, function(s) colMeans(matrix(vals[[s]][pool][idx], nrow = k))) })
  names(nulls) <- sizes
  for (s in stats) {
    mu <- mean(tab[[s]], na.rm = TRUE)
    obs[, (paste0("p_", s)) := mapply(function(v, k) { nd <- nulls[[as.character(k)]][, s]
      (1 + sum(abs(nd - mu) >= abs(v - mu) - 1e-12)) / (length(nd) + 1) }, get(paste0("mean_", s)), n)]
    obs[, (paste0("bh_", s)) := p.adjust(get(paste0("p_", s)), "BH")]
  }
  obs[, scan := label][] }
A <- scan(bal, c("rho_prot", "rho_rna", "gap"), "ballooning_within_disease")
B <- scan(dis, c("logFC_prot", "logFC_rna", "d"), "MASLD_vs_control")
fwrite(A[order(p_gap, mean_gap)], file.path(out, "complexes_ballooning_protein_vs_rna.tsv"), sep = "\t")
fwrite(B[order(p_d, mean_d)], file.path(out, "complexes_disease_protein_vs_rna.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
show <- function(x, cols) print(x[, ..cols], digits = 2)
show(A[order(mean_gap)][1:40], c("term", "n", "mean_rho_prot", "mean_rho_rna", "mean_gap", "bh_gap"))
show(A[order(-mean_gap)][1:15], c("term", "n", "mean_rho_prot", "mean_rho_rna", "mean_gap", "bh_gap"))
show(B[order(mean_d)][1:25], c("term", "n", "mean_logFC_prot", "mean_logFC_rna", "mean_d", "bh_d"))
show(B[order(-mean_d)][1:25], c("term", "n", "mean_logFC_prot", "mean_logFC_rna", "mean_d", "bh_d"))
