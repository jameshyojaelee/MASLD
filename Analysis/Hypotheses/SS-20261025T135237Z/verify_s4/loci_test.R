suppressPackageStartupMessages(library(data.table))
SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
cg <- fread(file.path(SS, "frozen/dig1/coloc_genes_lead_z.tsv")); lc <- fread(file.path(SS, "frozen/dig1/loci.tsv"))
cg <- merge(unique(cg[, .(ensembl, locus)]), lc[, .(locus, cls)], by = "locus")[cls %in% c("release", "disease")]
lo <- fread(file.path(SS, "out_s3fgh-20261026T233706Z/s3g/limma_stage_num_t.tsv.gz"))
theirs_genes <- fread(file.path(SS, "out_s3fgh-20261026T233706Z/s3g/loci_gene_t.tsv"))
for (f in c("limma_t_TMM_renormTRUE.tsv", "THEIRS")) {
  lt <- if (f == "THEIRS") lo[, .(ensembl, t)] else fread(f)[, .(ensembl = sub("\\..*$", "", gene), t)]
  g <- unique(merge(cg, lt, by = "ensembl"))
  w <- wilcox.test(abs(g[cls == "release"]$t), abs(g[cls == "disease"]$t), conf.int = TRUE, exact = FALSE)
  m <- merge(g, theirs_genes[, .(ensembl, locus, t_o = t)], by = c("ensembl", "locus"))
  cat(f, "n", nrow(g), sum(g$cls == "release"), sum(g$cls == "disease"), "W", w$statistic, "p", format(w$p.value, digits = 15),
      "est", format(w$estimate, digits = 15), "max |t - their gene t|", max(abs(m$t - m$t_o)), "\n")
}
