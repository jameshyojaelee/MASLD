suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma)})
method <- commandArgs(TRUE)[1]
SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
d <- readRDS(file.path(SS, "out_s2_qc/dge_included.rds")); cnt <- d$counts
st <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv", colClasses="character")
st[] <- lapply(st, function(x) gsub("\r","",x)); setnames(st, gsub("\r","",names(st)))
stg <- suppressWarnings(as.integer(st$nash_crn_kleiner_fibrosis_stage[match(colnames(cnt), st$sample_alias)]))
k <- stg %in% 1:4 & d$samples$inferred_sex %in% c("F","M")
cs <- cnt[, k]; s <- stg[k]; sm <- d$samples[k, ]
b <- factor(ifelse(s == 1, "F0F1", ifelse(s == 2, "F2", "F3F4")))
X <- model.matrix(~ factor(sm$batch) + factor(sm$inferred_sex) + s)
y <- DGEList(cs, lib.size = colSums(cs)); y <- calcNormFactors(y, method = method)
kf <- filterByExpr(y, group = b); y <- y[kf, , keep.lib.sizes = TRUE]
cat(method, "norm-then-filter genes", nrow(y), "\n")
for (renorm in c(FALSE, TRUE)) {
  yy <- y; if (renorm) yy <- calcNormFactors(yy, method = method)
  v <- voomWithQualityWeights(yy, X); fit <- eBayes(lmFit(v, X))
  lt <- data.table(gene = rownames(yy), t = fit$t[, "s"])
  lo <- fread(file.path(SS, "out_s3fgh-20261026T233706Z/s3g/limma_stage_num_t.tsv.gz"))
  mm <- merge(lt, lo, by = "gene")
  cat(method, "renorm", renorm, "n", nrow(lt), "overlap", nrow(mm), "max abs t diff", max(abs(mm$t.x - mm$t.y)), "\n")
  fwrite(lt, sprintf("%s/verify_s4/limma_t_%s_renorm%s.tsv", SS, method, renorm), sep = "\t")
}
