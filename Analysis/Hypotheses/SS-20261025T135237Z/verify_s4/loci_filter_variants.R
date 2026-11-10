suppressPackageStartupMessages({library(data.table); library(edgeR)})
SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
d <- readRDS(file.path(SS, "out_s2_qc/dge_included.rds")); cnt <- d$counts
st <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv", colClasses="character")
st[] <- lapply(st, function(x) gsub("\r","",x)); setnames(st, gsub("\r","",names(st)))
stg <- suppressWarnings(as.integer(st$nash_crn_kleiner_fibrosis_stage[match(colnames(cnt), st$sample_alias)]))
lo <- fread(file.path(SS, "out_s3fgh-20261026T233706Z/s3g/limma_stage_num_t.tsv.gz"))
theirs <- lo$gene
cat("lib.size equal colSums:", all.equal(d$samples$lib.size, colSums(cnt)), "\n")
staged <- stg %in% 1:4
binf <- function(s) factor(ifelse(is.na(s) | !(s %in% 1:4), "U", ifelse(s == 1, "F0F1", ifelse(s == 2, "F2", "F3F4"))))
rep <- function(lab, k) cat(lab, sum(k), "same set as theirs:", setequal(rownames(cnt)[k], theirs), "\n")
cs <- cnt[, staged]; b <- binf(stg[staged])
rep("staged group=bin colSums lib", filterByExpr(DGEList(cs), group = b))
rep("staged group=bin dge lib", filterByExpr(DGEList(cs, lib.size = d$samples$lib.size[staged]), group = b))
rep("staged group=bin assigned lib", filterByExpr(DGEList(cs, lib.size = d$samples$assigned_fragments[staged]), group = b))
rep("all515 group=bin+U", filterByExpr(DGEList(cnt), group = binf(stg)))
rep("staged no group", filterByExpr(DGEList(cs)))
X <- model.matrix(~ factor(d$samples$batch[staged]) + d$samples$inferred_sex[staged] + stg[staged])
rep("staged design", filterByExpr(DGEList(cs), design = X))
y <- DGEList(cs); y <- calcNormFactors(y, method="TMM")
rep("staged group=bin after TMM", filterByExpr(y, group = b))
