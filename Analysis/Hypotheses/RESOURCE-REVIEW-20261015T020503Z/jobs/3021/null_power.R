# Gate L02-02 (job 3021): permutation-null spread and per-gene detectable effects
# for promoter H3K27ac in GSE267119 (99 participants). Feasibility only.
# The observed disease-vs-NOR correlation is NOT computed here.
suppressMessages({library(edgeR); library(data.table)})
set.seed(20261015)
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out  <- file.path(root, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3021")
h3   <- file.path(root, "Analysis/MASLD_Model_Benchmark/executions/gse267145-matrix-audit-21064333/raw/GSE267119_H3K27ac.txt.gz")
join <- fread(file.path(root, "Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"))
gm   <- fread(file.path(root, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L02-02/gene_promoter_map.tsv"))

hdr <- scan(h3, what = "", nlines = 1, quiet = TRUE)
dt  <- fread(cmd = paste("zcat", h3), skip = 1, header = FALSE)
keys <- gsub('"', "", dt[[1]])
cnt <- as.matrix(dt[, -1, with = FALSE]); rownames(cnt) <- keys
pid <- sub("_.*", "", hdr); colnames(cnt) <- pid
stopifnot(ncol(cnt) == 99, nrow(cnt) == 96460, all(pid %in% join$participant_id))
st  <- join$stage[match(pid, join$participant_id)]
sex <- join$sex[match(pid, join$participant_id)]
dis <- st != "NOR"
cat("n disease", sum(dis), "n NOR", sum(!dis), "\n")

y <- DGEList(cnt); y <- calcNormFactors(y, method = "TMM")
lc <- cpm(y, log = TRUE, prior.count = 2)

hit <- gm[canonical == TRUE & n_reg_gene_tss > 0]
cat("hit genes", nrow(hit), "distinct regions", length(unique(hit$best_region)), "\n")
ri <- match(hit$best_region, rownames(lc)); stopifnot(!anyNA(ri))
X <- lc[unique(ri), , drop = FALSE]; gi <- match(ri, unique(ri))

# ---- 1. label-permutation null of the global Spearman (proxy: difference of mean logCPM)
B <- 2000
w_of <- function(d) ifelse(d, 1 / sum(d), -1 / sum(!d))
null_rho <- numeric(B)
for (b in seq_len(B)) {
  d <- sample(dis)
  diff <- as.vector(X %*% w_of(d))
  null_rho[b] <- cor(hit$logFC, diff[gi], method = "spearman")
}
q <- quantile(null_rho, c(0.025, 0.5, 0.975))
mde_rho <- (qnorm(0.975) + qnorm(0.8)) * sd(null_rho)
cat(sprintf("NULL global Spearman (B=%d, seed 20261015): mean %.4f sd %.4f q2.5 %.4f q50 %.4f q97.5 %.4f; MDE(alpha .05 two-sided, power .8) ~ %.3f\n",
            B, mean(null_rho), sd(null_rho), q[1], q[2], q[3], mde_rho))
cat(sprintf("Permutation p floor with 10,000 draws: %.6f\n", 1 / 10001))

# ---- 2. per-gene detectable log2FC from pooled within-group SD of logCPM
psd <- function(M, g) {
  v1 <- apply(M[, g, drop = FALSE], 1, var); v0 <- apply(M[, !g, drop = FALSE], 1, var)
  sqrt(((sum(g) - 1) * v1 + (sum(!g) - 1) * v0) / (length(g) - 2))
}
sdw <- psd(X, dis)[gi]
z_nom <- qnorm(0.975) + qnorm(0.8)
z_bon <- qnorm(1 - 0.025 / nrow(hit)) + qnorm(0.8)
k <- sqrt(1 / sum(dis) + 1 / sum(!dis))
mde_nom <- z_nom * sdw * k; mde_bon <- z_bon * sdw * k
cat("within-group SD of promoter logCPM, quantiles 10/25/50/75/90:", round(quantile(sdw, c(.1,.25,.5,.75,.9)), 3), "\n")
cat("per-gene MDE log2FC, disease(75) vs NOR(24), alpha .05 nominal: median", round(median(mde_nom), 3),
    "IQR", round(quantile(mde_nom, c(.25,.75)), 3), "\n")
cat("per-gene MDE log2FC, Bonferroni bound over", nrow(hit), "(alpha", signif(0.05/nrow(hit),3), "): median", round(median(mde_bon), 3),
    "IQR", round(quantile(mde_bon, c(.25,.75)), 3), "\n")
f23 <- st %in% c("NASH_F23", "NOR")
sd23 <- psd(X[, f23], st[f23] == "NASH_F23")[gi]
cat("for reference, NASH_F23(13) vs NOR(24) per-gene MDE, nominal: median",
    round(median(z_nom * sd23 * sqrt(1/13 + 1/24)), 3), "\n")

# ---- 3. positive control region (COL1A1, region within 1 kb of TSS)
pc <- hit[symbol == "COL1A1"]
pr <- pc$best_region
v  <- lc[pr, ]
cat("COL1A1 region", pr, "mean count", round(pc$best_region_mean_count, 1), "\n")
cat("COL1A1 within-group SD", round(sdw[hit$symbol == "COL1A1"], 3),
    "; MDE nominal disease vs NOR", round(mde_nom[hit$symbol == "COL1A1"], 3),
    "; MDE nominal F23 vs NOR", round(z_nom * sd23[hit$symbol == "COL1A1"] * sqrt(1/13 + 1/24), 3), "\n")
cat("COL1A1 region mean logCPM by stage:\n"); print(round(tapply(v, st, mean), 3))

# ---- 4. sex structure in the design
print(table(stage = st, sex = sex))

fwrite(data.table(null_rho = null_rho), file.path(out, "null_rho_proxy.tsv"), sep = "\t")
fwrite(data.table(gene = hit$gene, symbol = hit$symbol, bulk_logFC = hit$logFC, region = hit$best_region,
                  within_sd = sdw, mde_nominal = mde_nom, mde_bonferroni = mde_bon),
       file.path(out, "per_gene_mde.tsv"), sep = "\t")

# ---- 5. cost of ONE full edgeR QL fit on all 96,460 regions (what each of 10,000 draws would need)
design <- model.matrix(~ dis + sex)
tm <- system.time({
  keep <- filterByExpr(y, design)
  yk <- y[keep, , keep.lib.sizes = FALSE]; yk <- calcNormFactors(yk)
  fit <- glmQLFit(yk, design)
})
cat("filterByExpr keeps", sum(keep), "of", nrow(y), "regions; hit regions kept:",
    sum(unique(hit$best_region) %in% rownames(yk)), "of", length(unique(hit$best_region)), "\n")
cat(sprintf("one full glmQLFit elapsed %.1f s -> 10,000 draws ~ %.1f CPU-hours\n", tm[["elapsed"]], tm[["elapsed"]] * 10000 / 3600))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
