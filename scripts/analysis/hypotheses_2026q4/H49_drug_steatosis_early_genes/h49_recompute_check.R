# Independent recomputation of H49 headline numbers from output tables and one per-compound RDS.
args <- commandArgs(TRUE); sd <- args[1]; pcd <- args[2]
rd <- function(f) read.delim(file.path(sd, f), check.names = FALSE)
pt <- rd("primary_tests.tsv")
q <- p.adjust(pt$p_two, "BH")
cat("max |q_BH - recomputed|:", max(abs(q - pt$q_BH)), "\n")
conc <- ifelse(pt$direction == "up", pt$S - pt$null_median, pt$null_median - pt$S)
vq <- rd("vehicle_summary.tsv")
vq95 <- setNames(vq$q95_abs_dev, vq$direction)[pt$direction]
vfail <- !setNames(vq$pass, vq$direction)[pt$direction]
pass_v <- ifelse(vfail, abs(pt$S - pt$null_median) > vq95, TRUE)
sigC <- conc > 0 & q < 0.05 & pass_v; sigA <- conc < 0 & q < 0.05 & pass_v
cat("concordant sig:", paste(pt$drug[sigC], pt$direction[sigC]), "\n")
cat("anticoncordant sig:", paste(pt$drug[sigA], pt$direction[sigA]), "\n")
nt <- rd("negative_control_drugs.tsv")
ncon <- ifelse(nt$direction == "up", nt$S - nt$null_median, nt$null_median - nt$S)
negmax <- tapply(ncon, nt$direction, max)
spec <- sigC & conc > negmax[pt$direction]
cat("specific:", paste(pt$drug[spec], pt$direction[spec]), "\n")
cat("F q<0.05 concordant:", paste(pt$drug[p.adjust(pt$p_F, "BH") < 0.05], pt$direction[p.adjust(pt$p_F, "BH") < 0.05]), "\n")
# vehicle false-positive fraction
vv <- rd("vehicle_vs_vehicle.tsv")
cat("vehicle n:", nrow(vv), " frac p<0.05 up:", mean(vv$p_up < 0.05), " down:", mean(vv$p_down < 0.05), "\n")
# S for amiodarone early-up recomputed from the RDS, probe selection and gene sets
ps <- rd("probe_selection.tsv"); gs <- rd("gene_sets_used.tsv")
r <- readRDS(file.path(pcd, "amiodarone.rds"))
up <- as.character(gs$entrez[gs$set == "early_up"]); dn <- as.character(gs$entrez[gs$set == "early_down"])
pid <- setNames(ps$ID, as.character(ps$Gene.ID))
z <- function(g) { t <- r$t[pid[g], "High_24h"]; d <- r$df[pid[g]]
  ifelse(t > 0, qnorm(pt(t, d, lower.tail = FALSE), lower.tail = FALSE), qnorm(pt(t, d))) }
cat("amiodarone High_24h S_up recomputed:", mean(z(up)), " table:", pt$S[pt$drug == "amiodarone" & pt$direction == "up"], "\n")
cat("amiodarone High_24h S_down recomputed:", mean(z(dn)), " table:", pt$S[pt$drug == "amiodarone" & pt$direction == "down"], "\n")
cat("amiodarone arrays:", r$n_arrays, "; groups:", paste(names(r$groups), r$groups, collapse = ", "), "\n")
# replicate counts for the primary condition across the 10 drugs
for (k in c("amiodarone", "tamoxifen", "perhexiline", "tetracycline", "nefazodone", "carbamazepine",
            "enalapril", "isoniazid", "ketoconazole", "valproicacid")) {
  g <- readRDS(file.path(pcd, paste0(k, ".rds")))$groups
  cat(k, "High_24h n =", g[["High_24h"]], "Control_24h n =", g[["Control_24h"]], "\n")
}
