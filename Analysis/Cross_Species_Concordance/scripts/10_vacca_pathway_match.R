#!/usr/bin/env Rscript
# ============================================================================
# 10_vacca_pathway_match.R
#
# Script 09 showed Vacca's DHPS is reproduced at the PATHWAY level (pathArm:
# human KEGG-pathway signature enriched in a model's pathway-NES profile;
# Spearman 0.73 vs published DHPS on their data), NOT the gene level (0.22).
#
# DECISIVE TEST: apply that VALIDATED pathway construction to OUR 4 mouse diets
# (our per-diet KEGG NES from fgsea_mouse_results.csv) using the SAME human
# pathway reference (Vacca MOESM6), and ask: do our diets now rank CONSISTENTLY
# with Vacca's same-diet models? If yes ⇒ our cross-species MATCHES LITMUS once
# computed at the right (pathway) resolution ⇒ MATCH branch is publishable.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl); library(fgsea) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(CS, "vacca_benchmark")
set.seed(42)
nrm  <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
pnrm <- function(x) toupper(trimws(gsub("_", " ", sub("^KEGG_", "", x))))   # KEGG name → MOESM6 space
sp   <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))

# PEP proximity = correlation of a model's pathway-NES vector with the human
# reference pathway-NES vector over a SHARED pathway set (the validated DHPS
# construction from 09: pepP reproduced DHPS at Spearman ~0.64).

# ── Vacca human pathway reference + per-model PEPs (MOESM6) ───────────────────
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"), sheet = "NES"))
setnames(p6, 1, "pathway"); p6[, pname := pnrm(pathway)]
href <- rowMeans(cbind(suppressWarnings(as.numeric(p6[["UCAM/VCU: Severe vs Mild"]])),
                       suppressWarnings(as.numeric(p6[["EPoS: Severe vs Mild"]]))), na.rm = TRUE)
human_ref <- data.table(pname = p6$pname, h_nes = href)
model_cols <- setdiff(names(p6)[-c(1, ncol(p6))], grep("UCAM|EPoS|pname", names(p6), value = TRUE))

# ── OUR per-diet KEGG pathway NES (reuse fgsea_mouse_results.csv) ─────────────
mr <- fread(file.path(CS, "fgsea_mouse_results.csv"))[grepl("^KEGG_", pathway)]
mr[, pname := pnrm(pathway)]
shared <- intersect(unique(mr$pname), human_ref[is.finite(h_nes)]$pname)
cat(sprintf("Shared KEGG pathways (ours ∩ Vacca human ref): %d.\n", length(shared)))
href_v <- human_ref[match(shared, pname)]$h_nes

# our diets: PEP correlation over shared pathways
diets <- c("MCD", "HFD", "CDAHFD", "FPC")
our <- rbindlist(lapply(diets, function(d) {
  v <- mr[source == d][match(shared, pname)]
  data.table(our_diet = d,
             our_pepP = suppressWarnings(cor(v$NES, href_v, method = "pearson",  use = "complete.obs")),
             our_pepS = suppressWarnings(cor(v$NES, href_v, method = "spearman", use = "complete.obs")))
}))

# ── Vacca per-model PEP over the SAME shared pathways + published DHPS ────────
for (cc in model_cols) suppressWarnings(p6[[cc]] <- as.numeric(p6[[cc]]))
vm_dhps <- fread(file.path(OUT, "dhps_reconstruction_per_model.csv"))[, .(model, diet_group, dhps_mean)]
vmod <- rbindlist(lapply(model_cols, function(mc) {
  v <- p6[match(shared, pname)][[mc]]
  r <- data.table(model = mc,
             pepP = suppressWarnings(cor(v, href_v, method = "pearson", use = "complete.obs")))
  r[, key := nrm(mc)][]
}))
vm_dhps[, key := nrm(model)]
vmod <- merge(vmod, vm_dhps[, .(key, diet_group, dhps_mean)], by = "key")[!grepl("^R[-.]", model)]
cat(sprintf("Vacca per-model PEP (over %d shared paths) vs published DHPS: Spearman %.2f (sanity ≥0.5 expected)\n",
            length(shared), sp(vmod$pepP, vmod$dhps_mean)))

# map our diet → Vacca diet group(s); attach same-diet DHPS + same-diet PEP
grp <- list(MCD = "CDD", CDAHFD = "CDHFD", HFD = "HFD", FPC = "^WD")
sel_of <- function(d) if (d == "FPC") vmod[grepl("^WD", diet_group)] else vmod[diet_group == grp[[d]]]
our[, vacca_DHPS := sapply(our_diet, function(d) mean(sel_of(d)$dhps_mean, na.rm = TRUE))]
our[, vacca_pepP := sapply(our_diet, function(d) mean(sel_of(d)$pepP,      na.rm = TRUE))]

cat("\n=== OUR 4 diets — pathway PEP proximity vs Vacca same-diet ===\n")
print(our[order(-vacca_DHPS)])
cat(sprintf("\n  Spearman across 4 diets:\n    our_pepP vs Vacca same-diet DHPS = %.2f\n    our_pepP vs Vacca same-diet PEP  = %.2f\n",
            sp(our$our_pepP, our$vacca_DHPS), sp(our$our_pepP, our$vacca_pepP)))
cat("  (for contrast — script 07 gene-level was Spearman -1.00)\n")

# place our diets among Vacca's models on the PEP scale
allm <- rbind(vmod[, .(label = model, group = "Vacca", pepP)],
              our[, .(label = paste0("OUR_", our_diet), group = "Ours", pepP = our_pepP)], fill = TRUE)
setorder(allm, -pepP); allm[, rank := .I]
cat("\n  Ranking of all models by pathway PEP (our diets flagged):\n")
print(allm[group == "Ours" | rank <= 5 | rank > .N - 3, .(rank, label, pepP = round(pepP, 2))])

fwrite(our, file.path(OUT, "our_diets_pathway_dhps.csv"))
fwrite(vmod[order(-pepP)], file.path(OUT, "vacca_models_pathway_pep.csv"))
cat("\n", strrep("=", 78), "\n")
cat("VERDICT: positive our_pepP-vs-Vacca-DHPS ⇒ our cross-species MATCHES LITMUS at\n")
cat("  pathway resolution ⇒ MATCH publishable. Negative/flat ⇒ our mouse DATA differs.\n")
cat(strrep("=", 78), "\n")
