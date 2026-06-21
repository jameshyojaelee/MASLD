#!/usr/bin/env Rscript
# ============================================================================
# q1_metabolic_arm_dhps.R
#
# Q1: Vacca's published DHPS is TWO-ARMED — a metabolic arm (MOESM9 col5) and a
# fibrotic arm (MOESM9 col9). Script 10 reproduced the *mean/fibrotic-leaning*
# DHPS with a PEP-correlation against the human Severe-vs-Mild pathway reference.
# Here we test ALTERNATIVE constructions purpose-built to reproduce the METABOLIC
# arm specifically, then apply the winner to our 4 diets.
#
# Candidates (each = a per-model scalar, correlated across Vacca's 37 mouse models
# vs the published METABOLIC-arm DHPS, MOESM9 col5):
#   q1a  PEP-corr vs human "Mild vs Control" pathway ref (early/metabolic contrast)
#        + also vs Moderate-vs-Mild and Severe-vs-Mild as comparators.
#   q1b  PEP-corr restricted to METABOLIC KEGG pathways (lipid/glucose/insulin/
#        PPAR/fatty-acid/cholesterol/...), vs each human contrast.
#   q1c  GENE-arm DSEA-analog using ONLY metabolic-pathway genes (KEGG metabolic
#        leading-edge union), Vacca 951 sig restricted to those genes, ranked by
#        each model's per-gene L2FC (MOESM4).
# Baseline carried for reference: full-pathway PEP vs Severe-vs-Mild (= script 10).
#
# WINNER = candidate with max |Spearman| vs metabolic-arm DHPS. Apply to our 4
# diets (our per-diet KEGG NES, our per-diet DE) and report.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl); library(fgsea); library(msigdbr) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(CS, "vacca_benchmark")
ANNOT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation")
MPD   <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
strip <- function(x) gsub("\\..*", "", x)
nrm   <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
pnrm  <- function(x) toupper(trimws(gsub("_", " ", sub("^KEGG_", "", x))))
sp    <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
pe    <- function(a, b) suppressWarnings(cor(a, b, method = "pearson",  use = "complete.obs"))

# Metabolic pathway keyword filter (applied to UPPERCASE space-separated names)
METAB_KW <- paste0("LIPID|GLUCOS|INSULIN|PPAR|FATTY ACID|FATTY.?ACID|GLYCOLYSIS|",
  "GLUCONEO|CHOLESTEROL|BILE|STEROID|UNSATURATED|LINOLEIC|KETONE|",
  "CARBON METABOLISM|RETINOL|ARACHIDONIC|PROPANOATE|BUTANOATE|PYRUVATE|",
  "GLYCINE SERINE|FAT DIGESTION|ADIPO|TYPE II DIABETES")

# ── Published DHPS (MOESM9): col5 = metabolic arm, col9 = fibrotic arm ─────────
m9 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM9_ESM.xlsx"),
                               col_names = FALSE, skip = 2))
setnames(m9, 1:10, c("model","diet_group","mP","mH","m_DHPS","mM","fP","fH","f_DHPS","fM"))
m9 <- m9[!is.na(model)]
for (cc in c("m_DHPS","f_DHPS")) m9[[cc]] <- as.numeric(m9[[cc]])
m9[, key := nrm(model)]

# ── MOESM6 NES: human pathway reference + per-model PEPs ──────────────────────
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"), sheet = "NES"))
setnames(p6, 1, "pathway"); p6[, pname := toupper(trimws(pathway))]
human_cols <- grep("UCAM|EPoS", names(p6), value = TRUE)
model_cols <- setdiff(names(p6)[-1], c(human_cols, "pname"))
for (cc in c(human_cols, model_cols)) suppressWarnings(p6[[cc]] <- as.numeric(p6[[cc]]))

# human reference contrasts (in pathway space)
href <- list(
  MildCtrl = p6[["UCAM/VCU: Mild vs Control"]],
  ModMild  = rowMeans(cbind(p6[["UCAM/VCU: Moderate vs Mild"]], p6[["EPoS: Moderate vs Mild"]]), na.rm = TRUE),
  SevMild  = rowMeans(cbind(p6[["UCAM/VCU: Severe vs Mild"]],   p6[["EPoS: Severe vs Mild"]]),   na.rm = TRUE))

# metabolic-pathway mask over Vacca's 56 pathways
metab_mask <- grepl(METAB_KW, p6$pname)
cat(sprintf("Vacca pathway set: %d total ; %d metabolic\n", nrow(p6), sum(metab_mask)))
cat("  metabolic pathways:\n"); for (z in p6$pname[metab_mask]) cat("    -", z, "\n")

# ── q1a + q1b : per-model PEP correlations (full set and metabolic-only) ──────
recon <- rbindlist(lapply(model_cols, function(mc) {
  mv <- p6[[mc]]
  r <- data.table(model = mc); r[, key := nrm(mc)]
  for (h in names(href)) {
    r[[paste0("pepFull_", h)]]  <- pe(mv,             href[[h]])             # full-set PEP (q1a comparators)
    r[[paste0("pepMetab_", h)]] <- pe(mv[metab_mask], href[[h]][metab_mask]) # metabolic-only PEP (q1b)
  }
  r
}))

# ── q1c : metabolic-gene DSEA-analog from Vacca per-gene L2FC (MOESM4) ────────
# Build metabolic gene universe = leading genes of KEGG metabolic pathways (msigdbr),
# intersect Vacca 951 signature with them, then DSEA-analog on each model's L2FC.
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"), sheet = "Table S4"))
setnames(s4, c("Gene_used_in_DSEA(0:No/1:Yes)", "ENSGENEID_MOUSE"),
         c("in_dsea", "mouse_id"), skip_absent = TRUE)
s4[, mouse_id := strip(mouse_id)]
l2fc_cols  <- grep("^L2FC_", names(s4), value = TRUE)
hcol4      <- grep("UCAM|EPoS", l2fc_cols, value = TRUE)
mcol4      <- setdiff(l2fc_cols, hcol4)
for (cc in l2fc_cols) suppressWarnings(s4[, (cc) := as.numeric(get(cc))])
prog_cols <- c("L2FC_UCAM/VCU: Severe vs Mild", "L2FC_EPoS: Severe vs Mild")
s4[, vacca_prog_lfc := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = prog_cols]

# KEGG metabolic gene universe (mouse ENSEMBL) via msigdbr (>=10.0: collection/subcollection;
# KEGG legacy set renamed CP:KEGG_LEGACY)
kegg <- as.data.table(msigdbr(species = "Mus musculus", collection = "C2",
                              subcollection = "CP:KEGG_LEGACY"))
ens_col <- intersect(c("ensembl_gene","db_ensembl_gene","ensembl_gene_id"), names(kegg))[1]
kegg[, pname := pnrm(gs_name)]
metab_kegg_genes <- unique(kegg[grepl(METAB_KW, pname)][[ens_col]])
cat(sprintf("\nKEGG metabolic gene universe (mouse ENSEMBL): %d genes from %d metabolic KEGG sets\n",
            length(metab_kegg_genes), length(unique(kegg[grepl(METAB_KW, pname)]$pname))))

vsig_metab <- s4[in_dsea == 1 & !is.na(mouse_id) & mouse_id != "" & is.finite(vacca_prog_lfc) &
                 mouse_id %in% metab_kegg_genes]
metab_sets <- list(up   = unique(vsig_metab[vacca_prog_lfc > 0]$mouse_id),
                   down = unique(vsig_metab[vacca_prog_lfc < 0]$mouse_id))
cat(sprintf("Vacca metabolic-gene signature: %d up / %d down\n",
            length(metab_sets$up), length(metab_sets$down)))

dsea_analog <- function(stat, sets) {
  stat <- stat[is.finite(stat)]; stat <- stat[!duplicated(names(stat))]
  if (length(stat) < 200 || length(sets$up) < 5 || length(sets$down) < 5) return(NA_real_)
  fg <- as.data.table(suppressWarnings(fgsea(sets, stat, scoreType = "std", nPermSimple = 5000)))
  nu <- fg[pathway == "up", NES]; nd <- fg[pathway == "down", NES]
  if (length(nu) == 0) nu <- 0; if (length(nd) == 0) nd <- 0
  (nu - nd) / 2
}
geneArm <- rbindlist(lapply(mcol4, function(col) {
  v <- s4[is.finite(get(col)) & mouse_id != "", .(mouse_id, lfc = get(col))][!duplicated(mouse_id)]
  st <- setNames(v$lfc, v$mouse_id)
  rr <- data.table(metabGeneArm = dsea_analog(st, metab_sets))
  rr[, key := nrm(sub("^L2FC_", "", col))][]
}))
recon <- merge(recon, geneArm, by = "key", all.x = TRUE)

# ── Score every candidate vs the METABOLIC-arm DHPS (col5) across mouse models ─
mg <- merge(recon, m9[, .(key, diet_group, m_DHPS, f_DHPS)], by = "key")
mg <- mg[!grepl("^R[-.]", model)]   # drop rats → mouse models only
cat(sprintf("\nMouse models matched to published DHPS: %d\n", nrow(mg)))

cand <- setdiff(names(mg), c("key","model","diet_group","m_DHPS","f_DHPS"))
score <- rbindlist(lapply(cand, function(cc) data.table(
  candidate      = cc,
  spearman_metab = sp(mg[[cc]], mg$m_DHPS),    # vs METABOLIC arm (target)
  spearman_fibro = sp(mg[[cc]], mg$f_DHPS),    # vs fibrotic arm (contrast)
  n              = sum(is.finite(mg[[cc]])))))
score[, abs_metab := abs(spearman_metab)]
setorder(score, -abs_metab)
cat("\n=== CANDIDATE vs published METABOLIC-arm DHPS (Spearman, mouse models) ===\n")
print(score[, .(candidate, spearman_metab = round(spearman_metab,3),
                spearman_fibro = round(spearman_fibro,3), n)])
fwrite(score, file.path(OUT, "q1_metabolic_arm_candidate_scores.csv"))

winner <- score[which.max(abs_metab)]
cat(sprintf("\n>> WINNER (best metabolic-arm reproduction): %s  (Spearman = %.3f vs metab DHPS; %.3f vs fibro)\n",
            winner$candidate, winner$spearman_metab, winner$spearman_fibro))
# script-10 baseline for context
base_row <- score[candidate == "pepFull_SevMild"]
cat(sprintf("   baseline (full-set PEP vs Severe-vs-Mild, ~script10): Spearman %.3f vs metab DHPS\n",
            base_row$spearman_metab))

# ════════════════════════════════════════════════════════════════════════════
# APPLY THE WINNER TO OUR 4 DIETS
# ════════════════════════════════════════════════════════════════════════════
cat("\n=== Apply winner to OUR 4 diets ===\n")

# our per-diet KEGG NES (for PEP candidates)
mr <- fread(file.path(CS, "fgsea_mouse_results.csv"))[grepl("^KEGG_", pathway)]
mr[, pname := pnrm(pathway)]
diets <- c("MCD","HFD","CDAHFD","FPC")

our <- NULL
win <- winner$candidate

if (grepl("^pep", win)) {
  # which human contrast + whether metabolic-restricted
  hkey <- sub(".*_", "", win)                       # MildCtrl / ModMild / SevMild
  is_metab <- grepl("^pepMetab", win)
  href_dt <- data.table(pname = p6$pname, h_nes = href[[hkey]], is_metab = metab_mask)
  pool <- if (is_metab) href_dt[is_metab == TRUE] else href_dt
  shared <- intersect(unique(mr$pname), pool[is.finite(h_nes)]$pname)
  cat(sprintf("Winner is PEP (%s, metab_restricted=%s). Shared pathways ours∩Vacca: %d\n",
              hkey, is_metab, length(shared)))
  href_v <- pool[match(shared, pname)]$h_nes
  our <- rbindlist(lapply(diets, function(d) {
    v <- mr[source == d][match(shared, pname)]
    data.table(our_diet = d, our_score = pe(v$NES, href_v),
               our_score_spearman = sp(v$NES, href_v), n_shared = length(shared))
  }))
} else {
  # gene-arm winner: metabolic-gene DSEA-analog on our per-diet DE
  cat("Winner is the metabolic GENE-arm. Applying metab_sets to our per-diet DE t-stats.\n")
  our <- rbindlist(lapply(diets, function(d) {
    md <- fread(file.path(MPD, paste0(d, "_de_results.csv")))
    md[, mbase := strip(gene)]; md <- md[is.finite(t)]
    st <- setNames(md$t, md$mbase)
    data.table(our_diet = d, our_score = dsea_analog(st, metab_sets), n_shared = NA_integer_)
  }))
}

# attach Vacca same-diet metabolic-arm DHPS + same-diet winner score
grp <- list(MCD = "CDD", CDAHFD = "CDHFD", HFD = "HFD")
sel_of <- function(d) if (d == "FPC") mg[grepl("^WD", diet_group)] else mg[diet_group == grp[[d]]]
our[, vacca_metab_DHPS  := sapply(our_diet, function(d) mean(sel_of(d)$m_DHPS, na.rm = TRUE))]
our[, vacca_winnerScore := sapply(our_diet, function(d) mean(sel_of(d)[[win]], na.rm = TRUE))]
setorder(our, -our_score)
cat("\n  OUR 4 diets under the winning metabolic-arm construction:\n")
print(our)
cat(sprintf("\n  Spearman across 4 diets: our_score vs Vacca same-diet metabolic DHPS = %.2f\n",
            sp(our$our_score, our$vacca_metab_DHPS)))

# rank our diets among Vacca models on the winner scale
allm <- rbind(mg[, .(label = model, group = "Vacca", s = get(win))],
              our[, .(label = paste0("OUR_", our_diet), group = "Ours", s = our_score)])
setorder(allm, -s); allm[, rank := .I]
cat("\n  Our diets' rank among all", nrow(allm), "models on the winner score:\n")
print(allm[group == "Ours", .(rank, label, score = round(s,3))])

fwrite(our, file.path(OUT, "q1_metabolic_arm_our_diets.csv"))
fwrite(allm, file.path(OUT, "q1_metabolic_arm_ranking.csv"))
cat("\n", strrep("=", 78), "\n")
