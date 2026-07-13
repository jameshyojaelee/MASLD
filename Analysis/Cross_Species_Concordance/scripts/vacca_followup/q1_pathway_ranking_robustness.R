#!/usr/bin/env Rscript
# ============================================================================
# q1_pathway_ranking_robustness.R
#
# Q1 ROBUSTNESS: pathway-PEP ranking of our 4 mouse diets vs the Vacca human
# pathway reference. Script 10 fixed ONE setting (Pearson PEP; human ref =
# Severe-vs-Mild averaged over UCAM/VCU+EPoS; all shared KEGG pathways) and got
# FPC 0.924 (top) > MCD 0.690 > HFD 0.642 > CDAHFD 0.640.
#
# Here we sweep THREE knobs and ask whether two qualitative claims survive:
#   (A) FPC = top-ranked diet
#   (B) MCD > CDAHFD ordering
# Knobs:
#   (i)   correlation method:  Pearson vs Spearman
#   (ii)  human-reference contrast: Mild-vs-Control / Moderate-vs-Mild /
#         Severe-vs-Mild   (each = rowMeans over the UCAM/VCU + EPoS columns
#         that exist for that contrast)
#   (iii) pathway set: ALL shared KEGG  vs  LIVER-RELEVANT subset
#         (metabolic + fibrotic/ECM + liver-injury core; immune/cancer/generic
#          signalling dropped).
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(CS, "vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

pnrm <- function(x) toupper(trimws(gsub("_", " ", sub("^KEGG_", "", x))))  # KEGG -> MOESM6 space
cc   <- function(a, b, m) suppressWarnings(cor(a, b, method = m, use = "complete.obs"))

# ── Vacca human pathway reference (MOESM6 'NES') ──────────────────────────────
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"), sheet = "NES"))
setnames(p6, 1, "pathway"); p6[, pname := toupper(trimws(pathway))]

# Build the 3 human-reference contrasts; each averages whatever UCAM/VCU + EPoS
# columns exist for that contrast (Mild-vs-Control has NO EPoS column).
ref_cols <- list(
  MildVsControl  = c("UCAM/VCU: Mild vs Control"),
  ModerateVsMild = c("UCAM/VCU: Moderate vs Mild", "EPoS: Moderate vs Mild"),
  SevereVsMild   = c("UCAM/VCU: Severe vs Mild",   "EPoS: Severe vs Mild")
)
make_ref <- function(cols) {
  M <- sapply(cols, function(cn) suppressWarnings(as.numeric(p6[[cn]])))
  rowMeans(matrix(M, nrow = nrow(p6)), na.rm = TRUE)
}
human_refs <- lapply(ref_cols, function(cols)
  data.table(pname = p6$pname, h_nes = make_ref(cols)))

# ── OUR per-diet KEGG pathway NES (fgsea_mouse_results.csv) ───────────────────
mr <- fread(file.path(CS, "fgsea_mouse_results.csv"))[grepl("^KEGG_", pathway)]
mr[, pname := pnrm(pathway)]
diets <- c("MCD", "HFD", "CDAHFD", "FPC")

# Shared pathway universe (ours ∩ Vacca's 56 DSEA pathways)
shared_all <- intersect(unique(mr$pname),
                        p6[is.finite(make_ref(ref_cols$SevereVsMild))]$pname)
shared_all <- intersect(unique(mr$pname), p6$pname)   # full name-match set

# ── LIVER-RELEVANT subset: metabolic + fibrotic/ECM + liver-injury core ───────
# Keyword-defined over the Vacca pathway names; excludes generic immune/cancer/
# neuro signalling that is not liver-disease-axis-specific.
liver_kw <- c("LIPID","FATTY ACID","CHOLESTEROL","KETONE","GLYCOL","GLUCONEO",
              "PPAR","DIABETES","INSULIN","BILE","RETINOL","ARACHIDONIC",
              "LINOLEIC","STEROID","BIOSYNTHESIS OF UNSATURATED","DRUG METAB",
              "METABOLISM OF XENOBIOTICS","PEROXISOME","ADIPO",
              # fibrosis / ECM / stellate axis
              "ECM","FOCAL ADHESION","TGF","ACTIN","COMPLEMENT AND COAGULATION",
              "HEDGEHOG","WNT","HIPPO","AGE RAGE",
              # liver injury / hepatocyte stress
              "P53","APOPTOSIS","HEPATITIS","PRIMARY BILE")
is_liver <- function(pn) Reduce(`|`, lapply(liver_kw, function(k) grepl(k, pn)))
liver_set <- shared_all[is_liver(shared_all)]

cat(sprintf("Name-matched shared KEGG pathways (ours ∩ Vacca 56): %d\n", length(shared_all)))
cat(sprintf("Liver-relevant subset of those: %d\n", length(liver_set)))
cat("Liver-relevant pathways:\n"); print(sort(liver_set))
cat("\nDropped (non-liver) pathways:\n"); print(sort(setdiff(shared_all, liver_set)))

# ── Sweep: for every (method, ref-contrast, pathway-set), rank the 4 diets ────
methods   <- c(Pearson = "pearson", Spearman = "spearman")
pathsets  <- list(AllShared = shared_all, LiverRelevant = liver_set)

rows <- list()
for (rn in names(human_refs)) {
  href <- human_refs[[rn]]
  for (psn in names(pathsets)) {
    ps <- pathsets[[psn]]
    hv <- href[match(ps, pname)]$h_nes
    for (mn in names(methods)) {
      pep <- sapply(diets, function(d) {
        v <- mr[source == d][match(ps, pname)]$NES
        cc(v, hv, methods[[mn]])
      })
      ord  <- diets[order(-pep)]
      rk   <- rank(-pep)
      names(rk) <- diets
      rows[[length(rows) + 1]] <- data.table(
        ref_contrast = rn, pathway_set = psn, method = mn, n_paths = length(ps),
        pep_MCD = pep["MCD"], pep_HFD = pep["HFD"],
        pep_CDAHFD = pep["CDAHFD"], pep_FPC = pep["FPC"],
        top_diet = ord[1], rank_order = paste(ord, collapse = ">"),
        FPC_is_top = ord[1] == "FPC",
        MCD_gt_CDAHFD = pep["MCD"] > pep["CDAHFD"],
        rank_FPC = rk["FPC"], rank_MCD = rk["MCD"], rank_CDAHFD = rk["CDAHFD"])
    }
  }
}
res <- rbindlist(rows)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("FULL SWEEP (", nrow(res), " settings = 3 ref-contrasts x 2 path-sets x 2 methods)\n", sep = "")
cat(strrep("=", 78), "\n", sep = "")
print(res[, .(ref_contrast, pathway_set, method, n_paths,
              pep_FPC = round(pep_FPC,2), pep_MCD = round(pep_MCD,2),
              pep_CDAHFD = round(pep_CDAHFD,2), pep_HFD = round(pep_HFD,2),
              top_diet, FPC_top = FPC_is_top, MCD_gt_CDAHFD)])

cat("\n--- STABILITY SUMMARY (across", nrow(res), "settings) ---\n")
cat(sprintf("  FPC = top diet:        %d / %d (%.0f%%)\n",
            sum(res$FPC_is_top), nrow(res), 100*mean(res$FPC_is_top)))
cat(sprintf("  MCD > CDAHFD:          %d / %d (%.0f%%)\n",
            sum(res$MCD_gt_CDAHFD), nrow(res), 100*mean(res$MCD_gt_CDAHFD)))
cat(sprintf("  FPC mean rank: %.2f | MCD mean rank: %.2f | CDAHFD mean rank: %.2f | HFD mean rank: %.2f\n",
            mean(res$rank_FPC), mean(res$rank_MCD), mean(res$rank_CDAHFD),
            mean(rank_helper <- sapply(seq_len(nrow(res)), function(i){
              pep <- c(res$pep_MCD[i],res$pep_HFD[i],res$pep_CDAHFD[i],res$pep_FPC[i])
              rank(-pep)[2]}))))
cat("\n  Distinct rank orders observed:\n")
print(res[, .N, by = rank_order][order(-N)])

# Per-knob marginal stability
cat("\n--- MARGINAL: FPC-top & MCD>CDAHFD by each knob ---\n")
for (knob in c("method","ref_contrast","pathway_set")) {
  cat(sprintf("  by %s:\n", knob))
  agg <- res[, .(FPC_top = sprintf("%d/%d", sum(FPC_is_top), .N),
                 MCD_gt_CDAHFD = sprintf("%d/%d", sum(MCD_gt_CDAHFD), .N)),
             by = knob]
  print(agg)
}

fwrite(res, file.path(OUT, "q1_pathway_ranking_robustness.csv"))
cat("\nWrote:", file.path(OUT, "q1_pathway_ranking_robustness.csv"), "\n")
