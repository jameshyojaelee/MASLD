#!/usr/bin/env Rscript
# ============================================================================
# q1_metab_fibro_2x2.R   (Vacca follow-up Q1)
#
# RESOLVE Q1: in our prior single-score pathway match (script 10) MCD ranked
# ABOVE CDAHFD/HFD (MCD 0.69 > CDAHFD/HFD 0.64) using the FIBROSIS (Severe-vs-Mild)
# human reference. But Vacca's DHPS is TWO-ARMED — a METABOLIC arm and a FIBROTIC
# arm with DIFFERENT human references. MCD is fibrosis-strong but metabolically
# WEAK; its Vacca penalty is the METABOLIC arm. We therefore split our diets'
# pathway proximity into the SAME two arms and ask:
#   (i) is MCD LOW-metabolic / HIGH-fibrotic (matching Vacca MCD metab~0.55/fibro~0.75)?
#   (ii) does ranking diets by the METABOLIC arm put MCD BELOW CDAHFD/HFD?
#
# References (q1a/q1b):
#   q1a METABOLIC ref = human pathway-NES early/metabolic axis
#       = mean(UCAM/VCU: Mild vs Control, UCAM/VCU: Moderate vs Mild, EPoS: Moderate vs Mild)
#   q1b FIBROTIC ref  = human pathway-NES fibrosis axis
#       = mean(UCAM/VCU: Severe vs Mild, EPoS: Severe vs Mild)
# Proximity = Pearson corr of a model's pathway-NES vector with the arm reference
# over the SHARED KEGG-pathway set (the validated PEP/pathArm construction, script 09/10).
#
# Validation: per-arm proximity on VACCA's OWN models must track the published
# per-arm DHPS (m_DHPS metabolic / f_DHPS fibrotic, MOESM9). Then apply to OUR 4 diets.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
CS    <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
OUT   <- file.path(CS, "vacca_benchmark")
set.seed(42)
pnrm <- function(x) toupper(trimws(gsub("_", " ", sub("^KEGG_", "", x))))  # KEGG -> MOESM6 space
nrm  <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
sp   <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
pr   <- function(a, b) suppressWarnings(cor(a, b, method = "pearson",  use = "complete.obs"))

# ── Vacca human pathway reference + per-model pathway PEPs (MOESM6) ────────────
p6 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM6_ESM.xlsx"), sheet = "NES"))
setnames(p6, 1, "pathway"); p6[, pname := pnrm(pathway)]
numify <- function(col) suppressWarnings(as.numeric(p6[[col]]))

# q1a METABOLIC reference (early/metabolic axis)
metab_ref <- rowMeans(cbind(numify("UCAM/VCU: Mild vs Control"),
                            numify("UCAM/VCU: Moderate vs Mild"),
                            numify("EPoS: Moderate vs Mild")), na.rm = TRUE)
# q1b FIBROTIC reference (Severe-vs-Mild axis)
fibro_ref <- rowMeans(cbind(numify("UCAM/VCU: Severe vs Mild"),
                            numify("EPoS: Severe vs Mild")), na.rm = TRUE)
human_ref <- data.table(pname = p6$pname, metab = metab_ref, fibro = fibro_ref)

human_cols <- grep("UCAM|EPoS", names(p6), value = TRUE)
model_cols <- setdiff(names(p6)[-1], c(human_cols, "pname"))
for (cc in model_cols) suppressWarnings(p6[[cc]] <- as.numeric(p6[[cc]]))

# ── OUR per-diet KEGG pathway NES (fgsea_mouse_results.csv) ────────────────────
mr <- fread(file.path(CS, "fgsea_mouse_results.csv"))[grepl("^KEGG_", pathway)]
mr[, pname := pnrm(pathway)]

shared <- Reduce(intersect, list(unique(mr$pname),
                                 human_ref[is.finite(metab) & is.finite(fibro)]$pname))
cat(sprintf("Shared KEGG pathways (ours ∩ Vacca human ref, both arms finite): %d\n", length(shared)))
hr_m <- human_ref[match(shared, pname)]$metab
hr_f <- human_ref[match(shared, pname)]$fibro
cat(sprintf("  arm-reference independence: Pearson(metab_ref, fibro_ref) over shared = %.2f\n\n",
            pr(hr_m, hr_f)))

# ── (A) VALIDATE the two-arm construction on VACCA's OWN models ────────────────
#   per-model metabolic-proximity vs published m_DHPS ; fibrotic-prox vs f_DHPS
m9 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM9_ESM.xlsx"),
                               col_names = FALSE, skip = 2))
setnames(m9, 1:10, c("model","diet_group","mP","mH","m_DHPS","mM","fP","fH","f_DHPS","fM"))
m9 <- m9[!is.na(model)]; for (cc in c("m_DHPS","f_DHPS")) m9[[cc]] <- as.numeric(m9[[cc]])
m9[, key := nrm(model)]

vmod <- rbindlist(lapply(model_cols, function(mc) {
  v <- p6[match(shared, pname)][[mc]]
  r <- data.table(model = mc,
             metab_prox = pr(v, hr_m),
             fibro_prox = pr(v, hr_f))
  r[, key := nrm(mc)][]
}))
vmod <- merge(vmod, m9[, .(key, diet_group, m_DHPS, f_DHPS)], by = "key")[!grepl("^R[-.]", model)]
cat(sprintf("=== (A) Two-arm validation on Vacca's own mouse models (n=%d) ===\n", nrow(vmod)))
cat(sprintf("  metab_prox  vs published m_DHPS : Spearman %.2f\n", sp(vmod$metab_prox, vmod$m_DHPS)))
cat(sprintf("  fibro_prox  vs published f_DHPS : Spearman %.2f\n", sp(vmod$fibro_prox, vmod$f_DHPS)))
cat(sprintf("  cross-check metab_prox vs f_DHPS=%.2f ; fibro_prox vs m_DHPS=%.2f (should be lower/off-diagonal)\n\n",
            sp(vmod$metab_prox, vmod$f_DHPS), sp(vmod$fibro_prox, vmod$m_DHPS)))
fwrite(vmod[order(-m_DHPS)], file.path(OUT, "q1_vacca_model_arms.csv"))

# Vacca same-diet published arm targets + same-diet reconstructed arm proximity
grp <- list(MCD = "CDD", CDAHFD = "CDHFD", HFD = "HFD", FPC = "^WD")
sel_of <- function(d) if (d == "FPC") vmod[grepl("^WD", diet_group)] else vmod[diet_group == grp[[d]]]

# ── (B) apply two-arm construction to OUR 4 diets ─────────────────────────────
diets <- c("MCD","HFD","CDAHFD","FPC")
our <- rbindlist(lapply(diets, function(d) {
  v <- mr[source == d][match(shared, pname)]$NES
  sd <- sel_of(d)
  data.table(diet = d,
             our_metab_prox = pr(v, hr_m),
             our_fibro_prox = pr(v, hr_f),
             vacca_metab_DHPS = mean(sd$m_DHPS, na.rm = TRUE),   # published metabolic arm
             vacca_fibro_DHPS = mean(sd$f_DHPS, na.rm = TRUE),   # published fibrotic arm
             vacca_metab_prox = mean(sd$metab_prox, na.rm = TRUE),
             vacca_fibro_prox = mean(sd$fibro_prox, na.rm = TRUE))
}))

# rankings: lower rank-number = stronger on that arm
our[, our_metab_rank := frank(-our_metab_prox)]
our[, our_fibro_rank := frank(-our_fibro_prox)]
our[, vacca_metab_rank := frank(-vacca_metab_DHPS)]
our[, vacca_fibro_rank := frank(-vacca_fibro_DHPS)]

cat("=== (B) OUR 4 diets — per-diet (metabolic-proximity, fibrotic-proximity) 2x2 ===\n")
print(our[order(our_metab_rank),
          .(diet,
            our_metab_prox = round(our_metab_prox, 3),
            our_fibro_prox = round(our_fibro_prox, 3),
            our_metab_rank, our_fibro_rank,
            vacca_metab_DHPS = round(vacca_metab_DHPS, 3),
            vacca_fibro_DHPS = round(vacca_fibro_DHPS, 3))])

cat("\n=== METABOLIC-ARM ranking of our 4 diets (strongest metabolic first) ===\n")
print(our[order(our_metab_rank), .(rank = our_metab_rank, diet,
            our_metab_prox = round(our_metab_prox, 3))])

cat("\n=== FIBROTIC-ARM ranking of our 4 diets (strongest fibrotic first) ===\n")
print(our[order(our_fibro_rank), .(rank = our_fibro_rank, diet,
            our_fibro_prox = round(our_fibro_prox, 3))])

# Q1 resolution checks
mcd_metab_rank <- our[diet == "MCD", our_metab_rank]
cdahfd_metab_rank <- our[diet == "CDAHFD", our_metab_rank]
hfd_metab_rank <- our[diet == "HFD", our_metab_rank]
mcd_below_both <- (mcd_metab_rank > cdahfd_metab_rank) & (mcd_metab_rank > hfd_metab_rank)
# is MCD low-metab / high-fibro? (relative to its own diets)
mcd_low_metab  <- our[diet == "MCD", our_metab_rank] >= 3   # bottom half of 4
mcd_high_fibro <- our[diet == "MCD", our_fibro_rank] <= 2   # top half of 4

cat("\n=== Q1 RESOLUTION ===\n")
cat(sprintf("  MCD metabolic-arm rank = %d/4 ; fibrotic-arm rank = %d/4\n",
            mcd_metab_rank, our[diet=="MCD", our_fibro_rank]))
cat(sprintf("  MCD is LOW-metabolic (rank>=3): %s ; HIGH-fibrotic (rank<=2): %s\n",
            mcd_low_metab, mcd_high_fibro))
cat(sprintf("  Ranking by METABOLIC arm puts MCD BELOW both CDAHFD and HFD: %s\n", mcd_below_both))
cat(sprintf("    (MCD metab_rank %d vs CDAHFD %d, HFD %d)\n",
            mcd_metab_rank, cdahfd_metab_rank, hfd_metab_rank))
cat(sprintf("  Vacca MCD targets: metab m_DHPS=%.2f (expect ~0.55), fibro f_DHPS=%.2f (expect ~0.75)\n",
            our[diet=="MCD", vacca_metab_DHPS], our[diet=="MCD", vacca_fibro_DHPS]))

# Spearman of our per-arm proximity vs Vacca same-diet per-arm DHPS across 4 diets
cat(sprintf("\n  Across 4 diets: Spearman(our_metab_prox, vacca_metab_DHPS)=%.2f ; Spearman(our_fibro_prox, vacca_fibro_DHPS)=%.2f\n",
            sp(our$our_metab_prox, our$vacca_metab_DHPS),
            sp(our$our_fibro_prox, our$vacca_fibro_DHPS)))

fwrite(our, file.path(OUT, "q1_metab_fibro_2x2.csv"))
cat(sprintf("\nWritten: %s\n", file.path(OUT, "q1_metab_fibro_2x2.csv")))
cat(strrep("=", 78), "\n")
