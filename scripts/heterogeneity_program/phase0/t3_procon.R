#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# T3 — PROCON: program-level multimodal convergence over the frozen dictionary.
#
# Per-modality per-gene signed ranking statistic → cameraPR competitive
# enrichment per program → metafor REML (pooled signed effect + I², PRIMARY) +
# ACAT (omnibus, CO-PRIMARY). Genetic (COLOC PP.H4) and the T1 dispersion axis
# (|dv_t|) enter as one-sided enrichment arms. Heuristic convergence score, NOT
# a posterior. Permutation calibration (gene-label, decile-matched) gated by N_PERM.
#
# Modes: N_PERM=0 observed (fast); >0 adds the random-matched-set permutation null.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(data.table) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/heterogeneity_program/lib/procon.R"))
N_PERM  <- as.integer(Sys.getenv("N_PERM", "0"))
P  <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")

atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
dv    <- fread(file.path(P, "dv_atlas_columns.tsv"))
atlas <- merge(atlas, dv[, .(human_symbol, dv_t)], by = "human_symbol", all.x = TRUE)
memb  <- fread(file.path(P, "program_gene_membership.tsv"))
sets  <- split(memb$gene, memb$program_id)
cat(sprintf("atlas %d genes; %d programs\n", nrow(atlas), length(sets)))

sgnl <- function(lfc, padj) sign(lfc) * -log10(pmax(padj, 1e-300))
named <- function(v) { x <- v; names(x) <- atlas$human_symbol; x[is.finite(x)] }

# signed DE-type modalities (metafor directional convergence)
mods_signed <- list(
  bulk       = named(atlas$bulk_tstat),
  mouse      = named(sgnl(atlas$mouse_meta_logFC, atlas$mouse_meta_padj)),
  proteomics = named(sgnl(atlas$best_protein_logFC, atlas$best_protein_padj)),
  sc_hep     = named(sgnl(atlas$sc_hepatocyte_logFC, atlas$sc_hepatocyte_padj)),
  spatial    = named(sgnl(atlas$spatial_govaere2026_cosmx_hep_mash_logfc,
                          atlas$spatial_govaere2026_cosmx_hep_mash_padj)))
cat("signed-modality coverage:\n")
for (m in names(mods_signed)) cat(sprintf("  %-11s %d genes\n", m, length(mods_signed[[m]])))

# per-modality coverage → sei (review H1: down-weight thin modalities so a
# 848-gene spatial arm carries ∝ its coverage, not equal to bulk's 27k).
cov_m <- sapply(mods_signed, length); sei_m <- sqrt(max(cov_m) / cov_m)
cat("modality coverage / sei (metafor down-weight ∝ 1/sqrt coverage):\n")
for (m in names(sei_m)) cat(sprintf("  %-11s cov=%5d  sei=%.2f\n", m, cov_m[[m]], sei_m[[m]]))

long <- rbindlist(lapply(names(mods_signed), function(m) {
  r <- program_camera(mods_signed[[m]], sets); if (is.null(r)) return(NULL)
  r[, `:=`(modality = m, sei = sei_m[[m]])]; r }))
res <- procon_combine(long)   # metafor (coverage-weighted) = CONVERGENCE headline

# one-sided (upper-tail) enrichment arms (review M3): COLOC PP.H4 + dispersion |dv_t|
gen <- program_camera(named(atlas$coloc_susie_best_pp4), sets)
dis <- program_camera(named(abs(atlas$dv_t)), sets)
gen[, p1 := pnorm(z, lower.tail = FALSE)]; dis[, p1 := pnorm(z, lower.tail = FALSE)]
res <- merge(res, gen[, .(program_id, genetic_z = z, genetic_p1 = p1)], by = "program_id", all.x = TRUE)
res <- merge(res, dis[, .(program_id, dispersion_z = z, dispersion_p1 = p1)], by = "program_id", all.x = TRUE)

# omnibus ANY-SIGNAL ACAT (review H2: sign-AGNOSTIC, NOT convergence — convergence
# is metafor only) over signed two-sided p's + genetic/dispersion ONE-sided p's.
allp <- rbind(long[, .(program_id, p)], gen[, .(program_id, p = p1)], dis[, .(program_id, p = p1)])
anysig <- allp[, .(procon_p_anysignal = acat(p)), by = program_id]
res <- merge(res, anysig, by = "program_id", all.x = TRUE)

res[, source := sub("_.*", "", program_id)]
res[, q_meta := p.adjust(procon_p_meta, "BH")]                 # parametric metafor BH
res[, q_anysignal := p.adjust(procon_p_anysignal, "BH")]       # any-modality omnibus (NOT convergence)

# ── Gene-label permutation calibration (rigor layer) ──────────────────────────
# Null of the convergence effect from RANDOM gene sets of matched size: permute
# the gene→stat assignment within each modality, recompute cameraPR→metafor per
# program. Calibrates "more convergent than a random set" beyond the parametric BH.
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
if (N_PERM > 0) {
  cat(sprintf("\nGene-label permutation calibration: N_PERM=%d, cores=%d\n", N_PERM, N_CORES))
  obs_eff <- setNames(abs(res$procon_effect), res$program_id)
  set.seed(42); perm_seeds <- sample.int(1e7, N_PERM)
  one_perm <- function(sd) {
    set.seed(sd)
    lb <- rbindlist(lapply(names(mods_signed), function(m) {
      sm <- mods_signed[[m]]; names(sm) <- sample(names(sm))        # permute gene labels
      r <- program_camera(sm, sets); if (is.null(r)) return(NULL)
      r[, `:=`(modality = m, sei = sei_m[[m]])]; r }))
    rb <- procon_combine(lb)
    e <- setNames(abs(rb$procon_effect), rb$program_id)[names(obs_eff)]
    as.integer(!is.na(e) & e >= obs_eff)
  }
  exc <- Reduce(`+`, parallel::mclapply(perm_seeds, one_perm, mc.cores = N_CORES))
  res[, procon_p_perm := (exc[match(program_id, names(obs_eff))] + 1) / (N_PERM + 1)]
  res[, procon_q_perm := p.adjust(procon_p_perm, "BH")]
  cat(sprintf("  CONVERGENT at perm-q<0.05: %d  (vs parametric metafor-q<0.05: %d)\n",
              res[procon_q_perm < 0.05, .N], res[q_meta < 0.05, .N]))
}
res[, conv_q := if ("procon_q_perm" %in% names(res)) procon_q_perm else q_meta]   # headline convergence
res[, label := procon_label(procon_effect, procon_I2, conv_q, n_modalities, sign_frac)]
setorder(res, conv_q, procon_p_anysignal)
fwrite(res, file.path(P, sprintf("procon_convergence%s.tsv", if (N_PERM>0) "" else "_observed")), sep = "\t")

cat(sprintf("\n── PROCON: %d programs | CONVERGENT (conv-q<0.05): %d | any-signal (ACAT q<0.05): %d ──\n",
            nrow(res), res[conv_q<0.05,.N], res[q_anysignal<0.05,.N]))
cat("4-state labels:\n"); print(res[, .N, by = label][order(-N)])
cat("\nTop 15 convergent programs:\n")
print(res[order(conv_q)][1:15, .(program_id, n_mod = n_modalities, effect = round(procon_effect,2),
       I2 = round(procon_I2,0), conv_q = signif(conv_q,2), genetic_z = round(genetic_z,1),
       disp_z = round(dispersion_z,1), label)])
