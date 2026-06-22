#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# PROCON core calibration test. Synthetic: 8 modalities × 2000 genes × 24 programs.
#   concordant programs (signal in ALL modalities) → metafor pooled p sig, I² low;
#   single-modality programs (signal in 1)          → ACAT sig (omnibus), metafor weaker;
#   null programs                                    → both p ~ uniform (FPR≈0.05).
# Validates cameraPR roll-up + metafor + ACAT before wiring real 46d modality stats.
# Lightweight — safe off-node.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(data.table) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/heterogeneity_program/lib/procon.R"))

set.seed(11)
G <- 2000L; M <- 8L; P <- 24L; psize <- 40L
genes <- sprintf("g%04d", seq_len(G))
sets <- lapply(seq_len(P), function(i) genes[((i-1)*psize + 1):(i*psize)])
names(sets) <- sprintf("P%02d", seq_len(P))
# sparse, sign-BALANCED signal (mirrors real data: most genes null, up≈down) so the
# COMPETITIVE cameraPR background stays centred and the true-null programs are calibrated.
cls <- c(rep("conc_up",2), rep("conc_dn",2), rep("single",2), rep("null",18))   # 24
mod_stats <- lapply(seq_len(M), function(m) {
  s <- setNames(rnorm(G), genes)
  for (i in seq_len(P)) {
    if (cls[i] == "conc_up")                s[sets[[i]]] <- s[sets[[i]]] + 1.2
    if (cls[i] == "conc_dn")                s[sets[[i]]] <- s[sets[[i]]] - 1.2
    if (cls[i] == "single" && m == 1L)      s[sets[[i]]] <- s[sets[[i]]] + 1.6
  }
  s
})

# cameraPR per modality → long (program_id, modality, z, p)
long <- rbindlist(lapply(seq_len(M), function(m) {
  r <- program_camera(mod_stats[[m]], sets); r[, modality := m]; r
}))
res <- procon_combine(long)
res[, cls := cls[match(program_id, names(sets))]]
res[, q_meta := p.adjust(procon_p_meta, "BH")]
res[, q_acat := p.adjust(procon_p_acat, "BH")]
res[, label := procon_label(procon_effect, procon_I2, q_meta, n_modalities, sign_frac)]

cat("── PROCON per-class results (median over 8 programs each) ──\n")
print(res[, .(meta_p = round(median(procon_p_meta),4), acat_p = round(median(procon_p_acat),4),
              I2 = round(median(procon_I2, na.rm=TRUE),1),
              n_meta_sig = sum(q_meta < 0.05), n_acat_sig = sum(q_acat < 0.05)), by = cls])

conc <- res[cls %in% c("conc_up","conc_dn")]; sing <- res[cls=="single"]; nul <- res[cls=="null"]
pass <- mean(conc$q_meta < 0.05) >= 0.75 && median(conc$procon_I2, na.rm=TRUE) < 50 &&
        mean(sing$q_acat < 0.05) >= 0.75 &&
        mean(nul$procon_p_meta < 0.05) <= 0.15 && mean(nul$procon_p_acat < 0.05) <= 0.15
cat(sprintf("\nconcordant: metafor-sig %.0f%% (I²=%.0f) | single: ACAT-sig %.0f%% | null FPR meta %.0f%% acat %.0f%%\n",
            100*mean(conc$q_meta<0.05), median(conc$procon_I2,na.rm=TRUE),
            100*mean(sing$q_acat<0.05), 100*mean(nul$procon_p_meta<0.05), 100*mean(nul$procon_p_acat<0.05)))
cat("=== ", if (pass) "PASS" else "FAIL",
    " (metafor finds concordant; ACAT finds single-modality; null calibrated) ===\n", sep="")
if (!pass) quit(status = 1)
