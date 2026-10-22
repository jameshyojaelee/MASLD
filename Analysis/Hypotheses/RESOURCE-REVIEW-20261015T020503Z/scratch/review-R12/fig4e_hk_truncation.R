suppressPackageStartupMessages({library(data.table); library(metafor)})
cat("metafor", as.character(packageVersion("metafor")), "\n")
base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/"
out <- list()
for (arm in c("new774_plus_GSE174478", "new774_three_cohorts")) {
  ce <- fread(file.path(base, arm, "cohort_effects.tsv"))
  me <- if (file.exists(file.path(base, arm, "meta_effects.tsv"))) fread(file.path(base, arm, "meta_effects.tsv")) else NULL
  ce <- ce[estimable == TRUE & is.finite(beta) & is.finite(se_hc3) & se_hc3 > 0]
  g <- unique(ce[, .(feature_id, axis, model_kind)])
  res <- rbindlist(lapply(seq_len(nrow(g)), function(i) {
    d <- ce[feature_id == g$feature_id[i] & axis == g$axis[i] & model_kind == g$model_kind[i]]
    if (nrow(d) < 3) return(NULL)
    f1 <- rma(yi = d$beta, sei = d$se_hc3, method = "REML", test = "knha")
    f2 <- rma(yi = d$beta, sei = d$se_hc3, method = "REML", test = "adhoc")
    f0 <- rma(yi = d$beta, sei = d$se_hc3, method = "REML", test = "z")
    data.table(feature_id = g$feature_id[i], axis = g$axis[i], model_kind = g$model_kind[i], k = nrow(d),
               beta = as.numeric(f1$b), se_knha = f1$se, p_knha = f1$pval, se_adhoc = f2$se, p_adhoc = f2$pval,
               se_z = f0$se, p_z = f0$pval, tau2 = f1$tau2, QE = f1$QE)
  }))
  res[, family_size := .N, by = model_kind]
  res[, q_knha := p.adjust(p_knha, "BH"), by = model_kind]
  res[, q_adhoc := p.adjust(p_adhoc, "BH"), by = model_kind]
  res[, arm := arm]
  if (!is.null(me)) {
    chk <- merge(res, me[, .(feature_id, axis, model_kind, q_value, se_meta)], by = c("feature_id", "axis", "model_kind"))
    cat(arm, "max |q_knha - stored q|:", max(abs(chk$q_knha - chk$q_value)), " max |se diff|:", max(abs(chk$se_knha - chk$se_meta)), "\n")
  }
  cat(arm, "\n"); print(res[, .(family = .N, calls_knha = sum(q_knha < 0.05), calls_adhoc = sum(q_adhoc < 0.05),
      frac_se_knha_below_se_z = mean(se_knha < se_z)), by = .(model_kind, axis)])
  out[[arm]] <- res
}
all <- rbindlist(out)
fwrite(all, "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R12/fig4e_hk_truncation.tsv", sep = "\t")
duc <- "hotspot_hepatocytes_48f39dd4d817a10e"; ecm <- "hotspot_hepatocytes_f05c535ae5bbc0b9"
print(all[feature_id %in% c(duc, ecm) & axis == "fibrosis", .(arm, feature_id, model_kind, beta, se_knha, se_z, se_adhoc, q_knha, q_adhoc, QE, tau2)])
print(all[arm == "new774_plus_GSE174478" & model_kind == "composition_pc" & q_knha < 0.05, .(feature_id, axis, beta, se_knha, se_z, q_knha, q_adhoc, QE)])
