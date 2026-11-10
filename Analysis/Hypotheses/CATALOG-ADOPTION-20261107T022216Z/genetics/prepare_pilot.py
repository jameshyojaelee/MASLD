#!/usr/bin/env python3
"""Snapshot the unchanged fit producer and inject isolated diagnostic exports."""
import csv
import hashlib
import json
from pathlib import Path

here = Path(__file__).resolve().parent
root = here.parents[3]
fm = root / "GWAS/finemapping"
source = fm / "src/06_susie_coloc.R"
runner = source.read_text()
anchor = "          if (!is.null(s_gwas) && isTRUE(s_gwas$converged)) {"
assert runner.count(anchor) == 1
injection = '''
            # Diagnostic exports only: unchanged z, LD and fits, before colocalization.
            pilot_dir <- Sys.getenv("GENETICS_PILOT_EXPORT_DIR")
            stopifnot(nzchar(pilot_dir))
            dir.create(pilot_dir, recursive = TRUE, showWarnings = FALSE)
            pilot_tag <- paste(gwas_name, gene_id, sep = "__")
            saveRDS(list(gwas_fit = s_gwas, eqtl_fit = s_eqtl, z = gwas_z,
              R_regularized = R_reg, aligned_summary = ss_ld, snp_ids = snp_ids,
              eqtl_fit_path = eqtl_rds, n_eff = gwas_n_eff, L = SUSIE_L,
              ld_diagonal = LD_REGULARIZE, residual_variance = 1,
              p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12,
              build = "GRCh37/hg19", allele_harmonization = c("match_or_flip", "frequency_resolved_palindromes"),
              lambda_s_locus = susie_lambda_s, seed = 42),
              file.path(pilot_dir, paste0(pilot_tag, "__inputs_and_fits.rds")))
            data.table::fwrite(ss_ld, file.path(pilot_dir,
              paste0(pilot_tag, "__aligned_variants.tsv")), sep = "\\t")
            pilot_cs <- rbindlist(lapply(seq_along(s_gwas$sets$cs), function(k) {
              ix <- s_gwas$sets$cs[[k]]
              data.table(signal_index = s_gwas$sets$cs_index[k],
                variant = snp_ids[ix], alpha = s_gwas$alpha[s_gwas$sets$cs_index[k], ix])
            }))
            data.table::fwrite(pilot_cs, file.path(pilot_dir,
              paste0(pilot_tag, "__gwas_credible_sets.tsv")), sep = "\\t")
'''
with (here / "pilot_runner_snapshot.R").open("x") as handle:
    handle.write(runner.replace(anchor, anchor + injection))
targets = [dict(gwas_name="MVP_ALT_EUR", chr=17, ensembl="ENSG00000133030", gene="MPRIP", role="selected_LD_sensitive"),
           dict(gwas_name="MVP_NAFLD_EUR", chr=15, ensembl="ENSG00000069667", gene="RORA", role="reconstruction_control"),
           dict(gwas_name="UKBB_ALT", chr=2, ensembl="ENSG00000163586", gene="FABP1", role="reconstruction_control")]
with (here / "pilot_targets.tsv").open("x", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(targets[0]), delimiter="\t")
    writer.writeheader(); writer.writerows(targets)
paths = [source, fm / "src/coloc_signal_eligibility.R", fm / "src/finemapping_functions.R",
         fm / "src/perf/get_ld_per_locus_cached.R", here / "pilot_runner_snapshot.R",
         here / "pilot_targets.tsv", here / "SPEC.md"]
for target in targets:
    paths.append(fm / "results/eqtl_susie_polyfun" / f"chr{target['chr']}" / f"{target['ensembl']}_susie.rds")
with (here / "pilot_sources.sha256").open("x") as handle:
    for path in paths:
        handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path}\n")
print(json.dumps(dict(status="prepared", target_fits=3, original_runner_sha256=hashlib.sha256(source.read_bytes()).hexdigest()), indent=2))
