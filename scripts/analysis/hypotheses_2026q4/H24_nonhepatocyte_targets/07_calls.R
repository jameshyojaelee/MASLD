#!/usr/bin/env Rscript
# H24 step 7: apply the prespecified call rule (PRESPEC_H24.md) and summarize.
# Usage: Rscript 07_calls.R <OUT_DIR>
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
args <- commandArgs(trailingOnly = TRUE); OUT <- args[[1]]
res <- file.path(OUT, "results")
PIP_MIN <- 0.1; PP_CALL <- 0.8; PP_SENS <- 0.5
ZHU_SET <- c("APOC1", "APOC1P1", "APOE", "TOMM40", "APOC2", "APOC4-APOC2")

tests <- fread(file.path(res, "l1_l2_link_tests.tsv"))
diag <- fread(file.path(res, "l1_lineage_diagnostics.tsv"))
col_ok <- file.exists(file.path(res, "coloc_results.tsv")) && file.size(file.path(res, "coloc_results.tsv")) > 10
col <- if (col_ok) fread(file.path(res, "coloc_results.tsv")) else data.table()
if (nrow(col) && !"gse296875_peak_id" %in% names(col)) col[, gse296875_peak_id := NA_character_]
cstat <- fread(file.path(res, "coloc_status.tsv"))
cloc <- fread(file.path(res, "coloc_loci.tsv"))
phys <- fread(file.path(OUT, "inventory_v2/physical_loci.tsv"))
vpk <- fread(file.path(OUT, "inventory_v2/variant_peaks.tsv"))
cs <- jsonlite::fromJSON(file.path(res, "coloc_summary.json"))

# ---- S2: coloc check at control loci (L0) ----
ctrl <- cloc[locus_set %in% c("CTRL_L0", "REF_PNPLA3")]
mach <- rbindlist(lapply(seq_len(nrow(ctrl)), function(i) {
  rel <- strsplit(ctrl$release_direct_coloc_genes[i], ";")[[1]]
  x <- if (nrow(col)) col[layer == "L0" & physical_locus == ctrl$physical_locus[i]] else data.table()
  st <- cstat[layer == "L0" & physical_locus == ctrl$physical_locus[i]]$status
  hit <- if (nrow(x)) x[symbol %in% rel] else data.table()
  data.table(physical_locus = ctrl$physical_locus[i], release_genes = ctrl$release_direct_coloc_genes[i],
             l0_status = if (length(st)) st[1] else NA_character_,
             best_release_gene = if (nrow(hit)) hit[which.max(PP.H4)]$symbol else NA_character_,
             best_release_pp_h4 = if (nrow(hit)) max(hit$PP.H4) else NA_real_)
}))
mach[, recovered := !is.na(best_release_pp_h4) & best_release_pp_h4 > PP_SENS]
machinery_pass <- any(mach$recovered)
fwrite(mach, file.path(res, "control_l0_recovery.tsv"), sep = "\t")

# ---- per-pair layers ----
h <- copy(tests)
h[, call_eligible := l1_support & r > 0 & !is.na(max_cs_pip_in_peak) & max_cs_pip_in_peak >= PIP_MIN]
if (nrow(col)) {
  l3 <- col[layer == "L3", .(l3_max_pp_h4 = max(PP.H4), l3_best_dataset = dataset[which.max(PP.H4)]),
            by = .(physical_locus, ensembl = gene_id, lineage = proxy_lineage)]
  l0 <- col[layer == "L0", .(l0_max_pp_h4 = max(PP.H4)), by = .(physical_locus, ensembl = gene_id)]
  l5 <- col[layer == "L5", .(l5_max_pp_h4 = max(PP.H4)), by = .(physical_locus, ensembl = gene_id)]
  l4 <- col[layer == "L4", .(l4_max_pp_h4 = max(PP.H4)), by = .(physical_locus, peak_id = gse296875_peak_id)]
  h <- merge(h, l3, by = c("physical_locus", "ensembl", "lineage"), all.x = TRUE)
  h <- merge(h, l0, by = c("physical_locus", "ensembl"), all.x = TRUE)
  h <- merge(h, l5, by = c("physical_locus", "ensembl"), all.x = TRUE)
  h <- merge(h, l4, by = c("physical_locus", "peak_id"), all.x = TRUE)
} else h[, `:=`(l3_max_pp_h4 = NA_real_, l3_best_dataset = NA_character_, l0_max_pp_h4 = NA_real_,
                l5_max_pp_h4 = NA_real_, l4_max_pp_h4 = NA_real_)]
l3_usable <- isTRUE(cs$network_ok) && machinery_pass
h[, l2_support := !is.na(l2_support) & l2_support]
h[, l3_support := l3_usable & !is.na(l3_max_pp_h4) & l3_max_pp_h4 >= PP_CALL]
h[, l3_support_sens := l3_usable & !is.na(l3_max_pp_h4) & l3_max_pp_h4 >= PP_SENS]
h[, two_layer_call := call_eligible & (l2_support | l3_support)]
h[, two_layer_call_sens := call_eligible & (l2_support | l3_support_sens)]
h[, non_hep_call := two_layer_call & lineage != "hepatocyte"]
# S4: HSD17B13 non-hepatocyte call withheld
h[, s4_withheld := non_hep_call & symbol == "HSD17B13"]
h[s4_withheld == TRUE, `:=`(two_layer_call = FALSE, non_hep_call = FALSE)]
setorder(h, family, physical_locus, lineage, p_perm)
fwrite(h, file.path(res, "h24_pairs_all_layers.tsv"), sep = "\t")

# ---- locus table (H24) ----
inform <- diag[status == "tested" & l1_informative == TRUE, lineage]
hl <- h[family == "H24"]
loc_ids <- unique(phys[set == "H24", .(physical_locus, strict = h24_strict_no_any_coloc)])[, .(strict = any(strict)), by = physical_locus]
loc <- rbindlist(lapply(loc_ids$physical_locus, function(pl) {
  x <- hl[physical_locus == pl]
  st <- cstat[physical_locus == pl]
  fmt <- function(d) if (nrow(d)) paste(unique(sprintf("%s:%s", d$symbol, d$lineage)), collapse = ";") else ""
  data.table(physical_locus = pl, strict = loc_ids[physical_locus == pl]$strict,
             variant_peaks = nrow(vpk[set == "H24" & physical_locus == pl]),
             l1_tests = nrow(x), l1_lineages = paste(sort(unique(x$lineage)), collapse = ";"),
             l1_supported = sum(x$l1_support), call_eligible = sum(x$call_eligible),
             l1_supported_pairs = fmt(x[l1_support == TRUE]),
             l2_tested = sum(x$l2_status == "tested", na.rm = TRUE), l2_supported = sum(x$l2_support),
             l3_datasets_tested = sum(st$layer == "L3" & st$status == "tested"),
             l3_status = paste(unique(st[layer == "L3"]$status), collapse = ";"),
             calls = fmt(x[two_layer_call == TRUE]), non_hep_calls = fmt(x[non_hep_call == TRUE]),
             calls_sens = fmt(x[two_layer_call_sens == TRUE]),
             l3_any_ge_0.8 = if (nrow(col)) paste(unique(col[layer == "L3" & physical_locus == pl & PP.H4 >= PP_CALL,
                                                                 sprintf("%s:%s(%s)", symbol, proxy_lineage, dataset)]), collapse = ";") else "",
             l0_any_ge_0.8 = if (nrow(col)) paste(unique(col[layer == "L0" & physical_locus == pl & PP.H4 >= PP_CALL]$symbol), collapse = ";") else "",
             l5_any_ge_0.8 = if (nrow(col)) paste(unique(col[layer == "L5" & physical_locus == pl & PP.H4 >= PP_CALL]$symbol), collapse = ";") else "",
             l4_any_ge_0.8 = if (nrow(col)) paste(unique(col[layer == "L4" & physical_locus == pl & PP.H4 >= PP_CALL]$gse296875_peak_id), collapse = ";") else "")
}))
fwrite(loc, file.path(res, "h24_locus_calls.tsv"), sep = "\t")

# ---- descriptive: variant-peak accessibility per lineage (fragments per 1,000 nuclei, pooled) ----
units <- fread(file.path(OUT, "gse296875/units.tsv"))
dims <- fromJSON(file.path(OUT, "gse296875/dims.json"))
con <- file(file.path(OUT, "gse296875/atac_peak_by_unit.f32"), "rb")
acc <- rbindlist(lapply(unique(vpk$peak_index), function(p) {
  v <- vapply(units$unit_index, function(j) {
    seek(con, (as.numeric(j) * dims$atac[1] + (p - 1)) * 4); readBin(con, "numeric", 1, size = 4, endian = "little")
  }, numeric(1))
  d <- data.table(units[, .(lineage_id, nuclei, atac_in_peak_fragments)], count = v)
  d[, .(peak_index = p, fragments = sum(count), nuclei = sum(nuclei),
        per_1000_nuclei = 1000 * sum(count) / sum(nuclei),
        cpm_in_peaks = 1e6 * sum(count) / sum(atac_in_peak_fragments)), by = lineage_id]
}))
close(con)
acc <- merge(vpk[, .(peak_index, peak_id, set, physical_locus, max_cs_pip_in_peak)], acc, by = "peak_index")
fwrite(acc, file.path(res, "variant_peak_accessibility_by_lineage.tsv"), sep = "\t")

# ---- references ----
zhu <- h[family == "REF_ZHU"]
zhu[, in_candidate_set := symbol %in% ZHU_SET]
pn <- h[family == "REF_PNPLA3"]
hsd <- h[family == "H24" & physical_locus == "4:2+4:9"]
sm <- list(
  release_state = "COLOC r2 candidate susie_coloc_r2_20260924T053244Z; not adopted; corrected 50-study rerun active",
  h24_physical_loci = nrow(loc), strict_physical_loci = sum(loc$strict),
  loci_with_variant_peak = sum(loc$variant_peaks > 0), loci_l1_tested = sum(loc$l1_tests > 0),
  loci_with_l1_support = sum(loc$l1_supported > 0),
  loci_with_two_layer_call = sum(loc$calls != ""), loci_with_non_hep_call = sum(loc$non_hep_calls != ""),
  strict_loci_with_non_hep_call = sum(loc$strict & loc$non_hep_calls != ""),
  loci_with_two_layer_call_sens = sum(loc$calls_sens != ""),
  call_lineages = sort(unique(hl[two_layer_call == TRUE]$lineage)),
  h24_tests = nrow(hl), h24_tests_by_lineage = hl[, .N, by = lineage],
  h24_l1_supported = sum(hl$l1_support), h24_call_eligible = sum(hl$call_eligible),
  informative_lineages = inform, lineage_diagnostics = diag,
  network_ok = cs$network_ok, coloc_check_recovered = sprintf("%d/%d", sum(mach$recovered), nrow(mach)),
  coloc_check_pass = machinery_pass, l3_usable_for_calls = l3_usable,
  hsd17b13_locus = list(tests = nrow(hsd), l1_supported = hsd[l1_support == TRUE, .(symbol, lineage, r, q_bh)],
                        hsd17b13_rows = hsd[symbol == "HSD17B13", .(lineage, peak_id, r, p_perm, q_bh, l1_support)],
                        non_hep_call_withheld = any(h$s4_withheld)),
  ref_pnpla3 = list(tests = nrow(pn), l1_supported = pn[l1_support == TRUE, .(symbol, lineage, r, q_bh, l2_support)]),
  ref_zhu = list(tests = nrow(zhu), l1_supported = zhu[l1_support == TRUE, .(symbol, r, q_bh, in_candidate_set)],
                 candidate_genes_tested = sort(unique(zhu[in_candidate_set == TRUE]$symbol))))
write_json(sm, file.path(res, "h24_summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
print(loc)
writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_calls.txt"))
cat("calls done\n")
