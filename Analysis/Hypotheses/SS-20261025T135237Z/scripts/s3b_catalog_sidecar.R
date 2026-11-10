# SS job S3b (prespecified; PRESPEC_SS.md F-SS3 per-gene states; rules in PRESPEC_SS_addendum_S3b.md). Catalog sidecar:
# one row per r8 Gene Catalog gene (frozen passport_gene_index.parquet) with the SteatoSITE stage_num coefficient
# (`~ batch + inferred_sex + stage_num`, limma voomWithQualityWeights), its SE and BH q over kept genes, the H51-style
# F3F4 per-gene state, and the new column steatosite_stage_state. Reads S3a outputs; fits nothing. r8 classes are copied
# and checked unchanged; no class is edited.
# Unit: one library per genetic individual (as S3a). Family for q: all genes kept by S3a filterByExpr.
# Usage: Rscript s3b_catalog_sidecar.R --s3a <S3a out dir> --input <S2 dir> --stage <stage tsv> --out <new dir>
#        [--allow-unhashed-prespec]
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3a_common.R")
A <- parse_args(commandArgs(trailingOnly = TRUE), c("s3a", "input", "stage", "out"), "allow-unhashed-prespec")
PY <- "/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
OUT <- new_outdir(A$out)
log_msg("S3b start. S3a =", A$s3a, " OUT =", OUT)
verify_prespec(OUT, isTRUE(A$`allow-unhashed-prespec`))
verify_frozen(OUT)
# The S3a run must have used these exact inputs
ih <- fread(file.path(A$s3a, "input_sha256.tsv"))
now <- c(units = sha256(file.path(A$input, "units.tsv")), dge = sha256(file.path(A$input, "dge_included.rds")), stage = sha256(A$stage))
stopifnot(all(ih[match(names(now), id), sha256] == now))
sc <- list.files(file.path(SS, "scripts"), pattern = "^s3[ab]_.*\\.R$", full.names = TRUE)
fwrite(data.table(script = basename(sc), sha256 = vapply(sc, sha256, "")), file.path(OUT, "script_sha256.tsv"), sep = "\t")
res <- jsonlite::read_json(file.path(A$s3a, "result.json"), simplifyVector = TRUE)
f3f4 <- as.data.table(res$F_SS1)[contrast == "F3F4"]
adopt <- isTRUE(f3f4$holm_pass) && isTRUE(res$PC_pass) && !isTRUE(res$debug)
why <- if (isTRUE(res$debug)) "not assigned: S3a debug run" else if (!isTRUE(res$PC_pass)) "not assigned: PC failed" else
  if (!isTRUE(f3f4$holm_pass)) "not assigned: F-SS1 F3F4 did not pass Holm" else NA_character_

# Catalog (parquet -> tsv with pyarrow, isolated interpreter)
cat_tsv <- file.path(OUT, "catalog_r8_input.tsv")
py <- sprintf("import pyarrow.parquet as pq, pyarrow.csv as pc; t = pq.read_table('%s'); pc.write_csv(t, '%s', pc.WriteOptions(delimiter='\\t'))",
              file.path(FROZEN, "catalog/passport_gene_index.parquet"), cat_tsv)
stopifnot(system2(PY, c("-I", "-c", shQuote(py))) == 0L)
catg <- fread(cat_tsv, sep = "\t")
stopifnot(nrow(catg) == 28233L, !anyDuplicated(catg$ensembl_id))
cls0 <- catg[, .(ensembl_id, primary_evidence_class)]

ann <- load_ann()
num <- fread(file.path(A$s3a, "gene_num.tsv.gz"))[coef == "stage_num"]
pg <- fread(file.path(A$s3a, "per_gene_states.tsv.gz"))[contrast == "F3F4"]
sc_tab <- merge(catg[, .(ensembl_id, symbol, primary_evidence_class, role_hypothesis)], ann[, .(ensembl_id = gene_id_base, gene_id_versioned)], by = "ensembl_id", all.x = TRUE)
sc_tab[, in_release_universe := !is.na(gene_id_versioned)]
sc_tab <- merge(sc_tab, num[, .(gene_id_versioned = gene, steatosite_stage_num_coef = logFC, steatosite_stage_num_SE = SE,
  steatosite_stage_num_t = t, steatosite_stage_num_p = P, steatosite_stage_num_q = FDR)], by = "gene_id_versioned", all.x = TRUE)
sc_tab[, kept_steatosite := !is.na(steatosite_stage_num_coef)]
sc_tab <- merge(sc_tab, pg[, .(gene_id_versioned = gene, D0_F3F4_log2FC = D0_log2FC, steatosite_F3F4_log2FC = G_log2FC, steatosite_F3F4_BH_q = G_BH_q,
  steatosite_F3F4_MDE_g = MDE_g, steatosite_F3F4_e_g = e_g, steatosite_F3F4_state_computed = state)], by = "gene_id_versioned", all.x = TRUE)
sc_tab[, steatosite_stage_state := fifelse(!in_release_universe, "not in release gene universe",
  fifelse(is.na(steatosite_F3F4_state_computed), "not a discovery F3F4 hit",
  fifelse(steatosite_F3F4_state_computed == "not measurable in SteatoSITE" | adopt, steatosite_F3F4_state_computed, why)))]
setcolorder(sc_tab, c("ensembl_id", "symbol", "gene_id_versioned", "primary_evidence_class", "role_hypothesis", "steatosite_stage_state"))
setorder(sc_tab, ensembl_id)
# r8 classes unchanged
chk <- merge(cls0, sc_tab[, .(ensembl_id, cls1 = primary_evidence_class)], by = "ensembl_id")
stopifnot(nrow(chk) == 28233L, identical(chk$primary_evidence_class, chk$cls1))
fwrite(sc_tab, file.path(OUT, "steatosite_catalog_sidecar.tsv.gz"), sep = "\t")
cnt <- sc_tab[, .N, by = .(primary_evidence_class, steatosite_stage_state)][order(primary_evidence_class, steatosite_stage_state)]
fwrite(cnt, file.path(OUT, "state_by_class.tsv"), sep = "\t")
unlink(cat_tsv)

rd <- c("# SteatoSITE Catalog sidecar (S3b)", "",
  "One row per r8 Gene Catalog gene (28,233; `gene-catalog-corrected-r8-20261009T142027Z`). Adds columns; r8 classes are copied and verified unchanged.", "",
  paste0("S3a source: ", normalizePath(A$s3a), if (isTRUE(res$debug)) " (DEBUG run: not a result)" else ""),
  paste0("SteatoSITE units: ", res$n_units, " staged genetic individuals; genes kept: ", res$n_kept_genes, "."),
  paste0("F-SS1 F3F4 Holm pass: ", isTRUE(f3f4$holm_pass), "; PC pass: ", isTRUE(res$PC_pass), "; states adopted into steatosite_stage_state: ", adopt, "."), "",
  "Columns:", "",
  "- `steatosite_stage_num_coef`, `_SE`, `_t`, `_p`, `_q`: log2 expression change per fibrosis-stage code (1-4) in `~ batch + inferred_sex + stage_num`; q = BH over kept genes. NA = not kept by filterByExpr or not in the release universe.",
  "- `D0_F3F4_log2FC`, `steatosite_F3F4_log2FC`, `steatosite_F3F4_BH_q`, `_MDE_g`, `_e_g`: H51 per-gene quantities for frozen D0 F3F4 BH hits (q BH within kept hits; MDE_g = 3.65 x SE; e_g = s_c x |D0 log2FC|, s_c frozen from H51).",
  "- `steatosite_F3F4_state_computed`: H51 state as computed, regardless of the family result.",
  "- `steatosite_stage_state`: replicated / opposite sign / not replicated at adequate power / underpowered / not measurable in SteatoSITE for D0 F3F4 hits, only when F-SS1 F3F4 passes Holm and PC passes; otherwise `not assigned: <reason>`. Non-hits: `not a discovery F3F4 hit`; genes outside the 24,115 release universe: `not in release gene universe`.", "",
  "Cross-sectional, FFPE total RNA; a state is a replication label for a stage association, not a causal or class change.", "",
  "Counts by r8 class:", "", tbl(cnt))
writeLines(rd, file.path(OUT, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
log_msg("S3b DONE; adopted =", adopt)
