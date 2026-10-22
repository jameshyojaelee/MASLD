suppressPackageStartupMessages(library(data.table))
repo <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
meta <- as.data.table(readRDS(file.path(repo, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds")))
old <- fread(cmd = paste("zcat", file.path(repo, "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/sealed_inputs/nonholdout_composition.tsv.gz")), check.names = FALSE)
new <- fread(cmd = paste("zcat", file.path(repo, "figures/candidates/composition-774-20261011T151616Z/composition/accepted_composition_774_plus_unchanged4.tsv.gz")), check.names = FALSE)
ros <- fread(file.path(repo, "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_plus_GSE174478/model_roster.tsv"))
cat("old rows", nrow(old), "unique", uniqueN(old$sample_id), "; new rows", nrow(new), "unique", uniqueN(new$sample_id), "\n")
cat("old by dataset\n"); print(old[, .N, by = dataset][order(dataset)])
cat("new by dataset\n"); print(new[, .N, by = dataset][order(dataset)])
lin <- setdiff(names(new), c("sample_id", "dataset"))
stopifnot(identical(lin, setdiff(names(old), c("sample_id", "dataset"))), length(lin) == 16)
cat("PCA record set = 774 meta records, by dataset x group, with earlier-fraction coverage\n")
meta[, in_old := sample_id %in% old$sample_id]
meta[, in_new := sample_id %in% new$sample_id]
meta[, in_model := sample_id %in% ros$sample_id]
print(meta[, .(n = .N, in_new = sum(in_new), in_old = sum(in_old), in_model = sum(in_model)), by = .(dataset, group_binary)][order(dataset, group_binary)])
cat("model roster by dataset\n"); print(ros[, .N, by = dataset])
cat("model roster covered: old", sum(ros$sample_id %in% old$sample_id), "new", sum(ros$sample_id %in% new$sample_id), "\n")
# per-cohort difference over the full PCA record set (fractions only; no model outcome)
m <- merge(old[, c("sample_id", lin), with = FALSE], new[, c("sample_id", "dataset", lin), with = FALSE], by = "sample_id", suffixes = c(".o", ".n"))
m <- m[sample_id %in% meta$sample_id | dataset == "GSE174478"]
d <- abs(as.matrix(m[, paste0(lin, ".o"), with = FALSE]) - as.matrix(m[, paste0(lin, ".n"), with = FALSE]))
m[, maxdiff := apply(d, 1, max)]
print(m[, .(n_shared = .N, n_changed = sum(maxdiff > 1e-6), max_abs_diff = max(maxdiff)), by = dataset][order(dataset)])
cat("old fractions below 1e-6 in the 457 model records:", sum(as.matrix(old[sample_id %in% ros$sample_id, ..lin]) < 1e-6), "; exact zeros:", sum(as.matrix(old[sample_id %in% ros$sample_id, ..lin]) == 0), "\n")
cat("new fractions below 1e-6 in the 457 model records:", sum(as.matrix(new[sample_id %in% ros$sample_id, ..lin]) < 1e-6), "; exact zeros:", sum(as.matrix(new[sample_id %in% ros$sample_id, ..lin]) == 0), "\n")
