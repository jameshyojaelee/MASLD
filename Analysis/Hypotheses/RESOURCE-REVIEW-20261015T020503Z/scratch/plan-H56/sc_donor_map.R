# Header-only check: map hepatocyte pseudobulk runs to donors and match the corrected roster (no expression read).
suppressPackageStartupMessages(library(data.table))
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
m <- build_srr_to_donor_map(BASE)
F <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk")
hdr <- names(fread(cmd = paste("zcat", file.path(F, "pseudobulk_Hepatocytes.tsv.gz"), "| head -2")))[-1]
meta <- fread(file.path(F, "sample_celltype_metadata.tsv"))[cell_type == "Hepatocytes"]
cat("hepatocyte pseudobulk run columns", length(hdr), "; metadata rows", nrow(meta), "; header in metadata", sum(hdr %in% meta$sample), "\n")
don <- ifelse(hdr %in% names(m), m[hdr], NA)
meta[, donor := ifelse(sample %in% names(m), m[sample], NA_character_)]
ros <- fread(file.path(BASE, "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/current_source_diagnosis_donor_roster.tsv"))
cat("roster donors", nrow(ros), "; non-excluded", ros[exclude == FALSE | exclude == "False", .N], "\n")
cat("runs mapped by pairing files", sum(!is.na(don)), "\n")
# pass-through runs: show dataset and example IDs
pt <- meta[is.na(donor)]
print(pt[, .N, by = dataset])
print(head(ros[dataset %in% c("GSE174748", "GSE189600"), .(donor, dataset, disease_stage_coarse)], 10))
print(head(pt[dataset %in% c("GSE174748", "GSE189600"), .(sample, dataset, n_cells, disease_stage_coarse)], 10))
mapped <- unique(na.omit(meta$donor))
cat("distinct mapped donors", length(mapped), "; in roster", sum(mapped %in% ros$donor), "\n")
# donor-level hepatocyte cell totals for non-excluded roster donors
meta[is.na(donor), donor := sample]
dn <- meta[, .(n_cells = sum(n_cells)), by = .(donor, dataset)]
r <- merge(ros[exclude == FALSE | exclude == "False"], dn, by = c("donor", "dataset"), all.x = TRUE)
cat("non-excluded roster donors with hepatocyte pseudobulk", r[!is.na(n_cells), .N], "; with >=50 hepatocytes", r[n_cells >= 50, .N], "\n")
print(r[n_cells >= 50, .N, by = .(dataset, healthy = disease_stage_coarse == "Healthy")][order(dataset)])
print(r[is.na(n_cells) | n_cells < 50, .(donor, dataset, disease_stage_coarse, n_cells)])
