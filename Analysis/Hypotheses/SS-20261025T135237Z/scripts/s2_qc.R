#!/usr/bin/env Rscript
# SS job S2, QC part (PRESPEC_SS section 3). Builds out_s2_qc/ per S2_CONTRACT.md. The stage file is never opened here;
# no expression value other than library totals, protein-coding detection, XIST and the Y gene is looked at.
# Usage: s2_qc.R <assembled_dir> <identity_dir> <out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR); library(jsonlite); library(igraph)})
argv <- commandArgs(trailingOnly = TRUE)
stopifnot(length(argv) == 3L)
asm <- argv[[1]]; idd <- argv[[2]]; out <- argv[[3]]
if (file.exists(out)) stop("refusing to overwrite ", out)
SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
STAR <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/results/SteatoSITE/alignments/star")
META <- file.path(ROOT, "data/external/SteatoSite/fastq/PRJEB58625_run_metadata.tsv")
REL <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ANN <- file.path(ROOT, "figures/candidates/resource-corrections-20261005-v1/annotation_resolution_r8/annotation_correspondence.tsv")

counts <- readRDS(file.path(asm, "counts_raw.rds"))
libs <- colnames(counts)
ann <- fread(ANN, colClasses = "character")[match(rownames(counts), gene_id)]
stopifnot(!anyNA(ann$gene_id))

# Batch: sequencing project code in the submitted file name (e.g. ..._15110FJ0009L01_1.fastq.gz).
meta <- fread(META, colClasses = "character")
# Only the runs actually aligned count (mate_order.tsv in each alignment directory; one library per sample, amendment A2).
meta[, batch := sub(".*_([0-9]{5}[A-Z]{2})[0-9]{4}L0[0-9]_[12]\\.fastq\\.gz$", "\\1", sub(";.*", "", submitted_ftp))]
stopifnot(all(grepl("^[0-9]{5}[A-Z]{2}$", meta$batch)))
aligned <- rbindlist(lapply(libs, function(s) data.table(sample_alias = s,
  run_accession = fread(file.path(STAR, s, "mate_order.tsv"), header = FALSE)$V1)))
aligned[, batch := meta$batch[match(run_accession, meta$run_accession)]]
stopifnot(!anyNA(aligned$batch), all(meta$sample_alias[match(aligned$run_accession, meta$run_accession)] == aligned$sample_alias))
bt <- aligned[, .(batch = batch[1], n_batch_codes = uniqueN(batch)), by = sample_alias]
stopifnot(all(libs %in% bt$sample_alias), all(bt$n_batch_codes == 1L))

star_pct <- function(s, key) {
  x <- readLines(file.path(STAR, s, paste0(s, ".Log.final.out")))
  as.numeric(sub("%", "", trimws(sub(".*\\|", "", grep(key, x, fixed = TRUE, value = TRUE)[1]))))
}
rpg <- function(s) {  # column 4 = read 1 antisense = the strand counted (-s 2), after per-run mate order (amendment A1)
  x <- fread(file.path(STAR, s, paste0(s, ".ReadsPerGene.out.tab")), header = FALSE)
  c(nofeature = x$V4[x$V1 == "N_noFeature"], genes = sum(x$V4[!grepl("^N_", x$V1)]))
}
u <- data.table(sample_alias = libs, batch = bt$batch[match(libs, bt$sample_alias)],
                n_batch_codes = bt$n_batch_codes[match(libs, bt$sample_alias)],
                assigned_fragments = as.numeric(colSums(counts)))
u[, log10_assigned := log10(assigned_fragments)]
u[, unique_pct := vapply(sample_alias, star_pct, numeric(1), key = "Uniquely mapped reads %")]
r <- t(vapply(u$sample_alias, rpg, numeric(2)))
u[, intron_exon_ratio := r[, "nofeature"] / r[, "genes"]]  # proxy, see DEVIATIONS.md (D-S2-1)
pc <- ann$gene_biotype == "protein_coding"
u[, pc_genes_ge10 := as.integer(colSums(counts[pc, ] >= 10))]
u[, pc_floor := median(pc_genes_ge10) - 4 * mad(pc_genes_ge10), by = batch]

# Q1
u[, q1_pass := assigned_fragments >= 5e6 & unique_pct >= 40 & pc_genes_ge10 >= pc_floor]

# Q2: release sex rule (01_sample_qc.R:220-240): k-means (seed 42, 2 centres, nstart 25) on log1p CPM of XIST and DDX3Y
# (RPS4Y1 if DDX3Y absent); female cluster = higher mean XIST CPM. Possibly mixed -> excluded.
cp <- cpm(counts)
xist <- which(ann$gene_name == "XIST"); ygene <- which(ann$gene_name == "DDX3Y")
if (!length(ygene)) ygene <- which(ann$gene_name == "RPS4Y1")
stopifnot(length(xist) == 1L, length(ygene) == 1L)
u[, xist_cpm := cp[xist, ]]; u[, y_gene_cpm := cp[ygene, ]]
set.seed(42)
km <- kmeans(cbind(log1p(u$xist_cpm), log1p(u$y_gene_cpm)), centers = 2, nstart = 25)
fc <- which.max(tapply(u$xist_cpm, km$cluster, mean))
u[, sex_cluster := fifelse(km$cluster == fc, "F", "M")]
ho <- fread(file.path(idd, "ss_hom_other_allele.tsv"))
u[, other_allele_fraction := ho$other_allele_fraction[match(sample_alias, ho$sample_alias)]]
u[, mixed_expression := log2(xist_cpm + 0.5) > 5 & log2(y_gene_cpm + 0.5) > 3]
u[, mixed_allele := !is.na(other_allele_fraction) & other_allele_fraction > 0.02]
u[, q2_pass := !mixed_expression & !mixed_allele]
u[, inferred_sex := fifelse(q2_pass, sex_cluster, NA_character_)]

# Q3: identity. Same individual = KING >= 0.354 on >= 1,000 shared sites, or on 500-999 sites (amendment A3).
cs <- fread(file.path(idd, "ss_called_sites.tsv"))
u[, identity_sites_called := cs$identity_sites_called_dp10[match(sample_alias, cs$sample_alias)]]
u[, identity_screenable := !is.na(identity_sites_called) & identity_sites_called >= 1000]
kp <- fread(file.path(idd, "kinship_pairs.tsv"))
same <- kp[relation == "same_individual"]
# Amendment A3: pairs with KING >= 0.354 on 500-999 shared sites also count as the same individual.
k5 <- fread(file.path(idd, "kinship_same_individual_500_999_sites.tsv"))
if (nrow(k5)) same <- rbind(same, k5[, .(id1, id2, pair_type, nsnp, kinship, relation = "same_individual")], fill = TRUE)
ext <- same[pair_type == "ss_existing"]
ext_hit <- unique(c(ext$id1[ext$id1 %in% libs], ext$id2[ext$id2 %in% libs]))
u[, matches_existing_cohort := sample_alias %in% ext_hit]
g <- graph_from_data_frame(same[pair_type == "ss_ss", .(id1, id2)], directed = FALSE, vertices = data.frame(name = libs))
comp <- components(g)$membership
u[, identity_cluster := paste0("SSID", comp[sample_alias])]
u[, cluster_size := .N, by = identity_cluster]
u[, cluster_matches_existing := any(matches_existing_cohort), by = identity_cluster]
# One library per person: most assigned fragments among Q1- and Q2-passing members.
u[, rank_in_cluster := frank(-assigned_fragments, ties.method = "first"), by = .(identity_cluster, ok = q1_pass & q2_pass)]
u[, q3_pass := !cluster_matches_existing & (rank_in_cluster == 1L | !(q1_pass & q2_pass))]
u[, include := q1_pass & q2_pass & q3_pass]
u[, rank_in_cluster := NULL]
u[, exclusion_reason := fifelse(!q1_pass, "Q1_library",
                        fifelse(!q2_pass, "Q2_possibly_mixed",
                        fifelse(cluster_matches_existing, "Q3_matches_existing_cohort",
                        fifelse(!q3_pass, "Q3_same_individual_lower_depth", ""))))]
u[, analysis_unit_id := paste0("SS_", identity_cluster)]
setcolorder(u, c("sample_alias", "analysis_unit_id", "batch", "inferred_sex", "assigned_fragments", "log10_assigned",
                 "unique_pct", "intron_exon_ratio", "q1_pass", "q2_pass", "q3_pass", "include", "exclusion_reason"))
inc <- u[include == TRUE]
stopifnot(!anyDuplicated(inc$analysis_unit_id), !anyNA(inc$inferred_sex))

dir.create(file.path(out, "genotypes"), recursive = TRUE)
fwrite(u, file.path(out, "units.tsv"), sep = "\t")

rel_genes <- rownames(readRDS(file.path(REL, "merged_dge.rds")))
stopifnot(length(rel_genes) == 24115L)
m <- matrix(0L, length(rel_genes), nrow(inc), dimnames = list(rel_genes, inc$sample_alias))
hit <- intersect(rel_genes, rownames(counts))
m[hit, ] <- counts[hit, inc$sample_alias]
dge <- DGEList(counts = m, samples = as.data.frame(inc), genes = data.frame(gene_id_versioned = rel_genes))
dge$samples$lib.size <- colSums(m)
saveRDS(dge, file.path(out, "dge_included.rds"))
saveRDS(counts[, inc$sample_alias], file.path(out, "counts_all_genes.rds"))

for (s in inc$sample_alias) {
  for (ext in c(".vcf.gz", ".vcf.gz.tbi")) {
    src <- file.path(SS, "out_s1c_genotypes/vcf", paste0(s, ext))
    if (file.exists(src)) file.symlink(src, file.path(out, "genotypes", paste0(s, ext)))
  }
}
fwrite(data.table(sample_alias = inc$sample_alias, path = file.path(STAR, inc$sample_alias, paste0(inc$sample_alias, ".SJ.out.tab"))),
       file.path(out, "sj_paths.tsv"), sep = "\t")
fwrite(data.table(sample_alias = inc$sample_alias,
                  path = file.path(STAR, inc$sample_alias, paste0(inc$sample_alias, ".Aligned.sortedByCoord.out.bam"))),
       file.path(out, "bam_paths.tsv"), sep = "\t")

summ <- list(n_libraries = nrow(u), n_included = nrow(inc),
             excluded = as.list(table(factor(u$exclusion_reason[u$exclusion_reason != ""]))),
             q1_fail_parts = list(lt_5M = sum(u$assigned_fragments < 5e6), unique_lt_40 = sum(u$unique_pct < 40),
                                  pc_genes_below_batch_floor = sum(u$pc_genes_ge10 < u$pc_floor)),
             q2_parts = list(mixed_expression = sum(u$mixed_expression), mixed_allele = sum(u$mixed_allele)),
             identity_not_screenable_included = sum(!inc$identity_screenable),
             same_individual_clusters_gt1 = uniqueN(u[cluster_size > 1]$identity_cluster),
             libraries_matching_existing = sum(u$matches_existing_cohort),
             inferred_sex_included = as.list(table(inc$inferred_sex)),
             batches_included = as.list(table(inc$batch)),
             libraries_with_multiple_batch_codes = sum(u$n_batch_codes > 1),
             release_genes_absent = length(setdiff(rel_genes, rownames(counts))),
             strand = "r8 verified fragment counter -s 2 (read 1 antisense after per-run mate order; PRESPEC_SS amendment A1)", stage_read = FALSE)
writeLines(toJSON(summ, auto_unbox = TRUE, pretty = TRUE), file.path(out, "S2_summary.json"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat(toJSON(summ, auto_unbox = TRUE, pretty = TRUE), "\n")
