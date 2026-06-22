#!/usr/bin/env Rscript
# Precompute the fig5a "disrupted master-regulator targets" column.
# For each gene: how many DISEASE master-regulator TFs that are ALSO motif-disrupted
# by fine-mapped MASLD risk variants (the fig4 disease_master_regulators panel)
# transcriptionally regulate it, per the curated CollecTRI network. This carries
# fig4's variant->motif master-regulator story down to the convergence targets and
# is far denser than per-gene motif disruption (which only fills TF genes).
suppressPackageStartupMessages({ library(data.table) })
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT  <- file.path(BASE, "RNA-seq/results/multi_evidence/disrupted_mr_targets.csv")

dmr <- fread(file.path(BASE, "Analysis/ATAC/Human_Multiome/scenic_plus/disease_master_regulators.csv"))
mds <- fread(file.path(BASE, "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"))
mds_tfs <- unique(unlist(strsplit(mds$tf_name, "[;:]+")))          # TFs whose motif a risk variant disrupts
disrupted_mr <- dmr[tf_name %in% mds_tfs, tf_name]                 # disease master-reg AND motif-disrupted
cat(sprintf("disrupted disease master-regulators: %d  (COLOC-anchored: %s)\n",
            length(disrupted_mr), paste(dmr[is_coloc==TRUE & tf_name %in% disrupted_mr, tf_name], collapse=",")))

net <- as.data.table(decoupleR::get_collectri(organism = "human", split_complexes = FALSE))
setnames(net, names(net)[1:2], c("source", "target"))
sub <- net[source %in% disrupted_mr]
per_gene <- sub[, .(n_disrupted_mr_targets = uniqueN(source),
                    disrupted_mr_regulators = paste(sort(unique(source)), collapse = ";")),
                by = .(human_symbol = target)]
fwrite(per_gene, OUT)
cat(sprintf("wrote %s: %d genes are targets of >=1 disrupted master-regulator\n", OUT, nrow(per_gene)))
cat("range of n_disrupted_mr_targets:", paste(range(per_gene$n_disrupted_mr_targets), collapse="-"), "\n")
