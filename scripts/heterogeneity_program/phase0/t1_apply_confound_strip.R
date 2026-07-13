#!/usr/bin/env Rscript
# Post-hoc confound strip for the T1 DV results (sex/library-prep variance leak fix).
# Removes chrY/chrM/XIST + ribosomal (RP[SL]/MRP[SL]) + IG/TR + hemoglobin genes
# (repo convention; CLAUDE.md NMF strip) from the DV candidate set, re-BH over the
# cleaned universe, and re-derives dv_sig / partition / direction. DV is per-gene,
# so this is equivalent to excluding them from the universe (and conservative: the
# stripped genes' inflated variance was in the pooled null, so retained genes can
# only get LESS significant by keeping the old dv_p_pool). The canonical
# t1_dv_screen.R now also strips at the universe stage for future full re-runs.
suppressPackageStartupMessages({ library(data.table) })
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
P    <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")

r  <- fread(file.path(P, "dv_results.csv"))
gm <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
gms <- gm[match(r$gene_base, ensembl_base)]
hb <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
strip <- (gms$chromosome %in% c("chrY","chrM")) | (gms$gene_name == "XIST") |
         (grepl("^(RP[SL]|MRP[SL])", gms$gene_name) & !grepl("K[ABCL]?[0-9]", gms$gene_name)) |
         grepl("^IG_|^TR_", gms$gene_biotype) | (gms$gene_name %in% hb)
strip[is.na(strip)] <- FALSE
cat(sprintf("confound genes in universe: %d (of which DV-sig pre-strip: %d)\n",
            sum(strip), sum(strip & r$dv_sig == TRUE)))

# backup once, then strip + re-derive
bk <- file.path(P, "dv_results_preStrip.csv")
if (!file.exists(bk)) fwrite(r, bk)
r <- r[!strip]
r[, dv_q_pool := p.adjust(dv_p_pool, "BH")]                          # re-BH over cleaned universe
r[, dv_sig := dv_q_pool < 0.05 & loco_pass == TRUE]
r[, partition := fifelse(dv_sig & !is_mean_deg, "variance_only",
                  fifelse(dv_sig & is_mean_deg, "mean_and_variance",
                   fifelse(!dv_sig & is_mean_deg, "mean_only", "neither")))]
r[, dv_direction := fifelse(dv_sig & dv_logfc > 0, "fan_out",
                     fifelse(dv_sig & dv_logfc <= 0, "canalize", "ns"))]
fwrite(r, file.path(P, "dv_results.csv"))

# regenerate dv_atlas_columns.tsv (identical schema)
atl <- r[, .(human_symbol = symbol, gene_ensembl = gene_base, dv_t, dv_logfc,
             dv_p_perm, dv_q_perm, dv_p_pool, dv_q_pool,
             dv_nb_logratio = nb_logratio, dv_nb_concordant = nb_concordant,
             dv_loco_signfrac = loco_signfrac, dv_loco_pass = loco_pass,
             dv_sig, dv_direction, dv_partition = partition)]
if (!file.exists(file.path(P, "dv_atlas_columns_preStrip.tsv")))
  file.copy(file.path(P, "dv_atlas_columns.tsv"), file.path(P, "dv_atlas_columns_preStrip.tsv"))
fwrite(atl, file.path(P, "dv_atlas_columns.tsv"), sep = "\t")

cat(sprintf("\nPOST-STRIP universe %d genes\n", nrow(r)))
cat("partition:\n"); print(r[, .N, by = partition][order(-N)])
cat(sprintf("DV-sig %d | variance_only %d | fan_out %d | canalize %d\n",
            sum(r$dv_sig), sum(r$partition=="variance_only"),
            sum(r$dv_direction=="fan_out"), sum(r$dv_direction=="canalize")))
cat(sprintf("Jaccard(DV-sig, Tier-1 mean-DEG): %.3f\n",
            { a<-r[dv_sig==TRUE, gene_base]; b<-r[is_mean_deg==TRUE, gene_base]
              length(intersect(a,b))/length(union(a,b)) }))
