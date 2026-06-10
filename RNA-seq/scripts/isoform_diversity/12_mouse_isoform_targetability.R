#!/usr/bin/env Rscript
# Bridge to Cas13: map the GSE213621 MANE-switching genes to their MOUSE ortholog's
# isoform structure (mouse liver expression), to flag which have a designable
# isoform-SELECTIVE window vs only pan-isoform targeting. Human DTU = gene-level prior;
# actual targetability comes from mouse-empirical isoform structure (no isoform orthology).
suppressPackageStartupMessages({ library(data.table) })
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
HRES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/human")
MRES <- file.path(PROJ, "RNA-seq/results/isoform_diversity/mouse")
strip <- function(x) sub("\\..*$","",x)

# 67 human MANE-switchers (GSE213621)
tg <- fread(file.path(HRES, "gse213621_dtu_targets.tsv"))[switch_away_from_mane == TRUE]
tg[, h_ensg := strip(gene_id)]
cat(sprintf("[bridge] %d human MANE-switching targets\n", nrow(tg)))

# human ENSG -> mouse ENSMUSG (tier-H; flag 1:1)
pairs <- fread(cmd = paste0("zcat ", PROJ, "/data/external/orthologs/master_ortholog_table.tsv.gz"))
hc <- grep("^tier_H", names(pairs), value = TRUE)
pairs[, isH := Reduce(`|`, lapply(.SD, function(x) toupper(as.character(x)) %in% c("TRUE","T","1","YES"))), .SDcols = hc]
ph <- pairs[isH == TRUE, .(h = strip(human_ensembl), m = strip(mouse_ensembl), msym = mouse_symbol)]
ph <- ph[h != "" & m != ""]
n_h_per_m <- ph[, .N, by = m]; n_m_per_h <- ph[, .N, by = h]
ph[, one2one := h %in% n_m_per_h[N==1, h] & m %in% n_h_per_m[N==1, m]]

# mouse isoform structure from mouse liver expression (all Control+Disease samples)
mtpm <- readRDS(file.path(MRES, "tx_tpm.rds"))
msamp<- fread(file.path(MRES, "samples.tsv")); setkey(msamp, sample_id); msamp <- msamp[colnames(mtpm)]
keep <- msamp$group_binary %in% c("Control","Disease")
mtpm <- mtpm[, keep]
mt2g <- fread(file.path(PROJ, "RNA-seq/results/isoform_diversity/_index/tx2gene_vM38_primary.tsv.gz"))
mt2g <- mt2g[match(rownames(mtpm), txname)]
mt2g[, mensg := strip(geneid)]
meanTPM <- rowMeans(mtpm)                                  # mean across mouse liver samples

mouse_struct <- function(mensg) {
  idx <- which(mt2g$mensg == mensg)
  if (!length(idx)) return(NULL)
  tp <- meanTPM[idx]; gsum <- sum(tp)
  if (gsum < 1) return(data.table(mouse_n_tx=length(idx), mouse_n_expr_iso=0L,
                                  mouse_dom_iso=NA, mouse_dom_prop=NA, mouse_canon_dominant=NA))
  prop <- tp / gsum
  expr <- prop >= 0.05 & tp >= 1                           # meaningfully expressed isoforms (>=5% of gene)
  canon <- mt2g$is_ensembl_canonical[idx] == 1
  data.table(mouse_n_tx = length(idx), mouse_n_expr_iso = sum(expr),
             mouse_dom_iso = rownames(mtpm)[idx][which.max(prop)],
             mouse_dom_prop = round(max(prop), 3),
             mouse_canon_dominant = if (any(canon)) which.max(prop) == which(canon)[which.max(prop[canon])] else NA)
}

rows <- rbindlist(lapply(seq_len(nrow(tg)), function(i) {
  eg <- tg$h_ensg[i]
  orth <- ph[h == eg]
  if (!nrow(orth)) return(data.table(symbol=tg$symbol[i], h_ensg=eg, human_disease_dprop=tg$disease_dprop[i],
                                      mouse_ensg=NA, mouse_symbol=NA, ortholog="none"))
  o <- orth[order(-one2one)][1]                            # prefer 1:1
  st <- mouse_struct(o$m)
  cbind(data.table(symbol=tg$symbol[i], h_ensg=eg, human_disease_dprop=tg$disease_dprop[i],
                   mouse_ensg=o$m, mouse_symbol=o$msym,
                   ortholog=ifelse(o$one2one,"1to1","many")),
        if (is.null(st)) data.table(mouse_n_tx=NA,mouse_n_expr_iso=NA,mouse_dom_iso=NA,mouse_dom_prop=NA,mouse_canon_dominant=NA) else st)
}), fill = TRUE)

# targetability verdict
rows[, isoform_selective_possible := !is.na(mouse_n_expr_iso) & mouse_n_expr_iso >= 2 & mouse_dom_prop < 0.85]
rows[, targeting := fifelse(is.na(mouse_ensg), "no mouse ortholog",
                    fifelse(is.na(mouse_n_expr_iso) | mouse_n_expr_iso == 0, "mouse gene not expressed",
                    fifelse(isoform_selective_possible, "isoform-selective possible", "pan-only (single dominant)")))]
setorder(rows, -human_disease_dprop)
fwrite(rows, file.path(MRES, "gse213621_targets_mouse_isoform_structure.tsv"), sep = "\t")

cat(sprintf("[bridge] ortholog: %d mapped (%d 1:1) / %d none\n",
            sum(rows$ortholog!="none"), sum(rows$ortholog=="1to1"), sum(rows$ortholog=="none")))
cat("[bridge] targeting verdict:\n"); print(rows[, .N, by=targeting][order(-N)])
cat("[bridge] isoform-selective-possible genes (top by human Δprop):\n")
print(head(rows[isoform_selective_possible==TRUE, .(symbol, mouse_symbol, ortholog, mouse_n_expr_iso, mouse_dom_prop, human_disease_dprop)], 20))
cat("-> gse213621_targets_mouse_isoform_structure.tsv\n")
