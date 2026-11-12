#!/usr/bin/env Rscript
# MQTL2 step 04: ADDED control, not prespecified (requested by team lead after submission, 2026-11-11).
# GCKR rs1260326 has no ALT signal (z about 1.4 in Pan-UKBB and MVP ALT), so the GCKR-triglyceride control is not
# applicable to ALT. Pipeline check instead: Pan-UKBB EUR GGT (GRCh37, effect = alt, N 400,751, irnt) vs the
# triglyceride-class lipid species (Ottensmann TAG/DAG) at the locus-17 window, using the GRCh37 window extracts
# written by 02_coloc.R. Same harmonisation and default priors as 02. Also r2 of locus-17 lead with rs1260326
# in 1000 Genomes EUR (data/1kg_eur, 379 people), via PLINK 1.9.
# Usage: Rscript 04_ggt_gckr_control.R <out_dir> <loci_file>
suppressMessages({library(data.table); library(coloc)})
set.seed(20261110)
args <- commandArgs(TRUE); OUT <- args[1]; LOCI_FILE <- args[2]
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
TABIX <- "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/tabix"
GGT <- file.path(ROOT, "GWAS/MR_Data/PanUKBB/biomarkers-30730-GGT-both_sexes-irnt.tsv.bgz")
loci <- fread(LOCI_FILE, colClasses = list(character = "locus_id"))
L <- loci[locus_id == "17"]
g <- fread(cmd = sprintf("%s %s %d:%d-%d | cut -f1-4,19,25,31,37,43", TABIX, GGT, L$chromosome, L$win_start, L$win_end),
           header = FALSE, col.names = c("chr", "pos", "ref", "alt", "af", "beta", "se", "nlp", "lowconf"))
g <- g[lowconf %in% c("false", FALSE) & !is.na(beta) & se > 0]
cat("GGT window SNPs:", nrow(g), "| max -log10p", max(g$nlp, na.rm = TRUE), "at", g$pos[which.max(g$nlp)], "\n")
mf <- fread(file.path(OUT, "metabolite_files.tsv"))[source == "ottensmann" & grepl("^(Triacylglycerol|Diacylglycerol)", trait)]
res <- rbindlist(lapply(seq_len(nrow(mf)), function(k) {
  f <- file.path(OUT, "windows", paste0(mf$accession[k], ".win.tsv.gz"))
  if (!file.exists(f)) return(NULL)
  q <- fread(f)[chr == L$chromosome & pos >= L$win_start & pos <= L$win_end]
  m <- merge(q[, .(pos, q_ea = ea, q_oa = oa, q_beta = beta, q_se = se)],
             g[, .(pos, g_eff = alt, g_oth = ref, g_beta = beta, g_se = se)], by = "pos")
  m <- m[(q_ea == g_eff & q_oa == g_oth) | (q_ea == g_oth & q_oa == g_eff)]
  m[q_ea == g_oth, q_beta := -q_beta]
  m <- m[!((g_eff %in% c("A", "T") & g_oth %in% c("A", "T")) | (g_eff %in% c("C", "G") & g_oth %in% c("C", "G")))]
  m <- m[!(duplicated(pos) | duplicated(pos, fromLast = TRUE))]
  if (nrow(m) < 50) return(NULL)
  s <- as.character(m$pos)
  a <- suppressMessages(coloc.abf(list(beta = m$g_beta, varbeta = m$g_se^2, N = 400751, type = "quant", sdY = 1, snp = s),
                                  list(beta = m$q_beta, varbeta = m$q_se^2, N = mf$N[k], type = "quant", sdY = 1, snp = s)))
  ti <- which.max(a$results$SNP.PP.H4)
  data.table(accession = mf$accession[k], trait = mf$trait[k], n_snps = nrow(m),
             PP.H3 = a$summary[["PP.H3.abf"]], PP.H4 = a$summary[["PP.H4.abf"]], top_pos = a$results$snp[ti])
}))
setorder(res, -PP.H4)
fwrite(res, file.path(OUT, "control_added_GGT_GCKR_triglyceride.tsv"), sep = "\t")
cat("TAG/DAG traits tested:", nrow(res), "| PP.H4 >= 0.8:", sum(res$PP.H4 >= 0.8), "\n"); print(head(res, 8))
# r2 between locus-17 lead (2:27741237) and rs1260326 (2:27730940), 1000G EUR GRCh37
ld <- file.path(OUT, "ld_gckr")
system(sprintf("plink --bfile %s/data/1kg_eur/chr2_eur --ld 2:27741237:T:C 2:27730940:T:C --out %s > /dev/null 2>&1", ROOT, ld))
cat(grep("R-sq|D'", readLines(paste0(ld, ".log")), value = TRUE), sep = "\n")
