suppressPackageStartupMessages(library(data.table))
S <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R02-04-coverage/"
st <- readRDS(paste0(S, "state.rds")); base <- st$base; fm0 <- st$fm0; gene_y <- st$gene_y
pb <- base[plotted == TRUE]
cat("plotted:", nrow(pb), "\n")
pmx <- function(...) suppressWarnings(pmax(..., na.rm = TRUE))
pb[, sc := pmx(susie_pip, carma_pip)]; pb[is.infinite(sc), sc := NA]
pb[, joint := pmx(susiex_pip, mesusie_pip)]
cat("x from SuSiEx/MESuSiE above both SuSiE and CARMA:", sum(!is.na(pb$joint) & pb$causal_pip > fifelse(is.na(pb$sc), -Inf, pb$sc) & pb$joint == pb$causal_pip), "\n")
cat("  of which susiex is the max:", sum(!is.na(pb$susiex_pip) & pb$susiex_pip == pb$causal_pip & pb$causal_pip > fifelse(is.na(pb$sc), -Inf, pb$sc)),
    " mesusie is the max:", sum(!is.na(pb$mesusie_pip) & pb$mesusie_pip == pb$causal_pip & pb$causal_pip > fifelse(is.na(pb$sc), -Inf, pb$sc)), "\n")
cat("x from susie_pip where recommended_pip excluded it (non-EUR or non-converged), susie > recommended:",
    sum(!is.na(pb$susie_pip) & pb$susie_pip == pb$causal_pip & (is.na(pb$recommended_pip) | pb$susie_pip > pb$recommended_pip)), "\n")
print(pb[!is.na(susie_pip) & susie_pip == causal_pip & (is.na(recommended_pip) | susie_pip > recommended_pip), .N, by = .(ancestry, susie_converged)])
cat("x and y from different GWAS:", sum(pb$x_gwas != pb$susie_gwas), "\n")
cat("x and y from different lead SNP:", sum(pb$x_top_snp != pb$y_top_snp), "\n")
cat("x and y differ in GWAS or SNP:", sum(pb$x_gwas != pb$susie_gwas | pb$x_top_snp != pb$y_top_snp), "\n")
corner <- function(d) sum(d$causal_pip >= 0.9 & d$pp4 < 0.5)
cat("baseline corner PIP>=0.9 & PP.H4<0.5 (plotted):", corner(pb), "\n")

# Alternative A: SuSiE-only x, same gene-level max over all leads/GWAS
fm <- copy(fm0); fm[, cp := susie_pip]; fm <- fm[is.finite(cp)]
fm <- fm[order(study, chromosome, position, -cp, source_rank)]
fu <- unique(fm, by = c("study","chromosome","position"))
cm <- merge(st$coloc, fu[, .(study, chromosome, position, cp)], by.x = c("gwas_name","snp_chr","snp_pos"), by.y = c("study","chromosome","position"))
gx <- cm[order(gene, -cp), .SD[1], by = gene][, .(gene, xA = cp)]
A <- merge(base[, .(gene, causal_pip, pp4, plotted)], gx, by = "gene", all.x = TRUE)
roster <- c("HSD17B13","TM6SF2","MBOAT7","GCKR","RORA")
A[, plotA := !is.na(xA) & ((pp4 > 0.05) | (xA >= 0.99) | gene %in% roster)]
Ap <- A[plotted == TRUE]
cat("\nA SuSiE-only x: plotted baseline points with x changed (|dx|>1e-4):", sum(is.na(Ap$xA) | abs(Ap$xA - Ap$causal_pip) > 1e-4),
    " (no SuSiE x at all:", sum(is.na(Ap$xA)), ")\n")
cat("  baseline plotted but no longer pass gate:", sum(A$plotted & !A$plotA), "  newly pass:", sum(!A$plotted & A$plotA), "\n")
cat("  of baseline plotted: x>=0.99 base and xA<0.99 (or NA):", sum(Ap$causal_pip >= 0.99 & (is.na(Ap$xA) | Ap$xA < 0.99)), "\n")
cat("  plotted under A:", sum(A$plotA), "  corner under A:", sum(A$plotA & A$xA >= 0.9 & A$pp4 < 0.5), "\n")

# Alternative B: SuSiE PIP at the same GWAS and lead SNP as y
fmB <- fm0[is.finite(susie_pip), .(sp = max(susie_pip)), by = .(study, chromosome, position)]
gy <- copy(gene_y)[!is.na(pp4_susie)]
gy[, c("c","p") := tstrsplit(y_top_snp, ":", fixed = TRUE, type.convert = TRUE)]
Bt <- merge(gy, fmB, by.x = c("susie_gwas","c","p"), by.y = c("study","chromosome","position"), all.x = TRUE)
Bt[, pp4 := pmax(pp4_susie, 0)]
Bt[, plotB := !is.na(sp) & ((pp4 > 0.05) | (sp >= 0.99) | gene %in% roster)]
cat("\nB same-GWAS same-lead SuSiE x: genes with pp4_susie:", nrow(Bt), " with finite x:", sum(!is.na(Bt$sp)), " plotted:", sum(Bt$plotB),
    "  corner:", sum(Bt$plotB & Bt$sp >= 0.9 & Bt$pp4 < 0.5), "\n")
# also all-method x at same GWAS + lead (to separate the method effect from the GWAS/lead effect)
fmC <- copy(fm0); fmC[, cp := pmx(recommended_pip, susie_pip)]; fmC <- fmC[is.finite(cp), .(cp = max(cp)), by = .(study, chromosome, position)]
Ct <- merge(gy, fmC, by.x = c("susie_gwas","c","p"), by.y = c("study","chromosome","position"), all.x = TRUE)
Ct[, pp4 := pmax(pp4_susie, 0)]
Ct[, plotC := !is.na(cp) & ((pp4 > 0.05) | (cp >= 0.99) | gene %in% roster)]
cat("C same-GWAS same-lead producer x: plotted:", sum(Ct$plotC), " corner:", sum(Ct$plotC & Ct$cp >= 0.9 & Ct$pp4 < 0.5), "\n")
cat("\nRoster:\n")
r <- merge(base[gene %in% roster, .(gene, base_x = causal_pip, x_gwas, x_top_snp, y_gwas = susie_gwas, y_top_snp, pp4)],
           A[gene %in% roster, .(gene, xA)], by = "gene")
r <- merge(r, Bt[gene %in% roster, .(gene, xB = sp)], by = "gene", all.x = TRUE)
r <- merge(r, Ct[gene %in% roster, .(gene, xC = cp)], by = "gene", all.x = TRUE)
print(r)
