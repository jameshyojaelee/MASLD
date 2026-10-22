# Machado 2009 dichromat simulation (colorspace::deutan/protan, severity 1) and WCAG 2.1 contrast for Fig 2 colours.
suppressPackageStartupMessages(library(colorspace))
sets <- list(
  B_consequence = c(noncoding = "#1565C0", coding_this = "#C9265E", coding_other = "#7B1FA2"),
  C_classes = c(coding_this = "#C9265E", coding_other = "#E8A3BC", UTR = "#4DB6AC", promoter = "#1565C0", intron = "#64B5F6", intergenic = "#B0BEC5"),
  D_classes = c(protein_alt = "#C9265E", splice = "#8E3B76", syn_utr = "#4DB6AC", noncoding = "#B0BEC5"),
  E_blocks = c(EUR_only = "#D6DBDE", shared = "#A7B6BE", subthreshold = "#E9ECEE", GWS = "#F79268"),
  F_states = c(multi = "#1565C0", single_only = "#64B5F6", le_0.5 = "#B0BEC5"),
  A_ancestry = c(EUR = "#4C72B0", AFR = "#55A868", EAS = "#C44E52", AMR = "#DD8452", SAS = "#8172B3"),
  A_bands = c(gt0.9 = "#08519C", b0.8 = "#4292C6", b0.5 = "#9ECAE1"),
  GH_marks = c(CS = "#1565C0", eQTL = "#E07B39", lead = "#C2185B", coloc = "#FFD600", gwas_r2_high = "#1B2F5B", gwas_r2_low = "#E6E6E5"),
  GH_genes = c(up_sat = "#C9265E", up_mid = "#F5A0BD", down_sat = "#1565C0", down_mid = "#9DD0F9", ns = "#9E9E9E"))
lab <- function(h) as(hex2RGB(h), "LAB")@coords
de <- function(a, b) sqrt(sum((lab(a) - lab(b))^2))
rows <- list()
for (s in names(sets)) {
  v <- sets[[s]]
  for (cvd in c("normal", "deutan", "protan")) {
    sim <- switch(cvd, normal = v, deutan = deutan(v, 1), protan = protan(v, 1)); names(sim) <- names(v)
    cmb <- combn(names(v), 2)
    d <- apply(cmb, 2, function(x) de(sim[x[1]], sim[x[2]]))
    i <- which.min(d)
    rows[[length(rows) + 1]] <- data.frame(set = s, vision = cvd, closest_pair = paste(cmb[, i], collapse = "~"), deltaE76 = round(d[i], 1))
  }
}
cv <- do.call(rbind, rows); print(cv, row.names = FALSE)
lum <- function(h) { c <- hex2RGB(h)@coords; c <- ifelse(c <= 0.03928, c / 12.92, ((c + 0.055) / 1.055)^2.4); sum(c * c(0.2126, 0.7152, 0.0722)) }
cr <- function(f, b) { l <- sort(c(lum(f), lum(b)), decreasing = TRUE); (l[1] + 0.05) / (l[2] + 0.05) }
txt <- data.frame(text = c("D value labels", "D legend labels", "D 23.1 on magenta", "G/H annotation", "G/H text", "ICE2 gene", "FABP1 gene", "n.s. genes", "RORA gene", "ANXA2 gene"),
                  fg = c("#333333", "#4D4D4D", "#FFFFFF", "#8D6E00", "#333333", "#9DD0F9", "#F5A0BD", "#9E9E9E", "#1565C0", "#C9265E"),
                  bg = c("#FFFFFF", "#FFFFFF", "#C9265E", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF"))
txt$contrast <- round(mapply(cr, txt$fg, txt$bg), 2); txt$pass_AA_4.5 <- txt$contrast >= 4.5
print(txt, row.names = FALSE)
write.table(cv, "out_geom/cvd_deltaE.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(txt, "out_geom/wcag_contrast.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), "out_geom/cvd_sessionInfo.txt")
