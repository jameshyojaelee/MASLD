m <- read.delim("figures/main/fig4_singlecell_programs/source_tables/current_candidate/fig4e_all117_bulk_axis_map.tsv", check.names=FALSE)
cat("rows", nrow(m), " testable:", paste(names(table(m$testable)), table(m$testable), collapse=" "), "\n")
t <- m[m$testable %in% c(TRUE,"TRUE","True"),]
for (pre in c("", "composition_pc_")) {
  pf <- t[[paste0(pre,"p_value_meta_fibrosis")]]; pn <- t[[paste0(pre,"p_value_meta_nas")]]
  p <- c(pf, pn); ok <- !is.na(p)
  q <- rep(NA, length(p)); q[ok] <- p.adjust(p[ok], "BH")
  nf <- length(pf)
  cat(pre, "n tests", sum(ok), " BH224 q<0.05 fibrosis:", sum(q[1:nf] < 0.05, na.rm=TRUE), " nas:", sum(q[(nf+1):length(q)] < 0.05, na.rm=TRUE), "\n")
  cat(pre, "stored q<0.05 fibrosis:", sum(t[[paste0(pre,"q_value_fibrosis")]] < 0.05, na.rm=TRUE), " nas:", sum(t[[paste0(pre,"q_value_nas")]] < 0.05, na.rm=TRUE), "\n")
  cat(pre, "max abs diff recomputed vs stored q:", max(abs(q - c(t[[paste0(pre,"q_value_fibrosis")]], t[[paste0(pre,"q_value_nas")]])), na.rm=TRUE), "\n")
}
