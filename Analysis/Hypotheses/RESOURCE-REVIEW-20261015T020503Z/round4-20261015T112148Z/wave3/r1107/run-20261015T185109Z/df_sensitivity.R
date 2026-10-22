library(data.table); a <- commandArgs(TRUE)
w <- fread(a[1]); w[, p_df_n2 := 2*pt(-abs(rho*sqrt((n-2)/(1-rho^2))), n-2)]; w[, q_df_n2 := p.adjust(p_df_n2, "BH")]
fwrite(w, a[2], sep = "\t")
