#!/usr/bin/env Rscript

# Small hand-calculated fixture protecting the biological-unit/interface table.
fixture <- data.frame(
  gene = paste0("G", 1:8),
  genetic = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
  state = c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
  pp4 = c(0.9, 0.8, 0.2, 0.1, 0.7, 0.3, 0.0, 0.4),
  treat_t = c(4, 1, -3, 0.2, -1, 2, 0, -0.5)
)

a <- sum(fixture$genetic & fixture$state)
b <- sum(fixture$genetic & !fixture$state)
c <- sum(!fixture$genetic & fixture$state)
d <- sum(!fixture$genetic & !fixture$state)
stopifnot(identical(c(a, b, c, d), c(1L, 2L, 2L, 3L)))
ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2L, byrow = TRUE))
stopifnot(is.finite(unname(ft$estimate)), ft$p.value >= 0, ft$p.value <= 1)
rho_signed <- suppressWarnings(cor(fixture$pp4, fixture$treat_t, method = "spearman"))
rho_absolute <- suppressWarnings(cor(fixture$pp4, abs(fixture$treat_t), method = "spearman"))
stopifnot(is.finite(rho_signed), is.finite(rho_absolute))
cat("PASS: hand-calculated power fixture\n")

