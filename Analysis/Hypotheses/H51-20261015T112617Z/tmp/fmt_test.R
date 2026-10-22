for (f in c("h51_common.R","stage1_blinded.R","stage2_outcome.R")) invisible(parse(f))
suppressMessages(library(data.table))
l <- readLines("h51_common.R"); i <- grep("^fmt <- ", l); eval(parse(text = l[i:(i + 1)]))
print(fmt(data.table(a = c(23648L, 5L), b = c(0.123456, 2), c = c(TRUE, FALSE), d = c("x", "y"))))
