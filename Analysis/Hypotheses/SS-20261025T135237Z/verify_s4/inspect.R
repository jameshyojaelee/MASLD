d <- readRDS("out_s2_qc/dge_included.rds"); print(class(d)); print(dim(d)); str(d, max.level=1)
if (!is.null(d$genes)) print(head(d$genes)); print(head(colnames(d))); print(head(d$samples))
