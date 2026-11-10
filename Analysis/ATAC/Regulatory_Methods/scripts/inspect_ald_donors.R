args <- commandArgs(trailingOnly=TRUE)
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
# Serialized metadata slots can be inspected without loading the analysis package.
x <- readRDS(gzfile(args[1]))
cat('Class:', class(x), '\n')
cat('Assays:', names(x@assays), '\n')
m <- x@meta.data
cat('Metadata columns:', colnames(m), '\n')
for (name in colnames(m)) {
  v <- m[[name]]
  if (length(unique(v)) <= 50) {
    cat('\nFIELD', name, '\n'); print(table(v, useNA='ifany'))
  }
}
cat('\nCell identity examples:', head(rownames(m)), '\n')
write.table(m, args[2], sep='\t', quote=FALSE, col.names=NA)
print(sessionInfo())
