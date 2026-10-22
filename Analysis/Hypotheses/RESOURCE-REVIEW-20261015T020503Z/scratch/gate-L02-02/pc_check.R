suppressMessages({library(edgeR); library(data.table)})
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
h3 <- file.path(root, "Analysis/MASLD_Model_Benchmark/executions/gse267145-matrix-audit-21064333/raw/GSE267119_H3K27ac.txt.gz")
join <- fread(file.path(root, "Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"))
hdr <- scan(h3, what = "", nlines = 1, quiet = TRUE)
dt <- fread(cmd = paste("zcat", h3), skip = 1, header = FALSE)
cnt <- as.matrix(dt[, -1, with = FALSE]); rownames(cnt) <- gsub('"', "", dt[[1]])
pid <- sub("_.*", "", hdr)
st <- join$stage[match(pid, join$participant_id)]; sex <- join$sex[match(pid, join$participant_id)]
lc <- cpm(calcNormFactors(DGEList(cnt), method = "TMM"), log = TRUE, prior.count = 2)
for (r in c("chr17:50200049-50200840", "chr17:50210153-50210663")) {
  v <- lc[r, ]; d <- st != "NOR"; f <- sex == "F"
  cat(r, ": disease(75)-NOR(24)", round(mean(v[d]) - mean(v[!d]), 3),
      "| women disease(61)-NOR(24)", round(mean(v[d & f]) - mean(v[!d & f]), 3),
      "| F23(13)-NOR(24)", round(mean(v[st == "NASH_F23"]) - mean(v[st == "NOR"]), 3),
      "| Welch p disease vs NOR", signif(t.test(v[d], v[!d])$p.value, 3), "\n")
}
