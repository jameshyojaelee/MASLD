# DEBUG ONLY: fake exposure. Pilot values permuted across all 224 slides so no real exposure-expression link exists.
suppressMessages(library(data.table)); a <- commandArgs(TRUE)
RR <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/"
o <- fread(paste0(RR, "jobs/3111/order.tsv")); p <- fread(paste0(RR, "jobs/3111/morph_fat.tsv"))[status == "ok"]
set.seed(1); k <- sample.int(nrow(p), 224, replace = TRUE)
fk <- data.table(donor = o$donor, uid = o$uid, status = "ok", mpp = 1.977, tissue_px = sample(p$tissue_px[k]), white_frac = sample(p$white_frac[k]),
                 fat_frac = sample(p$fat_frac[k]), n_droplets = 0, secs = 0)
fwrite(fk, a[1], sep = "\t")
