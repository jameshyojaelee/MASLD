#!/bin/bash
#SBATCH --job-name=Rprof
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/profsum_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/profsum_%j.err
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
Rscript -e '
f <- "'"$FM"'/logs/coloc_perf/coloc_prof_19671233.out"
s <- summaryRprof(f)
tot <- s$sampling.time
cat(sprintf("total sampled: %.0f s (%.1f h)\n\n", tot, tot/3600))
bs <- head(s$by.self, 20); bs$pct <- round(100*bs$self.time/tot,1)
cat("--- BY SELF (where CPU sits) ---\n"); print(bs[,c("self.time","pct")])
bt <- s$by.total
fns <- c("get_ld_per_locus","fread","bdiag","as.matrix","susie_rss","estimate_s_rss",
         "eigen","coloc.susie","coloc.abf","find_LD_block","runsusie")
fq <- paste0(dQuote(fns, FALSE)); k <- fq[fq %in% rownames(bt)]
b2 <- bt[k,,drop=FALSE]; b2$pct <- round(100*b2$total.time/tot,1)
cat("\n--- BY TOTAL (stage owners) ---\n"); print(b2[order(-b2$total.time),c("total.time","pct","self.time")])
gt <- function(x) { q <- dQuote(x, FALSE); if (q %in% rownames(bt)) bt[q,"total.time"] else 0 }
ld <- gt("get_ld_per_locus"); es <- gt("estimate_s_rss"); sr <- gt("susie_rss")
cl <- gt("coloc.susie")+gt("coloc.abf")
cat("\n--- HEADLINE SPLIT ---\n")
cat(sprintf("  LD load + reshape : %8.0f s  %5.1f%%\n", ld, 100*ld/tot))
cat(sprintf("  estimate_s_rss    : %8.0f s  %5.1f%%\n", es, 100*es/tot))
cat(sprintf("  susie_rss         : %8.0f s  %5.1f%%\n", sr, 100*sr/tot))
cat(sprintf("  coloc.*           : %8.0f s  %5.1f%%\n", cl, 100*cl/tot))
cat(sprintf("  STATISTICS TOTAL  : %8.0f s  %5.1f%%\n", es+sr+cl, 100*(es+sr+cl)/tot))
'
