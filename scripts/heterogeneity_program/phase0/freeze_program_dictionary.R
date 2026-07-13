#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Phase 0f / T3 — freeze the program dictionary (program → gene SET) for PROCON.
#
# Union of pre-annotated, label-free programs: cNMF k16, Hotspot per-cell-type
# soft modules, bulk-NMF k6, WGCNA cross-species-preserved modules, Hallmark.
# Genes harmonised to human_symbol. Programs with Jaccard > 0.5 collapsed
# (keep the larger). Effective number of independent programs from the overlap
# eigenstructure (drives across-program multiplicity in PROCON). cameraPR is
# set-based, so membership (not weights) is what's needed.
# Output: program_dictionary.tsv (summary) + program_gene_membership.tsv (long).
# Lightweight-ish (reads + NMF basis) — run on a compute node via srun.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(data.table) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SC  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2")
OUT <- file.path(BASE, "RNA-seq/results/heterogeneity_program/phase0")
INT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
MIN_GENES <- 10L

# ENSG → symbol map (for bulk-NMF rownames)
deg <- fread(file.path(INT, "canonical_deg_results.csv"))
deg[, gene_base := sub("\\..*$", "", gene)]
ens2sym <- setNames(deg$symbol, deg$gene_base)
to_sym <- function(g) { gb <- sub("\\..*$", "", g); s <- ens2sym[gb]; ifelse(is.na(s) | s=="", g, s) }

mem <- list()   # program_id → character vector of symbols
add <- function(id, genes) { genes <- unique(genes[!is.na(genes) & genes != ""]); if (length(genes) >= MIN_GENES) mem[[id]] <<- genes }

# ── cNMF k16 ──────────────────────────────────────────────────────────────────
cn <- fread(file.path(SC, "mcp/cnmf_annot/global/program_topgenes.k16.tsv"))
for (p in unique(cn$program)) add(sprintf("cnmf_P%02d", p), to_sym(cn[program == p, gene_name]))

# ── Hotspot per-cell-type modules ─────────────────────────────────────────────
for (ct in c("hepatocytes","macrophages","fibroblasts","cholangiocytes","tcells")) {
  f <- file.path(SC, "hotspot_modules", ct, "module_genes.tsv")
  if (!file.exists(f)) next
  hm <- fread(f)
  for (m in unique(hm$module)) add(sprintf("hs_%s_M%02d", substr(ct,1,3), m), to_sym(hm[module == m, gene]))
}

# ── bulk-NMF k6 (top genes per program by basis loading) ──────────────────────
try({
  suppressPackageStartupMessages(library(NMF))
  x <- readRDS(file.path(BASE, "RNA-seq/results/subtypes/nmf_shard_clean_k6.rds"))
  W <- basis(x$fit)                       # gene × 6
  for (p in seq_len(ncol(W))) {
    top <- names(sort(W[, p], decreasing = TRUE))[1:120]
    add(sprintf("bnmf_P%d", p), to_sym(top))
  }
}, silent = TRUE)

# ── WGCNA modules ─────────────────────────────────────────────────────────────
try({
  wg <- fread(file.path(BASE, "Analysis/Cross_Species_Concordance/results/wgcna_module_assignments.csv"))
  setnames(wg, 1:2, c("gene","module"))
  for (m in unique(wg$module)) if (m != 0 && m != "grey") add(sprintf("wgcna_M%s", m), to_sym(wg[module == m, gene]))
}, silent = TRUE)

# ── Hallmark (msigdbr if available) ───────────────────────────────────────────
try({
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    h <- as.data.table(msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
    for (g in unique(h$gs_name)) add(g, h[gs_name == g, gene_symbol])
  }
}, silent = TRUE)

cat(sprintf("Assembled %d raw programs (sources: %s)\n", length(mem),
            paste(unique(sub("_.*","",names(mem))), collapse=", ")))

# ── Collapse Jaccard > 0.5 (keep larger) ──────────────────────────────────────
ids <- names(mem); n <- length(ids); sizes <- lengths(mem)
ord <- order(-sizes); ids <- ids[ord]; keep <- logical(length(ids)); names(keep) <- ids
kept <- character(0)
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
for (id in ids) {
  redundant <- FALSE
  for (k in kept) if (jacc(mem[[id]], mem[[k]]) > 0.5) { redundant <- TRUE; break }
  if (!redundant) kept <- c(kept, id)
}
cat(sprintf("After Jaccard>0.5 collapse: %d programs (from %d)\n", length(kept), length(mem)))

# ── Effective number of independent programs (overlap eigenstructure) ─────────
allg <- sort(unique(unlist(mem[kept])))
M <- sapply(kept, function(id) as.integer(allg %in% mem[[id]]))   # genes × programs
C <- suppressWarnings(cor(M)); C[is.na(C)] <- 0
lam <- pmax(eigen(C, symmetric = TRUE, only.values = TRUE)$values, 0)
keff <- sum(lam)^2 / sum(lam^2)
cat(sprintf("Effective # independent programs (K_eff) = %.1f of %d\n", keff, length(kept)))

# ── Write ─────────────────────────────────────────────────────────────────────
longdt <- rbindlist(lapply(kept, function(id) data.table(program_id = id, gene = mem[[id]])))
longdt[, source := sub("_.*", "", program_id)]
fwrite(longdt, file.path(OUT, "program_gene_membership.tsv"), sep = "\t")
summ <- longdt[, .(n_genes = .N, source = source[1]), by = program_id][order(program_id)]
attr(summ, "k_eff") <- keff
fwrite(summ, file.path(OUT, "program_dictionary.tsv"), sep = "\t")
fwrite(data.table(n_programs = length(kept), n_genes_total = length(allg), k_eff = round(keff,1)),
       file.path(OUT, "program_dictionary_meta.tsv"), sep = "\t")
cat(sprintf("\nWrote program_dictionary.tsv (%d programs), program_gene_membership.tsv (%d rows), %d unique genes\n",
            length(kept), nrow(longdt), length(allg)))
cat("Programs by source:\n"); print(summ[, .N, by = source][order(-N)])
