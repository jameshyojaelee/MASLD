#!/usr/bin/env Rscript
# =============================================================================
# 95d_program_pathway_naming.R — Pathway-driven naming of the k=6 bulk NMF programs
#
# Scores each NMF program against a liver-curated + MSigDB Hallmark gene-set
# panel by fgsea on the (per-program centered) gene loadings, writes an audit
# atlas of enrichments, proposes concise pathway-derived labels, and (DEFAULT)
# writes the CANONICAL program_labels.csv + remaps nmf_assignments.csv from
# the human-approved CURATED_LABELS below.
#
# Pipeline order: 44 -> 95 -> 95c -> 95d (95d is the FINAL labeling authority).
#
# Emits:
#   program_pathway_enrichment.csv       — program_code, set_name, source_db, NES, pval, padj, size, leadingEdge
#   liver_curated_genesets.tsv           — set_name, gene, source (panel audit)
#   program_labels_proposed.csv          — program_code, top_pathway, evidence, proposed_label (pre-filled)
#   program_labels.csv                   — CANONICAL labels (rewritten, default ON)
#   program_rename_mapping.csv           — program_code, old_label, new_label (audit trail)
#   program_labels_prepathway_backup.csv — one-time backup of pre-rename labels
#   nmf_assignments.csv                  — dominant_program remapped (backup first)
#   nmf_assignments_prepathway_backup.csv — one-time backup of pre-rename assignments
#
# Escape hatch (dry-run / re-curation): set NMF_LABEL_DRYRUN=TRUE to skip all
# canonical writes and only emit the enrichment atlas + proposed file.
#
# Run on a compute node (rnaseq env), NEVER the login node:
#   srun --partition=cpu --qos=interactive --mem=32G --cpus-per-task=4 --time=48:00:00 \
#     Rscript RNA-seq/95d_program_pathway_naming.R
# =============================================================================

# ---------------------------------------------------------------------------
# Human-approved pathway-derived program names (2026-06-24). These are the
# AUTHORITATIVE labels; 95d is the final labeling step (44 -> 95 -> 95c -> 95d).
# P2/P4/P5 are curatorial calls (not the literal top enriched set — see evidence
# table below); re-curate only if the NMF cache changes (new gene loadings).
#
# Evidence table (fgsea on per-program centered loadings, k=6 canonical cache):
#   P1 Inflammatory-EMT  : HALLMARK_EMT(NES=2.80) + INFLAMMATORY_RESPONSE(NES=2.73)
#   P2 Hepatocyte-Metabolic: top set IFN-γ(NES=1.96,weak); named for bile-acid+
#                            xenobiotic+complement context; hepatocyte-metabolic call
#   P3 Fibrotic-ECM      : EMT(NES=2.57) + LIVER_FIBROTIC_ECM(NES=2.31)
#   P4 Unresolved        : ns (no set padj<0.05) — curatorial
#   P5 Noncoding         : ns (lncRNA-dominated top genes) — curatorial
#   P6 Skeletal-muscle   : HALLMARK_MYOGENESIS(NES=3.20) + LIVER_SKELETAL_MUSCLE_CONTAM(NES=2.85)
# ---------------------------------------------------------------------------
CURATED_LABELS <- c(
  P1 = "Inflammatory-EMT",
  P2 = "Hepatocyte-Metabolic",
  P3 = "Fibrotic-ECM",
  P4 = "Unresolved",
  P5 = "Noncoding",
  P6 = "Skeletal-muscle"
)

suppressPackageStartupMessages({
  library(data.table)
  library(NMF)
  library(fgsea)
  library(msigdbr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

in_dir     <- file.path(BASE, "RNA-seq/results/subtypes")
out_dir    <- Sys.getenv("SUBTYPE_OUT_DIR", unset = in_dir)
atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

chosen_k <- as.integer(Sys.getenv("NMF_CHOSEN_K", unset = "6"))

# ---- Load NMF cache (same precedence as Script 44) ----
env_cache   <- Sys.getenv("NMF_CACHE_PATH", "")
clean_cache <- file.path(in_dir, "nmf_results_cache_clean.rds")
cache_path  <- if (nzchar(env_cache)) env_cache else clean_cache
if (!file.exists(cache_path)) stop("No NMF cache at ", cache_path)
cat("Loading NMF cache:", basename(cache_path), "\n")
cached <- readRDS(cache_path)
res_k  <- cached$nmf_results[[as.character(chosen_k)]]
if (is.null(res_k)) stop("k=", chosen_k, " not in cache keys: ",
                         paste(names(cached$nmf_results), collapse = ","))
W <- basis(res_k)                        # genes x k
rownames(W) <- rownames(cached$mat_nn)   # Ensembl IDs (versioned)

# ---- Ensembl -> symbol map (same as Script 44) ----
atlas     <- fread(atlas_path, select = c("ensembl_id", "human_symbol"))
id2sym    <- setNames(atlas$human_symbol, atlas$ensembl_id)
strip_ver <- function(x) sub("\\.\\d+$", "", x)
sym_W     <- id2sym[strip_ver(rownames(W))]; names(sym_W) <- rownames(W)

# ---- Per-program centering (same as Script 44) ----
W_centered <- sapply(seq_len(chosen_k), function(j) W[, j] - rowMeans(W[, -j, drop = FALSE]))
rownames(W_centered) <- rownames(W)

# ============================================================
# Reference panel: MSigDB Hallmark + curated liver sets
# ============================================================
# CANON sets copied verbatim from 44_molecular_subtyping.R:69-96 (provenance).
FIBROTIC_CANON <- c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A3","COL15A1",
  "LUM","DCN","FBN1","LOX","SPARC","ELN","FN1","ACTA2","DES","MYOCD","CNN1",
  "TAGLN","MYL9","MYH11","FLNC","PDGFRB","THY1","VIM","CDH11","LTBP2","CTSK",
  "TIMP1","PDPN","SMOC2","EDIL3","POSTN","RCN3","BGN","LRRC15","SFRP4","PTPRQ")
HEPATOCYTE_CANON <- c("ALB","APOB","APOA1","APOA2","TF","SERPINA1","TTR","HP",
  "HRG","HPX","AHSG","CP","C3","FGA","FGB","FGG","PCK1","G6PC","HNF4A","HNF1A",
  "CPS1","FABP1")
NEURAL_CANON <- c("NCAM1","NCAM2","ASCL1","NLGN4X","NLGN4Y","CNTN2","STMN2",
  "NRXN1","NEFL","GAP43","SYN1","GABRG2","RBFOX1","RBFOX3","MAP2","TUBB3")
METABOLIC_CANON <- c("CYP3A4","CYP2E1","ALB","APOB","APOA1","G6PC","PCK1",
  "HMGCS2","CPT1A","ACOX1","SLC10A1","ABCB11","FABP1","ALDOB","GLUD1","TAT",
  "SERPINA1","CPS1","HNF4A")
HCC_CANON <- c("AKR1B10","GPC3","AFP","HKDC1","IGF2BP3","SPINK1")
CELLCYCLE_CANON <- c("MKI67","TOP2A","PCNA","CCNB1","CDK1","AURKB","BIRC5")

# New curated liver cell-identity / zonation / contamination sets.
liver_curated <- list(
  LIVER_HEPATOCYTE_IDENTITY            = HEPATOCYTE_CANON,
  LIVER_HEPATOCYTE_METABOLIC           = METABOLIC_CANON,
  LIVER_HEPATIC_STELLATE_MYOFIBROBLAST = c("ACTA2","COL1A1","COL1A2","COL3A1",
    "PDGFRB","PDGFRA","TAGLN","MYH11","LUM","DCN","TIMP1","LOX","THY1","POSTN","RGS5"),
  LIVER_FIBROTIC_ECM                   = FIBROTIC_CANON,
  LIVER_CHOLANGIOCYTE                  = c("KRT19","KRT7","SOX9","EPCAM","CFTR",
    "HNF1B","ANXA4","SPP1","MUC1","CLDN4","CLDN7"),
  LIVER_ENDOTHELIAL_LSEC               = c("PECAM1","CLEC4G","STAB2","OIT3","LYVE1",
    "FCGR2B","VWF","CD34","ENG","KDR","FLT1"),
  LIVER_KUPFFER_MACROPHAGE             = c("CD68","CD163","MARCO","VSIG4","TIMD4",
    "CD5L","TREM2","CSF1R","MRC1","LYZ","C1QA","C1QB"),
  LIVER_ZONATION_PERIPORTAL            = c("ARG1","ASS1","PCK1","CPS1","SDS","HAL",
    "AGXT","GLS2","ALDOB"),
  LIVER_ZONATION_PERICENTRAL           = c("GLUL","CYP2E1","CYP1A2","OAT","ADH4",
    "CYP1A1","SLC1A2","RHBG","GSTA1"),
  LIVER_NEURAL_CONTAM                  = NEURAL_CANON,
  LIVER_SKELETAL_MUSCLE_CONTAM         = c("MYH1","MYH2","MYH7","CKM","MYL2","ACTN2",
    "TNNT3","TNNT1","MYBPC1","MYBPC2","TRDN","XIRP2","NEB","TTN","DES","MYBPH"),
  LIVER_CELLCYCLE_PROLIFERATION        = CELLCYCLE_CANON,
  LIVER_HCC_TUMOR                      = HCC_CANON
)

# Hallmark
hall      <- msigdbr(species = "Homo sapiens", collection = "H")
hall_list <- split(hall$gene_symbol, hall$gs_name)

panel  <- c(hall_list, liver_curated)
src_db <- c(setNames(rep("hallmark", length(hall_list)), names(hall_list)),
            setNames(rep("liver_curated", length(liver_curated)), names(liver_curated)))

# Panel sanity check (closest thing to a unit test here)
stopifnot(all(lengths(panel) > 0))
stopifnot(length(liver_curated) == 13L)
cat(sprintf("Panel: %d Hallmark + %d liver_curated = %d sets\n",
            length(hall_list), length(liver_curated), length(panel)))

# Emit panel audit TSV
panel_tsv <- rbindlist(lapply(names(panel), function(s)
  data.table(set_name = s, gene = panel[[s]], source = src_db[s])))
fwrite(panel_tsv, file.path(out_dir, "liver_curated_genesets.tsv"), sep = "\t")
cat("Wrote liver_curated_genesets.tsv (", nrow(panel_tsv), " rows)\n", sep = "")

# ============================================================
# fgsea per program against the combined panel
# ============================================================
# minSize lowered 10 -> 3 vs Script 44 to admit the deliberately small curated
# marker panels (e.g. HCC n=6, zonation n~9). This only ADDS sets relative to a
# Hallmark-only run; large Hallmark sets are unaffected. 'size' is retained so
# tiny-set NES can be read as supportive, not definitive.
all_fg <- list()
for (j in seq_len(chosen_k)) {
  w_j <- W_centered[, j]
  has_sym <- !is.na(sym_W) & sym_W != ""
  ranked <- w_j[has_sym]; names(ranked) <- sym_W[has_sym]
  ranked <- ranked[!duplicated(names(ranked))]
  ranked <- sort(ranked, decreasing = TRUE)
  set.seed(42)
  fg <- suppressWarnings(fgsea(pathways = panel, stats = ranked,
                               minSize = 3, maxSize = 500,
                               nPermSimple = 10000, scoreType = "std"))
  fg$program_code <- paste0("P", j)
  fg$source_db    <- src_db[fg$pathway]
  fg$leadingEdge  <- vapply(fg$leadingEdge, paste, character(1), collapse = ";")
  all_fg[[j]] <- as.data.table(fg)
}
enrich <- rbindlist(all_fg)
setnames(enrich, "pathway", "set_name")
setcolorder(enrich, c("program_code","set_name","source_db","NES","pval","padj","size","leadingEdge"))
enrich <- enrich[order(program_code, -NES)]
fwrite(enrich, file.path(out_dir, "program_pathway_enrichment.csv"))
cat("Wrote program_pathway_enrichment.csv (", nrow(enrich), " rows)\n", sep = "")

# ============================================================
# Propose labels: top significant set per program (by NES)
# ============================================================
sig <- enrich[padj < 0.05 & NES > 0][order(program_code, -NES)]
prop_tbl <- sig[, .(
  top_pathway = sprintf("%s(NES=%.2f,padj=%.1e)", set_name[1], NES[1], padj[1]),
  evidence    = paste(sprintf("%s(NES=%.2f)", set_name[seq_len(min(5, .N))],
                              NES[seq_len(min(5, .N))]), collapse = "; ")
), by = program_code]
allp <- data.table(program_code = paste0("P", seq_len(chosen_k)))
prop_tbl <- merge(allp, prop_tbl, by = "program_code", all.x = TRUE)
prop_tbl[is.na(top_pathway),
         `:=`(top_pathway = "ns", evidence = "no set padj<0.05 — name from top genes")]
# Pre-fill proposed_label from CURATED_LABELS so a fresh run never loses curation.
prop_tbl[, proposed_label := unname(CURATED_LABELS[program_code])]
fwrite(prop_tbl, file.path(out_dir, "program_labels_proposed.csv"))

cat("\n=== TOP 5 SETS PER PROGRAM (padj<0.05, NES>0) ===\n")
for (p in paste0("P", seq_len(chosen_k))) {
  cat("\n", p, ":\n", sep = "")
  sub <- enrich[program_code == p & padj < 0.05 & NES > 0][order(-NES)][seq_len(min(5, .N))]
  if (nrow(sub) == 0 || is.na(sub$set_name[1])) cat("  (none significant)\n")
  else for (i in seq_len(nrow(sub)))
    cat(sprintf("  %-45s NES=%.2f padj=%.1e [%s]\n",
        sub$set_name[i], sub$NES[i], sub$padj[i], sub$source_db[i]))
}
# ============================================================
# Canonical-write block (DEFAULT ON; set NMF_LABEL_DRYRUN=TRUE to skip)
# ============================================================
dry_run <- toupper(Sys.getenv("NMF_LABEL_DRYRUN", "FALSE")) == "TRUE"

if (dry_run) {
  cat("\n[DRY-RUN] NMF_LABEL_DRYRUN=TRUE — skipping all canonical writes.\n")
} else {
  cat("\n=== CANONICAL LABEL WRITE (NMF_LABEL_DRYRUN=FALSE) ===\n")

  labels_path      <- file.path(out_dir, "program_labels.csv")
  labels_bak_path  <- file.path(out_dir, "program_labels_prepathway_backup.csv")
  assign_path      <- file.path(out_dir, "nmf_assignments.csv")
  assign_bak_path  <- file.path(out_dir, "nmf_assignments_prepathway_backup.csv")
  rename_path      <- file.path(out_dir, "program_rename_mapping.csv")

  # -- Require program_labels.csv (written by Script 44 upstream) --
  if (!file.exists(labels_path)) {
    stop("program_labels.csv not found at ", labels_path,
         "\n  Run Script 44 first to generate the canonical label file.")
  }

  # -- Load existing labels (skip comment lines starting with #) --
  existing_labels <- fread(cmd = paste0("grep -v '^#' ", shQuote(labels_path)))

  # -- Validate: all program_codes present in CURATED_LABELS --
  missing_codes <- setdiff(existing_labels$program_code, names(CURATED_LABELS))
  if (length(missing_codes) > 0) {
    stop("program_labels.csv contains program codes not in CURATED_LABELS: ",
         paste(missing_codes, collapse = ", "),
         "\n  The NMF cache may have changed — re-curate CURATED_LABELS in 95d.")
  }

  # -- Drift guard: build and print reconciliation table --
  # top_sig_set per program from freshly computed enrichment
  top_sig <- enrich[padj < 0.05 & NES > 0][order(program_code, -NES)]
  top_sig_per_prog <- top_sig[, .(current_top_sig_set = if (.N > 0) set_name[1] else "ns",
                                   current_top_NES     = if (.N > 0) NES[1] else NA_real_),
                               by = program_code]
  recon <- merge(
    data.table(program_code = names(CURATED_LABELS),
               assigned_curated_label = unname(CURATED_LABELS)),
    top_sig_per_prog, by = "program_code", all.x = TRUE
  )
  recon[is.na(current_top_sig_set), current_top_sig_set := "ns"]

  cat("\n--- DRIFT RECONCILIATION TABLE ---\n")
  cat(sprintf("%-6s %-22s %-50s %s\n",
              "Code", "Curated Label", "Current Top Sig Set", "NES"))
  cat(strrep("-", 95), "\n")
  for (i in seq_len(nrow(recon))) {
    cat(sprintf("%-6s %-22s %-50s %s\n",
                recon$program_code[i],
                recon$assigned_curated_label[i],
                recon$current_top_sig_set[i],
                ifelse(is.na(recon$current_top_NES[i]), "NA",
                       sprintf("%.2f", recon$current_top_NES[i]))))
  }
  cat(strrep("-", 95), "\n")

  # Drift warnings
  pathway_labeled <- c("P1","P2","P3","P6")
  curatorial_ns   <- c("P4","P5")
  for (i in seq_len(nrow(recon))) {
    code   <- recon$program_code[i]
    top_s  <- recon$current_top_sig_set[i]
    top_n  <- recon$current_top_NES[i]
    has_strong_sig <- !is.na(top_n) && top_n >= 2.0 && top_s != "ns"
    if (code %in% pathway_labeled && top_s == "ns") {
      warning(sprintf(
        "DRIFT: %s assigned '%s' (pathway-based label) but now shows NO significant set. Curation may be stale.",
        code, recon$assigned_curated_label[i]))
    }
    if (code %in% curatorial_ns && has_strong_sig) {
      warning(sprintf(
        "DRIFT: %s assigned '%s' (curatorial no-sig label) but now has strong set %s (NES=%.2f). Curation may be stale.",
        code, recon$assigned_curated_label[i], top_s, top_n))
    }
  }

  # -- Backup program_labels.csv ONLY if backup does not already exist --
  # NOTE: this is a ONE-TIME pre-rename snapshot. On a from-scratch refit
  # (44 -> 95 -> 95c -> 95d), Script 44 rewrites program_labels.csv with fresh
  # old-vocabulary labels but does NOT touch this backup; an existing backup is
  # preserved (not re-snapshotted). To capture a new pre-rename snapshot, delete
  # *_prepathway_backup.csv before re-running.
  if (!file.exists(labels_bak_path)) {
    file.copy(labels_path, labels_bak_path)
    cat("Created backup: program_labels_prepathway_backup.csv\n")
  } else {
    cat("Backup already exists (skipping): program_labels_prepathway_backup.csv\n")
  }

  # -- Build and write rename mapping (old_label = current file's biological_label) --
  rename_map <- data.table(
    program_code = existing_labels$program_code,
    old_label    = existing_labels$biological_label,
    new_label    = unname(CURATED_LABELS[existing_labels$program_code])
  )
  fwrite(rename_map, rename_path)
  cat("Wrote program_rename_mapping.csv\n")

  # -- top_pathway per program from freshly computed enrichment --
  top_path_per_prog <- top_sig[, .(
    top_pathway_str = if (.N > 0)
      sprintf("%s(NES=%.2f,padj=%.1e)", set_name[1], NES[1], padj[1])
    else "ns"
  ), by = program_code]
  # ensure all 6 programs represented
  all_progs <- data.table(program_code = paste0("P", seq_len(chosen_k)))
  top_path_per_prog <- merge(all_progs, top_path_per_prog, by = "program_code", all.x = TRUE)
  top_path_per_prog[is.na(top_pathway_str), top_pathway_str := "ns"]

  # -- Rewrite program_labels.csv (preserve schema + top_genes) --
  new_labels <- data.table(
    program_code     = existing_labels$program_code,
    biological_label = unname(CURATED_LABELS[existing_labels$program_code]),
    label_category   = unname(CURATED_LABELS[existing_labels$program_code]),
    top_pathway      = top_path_per_prog$top_pathway_str[
                         match(existing_labels$program_code, top_path_per_prog$program_code)],
    top_genes        = existing_labels$top_genes
  )
  fwrite(new_labels, labels_path)
  cat("Wrote program_labels.csv (", nrow(new_labels), " programs)\n", sep = "")

  # -- Remap nmf_assignments.csv dominant_program via CURATED_LABELS --
  if (!file.exists(assign_path)) {
    cat("nmf_assignments.csv not found — skipping remap (run Script 44 to generate).\n")
  } else {
    assign_dt <- fread(assign_path)
    # Backup only if absent
    if (!file.exists(assign_bak_path)) {
      file.copy(assign_path, assign_bak_path)
      cat("Created backup: nmf_assignments_prepathway_backup.csv\n")
    } else {
      cat("Backup already exists (skipping): nmf_assignments_prepathway_backup.csv\n")
    }
    # Remap dominant_program from dominant_program_code
    if ("dominant_program_code" %in% names(assign_dt)) {
      assign_dt[, dominant_program := unname(CURATED_LABELS[dominant_program_code])]
    } else {
      warning("nmf_assignments.csv has no dominant_program_code column — dominant_program not remapped.")
    }
    fwrite(assign_dt, assign_path)
    cat("Wrote nmf_assignments.csv (", nrow(assign_dt), " samples remapped)\n", sep = "")
  }

  cat("\n=== CANONICAL WRITE COMPLETE ===\n")
  cat("  program_labels.csv        : ", nrow(new_labels), " programs with curated names\n", sep = "")
  cat("  program_rename_mapping.csv: old->new label audit trail\n")
  cat("  Backups: skipped if already present (idempotent).\n")
}

cat("\nDone:", format(Sys.time()), "\n")
