suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
  library(msigdbr)
  library(ggplot2)
})

pdf.options(useDingbats = FALSE)

cat("=== Phase 3b: TF Activity + PROGENy Concordance ===\n\n")

BASE    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
H_INT   <- file.path(BASE, "Human/Patient_Cohorts/analysis/integration")
ANNOT   <- file.path(H_INT, "results/gene_annotation")
DS_DIR  <- file.path(H_INT, "results/disease_signatures")
INT_DIR <- file.path(H_INT, "results/integration")
MOUSE_PD <- file.path(BASE, "Mouse/Unified_Integration/results/per_diet")
WD      <- file.path(BASE, "Analysis/Cross_Species_Concordance")
RES     <- file.path(WD, "results")

DIETS <- c("MCD", "HFD", "CDAHFD", "FPC")  # 2026-05-29: LIDPAD archived/dropped

# ============================================================
#  Part A: TF Activity via MSigDB TF target gene sets + decoupleR
# ============================================================
cat("=== Part A: TF Activity ===\n")

# Build TF regulon network from MSigDB C3:TFT (Regulatory Target Gene Sets)
# This avoids needing OmniPath/CollecTRI/DoRothEA (unreachable from HPC)
cat("Building TF regulon network from MSigDB C3:TFT...\n")
h_tft <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C3", subcollection = "TFT:GTRD"))
m_tft <- as.data.table(msigdbr(species = "Mus musculus", collection = "C3", subcollection = "TFT:GTRD"))

# Extract TF name from gene set name (format: GTRD_TF_TARGETS)
# Convert to decoupleR format: source (TF), target (gene), mor (1 = activation)
#
# LIMITATION: MSigDB TFT gene sets do not distinguish activation vs repression.
# All mor=1 (activation assumed). This means:
# - TF activity signs may be inverted for repressive TFs
# - Concordance of activity *direction* is interpretable only for activators
# - The magnitude and rank of concordance (rho) remains informative
# - We report rho_abs (magnitude-only concordance) to address this
# For signed regulons, switch to CollecTRI/DoRothEA when server is accessible.
build_tf_net <- function(tft_dt) {
  net <- tft_dt[, .(source = gsub("_.*", "", gs_name), target = gene_symbol, mor = 1)]
  # Remove duplicate edges (same TF-target pair from overlapping gene sets)
  net <- unique(net, by = c("source", "target"))
  # Keep only TFs with at least 10 targets
  tf_counts <- net[, .N, by = source]
  good_tfs <- tf_counts[N >= 10 & N <= 2000, source]
  net <- net[source %in% good_tfs]
  net
}

h_net_tf <- build_tf_net(h_tft)
m_net_tf <- build_tf_net(m_tft)
cat(sprintf("  Human TF network: %d interactions, %d TFs\n",
  nrow(h_net_tf), uniqueN(h_net_tf$source)))
cat(sprintf("  Mouse TF network: %d interactions, %d TFs\n",
  nrow(m_net_tf), uniqueN(m_net_tf$source)))

# ============================================================
#  Run decoupleR on DE results (stat-based approach)
# ============================================================
# Helper: extract t-statistics as a named matrix for decoupleR
make_stat_mat <- function(dt, t_col = "t", symbol_col = "symbol") {
  dt <- dt[!is.na(get(symbol_col)) & get(symbol_col) != ""]
  # Average duplicates
  dt_avg <- dt[, .(t_val = mean(get(t_col), na.rm = TRUE)), by = symbol_col]
  mat <- matrix(dt_avg$t_val, ncol = 1, dimnames = list(dt_avg[[symbol_col]], "condition"))
  return(mat)
}

cat("\n--- Human TF activity ---\n")
human_tf_results <- list()

human_files <- list(
  disease_vs_ctrl = list(path = file.path(INT_DIR, "dream_results.csv"), padj_col = "padj"),
  nafl_vs_nash    = list(path = file.path(DS_DIR, "nafl_vs_nash_dream.csv"), padj_col = "adj.P.Val"),
  fibrosis        = list(path = file.path(DS_DIR, "fibrosis_dream.csv"), padj_col = "adj.P.Val")
)

h_annot <- fread(file.path(ANNOT, "human_ensg_to_symbol.tsv"))

for (sig_name in names(human_files)) {
  info <- human_files[[sig_name]]
  dt <- fread(info$path)
  if (!"symbol" %in% names(dt)) {
    dt[, gene_base := gsub("\\..*", "", gene)]
    dt <- merge(dt, h_annot[, .(gene_base, symbol)], by = "gene_base", all.x = TRUE)
  }

  mat <- make_stat_mat(dt, "t", "symbol")
  cat(sprintf("  %s: %d genes\n", sig_name, nrow(mat)))

  tf_act <- run_wmean(mat, h_net_tf, .source = "source", .target = "target",
                       .mor = "mor", times = 1000, minsize = 5)
  tf_act <- as.data.table(tf_act)
  tf_act[, source_name := sig_name]
  human_tf_results[[sig_name]] <- tf_act[statistic == "norm_wmean",
    .(source, score, p_value, source_name)]
}

human_tf_all <- rbindlist(human_tf_results)
cat(sprintf("  Total TFs scored: %d\n", uniqueN(human_tf_all$source)))

cat("\n--- Mouse TF activity ---\n")
mouse_tf_results <- list()
m_annot <- fread(file.path(ANNOT, "mouse_ensmusg_to_symbol.tsv"))

for (diet in DIETS) {
  dt <- fread(file.path(MOUSE_PD, paste0(diet, "_de_results.csv")))
  if (!"symbol" %in% names(dt)) {
    dt[, gene_base := gsub("\\..*", "", gene)]
    dt <- merge(dt, m_annot[, .(gene_base, symbol)], by = "gene_base", all.x = TRUE)
  }

  mat <- make_stat_mat(dt, "t", "symbol")
  cat(sprintf("  %s: %d genes\n", diet, nrow(mat)))

  tf_act <- run_wmean(mat, m_net_tf, .source = "source", .target = "target",
                       .mor = "mor", times = 1000, minsize = 5)
  tf_act <- as.data.table(tf_act)
  tf_act[, source_name := diet]
  mouse_tf_results[[diet]] <- tf_act[statistic == "norm_wmean",
    .(source, score, p_value, source_name)]
}

mouse_tf_all <- rbindlist(mouse_tf_results)

# Save
fwrite(human_tf_all, file.path(RES, "tf_activity_human.csv"))
fwrite(mouse_tf_all, file.path(RES, "tf_activity_mouse.csv"))

# ============================================================
#  TF concordance
# ============================================================
cat("\n=== TF Concordance ===\n")

tf_conc <- data.table()
for (sig_name in names(human_files)) {
  h_tf <- human_tf_all[source_name == sig_name, .(source, h_score = score, h_pval = p_value)]

  for (diet in DIETS) {
    m_tf <- mouse_tf_all[source_name == diet, .(source, m_score = score, m_pval = p_value)]

    # TF names may differ (case); normalize
    h_tf[, tf_name := toupper(source)]
    m_tf[, tf_name := toupper(source)]
    paired <- merge(h_tf, m_tf, by = "tf_name")

    if (nrow(paired) < 5) next

    rho <- cor(paired$h_score, paired$m_score, method = "spearman", use = "complete.obs")
    # Unsigned concordance (immune to mor=1 sign assumption)
    rho_abs <- cor(abs(paired$h_score), abs(paired$m_score),
                   method = "spearman", use = "complete.obs")
    both_active <- paired[h_pval < 0.05 & m_pval < 0.05]
    conc <- if (nrow(both_active) > 0) sum(sign(both_active$h_score) == sign(both_active$m_score)) else 0
    disc <- if (nrow(both_active) > 0) sum(sign(both_active$h_score) != sign(both_active$m_score)) else 0

    row <- data.table(
      human_signature = sig_name,
      diet = diet,
      n_common_tfs = nrow(paired),
      rho_tf = round(rho, 4),
      rho_tf_abs = round(rho_abs, 4),
      n_both_active = nrow(both_active),
      n_concordant = conc,
      n_discordant = disc
    )
    tf_conc <- rbindlist(list(tf_conc, row))

    cat(sprintf("  %s × %s: ρ_TF=%.3f (abs=%.3f), active_conc=%d/%d\n",
      sig_name, diet, rho, rho_abs, conc, conc + disc))
  }
}

fwrite(tf_conc, file.path(RES, "tf_concordance.csv"))

# ============================================================
#  Part B: PROGENy Pathway Activity
# ============================================================
cat("\n=== Part B: PROGENy Activity ===\n")

# Get PROGENy model
h_progeny <- tryCatch(
  get_progeny(organism = "human", top = 500),
  error = function(e) {
    cat("PROGENy download failed:", e$message, "\n")
    NULL
  }
)
m_progeny <- tryCatch(
  get_progeny(organism = "mouse", top = 500),
  error = function(e) {
    cat("PROGENy download failed:", e$message, "\n")
    NULL
  }
)

if (!is.null(h_progeny) && !is.null(m_progeny)) {
  cat(sprintf("  Human PROGENy: %d footprints, %d pathways\n",
    nrow(h_progeny), uniqueN(h_progeny$source)))
  cat(sprintf("  Mouse PROGENy: %d footprints, %d pathways\n",
    nrow(m_progeny), uniqueN(m_progeny$source)))

  # Run on human signatures
  human_prog <- list()
  for (sig_name in names(human_files)) {
    info <- human_files[[sig_name]]
    dt <- fread(info$path)
    if (!"symbol" %in% names(dt)) {
      dt[, gene_base := gsub("\\..*", "", gene)]
      dt <- merge(dt, h_annot[, .(gene_base, symbol)], by = "gene_base", all.x = TRUE)
    }
    mat <- make_stat_mat(dt, "t", "symbol")
    pa <- run_wmean(mat, h_progeny, .source = "source", .target = "target",
                     .mor = "weight", times = 1000, minsize = 5)
    pa <- as.data.table(pa)
    pa[, source_name := sig_name]
    human_prog[[sig_name]] <- pa[statistic == "norm_wmean"]
  }

  # Run on mouse diets
  mouse_prog <- list()
  for (diet in DIETS) {
    dt <- fread(file.path(MOUSE_PD, paste0(diet, "_de_results.csv")))
    if (!"symbol" %in% names(dt)) {
      dt[, gene_base := gsub("\\..*", "", gene)]
      dt <- merge(dt, m_annot[, .(gene_base, symbol)], by = "gene_base", all.x = TRUE)
    }
    mat <- make_stat_mat(dt, "t", "symbol")
    pa <- run_wmean(mat, m_progeny, .source = "source", .target = "target",
                     .mor = "weight", times = 1000, minsize = 5)
    pa <- as.data.table(pa)
    pa[, source_name := diet]
    mouse_prog[[diet]] <- pa[statistic == "norm_wmean"]
  }

  h_prog_all <- rbindlist(human_prog)
  m_prog_all <- rbindlist(mouse_prog)
  fwrite(h_prog_all[, .(source, score, p_value, source_name)],
    file.path(RES, "progeny_scores_human.csv"))
  fwrite(m_prog_all[, .(source, score, p_value, source_name)],
    file.path(RES, "progeny_scores_mouse.csv"))

  # PROGENy concordance
  prog_conc <- data.table()
  for (sig_name in names(human_files)) {
    h_p <- h_prog_all[source_name == sig_name, .(source, h_score = score, h_pval = p_value)]
    for (diet in DIETS) {
      m_p <- m_prog_all[source_name == diet, .(source, m_score = score, m_pval = p_value)]
      paired <- merge(h_p, m_p, by = "source")
      if (nrow(paired) < 3) next
      rho <- cor(paired$h_score, paired$m_score, method = "spearman", use = "complete.obs")
      row <- data.table(
        human_signature = sig_name, diet = diet,
        n_pathways = nrow(paired), rho_progeny = round(rho, 4)
      )
      prog_conc <- rbindlist(list(prog_conc, row))
      cat(sprintf("  PROGENy %s × %s: ρ=%.3f\n", sig_name, diet, rho))
    }
  }
  fwrite(prog_conc, file.path(RES, "progeny_concordance.csv"))
} else {
  cat("  PROGENy skipped (download failed)\n")
  # Write empty CSVs so downstream scripts can load without error
  fwrite(data.table(source = character(), score = numeric(),
    p_value = numeric(), source_name = character()),
    file.path(RES, "progeny_scores_human.csv"))
  fwrite(data.table(source = character(), score = numeric(),
    p_value = numeric(), source_name = character()),
    file.path(RES, "progeny_scores_mouse.csv"))
  fwrite(data.table(human_signature = character(), diet = character(),
    n_pathways = integer(), rho_progeny = numeric()),
    file.path(RES, "progeny_concordance.csv"))
}

cat("\nSaved: tf_activity_*.csv, tf_concordance.csv, progeny_*.csv\n")
cat("=== Phase 3b complete ===\n")
