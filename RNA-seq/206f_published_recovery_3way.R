#!/usr/bin/env Rscript
# 206f_published_recovery_3way.R
# Three-way DEG recovery comparison now that the canonical method is
# limma-voom quality-weighted C2 (canonical_deg_results.csv, 2026-06-08).
#
# For every human cohort with a usable published DEG list, measure how much
# of that published list is recovered by:
#   (A) our CANONICAL pooled (cohort-adjusted) limma-voom-qw C2 disease signature
#   (B) our MATCHED per-study DE on that cohort's own data, using the contrast
#       that MATCHES the published table:
#         disease-vs-control studies -> per_study/{GSE}_de_results.csv
#         stage/ordinal studies      -> results/audit_sensitivity/published_recovery/
#                                       matched_perstudy_{GSE}.csv (rebuilt 2026-06-11)
# reported side by side: recovery %, LFC Spearman rho (rank-based, robust to the
# extreme low-count LFC outliers in older published lists), % sign agreement, and
# a threshold-free median percentile of the published genes in our ranking.
#
# Honest annotations per study:
#   in_canonical_pool : TRUE for the 5 control-bearing cohorts pooled into (A);
#                       for the rest, (A) is genuine external replication.
#   pub_source        : "published_table" (their reported genes) vs
#                       "reproduced" (we re-ran their method; no table available)
#   contrast_type     : disease_vs_control / stage / ordinal / signature
#
# Output: figures/supplementary/figS_methods_validation/sensitivity/
#   published_recovery_3way.csv   (the side-by-side table)
#   figS_published_recovery_3way.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
PER  <- file.path(INT, "per_study")
PUB  <- file.path(BASE, "data/published_degs")
OUTDIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# optional publication theme (guarded)
theme_pub <- theme_bw
suppressWarnings(tryCatch(
  source(file.path(BASE, "scripts/figures/publication_theme.R")),
  error = function(e) NULL))

OUR_PADJ <- 0.1           # DEG threshold for OUR results when checking recovery of
                          # published lists (loosened 0.05->0.1 2026-06-11; matches the
                          # pipeline's exploratory-flagging threshold. Published-set
                          # thresholds are each study's own cutoff, unchanged.)
cat("=== 206f: Published DEG recovery — canonical (LVQW C2) vs per-study ===\n\n")

# ── Gene symbol -> Ensembl base map ─────────────────────────────────────────
annot <- fread(file.path(INT, "gene_annotation/human_ensg_to_symbol.tsv"))
annot_pc <- annot[gene_type == "protein_coding"]
sym2ens <- setNames(annot_pc$gene_base, annot_pc$symbol)
extra <- annot[!gene_base %in% annot_pc$gene_base & !symbol %in% names(sym2ens)]
extra <- extra[!duplicated(symbol)]
sym2ens <- c(sym2ens, setNames(extra$gene_base, extra$symbol))
sym2ens_up <- setNames(sym2ens, toupper(names(sym2ens)))
map_sym <- function(s) {
  s <- as.character(s); m <- sym2ens[s]
  miss <- is.na(m); if (any(miss)) m[miss] <- sym2ens_up[toupper(s[miss])]
  unname(m)
}
strip_ver <- function(x) sub("\\.[0-9]+$", "", as.character(x))

# ── Canonical pooled limma-voom-qw C2 ───────────────────────────────────────
can <- fread(file.path(INT, "integration/canonical_deg_results.csv"))
can[, gene_base := strip_ver(gene)]
can <- can[!is.na(padj)]
can[, neglogp := -log10(pmax(P.Value, .Machine$double.xmin))]
can[, pctile := frank(neglogp) / .N]          # threshold-free rank (1 = most sig)
can_sig <- can[padj < OUR_PADJ, gene_base]
cat(sprintf("Canonical: %d genes tested, %d DEG at padj<%.2f\n",
            nrow(can), length(can_sig), OUR_PADJ))

# ── Matched per-study loader ────────────────────────────────────────────────
# "Matched" = our DE on that cohort's own data, using the contrast that MATCHES
# the published table (not our generic disease-vs-control) — so the per-study
# arm is a fair head-to-head with the published list.
#   disease-vs-control studies -> per_study/{GSE}_de_results.csv (already matched)
#   stage / ordinal studies    -> matched contrasts rebuilt by subagents 2026-06-11
#                                 (RNA-seq/results/audit_sensitivity/published_recovery/)
MATCHED_DIR <- file.path(BASE, "RNA-seq/results/audit_sensitivity/published_recovery")
matched_map <- list(
  GSE126848   = list(type="perstudy", contrast="disease vs control"),
  GSE167523   = list(type="perstudy", contrast="NASH vs NAFL"),
  PRJNA512027 = list(type="perstudy", contrast="inflammation/fibrosis vs normal"),
  GSE130970   = list(type="matched", file="matched_perstudy_GSE130970.csv", contrast="NAS + fibrosis ordinal regression (matches Hoang design)"),
  GSE135251   = list(type="matched", file="matched_perstudy_GSE135251.csv", contrast="NASH F2/F3/F4 vs NAFL + NAS>=4"),
  GSE162694   = list(type="matched", file="matched_perstudy_GSE162694.csv", contrast="fibrosis F*vsF0 union"),
  GSE240729   = list(type="matched", file="matched_perstudy_GSE240729.csv", contrast="fibrosis-stage ordinal regression (well-powered analog of pairwise Fx-vs-F0)")
)
load_matched <- function(gse) {
  m <- matched_map[[gse]]; if (is.null(m)) return(NULL)
  if (m$type == "perstudy") {
    f <- file.path(PER, paste0(gse, "_de_results.csv")); if (!file.exists(f)) return(NULL)
    d <- fread(f)[!is.na(adj.P.Val)]
    d <- d[, .(gene_base = strip_ver(gene), logFC, padj = adj.P.Val)]
  } else {
    f <- file.path(MATCHED_DIR, m$file); if (!file.exists(f)) return(NULL)
    d <- fread(f)[!is.na(our_padj)]
    d <- d[, .(gene_base, logFC = our_logFC, padj = our_padj)]
  }
  d[, neglogp := -log10(pmax(padj, .Machine$double.xmin))]
  d[, pctile  := frank(neglogp) / .N]
  d[]
}

# ── Published-DEG loaders: each returns data.table(gene_base, pub_dir, pub_lfc)
# pub_dir in {"up","down",NA}; pub_lfc numeric (log2FC or coefficient) or NA.
# one row per published DEG (deduped, keeping most significant).
dedup_pub <- function(dt) {
  if (!"pub_lfc" %in% names(dt)) dt[, pub_lfc := NA_real_]
  dt <- dt[!is.na(gene_base) & gene_base != ""]
  dt <- dt[order(rank_key)]
  dt[!duplicated(gene_base), .(gene_base, pub_dir, pub_lfc)]
}

load_published <- function(study) switch(study,

  # GSE126848 Suppli — reproduced DESeq2, NAFL/NASH vs normal-weight controls
  "GSE126848" = {
    d <- fread(file.path(PUB, "GSE126848/deseq2_NAFL_NASH_vs_Normalweight.csv"))
    setnames(d, 1, "gene")
    d <- d[!is.na(padj) & padj < 0.05]
    d[, `:=`(gene_base = strip_ver(gene),
             pub_dir = fifelse(log2FoldChange > 0, "up", "down"),
             pub_lfc = log2FoldChange,
             rank_key = padj)]
    dedup_pub(d)
  },

  # GSE130970 Hoang — MOESM2 NAS + fibrosis ordinal regression, adj_P<0.01
  "GSE130970" = {
    xl <- file.path(PUB, "GSE130970/MOESM2.xlsx")
    rd <- function(sheet) {
      r <- as.data.table(read_excel(xl, sheet = sheet))
      setnames(r, names(r), tolower(gsub("\\s+", "_", names(r))))
      r <- r[adj_p < 0.01]
      r[, `:=`(gene_base = map_sym(gene_symbol),
               pub_dir = fifelse(coefficient > 0, "up", "down"),
               pub_lfc = as.numeric(coefficient),
               rank_key = adj_p)]
      r
    }
    dedup_pub(rbindlist(list(rd("NAS ordinal regression"),
                             rd("fibrosis ordinal regression")), use.names = TRUE, fill = TRUE))
  },

  # GSE135251 Govaere — Sci Transl Med S-tables (already q<0.05 & |FC|>1.5)
  "GSE135251" = {
    xl <- file.path(PUB, "GSE135251/aba4448_supplementary_tables.xlsx")
    sheets <- c("Table S3", "Table S4", "Table S5", "Table S8", "Table S9")
    rd <- function(sn) {
      r <- as.data.table(read_excel(xl, sheet = sn, skip = 1))
      setnames(r, names(r), tolower(gsub("\\s+", "_", names(r))))
      ens <- grep("ensembl|gene_id", names(r), value = TRUE)[1]
      lfc <- grep("log2?fc", names(r), value = TRUE)[1]
      q   <- grep("q_value|qvalue|q.value", names(r), value = TRUE)[1]
      r[, `:=`(gene_base = strip_ver(r[[ens]]),
               pub_dir = fifelse(as.numeric(r[[lfc]]) > 0, "up", "down"),
               pub_lfc = as.numeric(r[[lfc]]),
               rank_key = as.numeric(r[[q]]))]
      r[, .(gene_base, pub_dir, pub_lfc, rank_key)]
    }
    dedup_pub(rbindlist(lapply(sheets, rd)))
  },

  # GSE162694 Pantano — 98-gene fibrosis signature (Table S3)
  "GSE162694" = {
    r <- as.data.table(read_excel(file.path(PUB, "GSE162694/supp_info_3.xlsx"),
                                  sheet = "Table S3", skip = 1))
    # columns: Column1, gene(ENSG), coefficient, entrezgene_id, external_gene_name, ...
    setnames(r, 2, "gene"); setnames(r, 3, "coef")
    r <- r[grepl("^ENSG", gene)]
    r[, `:=`(gene_base = strip_ver(gene),
             pub_dir = fifelse(as.numeric(coef) > 0, "up", "down"),
             pub_lfc = as.numeric(coef),
             rank_key = seq_len(.N))]
    dedup_pub(r)
  },

  # GSE167523 Kozumi — 137 upregulated genes, NASH vs NAFL (up-only; no LFC spread)
  "GSE167523" = {
    r <- fread(file.path(PUB, "GSE167523/table_s2_genes.csv"))
    r[, `:=`(gene_base = map_sym(symbol), pub_dir = "up",
             pub_lfc = log2(as.numeric(fold_change)), rank_key = adj_pval)]
    dedup_pub(r)
  },

  # GSE240729 Verschuren — REAL published human log2FC from Source Data Fig 4
  # (source_data.xlsx "Figure 4": F3 block cols 8-9, F4 block cols 14-15 =
  # human-mouse concordant gene subset, n~528). The deseq2_F*_vs_F0_degs.csv
  # files are an internal reproduction the GEO README flags as unreproducible
  # (muscle-contamination + F0/F4 orientation artifacts) — DO NOT use them.
  "GSE240729" = {
    r <- as.data.table(read_excel(file.path(PUB, "GSE240729/source_data.xlsx"),
                                  sheet = "Figure 4", col_names = FALSE, skip = 2))
    blk <- function(symc, lfcc)
      data.table(sym = as.character(r[[symc]]),
                 lfc = suppressWarnings(as.numeric(r[[lfcc]])))[!is.na(sym) & !is.na(lfc)]
    d <- rbindlist(list(blk(8, 9), blk(14, 15)))      # F3 (8-9) + F4 (14-15)
    d[, `:=`(gene_base = map_sym(sym),
             pub_dir = fifelse(lfc > 0, "up", "down"),
             pub_lfc = lfc,
             rank_key = -abs(lfc))]                    # keep largest |logFC| per gene
    dedup_pub(d)
  },

  # PRJNA512027 Gerhard — inflammation + fibrosis vs normal (q<0.05, |log2FC|>=1)
  "PRJNA512027" = {
    rd <- function(tab) {
      d <- fread(file.path(PUB, sprintf("PRJNA512027/%s", tab)))
      setnames(d, c("Gene_ID", "log2FC2", "q-value3"), c("sym","lfc","q"), skip_absent = TRUE)
      d[, `:=`(gene_base = map_sym(sym),
               pub_dir = fifelse(lfc > 0, "up", "down"),
               pub_lfc = as.numeric(lfc),
               rank_key = q)]
      d[, .(gene_base, pub_dir, pub_lfc, rank_key)]
    }
    dedup_pub(rbindlist(list(rd("table_s2_inflammation_vs_normal.csv"),
                             rd("table_s3_fibrosis_vs_normal.csv"))))
  },

  NULL
)

# ── Study configuration ─────────────────────────────────────────────────────
cfg <- data.table(
  study      = c("GSE126848","GSE130970","GSE135251","GSE162694","GSE167523",
                 "GSE213621","GSE240729","PRJNA512027","GSE174478","GSE193066"),
  label      = c("Suppli 2019","Hoang 2019","Govaere 2020","Pantano 2021","Kozumi 2021",
                 "Chen 2023","Verschuren 2024","Gerhard 2018","Kawamura 2022","Fujiwara 2022"),
  in_pool    = c(TRUE,TRUE,TRUE,TRUE,FALSE, TRUE,FALSE,FALSE,FALSE,FALSE),
  pub_source = c("reproduced","published_table","published_table","published_table",
                 "published_table","none","published_table","published_table","none","signature"),
  contrast   = c("disease_vs_control","ordinal","stage","signature","stage_up_only",
                 "none","stage","disease_vs_control","none","prognostic"),
  note       = c("", "", "per-study=disease-vs-ctrl (published=stage; expected under-recovery)",
                 "", "no healthy controls; per-study=NASH-vs-NAFL",
                 "no published DEG list (validation-only cohort)",
                 "real published log2FC from Source Data Fig 4 (human-mouse concordant subset, n~528)",
                 "DROPPED for L0/S0 batch confound; published list artifact-dominated",
                 "no human DE performed (validation-only cohort)",
                 "prognostic PLS signature, not a DEG contrast")
)

# ── Recovery computation ────────────────────────────────────────────────────
# recovery % : fraction of published DEGs (present in our tested universe) that
#              are DEG in ours at OUR_PADJ.
# lfc_rho    : Spearman corr of published LFC vs our LFC over ALL shared genes
#              (rank-based; immune to extreme low-count LFC outliers).
# sign_agree : % sign agreement of published vs our LFC over all shared genes.
# med_pctile : median significance-percentile of published genes in our ranking.
recovery_vs <- function(pub, res, sig_genes) {
  if ("logFC" %in% names(res)) lfc_col <- "logFC" else lfc_col <- "shrunk_logFC"
  in_univ <- pub[gene_base %in% res$gene_base]
  n_univ  <- nrow(in_univ)
  if (n_univ == 0)
    return(list(n_univ = 0, rec_pct = NA, lfc_rho = NA, sign_agree = NA, med_pctile = NA))
  rec     <- in_univ[gene_base %in% sig_genes]
  our_lfc <- res[match(in_univ$gene_base, gene_base), get(lfc_col)]
  med_pctile <- median(res[match(in_univ$gene_base, gene_base), pctile], na.rm = TRUE)
  # LFC concordance (only if published LFC available and has spread)
  have_lfc <- !is.na(in_univ$pub_lfc) & is.finite(our_lfc)
  if (sum(have_lfc) >= 10 && length(unique(in_univ$pub_lfc[have_lfc])) > 5) {
    lfc_rho    <- cor(in_univ$pub_lfc[have_lfc], our_lfc[have_lfc], method = "spearman")
    sign_agree <- mean(sign(in_univ$pub_lfc[have_lfc]) == sign(our_lfc[have_lfc]))
  } else { lfc_rho <- NA; sign_agree <- NA }
  list(n_univ = n_univ, rec_pct = 100 * nrow(rec) / n_univ,
       lfc_rho = lfc_rho, sign_agree = sign_agree, med_pctile = med_pctile)
}

rows <- list(); pub_store <- list()
for (i in seq_len(nrow(cfg))) {
  st <- cfg$study[i]
  cat(sprintf("[%s] %s ", st, cfg$label[i]))
  pub <- tryCatch(load_published(st), error = function(e) {
    cat(sprintf("(published load FAILED: %s) ", conditionMessage(e))); NULL })
  if (!is.null(pub) && nrow(pub) > 0) pub_store[[st]] <- pub
  ps  <- load_matched(st)
  n_ps_deg <- if (!is.null(ps)) sum(ps$padj < OUR_PADJ) else NA
  matched_contrast <- if (!is.null(matched_map[[st]])) matched_map[[st]]$contrast else NA

  base <- data.table(study = st, label = cfg$label[i], in_canonical_pool = cfg$in_pool[i],
                     pub_source = cfg$pub_source[i], contrast_type = cfg$contrast[i],
                     matched_contrast = matched_contrast,
                     n_published = if (is.null(pub)) NA_integer_ else nrow(pub),
                     perstudy_n_deg = n_ps_deg, note = cfg$note[i])
  if (is.null(pub) || nrow(pub) == 0) {
    cat("-> no published DEG list\n")
    rows[[st]] <- cbind(base, data.table(
      can_n_in_universe=NA, can_recovery_pct=NA, can_lfc_rho=NA, can_sign_agree=NA, can_med_pctile=NA,
      ps_n_in_universe=NA,  ps_recovery_pct=NA,  ps_lfc_rho=NA,  ps_sign_agree=NA,  ps_med_pctile=NA))
    next
  }
  canR <- recovery_vs(pub, can, can_sig)
  if (!is.null(ps)) {
    psR <- recovery_vs(pub, ps, ps[padj < OUR_PADJ, gene_base])
  } else psR <- list(n_univ=NA, rec_pct=NA, lfc_rho=NA, sign_agree=NA, med_pctile=NA)

  cat(sprintf("-> %d pub | canon: rec %.0f%% rho %+.2f | perstudy: rec %s%% rho %s\n",
              nrow(pub), canR$rec_pct, ifelse(is.na(canR$lfc_rho),0,canR$lfc_rho),
              ifelse(is.na(psR$rec_pct),"NA",sprintf("%.0f",psR$rec_pct)),
              ifelse(is.na(psR$lfc_rho),"NA",sprintf("%+.2f",psR$lfc_rho))))
  rows[[st]] <- cbind(base, data.table(
    can_n_in_universe=canR$n_univ, can_recovery_pct=round(canR$rec_pct,1),
    can_lfc_rho=round(canR$lfc_rho,3), can_sign_agree=round(100*canR$sign_agree,1),
    can_med_pctile=round(100*canR$med_pctile,1),
    ps_n_in_universe=psR$n_univ, ps_recovery_pct=round(ifelse(is.na(psR$rec_pct),NA,psR$rec_pct),1),
    ps_lfc_rho=round(psR$lfc_rho,3), ps_sign_agree=round(100*psR$sign_agree,1),
    ps_med_pctile=round(100*ifelse(is.na(psR$med_pctile),NA,psR$med_pctile),1)))
}
res_tbl <- rbindlist(rows, use.names = TRUE, fill = TRUE)
fwrite(res_tbl, file.path(OUTDIR, "published_recovery_3way.csv"))
cat(sprintf("\nWrote %s\n", file.path(OUTDIR, "published_recovery_3way.csv")))
print(res_tbl[!is.na(n_published),
      .(label, in_pool=in_canonical_pool, src=pub_source, n_pub=n_published,
        can_rec=can_recovery_pct, can_rho=can_lfc_rho,
        ps_rec=ps_recovery_pct, ps_rho=ps_lfc_rho)])

# ════════════════════════════════════════════════════════════════════════════
# FIGURE: systematic reproduction of published DEGs by our canonical signature
# ════════════════════════════════════════════════════════════════════════════
DV <- "#D41159"; PV <- "#1A85FF"; GREY <- "#9E9E9E"

# ── Panel A: TWO-ARM LFC concordance scatter (published vs canonical AND vs
#    matched per-study), faceted by study. Only studies with a SIGNED published
#    effect size (Pantano=unsigned importance and Kozumi=up-only are excluded;
#    they are assessed by recovery % in panel B).
SIGNED <- c("GSE126848","GSE130970","GSE135251","PRJNA512027","GSE240729")
arm_lvls <- c("Integrated","Matched per-study")
scatter_dat <- rbindlist(lapply(intersect(SIGNED, names(pub_store)), function(st) {
  p <- pub_store[[st]][!is.na(pub_lfc), .(gene_base, pub_lfc)]
  if (nrow(p) < 10) return(NULL)
  out <- list()
  mc <- merge(p, can[, .(gene_base, our_lfc = logFC)], by = "gene_base")
  if (nrow(mc) >= 10) { mc[, arm := "Integrated"]; out[["c"]] <- mc }
  ps <- load_matched(st)
  if (!is.null(ps)) {
    mm <- merge(p, ps[, .(gene_base, our_lfc = logFC)], by = "gene_base")
    if (nrow(mm) >= 10) { mm[, arm := "Matched per-study"]; out[["p"]] <- mm }
  }
  if (!length(out)) return(NULL)
  r <- rbindlist(out); r[, study := st]; r
}), use.names = TRUE)
scatter_dat[, pub_lfc_w := pmax(pmin(pub_lfc, 6), -6)]   # winsorize display only
scatter_dat[, arm := factor(arm, levels = arm_lvls)]
# facet label carries BOTH rho values
scatter_dat[, facet := sprintf("%s\nintegrated rho=%+.2f\nmatched rho=%+.2f",
              cfg$label[match(study, cfg$study)],
              res_tbl$can_lfc_rho[match(study, res_tbl$study)],
              res_tbl$ps_lfc_rho[match(study, res_tbl$study)])]
ford <- res_tbl[study %in% scatter_dat$study][order(-ps_lfc_rho)]
scatter_dat[, facet := factor(facet, levels = sapply(ford$study, function(s) scatter_dat[study==s, facet[1]]))]

pA <- ggplot(scatter_dat, aes(pub_lfc_w, our_lfc, color = arm)) +
  geom_hline(yintercept = 0, color = GREY, linewidth = 0.3) +
  geom_vline(xintercept = 0, color = GREY, linewidth = 0.3) +
  geom_point(alpha = 0.2, size = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = c("Integrated"=DV, "Matched per-study"=PV)) +
  facet_wrap(~ facet, ncol = 5, scales = "free_x") +
  labs(title = "Published vs our DEG effect sizes",
       x = "Published log2 fold-change", y = "Our logFC", color = NULL) +
  theme_bw(base_size = 16) + theme_pub() +
  theme(strip.text  = element_text(size = 10.5, lineheight = 0.95),
        plot.title  = element_text(size = 18, face = "bold"),
        axis.title  = element_text(size = 15),
        axis.text   = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.position = "top")

# ── Panel B: recovery % side by side (canonical vs matched per-study) ────────
plt <- res_tbl[!is.na(n_published) & !is.na(can_recovery_pct)][order(can_recovery_pct)]
plt[, label2 := sprintf("%s (n=%d%s)", label, n_published,
                        ifelse(pub_source=="reproduced"," *","" ))]
plt[, label2 := factor(label2, levels = label2)]
mlt <- melt(plt, id.vars = "label2",
            measure.vars = c("ps_recovery_pct","can_recovery_pct"),
            variable.name = "arm", value.name = "recovery")
mlt[, arm := factor(arm, levels=c("ps_recovery_pct","can_recovery_pct"),
                    labels=c("Matched per-study","Integrated"))]
pB <- ggplot(mlt, aes(label2, recovery, fill = arm)) +
  geom_col(position = position_dodge(width=0.78), width=0.72) +
  geom_text(aes(label = ifelse(is.na(recovery),"n/a",sprintf("%.0f",recovery))),
            position = position_dodge(width=0.78), hjust=-0.2, size=3.6) +
  coord_flip(clip="off") + scale_y_continuous(limits=c(0,108), breaks=seq(0,100,25)) +
  scale_fill_manual(values = c("Matched per-study"=PV, "Integrated"=DV)) +
  labs(title="DEG recovery: integrated vs matched per-study",
       x=NULL, y="Recovery (%)", fill=NULL) +
  theme_bw(base_size = 13) + theme_pub() + theme(legend.position="top")

# ── Panel C: threshold-free — median significance percentile ────────────────
pctd <- res_tbl[!is.na(can_med_pctile)][order(can_med_pctile)]
pctd[, label2 := factor(label, levels = label)]
pC <- ggplot(pctd, aes(label2, can_med_pctile)) +
  geom_segment(aes(xend=label2, y=50, yend=can_med_pctile), color=GREY, linewidth=0.4) +
  geom_point(color=DV, size=3.6) +
  geom_hline(yintercept=50, linetype="dashed", color=GREY) +
  geom_text(aes(label=sprintf("%.0f",can_med_pctile)), hjust=-0.45, size=3.6) +
  coord_flip(clip="off") + scale_y_continuous(limits=c(40,103)) +
  labs(title="Where published DEGs rank in our signature",
       x=NULL, y="Median percentile") +
  theme_bw(base_size = 13) + theme_pub()

# ── Individual panels (proportionally sized, larger fonts) ──────────────────
ggsave(file.path(OUTDIR, "figS_published_recovery_panelA_lfc_scatter.pdf"),
       pA, width = 13.5, height = 6.2, useDingbats = FALSE)
ggsave(file.path(OUTDIR, "figS_published_recovery_panelB_recovery.pdf"),
       pB, width = 6.2, height = 3.9, useDingbats = FALSE)
ggsave(file.path(OUTDIR, "figS_published_recovery_panelC_percentile.pdf"),
       pC, width = 6.2, height = 3.9, useDingbats = FALSE)
# combined (reference) — kept, but the individual panels above are the primary deliverable
fig <- pA / (pB | pC) + plot_layout(heights = c(1.1, 1.2))
ggsave(file.path(OUTDIR, "figS_published_recovery_3way.pdf"),
       fig, width = 15, height = 11, useDingbats = FALSE)
cat(sprintf("Wrote individual panels A/B/C + combined to %s\n", OUTDIR))
cat("\n=== done ===\n")
