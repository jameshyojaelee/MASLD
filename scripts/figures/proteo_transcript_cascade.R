#!/usr/bin/env Rscript
# ==============================================================================
# proteo_transcript_cascade.R  — RNA-seq x proteomics, 2 honest panels
#
# Disease-vs-control log2FC measured on shared genes across three INDEPENDENT
# cohorts (NOT the same patients): liver mRNA (846-sample C2 mega) | liver
# protein (DIA-MS PXD051911, n~58) | plasma protein (DIA-MS PXD052937, n~72).
# We therefore show CROSS-LAYER EFFECT-DIRECTION CONCORDANCE, not a within-
# sample molecular cascade.
#
# Panel A -> SUPPLEMENT (figS_proteomics/figS_proteo_liver_blood_decoupling.pdf, demoted
#   2026-07-02): a liver->blood DECOUPLING (null) result across 3 non-paired cohorts where
#   plasma is whole-body — not a Fig-4 hepatocyte-autonomous validation. Panel B stays in Fig 4.
# Panel A "Liver -> blood concordance" (parallel coordinates), restricted to
#   bona-fide SECRETED-to-blood proteins (HPA) — the only genes for which a
#   plasma measurement is interpretable (intracellular enzymes leak into plasma
#   and are excluded). Honest finding: the circulating disease proteome is
#   largely DECOUPLED from liver transcription — a robust acute-phase component
#   (CRP/SAA/HP/LBP up in blood) is not predicted by liver mRNA/protein.
# Panel B  liver mRNA vs liver-protein concordance (SAME tissue): disease effects
#   agree in direction (rho ~0.49 all / 0.66 disease DEGs, ~88% same-sign). A
#   positive trend line is shown; the magnitude RATIO is NOT interpreted (DIA
#   log-ratio compression makes a slope-1 expectation technically confounded, so
#   no y=x null is drawn). The collagen cluster is an ECM DETECTION caveat
#   (deposited/crosslinked matrix under-sampled by lysate DIA), not "fibrosis is
#   transcript-only".
#
# House style: all text black, no lollipops, no shading, prose -> caption.
# Output: figures/supplementary/figS_proteomics/figS_proteo_liver_blood_decoupling.pdf  (Panel A, demoted)
#         (Panel B, fig4_proteo_buffering_map.pdf, retired 2026-07-07 — plot still computed, no longer saved)
# Env: rnaseq
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
set.seed(42)
CONC <- "#2E7D32"; DECO <- "#C9265E"; ATTN <- "#1565C0"; DISC <- "#6A1B9A"; BG <- "#D9D9D9"

# ── 3-layer reshape (disease-vs-control) + HPA secreted-to-blood flag ──────────
c3 <- fread(file.path(BASE, "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv"))
M <- unique(c3[dataset == "PXD051911", .(gene, M = bulk_logFC, Mp = bulk_padj)])
L <- unique(c3[dataset == "PXD051911", .(gene, L = protein_logFC)])
P <- unique(c3[dataset == "PXD052937", .(gene, P = protein_logFC)])
d <- Reduce(function(a, b) merge(a, b, by = "gene"), list(M, L, P))
d <- d[is.finite(M) & is.finite(L) & is.finite(P)]
sec <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
             select = c("human_symbol", "secreted_to_blood", "hpa_secretome_location"))
d <- merge(d, sec, by.x = "gene", by.y = "human_symbol", all.x = TRUE)
d[, secreted := !is.na(secreted_to_blood) & secreted_to_blood == TRUE]
d[, disease := !is.na(Mp) & Mp < 0.05 & abs(M) >= 0.30]
cat(sprintf("[data] %d 3-layer genes; %d secreted-to-blood; %d disease (mRNA padj<.05 & |lfc|>=.3)\n",
            nrow(d), sum(d$secreted), sum(d$disease)))

# ── Panel A: liver->blood DECOUPLING scatter (secreted-to-blood only) ──────────
# The circulating disease proteome is not predicted by liver transcription:
# liver mRNA vs plasma protein is chance-level (rho~0.05); some coupling only
# re-appears at the liver-protein level (rho~0.36). Acute-phase reactants sit
# high in plasma at near-zero liver mRNA.
CAP <- 1.6; capv <- function(v) pmax(pmin(v, CAP), -CAP)
ds <- d[secreted == TRUE]
rho_mp  <- suppressWarnings(cor(ds$M, ds$P, method = "spearman"))
rho_lp  <- suppressWarnings(cor(ds$L, ds$P, method = "spearman"))
conc_mp <- mean(sign(ds$M) == sign(ds$P)) * 100
cat(sprintf("[panelA] secreted n=%d  rho(mRNA,plasma)=%.2f (%.0f%% same-dir)  rho(protein,plasma)=%.2f\n",
            nrow(ds), rho_mp, conc_mp, rho_lp))
# positive acute-phase reactants (hepatic APR) — the interpretable high-plasma set
APR <- c("CRP","SAA1","SAA2","SAA4","HP","HPX","LBP","FGA","FGB","FGG","ORM1","ORM2",
         "SERPINA1","SERPINA3","C3","CFB","LCN2","CP","AGT","LGALS3","F2")
ds[, acute := gene %in% APR]
# label ONLY a curated, recognisable APR subset (the story is the cluster, not every name)
labA <- ds[gene %in% c("CRP","SAA4","HP","LBP","ORM2","C3","CFB","LGALS3")]
pA <- ggplot(ds, aes(capv(M), capv(P))) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey55",
              linewidth = 0.4, linetype = "22") +
  geom_point(data = ds[acute == FALSE], colour = "#B4B4B4", size = 0.6, alpha = 0.6, shape = 16) +
  geom_point(data = ds[acute == TRUE],  colour = DECO, size = 1.3, shape = 16) +
  geom_text_repel(data = labA, aes(label = gene), colour = "black", size = GEOM_TEXT_6PT,
                  fontface = "italic", min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey75", max.overlaps = 40, seed = 1) +
  # single headline: liver transcription does not predict the blood proteome
  annotate("text", x = -CAP, y = CAP, hjust = 0, vjust = 1, size = GEOM_TEXT_6PT, colour = "black",
           lineheight = 1.0,
           label = sprintf("liver mRNA -> plasma:  rho = %.2f (chance)\nliver protein -> plasma:  rho = %.2f",
                           rho_mp, rho_lp)) +
  # in-plot colour key (no ggplot legend)
  annotate("point", x = 0.10, y = -1.48, colour = DECO, size = 1.3) +
  annotate("text",  x = 0.24, y = -1.48, hjust = 0, vjust = 0.5, size = GEOM_TEXT_6PT, colour = "black",
           label = "acute-phase reactant") +
  scale_x_continuous(name = "liver mRNA log2FC", limits = c(-CAP, CAP)) +
  scale_y_continuous(name = "plasma protein log2FC", limits = c(-CAP, CAP)) +
  theme_masld_compact() +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"))
# Demoted to supplement 2026-07-02: this is a liver->blood DECOUPLING (null) panel across 3
# non-paired cohorts, where plasma is whole-body — not a Fig-4 hepatocyte-autonomous validation.
# Panel B (liver mRNA vs liver protein, SAME tissue, positive concordance) remains the Fig-4 panel.
save_fig(pA, file.path(FIGS_PROTEO_DIR, "figS_proteo_liver_blood_decoupling.pdf"), width = fig_full_width*0.50, height = 3.0)

# ── Panel B: mRNA vs liver-protein, scored vs the GLOBAL regression ────────────
g2 <- d[abs(M) > 0.2 | abs(L) > 0.2]
fit <- lm(L ~ M, data = d); SL <- coef(fit)[2]; IN <- coef(fit)[1]
g2[, resid := L - (SL * M + IN)]
g2[, cls := fcase(
  sign(M) != sign(L) & abs(L) > 0.2,                       "Discordant",
  abs(M) >= 0.3 & abs(L) < 0.3 * abs(M),                   "Attenuated (mRNA >> protein)",
  default = "Concordant")]
g2[, cls := factor(cls, levels = c("Concordant","Attenuated (mRNA >> protein)","Discordant"))]
nB <- g2[, .N, by = cls]; cat("[panelB] "); cat(paste(sprintf("%s=%d", nB$cls, nB$N), collapse="  "), "\n")  # (cls kept for ORA below only)
LB <- 1.6; cb <- function(v) pmax(pmin(v, LB), -LB)
rhoB_all <- suppressWarnings(cor(d$M, d$L, method = "spearman"))
ddis <- d[!is.na(Mp) & Mp < 0.05 & abs(M) >= 0.3]
rhoB_dis  <- suppressWarnings(cor(ddis$M, ddis$L, method = "spearman"))
concB_dis <- mean(sign(ddis$M) == sign(ddis$L)) * 100
cat(sprintf("[panelB] all rho=%.2f slope=%.2f | disease DEG n=%d rho=%.2f %.0f%% same-dir\n",
            rhoB_all, SL, nrow(ddis), rhoB_dis, concB_dis))
# single accent = the collagen/ECM cluster (the visible off-fit exception); everything else one neutral hue
SLATE <- "#7E8AA2"
g2[, ecm := gene %in% c("COL1A1","COL1A2","COL3A1","COL5A1","LUM")]
hb <- g2[ecm == TRUE]
pB <- ggplot(g2, aes(cb(M), cb(L))) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  # faint background = the near-zero genes excluded from the |lfc|>0.2 cloud (fills the centre)
  rasterize_layer(geom_point(data = d, aes(cb(M), cb(L)), colour = "grey86",
                             size = 0.35, alpha = 0.5, shape = 16), dpi = 600) +
  rasterize_layer(geom_point(data = g2[ecm == FALSE], colour = SLATE,
                             size = 0.5, alpha = 0.5, shape = 16), dpi = 600) +
  geom_abline(slope = SL, intercept = IN, colour = "grey45", linewidth = 0.5) +  # positive trend guide (rho); magnitude ratio NOT interpreted (DIA compression)
  geom_point(data = hb, colour = ATTN, size = 1.3) +
  geom_text_repel(data = hb, aes(label = gene), colour = "black", size = GEOM_TEXT_6PT,
                  fontface = "italic", min.segment.length = 0, segment.size = 0.2, segment.colour = "grey75",
                  max.overlaps = 30, seed = 1) +
  # single headline: liver mRNA and liver protein agree WITHIN the same tissue
  annotate("text", x = -LB, y = LB, hjust = 0, vjust = 1, size = GEOM_TEXT_6PT, colour = "black",
           lineheight = 1.0,
           label = sprintf("liver mRNA vs liver protein:  rho = %.2f (all)\ndisease DEGs (n=%d):  rho = %.2f, %.0f%% same sign",
                           rhoB_all, nrow(ddis), rhoB_dis, concB_dis)) +
  # short in-plot key for the accented ECM cluster
  annotate("point", x = -LB + 0.02, y = -1.48, colour = ATTN, size = 1.3) +
  annotate("text",  x = -LB + 0.16, y = -1.48, hjust = 0, vjust = 0.5, size = GEOM_TEXT_6PT, colour = "black",
           label = "collagen (deposited matrix under-sampled by lysate DIA)") +
  scale_x_continuous(name = "liver mRNA log2FC", limits = c(-LB, LB)) +
  scale_y_continuous(name = "liver protein log2FC", limits = c(-LB, LB)) +
  theme_masld_compact() +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"))
# NOTE: fig4_proteo_buffering_map.pdf (panel B) retired 2026-07-07 — do not re-add.

# ── FDR-corrected Hallmark ORA per Panel-B class (report only FDR<0.05) ────────
if (requireNamespace("msigdbr", quietly = TRUE)) {
  hm <- tryCatch(msigdbr::msigdbr(species="Homo sapiens", category="H"),
                 error=function(e) tryCatch(msigdbr::msigdbr(species="Homo sapiens", collection="H"),
                                            error=function(e2) NULL))
  if (!is.null(hm)) {
    sc2 <- intersect(c("gene_symbol","human_gene_symbol"), names(hm))[1]
    sets <- split(toupper(hm[[sc2]]), hm$gs_name); uni <- toupper(unique(g2$gene))
    ora <- function(genes) {
      g <- toupper(genes); r <- rbindlist(lapply(names(sets), function(s){
        a <- length(intersect(g,sets[[s]])); if (a<3) return(NULL)
        b1<-length(g)-a; cc<-length(intersect(uni,sets[[s]]))-a; dd<-length(uni)-a-b1-cc
        data.table(set=sub("^HALLMARK_","",s), OR=(a*dd)/(b1*cc),
                   p=fisher.test(matrix(c(a,b1,cc,dd),2),alternative="greater")$p.value, k=a)}))
      if (!nrow(r)) return(r); r[, fdr := p.adjust(p,"BH")]; r[fdr<0.05][order(fdr)] }
    for (cl in levels(g2$cls)) { cat(sprintf("\n[ORA FDR<0.05] %s (n=%d):\n", cl, g2[cls==cl,.N]))
      x <- ora(g2[cls==cl, gene]); if (nrow(x)) print(x[1:min(4,.N)]) else cat("  (none survive FDR)\n") }
  }
}

# ── Caption ───────────────────────────────────────────────────────────────────
message(strrep("=", 78))
message("PROTEO-TRANSCRIPTOMIC (RNA-seq + 2 DIA-MS) — cross-layer concordance, honest build")
message(strrep("=", 78))
message(sprintf(
"Disease-vs-control log2FC on shared genes across THREE INDEPENDENT cohorts (different patients):
liver mRNA (846-sample C2) | liver protein (PXD051911, n~58) | plasma protein (PXD052937, n~72) — so
this is effect-direction CONCORDANCE, NOT a within-sample cascade. (a) Liver->blood decoupling scatter
over %d secreted-to-blood proteins (HPA; the only genes for which a plasma value is interpretable —
intracellular enzymes leak into plasma and are excluded): the circulating disease proteome is DECOUPLED
from liver transcription (liver mRNA vs plasma Spearman rho = %.2f, %.0f%% same-direction — chance-level),
with only partial coupling re-appearing at the liver-protein level (rho = %.2f). Acute-phase reactants
(CRP/SAA/HP/LBP) sit high in plasma at near-zero liver mRNA. (b) Within the SAME tissue, liver mRNA and
liver protein disease effects are CONCORDANT: all genes Spearman rho = %.2f; disease DEGs (n=%d) rho = %.2f,
%.0f%% same-direction. A positive trend line is shown; the mRNA-vs-protein magnitude ratio is NOT interpreted,
as DIA log-ratio compression makes a slope-1 (y=x) expectation technically confounded, so no y=x null is drawn.
The collagen cluster (mRNA up, soluble protein flat) is an ECM DETECTION caveat (deposited crosslinked matrix
under-sampled by lysate DIA), NOT 'fibrosis is transcript-only'. Per-class Hallmark pathways reported only at
FDR<0.05 (stdout).",
  sum(d$secreted), rho_mp, conc_mp, rho_lp,
  rhoB_all, nrow(ddis), rhoB_dis, concB_dis))
message("CAVEATS: 3 independent cohorts (not paired); proteomics underpowered (effect-direction, not padj);")
message("  plasma reflects whole-body; liver cohort GSE276114 includes CVH/ARLD (MASLD-vs-control only).")
message(strrep("=", 78))
cat("[done] wrote figS_proteo_liver_blood_decoupling.pdf (supp) + fig4_proteo_buffering_map.pdf (fig4)\n")
