#!/usr/bin/env Rscript
# figS_multimethod_power_design.R
# Panel B0 (VISUAL methods primer) for the multi-method power analysis: SHOWS how
# the known-truth NB simulation works, by actually simulating from the cached real
# NB params and plotting the data — so panels B/B2 (the grid sweep) make sense.
#   power_design.pdf  (was panelB0_power_design.pdf)
#   A  the simulated cohort data: 10% of genes carry a planted Disease effect
#   B  what tau2 does: cross-cohort consistency of the planted effect
#   C  one run scored vs the known truth (volcano; planted DE vs null)
# Simulator mirrors power_multimethod.R::simulate_rep EXACTLY (same constants).
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(edgeR); library(limma)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
st <- file.path(BASE, "scripts/figures/publication_theme.R")
if (file.exists(st)) try(source(st), silent = TRUE)
OUT <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

th <- if (exists("theme_masld")) function(...) theme_masld(...) else function(...) theme_bw(...)
CTRL <- "#9E9E9E"; DIS <- "#C2185B"
UP <- "#C2185B"; DOWN <- "#1565C0"; NULLC <- "#BDBDBD"

# --- real NB params used by the actual power harness ------------------------
nb <- readRDS(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/multimethod_validation/power/nb_params.rds"))
MU0 <- nb$mu; PHI0 <- nb$phi; ALL_GENES <- nb$genes
G <- length(ALL_GENES)

# --- constants (verbatim from power_multimethod.R) --------------------------
K <- 5L; PI_DE <- 0.10; PI_SEX <- 0.05; SEX_LFC <- 0.5; COHORT_OFFSET_SD <- 0.2

# --- simulator (mirrors simulate_rep) ; returns counts, meta, truth, delta ---
simulate_rep <- function(n_per_group, true_lfc, tau2, seed) {
  set.seed(seed)
  truth_g <- rbinom(G, 1, PI_DE)
  sign_g  <- sample(c(-1, 1), G, TRUE)
  lfc_g   <- truth_g * sign_g * true_lfc
  sex_gene <- rbinom(G, 1, PI_SEX)
  sexeff_g <- sex_gene * sample(c(-1, 1), G, TRUE) * SEX_LFC
  offset_k <- rnorm(K, 0, COHORT_OFFSET_SD)
  eps   <- matrix(if (tau2 > 0) rnorm(G * K, 0, sqrt(tau2)) else 0, G, K)
  delta <- (lfc_g + eps) * truth_g                       # genes x cohorts
  n_tot <- K * 2L * n_per_group
  counts <- matrix(0L, G, n_tot, dimnames = list(ALL_GENES, NULL))
  dataset <- character(n_tot); arm <- character(n_tot); col <- 0L
  for (k in seq_len(K)) for (grp in c(0L, 1L)) for (r in seq_len(n_per_group)) {
    col <- col + 1L; sex01 <- rbinom(1, 1, 0.5)
    log2mean <- log2(pmax(MU0, 1e-8)) + offset_k[k] + delta[, k] * grp + sexeff_g * sex01
    counts[, col] <- rnbinom(G, mu = 2^log2mean, size = 1 / PHI0)
    dataset[col] <- paste0("simC", k); arm[col] <- if (grp) "Disease" else "Control"
  }
  storage.mode(counts) <- "integer"
  list(counts = counts, dataset = dataset, arm = arm,
       truth = truth_g, sign = sign_g, delta = delta, lfc_g = lfc_g)
}

# ===========================================================================
# Panel A — the simulated data (one representative cell: n=20, LFC=1, tau2=0.0375)
# ===========================================================================
s <- simulate_rep(n_per_group = 20L, true_lfc = 1.0, tau2 = 0.0375, seed = 101L)
cpm_l <- edgeR::cpm(s$counts, log = TRUE, prior.count = 1)
expr  <- rowMeans(cpm_l)
de_up   <- which(s$truth == 1L & s$sign ==  1L & expr > median(expr))
de_dn   <- which(s$truth == 1L & s$sign == -1L & expr > median(expr))
nullg   <- which(s$truth == 0L & expr > median(expr))
set.seed(1)
sel <- c(head(de_up[order(-expr[de_up])], 30),
         head(de_dn[order(-expr[de_dn])], 30),
         sample(nullg, 40))
blk <- factor(c(rep("Planted UP", 30), rep("Planted DOWN", 30), rep("Null (90%)", 40)),
              levels = c("Planted UP", "Planted DOWN", "Null (90%)"))
z <- t(scale(t(cpm_l[sel, , drop = FALSE])))            # z-score per gene
z[z >  2] <- 2; z[z < -2] <- -2
# column order: cohort then sample; facet by arm
ord <- order(s$arm, s$dataset)
dA <- data.table(
  gene = factor(rep(seq_along(sel), times = ncol(z)), levels = rev(seq_along(sel))),
  block = rep(blk, times = ncol(z)),
  samp  = factor(rep(seq_len(ncol(z)), each = length(sel))),
  arm   = rep(s$arm, each = length(sel)),
  cohort = rep(s$dataset, each = length(sel)),
  z = as.vector(z))
# within-arm sample order by cohort
samp_ord <- data.table(samp = factor(seq_len(ncol(z))), arm = s$arm, cohort = s$dataset)
setorder(samp_ord, arm, cohort)
dA[, samp := factor(samp, levels = unique(samp_ord$samp))]

pA <- ggplot(dA, aes(samp, gene, fill = z)) +
  geom_raster() +
  facet_grid(block ~ arm, scales = "free", space = "free", switch = "y") +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                       midpoint = 0, limits = c(-2, 2),
                       name = "expr (z)", breaks = c(-2, 0, 2)) +
  labs(x = "samples  (5 cohorts x 20 per arm)", y = NULL) +
  th(base_size = 7) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "right",
        legend.key.height = unit(0.35, "cm"), legend.key.width = unit(0.25, "cm"),
        strip.text.y.left = element_text(angle = 0, size = 6),
        strip.text.x = element_text(face = "plain", size = 6),
        panel.spacing = unit(0.06, "lines"))

# ===========================================================================
# Panel B — what tau2 does: cross-cohort consistency of the planted effect
#   For DE genes at true_lfc=1, show the realised per-cohort effect delta_{g,k}
#   at tau2 = 0 / 0.0375 / 0.15. Tight at 0, spreads as tau2 grows.
# ===========================================================================
tau_levels <- c(0, 0.0375, 0.15)
dB <- rbindlist(lapply(tau_levels, function(t2) {
  ss <- simulate_rep(20L, 1.0, t2, seed = 202L)
  de <- which(ss$truth == 1L)
  de <- sample(de, min(120, length(de)))
  dt <- data.table(
    eff = as.vector(ss$delta[de, ]),                    # realised log2FC per cohort
    dir = rep(ifelse(ss$sign[de] > 0, "up", "down"), times = K),
    tau2 = t2)
  dt[eff != 0]
}))
dB[, tau_lab := factor(tau2, levels = tau_levels,
     labels = c("tau2 = 0\n(identical\nacross cohorts)",
                "tau2 = 0.0375\n(moderate)",
                "tau2 = 0.15\n(heterogeneous)"))]
pB <- ggplot(dB, aes(tau_lab, eff, colour = dir)) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_jitter(width = 0.22, height = 0, size = 0.35, alpha = 0.45) +
  annotate("text", x = 0.62, y = 1, label = "planted\n+1", size = GEOM_TEXT_6PT, colour = "grey35", hjust = 1, lineheight = 0.8) +
  annotate("text", x = 0.62, y = -1, label = "planted\n-1", size = GEOM_TEXT_6PT, colour = "grey35", hjust = 1, lineheight = 0.8) +
  scale_colour_manual(values = c(up = UP, down = DOWN), guide = "none") +
  coord_cartesian(ylim = c(-2.6, 2.6), clip = "off") +
  labs(x = NULL, y = "per-cohort effect (log2FC)") +
  th(base_size = 7) +
  theme(axis.text.x = element_text(size = 6, lineheight = 0.8),
        plot.margin = margin(4, 6, 4, 14))

# ===========================================================================
# Panel C — one run scored vs the known truth (fast voom+limma on the rep)
# ===========================================================================
ct <- s$counts
keep <- edgeR::filterByExpr(ct, group = factor(s$arm))
ct <- ct[keep, ]; truth_k <- s$truth[keep]
dge <- edgeR::calcNormFactors(edgeR::DGEList(ct))
des <- model.matrix(~ factor(s$dataset) + factor(s$arm, levels = c("Control", "Disease")))
v   <- limma::voom(dge, des)
fit <- limma::eBayes(limma::lmFit(v, des))
tt  <- limma::topTable(fit, coef = ncol(des), n = Inf, sort.by = "none")
dC <- data.table(logFC = tt$logFC, p = tt$P.Value, padj = tt$adj.P.Val,
                 truth = ifelse(truth_k[match(rownames(tt), rownames(ct))] == 1L,
                                "planted DE", "null"))
dC[, truth := factor(truth, levels = c("null", "planted DE"))]
called <- dC$padj < 0.05
pw  <- round(100 * mean(called[dC$truth == "planted DE"]), 0)
fdr <- if (sum(called) > 0) round(100 * sum(called & dC$truth == "null") / sum(called), 0) else NA
pthr <- max(dC$p[called], na.rm = TRUE)                  # p at BH<0.05 boundary

pC <- ggplot(dC[order(truth)], aes(logFC, -log10(p), colour = truth)) +
  geom_hline(yintercept = -log10(pthr), linetype = "dashed", colour = "#D6604D", linewidth = 0.35) +
  geom_point(size = 0.5, alpha = 0.55) +
  scale_colour_manual(values = c("null" = NULLC, "planted DE" = UP), name = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, size = GEOM_TEXT_6PT, colour = "grey20",
           label = sprintf("power = %d%% of planted DE\nFDR = %d%% of calls", pw, fdr)) +
  annotate("text", x = Inf, y = -log10(pthr), hjust = 1.05, vjust = -0.5, size = GEOM_TEXT_6PT,
           colour = "#D6604D", label = "FDR < 0.05") +
  labs(x = "estimated log2FC", y = "-log10(p)") +
  th(base_size = 7) +
  theme(legend.position = c(0.99, 0.5), legend.justification = c(1, 0.5),
        legend.text = element_text(size = 6), legend.key.size = unit(0.3, "cm"))

message("[caption] A: Simulated cohorts - 10% of genes carry a planted Disease effect; counts ~ NB anchored on real Control mean/dispersion, planted genes shift in Disease, null genes do not.")
message("[caption] B: tau2 sets cross-cohort consistency; each point = one DE gene in one cohort, spread is what metafor's RE model must absorb.")
message("[caption] C: Each run scored vs the known truth; planted-DE points recovered above the line = power, null points above = false discoveries.")

final <- pA / (pB | pC) + plot_layout(heights = c(1, 0.95))
ggsave(file.path(OUT, "power_design.pdf"), final,
       width = fig_full_width, height = 6.6 * fig_full_width / 8.8, device = cairo_pdf)
cat(sprintf("Wrote power_design.pdf  (run power=%d%% FDR=%d%%)\n", pw, fdr))
