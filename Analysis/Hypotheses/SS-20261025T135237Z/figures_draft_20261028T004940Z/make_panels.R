# Draft supplementary panels for SteatoSITE replication (S3c, S3d, S3e, S3g outputs).
# Reads finished outputs only; writes PDFs + source TSVs into this directory.
# Unit: donor (genetic individual after S2 QC). Staged donors only (n 493). Cross-sectional.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(edgeR); library(ggrepel)})
set.seed(20261028)
SS <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z"
OUT <- file.path(SS, "figures_draft_20261028T004940Z")
source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))

FS <- 6; GT <- FS / .pt
STAGE_LAB <- c("1" = "F0–F1", "2" = "F2", "3" = "F3", "4" = "F4")
STAGE_COL <- setNames(colorRampPalette(c("#F4A674", "#C9265E", "#650a33"))(4), STAGE_LAB)
UPC <- "#C9265E"; DNC <- "#1565C0"; NSC <- "#9E9E9E"

theme_ss <- function() theme_classic(base_size = FS, base_family = "Helvetica") %+replace% theme(
  text = element_text(size = FS, face = "plain", colour = "black"),
  plot.title = element_text(size = FS, face = "plain", hjust = 0, margin = margin(b = 3)),
  axis.text = element_text(size = FS, colour = "black"), axis.title = element_text(size = FS),
  legend.text = element_text(size = FS), legend.title = element_text(size = FS),
  axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
  legend.key.size = unit(0.25, "cm"), plot.margin = margin(4, 6, 4, 4))

save_pdf <- function(p, name, w, h) {
  f <- file.path(OUT, paste0(name, ".pdf"))
  # grDevices::pdf() cannot load encodings in this R build (share/encodings missing), so cairo_pdf is used;
  # cairo draws points as paths, so the Dingbats font issue that useDingbats = FALSE avoids does not arise.
  cairo_pdf(f, width = w, height = h, family = "Helvetica")
  print(p); invisible(dev.off()); f }

# Per-donor stage panel: box + jitter, n per stage under each box.
stage_panel <- function(dt, ycol, ylab, title, note, name, italic_title = FALSE) {
  z <- dt[!is.na(get(ycol)) & !is.na(stage_num)]
  z[, stage := factor(STAGE_LAB[as.character(stage_num)], levels = STAGE_LAB)]
  ns <- z[, .N, by = stage]; setorder(ns, stage)
  yr <- range(z[[ycol]]); pad <- diff(yr)
  ns[, y := yr[1] - 0.10 * pad]; ns[, lab := paste0("n ", N)]
  head <- (0.075 * (lengths(regmatches(note, gregexpr("\n", note))) + 1) + 0.03) * pad
  p <- ggplot(z, aes(stage, .data[[ycol]])) +
    geom_jitter(aes(colour = stage), width = 0.22, height = 0, size = 0.45, alpha = 0.6, stroke = 0) +
    geom_boxplot(fill = NA, outlier.shape = NA, width = 0.55, linewidth = 0.3, colour = "black") +
    geom_text(data = ns, aes(stage, y, label = lab), size = GT, inherit.aes = FALSE) +
    annotate("text", x = 0.5, y = yr[2] + head, label = note, hjust = 0, vjust = 1, size = GT) +
    scale_colour_manual(values = STAGE_COL, guide = "none") +
    scale_y_continuous(limits = c(yr[1] - 0.15 * pad, yr[2] + head + 0.02 * pad)) +
    labs(x = "Fibrosis stage (SteatoSITE, cross-sectional)", y = ylab, title = title) + theme_ss()
  if (italic_title) p <- p + theme(plot.title = element_text(size = FS, face = "italic", hjust = 0, margin = margin(b = 3)))
  src <- z[, c("sample_alias", "stage_num", "stage", ycol), with = FALSE]
  fwrite(src, file.path(OUT, paste0(name, "_source.tsv")), sep = "\t")
  fwrite(ns[, .(stage, n = N)], file.path(OUT, paste0(name, "_n_per_stage.tsv")), sep = "\t")
  save_pdf(p, name, 2.3, 2.4)
}
f3 <- function(x) formatC(x, format = "f", digits = 3)

## ---- 1a continuum transport ----
pu <- fread(file.path(SS, "out_s3c-20261026T233706Z/per_unit_scores.tsv"))
pu <- pu[!is.na(stage_num)]
rho1 <- cor(pu$fixed_projection, pu$stage_num, method = "spearman")
cat("1a rho", rho1, "n", nrow(pu), "\n")
stage_panel(pu, "fixed_projection", "Frozen 139-gene continuum score\n(fixed discovery loading)",
            "Discovery continuum rises with fibrosis stage",
            paste0("Spearman rho ", f3(rho1), ", n ", nrow(pu), " donors"), "p1a_continuum_by_stage")

## ---- 1b matched-gene-set null ----
nv <- fread(file.path(SS, "out_s3c-20261026T233706Z/nullvec_SteatoSITE_S1_all_staged_nullA.tsv"))
mn <- fread(file.path(SS, "out_s3c-20261026T233706Z/matched_null_results.tsv"))[stratum == "S1_all_staged" & null == "A"]
tc <- fread(file.path(SS, "out_s3c-20261026T233706Z/three_cohort_descriptive_holm.tsv"))
stopifnot(nrow(mn) == 1, abs(mn$observed_rho - rho1) < 1e-9)
ref <- tc[cohort != "SteatoSITE"]
reftxt <- paste0("Discovery-era cohorts, same score:\n",
                 paste0(ref$cohort, " rho ", f3(ref$observed_rho), " (n ", ref$n, ")", collapse = "\n"))
p1b <- ggplot(nv, aes(rho_signed)) +
  geom_histogram(bins = 40, fill = NSC, colour = NA) +
  geom_vline(xintercept = rho1, colour = UPC, linewidth = 0.5) +
  geom_vline(xintercept = mn$null_q95, colour = "black", linewidth = 0.3, linetype = "22") +
  annotate("text", x = rho1 + 0.02, y = Inf, label = paste0("observed\n", f3(rho1)), hjust = 0, vjust = 1.2, size = GT, colour = UPC) +
  annotate("text", x = rho1 + 0.02, y = Inf, label = paste0("null 95th\npercentile\n", f3(mn$null_q95)), hjust = 0, vjust = 2.6, size = GT) +
  scale_x_continuous(breaks = seq(-0.5, 1, 0.5)) + coord_cartesian(xlim = c(min(nv$rho_signed) - 0.02, 1.02)) +
  annotate("text", x = -Inf, y = Inf, label = paste0(nrow(nv), " expression-matched random gene sets\nn ", mn$n, " donors\n\n", reftxt),
           hjust = -0.03, vjust = 1.1, size = GT) +
  labs(x = "Spearman rho with fibrosis stage", y = "Random gene sets",
       title = "Frozen continuum genes beat matched random gene sets") + theme_ss()
fwrite(nv[, .(draw_index, rho_signed)], file.path(OUT, "p1b_continuum_null_source.tsv"), sep = "\t")
save_pdf(p1b, "p1b_continuum_matched_null", 3.0, 2.4)

## ---- 2 program stage slopes vs discovery ----
pr <- fread(file.path(SS, "out_s3d-20261026T233706Z/FSS5_all_programs.tsv"))[estimable == TRUE]
pr_nofrozen <- pr[!is.finite(frozen_stage_beta)]   # 1 program (Activated stellate PDGFRA) has no frozen discovery beta
pr <- pr[is.finite(frozen_stage_beta)]
r2 <- cor(pr$frozen_stage_beta, pr$beta)
cat("2 pearson", r2, "n", nrow(pr), "q<0.05", sum(pr$q_bh < 0.05), "\n")
pr[, cls := fifelse(q_bh < 0.05 & beta > 0, "Higher at later stage (BH q < 0.05)",
              fifelse(q_bh < 0.05 & beta < 0, "Lower at later stage (BH q < 0.05)", "BH q ≥ 0.05"))]
pr[, cls := factor(cls, levels = c("Higher at later stage (BH q < 0.05)", "Lower at later stage (BH q < 0.05)", "BH q ≥ 0.05"))]
pr[, lab := paste0(module_name, " [", cell_type, "]")]
top <- rbind(pr[order(-beta)][1:3], pr[order(beta)][1:3])
lim <- range(c(pr$beta, pr$frozen_stage_beta)); lim <- lim + c(-0.05, 0.05) * diff(lim)
cnt <- pr[, .N, by = cls]
p2 <- ggplot(pr, aes(frozen_stage_beta, beta)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.2, linetype = "22", colour = "grey50") +
  geom_point(aes(colour = cls, shape = cls), size = 0.9, alpha = 0.85, stroke = 0.3) +
  geom_text_repel(data = top[beta > 0], aes(label = lab), size = GT, segment.size = 0.2, min.segment.length = 0,
                  direction = "y", hjust = 0, xlim = c(lim[1], -0.02), ylim = c(0.15, lim[2]), max.overlaps = Inf,
                  box.padding = 0.2, seed = 1) +
  geom_text_repel(data = top[beta < 0], aes(label = lab), size = GT, segment.size = 0.2, min.segment.length = 0,
                  direction = "y", hjust = 0, xlim = c(0.05, lim[2]), ylim = c(lim[1], -0.15), max.overlaps = Inf,
                  box.padding = 0.2, seed = 1) +
  scale_colour_manual(values = setNames(c(UPC, DNC, NSC), levels(pr$cls)),
                      labels = paste0(levels(pr$cls), ", ", cnt$N[match(levels(pr$cls), cnt$cls)]), name = NULL) +
  scale_shape_manual(values = setNames(c(16, 17, 1), levels(pr$cls)),
                     labels = paste0(levels(pr$cls), ", ", cnt$N[match(levels(pr$cls), cnt$cls)]), name = NULL) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Discovery stage slope (frozen, per stage)", y = "SteatoSITE stage slope (per stage code)",
       title = paste0("Cell-program stage trends match discovery\nPearson r ", f3(r2), ", ", nrow(pr),
                      " programs; n ", pr$n[1], " donors")) +
  theme_ss() + theme(legend.position = "bottom", legend.direction = "vertical")
fwrite(pr[, .(feature_id, cell_type, module_name, steatosite_beta = beta, se_hc3, q_bh, frozen_stage_beta, cls)],
       file.path(OUT, "p2_program_slopes_source.tsv"), sep = "\t")
save_pdf(p2, "p2_program_stage_slopes_vs_discovery", 3.4, 4.1)

## ---- 3 IgA share vs stage ----
h <- fread(file.path(SS, "out_s3e-20261026T233706Z/h38_person_table.tsv"))
h <- h[!is.na(stage_num) & n >= 10]
h[, iga_share := k / n]
res <- jsonlite::fromJSON(file.path(SS, "out_s3e-20261026T233706Z/results.json"))$H38$T2
cat("3 b", res$b, "persons", res$persons, "nrow", nrow(h), "\n")
stage_panel(h, "iga_share", "IgA share of heavy-chain constant-region\nfragments, IgA / (IgA + IgG + IgM)",
            "IgA share of liver immunoglobulin rises with stage",
            paste0("beta-binomial b ", f3(res$b), " (SE ", f3(res$se), ")\n",
                   "per stage code | batch, sex, B-cell score\nn ", res$persons, " donors"),
            "p3_iga_share_by_stage")

## ---- 4 DIG1 genes by stage ----
us <- fread(file.path(SS, "out_s3fgh-20261026T233706Z/s3g/unit_scores.tsv"))
U <- fread(file.path(SS, "out_s2_qc/units.tsv"))[toupper(as.character(include)) == "TRUE"]
d <- readRDS(file.path(SS, "out_s2_qc/dge_included.rds"))
stopifnot(identical(colnames(d$counts), as.character(U$sample_alias)), identical(as.character(us$sample_alias), colnames(d$counts)))
deg <- fread(file.path(SS, "frozen/release/deg_results.csv"))
sym <- deg$symbol[match(rownames(d$counts), deg$gene)]
d <- calcNormFactors(d); lc1 <- cpm(d, log = TRUE, prior.count = 1)
genes <- c("TNFSF10", "TNFRSF10A", "TNFRSF10B", "GNMT", "CDKN1A")
for (g in genes) { i <- which(sym == g); stopifnot(length(i) == 1); set(us, j = paste0("l2_", g), value = lc1[i, ]) }
# check: GNMT / CDKN1A z-scores in unit_scores equal z of these log CPM
zc <- function(v) (v - mean(v)) / sd(v)
stopifnot(max(abs(zc(us$l2_GNMT) - us$j10_GNMT)) < 1e-8, max(abs(zc(us$l2_CDKN1A) - us$j10_p21)) < 1e-8)
ft <- fread(file.path(SS, "out_s3fgh-20261026T233706Z/s3g/f_ss7_tests.tsv"))
fid <- c(TNFSF10 = "s45a", TNFRSF10A = "s45c", TNFRSF10B = "s45b", GNMT = "s11c", CDKN1A = "s11b", biliary = "s11a")
note <- function(k) { r <- ft[id == fid[[k]]]
  paste0("partial Spearman rho ", f3(r$est), " | batch, sex\nn ", r$n, " donors") }
for (g in genes) stage_panel(us, paste0("l2_", g), "log2 CPM (TMM, prior 1)", g, note(g),
                             paste0("p4_", g, "_by_stage"), italic_title = TRUE)
stage_panel(us, "j10_biliary", "Biliary score (mean z: SOX9, EPCAM, KRT7,\nKRT19, CD24, SPP1, KRT23)",
            "Biliary score", note("biliary"), "p4_biliary_score_by_stage")

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
