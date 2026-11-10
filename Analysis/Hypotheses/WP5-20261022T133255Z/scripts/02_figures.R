# WP5 candidate panels (not adopted). Individual PDFs, 6 pt, no bold.
# KEY MESSAGE (A): whole-transcriptome assays observe 112-117 of the 117 fixed programs; targeted
#   CosMx panels observe 12-14 and the Xenium panel none.
# KEY MESSAGE (B): probe-based whole-transcriptome assays (Flex, GeoMx WTA) carry none of the
#   lncRNA member weight, while carrying 97-99% of the remaining weight.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))
grDevices::pdf.options(useDingbats = FALSE)
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/WP5-20261022T133255Z"
OUT <- file.path(RUN, "figures")
if (dir.exists(OUT)) stop("refusing to overwrite")
dir.create(OUT)
FS <- 6
th <- theme_classic(base_size = FS) +
  theme(text = element_text(size = FS), axis.text = element_text(size = FS, colour = "black"),
        legend.text = element_text(size = FS), legend.title = element_text(size = FS),
        plot.title = element_text(size = FS, face = "plain"), strip.text = element_text(size = FS),
        strip.background = element_blank(), axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25))

p <- fread(file.path(RUN, "outputs/observability_panel_117.tsv"))
lab <- c(GSE344087 = "GSE344087 CITE-seq", GSE298719 = "GSE298719 Flex", GSE255312 = "GSE255312 GeoMx WTA",
         GSE192741 = "GSE192741 Visium", HRA007511_HMSMA = "HMSMA Visium", Vu_et_al_2025 = "Vu 2025 Visium",
         Yakubovsky2026 = "Yakubovsky ST", Govaere2026_CosMx = "Govaere CosMx", Fukumoto2025_CosMx = "Fukumoto CosMx",
         GSE292268 = "GSE292268 CosMx", SCP3359 = "SCP3359 Xenium", Govaere2026_GeoMx = "Govaere GeoMx",
         GSE287826 = "GSE287826 GeoMx")
p <- p[dataset_id %in% names(lab)]
p[, retained := pmin(as.numeric(retained_l1_weight), 1)]  # 1 + 2e-16 rounding
p[dataset_id %in% c("Govaere2026_GeoMx", "GSE287826"), retained := NA]  # no gene axis loaded
n_obs <- p[, .(n = sum(coverage_status == "observable")), by = dataset_id]
ord_ds <- n_obs[order(-n, dataset_id)]$dataset_id
ds_lab <- setNames(paste0(lab[n_obs$dataset_id], " (", n_obs$n, ")"), n_obs$dataset_id)
p[, ds := factor(ds_lab[dataset_id], levels = ds_lab[ord_ds])]
stopifnot(!anyNA(p$ds))
po <- p[, .(m = mean(retained, na.rm = TRUE)), by = .(program_uid, cell_type)][order(cell_type, m)]
p[, prog := factor(program_uid, levels = po$program_uid)]
p[, cell_type := factor(cell_type, levels = unique(po$cell_type))]

gA <- ggplot(p, aes(ds, prog, fill = retained)) +
  geom_tile() +
  geom_point(data = p[coverage_status == "observable"], shape = 16, size = 0.15, colour = "black") +
  scale_fill_gradientn(colours = blue_gradient, limits = c(0, 1), na.value = "#E6E6E5",
                       name = "Fixed weight on\ndetected genes") +
  facet_grid(cell_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = "117 fixed programs by lineage",
       title = "Dot: observable (>= 8 genes and >= 0.20 of weight detected\nin >= 1% of cells or spots). Grey: no gene axis loaded") +
  th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
             axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y.left = element_text(angle = 0),
             strip.placement = "outside", panel.spacing = unit(0.5, "pt"))
ggsave(file.path(OUT, "panel_A_observability_117_programs.pdf"), gA, width = 4.2, height = 4.6, units = "in")

l <- fread(file.path(RUN, "outputs/lncrna_axis_presence_summary.tsv"))
l <- melt(l[, .(dataset_id, lncRNA = lncrna_on_axis_frac, other = other_on_axis_frac)], id.vars = "dataset_id",
          variable.name = "members", value.name = "frac")
l[, assay := factor(lab[dataset_id], levels = lab[c("GSE344087", "GSE298719", "GSE255312", "Fukumoto2025_CosMx", "GSE292268", "SCP3359")])]
l[, members := factor(ifelse(members == "lncRNA", "lncRNA members", "Other members"), levels = c("Other members", "lncRNA members"))]
gB <- ggplot(l, aes(frac, assay, colour = members)) +
  geom_line(aes(group = assay), colour = "#9E9E9E", linewidth = 0.25) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = c("Other members" = "#9E9E9E", "lncRNA members" = palette1[1]), name = NULL) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_y_discrete(limits = rev) +
  labs(x = "Fixed weight on the assay gene axis\n(84 programs with lncRNA members)", y = NULL,
       title = "84 programs have lncRNA members\n(lncRNA = 6.9% of their weight)") +
  th + theme(legend.position = "top")
ggsave(file.path(OUT, "panel_B_lncrna_axis_presence.pdf"), gB, width = 3.2, height = 1.9, units = "in")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
