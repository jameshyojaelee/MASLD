# WP5 candidate panel (not adopted). Individual PDF, 6 pt, no bold.
# KEY MESSAGE: whole-transcriptome assays observe 112-117 of the 117 fixed programs; targeted
#   CosMx panels observe 12-14 and the Xenium panel none.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))
grDevices::pdf.options(useDingbats = FALSE)
RUN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/WP5-20261022T133255Z"
OUT <- file.path(RUN, "figures_bar")
if (dir.exists(OUT)) stop("refusing to overwrite")
dir.create(OUT)
FS <- 6
th <- theme_classic(base_size = FS) +
  theme(text = element_text(size = FS), axis.text = element_text(size = FS, colour = "black"),
        legend.text = element_text(size = FS), legend.title = element_text(size = FS),
        plot.title = element_text(size = FS, face = "plain"), plot.title.position = "plot",
        axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25))

s <- fread(file.path(RUN, "outputs/observability_summary.tsv"))
lab <- c(GSE344087 = "GSE344087 CITE-seq", GSE298719 = "GSE298719 Flex", GSE255312 = "GSE255312 GeoMx WTA",
         GSE192741 = "GSE192741 Visium", HRA007511_HMSMA = "HMSMA Visium", Vu_et_al_2025 = "Vu 2025 Visium",
         Yakubovsky2026 = "Yakubovsky ST", Govaere2026_CosMx = "Govaere CosMx", Fukumoto2025_CosMx = "Fukumoto CosMx",
         GSE292268 = "GSE292268 CosMx", SCP3359 = "SCP3359 Xenium",
         Govaere2026_GeoMx = "Govaere GeoMx (no gene list)", GSE287826 = "GSE287826 GeoMx (no gene list)")
s <- s[dataset_id %in% names(lab)]
stopifnot(nrow(s) == 13, all(s$observable + s$partial + s$untestable == 117))
s[, ds := factor(lab[dataset_id], levels = lab[s[order(observable, partial), dataset_id]])]
l <- melt(s[, .(ds, observable, partial, untestable)], id.vars = "ds", variable.name = "state", value.name = "n")
l[, state := factor(state, levels = c("untestable", "partial", "observable"),
                    labels = c("Untestable", "Partial", "Observable"))]

g <- ggplot(l, aes(n, ds, fill = state)) +
  geom_col(width = 0.75) +
  geom_text(data = s, aes(x = 117, y = ds, label = observable), inherit.aes = FALSE,
            hjust = -0.4, size = FS / .pt) +
  scale_fill_manual(values = c(Observable = blue_gradient[10], Partial = blue_gradient[3], Untestable = "#E6E6E5"),
                    breaks = c("Observable", "Partial", "Untestable"), name = NULL) +
  scale_x_continuous(limits = c(0, 132), breaks = c(0, 39, 78, 117), expand = c(0, 0)) +
  labs(x = "Fixed programs (of 117); number = observable", y = NULL,
       title = "Observable: >= 8 genes and >= 0.20 of\nprogram weight detected in >= 1% of cells") +
  th + theme(legend.position = "top", legend.key.size = unit(6, "pt"))
ggsave(file.path(OUT, "panel_observable_programs_bar.pdf"), g, width = 3.0, height = 2.4, units = "in")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
