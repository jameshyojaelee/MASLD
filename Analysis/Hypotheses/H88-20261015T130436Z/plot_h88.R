# H88 descriptive panels, read from the run outputs. Individual PDF panels, font 6, no bold, black text,
# inferred controls in gray #9E9E9E. Usage: Rscript plot_h88.R <out dir> <figure dir>
# KEY MESSAGE: per-donor excess of coordinated over non-coordinated pair co-variation across spots, against
# what random reassignment of the bulk class label within lineage pairs gives, with the zonation and portal controls beside it.
suppressPackageStartupMessages({ library(ggplot2); library(ggbeeswarm) })
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]; figdir <- args[2]
if (dir.exists(figdir) && length(list.files(figdir)) > 0) stop("figure dir not empty; nothing is overwritten")
dir.create(figdir, showWarnings = FALSE)
set.seed(20260923)
FS <- 6
CTRL_GRAY <- "#9E9E9E"
th <- theme_classic(base_size = FS, base_family = "Helvetica") +
  theme(text = element_text(size = FS, colour = "black", face = "plain"),
        axis.text = element_text(size = FS, colour = "black"), axis.title = element_text(size = FS, colour = "black"),
        legend.text = element_text(size = FS, colour = "black"), legend.title = element_text(size = FS, colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3), legend.position = "none")
save_pdf <- function(p, f, w, h) {
  grDevices::pdf(file.path(figdir, f), width = w, height = h, useDingbats = FALSE)
  print(p); invisible(grDevices::dev.off())
}
pdn <- read.delim(file.path(out, "per_donor.tsv"))
f1 <- read.delim(file.path(out, "f1_results.tsv"))
nul <- read.delim(file.path(out, "null_T1.tsv.gz"))
t1 <- f1[f1$test == "T1", ]
pdn$set <- ifelse(pdn$group == "primary", "Primary (32)", "Inferred control (3)")
cols <- c("Primary (32)" = "black", "Inferred control (3)" = CTRL_GRAY)

# Panel 1: T1 label-null distribution of D minus its mean, with the observed D_excess
nd <- data.frame(x = nul$T1_D_null - mean(nul$T1_D_null))
p1 <- ggplot(nd, aes(x)) + geom_histogram(bins = 60, fill = CTRL_GRAY, colour = NA) +
  geom_vline(xintercept = t1$D_excess, linewidth = 0.4, colour = "black") +
  labs(x = "Coordinated minus non-coordinated mean Fisher z,\nbeyond the within-stratum label null", y = "Label relabelings (of 10,000)") + th
save_pdf(p1, "H88_T1_label_null.pdf", 2.4, 1.9)

# Panel 2: per-donor e_i
p2 <- ggplot(pdn, aes(x = set, y = e_T1, colour = set)) + geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2, colour = CTRL_GRAY) +
  geom_quasirandom(size = 0.8, width = 0.25) +
  stat_summary(data = pdn[pdn$group == "primary", ], fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "black") +
  scale_colour_manual(values = cols) + labs(x = NULL, y = "Per-donor excess e_i (Fisher z)") + th
save_pdf(p2, "H88_per_donor_excess.pdf", 1.9, 2.0)

# Panel 3: positive controls per donor (PC1 zonation, PC3 portal tract)
pc <- rbind(data.frame(control = "PC1 pericentral vs periportal", z = pdn$pc1_z, set = pdn$set),
            data.frame(control = "PC3 cholangiocyte vs portal fibroblast", z = pdn$pc3_z, set = pdn$set))
pc <- pc[!is.na(pc$z), ]
p3 <- ggplot(pc, aes(x = control, y = z, colour = set)) + geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2, colour = CTRL_GRAY) +
  geom_quasirandom(size = 0.8, width = 0.25) + scale_colour_manual(values = cols) +
  scale_x_discrete(labels = function(x) gsub(" vs ", "\nvs ", x)) + labs(x = NULL, y = "Spot-level Fisher z per donor") + th
save_pdf(p3, "H88_positive_controls.pdf", 2.2, 2.0)

# Panel 4: e_i against spots per donor
p4 <- ggplot(pdn, aes(x = n_spots, y = e_T1, colour = set)) + geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2, colour = CTRL_GRAY) +
  geom_point(size = 0.8) + scale_colour_manual(values = cols) + labs(x = "Spots per array", y = "Per-donor excess e_i (Fisher z)") + th
save_pdf(p4, "H88_excess_vs_spots.pdf", 2.2, 2.0)
message("Caption: e_i = donor i's coordinated minus non-coordinated mean Fisher z across spots, minus its mean over 10,000 within-stratum ",
        "relabelings of the bulk class (T1). Black, 32 primary arrays; gray, 3 inferred controls (not in the test).")
writeLines(capture.output(sessionInfo()), file.path(figdir, "sessionInfo.txt"))
